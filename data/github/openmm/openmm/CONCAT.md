[![GH Actions Status](https://github.com/openmm/openmm/workflows/CI/badge.svg)](https://github.com/openmm/openmm/actions?query=branch%3Amaster+workflow%3ACI)
[![Conda](https://img.shields.io/conda/v/conda-forge/openmm.svg)](https://anaconda.org/conda-forge/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/conda-forge/openmm/badges/downloads.svg)](https://anaconda.org/conda-forge/openmm)

## OpenMM: A High Performance Molecular Dynamics Library

Introduction
------------

[OpenMM](http://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. It
provides a combination of extreme flexibility (through custom forces and integrators), openness, and high performance (especially on recent GPUs) that make it truly unique among simulation codes.  

Getting Help
------------

Need Help? Check out the [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).
## How to Get Support for OpenMM

There are two main venues for getting support for OpenMM: the [discussion forum](https://simtk.org/forums/viewforum.php?f=161)
and the [Github repository](https://github.com/openmm/openmm).  There is some overlap
between the two, but generally speaking the forum is for user oriented issues while the
repository is for developer oriented issues.  If you have a question about how to use OpenMM
(including writing programs that access it through its public API), post on the forum.  If
you want to suggest a change to the code, or if you think you have found a bug,
open an issue on Github.  The core developers monitor both, so don't worry if you aren't
sure which one is most appropriate for your question.  We will see it either way.

You also may want to consult the [documentation](http://docs.openmm.org/).  It is quite
thorough, and you may be able to find the answer to your question.# Contributor Covenant Code of Conduct

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
reported by contacting the project team at 'peastman@stanford.edu'. The project team will
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
## How to Contribute to OpenMM Development

We welcome anyone who wants to contribute to the project, whether by adding a feature,
fixing a bug, or improving documentation.  The process is quite simple.

First, it is always best to begin by opening an issue on Github that describes the change you
want to make.  This gives everyone a chance to discuss it before you put in a lot of work.
For bug fixes, we will confirm that the behavior is actually a bug and that the proposed fix
is correct.  For new features, we will decide whether the proposed feature is something we
want and discuss possible designs for it.

Once everyone is in agreement, the next step is to
[create a pull request](https://help.github.com/en/articles/about-pull-requests) with the code changes.
For larger features, feel free to create the pull request even before the implementaton is
finished so as to get early feedback on the code.  When doing this, put the letters "WIP" at
the start of the title of the pull request to indicate it is still a work in progress.

For new features, consult the [New Feature Checklist](https://github.com/openmm/openmm/wiki/Checklist-when-adding-a-new-feature),
which lists various items that need to be included before the feature can be merged (documentation,
tests, serialization, support for all APIs, etc.).  Not every item is necessarily applicable to
every new feature, but usually at least some of them are.

The core developers will review the pull request and may suggest changes.  Simply push the
changes to the branch that is being pulled from, and they will automatically be added to the
pull request.  In addition, the full test suite is automatically run on every pull request,
and rerun every time a change is added.  Once the tests are passing and everyone is satisfied
with the code, the pull request will be merged.  Congratulations on a succesful contribution!# Conda Forge releases for OpenMM

## Final releases

Create a new version tag on `openmm/openmm` and publish a new GitHub release. The Conda Forge bots would auto submit a PR with the new version within a couple hours _if_ we were using the `source.url` field pointing to a GH release. Since we build off the full repo, I am not sure the automation will work, so we need to switch to manual mode.

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. On your `conda-forge/openmm-feedstock` **fork**, create a new branch off most recent upstream branch, be it `master` or `rc`.
3. Edit `meta.yaml`:
   - [ ] Update `package.version`.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.git_rev` points to the release tag you want to publish. This could be a git commit too, but a tag is preferred.
4. Commit and push to **your fork**. Do NOT push to `upstream`.
5. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targetting `conda-forge/openmm-feedstock`'s `master`, from **your fork**.
6. Review the checklist and open the PR.
7. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
8. Wait for all green, reviews and then merge. Always make sure you are merging to `master`.
9. Once merged, check the CI status on the `master` branch.

   - It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
   - If a CI provider was not triggered (for whatever reason), `master` might need a little _push_ (no pun intended). An empty commit to `master` will do:

   ```
   # make sure conda-forge/openmm-feedstock is configured as `upstream`
   git remote -v
   git checkout master
   git fetch upstream master
   git merge upstream/master
   git commit --allow-empty -m "Trigger CI"
   git push upstream master
   ```

## Release candidates

> Technically, once you publish an RC tag on GitHub Releases, the bots will pick it up, but we haven't tested this yet.

Manual instructions:

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. Create a new branch from the most recent upstream branch, be it `master` or `rc`.
3. Edit `meta.yaml`:
   - [ ] Update `package.version`. It should be the new version number plus `rcX`, `X` being a number. Check [CFEP05](https://github.com/conda-forge/cfep/blob/master/cfep-05.md) in case of doubt.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.git_rev` points to the release tag you want to publish. This could be a git commit too, but a tag is preferred.
4. Edit `conda_build_config.yaml`. This a KEY difference: RC releases are published to a different label!

   - [ ] `channel_targets` should be set to `- conda-forge openmm_rc`:

   ```yaml
   channel_targets:
     - conda-forge openmm_rc
   ```

5. Commit and push to **your fork**. Do NOT push to `upstream`.
6. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targetting `conda-forge/openmm-feedstock`'s `rc`, from **your fork**. Again, we are **targetting** the `rc` branch, NOT master. This is a KEY difference. RC candidates stay on `rc`.
7. Review the checklist and open the PR.
8. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
9. Wait for all green, reviews and then merge. Always make sure you are merging to `rc`.
10. Once merged, check the CI status on the `rc` branch.

- It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
- If a CI provider was not triggered (for whatever reason), `rc` might need a little _push_ (no pun intended). An empty commit to `rc` will do:

```
# make sure conda-forge/openmm-feedstock is configured as `upstream`
git remote -v
git checkout rc
git fetch upstream rc
git merge upstream/rc
git commit --allow-empty -m "Trigger CI"
git push upstream rc
```
<!-- Authored by Jaime Rodríguez-Guerra, Chodera Lab. December 2020 -->

# Our Continuous Integration setup

OpenMM can be described as a C++ library with wrappers available in different programming languages (Python, C, Fortran). The heavy lifting is performed by the backend platforms, which can be based on CPU, CUDA and/or OpenCL (and possibly more in the future). All of this is supported for different operating systems and architectures. As a result, the CI setup can get a bit involved, but this document will try to clarify how it works and what we support.

## Implementation overview

OpenMM's CI runs mainly on GitHub Actions, with one separate Jenkins box running the GPU tests (generously provided by Jason Swails).

The build matrix covers:

- Operating systems and architecture:
  - Linux x64
  - MacOS Intel
  - Windows
  - Linux ppc64le (PowerPC)
  - Linux aarch64 (ARM)
- Python
  - CPython 3.6, 3.7, 3.8, 3.9
- CUDA versions
  - 10.0 and above (Linux x64, Linux ppc64le, Windows)
- OpenCL implementations
  - Nvidia (tested along CUDA)
  - AMD 3.0
- Sysroots and C++ Compilers
  - Linux: System's GCC 7 and whatever conda-forge is pinning (GCC 9 as of writing)
  - MacOS: System's, targetting 10.9 SDK
  - Windows: VS2019

Before I describe the pipelines, I will clarify some concepts and idiosyncrasies in GitHub Actions

- The configuration file lives on `.github/workflows/CI.yml`. This directory can host more than one YML _workflow_, each describing a set of event that will trigger a run.
- The workflow specifies a set of triggers (key `on`) and a list of `jobs` to run. We run the `CI` workflow for:
  - Pushes to `master`
  - Pull requests targetting `master`
  - Nightlies
- Currently, the workflow contains four jobs: `unix`, `windows`, `docker`, `docs`. Each job can be run several times, depending on the configuration of `jobs.*.strategy.matrix`. All those jobs replicas will run in parallel and individually. The [`Actions > Summary`](https://github.com/openmm/openmm/actions/runs/451301350) overview can help visualize this.
- Within each job, you find `steps`. A step can either run a script on a `shell` or use a GitHub `action` to perform a task.
  - For example, cloning the repo or setting up Miniconda are both independent GitHub _actions_. You will recognize this because they contain the keyword `uses:`.
  - Running CMake is a shell step, which uses `run:`.
  - Note 1: Each step is run a new shell session. Environment variables won't survive across steps, unless you add them to the `$GITHUB_ENV` file: `echo "VARIABLE=VALUE" >> ${GITHUB_ENV}`. You can also use step `outputs` but that's more involved and rarely needed.
  - Note 2: Due to the design of `conda-incubator/setup-miniconda`, all subsequent steps that rely on a conda environment require us to specify an OS-dependent custom shell. Do remember this if you need to add more steps in the job!
- Steps can be run or skipped based on conditions expressed inside an `if:` key. This is how we control whether we need to install CUDA or not, for example. Jobs can have `if` check, if needed.
- Steps can define environment variables in their `env:` key, but they will only be available in that step. A `job` can do it too, and these will be available for all steps.

## Details per operating system

The different implementations are very similar to what we do on Linux x64, so I will explain this one on detail and the rest will only comment on the relevant differences.

### Linux x64

- Part of the `unix` pipeline.
- Runs on `ubuntu-latest`, as provided by GitHub Actions.
- Uses `conda-incubator/setup-miniconda` to setup the bundled Miniconda and install a conda environment available that provides the building and testing dependencies (CMake, Swig, the adequate Python version, etc). These environment files are located under `devtools/ci/gh-actions/conda-envs`, per operating system.
- Depending on the matrix configuration, we also install CUDA and/or AMD's OpenCL. These conditional steps are evaluated using GHA's builtin `if` mechanism. Ideally we would install this within the conda environment, but sometimes they are not available (licensing issues, etc(), so we delegate that to the system packages or vendor installers.
  - For CUDA, we check whether `cuda-version` is not empty, and pass it to `devtools/ci/gh-actions/scripts/install_cuda.sh` as an environment variable.
  - For OpenCL, we check whether `OPENCL` is `true` and run `devtools/ci/gh-actions/scripts/install_amd_opencl.sh`. This relies on a installer located in a S3 bucket. This could be refactored to install different OpenCL implementations (ROCm, Intel, etc).
- Some matrix entries require us to install the conda forge compilers, which are used instead of the system's if present.
- Now we need to configure and download the CCache contents. The keys are built off the matrix name, and a `YYYYDDMM-HHMMSS` timestamp. A secret `CACHE_VERSION` is also included so one can bump the cache by modifying this secret in the repository settings. The configuration is done through environment variables defined at the beginning of the job (key `jobs.unix.env`).
- CMake is finally invoked, targetting the conda environment as destination (`CONDA_PREFIX`). Additional flags are passed from the matrix configuration. This is how we enable or disable features per matrix entry.
- CCache performance is assessed.
- Then we build the C++ libraries and Python wrappers, but separately. This way we can visually check which part failed more easily. Tests are also run separately for the same reason. Whether Python is built and/or tested is checked through the contents of `CMAKE_FLAGS`.

### MacOS Intel

- Part of the `unix` pipeline.
- Runs on `macos-latest`.
- Uses `conda-incubator/setup-miniconda`, pointing to the relevant environment file.
- Neither CUDA nor OpenCL installation scripts are run. Instead, we download and install the 10.9 SDK using `devtools/ci/gh-actions/scripts/install_macos_sdk.sh`. This is done so we can mimic what Conda Forge does in their feedstocks. Check the scripts comments for more info.
- Everything else is the same.

### Windows

- Sole member of the `windows` pipeline.
- Runs on `windows-latest`.
- Uses `conda-incubator/setup-miniconda`, pointing to the relevant environment file.
- Installs CUDA with the Nvidia installers using `devtools/ci/gh-actions/scripts/install_cuda.bat`, which requires an environment variable `CUDA_VERSION`, exported from the corresponding matrix entry. Again, this only runs if `matrix.cuda-version` is not empty.
- Everything else is the same.

### PowerPC & ARM

- Part of the `docker` pipeline.
- These run on a Docker image on top of `ubuntu-latest`. The Docker image itself depends on the architecture chosen (ppc64le, aarch64) and what CUDA version we want. These are provided by Conda Forge, so they have `conda` preinstalled and ready to go.
- Since it's a different architecture, we need to configure QEMU first. This is done automatically with a Docker image, mimicking what Conda Forge does.
- We start the Docker image. The working directory (`$GITHUB_WORKSPACE`) is mounted with read/write permissions on `/home/conda/workspace`, so we can communicate back with the host using files, and also use CCache.
- The Docker image will run `devtools/ci/gh-actions/scripts/run_steps_inside_docker_image.sh`. This script mostly does what you saw for Linux x64, with some differences:
  - We don't need to install CUDA or setup Miniconda, because they are preinstalled in the Docker image.
  - We patch some dependencies from the environment file because they are not available for this architecture. To save one conda environment solve, we also patch the Python version in the environment file.
  - These images don't come with a system compiler, so we specify one in the matrix configuration:
    - If `compilers` contains a value that starts with `devtoolset-`, we understand we want a CentOS devtoolse. So far, we specify `devtoolset-7`.
    - If `compilers` is any other thing, we understand that's a (space-separated series of) Conda packages. Since Conda Forge provides a metapackage named `compilers` that will install all of them for the current platform, we use that one. That's why some entries have a `compilers: compilers` entry.
  - Everything else runs as usual.
- Do note that the whole Docker run is a single GitHub Actions step, so it's not as visually appealing. I tried my best to group the commands with the `::group::` syntax so it's easier to follow, but it's not the same.
- If the script runs successfully, it will create an empty file. We test for existence after the Docker run to make sure.

> Note: Since these use software emulation, they are really slow. Still, they can run successfully within the 6h GHA provides. If GHA upgrades to better CI machines with hardware based virtualization, they might be able to run with close-to-native performance.

### Docs

This is a Linux-x64 pipeline optimized for building the documentation only. It's provided as a separate entry because I didn't want to overcomplicate the `if:` logic in the `unix` pipeline. It's essentially the same, but:

- It uses a different environment file in `setup-miniconda`.
- It only builds the docs, and their dependencies. No tests, for example.
- It contains a deployment step, which will copy the contents to the S3 bucket _only_ when run on `master`, ignoring cron jobs. The required secrets must be defined in the repository settings with the following exact key names. Just copy paste the values there. GitHub will encrypt and mask them.
  - `AWS_S3_BUCKET`
  - `AWS_ACCESS_KEY_ID`
  - `AWS_SECRET_ACCESS_KEY`
- It will also check for dead links using a Node package. This is run _after_ deployment so it won't prevent that, but it will still signal the job as failed if the docs contain broken links.

## Shortcomings

There are some limitations when compared to other CI services, but I guess this list will be shorter over time:

- Cache cannot be invalidated directly. Instead, I included a secret `CACHE_VERSION` that is part of the cache key. If you change the value of this secret, it will functionally prevent access to the previous cache. It also expires every 7 days. Note that since this trick uses a secret, the value of `CACHE_VERSION` will be masked in the log output. As a result, make sure to use something short but meaningless and difficult to find in the wild (e.g. `pqgbhl` instead of `0`).
- There's no `ci skip` functionality (yet).

## Extra content

### How to debug PowerPC / ARM locally

From the root of the repository, run the following script. There are
some variables you might want to edit (PPC vs ARM, Python version, etc).
Take a look to the script first in that case.

```bash
bash devtools/ci/gh-actions/start_docker_locally.sh
```

You will be inside the Docker image after a few moments. The repo root has
been mounted to `/home/conda/workspace`.

Run this other script to reproduce the CI steps exactly. Do NOT `source` scripts,
since a failure will exit Docker altogether. Always use new `bash` processes
to avoid starting from scratch.

```bash
bash /home/conda/workspace/devtools/ci/gh-actions/scripts/run_steps_inside_docker_image.sh
```
# Packaging OpenMM into ZIP installers

Set your environment variable `TAG` to the git tag for the release:
```bash
# OpenMM 7.1.1
export TAG="c1a64aa"
```

## Source

Start the docker container:
```bash
docker run -i -t --rm -e TAG -v `pwd`:/io jchodera/omnia-build-box:cuda80-amd30-clang38 bash
```
Patch the docker container for missing LaTeX files:
```
tlmgr install fncychap tabulary capt-of eqparbox environ trimspaces
```
Build the installer inside the docker container:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/source/prepare.sh       
source openmm/devtools/packaging/scripts/source/build.sh
source openmm/devtools/packaging/scripts/source/package.sh
# Recover the packages to host directory
cp packaging/compressed/* /io
```

## Linux

Start the docker container:
```bash
docker run -i -t --rm -e TAG -v `pwd`:/io jchodera/omnia-build-box:cuda80-amd30-clang38 bash
```
Patch the docker container for missing LaTeX files:
```
tlmgr install fncychap tabulary capt-of eqparbox environ trimspaces
```
Build the installer inside the docker container:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/linux/prepare.sh
source openmm/devtools/packaging/scripts/linux/build.sh
source openmm/devtools/packaging/scripts/linux/package.sh
# Recover the packages to host directory
cp packaging/compressed/* /io
```

## OS X

On an `osx` machine with XCode and the OS X 10.9 frameworks installed:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/osx/prepare.sh
source openmm/devtools/packaging/scripts/osx/build.sh
source openmm/devtools/packaging/scripts/osx/package.sh
```
# Manifests for automated packaging of source and binary distributions

A detailed explanation of packaging protocols can be found on the developer wiki:
https://github.com/pandegroup/openmm/wiki/Packaging-OpenMM-installers

## Contents
* `source/` - directories and files to be copied from the GitHub repo for a source distribution
* `binary/` - directories and files to be copied from install directory after build for a binary distribution

=============
Extra classes
=============

Tabulated functions
===================

These classes use table of values to define a mathematical function and can be
used by various :ref:`custom forces <custom-forces>`.
The :ref:`OpenMM::TabulatedFunction` class is an abstract class that the other classes
extend.

.. toctree::
    :maxdepth: 2

    generated/TabulatedFunction

    generated/Continuous1DFunction
    generated/Continuous2DFunction
    generated/Continuous3DFunction
    generated/Discrete1DFunction
    generated/Discrete2DFunction
    generated/Discrete3DFunction

Virtual Sites
=============

A virtual site is a particle whose position is computed directly from the
positions of other particles. The :ref:`OpenMM::VirtualSite` class is an abstract
class that the other classes extend.

.. toctree::
    :maxdepth: 2

    generated/VirtualSite

    generated/LocalCoordinatesSite
    generated/OutOfPlaneSite
    generated/ThreeParticleAverageSite
    generated/TwoParticleAverageSite

Serialization
=============

These classes are used to serialize other objects, allowing them to be stored on
disk.

.. toctree::
    :maxdepth: 2

    generated/SerializationNode
    generated/SerializationProxy
    generated/XmlSerializer

Other classes
=============

These classes don't fit neatly into the other categories, but that is not to say
that they aren't important!

.. toctree::
    :maxdepth: 2

    generated/LocalEnergyMinimizer
    generated/NoseHooverChain
    generated/OpenMMException
    generated/Vec3
==============
OpenMM C++ API
==============

The C++ API provides information about the classes and methods available in OpenMM for C++ developers. OpenMM uses an object-oriented API that makes all its functionality available through a small number of classes.

Core classes
============

.. toctree::
    :maxdepth: 1
    :hidden:

    generated/System
    generated/Context
    generated/State
    generated/Platform

:cpp:class:`OpenMM::System`
---------------------------

A ``System`` specifies generic properties of the molecular system to be
simulated: the number of particles it contains, the mass of each one, the size
of the periodic box, and so on.  The interactions between the particles are
specified through a set of :ref:`Force <forces>` objects that are added to the
``System``.  Force field specific parameters, such as particle charges, are
stored in these ``Force`` objects, not as direct properties of the ``System``.

:cpp:class:`OpenMM::Context`
----------------------------

A ``Context`` stores all of the state information for a simulation: particle
positions and velocities, as well as arbitrary parameters defined by the
``Forces`` in the System.  It is possible to create multiple ``Contexts`` for a
single ``System``, and thus have multiple simulations of that ``System`` in
progress at the same time. ``Context`` does not provide methods for accessing
state variables directly; they must be read via a ``State`` object.


:cpp:class:`OpenMM::State`
--------------------------

A ``State`` object must be constructed before data can be read from a
simulation. State variables are not accessible directly via a ``Context`` in
order to make explicit the precise time that a variable reflects. A ``State``
is created by calling a method on a ``Context`` and stores only the information
requested at invocation.

:cpp:class:`OpenMM::Platform`
-----------------------------

A ``Platform`` is a single implementation of OpenMM at a low level. This allows
the same high level API documented here to be used on all sorts of compute
hardware, from GPUs to supercomputers. A ``Platform`` implements some set of
kernels, which define which operations it supports. Writing a new ``Platform``
allows OpenMM to be ported to new hardware or to be implemented in a new way
without rewriting the entire application.

Forces
======

``Force`` objects define the behavior of the particles in a ``System``. The
``Force`` class is actually slightly more general than its name suggests.  A
``Force`` can, indeed, apply forces to particles, but it can also directly
modify particle positions and velocities in arbitrary ways.  Some thermostats
and barostats, for example, can be implemented as ``Force`` classes.  Examples
of Force subclasses include :cpp:class:`HarmonicBondForce
<OpenMM::HarmonicBondForce>`, :cpp:class:`NonbondedForce
<OpenMM::NonbondedForce>`, and :cpp:class:`MonteCarloBarostat
<OpenMM::MonteCarloBarostat>`.

.. toctree::
    :maxdepth: 2

    forces

Integrators
===========

An ``Integrator`` implements an algorithm for advancing the simulation through
time.  They provide a ``Context`` a means of stepping the simulation forward,
and must be coupled to a ``Context`` to function. Examples of Integrator
subclasses include :cpp:class:`LangevinIntegrator <OpenMM::LangevinIntegrator>`,
:cpp:class:`VerletIntegrator <OpenMM::VerletIntegrator>`, and :cpp:class:`BrownianIntegrator <OpenMM::BrownianIntegrator>`.

.. toctree::
    :maxdepth: 2

    integrators

Extras
======

OpenMM's public API includes a few more classes that support the above.

.. toctree::
    :maxdepth: 2

    extras
.. _forces:

======
Forces
======

The ``Force`` abstract class
============================

The ``Force`` objects added to a ``System`` define the behavior of the
particles. ``Force`` is an abstract class; subclasses implement specific behaviors. Classes that extend ``Force`` may implement actual physical forces, or any number of processes that either actually apply forces to particles or directly modify their positions or momenta.

.. toctree::
    :maxdepth: 2

    generated/Force


Common bonded and non-bonded forces
===================================

These classes implement forces that are widely used in biomolecular simulation.

.. toctree::
    :maxdepth: 2

    generated/CMAPTorsionForce
    generated/DrudeForce
    generated/GBSAOBCForce
    generated/GayBerneForce
    generated/HarmonicAngleForce
    generated/HarmonicBondForce
    generated/NonbondedForce
    generated/PeriodicTorsionForce
    generated/RBTorsionForce


AMOEBA forces
=============

These forces are used to implement the polarizable AMOEBA force fields.

.. toctree::
    :maxdepth: 2

    generated/AmoebaGeneralizedKirkwoodForce
    generated/AmoebaMultipoleForce
    generated/AmoebaTorsionTorsionForce
    generated/AmoebaVdwForce
    generated/AmoebaWcaDispersionForce
    generated/HippoNonbondedForce


Pseudo-forces
=============

These inherit from ``Force``, but do not describe physical forces. They are used
to implement thermostats or barostats, or otherwise modify the simulation from
step to step. They are conceptually closer to modifications to the integrator,
but providing them as a ``Force`` simplifies implementation and allows them to
be combined in arbitrary ways.

.. toctree::
    :maxdepth: 2

    generated/AndersenThermostat
    generated/CMMotionRemover
    generated/MonteCarloAnisotropicBarostat
    generated/MonteCarloBarostat
    generated/MonteCarloFlexibleBarostat
    generated/MonteCarloMembraneBarostat
    generated/RMSDForce
    generated/RPMDMonteCarloBarostat


.. _custom-forces:

Customizing ``Force``
=====================

OpenMM provides a number of classes that make it easier to implement custom
forces for common scenarios. These classes implement constructors that take an
algebraic expression as a string. The class is instantiated (not extended) to
provide a ``Force`` object that efficiently implements the provided
expression.

.. toctree::
    :maxdepth: 2

    generated/CustomAngleForce
    generated/CustomBondForce
    generated/CustomCVForce
    generated/CustomCentroidBondForce
    generated/CustomCompoundBondForce
    generated/CustomExternalForce
    generated/CustomGBForce
    generated/CustomHbondForce
    generated/CustomManyParticleForce
    generated/CustomNonbondedForce
    generated/CustomTorsionForce
===========
Integrators
===========

The ``Integrator`` abstract class
=================================

An ``Integrator`` implements an algorithm for advancing the simulation through
time.  ``Integrator`` is an abstract class; subclasses implement specific
algorithms.

.. toctree::
    :maxdepth: 2

    generated/Integrator


General purpose integrators
===========================

These are integrators appropriate for traditional MD and BD simulations.

.. toctree::
    :maxdepth: 2

    generated/BrownianIntegrator
    generated/LangevinIntegrator
    generated/LangevinMiddleIntegrator
    generated/NoseHooverIntegrator
    generated/VariableLangevinIntegrator
    generated/VariableVerletIntegrator
    generated/VerletIntegrator


Drude integrators
=================

These integrators permit modelling polarization with a Drude particle.

.. toctree::
    :maxdepth: 2

    generated/DrudeIntegrator
    generated/DrudeLangevinIntegrator
    generated/DrudeNoseHooverIntegrator
    generated/DrudeSCFIntegrator


Ring Polymer Molecular Dynamics integrators
===========================================

The RPMD integrator implements Ring Polymer MD.

.. toctree::
    :maxdepth: 2

    generated/RPMDIntegrator


Customizing ``Integrator``
==========================

These classes facilitate customisation of the integrator. ``CustomIntegrator``
allows a wide variety of integration algorithms to be implemented efficiently
without writing any low-level code. The integrator is built up as a series of
steps, each defined as an algebraic expression. ``CompoundIntegrator`` allows
different integrators to be combined by making it possible to switch the active
integrator in the middle of a simulation.

.. toctree::
    :maxdepth: 2

    generated/CustomIntegrator
    generated/CompoundIntegrator
.. only:: html

   Bibliography
   ############

.. bibliography:: references.bib
   :style: unsrt
.. _the-openmm-library:

###########################
Part II: The OpenMM Library
###########################


.. toctree::
   :numbered: 3
   :maxdepth: 3

   library/01_introduction
   library/02_compiling
   library/03_tutorials
   library/04_platform_specifics
   library/05_languages_not_cpp
   library/06_integration_examples
   library/07_testing_validation
   library/08_amoeba_plugin
   library/09_rpmd_plugin
   library/10_drude_plugin
#####################################
OpenMM User's Manual and Theory Guide
#####################################

.. only:: latex

   .. include:: license.rst

.. toctree::
   :numbered:
   :maxdepth: 3

   introduction

.. toctree::
   :maxdepth: 3

   application
   library
   theory

.. toctree::
   :numbered:
   :maxdepth: 3

   zbibliography

.. only:: html

   .. include:: license.rst
Introduction
############

OpenMM consists of two parts:

#. A set of libraries that lets programmers easily add molecular simulation
   features to their programs
#. An “application layer” that exposes those features to end users who just want
   to run simulations


This guide is divided into three sections:

* :ref:`Part I <the-openmm-application-layer>`
  describes the application layer.  It is relevant to all users, but especially relevant to people
  who want to use OpenMM as a stand-alone application for running simulations.
* :ref:`Part II <the-openmm-library>`
  describes how to use the OpenMM libraries within your own applications.  It is primarily
  relevant to programmers who want to write simulation applications.
* :ref:`Part III <the-theory-behind-openmm>`
  describes the mathematical theory behind the features found in OpenMM.  It is relevant to all users.


Online Resources
****************

You can find more documentation and other material at our website
http://openmm.org.   Among other things there is a discussion forum,
a wiki, and videos of lectures on using OpenMM.


Referencing OpenMM
******************

Any work that uses OpenMM should cite the following publication:

P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp,
L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora,
B. R. Brooks, and V. S. Pande. "OpenMM 7: Rapid development of high performance
algorithms for molecular dynamics." PLOS Comp. Biol. 13(7): e1005659. (2017)

We depend on academic research grants to fund the OpenMM development efforts;
citations of our publication will help demonstrate the value of OpenMM.


Acknowledgments
***************

OpenMM software and all related activities, such as this manual, are funded by
the Simbios National Center for Biomedical Computing through the National
Institutes of Health Roadmap for Medical Research, Grant U54 GM072970, and by
National Institutes of Health grant R01-GM062868.

.. default-domain:: py

.. _the-openmm-application-layer:

####################################
Part I: The OpenMM Application Layer
####################################

.. toctree::
   :numbered: 3
   :maxdepth: 3

   application/01_getting_started
   application/02_running_sims
   application/03_model_building_editing
   application/04_advanced_sim_examples
   application/05_creating_ffs
.. _the-theory-behind-openmm:

##################################
Part III: The Theory Behind OpenMM
##################################

.. toctree::
   :numbered: 3
   :maxdepth: 3

   theory/01_introduction
   theory/02_standard_forces
   theory/03_custom_forces
   theory/04_integrators
   theory/05_other_features

.. _other-features:

Other Features
##############


Periodic Boundary Conditions
****************************

Many Force objects support periodic boundary conditions.  They act as if space
were tiled with infinitely repeating copies of the system, then compute the
forces acting on a single copy based on the infinite periodic copies.  In most
(but not all) cases, they apply a cutoff so that each particle only interacts
with a single copy of each other particle.

OpenMM supports triclinic periodic boxes.  This means the periodicity is defined
by three vectors, :math:`\mathbf{a}`\ , :math:`\mathbf{b}`\ , and
:math:`\mathbf{c}`\ .  Given a particle position, the infinite periodic copies
of that particle are generated by adding vectors of the form
:math:`i \mathbf{a}+j \mathbf{b}+k \mathbf{c}`\ , where :math:`i`\ ,
:math:`j`\ , and :math:`k` are arbitrary integers.

The periodic box vectors must be chosen to satisfy certain requirements.
Roughly speaking, :math:`\mathbf{a}`\ , :math:`\mathbf{b}`\ , and
:math:`\mathbf{c}` need to "mostly" correspond to the x, y, and z axes.  They
must have the form

.. math::
   \mathbf{a} = (a_x, 0, 0)

   \mathbf{b} = (b_x, b_y, 0)

   \mathbf{c} = (c_x, c_y, c_z)

It is always possible to put the box vectors into this form by rotating the
system until :math:`\mathbf{a}` is parallel to x and :math:`\mathbf{b}` lies in
the xy plane.

Furthermore, they must obey the following constraints:

.. math::
   a_x > 0, b_y > 0, c_z > 0

   a_x \ge 2 |b_x|

   a_x \ge 2 |c_x|

   b_y \ge 2 |c_y|

This effectively requires the box vectors to be specified in a particular
reduced form.  By forming combinations of box vectors (a process known as
"lattice reduction"), it is always possible to put them in this form without
changing the periodic system they represent.

These requirements have an important consequence: the periodic unit cell can
always be treated as an axis-aligned rectangular box of size
:math:`(a_x, b_y, c_z)`\ .  The remaining non-zero elements of the box vectors
cause the repeating copies of the system to be staggered relative to each other,
but they do not affect the shape or size of each copy.  The volume of the unit
cell is simply given by :math:`a_x b_y c_z`\ .

LocalEnergyMinimizer
********************

This provides an implementation of the L-BFGS optimization algorithm.
:cite:`Liu1989`  Given a Context specifying initial particle positions, it
searches for a nearby set of positions that represent a local minimum of the
potential energy.  Distance constraints are enforced during minimization by
adding a harmonic restraining force to the potential function.  The strength of
the restraining force is steadily increased until the minimum energy
configuration satisfies all constraints to within the tolerance specified by the
Context's Integrator.

XMLSerializer
*************

This provides the ability to “serialize” a System, Force, Integrator, or State
object to a portable XML format, then reconstruct it again later.  When
serializing a System, the XML data contains a complete copy of the entire system
definition, including all Forces that have been added to it.

Here are some examples of uses for this class:

#. A model building utility could generate a System in memory, then serialize it
   to a file on disk.  Other programs that perform simulation or analysis could
   then reconstruct the model by simply loading the XML file.
#. When running simulations on a cluster, all model construction could be done
   on a single node.  The Systems and Integrators could then be encoded as XML,
   allowing them to be easily transmitted to other nodes.


XMLSerializer is a templatized class that, in principle, can be used to
serialize any type of object.  At present, however, only System, Force,
Integrator, and State are supported.

Force Groups
************

It is possible to split the Force objects in a System into groups.  Those groups
can then be evaluated independently of each other.  This is done by calling
:code:`setForceGroup()` on the Force.  Some Force classes also
provide finer grained control over grouping.  For example, NonbondedForce allows
direct space computations to be in one group and reciprocal space computations
in a different group.

The most important use of force groups is for implementing multiple time step
algorithms with CustomIntegrator.  For example, you might evaluate the slowly
changing nonbonded interactions less frequently than the quickly changing bonded
ones.  This can be done by putting the slow and fast forces into separate
groups, then using a :class:`MTSIntegrator` or :class:`MTSLangevinIntegrator`
that evaluates the groups at different frequencies.

Another important use is to define forces that are not used when integrating
the equations of motion, but can still be queried efficiently.  To do this,
call :code:`setIntegrationForceGroups()` on the :class:`Integrator`.  Any groups
omitted will be ignored during simulation, but can be queried at any time by
calling :code:`getState()`.

Virtual Sites
*************

A virtual site is a particle whose position is computed directly from the
positions of other particles, not by integrating the equations of motion.  An
important example is the “extra sites” present in 4 and 5 site water models.
These particles are massless, and therefore cannot be integrated.  Instead,
their positions are computed from the positions of the massive particles in the
water molecule.

Virtual sites are specified by creating a VirtualSite object, then telling the
System to use it for a particular particle.  The VirtualSite defines the rules
for computing its position.  It is an abstract class with subclasses for
specific types of rules.  They are:

* TwoParticleAverageSite: The virtual site location is computed as a weighted
  average of the positions of two particles:

.. math::
   \mathbf{r}={w}_{1}\mathbf{r}_{1}+{w}_{2}\mathbf{r}_{2}

* ThreeParticleAverageSite: The virtual site location is computed as a weighted
  average of the positions of three particles:

.. math::
   \mathbf{r}={w}_{1}\mathbf{r}_{1}+{w}_{2}\mathbf{r}_{2}+{w}_{3}\mathbf{r}_{3}

* OutOfPlaneSite: The virtual site location is computed as a weighted average
  of the positions of three particles and the cross product of their relative
  displacements:

.. math::
   \mathbf{r}=\mathbf{r}_{1}+{w}_{12}\mathbf{r}_{12}+{w}_{13}\mathbf{r}_{13}+{w}_\mathit{cross}\left(\mathbf{r}_{12}\times \mathbf{r}_{13}\right)
..

  where :math:`\mathbf{r}_{12} = \mathbf{r}_{2}-\mathbf{r}_{1}` and
  :math:`\mathbf{r}_{13} = \mathbf{r}_{3}-\mathbf{r}_{1}`\ .  This allows
  the virtual site to be located outside the plane of the three particles.

* LocalCoordinatesSite: The locations of several other particles are used to compute a local
  coordinate system, and the virtual site is placed at a fixed location in that coordinate
  system.  The number of particles used to define the coordinate system is user defined.
  The origin of the coordinate system and the directions of its x and y axes
  are each specified as a weighted sum of the locations of the other particles:

.. math::
   \mathbf{o}={w}^{o}_{1}\mathbf{r}_{1} + {w}^{o}_{2}\mathbf{r}_{2} + ...

   \mathbf{dx}={w}^{x}_{1}\mathbf{r}_{1} + {w}^{x}_{2}\mathbf{r}_{2} + ...

   \mathbf{dy}={w}^{y}_{1}\mathbf{r}_{1} + {w}^{y}_{2}\mathbf{r}_{2} + ...

   \mathbf{dz}=\mathbf{dx}\times \mathbf{dy}
..

   These vectors are then used to construct a set of orthonormal coordinate axes as follows:

.. math::
   \mathbf{\hat{x}}=\mathbf{dx}/|\mathbf{dx}|

   \mathbf{\hat{z}}=\mathbf{dz}/|\mathbf{dz}|

   \mathbf{\hat{y}}=\mathbf{\hat{z}}\times \mathbf{\hat{x}}
..

   Finally, the position of the virtual site is set to

.. math::
   \mathbf{r}=\mathbf{o}+p_1\mathbf{\hat{x}}+p_2\mathbf{\hat{y}}+p_3\mathbf{\hat{z}}
..

Random Numbers with Stochastic Integrators and Forces
*****************************************************

OpenMM includes many stochastic integrators and forces that make extensive use
of random numbers. It is impossible to generate truly random numbers on a
computer like you would with a dice roll or coin flip in real life---instead
programs rely on pseudo-random number generators (PRNGs) that take some sort of
initial "seed" value and steps through a sequence of seemingly random numbers.

The exact implementation of the PRNGs is not important (in fact, each platform
may have its own PRNG whose performance is optimized for that hardware).  What
*is* important, however, is that the PRNG must generate a uniform distribution
of random numbers between 0 and 1. Random numbers drawn from this distribution
can be manipulated to yield random integers in a desired range or even a random
number from a different type of probability distribution function (e.g., a
normal distribution).

What this means is that the random numbers used by integrators and forces within
OpenMM cannot have any discernible pattern to them.  Patterns can be induced in
PRNGs in two principal ways:

1. The PRNG uses a bad algorithm with a short period.
2. Two PRNGs are started using the same seed

All PRNG algorithms in common use are periodic---that is their sequence of
random numbers repeats after a given *period*, defined by the number of "unique"
random numbers in the repeating pattern.  As long as this period is longer than
the total number of random numbers your application requires (preferably by
orders of magnitude), the first problem described above is avoided. All PRNGs
employed by OpenMM have periods far longer than any current simulation can cycle
through.

Point two is far more common in biomolecular simulation, and can result in very
strange artifacts that may be difficult to detect. For example, with Langevin
dynamics, two simulations that use the same sequence of random numbers appear to
synchronize in their global movements.\ :cite:`Uberuaga2004`\
:cite:`Sindhikara2009` It is therefore very important that the stochastic forces
and integrators in OpenMM generate unique sequences of pseudo-random numbers not
only within a single simulation, but between two different simulations of the
same system as well (including any restarts of previous simulations).

Every stochastic force and integrator that does (or could) make use of random
numbers has two instance methods attached to it: :meth:`getRandomNumberSeed()`
and :meth:`setRandomNumberSeed(int seed)`. If you set a unique random seed for
two different simulations (or different forces/integrators if applicable),
OpenMM guarantees that the generated sequences of random numbers will be
different (by contrast, no guarantee is made that the same seed will result in
identical random number sequences).

Since breaking simulations up into pieces and/or running multiple replicates of
a system to obtain more complete statistics is common practice, a new strategy
has been employed for OpenMM versions 6.3 and later with the aim of trying to
ensure that each simulation will be started with a unique random seed. A random
seed value of 0 (the default) will cause a unique random seed to be generated
when a new :class:`Context` is instantiated.

Prior to the introduction of this feature, deserializing a serialized
:class:`System` XML file would result in each stochastic force or integrator
being assigned the same random seed as the original instance that was
serialized. If you use a :class:`System` XML file generated by a version of
OpenMM older than 6.3 to start a new simulation, you should manually set the
random number seed of each stochastic force or integrator to 0 (or another
unique value).
.. _standard-forces:

Standard Forces
###############

The following classes implement standard force field terms that are widely used
in molecular simulations.

HarmonicBondForce
*****************

Each harmonic bond is represented by an energy term of the form



.. math::
   E=\frac{1}{2}k{\left(x-{x}_{0}\right)}^{2}


where *x* is the distance between the two particles, *x*\ :sub:`0` is
the equilibrium distance, and *k* is the force constant.  This produces a
force of magnitude *k*\ (\ *x*\ -\ *x*\ :sub:`0`\ ).

Be aware that some force fields define their harmonic bond parameters in a
slightly different way: *E* = *k*\ ´(\ *x*\ -\ *x*\ :sub:`0`\ )\
:sup:`2`\ , leading to a force of magnitude 2\ *k*\ ´(\ *x*\ -\ *x*\ :sub:`0`\ ).
Comparing these two forms, you can see that *k* = 2\ *k*\ ´.  Be sure to
check which form a particular force field uses, and if necessary multiply the
force constant by 2.

HarmonicAngleForce
******************

Each harmonic angle is represented by an energy term of the form


.. math::
   E=\frac{1}{2}k{\left(\theta-\theta_0\right)}^{2}


where :math:`\theta` is the angle formed by the three particles, :math:`\theta_0` is
the equilibrium angle, and *k* is the force constant.

As with HarmonicBondForce, be aware that some force fields define their harmonic
angle parameters as *E* = *k*\ ´(\ :math:`\theta`\ -\ :math:`\theta`\ :sub:`0`\ )\ :sup:`2`\ .
Be sure to check which form a particular force field uses, and if necessary
multiply the force constant by 2.

PeriodicTorsionForce
********************

Each torsion is represented by an energy term of the form


.. math::
   E=k\left(1+\text{cos}\left(n\theta-\theta_0\right)\right)


where :math:`\theta` is the dihedral angle formed by the four particles, :math:`\theta_0`
is the phase offset, *n* is the periodicity, and *k* is
the force constant.

RBTorsionForce
**************

Each torsion is represented by an energy term of the form


.. math::
   E=\sum _{i=0}^{5}{C}_{i}{\left(\text{cos}\phi\right)}^{i}


where :math:`\phi` is the dihedral angle formed by the four particles and
*C*\ :sub:`0` through *C*\ :sub:`5` are constant coefficients.

For reason of convention, PeriodicTorsionForce and RBTorsonForce define the
torsion angle differently. :math:`\theta` is zero when the first and last particles are
on the *same* side of the bond formed by the middle two particles (the
*cis* configuration), whereas :math:`\phi` is zero when they are on *opposite*
sides (the *trans* configuration).  This means that :math:`\theta` = :math:`\phi` - :math:`\pi`.

CMAPTorsionForce
****************

Each torsion pair is represented by an energy term of the form


.. math::
   E=f\left(\theta_1,\theta_2\right)


where :math:`\theta_1` and :math:`\theta_2` are the two dihedral angles
coupled by the term, and *f*\ (\ *x*\ ,\ *y*\ ) is defined by a user-supplied
grid of tabulated values.  A natural cubic spline surface is fit through the
tabulated values, then evaluated to determine the energy for arbitrary (\ :math:`\theta_1`\ ,
:math:`\theta_2`\ ) pairs.

NonbondedForce
**************

.. _lennard-jones-interaction:

Lennard-Jones Interaction
=========================

The Lennard-Jones interaction between each pair of particles is represented by
an energy term of the form


.. math::
   E=4\epsilon\left({\left(\frac{\sigma}{r}\right)}^{12}-{\left(\frac{\sigma}{r}\right)}^{6}\right)


where *r* is the distance between the two particles, :math:`\sigma` is the distance
at which the energy equals zero, and :math:`\epsilon` sets the strength of the
interaction.  If the NonbondedMethod in use is anything other than NoCutoff and
\ *r* is greater than the cutoff distance, the energy and force are both set
to zero.  Because the interaction decreases very quickly with distance, the
cutoff usually has little effect on the accuracy of simulations.

Optionally you can use a switching function to make the energy go smoothly to 0
at the cutoff distance.  When :math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by

.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}


where :math:`x = (r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`. This function decreases smoothly from 1 at
:math:`r = r_\mathit{switch}` to 0 at :math:`r = r_\mathit{cutoff}`, and has continuous first and
second derivatives at both ends.

When an exception has been added for a pair of particles, :math:`\sigma` and :math:`\epsilon`
are the parameters specified by the exception.  Otherwise they are determined
from the parameters of the individual particles using the Lorentz-Berthelot
combining rule:

.. math::
   \sigma=\frac{\sigma_1+\sigma_2}{2}

.. math::
   \epsilon=\sqrt{\epsilon_1 \epsilon_2}

When using periodic boundary conditions, NonbondedForce can optionally add a
term (known as a *long range dispersion correction*\ ) to the energy that
approximately represents the contribution from all interactions beyond the
cutoff distance:\ :cite:`Shirts2007`\

.. math::
   {E}_{\text{cor}}=\frac{{8\pi N}^{2}}{V}\left(\frac{\langle \epsilon_{ij}\sigma_{ij}^{12}\rangle}{{9r_c}^9}-\frac{\langle \epsilon_{ij}\sigma_{ij}^{6}\rangle}{{3r_c}^3}\right)

where *N* is the number of particles in the system, *V* is the volume of
the periodic box, :math:`r_c` is the cutoff distance, :math:`\sigma_{ij}` and
:math:`\epsilon_{ij}` are the interaction parameters between particle *i* and
particle *j*\ , and :math:`\langle \text{...} \rangle` represents an average over all pairs of particles in
the system.  When a switching function is in use, there is also a contribution
to the correction that depends on the integral of *E*\ ·(1-\ *S*\ ) over the
switching interval.  The long range dispersion correction is primarily useful
when running simulations at constant pressure, since it produces a more accurate
variation in system energy with respect to volume.

The Lennard-Jones interaction is often parameterized in two other equivalent
ways.  One is


.. math::
   E=\epsilon\left({\left(\frac{{r}_{\mathit{min}}}{r}\right)}^{\text{12}}-2{\left(\frac{{r}_{\mathit{min}}}{r}\right)}^{6}\right)


where :math:`r_\mathit{min}` (sometimes known as :math:`d_\mathit{min}`; this is not a
radius) is the center-to-center distance at which the energy is minimum.  It is
related to :math:`\sigma` by


.. math::
   \sigma=\frac{{r}_{\mathit{min}}}{{2}^{1/6}}


In turn, :math:`r_\mathit{min}` is related to the van der Waals radius by :math:`r_\mathit{min} = 2r_\mathit{vdw}`\ .

Another common form is



.. math::
   E=\frac{A}{{r}^{\text{12}}}-\frac{B}{{r}^{6}}


The coefficients A and B are related to :math:`\sigma` and :math:`\epsilon` by



.. math::
   \sigma={\left(\frac{A}{B}\right)}^{1/6}



.. math::
   \epsilon=\frac{{B}^{2}}{4A}


Coulomb Interaction Without Cutoff
==================================

The form of the Coulomb interaction between each pair of particles depends on
the NonbondedMethod in use.  For NoCutoff, it is given by


.. math::
   E=\frac{1}{4{\pi}{\epsilon}_{0}}\frac{{q}_{1}{q}_{2}}{r}


where *q*\ :sub:`1` and *q*\ :sub:`2` are the charges of the two
particles, and *r* is the distance between them.

Coulomb Interaction With Cutoff
===============================

For CutoffNonPeriodic or CutoffPeriodic, it is modified using the reaction field
approximation.  This is derived by assuming everything beyond the cutoff
distance is a solvent with a uniform dielectric constant.\ :cite:`Tironi1995`


.. math::
   E=\frac{{q}_{1}{q}_{2}}{4\pi\epsilon_0}\left(\frac{1}{r}+{k}_{\mathit{rf}}{r}^{2}-{c}_{\mathit{rf}}\right)


.. math::
   {k}_{\mathit{rf}}=\left(\frac{1}{{r_\mathit{cutoff}}^3}\right)\left(\frac{{\epsilon}_{\mathit{solvent}}-1}{2{\epsilon}_{\mathit{solvent}}+1}\right)


.. math::
   {c}_{\mathit{rf}}=\left(\frac{1}{{r}_{\mathit{cutoff}}}\right)\left(\frac{3{\epsilon}_{\mathit{solvent}}}{2{\epsilon}_{\mathit{solvent}}+1}\right)


where :math:`r_\mathit{cutoff}` is the cutoff distance and :math:`\epsilon_\mathit{solvent}` is
the dielectric constant of the solvent.  In the limit :math:`\epsilon_\mathit{solvent}` >> 1,
this causes the force to go to zero at the cutoff.

The reaction field approximation is not applied to nonbonded exceptions.  They
are always evaluated at full strength, regardless of the cutoff distance.  That
is because exceptions are primarily used to model 1-4 interactions, which are
really a type of bonded interaction.  They are typically parametrized without a
cutoff together with the other bonded interactions.

Coulomb Interaction With Ewald Summation
========================================

For Ewald, the total Coulomb energy is the sum of three terms: the *direct
space sum*\ , the *reciprocal space sum*\ , and the *self-energy term*\ .\
:cite:`Toukmaji1996`


.. math::
   E=E_{\mathit{dir}}+{E}_{\mathit{rec}}+{E}_{\mathit{self}}


.. math::
   E_{\mathit{dir}}=\frac{1}{2}\sum _{i,j}\sum_\mathbf{n}{q}_{i}{q}_{j}\frac{\text{erfc}\left({\mathit{\alpha r}}_{ij,\mathbf{n}}\right)}{r_{ij,\mathbf{n}}}


.. math::
   E_{\mathit{rec}}=\frac{1}{2{\pi}V}\sum _{i,j}q_i q_j\sum _{\mathbf{k}{\neq}0}\frac{\text{exp}(-(\pi \mathbf{k}/\alpha)^2+2\pi i \mathbf{k} \cdot (\mathbf{r}_{i}-\mathbf{r}_{j}))}{\mathbf{m}^2}


.. math::
   E_{\mathit{self}}=-\frac{\alpha}{\sqrt{\pi}}\sum _{i}{q}_{{i}^{2}}


In the above expressions, the indices *i* and *j* run over all
particles, **n** = (n\ :sub:`1`\ , n\ :sub:`2`\ , n\ :sub:`3`\ ) runs over
all copies of the periodic cell, and **k** = (k\ :sub:`1`\ , k\ :sub:`2`\ ,
k\ :sub:`3`\ ) runs over all integer wave vectors from (-k\ :sub:`max`\ ,
-k\ :sub:`max`\ , -k\ :sub:`max`\ ) to (k\ :sub:`max`\ , k\ :sub:`max`\ ,
k\ :sub:`max`\ ) excluding (0, 0, 0).  :math:`\mathbf{r}_i` is the position of
particle i , while :math:`r_{ij}` is the distance between particles *i* and *j*\ .
*V* is the volume of the periodic cell, and :math:`\alpha` is an internal parameter.

In the direct space sum, all pairs that are further apart than the cutoff
distance are ignored.  Because the cutoff is required to be less than half the
width of the periodic cell, the number of terms in this sum is never greater
than the square of the number of particles.

The error made by applying the direct space cutoff depends on the magnitude of
:math:`\text{erfc}({\alpha}r_\mathit{cutoff})`\ .  Similarly, the error made in the reciprocal space
sum by ignoring wave numbers beyond k\ :sub:`max` depends on the magnitude
of :math:`\text{exp}(-({\pi}k_{max}/{\alpha})^2`\ ).  By changing :math:`\alpha`, one can decrease the
error in either term while increasing the error in the other one.

Instead of having the user specify :math:`\alpha` and -k\ :sub:`max`\ , NonbondedForce
instead asks the user to choose an error tolerance :math:`\delta`.  It then calculates :math:`\alpha` as


.. math::
   \alpha =\sqrt{-\text{log}\left(2{\delta}\right)}/{r}_{\mathit{cutoff}}


Finally, it estimates the error in the reciprocal space sum as


.. math::
   \mathit{error}=\frac{k_{\mathit{max}}\sqrt{d\alpha}}{20}\text{exp}(-(\pi k_\mathit{max}/d\alpha)^2)


where *d* is the width of the periodic box, and selects the smallest value
for k\ :sub:`max` which gives *error* < :math:`\delta`\ .  (If the box is not square,
k\ :sub:`max` will have a different value along each axis.)

This means that the accuracy of the calculation is determined by :math:`\delta`\ .
:math:`r_\mathit{cutoff}` does not affect the accuracy of the result, but does affect the speed
of the calculation by changing the relative costs of the direct space and
reciprocal space sums.  You therefore should test different cutoffs to find the
value that gives best performance; this will in general vary both with the size
of the system and with the Platform being used for the calculation.  When the
optimal cutoff is used for every simulation, the overall cost of evaluating the
nonbonded forces scales as O(N\ :sup:`3/2`\ ) in the number of particles.

Be aware that the error tolerance :math:`\delta` is not a rigorous upper bound on the errors.
The formulas given above are empirically found to produce average relative
errors in the forces that are less than or similar to :math:`\delta` across a variety of
systems and parameter values, but no guarantees are made.  It is important to
validate your own simulations, and identify parameter values that produce
acceptable accuracy for each system.

Coulomb Interaction With Particle Mesh Ewald
============================================

The Particle Mesh Ewald (PME) algorithm\ :cite:`Essmann1995` is similar to
Ewald summation, but instead of calculating the reciprocal space sum directly,
it first distributes the particle charges onto nodes of a rectangular mesh using
5th order B-splines.  By using a Fast Fourier Transform, the sum can then be
computed very quickly, giving performance that scales as O(N log N) in the
number of particles (assuming the volume of the periodic box is proportional to
the number of particles).

As with Ewald summation, the user specifies the direct space cutoff :math:`r_\mathit{cutoff}`
and error tolerance :math:`\delta`\ .  NonbondedForce then selects :math:`\alpha` as


.. math::
   \alpha =\sqrt{-\text{log}\left(2\delta\right)}/{r}_\mathit{cutoff}


and the number of nodes in the mesh along each dimension as


.. math::
   n_\mathit{mesh}=\frac{2\alpha d}{{3\delta}^{1/5}}


where *d* is the width of the periodic box along that dimension.  Alternatively,
the user may choose to explicitly set values for these parameters.  (Note that
some Platforms may choose to use a larger value of :math:`n_\mathit{mesh}` than that
given by this equation.  For example, some FFT implementations require the mesh
size to be a multiple of certain small prime numbers, so a Platform might round
it up to the nearest permitted value.  It is guaranteed that :math:`n_\mathit{mesh}`
will never be smaller than the value given above.)

The comments in the previous section regarding the interpretation of :math:`\delta` for Ewald
summation also apply to PME, but even more so.  The behavior of the error for
PME is more complicated than for simple Ewald summation, and while the above
formulas will usually produce an average relative error in the forces less than
or similar to :math:`\delta`\ , this is not a rigorous guarantee.  PME is also more sensitive
to numerical round-off error than Ewald summation.  For Platforms that do
calculations in single precision, making :math:`\delta` too small (typically below about
5·10\ :sup:`-5`\ ) can actually cause the error to increase.

Lennard-Jones Interaction With Particle Mesh Ewald
==================================================

The PME algorithm can also be used for Lennard-Jones interactions.  Usually this
is not necessary, since Lennard-Jones forces are short ranged, but there are
situations (such as membrane simulations) where neglecting interactions beyond
the cutoff can measurably affect results.

For computational efficiency, certain approximations are made\ :cite:`Wennberg2015`.
Interactions beyond the cutoff distance include only the attractive :math:`1/r^6`
term, not the repulsive :math:`1/r^{12}` term.  Since the latter is much smaller
than the former at long distances, this usually has negligible effect.  Also,
the interaction between particles farther apart than the cutoff distance is
computed using geometric combination rules:

.. math::
   \sigma=\sqrt{\sigma_1 \sigma_2}

The effect of this approximation is also quite small, and it is still far more
accurate than ignoring the interactions altogether (which is what would happen
with PME).

The formula used to compute the number of nodes along each dimension of the mesh
is slightly different from the one used for Coulomb interactions:

.. math::
   n_\mathit{mesh}=\frac{\alpha d}{{3\delta}^{1/5}}

As before, this is an empirical formula.  It will usually produce an average
relative error in the forces less than or similar to :math:`\delta`\ , but that
is not guaranteed.

.. _gbsaobcforce:

GBSAOBCForce
************


Generalized Born Term
=====================

GBSAOBCForce consists of two energy terms: a Generalized Born Approximation term
to represent the electrostatic interaction between the solute and solvent, and a
surface area term to represent the free energy cost of solvating a neutral
molecule.  The Generalized Born energy is given by\ :cite:`Onufriev2004`


.. math::
   E\text{=-}\frac{1}{2}\left(\frac{1}{\epsilon_{\mathit{solute}}}-\frac{1}{\epsilon_{\mathit{solvent}}}\right)\sum _{i,j}\frac{{q}_{i}{q}_{j}}{{f}_{\text{GB}}\left({d}_{ij},{R}_{i},{R}_{j}\right)}


where the indices *i* and *j* run over all particles, :math:`\epsilon_\mathit{solute}`
and :math:`\epsilon_\mathit{solvent}` are the dielectric constants of the solute and solvent
respectively, :math:`q_i` is the charge of particle *i*\ , and :math:`d_{ij}` is the distance
between particles *i* and *j*\ .  :math:`f_\text{GB}(d_{ij}, R_i, R_j)` is defined as


.. math::
   {f}_{\text{GB}}\left({d}_{ij},{R}_{i},{R}_{j}\right)={\left[{d}_{ij}^2+{R}_{i}{R}_{j}\text{exp}\left(\frac{-{d}_{ij}^2}{{4R}_{i}{R}_{j}}\right)\right]}^{1/2}


:math:`R_i` is the Born radius of particle *i*\ , which calculated as


.. math::
   R_i=\frac{1}{\rho_i^{-1}-r_i^{-1}\text{tanh}\left(\alpha \Psi_{i}-{\beta \Psi}_i^2+{\gamma \Psi}_i^3\right)}


where :math:`\alpha`, :math:`\beta`, and :math:`\gamma` are the GB\ :sup:`OBC`\ II parameters :math:`\alpha` = 1, :math:`\beta` = 0.8, and :math:`\gamma` =
4.85.  :math:`\rho_i` is the adjusted atomic radius of particle *i*\ , which
is calculated from the atomic radius :math:`r_i` as :math:`\rho_i = r_i-0.009` nm.
:math:`\Psi_i` is calculated as an integral over the van der Waals
spheres of all particles outside particle *i*\ :


.. math::
   \Psi_i=\frac{\rho_i}{4\pi}\int_{\text{VDW}}\theta\left(|\mathbf{r}|-{\rho }_{i}\right)\frac{1}{{|\mathbf{r}|}^{4}}{d}^{3}\mathbf{r}


where :math:`\theta`\ (\ *r*\ ) is a step function that excludes the interior of particle
\ *i* from the integral.

Surface Area Term
=================

The surface area term is given by\ :cite:`Schaefer1998`\ :cite:`Ponder`


.. math::
   E=E_{SA} \cdot 4\pi \sum _{i}{\left({r}_{i}+{r}_{\mathit{solvent}}\right)}^{2}{\left(\frac{{r}_{i}}{{R}_{i}}\right)}^{6}


where :math:`r_i` is the atomic radius of particle *i*\ , :math:`r_i` is
its atomic radius, and :math:`r_\mathit{solvent}` is the solvent radius, which is taken
to be 0.14 nm.  The default value for the energy scale :math:`E_{SA}` is 2.25936 kJ/mol/nm\ :sup:`2`\ .


GayBerneForce
*************

This is similar to the Lennard-Jones interaction described in section :numref:`lennard-jones-interaction`,
but instead of being based on the distance between two point particles, it is based
on the distance of closest approach between two ellipsoids.\ :cite:`Everaers2003`
Let :math:`\mathbf{A}_1` and :math:`\mathbf{A}_2` be rotation matrices that transform
from the lab frame to the body frames of two interacting ellipsoids.  These rotations
are determined from the positions of other particles, as described in the API documentation.
Let :math:`\mathbf{r}_{12}` be the vector pointing from particle 1 to particle 2, and
:math:`\hat{\mathbf{r}}_{12}=\mathbf{r}_{12}/|\mathbf{r}_{12}|`.  Let :math:`\mathbf{S}_1`
and :math:`\mathbf{S}_2` be diagonal matrices containing the three radii of each particle:

.. math::
   \mathbf{S}_i=\begin{bmatrix}
   a_i & 0 & 0 \\
   0 & b_i & 0 \\
   0 & 0 & c_i
   \end{bmatrix}

The energy is computed as a product of three terms:

.. math::
   E=U_r(\mathbf{A}_1, \mathbf{A}_2, \mathbf{r}_{12}) \cdot \eta_{12}(\mathbf{A}_1, \mathbf{A}_2) \cdot \chi_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})

The first term describes the distance dependence, and is very similar in form to
the Lennard-Jones interaction:

.. math::
   U_r=4\epsilon\left({\left(\frac{\sigma}{h_{12}+\sigma}\right)}^{12}-{\left(\frac{\sigma}{h_{12}+\sigma}\right)}^{6}\right)

where :math:`h_{12}` is an approximation to the distance of closest approach between
the two ellipsoids:

.. math::
   h_{12}=|\mathbf{r}_{12}|-\sigma_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})

.. math::
   \sigma_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})=\left[ \frac{1}{2} \hat{\mathbf{r}}_{12}^T \mathbf{G}_{12}^{-1} \hat{\mathbf{r}}_{12} \right]^{-1/2}

.. math::
   \mathbf{G}_{12}=\mathbf{A}_1^T \mathbf{S}_1^2 \mathbf{A}_1 + \mathbf{A}_2^T \mathbf{S}_2^2 \mathbf{A}_2

The second term adjusts the energy based on the relative orientations of the two ellipsoids:

.. math::
   \eta_{12}(\mathbf{A}_1, \mathbf{A}_2)=\left[ \frac{2 s_1 s_2}{\text{det}(\mathbf{G}_{12})} \right]^{1/2}

.. math::
   s_i=(a_i b_i + c_i^2)\sqrt{a_i b_i}

The third term applies the user-defined scale factors :math:`e_a`, :math:`e_b`,
and :math:`e_c` that adjust the strength of the interaction along each axis:

.. math::
   \chi_{12}(\mathbf{A}_1, \mathbf{A}_2, \hat{\mathbf{r}}_{12})=(2 \hat{\mathbf{r}}_{12}^T \mathbf{B}_{12}^{-1} \hat{\mathbf{r}}_{12})^2

.. math::
   \mathbf{B}_{12}=\mathbf{A}_1^T \mathbf{E}_1 \mathbf{A}_1 + \mathbf{A}_2^T \mathbf{E}_2 \mathbf{A}_2

.. math::
   \mathbf{E}_i=\begin{bmatrix}
   e_{ai}^{-1/2} & 0 & 0 \\
   0 & e_{bi}^{-1/2} & 0 \\
   0 & 0 & e_{ci}^{-1/2}
   \end{bmatrix}

When using a cutoff, you can optionally use a switching function to make the energy go smoothly to 0
at the cutoff distance.  When :math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by

.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}

where :math:`x = (r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`. This function decreases smoothly from 1 at
:math:`r = r_\mathit{switch}` to 0 at :math:`r = r_\mathit{cutoff}`, and has continuous first and
second derivatives at both ends.


AndersenThermostat
******************

AndersenThermostat couples the system to a heat bath by randomly selecting a
subset of particles at the start of each time step, then setting their
velocities to new values chosen from a Boltzmann distribution.  This represents
the effect of random collisions between particles in the system and particles in
the heat bath.\ :cite:`Andersen1980`

The probability that a given particle will experience a collision in a given
time step is


.. math::
   P=1-{e}^{-f\Delta t}


where *f* is the collision frequency and :math:`\Delta t` is the step size.
Each component of its velocity is then set to


.. math::
   {v}_{i}=\sqrt{\frac{{k}_{B}T}{m}}R


where *T* is the thermostat temperature, *m* is the particle mass, and
*R* is a random number chosen from a normal distribution with mean of zero and
variance of one.

MonteCarloBarostat
******************

MonteCarloBarostat models the effect of constant pressure by allowing the size
of the periodic box to vary with time.\ :cite:`Chow1995`\ :cite:`Aqvist2004`
At regular intervals, it attempts a Monte Carlo step by scaling the box vectors
and the coordinates of each molecule’s center by a factor *s*\ .  The scale
factor *s* is chosen to change the volume of the periodic box from *V*
to *V*\ +\ :math:`\Delta`\ *V*\ :


.. math::
   s={\left(\frac{V+\Delta V}{V}\right)}^{1/3}


The change in volume is chosen randomly as


.. math::
   \Delta V=A\cdot r


where *A* is a scale factor and *r* is a random number uniformly
distributed between -1 and 1.  The step is accepted or rejected based on the
weight function


.. math::
   \Delta W=\Delta E+P\Delta V-Nk_{B}T \text{ln}\left(\frac{V+\Delta V}{V}\right)


where :math:`\Delta E` is the change in potential energy resulting from the step,
\ *P* is the pressure being applied to the system, *N* is the number of molecules in the
system, :math:`k_B` is Boltzmann’s constant, and *T* is the system
temperature.  In particular, if :math:`\Delta W\le 0` the step is always accepted.
If :math:`\Delta W > 0`\ , the step is accepted with probability
:math:`\text{exp}(-\Delta W/k_B T)`\ .

This algorithm tends to be more efficient than deterministic barostats such as
the Berendsen or Parrinello-Rahman algorithms, since it does not require an
expensive virial calculation at every time step.  Each Monte Carlo step involves
two energy evaluations, but this can be done much less often than every time
step.  It also does not require you to specify the compressibility of the
system, which usually is not known in advance.

The scale factor *A* that determines the size of the steps is chosen
automatically to produce an acceptance rate of approximately 50%.  It is
initially set to 1% of the periodic box volume.  The acceptance rate is then
monitored, and if it varies too much from 50% then *A* is modified
accordingly.

Each Monte Carlo step modifies particle positions by scaling the centroid of
each molecule, then applying the resulting displacement to each particle in the
molecule.  This ensures that each molecule is translated as a unit, so bond
lengths and constrained distances are unaffected.

MonteCarloBarostat assumes the simulation is being run at constant temperature
as well as pressure, and the simulation temperature affects the step acceptance
probability.  It does not itself perform temperature regulation, however.  You
must use another mechanism along with it to maintain the temperature, such as
LangevinIntegrator or AndersenThermostat.

MonteCarloAnisotropicBarostat
*****************************

MonteCarloAnisotropicBarostat is very similar to MonteCarloBarostat, but instead
of scaling the entire periodic box uniformly, each Monte Carlo step scales only
one axis of the box.  This allows the box to change shape, and is useful for
simulating anisotropic systems whose compressibility is different along
different directions.  It also allows a different pressure to be specified for
each axis.

You can specify that the barostat should only be applied to certain axes of the
box, keeping the other axes fixed.  This is useful, for example, when doing
constant surface area simulations of membranes.

MonteCarloMembraneBarostat
**************************

MonteCarloMembraneBarostat is very similar to MonteCarloBarostat, but it is
specialized for simulations of membranes.  It assumes the membrane lies in the
XY plane.  In addition to applying a uniform pressure to regulate the volume of
the periodic box, it also applies a uniform surface tension to regulate the
cross sectional area of the periodic box in the XY plane.  The weight function
for deciding whether to accept a step is

.. math::
   \Delta W=\Delta E+P\Delta V-S\Delta A-Nk_{B}T \text{ln}\left(\frac{V+\Delta V}{V}\right)

where *S* is the surface tension and :math:`\Delta`\ *A* is the change in cross
sectional area.  Notice that pressure and surface tension are defined with
opposite senses: a larger pressure tends to make the box smaller, but a larger
surface tension tends to make the box larger.

MonteCarloMembraneBarostat offers some additional options to customize the
behavior of the periodic box:

* The X and Y axes can be either

  * isotropic (they are always scaled by the same amount, so their ratio remains fixed)
  * anisotropic (they can change size independently)

* The Z axis can be either

  * free (its size changes independently of the X and Y axes)
  * fixed (its size does not change)
  * inversely varying with the X and Y axes (so the total box volume does not
    change)

MonteCarloFlexibleBarostat
**************************

MonteCarloFlexibleBarostat is very similar to MonteCarloBarostat, but it allows
the periodic box to be fully flexible.\ :cite:`Vandenhaute2021`  Monte Carlo
steps can change not just the lengths of the box sides, but also the angles.  It
is especially useful for simulations of bulk materials where the shape of a
crystal's unit cell may not be known in advance, or could even change with time
as it transitions between phases.

CMMotionRemover
***************

CMMotionRemover prevents the system from drifting in space by periodically
removing all center of mass motion.  At the start of every *n*\ ’th time step
(where *n* is set by the user), it calculates the total center of mass
velocity of the system:


.. math::
   \mathbf{v}_\text{CM}=\frac{\sum _{i}{m}_{i}\mathbf{v}_{i}}{\sum _{i}{m}_{i}}


where :math:`m_i` and :math:`\mathbf{v}_i` are the mass and velocity of particle
\ *i*\ .  It then subtracts :math:`\mathbf{v}_\text{CM}` from the velocity of every
particle.

RMSDForce
*********

RMSDForce computes the root-mean-squared deviation (RMSD) between the current
particle positions :math:`\mathbf{x}_i` and a set of reference positions
:math:`\mathbf{x}_i^\text{ref}`:

.. math::
   \text{RMSD} = \sqrt{\frac{\sum_{i} \| \mathbf{x}_i - \mathbf{x}_i^\text{ref} \|^2}{N}}

Before computing this, the reference positions are first translated and rotated
so as to minimize the RMSD.  The computed value is therefore :math:`argmin(\text{RMSD})`,
where the :math:`argmin` is taken over all possible translations and rotations.

This force is normally used with a CustomCVForce (see Section :numref:`customcvforce`).
One rarely wants a force whose energy exactly equals the RMSD, but there are many
situations where it is useful to have a restraining or biasing force that depends
on the RMSD in some way.
.. _custom-forces:

Custom Forces
#############

In addition to the standard forces described in the previous chapter, OpenMM
provides a number of “custom” force classes.   These classes provide detailed
control over the mathematical form of the force by allowing the user to specify
one or more arbitrary algebraic expressions.  The details of how to write these
custom expressions are described in section :numref:`writing-custom-expressions`\ .

CustomBondForce
***************

CustomBondForce is similar to HarmonicBondForce in that it represents an
interaction between certain pairs of particles as a function of the distance
between them, but it allows the precise form of the interaction to be specified
by the user.  That is, the interaction energy of each bond is given by


.. math::
   E=f\left(r\right)


where *f*\ (\ *r*\ ) is a user defined mathematical expression.

In addition to depending on the inter-particle distance *r*\ , the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in two ways:

* Global parameters have a single, fixed value.
* Per-bond parameters are defined by specifying a value for each bond.


CustomAngleForce
****************

CustomAngleForce is similar to HarmonicAngleForce in that it represents an
interaction between sets of three particles as a function of the angle between
them, but it allows the precise form of the interaction to be specified by the
user.  That is, the interaction energy of each angle is given by


.. math::
   E=f\left(\theta\right)


where :math:`f(\theta)` is a user defined mathematical expression.

In addition to depending on the angle :math:`\theta`\ , the energy may also depend on an
arbitrary set of user defined parameters.  Parameters may be specified in two
ways:

* Global parameters have a single, fixed value.
* Per-angle parameters are defined by specifying a value for each angle.


CustomTorsionForce
******************

CustomTorsionForce is similar to PeriodicTorsionForce in that it represents an
interaction between sets of four particles as a function of the dihedral angle
between them, but it allows the precise form of the interaction to be specified
by the user.  That is, the interaction energy of each angle is given by


.. math::
   E=f(\theta)


where :math:`f(\theta)` is a user defined mathematical expression.  The angle
:math:`\theta` is guaranteed to be in the range :math:`[-\pi, +\pi]`\ .  Like PeriodicTorsionForce, it
is defined to be zero when the first and last particles are on the same side of
the bond formed by the middle two particles (the *cis* configuration).

In addition to depending on the angle :math:`\theta`\ , the energy may also depend on an
arbitrary set of user defined parameters.  Parameters may be specified in two
ways:

* Global parameters have a single, fixed value.
* Per-torsion parameters are defined by specifying a value for each torsion.


.. _customnonbondedforce:

CustomNonbondedForce
********************

CustomNonbondedForce is similar to NonbondedForce in that it represents a
pairwise interaction between all particles in the System, but it allows the
precise form of the interaction to be specified by the user.  That is, the
interaction energy between each pair of particles is given by


.. math::
   E=f(r)


where *f*\ (\ *r*\ ) is a user defined mathematical expression.

In addition to depending on the inter-particle distance *r*\ , the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.


A CustomNonbondedForce can optionally be restricted to only a subset of particle
pairs in the System.  This is done by defining “interaction groups”.  See the
API documentation for details.

When using a cutoff, a switching function can optionally be applied to make the
energy go smoothly to 0 at the cutoff distance.  When
:math:`r_\mathit{switch} < r < r_\mathit{cutoff}`\ , the energy is multiplied by



.. math::
   S=1-{6x}^{5}+15{x}^{4}-10{x}^{3}


where :math:`x=(r-r_\mathit{switch})/(r_\mathit{cutoff}-r_\mathit{switch})`\ .
This function decreases smoothly from 1 at :math:`r=r_\mathit{switch}`
to 0 at :math:`r=r_\mathit{cutoff}`\ , and has continuous first and
second derivatives at both ends.

When using periodic boundary conditions, CustomNonbondedForce can optionally add
a term (known as a *long range truncation correction*\ ) to the energy that
approximately represents the contribution from all interactions beyond the
cutoff distance:\ :cite:`Shirts2007`


.. math::
   {E}_{cor}=\frac{2\pi N^2}{V}\left\langle\underset{{r}_\mathit{cutoff}}{\overset{\infty}{\int}}E(r)r^{2}dr\right\rangle


where *N* is the number of particles in the system, *V* is the volume of
the periodic box, and :math:`\langle \text{...} \rangle` represents an average over all pairs of particles in
the system.  When a switching function is in use, there is an additional
contribution to the correction given by


.. math::
   E_{cor}^\prime=\frac{2\pi N^2}{V}\left\langle\underset{{r}_\mathit{switch}}{\overset{{r}_\mathit{cutoff}}{\int }}E(r)(1-S(r))r^{2}dr\right\rangle


The long range dispersion correction is primarily useful when running
simulations at constant pressure, since it produces a more accurate variation in
system energy with respect to volume.

CustomExternalForce
*******************

CustomExternalForce represents a force that is applied independently to each
particle as a function of its position.   That is, the energy of each particle
is given by


.. math::
   E=f(x,y,z)


where *f*\ (\ *x*\ , *y*\ , *z*\ ) is a user defined mathematical
expression.

In addition to depending on the particle’s (\ *x*\ , *y*\ , *z*\ )
coordinates, the energy may also depend on an arbitrary set of user defined
parameters.  Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.


CustomCompoundBondForce
***********************

CustomCompoundBondForce supports a wide variety of bonded interactions.  It
defines a “bond” as a single energy term that depends on the positions of a
fixed set of particles.  The number of particles involved in a bond, and how the
energy depends on their positions, is configurable.  It may depend on the
positions of individual particles, the distances between pairs of particles, the
angles formed by sets of three particles, and the dihedral angles formed by sets
of four particles.  That is, the interaction energy of each bond is given by


.. math::
   E=f(\{x_i\},\{r_i\},\{\theta_i\},\{\phi_i\})


where *f*\ (\ *...*\ ) is a user defined mathematical expression.  It may
depend on an arbitrary set of positions {\ :math:`x_i`\ }, distances {\ :math:`r_i`\ },
angles {\ :math:`\theta_i`\ }, and dihedral angles {\ :math:`\phi_i`\ }
guaranteed to be in the range :math:`[-\pi, +\pi]`\ .

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from among the particles that make up the bond.  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on positions, distances, angles, and dihedrals, the
energy may also depend on an arbitrary set of user defined parameters.
Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-bond parameters are defined by specifying a value for each bond.


CustomCentroidBondForce
***********************

CustomCentroidBondForce is very similar to CustomCompoundBondForce, but instead
of creating bonds between individual particles, the bonds are between the
centers of groups of particles.  This is useful for purposes such as restraining
the distance between two molecules or pinning the center of mass of a single
molecule.

The first step in computing this force is to calculate the center position of
each defined group of particles.  This is calculated as a weighted average of
the positions of all the particles in the group, with the weights being user
defined.  The computation then proceeds exactly as with CustomCompoundBondForce,
but the energy of each "bond" is now calculated based on the centers of a set
of groups, rather than on the positions of individual particles.

This class supports all the same function types and features as
CustomCompoundBondForce.  In fact, any interaction that could be implemented
with CustomCompoundBondForce can also be implemented with this class, simply by
defining each group to contain only a single atom.


CustomManyParticleForce
***********************

CustomManyParticleForce is similar to CustomNonbondedForce in that it represents
a custom nonbonded interaction between particles, but it allows the interaction
to depend on more than two particles.  This allows it to represent a wide range
of non-pairwise interactions.  It is defined by specifying the number of
particles :math:`N` involved in the interaction and how the energy depends on
their positions.  More specifically, it takes a user specified energy function

.. math::
   E=f(\{x_i\},\{r_i\},\{\theta_i\},\{\phi_i\})

that may depend on an arbitrary set of positions {\ :math:`x_i`\ }, distances
{\ :math:`r_i`\ }, angles {\ :math:`\theta_i`\ }, and dihedral angles
{\ :math:`\phi_i`\ } from a particular set of :math:`N` particles.

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from among the particles in the set.  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on positions, distances, angles, and dihedrals, the
energy may also depend on an arbitrary set of user defined parameters.
Parameters may be specified in two ways:

* Global parameters have a single, fixed value.
* Per-particle parameters are defined by specifying a value for each particle.

The energy function is evaluated one or more times for every unique set of
:math:`N` particles in the system.  The exact number of times depends on the
*permutation mode*\ .  A set of :math:`N` particles has :math:`N!` possible
permutations.  In :code:`SinglePermutation` mode, the function is evaluated
for a single arbitrarily chosen one of those permutations.  In
:code:`UniqueCentralParticle` mode, the function is evaluated for :math:`N` of
those permutations, once with each particle as the "central particle".

The number of times the energy function is evaluated can be further restricted
by specifying *type filters*\ .  Each particle may have a "type" assigned to it,
and then each of the :math:`N` particles involved in an interaction may be
restricted to only a specified set of types.  This provides a great deal of
flexibility in controlling which particles interact with each other.


CustomGBForce
*************

CustomGBForce implements complex, multiple stage nonbonded interactions between
particles.  It is designed primarily for implementing Generalized Born implicit
solvation models, although it is not strictly limited to that purpose.

The interaction is specified as a series of computations, each defined by an
arbitrary algebraic expression.  These computations consist of some number of
per-particle *computed values*\ , followed by one or more *energy terms*\ .
A computed value is a scalar value that is computed for each particle in the
system.  It may depend on an arbitrary set of global and per-particle
parameters, and well as on other computed values that have been calculated
before it.  Once all computed values have been calculated, the energy terms and
their derivatives are evaluated to determine the system energy and particle
forces.  The energy terms may depend on global parameters, per-particle
parameters, and per-particle computed values.

Computed values can be calculated in two different ways:

* *Single particle* values are calculated by evaluating a user defined
  expression for each particle:

.. math::
  {value}_{i}=f\left(\text{.}\text{.}\text{.}\right)
..

  where *f*\ (...) may depend only on properties of particle *i* (its
  coordinates and parameters, as well as other computed values that have already
  been calculated).

* *Particle pair* values are calculated as a sum over pairs of particles:

.. math::
  {value}_{i}=\sum _{j\ne i}f\left(r,\text{...}\right)
..

  where the sum is over all other particles in the System, and *f*\ (\ *r*\ ,
  ...) is a function of the distance *r* between particles *i* and *j*\,
  as well as their parameters and computed values.

Energy terms may similarly be calculated per-particle or per-particle-pair.

* *Single particle* energy terms are calculated by evaluating a user
  defined expression for each particle:

.. math::
  E=f\left(\text{.}\text{.}\text{.}\right)
..

  where *f*\ (...) may depend only on properties of that particle (its
  coordinates, parameters, and computed values).

* *Particle pair* energy terms are calculated by evaluating a user defined
  expression once for every pair of particles in the System:

.. math::
  E=\sum _{i,j}f\left(r,\text{.}\text{.}\text{.}\right)
..

  where the sum is over all particle pairs *i* *< j*\ , and *f*\ (\ *r*\ ,
  ...) is a function of the distance *r* between particles *i* and *j*\,
  as well as their parameters and computed values.

Note that energy terms are assumed to be symmetric with respect to the two
interacting particles, and therefore are evaluated only once per pair.  In
contrast, expressions for computed values need not be symmetric and therefore
are calculated twice for each pair: once when calculating the value for the
first particle, and again when calculating the value for the second particle.

Be aware that, although this class is extremely general in the computations it
can define, particular Platforms may only support more restricted types of
computations.  In particular, all currently existing Platforms require that the
first computed value *must* be a particle pair computation, and all computed
values after the first *must* be single particle computations. This is
sufficient for most Generalized Born models, but might not permit some other
types of calculations to be implemented.

CustomHbondForce
****************

CustomHbondForce supports a wide variety of energy functions used to represent
hydrogen bonding.  It computes interactions between "donor" particle groups and
"acceptor" particle groups, where each group may include up to three particles.
Typically a donor group consists of a hydrogen atom and the atoms it is bonded
to, and an acceptor group consists of a negatively charged atom and the atoms it
is bonded to.  The interaction energy between each donor group and each acceptor
group is given by


.. math::
   E=f(\{r_i\},\{\theta_i\},\{\phi_i\})


where *f*\ (\ *...*\ ) is a user defined mathematical expression.  It may
depend on an arbitrary set of distances {\ :math:`r_i`\ }, angles {\ :math:`\theta_i`\ },
and dihedral angles {\ :math:`\phi_i`\ }.

Each distance, angle, or dihedral is defined by specifying a sequence of
particles chosen from the interacting donor and acceptor groups (up to six atoms
to choose from, since each group may contain up to three atoms).  A distance
variable is defined by two particles, and equals the distance between them.  An
angle variable is defined by three particles, and equals the angle between them.
A dihedral variable is defined by four particles, and equals the angle between
the first and last particles about the axis formed by the middle two particles.
It is equal to zero when the first and last particles are on the same side of
the axis.

In addition to depending on distances, angles, and dihedrals, the energy may
also depend on an arbitrary set of user defined parameters.  Parameters may be
specified in three ways:

* Global parameters have a single, fixed value.
* Per-donor parameters are defined by specifying a value for each donor group.
* Per-acceptor parameters are defined by specifying a value for each acceptor group.

.. _customcvforce:

CustomCVForce
*************

CustomCVForce computes an energy as a function of "collective variables".  A
collective variable may be any scalar valued function of the particle positions
and other parameters.  Each one is defined by a :code:`Force` object, so any
function that can be defined via any force class (either standard or custom) can
be used as a collective variable.  The energy is then computed as

.. math::
   E=f(...)

where *f*\ (...) is a user supplied mathematical expression of the collective
variables.  It also may depend on user defined global parameters.


.. _writing-custom-expressions:

Writing Custom Expressions
**************************

The custom forces described in this chapter involve user defined algebraic
expressions.  These expressions are specified as character strings, and may
involve a variety of standard operators and mathematical functions.

The following operators are supported: + (add), - (subtract), * (multiply), /
(divide), and ^ (power).  Parentheses “(“ and “)” may be used for grouping.

The following standard functions are supported: sqrt, exp, log, sin, cos, sec,
csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs,
floor, ceil, step, delta, select. step(x) = 0 if x < 0, 1 otherwise.
delta(x) = 1 if x is 0, 0 otherwise.  select(x,y,z) = z if x = 0, y otherwise.
Some custom forces allow additional functions to be defined from tabulated values.

Numbers may be given in either decimal or exponential form.  All of the
following are valid numbers: 5, -3.1, 1e6, and 3.12e-2.

The variables that may appear in expressions are specified in the API
documentation for each force class.  In addition, an expression may be followed
by definitions for intermediate values that appear in the expression.  A
semicolon “;” is used as a delimiter between value definitions.  For example,
the expression
::

    a^2+a*b+b^2; a=a1+a2; b=b1+b2

is exactly equivalent to
::

    (a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2

The definition of an intermediate value may itself involve other intermediate
values.  All uses of a value must appear *before* that value’s definition.

Setting Parameters
******************

Most custom forces have two types of parameters you can define.  The simplest type
are global parameters, which represent a single number.  The value is stored in
the :class:`Context`, and can be changed at any time by calling :meth:`setParameter`
on it.  Global parameters are designed to be very inexpensive to change.  Even if
you set a new value for a global parameter on every time step, the overhead will
usually be quite small.  There can be exceptions to this rule, however.  For
example, if a :class:`CustomNonbondedForce` uses a long range correction, changing
a global parameter may require the correction coefficient to be recalculated,
which is expensive.

The other type of parameter is ones that record many values, one for each element
of the force, such as per-particle or per-bond parameters.  These values are stored
directly in the force object itself, and hence are part of the system definition.
When a :class:`Context` is created, the values are copied over to it, and thereafter
the two are disconnected.  Modifying the force will have no effect on any
:class:`Context` that already exists.

Some forces do provide a way to modify these parameters via an :meth:`updateParametersInContext`
method.  These methods tend to be somewhat expensive, so it is best not to call
them too often.  On the other hand, they are still much less expensive than calling
:meth:`reinitialize` on the :class:`Context`, which is the other way of updating
the system definition for a running simulation.

Parameter Derivatives
*********************

Many custom forces have the ability to compute derivatives of the potential energy
with respect to global parameters.  To use this feature, first define a global
parameter that the energy depends on.  Then instruct the custom force to compute
the derivative with respect to that parameter by calling :meth:`addEnergyParameterDerivative()`
on it.  Whenever forces and energies are computed, the specified derivative will
then also be computed at the same time.  You can query it by calling :meth:`getState()`
on a :class:`Context`, just as you would query forces or energies.

An important application of this feature is to use it in combination with a
:class:`CustomIntegrator` (described in section :numref:`custom-integrator`\ ).  The
derivative can appear directly in expressions that define the integration
algorithm.  This can be used to implement algorithms such as lambda-dynamics,
where a global parameter is integrated as a dynamic variable.

.. _integrators-theory:

Integrators
###########


VerletIntegrator
****************

VerletIntegrator implements the leap-frog Verlet integration method.  The
positions and velocities stored in the context are offset from each other by
half a time step.  In each step, they are updated as follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2)=\mathbf{v}_{i}(t-\Delta t/2)+\mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\mathbf{v}_{i}(t+\Delta t/2)\Delta t


where :math:`\mathbf{v}_i` is the velocity of particle *i*\ , :math:`\mathbf{r}_i` is
its position, :math:`\mathbf{f}_i` is the force acting on it, :math:`m_i` is its
mass, and :math:`\Delta t` is the time step.

Because the positions are always half a time step later than the velocities,
care must be used when calculating the energy of the system.  In particular, the
potential energy and kinetic energy in a State correspond to different times,
and you cannot simply add them to get the total energy of the system.  Instead,
it is better to retrieve States after two successive time steps, calculate the
on-step velocities as


.. math::
   \mathbf{v}_{i}(t)=\frac{\mathbf{v}_{i}\left(t-\Delta t/2\right)+\mathbf{v}_{i}\left(t+\Delta t/2\right)}{2}


then use those velocities to calculate the kinetic energy at time *t*\ .

LangevinIntegator
*****************

LangevinIntegator simulates a system in contact with a heat bath by integrating
the Langevin equation of motion:


.. math::
   m_i\frac{d\mathbf{v}_i}{dt}=\mathbf{f}_i-\gamma m_i \mathbf{v}_i+\mathbf{R}_i


where :math:`\mathbf{v}_i` is the velocity of particle *i*\ , :math:`\mathbf{f}_i` is
the force acting on it, :math:`m_i` is its mass, :math:`\gamma` is the friction
coefficient, and :math:`\mathbf{R}_i` is an uncorrelated random force whose
components are chosen from a normal distribution with mean zero and variance
:math:`2m_i \gamma k_B T`\ , where *T* is the temperature of
the heat bath.

The integration is done using the Langevin leap-frog method. :cite:`Izaguirre2010`
In each step, the positions and velocities are updated as follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2)=\mathbf{v}_{i}(t-\Delta t/2)\alpha+\mathbf{f}_{i}(t)(1-\alpha)/\gamma{m}_{i} + \sqrt{kT(1-\alpha^2)/m}R


.. math::
   \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\mathbf{v}_{i}(t+\Delta t/2)\Delta t


where :math:`k` is Boltzmann's constant, :math:`T` is the temperature,
:math:`\gamma` is the friction coefficient, :math:`R` is a normally distributed
random number, and :math:`\alpha=\exp(-\gamma\Delta t)`.

The same comments about the offset between positions and velocities apply to
this integrator as to VerletIntegrator.

LangevinMiddleIntegrator
************************

This integrator is similar to LangevinIntegerator, but it instead uses the LFMiddle
discretization. :cite:`Zhang2019` In each step, the positions and velocities
are updated as follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t-\Delta t/2) + \mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t/2) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t+\Delta t/2)\Delta t/2


.. math::
   \mathbf{v'}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t+\Delta t/2)\alpha + \sqrt{kT(1-\alpha^2)/m}R


.. math::
   \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t+\Delta t/2) + \mathbf{v'}_{i}(t+\Delta t/2)\Delta t/2


This tends to produce more accurate sampling of configurational properties (such
as free energies), but less accurate sampling of kinetic properties.  Because
configurational properties are much more important than kinetic ones in most
simulations, this integrator is generally preferred over LangevinIntegrator.  It
often allows one to use a larger time step while still maintaining similar or
better accuracy.

One disadvantage of this integrator is that it requires applying constraints
twice per time step, compared to only once for LangevinIntegrator.  This
can make it slightly slower for systems that involve constraints.  However, this
usually is more than compensated by allowing you to use a larger time step.

.. _nosehoover-integrators-theory:

NoseHooverIntegrator
********************

Like LangevinMiddleIntegerator, this uses the LFMiddle discretization.
:cite:`Zhang2019` In each step, the positions and velocities are updated as
follows:


.. math::
   \mathbf{v}_{i}(t+\Delta t/2) = \mathbf{v}_{i}(t-\Delta t/2) + \mathbf{f}_{i}(t)\Delta t/{m}_{i}


.. math::
   \mathbf{r}_{i}(t+\Delta t/2) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t+\Delta t/2)\Delta t/2


.. math::
   \mathbf{v'}_{i}(t+\Delta t/2) = \mathrm{scale}\times\mathbf{v}_{i}(t+\Delta t/2)


.. math::
   \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t+\Delta t/2) + \mathbf{v'}_{i}(t+\Delta t/2)\Delta t/2


The universal scale factor used in the third step is determined by propagating
auxilliary degrees of freedom alongside the regular particles.  The original
Nosé-Hoover formulation used a single harmonic oscillator for the heat bath,
but this is problematic in small or stiff systems, which are non-ergodic, so
the chain formulation extends this by replacing the single oscillator
thermostat with a chain of connected oscillators.  :cite:`Martyna1992`  For
large systems a single oscillator (*i.e.* a chain length of one) will suffice,
but longer chains are necessary to properly thermostat non-ergodic systems.
The OpenMM default is to use a chain length of three to cover the latter case,
but this can be safely reduced to increase efficiency in large systems.

The heat bath propagation is performed using a multi-timestep algorithm.  Each
propagation step is discretized into substeps using a factorization from
Yoshida and Suzuki; the default discretization uses a :math:`\mathcal{O}(\Delta
t^6)` approach that uses 7 points, but 1, 3 or 5 points may also be used to
increase performace, at the expense of accuracy.  Each step is further
subdivided into multi-timesteps with a default of 3 multi time steps per
propagation; as with the number of Yoshida-Suziki points this value may be
increase to increase accuracy but with additional computational expense.

BrownianIntegrator
******************

BrownianIntegrator simulates a system in contact with a heat bath by integrating
the Brownian equation of motion:


.. math::
   \frac{d\mathbf{r}_i}{dt}=\frac{1}{\gamma m_i}\mathbf{f}_i+\mathbf{R}_i


where :math:`\mathbf{r}_i` is the position of particle *i*\ , :math:`\mathbf{f}_i` is
the force acting on it, :math:`\gamma` is the friction coefficient, and :math:`\mathbf{R}_i`
is an uncorrelated random force whose components are chosen from a normal
distribution with mean zero and variance :math:`2 k_B T/m_i  \gamma`,
where *T* is the temperature of the heat bath.

The Brownian equation of motion is derived from the Langevin equation of motion
in the limit of large :math:`\gamma`\ .  In that case, the velocity of a particle is
determined entirely by the instantaneous force acting on it, and kinetic energy
ceases to have much meaning, since it disappears as soon as the applied force is
removed.


VariableVerletIntegrator
************************

This is very similar to VerletIntegrator, but instead of using the same step
size for every time step, it continuously adjusts the step size to keep the
integration error below a user-specified tolerance.  It compares the positions
generated by Verlet integration with those that would be generated by an
explicit Euler integrator, and takes the difference between them as an estimate
of the integration error:


.. math::
   error={\left(\Delta t\right)}^{2}\sum _{i}\frac{|\mathbf{f}_{i}|}{m_i}


where :math:`\mathbf{f}_i` is the force acting on particle *i* and :math:`m_i`
is its mass.  (In practice, the error made by the Euler integrator is usually
larger than that made by the Verlet integrator, so this tends to overestimate
the true error.  Even so, it can provide a useful mechanism for step size
control.)

It then selects the value of :math:`\Delta t` that makes the error exactly equal the
specified error tolerance:


.. math::
   \Delta t=\sqrt{\frac{\delta}{\sum _{i}\frac{|\mathbf{f}_i|}{m_i}}}


where :math:`\delta` is the error tolerance.  This is the largest step that may be
taken consistent with the user-specified accuracy requirement.

(Note that the integrator may sometimes choose to use a smaller value for :math:`\Delta t`
than given above.  For example, it might restrict how much the step size
can grow from one step to the next, or keep the step size constant rather than
increasing it by a very small amount.  This behavior is not specified and may
vary between Platforms.  It is required, however, that :math:`\Delta t` never be larger
than the value given above.)

A variable time step integrator is generally superior to a fixed time step one
in both stability and efficiency.  It can take larger steps on average, but will
automatically reduce the step size to preserve accuracy and avoid instability
when unusually large forces occur.  Conversely, when each uses the same step
size on average, the variable time step one will usually be more accurate since
the time steps are concentrated in the most difficult areas of the trajectory.

Unlike a fixed step size Verlet integrator, variable step size Verlet is not
symplectic.  This means that for a given average step size, it will not conserve
energy as precisely over long time periods, even though each local region of the
trajectory is more accurate.  For this reason, it is most appropriate when
precise energy conservation is not important, such as when simulating a system
at constant temperature.  For constant energy simulations that must maintain the
energy accurately over long time periods, the fixed step size Verlet may be more
appropriate.

VariableLangevinIntegrator
**************************

This is similar to LangevinIntegrator, but it continuously adjusts the step size
using the same method as VariableVerletIntegrator.  It is usually preferred over
the fixed step size Langevin integrator for the reasons given above.
Furthermore, because Langevin dynamics involves a random force, it can never be
symplectic and therefore the fixed step size Verlet integrator’s advantages do
not apply to the Langevin integrator.

.. _custom-integrator:

CustomIntegrator
****************

CustomIntegrator is a very flexible class that can be used to implement a wide
range of integration methods.  This includes both deterministic and stochastic
integrators; Metropolized integrators; multiple time step integrators; and
algorithms that must integrate additional quantities along with the particle
positions and momenta.

The algorithm is specified as a series of computations that are executed in
order to perform a single time step.  Each computation computes the value (or
values) of a *variable*\ .  There are two types of variables: *global
variables* have a single value, while *per-DOF variables* have a separate
value for every degree of freedom (that is, every *x*\ , *y*\ , or *z*
component of a particle).  CustomIntegrator defines lots of variables you can
compute and/or use in computing other variables.  Some examples include the step
size (global), the particle positions (per-DOF), and the force acting on each
particle (per-DOF).  In addition, you can define as many variables as you want
for your own use.

The actual computations are defined by mathematical expressions as described in
section :numref:`writing-custom-expressions`\ .  Several types of computations are supported:

* *Global*\ : the expression is evaluated once, and the result is stored into
  a global variable.
* *Per-DOF*\ : the expression is evaluated once for every degree of freedom,
  and the results are stored into a per-DOF variable.
* *Sum*\ : the expression is evaluated once for every degree of freedom.  The
  results for all degrees of freedom are added together, and the sum is stored
  into a global variable.


There also are other, more specialized types of computations that do not involve
mathematical expressions.  For example, there are computations that apply
distance constraints, modifying the particle positions or velocities
accordingly.

CustomIntegrator is a very powerful tool, and this description only gives a
vague idea of the scope of its capabilities.  For full details and examples,
consult the API documentation.

.. _the-theory-behind-openmm-introduction:

Introduction
############

Overview
********

This guide describes the mathematical theory behind OpenMM.  For each
computational class, it describes what computations the class performs and how
it should be used.  This serves two purposes.  If you are using OpenMM within an
application, this guide teaches you how to use it correctly.  If you are
implementing the OpenMM API for a new Platform, it teaches you how to correctly
implement the required kernels.

On the other hand, many details are intentionally left unspecified.  Any
behavior that is not specified either in this guide or in the API documentation
is left up to the Platform, and may be implemented in different ways by
different Platforms.  For example, an Integrator is required to produce a
trajectory that satisfies constraints to within the user-specified tolerance,
but the algorithm used to enforce those constraints is left up to the Platform.
Similarly, this guide provides the functional form of each Force, but does not
specify what level of numerical precision it must be calculated to.

This is an essential feature of the design of OpenMM, because it allows the API
to be implemented efficiently on a wide variety of hardware and software
platforms, using whatever methods are most appropriate for each platform.  On
the other hand, it means that a single program may produce meaningfully
different results depending on which Platform it uses.  For example, different
constraint algorithms may have different regions of convergence, and thus a time
step that is stable on one platform may be unstable on a different one.  It is
essential that you validate your simulation methodology on each Platform you
intend to use, and do not assume that good results on one Platform will
guarantee good results on another Platform when using identical parameters.


.. _units:

Units
*****

There are several different sets of units widely used in molecular simulations.
For example, energies may be measured in kcal/mol or kJ/mol, distances may be in
Angstroms or nm, and angles may be in degrees or radians.  OpenMM uses the
following units everywhere.

===========  =================
Quantity     Units
===========  =================
distance     nm
time         ps
mass         atomic mass units
charge       proton charge
temperature  Kelvin
angle        radians
energy       kJ/mol
===========  =================

These units have the important feature that they form an internally consistent
set.  For example, a force always has the same units (kJ/mol/nm) whether it is
calculated as the gradient of an energy or as the product of a mass and an
acceleration.  This is not true in some other widely used unit systems, such as
those that express energy in kcal/mol.

The header file Units.h contains predefined constants for converting between the
OpenMM units and some other common units.  For example, if your application
expresses distances in Angstroms, you should multiply them by
OpenMM::NmPerAngstrom before passing them to OpenMM, and positions calculated by
OpenMM should be multiplied by OpenMM::AngstromsPerNm before passing them back
to your application.

.. _physical-constants:

Physical Constants
******************

OpenMM uses the CODATA 2018 values for all physical constants.  Here are the
specific values it uses for the constants that frequently come up in molecular
simulations.

========================================  =================================
Quantity                                  Value
========================================  =================================
Elementary Charge (:math:`e`)             1.602176634·10\ :sup:`-19` C
Boltzmann's Constant (:math:`k_B`)        1.380649·10\ :sup:`-23` J/K
Avogadro's Number (:math:`N_A`)           6.02214076·10\ :sup:`23`
Vacuum Permittivity (:math:`\epsilon_0`)  8.8541878128·10\ :sup:`-12` F/m
========================================  =================================

.. _openmm-tutorials:

OpenMM Tutorials
################


Example Files Overview
**********************

Four example files are provided in the examples folder, each designed with
a specific objective.

* **HelloArgon:**  A very simple example intended for verifying that you
  have installed OpenMM correctly.  It also introduces you to the basic classes
  within OpenMM.
* **HelloSodiumChloride:**  This example shows you our recommended strategy
  for integrating OpenMM into an existing molecular dynamics code.
* **HelloEthane:** The main purpose of this example is to demonstrate how
  to tell OpenMM about bonded forces (bond stretch, bond angle bend, dihedral
  torsion).
* **HelloWaterBox:**  This example shows you how to use OpenMM to model
  explicit solvation, including setting up periodic boundary conditions.  It runs
  extremely fast on a GPU but very, very slowly on a CPU, so it is an excellent
  example to use to compare performance on the GPU versus the CPU.  The other
  examples provided use systems where the performance difference would be too
  small to notice.


The two fundamental examples—HelloArgon and HelloSodiumChloride—are provided in
C++, C, and Fortran, as indicated in the table below.  The other two
examples—HelloEthane and HelloWaterBox—follow the same structure as
HelloSodiumChloride but demonstrate more calls within the OpenMM API.  They are
only provided in C++ but can be adapted to run in C and Fortran by following the
mappings described in Chapter :numref:`using-openmm-with-software-written-in-languages-other-than-c++`\ .
HelloArgon and HelloSodiumChloride also serve as examples of how to do these mappings.  The
sections below describe the HelloArgon, HelloSodiumChloride, and HelloEthane programs in more detail.

===============  ==============  ==========  ========  ========================================  ===============
Example          Solvent         Thermostat  Boundary  Forces & Constraints                      API
===============  ==============  ==========  ========  ========================================  ===============
Argon            Vacuum          None        None      Non-bonded\*                              C++, C, Fortran
Sodium Chloride  Implicit water  Langevin    None      Non-bonded\*                              C++, C, Fortran
Ethane           Vacuum          None        None      Non-bonded\*, stretch, bend, torsion      C++
Water Box        Explicit water  Andersen    Periodic  Non-bonded\*, stretch, bend, constraints  C++
===============  ==============  ==========  ========  ========================================  ===============

\*van der Waals and Coulomb forces

.. _running-example-files:

Running Example Files
**********************

The instructions below are for running the HelloArgon program.  A similar
process would be used to run the other examples.

Visual Studio
=============

Navigate to wherever you saved the example files.  Descend into the directory
folder VisualStudio. Double-click the file HelloArgon.sln (a Microsoft Visual
Studio Solution file).  Visual Studio will launch.

Note: These files were created using Visual Studio 8.  If you are using a more
recent version, it will ask if you want to convert the files to the new version.
Agree and continue through the conversion process.

In Visual Studio, make sure the "Solution Configuration" is set to "Release" and
not "Debug".  The “Solution Configuration” can be set using the drop-down menu
in the top toolbar, next to the green arrow (see :autonumref:`Figure,Visual Studio configuration`
below).  Due to incompatibilities among Visual Studio versions, we do not provide pre-compiled
debug binaries.



.. figure:: ../../images/VisualStudioSetConfiguration.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,Visual Studio configuration`:  Setting "Solution Configuration" to "Release" mode in Visual Studio




From the command options select Debug -> Start Without Debugging (or CTRL-F5).
See :autonumref:`Figure,run in Visual Studio`.  This will also compile the program, if it has not
previously been compiled.



.. figure:: ../../images/VisualStudioLaunch.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,run in Visual Studio`:  Run a program in Visual Studio

You should see a series of lines like the following output on your screen:
::

    REMARK  Using OpenMM platform Reference
    MODEL     1
    ATOM      1  AR   AR     1       0.000   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.000   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       10.000  0.000   0.000  1.00  0.00
    ENDMDL

    …

    MODEL     250
    ATOM      1  AR   AR     1       0.233   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.068   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.678   0.000   0.000  1.00  0.00
    ENDMDL
    MODEL     251
    ATOM      1  AR   AR     1       0.198   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.082   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.698   0.000   0.000  1.00  0.00
    ENDMDL
    MODEL     252
    ATOM      1  AR   AR     1       0.165   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.097   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.717   0.000   0.000  1.00  0.00
    ENDMDL


Determining the platform being used
-----------------------------------

The very first line of the output will indicate whether you are running on the
CPU (Reference platform) or a GPU (CUDA or OpenCL platform).  It will say one of
the following:
::

    REMARK  Using OpenMM platform Reference
    REMARK  Using OpenMM platform Cuda
    REMARK  Using OpenMM platform OpenCL

If you have a supported GPU, the program should, by default, run on the GPU.

Visualizing the results
------------------------

You can output the results to a PDB file that could be visualized using programs
like VMD (http://www.ks.uiuc.edu/Research/vmd/) or PyMol
(http://pymol.sourceforge.net/).  To do this within Visual Studios:

#. Right-click on the project name HelloArgon (not one of the files) and select
   the “Properties” option.
#. On the “Property Pages” form, select “Debugging” under the “Configuration
   Properties” node.
#. In the “Command Arguments” field, type:

   ::

       > argon.pdb

   This will save the output to a file called argon.pdb in the current working
   directory (default is the VisualStudio directory).  If you want to save it to
   another directory, you will need to specify the full path.

#. Select “OK”


Now, when you run the program in Visual Studio, no text will appear.  After a
short time, you should see the message “\ :code:`Press any key to continue…`\ ,”
indicating that the program is complete and that the PDB file has been
completely written.

Mac OS X/Linux
==============

Navigate to wherever you saved the example files.

Verify your makefile by consulting the MakefileNotes file in this directory, if
necessary.

Type:::

    make


Then run the program by typing:
::

    ./HelloArgon

You should see a series of lines like the following output on your screen:
::

    REMARK  Using OpenMM platform Reference
    MODEL     1
    ATOM      1  AR   AR     1       0.000   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.000   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       10.000  0.000   0.000  1.00  0.00
    ENDMDL

    ...

    MODEL     250
    ATOM      1  AR   AR     1       0.233   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.068   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.678   0.000   0.000  1.00  0.00
    ENDMDL
    MODEL     251
    ATOM      1  AR   AR     1       0.198   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.082   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.698   0.000   0.000  1.00  0.00
    ENDMDL
    MODEL     252
    ATOM      1  AR   AR     1       0.165   0.000   0.000  1.00  0.00
    ATOM      2  AR   AR     1       5.097   0.000   0.000  1.00  0.00
    ATOM      3  AR   AR     1       9.717   0.000   0.000  1.00  0.00
    ENDMDL


Determining the platform being used
-----------------------------------

The very first line of the output will indicate whether you are running on the
CPU (Reference platform) or a GPU (CUDA or OpenCL platform).  It will say one of
the following:
::

    REMARK  Using OpenMM platform Reference
    REMARK  Using OpenMM platform Cuda
    REMARK  Using OpenMM platform OpenCL

If you have a supported GPU, the program should, by default, run on the GPU.

Visualizing the results
------------------------

You can output the results to a PDB file that could be visualized using programs
like VMD (http://www.ks.uiuc.edu/Research/vmd/) or PyMol
(http://pymol.sourceforge.net/) by typing:
::

    ./HelloArgon > argon.pdb

Compiling Fortran and C examples
--------------------------------

The Makefile provided with the examples can also be used to compile the Fortran
and C examples.

The Fortran compiler needs to load a version of the libstdc++.dylib library that
is compatible with the version of gcc used to build OpenMM;   OpenMM for Mac is
compiled using gcc 4.2.  If you are compiling with a different version, edit the
Makefile and add the following flag to FCPPLIBS: :code:`–L/usr/lib/gcc/i686
-apple-darwin10/4.2.1`\ .

When the Makefile has been updated, type:
::

    make all

HelloArgon Program
******************

The HelloArgon program simulates three argon atoms in a vacuum.  It is a simple
program primarily intended for you to verify that you are able to compile, link,
and run with OpenMM.  It also demonstrates the basic calls needed to run a
simulation using OpenMM.

Including OpenMM-defined functions
==================================

The OpenMM header file *OpenMM.h* instructs the program to include
everything defined by the OpenMM libraries.  Include the header file by adding
the following line at the top of your program:  ::


    #include "OpenMM.h"

Running a program on GPU platforms
==================================

By default, a program will run on the Reference platform.  In order to run a
program on another platform (e.g., an NVIDIA or AMD GPU), you need to load the
required shared libraries for that other platform (e.g., Cuda, OpenCL).  The
easy way to do this is to call:

.. code-block:: c

    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());

This will load all the shared libraries (plug-ins) that can be found, so you do
not need to explicitly know which libraries are available on a given machine.
In this way, the program will be able to run on another platform, if it is
available.

Running a simulation using the OpenMM public API
================================================

The OpenMM public API was described in Section :numref:`the-openmm-public-api`\ .  Here you will
see how to use those classes to create a simple system of three argon atoms and run a short
simulation.  The main components of the simulation are within the function
:code:`simulateArgon()`\ :

#. **System** – We first establish a system and add a non-bonded force to
   it.  At this point, there are no particles in the system.

   .. code-block:: c

        // Create a system with nonbonded forces.
        OpenMM::System system;
        OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
        system.addForce(nonbond);

   We then add the three argon atoms to the system.  For this system, all the data
   for the particles are hard-coded into the program.  While not a realistic
   scenario, it makes the example simpler and clearer.  The
   :code:`std::vector<OpenMM::Vec3>` is an array of vectors of 3.

   .. code-block:: c

        // Create three atoms.
        std::vector<OpenMM::Vec3> initPosInNm(3);
        for (int a = 0; a < 3; ++a)
        {
            initPosInNm[a] = OpenMM::Vec3(0.5*a,0,0); // location, nm

            system.addParticle(39.95); // mass of Ar, grams per mole

            // charge, L-J sigma (nm), well depth (kJ)
            nonbond->addParticle(0.0, 0.3350, 0.996); // vdWRad(Ar)=.188 nm
        }

   **Units:** Be very careful with the units in your program.  It is very easy
   to make mistakes with the units, so we recommend including them in your variable
   names, as we have done here :code:`initPosInNm` (position in nanometers).
   OpenMM provides conversion constants that should be used whenever there are
   conversions to be done; for simplicity, we did not do that in HelloArgon, but
   all the other examples show the use of these constants.

   It is hard to overemphasize the importance of careful units handling—it is very
   easy to make a mistake despite, or perhaps because of, the trivial nature of
   units conversion.  For more information about the units used in OpenMM, see
   Section :numref:`units`.

   **Adding Particle Information:** Both the system and the non-bonded
   force require information about the particles.  The system just needs to know
   the mass of the particle.  The non-bonded force requires information about the
   charge (in this case, argon is uncharged), and the Lennard-Jones parameters
   sigma (zero-energy separation distance) and well depth (see Section :numref:`lennard-jones-interaction`
   for more details).

   Note that the van der Waals radius for argon is 0.188 nm and that it has already
   been converted to sigma (0.335 nm) in the example above where it is added to the
   non-bonded force;  in your code, you should make use of the appropriate
   conversion factor supplied with OpenMM as discussed in Section :numref:`units`\ .

#. **Integrator** – We next specify the integrator to use to perform the
   calculations.  In this case, we choose a Verlet integrator to run a constant
   energy simulation.  The only argument required is the step size in picoseconds.

   .. code-block:: c

        OpenMM::VerletIntegrator integrator(0.004); // step size in ps

   We have chosen to use 0.004 picoseconds, or 4 femtoseconds, which is larger than
   that used in a typical molecular dynamics simulation.  However, since this
   example does not have any bonds with higher frequency components, like most
   molecular dynamics simulations do, this is an acceptable value.

#. **Context** – The context is an object that consists of an integrator and
   a system.  It manages the state of the simulation.  The code below initializes
   the context.  We then let the context select the best platform available to run
   on, since this is not specifically specified, and print out the chosen platform.
   This is useful information, especially when debugging.

   .. code-block:: c

        // Let OpenMM Context choose best platform.
        OpenMM::Context context(system, integrator);
        printf("REMARK  Using OpenMM platform %s\n", context.getPlatform().getName().c_str());

   We then initialize the system, setting the initial time, as well as the initial
   positions and velocities of the atoms.  In this example, we leave time and
   velocity at their default values of zero.

   .. code-block:: c

        // Set starting positions of the atoms. Leave time and velocity zero.
        context.setPositions(initPosInNm);

#. **Initialize and run the simulation** – The next block of code runs the
   simulation and saves its output.  For each frame of the simulation (in this
   example, a frame is defined by the advancement interval of the integrator; see
   below), the current state of the simulation is obtained and written out to a
   PDB-formatted file.

   .. code-block:: c

        // Simulate.
        for (int frameNum=1; ;++frameNum) {
            // Output current state information.
            OpenMM::State state = context.getState(OpenMM::State::Positions);
            const double  timeInPs = state.getTime();
            writePdbFrame(frameNum, state); // output coordinates

   *Getting state information has to be done in bulk, asking for information for
   all the particles at once.*  This is computationally expensive since this
   information can reside on the GPUs and requires communication overhead to
   retrieve, so you do not want to do it very often.  In the above code, we only
   request the positions, since that is all that is needed, and time from the
   state.

   The simulation stops after 10 ps; otherwise we ask the integrator to take 10
   steps (so one frame is equivalent to 10 time steps).   Normally, we would want
   to take more than 10 steps at a time, but to get a reasonable-looking animation,
   we use 10.

   .. code-block:: c

         if (timeInPs >= 10.)
             break;

         // Advance state many steps at a time, for efficient use of OpenMM.
         integrator.step(10); // (use a lot more than this normally)

Error handling for OpenMM
=========================

Error handling for OpenMM is explicitly designed so you do not have to check the
status after every call.  If anything goes wrong, OpenMM throws an exception.
It uses standard exceptions, so on many platforms, you will get the exception
message automatically.  However, we recommend using :code:`try-catch` blocks
to ensure you do catch the exception.

.. code-block:: c

    int main()
    {
        try {
            simulateArgon();
            return 0; // success!
        }
        // Catch and report usage and runtime errors detected by OpenMM and fail.
        catch(const std::exception& e) {
            printf("EXCEPTION: %s\n", e.what());
            return 1; // failure!
        }
    }

Writing out PDB files
=====================

For the HelloArgon program, we provide a simple PDB file writing function
:code:`writePdbFrame` that *only* writes out argon atoms.  The function
has nothing to do with OpenMM except for using the OpenMM State.  The function
extracts the positions from the State in nanometers (10\ :sup:`-9` m) and
converts them to Angstroms (10\ :sup:`-10` m) to be compatible with the PDB
format.   Again, we emphasize how important it is to track the units being used!

.. code-block:: c

    void writePdbFrame(int frameNum, const OpenMM::State& state)
    {
        // Reference atomic positions in the OpenMM State.
        const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

        // Use PDB MODEL cards to number trajectory frames
        printf("MODEL     %d\n", frameNum); // start of frame
        for (int a = 0; a < (int)posInNm.size(); ++a)
        {
            printf("ATOM  %5d  AR   AR     1    ", a+1); // atom number
            printf("%8.3f%8.3f%8.3f  1.00  0.00\n",      // coordinates
            // "*10" converts nanometers to Angstroms
            posInNm[a][0]*10, posInNm[a][1]*10, posInNm[a][2]*10);
        }
        printf("ENDMDL\n"); // end of frame
    }

:code:`MODEL` and :code:`ENDMDL` are used to mark the beginning and end
of a frame, respectively.  By including multiple frames in a PDB file, you can
visualize the simulation trajectory.

HelloArgon output
=================

The output of the HelloArgon program can be saved to a *.pdb* file and
visualized using programs like VMD or PyMol (see Section :numref:`running-example-files`).
You should see three atoms moving linearly away and towards one another:


.. figure:: ../../images/Argon.png
   :align: center


You may need to adjust the van der Waals radius in your visualization program to
see the atoms colliding.

HelloSodiumChloride Program
***************************

The HelloSodiumChloride models several sodium (Na\ :sup:`+`\ ) and chloride
(Cl\ :sup:`-`\ ) ions in implicit solvent (using a Generalized Born/Surface Area, or
GBSA, OBC model).  As with the HelloArgon program, only non-bonded forces are
simulated.

The main purpose of this example is to illustrate our recommended strategy for
integrating OpenMM into an existing molecular dynamics (MD) code:

#. **Write a few, high-level interface routines containing all your OpenMM
   calls**\ :  Rather than make OpenMM calls throughout your program, we
   recommend writing a handful of interface routines that understand both your MD
   code’s data structures and OpenMM.  Organize these routines into a separate
   compilation unit so you do not have to make huge changes to your existing MD
   code.  These routines could be written in any language that is callable from the
   existing MD code.  We recommend writing them in C++ since that is what OpenMM is
   written in, but you can also write them in C or Fortran; see Chapter
   :numref:`using-openmm-with-software-written-in-languages-other-than-c++`\ .


#. **Call only these high-level interface routines from your existing MD
   code:**  This provides a clean separation between the existing MD code and
   OpenMM, so that changes to OpenMM will not directly impact the existing MD code.
   One way to implement this is to use opaque handles, a standard trick used (for
   example) for opening files in Linux.  An existing MD code can communicate with
   OpenMM via the handle, but knows none of the details of the handle.  It only has
   to hold on to the handle and give it back to OpenMM.


In the example described below, you will see how this strategy can be
implemented for a very simple MD code.  Chapter :numref:`examples-of-openmm-integration`
describes the strategies used in integrating OpenMM into real MD codes.

.. _simple-molecular-dynamics-system:

Simple molecular dynamics system
================================

The initial sections of HelloSodiumChloride.cpp represent a very simple
molecular dynamics system.  The system includes modeling and simulation
parameters and the atom and force field data.  It also provides a data structure
\ :code:`posInAng[3]` for storing the current state.  These sections represent
(in highly simplified form) information that would be available from an existing
MD code, and will be used to demonstrate how to integrate OpenMM with an
existing MD program.

.. code-block:: c

    // -----------------------------------------------------------------
    //                   MODELING AND SIMULATION PARAMETERS
    // -----------------------------------------------------------------
    static const double Temperature         = 300;     // Kelvins
    static const double FrictionInPerPs     = 91.;     // collisions per picosecond
    static const double SolventDielectric   = 80.;     // typical for water
    static const double SoluteDielectric    = 2.;      // typical for protein

    static const double StepSizeInFs        = 4;       // integration step size (fs)
    static const double ReportIntervalInFs  = 50;      // how often to issue PDB frame (fs)
    static const double SimulationTimeInPs  = 100;     // total simulation time (ps)

    // Decide whether to request energy calculations.
    static const bool   WantEnergy          = true;


    // -----------------------------------------------------------------
    //                          ATOM AND FORCE FIELD DATA
    // -----------------------------------------------------------------
    // This is not part of OpenMM; just a struct we can use to collect atom
    // parameters for this example. Normally atom parameters would come from the
    // force field's parameterization file. We're going to use data in Angstrom and
    // Kilocalorie units and show how to safely convert to OpenMM's internal unit
    // system which uses nanometers and kilojoules.
    static struct MyAtomInfo {
        const char* pdb;
        double      mass, charge, vdwRadiusInAng, vdwEnergyInKcal,
                    gbsaRadiusInAng, gbsaScaleFactor;
        double      initPosInAng[3];
        double      posInAng[3]; // leave room for runtime state info
    } atoms[] = {
    // pdb   mass  charge  vdwRad vdwEnergy   gbsaRad gbsaScale  initPos
    {" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     8, 0,  0},
    {" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,    -8, 0,  0},
    {" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 9,  0},
    {" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0,-9,  0},
    {" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 0,-10},
    {" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0, 0, 10},
    {""} // end of list
    };


Interface routines
==================

The key to our recommended integration strategy is the interface routines.  You
will need to decide what interface routines are required for effective
communication between your existing MD program and OpenMM, but typically there
will only be six or seven.  In our example, the following four routines suffice:

* **Initialize:** Data structures that already exist in your MD program
  (i.e., force fields, constraints, atoms in the system) are passed to the
  :code:`Initialize` routine, which makes appropriate calls to OpenMM and then
  returns a handle to the OpenMM object that can be used by the existing MD
  program.
* **Terminate:** Clean up the heap space allocated by :code:`Initialize`
  by passing the handle to the :code:`Terminate` routine.
* **Advance State:** The :code:`AdvanceState` routine advances the
  simulation.  It requires that the calling function, the existing MD code, gives
  it a handle.
* **Retrieve State:** When you want to do an analysis or generate some kind
  of report, you call the :code:`RetrieveState` routine.  You have to give it
  a handle.  It then fills in a data structure that is defined in the existing MD
  code, allowing the MD program to use it in its existing routines without further
  modification.

Note that these are just descriptions of the routines’ functions—you can call
them anything you like and implement them in whatever way makes sense for your
MD code.

In the example code, the four routines performing these functions, plus an
opaque data structure (the handle), would be declared, as shown below.  Then,
the main program, which sets up, runs, and reports on the simulation, accesses
these routines and the opaque data structure (in this case, the variable
:code:`omm`\ ).  As you can see, it does not have access to any OpenMM
declarations, only to the interface routines that you write so there is no need
to change the build environment.

.. code-block:: c

    struct MyOpenMMData;
    static MyOpenMMData* myInitializeOpenMM(const MyAtomInfo atoms[],
                                            double temperature,
                                            double frictionInPs,
                                            double solventDielectric,
                                            double soluteDielectric,
                                            double stepSizeInFs,
                                            std::string& platformName);
    static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
    static void          myGetOpenMMState(MyOpenMMData*,
                                          bool wantEnergy,
                                          double& time,
                                          double& energy,
                                          MyAtomInfo atoms[]);
    static void          myTerminateOpenMM(MyOpenMMData*);


    // -----------------------------------------------------------------
    //                                MAIN PROGRAM
    // -----------------------------------------------------------------
    int main() {
        const int NumReports     = (int)(SimulationTimeInPs*1000 / ReportIntervalInFs + 0.5);
        const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);

        // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
        // usage and runtime errors are caught and reported.
        try {
            double        time, energy;
            std::string   platformName;

            // Set up OpenMM data structures; returns OpenMM Platform name.
            MyOpenMMData* omm = myInitializeOpenMM(atoms, Temperature, FrictionInPerPs,
                 SolventDielectric, SoluteDielectric, StepSizeInFs, platformName);

            // Run the simulation:
            //  (1) Write the first line of the PDB file and the initial configuration.
            //  (2) Run silently entirely within OpenMM between reporting intervals.
            //  (3) Write a PDB frame when the time comes.
            printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
            myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
            myWritePDBFrame(1, time, energy, atoms);

            for (int frame=2; frame <= NumReports; ++frame) {
                myStepWithOpenMM(omm, NumSilentSteps);
                myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
                myWritePDBFrame(frame, time, energy, atoms);
            }

            // Clean up OpenMM data structures.
            myTerminateOpenMM(omm);

            return 0; // Normal return from main.
        }

        // Catch and report usage and runtime errors detected by OpenMM and fail.
        catch(const std::exception& e) {
            printf("EXCEPTION: %s\n", e.what());
            return 1;
        }
    }

We will examine the implementation of each of the four interface routines and
the opaque data structure (handle) in the sections below.

Units
-----

The simple molecular dynamics system described in Section :numref:`simple-molecular-dynamics-system`
employs the commonly used units of angstroms and kcals.  These differ from the units and
parameters used within OpenMM (see Section :numref:`units`\ ): nanometers and kilojoules.
These differences may be small but they are critical and must be carefully
accounted for in the interface routines.

Lennard-Jones potential
-----------------------

The Lennard-Jones potential describes the energy between two identical atoms as
the distance between them varies.

The van der Waals “size” parameter is used to identify the distance at which the
energy between these two atoms is at a minimum (that is, where the van der Waals
force is most attractive).  There are several ways to specify this parameter,
typically, either as the van der Waals radius r\ :sub:`vdw` or as the actual
distance between the two atoms d\ :sub:`min` (also called r\ :sub:`min`\ ),
which is twice the van der Waals radius r\ :sub:`vdw`\ .  A third way to
describe the potential is through sigma :math:`\sigma`, which identifies the distance at
which the energy function crosses zero as the atoms move closer together than
d\ :sub:`min`\ .  (See Section :numref:`lennard-jones-interaction` for more details about the
relationship between these).

:math:`\sigma` turns out to be about 0.89*d\ :sub:`min`\ , which is close enough to
d\ :sub:`min` that it makes it hard to distinguish the two.  Be very careful that
you use the correct value.  In the example below, we will show you how to use
the built-in OpenMM conversion constants to avoid errors.

Lennard-Jones parameters are defined for pairs of identical atoms, but must also
be applied to pairs of dissimilar atoms. That is done by “combining rules” that
differ among popular MD codes. Two of the most common are:

* Lorentz-Berthelot (used by AMBER, CHARMM):

.. math::
    r=\frac{r_i+r_j}{2}, \epsilon=\sqrt{\epsilon_i \epsilon_j}

* Jorgensen (used by OPLS):

.. math::
    r=\sqrt{r_i r_j}, \epsilon=\sqrt{\epsilon_i \epsilon_j}


where *r* = the effective van der Waals “size” parameter (minimum radius,
minimum distance, or zero crossing (sigma)), and :math:`\epsilon` = the effective van
der Waals energy well depth parameter, for the dissimilar pair of atoms *i*
and *j*\ .

OpenMM only implements Lorentz-Berthelot directly, but others can be implemented
using the CustomNonbondedForce class.  (See Section :numref:`customnonbondedforce` for details.)

Opaque handle MyOpenMMData
--------------------------

In this example, the handle used by the interface to OpenMM is a pointer to a
struct called :code:`MyOpenMMData.`  The pointer itself is opaque, meaning
the calling program has no knowledge of what the layout of the object it points
to is, or how to use it to directly interface with OpenMM.  The calling program
will simply pass this opaque handle from one interface routine to another.

There are many different ways to implement the handle.  The code below shows
just one example.  A simulation requires three OpenMM objects (a System, a
Context, and an Integrator) and so these must exist within the handle.  If other
objects were required for a simulation, you would just add them to your handle;
there would be no change in the main program using the handle.

.. code-block:: c

    struct MyOpenMMData {
        MyOpenMMData() : system(0), context(0), integrator(0) {}
        ~MyOpenMMData() {delete system; delete context; delete integrator;}
        OpenMM::System*         system;
        OpenMM::Context*        context;
        OpenMM::Integrator*     integrator;
    };

In addition to establishing pointers to the required three OpenMM objects,
:code:`MyOpenMMData` has a constructor :code:`MyOpenMMData()` that sets
the pointers for the three OpenMM objects to zero and a destructor
:code:`~MyOpenMMData()` that (in C++) gives the heap space back.  This was
done in-line in the HelloArgon program, but we recommend you use something like
the method here instead.

myInitializeOpenMM
-------------------

The :code:`myInitializeOpenMM` function takes the data structures and
simulation parameters from the existing MD code and returns a new handle that
can be used to do efficient computations with OpenMM.  It also returns the
:code:`platformName` so the calling program knows what platform (e.g., CUDA,
OpenCL, Reference) was used.

.. code-block:: c

    static MyOpenMMData*
    myInitializeOpenMM( const MyAtomInfo    atoms[],
                        double              temperature,
                        double              frictionInPs,
                        double              solventDielectric,
                        double              soluteDielectric,
                        double              stepSizeInFs,
                        std::string&        platformName)


This initialization routine is very similar to the HelloArgon example program,
except that objects are created and put in the handle.  For instance, just as in
the HelloArgon program, the first step is to load the OpenMM plug-ins, so that
the program will run on the best performing platform that is available.   Then,
a System is created **and** assigned to the handle :code:`omm`\ .
Similarly, forces are added to the System which is already in the handle.

.. code-block:: c

    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
           (OpenMM::Platform::getDefaultPluginsDirectory());

    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the OpenMM
    // System takes ownership of the force objects;don't delete them yourself.
    omm->system = new OpenMM::System();
    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
    OpenMM::GBSAOBCForce*   gbsa    = new OpenMM::GBSAOBCForce();
    omm->system->addForce(nonbond);
    omm->system->addForce(gbsa);

    // Specify dielectrics for GBSA implicit solvation.
    gbsa->setSolventDielectric(solventDielectric);
    gbsa->setSoluteDielectric(soluteDielectric);


In the next step, atoms are added to the System within the handle, with
information about each atom coming from the data structure that was passed into
the initialization function from the existing MD code.  As shown in the
HelloArgon program, both the System and the forces need information about the
atoms.  For those unfamiliar with the C++ Standard Template Library, the
:code:`push_back` function called at the end of this code snippet just adds
the given argument to the end of a C++ “vector” container.

.. code-block:: c

    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties(in MD units!).
    //  (3) GBSA needs charge, radius, and scale factor.
    //  (4) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; *atoms[n].pdb; ++n) {
         const MyAtomInfo& atom = atoms[n];

         omm->system->addParticle(atom.mass);

         nonbond->addParticle(atom.charge,
                             atom.vdwRadiusInAng * OpenMM::NmPerAngstrom
                                                 * OpenMM::SigmaPerVdwRadius,
                             atom.vdwEnergyInKcal * OpenMM::KJPerKcal);

         gbsa->addParticle(atom.charge,
                           atom.gbsaRadiusInAng * OpenMM::NmPerAngstrom,
                           atom.gbsaScaleFactor);

         // Convert the initial position to nm and append to the array.
         const Vec3 posInNm(atom.initPosInAng[0] * OpenMM::NmPerAngstrom,
                      atom.initPosInAng[1] * OpenMM::NmPerAngstrom,
                      atom.initPosInAng[2] * OpenMM::NmPerAngstrom);
         initialPosInNm.push_back(posInNm);


**Units:**  Here we emphasize the need to pay special attention to the
units.   As mentioned earlier, the existing MD code in this example uses units
of angstroms and kcals, but OpenMM uses nanometers and kilojoules.  So the
initialization routine will need to convert the values from the existing MD code
into the OpenMM units before assigning them to the OpenMM objects.

In the code above, we have used the unit conversion constants that come with
OpenMM (e.g., :code:`OpenMM::NmPerAngstrom`\ ) to perform these conversions.
Combined with the naming convention of including the units in the variable name
(e.g., :code:`initPosInAng`\ ), the unit conversion constants are useful
reminders to pay attention to units and minimize errors.

Finally, the initialization routine creates the Integrator and Context for the
simulation.  Again, note the change in units for the arguments!   The routine
then gets the platform that will be used to run the simulation and returns that,
along with the handle :code:`omm`\ , back to the calling function.

.. code-block:: c

    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero but could
    // have been set here.
    omm->integrator = new OpenMM::LangevinMiddleIntegrator(temperature,
    frictionInPs,
    stepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initialPosInNm);

    platformName = omm->context->getPlatform().getName();
    return omm;


myGetOpenMMState
----------------

The :code:`myGetOpenMMState` function takes the handle and returns the time,
energy, and data structure for the atoms in a way that the existing MD code can
use them without modification.

.. code-block:: c

    static void
    myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy,
                     double& timeInPs, double& energyInKcal, MyAtomInfo atoms[])

Again, this is another interface routine in which you need to be very careful of
your units!  Note the conversion from the OpenMM units back to the units used in
the existing MD code.

.. code-block:: c

    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
       infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
       infoMask += OpenMM::State::Energy;     // for pot. energy (more expensive)
    }
    // Forces are also available (and cheap).

    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
             atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;

    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;
    if (wantEnergy)
       energyInKcal = (state.getPotentialEnergy() + state.getKineticEnergy())
                      * OpenMM::KcalPerKJ;

myStepWithOpenMM
----------------

The :code:`myStepWithOpenMM` routine takes the handle, uses it to find the
Integrator, and then sets the number of steps for the Integrator to take.  It
does not return any values.

.. code-block:: c

    static void
    myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
        omm->integrator->step(numSteps);
    }

myTerminateOpenMM
-----------------

The :code:`myTerminateOpenMM` routine takes the handle and deletes all the
components, e.g., the Context and System, cleaning up the heap space.

.. code-block:: c

    static void
    myTerminateOpenMM(MyOpenMMData* omm) {
        delete omm;
    }


HelloEthane Program
*******************

The HelloEthane program simulates ethane (H3-C-C-H3) in a vacuum.  It is
structured similarly to the HelloSodiumChloride example, but includes bonded
forces (bond stretch, bond angle bend, dihedral torsion).  In setting up these
bonded forces, the program illustrates some of the other inconsistencies in
definitions and units that you should watch out for.

The bonded forces are added to the system within the initialization interface
routine, similar to how the non-bonded forces were added in the
HelloSodiumChloride example:

.. code-block:: c

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    OpenMM::NonbondedForce&         nonbond     = *new OpenMM::NonbondedForce();
    OpenMM::HarmonicBondForce&      bondStretch = *new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicAngleForce&     bondBend    = *new OpenMM::HarmonicAngleForce();
    OpenMM::PeriodicTorsionForce&   bondTorsion = *new OpenMM::PeriodicTorsionForce();
    system.addForce(&nonbond);
    system.addForce(&bondStretch);
    system.addForce(&bondBend);
    system.addForce(&bondTorsion);

\ **Constrainable and non-constrainable bonds:**  In the initialization
routine, we also set up the bonds.  If constraints are being used, then we tell
the System about the constrainable bonds:

.. code-block:: c

    std::vector< std::pair<int,int> > bondPairs;
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atom = bonds[i].atoms;
        const BondType& bond = bondType[bonds[i].type];

        if (UseConstraints && bond.canConstrain) {
            system.addConstraint(atom[0], atom[1],
                    bond.nominalLengthInAngstroms * OpenMM::NmPerAngstrom);
        }

Otherwise, we need to give the HarmonicBondForce the bond stretch parameters.

\ **Warning**\ *:* The constant used to specify the stiffness may be defined
differently between the existing MD code and OpenMM.  For instance, AMBER uses
the constant, as given in the harmonic *energy* term kx\ :sup:`2`\ , where
the force is 2kx (k = constant and x = distance).  OpenMM wants the constant, as
used in the *force* term kx (with energy 0.5 * kx\ :sup:`2`\ ).  So a factor
of 2 must be introduced when setting the bond stretch parameters in an OpenMM
system using data from an AMBER system.

.. code-block:: c

    bondStretch.addBond(atom[0], atom[1], bond.nominalLengthInAngstroms * OpenMM::NmPerAngstrom,
                        bond.stiffnessInKcalPerAngstrom2 * 2 * OpenMM::KJPerKcal *
                        OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);


**Non-bond exclusions:** Next, we deal with non-bond exclusions. These are
used for pairs of atoms that appear close to one another in the network of bonds
in a molecule. For atoms that close, normal non-bonded forces do not apply or
are reduced in magnitude.  First, we create a list of bonds to generate the non-
bond exclusions:

.. code-block:: c

    bondPairs.push_back(std::make_pair(atom[0], atom[1]));

OpenMM’s non-bonded force provides a convenient routine for creating the common
exceptions. These are: (1) for atoms connected by one bond (1-2) or connected by
just one additional bond (1-3), Coulomb and van der Waals terms do not apply;
and (2) for atoms connected by three bonds (1-4), Coulomb and van der Waals
terms apply but are reduced by a force-field dependent scale factor.  In
general, you may introduce additional exceptions, but the standard ones suffice
here and in many other circumstances.

.. code-block:: c

    // Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    // Create the 1-2-3 bond angle harmonic terms.
    for (int i=0; angles[i].type != EndOfList; ++i) {
         const int*       atom  = angles[i].atoms;
         const AngleType& angle = angleType[angles[i].type];

    // See note under bond stretch above regarding the factor of 2 here.
    bondBend.addAngle(atom[0],atom[1],atom[2],
    angle.nominalAngleInDegrees     * OpenMM::RadiansPerDegree,
    angle.stiffnessInKcalPerRadian2 * 2 *
    OpenMM::KJPerKcal);
    }

    // Create the 1-2-3-4 bond torsion (dihedral) terms.
    for (int i=0; torsions[i].type != EndOfList; ++i) {
         const int*         atom = torsions[i].atoms;
        const TorsionType& torsion = torsionType[torsions[i].type];
        bondTorsion.addTorsion(atom[0],atom[1],atom[2],atom[3],
                torsion.periodicity,
                torsion.phaseInDegrees  * OpenMM::RadiansPerDegree,
                torsion.amplitudeInKcal * OpenMM::KJPerKcal);
    }

The rest of the code is similar to the HelloSodiumChloride example and will not
be covered in detail here.  Please refer to the program HelloEthane.cpp itself,
which is well-commented, for additional details.
.. _examples-of-openmm-integration:

Examples of OpenMM Integration
###############################


GROMACS
*******

GROMACS is a large, complex application written primarily in C.  The
considerations involved in adapting it to use OpenMM are likely to be similar to
those faced by developers of other existing applications.

The first principle we followed in adapting GROMACS was to keep all OpenMM-related
code isolated to just a few files, while modifying as little of the
existing GROMACS code as possible.  This minimized the risk of breaking existing
parts of the code, while making the OpenMM-related parts as easy to work with as
possible.  It also minimized the need for C code to invoke the C++ API.  (This
would not be an issue if we used the OpenMM C API wrapper, but that is less
convenient than the C++ API, and placing all of the OpenMM calls into separate
C++ files solves the problem equally well.)  Nearly all of the OpenMM-specific
code is contained in a single file, openmm_wrapper.cpp.  It defines four
functions which encapsulate all of the interaction between OpenMM and the rest
of GROMACS:

\ :code:`openmm_init()`\ : As arguments, this function takes pointers to lots of
internal GROMACS data structures that describe the simulation to be run.  It
creates a System, Integrator, and Context based on them, then returns an opaque
reference to an object containing them.  That reference is an input argument to
all of the other functions defined in openmm_wrapper.cpp.  This allows
information to be passed between those functions without exposing it to the rest
of GROMACS.

\ :code:`openmm_take_one_step()`\ : This calls :code:`step(1)` on the
Integrator that was created by :code:`openmm_init()`\ .

\ :code:`openmm_copy_state()`\ : This calls :code:`getState()` on the
Context that was created by :code:`openmm_init()`\ , and then copies
information from the resulting State into various GROMACS data structures.  This
function is how state data generated by OpenMM is passed back to GROMACS for
output, analysis, etc.

\ :code:`openmm_cleanup()`\ : This is called at the end of the simulation.  It
deletes all the objects that were created by :code:`openmm_init()`\ .

This set of functions defines the interactions between GROMACS and OpenMM:
copying information from the application to OpenMM, performing integration,
copying information from OpenMM back to the application, and freeing resources
at the end of the simulation.  While the details of their implementations are
specific to GROMACS, this overall pattern is fairly generic.  A similar set of
functions can be used for many other applications as well.

TINKER-OpenMM
*************

TINKER is written primarily in Fortran, and uses common blocks extensively to
store application-wide parameters.  Rather than modify the TINKER build scripts
to allow C++ code, it was decided to use the OpenMM C API instead.  Despite
these differences, the overall approach used to add OpenMM support was very
similar to that used for GROMACS.

TINKER-OpenMM allows OpenMM to be used to calculate forces and energies and to
perform the integration in the main molecular dynamics loop. The only changes to
the TINKER source code are in the file :code:`dynamic.f` for the setup and
running of a simulation.  An added file, :code:`dynamic_openmm.c`\ , contains
the interface C code between TINKER and OpenMM.

The flow of the molecular dynamics simulation using OpenMM is as follows:

#. The TINKER code is used to read the AMOEBA parameter file,  the
   :code:`*.xyz` and :code:`*.key` files.  It then parses the command-line
   options.

#. The routine :code:`map_common_blocks_to_c_data_structs()` is
   called to map the FORTRAN common blocks to C data structures used in setting the
   parameters used by OpenMM.

#. The routine :code:`openmm_validate()` is called from
   :code:`dynamic.f` before the main loop.  This routine checks that all required
   options and settings obtained from the input in step (1) and common blocks in
   step (2) are available.  If an option or setting is unsupported, the program
   exits with an appropriate message.  The routine :code:`openmm_validate()`
   and the other OpenMM interface methods are in the file
   :code:`dynamic_openmm.c`\ .

#. :code:`openmm_init()` is called to create the OpenMM System,
   Integrator and Context objects..

#. :code:`openmm_take_steps()` is called to take a specified number
   of time steps.

#. :code:`openmm_update()` is then called to retrieve the state
   (energies/positions/velocities) and populate the appropriate TINKER data
   structures.  These values are converted from the OpenMM units of kJ/nm to kcal/Å
   when populating the TINKER arrays.

#. Once the main loop has completed, the routine
   :code:`openmm_cleanup()` is called to delete the OpenMM objects and release
   resources being used on the GPU.


.. _ring-polymer-molecular-dynamics-plugin:

Ring Polymer Molecular Dynamics (RPMD) Plugin
#############################################

Ring Polymer Molecular Dynamics (RPMD) provides an efficient approach to include
nuclear quantum effects in molecular simulations.\ :cite:`Craig2004`  When
used to calculate static equilibrium properties, RPMD reduces to path integral
molecular dynamics and gives an exact description of the effect of quantum
fluctuations for a given potential energy model.\ :cite:`Parrinello1984`  For
dynamical properties RPMD is no longer exact but has shown to be a good
approximation in many cases.

For a system with a classical potential energy *E*\ (\ *q*\ ), the RPMD
Hamiltonian is given by


.. math::
   H=\sum _{k=1}^{n}\left(\frac{{p}_{{k}^{2}}}{2m}+E({q}_{k})+\frac{m({k}_{B}Tn)^{2}}{2\hbar^{2}}({q}_{k}-{q}_{k-1})^{2}\right)


This Hamiltonian resembles that of a system of classical ring polymers where
different copies of the system are connected by harmonic springs.  Hence each
copy of the classical system is commonly referred to as a “bead”.  The spread of
the ring polymer representing each particle is directly related to its De
Broglie thermal wavelength (uncertainty in its position).

RPMD calculations must be converged with respect to the number *n* of beads
used.  Each bead is evolved at the effective temperature *nT*\ , where *T*
is the temperature for which properties are required.  The number of beads
needed to converge a calculation can be estimated using\ :cite:`Markland2008`


.. math::
   n>\frac{\hbar\omega_{max}}{{k}_{B}T}


where :math:`\omega_{max}` is the highest frequency in the problem.  For example, for
flexible liquid water the highest frequency is the OH stretch at around 3000
cm\ :sup:`-1`\ , so around 24 to 32 beads are needed depending on the accuracy
required.  For rigid water where the highest frequency is only around 1000
cm\ :sup:`-1`\ , only 6 beads are typically needed.  Due to the replication needed
of the classical system, the extra cost of the calculation compared to a
classical simulation increases linearly with the number of beads used.

This cost can be reduced by “contracting” the ring polymer to a smaller number
of beads.\ :cite:`Markland2008`  The rapidly changing forces are then computed
for the full number of beads, while slower changing forces are computed on a
smaller set.  In the case of flexible water, for example, a common arrangement
would be to compute the high frequency bonded forces on all 32 beads, the direct
space nonbonded forces on only 6 beads, and the reciprocal space nonbonded
forces on only a single bead.

Due to the stiff spring terms between the beads, NVE RPMD trajectories can
suffer from ergodicity problems and hence thermostatting is highly recommended,
especially when dynamical properties are not required.\ :cite:`Hall1984`  The
thermostat implemented here is the path integral Langevin equation (PILE)
approach.\ :cite:`Ceriotti2010`  This method couples an optimal white noise
Langevin thermostat to the normal modes of each polymer, leaving only one
parameter to be chosen by the user which controls the friction applied to the
center of mass of each ring polymer.  A good choice for this is to use a value
similar to that used in a classical calculation of the same system.

.. _testing-and-validation-of-openmm:

Testing and Validation of OpenMM
################################

The goal of testing and validation is to make sure that OpenMM works correctly.
That means that it runs without crashing or otherwise failing, and that it
produces correct results.  Furthermore, it must work correctly on a variety of
hardware platforms (e.g. different models of GPU), software platforms (e.g.
operating systems and OpenCL implementations), and types of simulations.

Three types of tests are used to validate OpenMM:

* **Unit tests:** These are small tests designed to test specific features
  or pieces of code in isolation.  For example, a test of HarmonicBondForce might
  create a System with just a few particles and bonds, compute the forces and
  energy, and compare them to the analytically expected values.  There are
  thousands of unit tests that collectively cover all of OpenMM.

* **System tests:** Whereas unit tests validate small features in
  isolation, system tests are designed to validate the entire library as a whole.
  They simulate realistic models of biomolecules and perform tests that are likely
  to fail if any problem exists anywhere in the library.

* **Direct comparison between OpenMM and other programs:**  The third type
  of validation performed is a direct comparison of the individual forces computed
  by OpenMM to those computed by other programs for a collection of biomolecules.


Each type of test is outlined in greater detail below; a discussion of the
current status of the tests is then given.


Description of Tests
********************


Unit tests
===========

The unit tests are with the source code, so if you build from source you can run
them yourself.  See Section :numref:`test-your-build` for details.  When you run the tests
(for example, by typing “make test” on Linux or Mac), it should produce output
something like this:
::

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

Each line represents a test suite, which may contain multiple unit tests.  If
all tests within a suite passed, it prints the word “Passed” and how long the
suite took to execute.  Otherwise it prints an error message.  If any tests
failed, you can then run them individually (each one is a separate executable)
to get more details on what went wrong.

System tests
============

Several different types of system tests are performed.  Each type is run for a
variety of systems, including both proteins and nucleic acids, and involving
both implicit and explicit solvent.  The full suite of tests is repeated for
both the CUDA and OpenCL platforms, using both single and double precision (and
for the integration tests, mixed precision as well), on a variety of operating
systems and hardware.  There are four types of tests:

* **Consistency between platforms:** The forces and energy are computed
  using the platform being tested, then compared to ones computed with the
  Reference platform.  The results are required to agree to within a small
  tolerance.
* **Energy-force consistency:** This verifies that the force really is the
  gradient of the energy.   It first computes the vector of forces for a given
  conformation.  It then generates four other conformations by displacing the
  particle positions by small amounts along the force direction.  It computes the
  energy of each one, uses those to calculate a fourth order finite difference
  approximation to the derivative along that direction, and compares it to the
  actual forces.  They are required to agree to within a small tolerance.
* **Energy conservation:** The system is simulated at constant energy using
  a Verlet integrator, and the total energy is periodically recorded.  A linear
  regression is used to estimate the rate of energy drift.  In addition, all
  constrained distances are monitored during the simulation to make sure they
  never differ from the expected values by more than the constraint tolerance.
* **Thermostability:** The system is simulated at constant temperature
  using a Langevin integrator.  The mean kinetic energy over the course of the
  simulation is computed and compared to the expected value based on the
  temperature.  In addition, all constrained distances are monitored during the
  simulation to make sure they never differ from the expected values by more than
  the constraint tolerance.


If you want to run the system tests yourself, they can be found in the
Subversion repository at https://simtk.org/svn/pyopenmm/trunk/test/system-tests.
Check out that directory, then execute the runAllTests.sh shell script.  It will
create a series of files with detailed information about the results of the
tests.  Be aware that running the full test suite may take a long time (possibly
several days) depending on the speed of your GPU.

Direct comparisons between OpenMM and other programs
====================================================

As a final check, identical systems are set up in OpenMM and in another program
(Gromacs 4.5 or Tinker 6.1), each one is used to compute the forces on atoms,
and the results are directly compared to each other.

Test Results
************

In this section, we highlight the major results obtained from the tests
described above.  They are not exhaustive, but should give a reasonable idea of
the level of accuracy you can expect from OpenMM.

Comparison to Reference Platform
================================

The differences between forces computed with the Reference platform and those
computed with the OpenCL or CUDA platform are shown in
:autonumref:`Table,force comparison between platforms`\ .  For every
atom, the relative difference between platforms was computed as
2·\|F\ :sub:`ref`\ –F\ :sub:`test`\ \|/(\|F\ :sub:`ref`\ \|+|F\ :sub:`test`\ \|), where
F\ :sub:`ref` is the force computed by the Reference platform and F\ :sub:`test`
is the force computed by the platform being tested (OpenCL or CUDA).  The median
over all atoms in a given system was computed to estimate the typical force
errors for that system.  Finally, the median of those values for all test
systems was computed to give the value shown in the table.

====================================  ========================  ====================  ===================  =====================
Force                                 OpenCL (single)           OpenCL (double)       CUDA (single)        CUDA (double)
====================================  ========================  ====================  ===================  =====================
Total Force                           2.53·10\ :sup:`-6`        1.44·10\ :sup:`-7`    2.56·10\ :sup:`-6`   8.78·10\ :sup:`-8`
HarmonicBondForce                     2.88·10\ :sup:`-6`        1.57·10\ :sup:`-13`   2.88·10\ :sup:`-6`   1.57·10\ :sup:`-13`
HarmonicAngleForce                    2.25·10\ :sup:`-5`        4.21·10\ :sup:`-7`    2.27·10\ :sup:`-5`   4.21·10\ :sup:`-7`
PeriodicTorsionForce                  8.23·10\ :sup:`-7`        2.44·10\ :sup:`-7`    9.27·10\ :sup:`-7`   2.56·10\ :sup:`-7`
RBTorsionForce                        4.86·10\ :sup:`-6`        1.46·10\ :sup:`-7`    4.72·10\ :sup:`-6`   1.4·10\ :sup:`-8`
NonbondedForce (no cutoff)            1.49·10\ :sup:`-6`        6.49·10\ :sup:`-8`    1.49·10\ :sup:`-6`   6.49·10\ :sup:`-8`
NonbondedForce (cutoff, nonperiodic)  9.74·10\ :sup:`-7`        4.88·10\ :sup:`-9`    9.73·10\ :sup:`-7`   4.88·10\ :sup:`-9`
NonbondedForce (cutoff, periodic)     9.82·10\ :sup:`-7`        4.88·10\ :sup:`-9`    9.8·10\ :sup:`-7`    4.88·10\ :sup:`-9`
NonbondedForce (Ewald)                1.33·10\ :sup:`-6`        5.22·10\ :sup:`-9`    1.33·10\ :sup:`-6`   5.22·10\ :sup:`-9`
NonbondedForce (PME)                  3.99·10\ :sup:`-5`        4.08·10\ :sup:`-6`    3.99·10\ :sup:`-5`   4.08·10\ :sup:`-6`
GBSAOBCForce (no cutoff)              3.0·10\ :sup:`-6`         1.76·10\ :sup:`-7`    3.09·10\ :sup:`-6`   9.4·10\ :sup:`-8`
GBSAOBCForce (cutoff, nonperiodic)    2.77·10\ :sup:`-6`        1.76·10\ :sup:`-7`    2.95·10\ :sup:`-6`   9.33·10\ :sup:`-8`
GBSAOBCForce (cutoff, periodic)       2.61·10\ :sup:`-6`        1.78·10\ :sup:`-7`    2.77·10\ :sup:`-6`   9.24·10\ :sup:`-8`
====================================  ========================  ====================  ===================  =====================

:autonumber:`Table,force comparison between platforms`\ :  Median relative difference in forces between Reference platform and
OpenCL/CUDA platform


Energy Conservation
===================

:autonumref:`Figure,energy drift` shows the total system energy versus time for three simulations of
ubiquitin in OBC implicit solvent.  All three simulations used the CUDA
platform, a Verlet integrator, a time step of 0.5 fs, no constraints, and no
cutoff on the nonbonded interactions.  They differ only in the level of numeric
precision that was used for calculations (see Chapter :numref:`platform-specific-properties`\ ).


.. figure:: ../../images/EnergyDrift.png
   :align: center

   :autonumber:`Figure,energy drift`: Total energy versus time for simulations run in three different
   precision modes.

For the mixed and double precision simulations, the drift in energy is almost
entirely diffusive with negligible systematic drift.  The single precision
simulation has a more significant upward drift with time, though the rate of
drift is still small compared to the rate of short term fluctuations.  Fitting a
straight line to each curve gives a long term rate of energy drift of 3.98
kJ/mole/ns for single precision, 0.217 kJ/mole/ns for mixed precision, and
0.00100 kJ/mole/ns for double precision.  In the more commonly reported units of
kT/ns/dof, these correspond to 4.3·10\ :sup:`-4` for single precision,
2.3·10\ :sup:`-5` for mixed precision, and 1.1·10\ :sup:`-7` for double precision.

Be aware that different simulation parameters will give different results.
These simulations were designed to minimize all sources of error except those
inherent in OpenMM.  There are other sources of error that may be significant in
other situations.  In particular:

* Using a larger time step increases the integration error (roughly
  proportional to *dt*\ :sup:`2`\ ).
* If a system involves constraints, the level of error will depend strongly on
  the constraint tolerance specified by the Integrator.
* When using Ewald summation or Particle Mesh Ewald, the accuracy will depend
  strongly on the Ewald error tolerance.
* Applying a distance cutoff to implicit solvent calculations will increase the
  error, and the shorter the cutoff is, the greater the error will be.


As a result, the rate of energy drift may be much greater in some simulations
than in the ones shown above.

Comparison to Gromacs
=====================

OpenMM and Gromacs 4.5.5 were each used to compute the atomic forces for
dihydrofolate reductase (DHFR) in implicit and explicit solvent.  The implicit
solvent calculations used the OBC solvent model and no cutoff on nonbonded
interactions.  The explicit solvent calculations used Particle Mesh Ewald and a
1 nm cutoff on direct space interactions.  For OpenMM, the Ewald error tolerance
was set to 10\ :sup:`-6`\ .  For Gromacs, :code:`fourierspacing` was set to
0.07 and :code:`ewald_rtol` to 10\ :sup:`-6`\ .  No constraints were applied
to any degrees of freedom.  Both programs used single precision.  The test was
repeated for OpenCL, CUDA, and CPU platforms.

For every atom, the relative difference between OpenMM and Gromacs was computed
as 2·\|F\ :sub:`MM`\ –F\ :sub:`Gro`\ \|/(\|F\ :sub:`MM`\ \|+\|F\ :sub:`Gro`\ \|),
where F\ :sub:`MM` is the force computed by OpenMM and F\ :sub:`Gro` is the
force computed by Gromacs.  The median over all atoms is shown in :autonumref:`Table,comparison to Gromacs`\ .

=============   ===================  ===================  ===================
Solvent Model   OpenCL               CUDA                 CPU
=============   ===================  ===================  ===================
Implicit        7.66·10\ :sup:`-6`   7.68·10\ :sup:`-6`   1.94·10\ :sup:`-5`
Explicit        6.77·10\ :sup:`-5`   6.78·10\ :sup:`-5`   9.89·10\ :sup:`-5`
=============   ===================  ===================  ===================

:autonumber:`Table,comparison to Gromacs`\ :  Median relative difference in forces between OpenMM and Gromacs

.. _drude-plugin:

Drude Plugin
############

Drude oscillators are a method for incorporating electronic polarizability into
a model.\ :cite:`Lamoureux2003`  For each polarizable particle, a second
particle (the “Drude particle”) is attached to it by an anisotropic harmonic
spring.  When both particles are at the same location, they are equivalent to an
ordinary point particle.  Applying an electric field causes the Drude particle
to move a short distance away from its parent particle, creating an induced
dipole moment.  The polarizability :math:`\alpha` is related to the charge *q* on
the Drude particle and the spring constant *k* by



.. math::
   \alpha =\frac{{q}^{2}}{k}


A damped interaction\ :cite:`Thole1981` is used between dipoles that are
bonded to each other.

The equations of motion can be integrated with three different methods:

#. In the Self Consistent Field (SCF) method, the ordinary particles are first
   updated as usual.  A local energy minimization is then performed to select new
   positions for the Drude particles.  This ensures that the induced dipole moments
   respond instantly to changes in their environments.  This method is accurate but
   computationally expensive.
#. In the extended Lagrangian method, the positions of the Drude particles are
   treated as dynamical variables, just like any other particles.  A small amount
   of mass is transferred from the parent particles to the Drude particles,
   allowing them to be integrated normally.  A dual Langevin or Nose-Hoover integrator is used to
   maintain the center of mass of each Drude particle pair at the system
   temperature, while using a much lower temperature for their relative internal
   motion.  In practice, this produces dipole moments very close to those from the
   SCF solution while being much faster to compute.
#. The Nosé-Hoover dual thermostat method.  In this approach the motion of
   non-Drude sites and center of mass motion of Drude pairs are thermostated to
   the target temperature with one thermostat.  Another thermostat is used to keep
   relative motion of Drude pairs to a different, typically much lower,
   temperature to maintain separation of nuclear and electronic degrees of
   freedom.  The minimal specification is as follows::

      DrudeNoseHooverIntegrator integrator(temperature, frequency,
                                           temperatureDrude, frequencyDrude,
                                           1*femtoseconds)

   Where the first and third arguments specify the center-of-mass temperature and
   relative temperature for each Drude pair, respecitvely.  The second and fourth
   arguments describe the frequency of interaction with the center-of-mass and
   relative heat baths, respectively, and should be specified with inverse time
   units.  The fifth argument is the timestep.  The multi-timestep and Nosé-Hoover
   chain length may also be specified, but sensible defaults are provided.
.. _compiling-openmm-from-source-code:

Compiling OpenMM from Source Code
#################################

This chapter describes the procedure for building and installing OpenMM from
source code.  In most cases, it is best to use the pre-built versions installed
with conda.  Sometimes you might need to build from source, though, such as if
you want to modify OpenMM, or if conda does not provide packages compatible with
your environment.

We first describe how to build on Linux or Mac.  We then describe how to build
on Windows, where the process is slightly different.

.. _compiling-openmm-from-source-linux:

Compiling on Linux and Mac
**************************

Prerequisites
=============

Before building OpenMM from source, you will need certain tools.

C++ compiler
------------

You must have a C++ compiler installed before attempting to build OpenMM from
source.  All recent versions of clang or gcc should work correctly.  On Linux,
you can install the compiler with your system's standard package manager (such
as apt or yum).  On MacOS, you can get a C++ compiler by installing the Xcode
developer tools from the App Store.  Alternatively you can use a package manager
such as Homebrew to install clang or gcc.

Python
------

You will need a 64-bit Python 3.x environment.  We recommend using Miniconda
(https://docs.conda.io/en/latest/miniconda.html), which includes the conda
package manager.

OpenMM Source Code
------------------

You will also need the OpenMM source code.  To download it:

#. Browse to https://github.com/openmm/openmm/releases.
#. Find the latest release and click the link to download the source as either
   a .zip or .tar.gz file.
#. Unpack the file.  Note the location where you unpacked the OpenMM source code.

Alternatively, if you want the most recent development version of the code rather
than the version corresponding to a particular release, you can get it from
https://github.com/openmm/openmm.  Be aware that the development code is constantly
changing, may contain bugs, and should never be used for production work.  If
you want a stable, well tested version of OpenMM, you should download the source
code for the latest release as described above.

CUDA or OpenCL Support
----------------------

If you want to compile OpenMM with support for running on GPUs, you will need
CUDA and/or OpenCL.  MacOS comes with OpenCL built in, so nothing else needs to
be installed.  For Linux, you need an appropriate SDK.

The easiest way is to install the most recent CUDA Toolkit from https://developer.nvidia.com/cuda-downloads.
It includes the headers and libraries needed to compile both CUDA and OpenCL
applications.  In addition, it has runtime libraries that are needed for running
CUDA applications.  The runtime components for OpenCL applications are included
with the GPU drivers from NVIDIA, AMD, and Intel, so make sure you have an
up-to-date driver.

Other Required Software
-----------------------

Several other tools are required to build OpenMM.  The easiest way to install
them is with conda.  The following command will install everything needed to
build OpenMM.
::

    conda install -c conda-forge cmake make cython swig fftw doxygen numpy

Step 1: Configure with CMake
============================

First, create a directory in which to build OpenMM.  Starting from the root
level of the OpenMM source tree (the directory containing the top CMakeLists.txt
file), execute the following commands:

.. code-block:: none

    mkdir build
    cd build
    ccmake ..

That is not a typo.  :code:`ccmake` has two c’s.  CCMake is the visual CMake
configuration tool.  Press :code:`c` within the CCMake interface to load the
OpenMM build scripts and begin configuring CMake.

There are several variables that can be adjusted in the CMake interface:

* Set the variable CMAKE_INSTALL_PREFIX to the location where you want to
  install OpenMM. If you are using conda environments this variable should point to
  the full path of the root directory of your environment.
* Set the variable PYTHON_EXECUTABLE to the Python interpreter you plan to use
  OpenMM with.  Usually this will be detected automatically.
* There are lots of options starting with OPENMM_BUILD that control
  whether to build particular features of OpenMM, such as plugins, API wrappers,
  and documentation.
* Usually the OpenCL library and headers will be detected automatically.  If for
  any reason CMake is unable to find them, set OPENCL_INCLUDE_DIR to point to
  the directory containing the headers (usually /usr/local/cuda/include on Linux)
  and OPENCL_LIBRARY to point to the library (usually /usr/local/cuda/lib64/libOpenCL.so
  on Linux).  

Configure (press “c”) again.  Adjust any variables that cause an error.

Continue to configure (press “c”) until no starred CMake variables are
displayed, then press “g” to generate the makefiles for building the project.

Step 2: Build
=============

Build OpenMM with the command::

    make

.. _test-your-build:

Step 3: Test your build
=======================

This step is optional but recommended. Tests can take up to several minutes depending on your
hardware configuration.

It is recommended that you make sure your local build of OpenMM works before trying
to install.

After OpenMM has been built, you should run the unit tests to make sure it
works.  Enter the command::

    make test

You should see a series of test results like this:

.. code-block:: none

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

:code:`Passed` is good.  :code:`FAILED` is bad.  If any tests fail, you
can run them individually to get more detailed error information.  For example,
::

    ./TestReferenceLangevinIntegrator

Note that some tests are stochastic, and therefore are expected to fail a small
fraction of the time.  These tests will say so in the error message:

.. code-block:: none

    exception: Assertion failure at TestReferenceLangevinIntegrator.cpp:129.  Expected 9.97741,
        found 10.7884 (This test is stochastic and may occasionally fail)

If you get an error message such as :code:`exception: Error launching CUDA compiler: 32512` you need
to specify the path to the CUDA compiler (nvcc) using the :code:`OPENMM_CUDA_COMPILER` environment
variable, for example using something like the following::

    OPENMM_CUDA_COMPILER=/<path_to_custom_cuda_dir>/nvcc

Step 3: Install
===============
Install your local build of OpenMM using the following command::

    make install

If you are installing to a system directory, such as /usr/local/openmm/, you will
need admin capabilities to install, in this case use::

    sudo make install

Step 3: Install the Python API
==============================

Build and install the Python API with the command::

    make PythonInstall

If you are installing into the system Python, such as /usr/bin/python, you will
need to type::

    sudo make PythonInstall

You can test the Python API installation using::

    python -m openmm.testInstallation

Congratulations! You have successfully built and installed OpenMM from source.

Compiling on Windows
********************

Prerequisites
=============

Before building OpenMM from source, you will need certain tools.

C++ compiler
------------

On Windows systems, use the C++ compiler in Visual Studio 2015 or later.  You
can download a free version of Visual Studio from https://visualstudio.microsoft.com.

Python
------

You will need a 64-bit Python 3.x environment.  We recommend using Miniconda
(https://docs.conda.io/en/latest/miniconda.html), which includes the conda
package manager.

CMake
-----

CMake is the build system used for OpenMM.  You must install CMake version 3.17
or higher before attempting to build OpenMM from source.  You can get CMake from
http://www.cmake.org/.

OpenMM Source Code
------------------

You will also need the OpenMM source code.  To download it:

#. Browse to https://github.com/openmm/openmm/releases.
#. Find the latest release and click the link to download the source as either
   a .zip or .tar.gz file.
#. Unpack the file.  Note the location where you unpacked the OpenMM source code.

Alternatively, if you want the most recent development version of the code rather
than the version corresponding to a particular release, you can get it from
https://github.com/openmm/openmm.  Be aware that the development code is constantly
changing, may contain bugs, and should never be used for production work.  If
you want a stable, well tested version of OpenMM, you should download the source
code for the latest release as described above.

CUDA or OpenCL Support
----------------------

If you want to compile OpenMM with support for running on GPUs, you will need
CUDA and/or OpenCL.  Install the most recent CUDA Toolkit from https://developer.nvidia.com/cuda-downloads.
It includes the headers and libraries needed to compile both CUDA and OpenCL
applications.  In addition, it has runtime libraries that are needed for running
CUDA applications.  The runtime components for OpenCL applications are included
with the GPU drivers from NVIDIA, AMD, and Intel, so make sure you have an
up-to-date driver.

Other Required Software
-----------------------

Several other tools are required to build OpenMM.  The easiest way to install
them is with conda.  From the Windows Start menu, select "Anaconda Prompt (Miniconda3)".
It will open a command window that is preconfigured for conda.  Enter the
following command to install everything needed to build OpenMM.
::

    conda install -c conda-forge cython swig fftw doxygen numpy

Step 1: Configure with CMake
============================

First, create a directory in which to build OpenMM.  In the "Anaconda Prompt"
window opened above, cd to the root level of the OpenMM source tree (the
directory containing the top CMakeLists.txt file).  Execute the following commands:

.. code-block:: none

    mkdir build
    cd build
    "C:\Program Files\CMake\bin\cmake-gui.exe"

This will launch the CMake GUI configuration tool.  It is critical that you
launch it from the "Anaconda Prompt" window as shown above.  Do *not* launch
it from the Start menu.  If you do, it will not be able to find the tools you
installed with conda.

#. In the box labeled "Where is the source code" browse to the OpenMM source directory
   (containing the top CMakeLists.txt file).
#. In the box labeled "Where to build the binaries" browse to the build
   directory you just created.
#. Click the "Configure" button at the bottom of the CMake window.
#. Select "Visual Studio 16 2019" from the  list of Generators (or whichever
   version you have installed) and click "Finish".

There are several variables that can be adjusted in the CMake interface:

* Set the variable CMAKE_INSTALL_PREFIX to the location where you want to
  install OpenMM.
* Set the variable PYTHON_EXECUTABLE to the Python interpreter you plan to use
  OpenMM with.  Usually this will be detected automatically.
* There are lots of options starting with OPENMM_BUILD that control
  whether to build particular features of OpenMM, such as plugins, API wrappers,
  and documentation.
* Usually the OpenCL library and headers will be detected automatically.  If for
  any reason CMake is unable to find them, set OPENCL_INCLUDE_DIR to point to
  the directory containing the headers (usually
  "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/include", except
  with the correct version number for the toolkit you installed) and
  OPENCL_LIBRARY to point to the library (usually "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/lib/x64/OpenCL.lib").  

Press "Configure" again.  Adjust any variables that cause an error.

Continue to press "Configure" until no red CMake variables are displayed, then
press "Generate" to create the Visual Studio project files for building OpenMM.

Step 2: Build and Install
=========================

#. Open the file :file:`OpenMM.sln` in your build directory in Visual Studio.
   Note that this file will appear as just :file:`OpenMM` if you have configured
   Explorer to hide file name extensions.
#. Set the configuration type to "Release" (not "Debug") in the toolbar.
#. From the Build menu, select "Build Solution".  This takes some time.
#. In the Solution Explorer, right-click on "INSTALL" and select "Build".

Step 3: Install the Python API
==============================

In the Solution Explorer, right-click on "PythonInstall" and select "Build".

Step 4: Test your build
=======================

After OpenMM has been built, you should run the unit tests to make sure it
works.  In the Solution Explorer, right-click on "RUN_TESTS" and select "Build".
You should see a series of test results like this:

.. code-block:: none

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

:code:`Passed` is good.  :code:`FAILED` is bad.  If any tests fail, you
can run them individually to get more detailed error information.  Right-click
on a test in the Solution Explorer and select "Debug > Start New Instance".

Note that some tests are stochastic, and therefore are expected to fail a small
fraction of the time.  These tests will say so in the error message:

.. code-block:: none

    exception: Assertion failure at TestReferenceLangevinIntegrator.cpp:129.  Expected 9.97741,
        found 10.7884 (This test is stochastic and may occasionally fail)

Congratulations! You have successfully built and installed OpenMM from source.


Building the Documentation (Optional)
*************************************

The documentation that you're currently reading, as well as the Developer Guide and API
documentation, can be built through CMake.  To do that, you need to install a few
additional tools.  The easiest way is to use :code:`conda` to install them into
your Python environment.  The following command installs everything needed to
build documentation in HTML format.
::

    conda install -c conda-forge sphinx sphinxcontrib-bibtex breathe jinja2

To build documentation in PDF format, you also must have a functional LaTeX
installation.  It can be obtained from https://www.latex-project.org/get.

If you want to build documentation, make sure that OPENMM_GENERATE_API_DOCS is
set to ON when configuring the build in CMake.

To build the documentation, use the following build targets.

* :code:`sphinxhtml`: Build the User Guide and Developer Guide in HTML format.

* :code:`sphinxpdf`: Build the User Guide and Developer Guide in PDF format.

* :code:`C++ApiDocs`: Build the C++ API documentation.

* :code:`PythonApiDocs`: Build the Python API documentation.  This target
  requires that you have already built the :code:`install` target, such as with
  :code:`make install`.

On Linux or Mac, build a target using the :code:`make` command.  For example,
::

    make sphinxhtml

On Windows, right-click on the target in the Solution Explorer and select "Build".

After building the documentation, build the :code:`install` target to install
the documentation into the installation directory (the one you specified with
CMAKE_INSTALL_PREFIX).

Using local build of OpenMM alongside conda tools that depend on it
*******************************************************************

A common case is to have a local build of OpenMM in the same environment as other tools
that depend on it. This can be achieved by forcing a remove of OpenMM when you install
your tools using conda.

We will use :code:`openmmtools` as an example here, but it can be replaced with any
other software package that requires OpenMM.

Step 1: Install your tools as usual
===================================

Install your tools using conda as you commonly do, for example using::

    conda install -c conda-forge  openmmtools

This will pull the conda-forge package of :code:`openmm` which we don't want since we want
to use our local build.

Step 2: Remove conda openmm package
===================================

To remove the openmm package that was installed in the previous step, we can use::

    conda remove --force openmm

This will remove the :code:`openmm` package without changing or removing dependencies.

Step 3: Install local build of openmm
=====================================

Now we just install our local build of :code:`openmm` as instructed in
:ref:`_compiling-openmm-from-source-linux`.. _using-openmm-with-software-written-in-languages-other-than-c++:

Using OpenMM with Software Written in Languages Other than C++
##############################################################

Although the native OpenMM API is object-oriented C++ code, it is possible to
directly translate the interface so that it is callable from C, Fortran 95, and
Python with no substantial conceptual changes. We have developed a
straightforward mapping for these languages that, while perhaps not the most
elegant possible, has several advantages:

* Almost all documentation, training, forum discussions, and so on are equally
  useful to users of all these languages. There are syntactic differences of
  course, but all the important concepts remain unchanged.
* We are able to generate the C, Fortran, and Python APIs from the C++ API.
  Obviously, this reduces development effort, but more importantly it means that
  the APIs are likely to be error-free and are always available immediately when
  the native API is updated.
* Because OpenMM performs expensive operations “in bulk” there is no noticeable
  overhead in accessing these operations through the C, Fortran, or Python APIs.
* All symbols introduced to a C or Fortran program begin with the prefix
  “\ :code:`OpenMM_`\ ” so will not interfere with symbols already in use.


*Availability of APIs in other languages:*  All necessary C and Fortran
bindings are built in to the main OpenMM library; no separate library is
required.  The Python wrappers are contained in a module that is distributed
with OpenMM and that can be installed by executing its setup.py script in the
standard way.

(This doesn’t apply to most users: if you are building your own OpenMM from
source using CMake and want the API bindings generated, be sure to enable the
:code:`OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS` option for C and Fortran, or
:code:`OPENMM_BUILD_PYTHON_WRAPPERS` option for Python.  The Python module
will be placed in a subdirectory of your main build directory called “python”)

*Documentation for APIs in other languages:*  While there is extensive
Doxygen documentation available for the C++ and Python APIs, there is no
separate on-line documentation for the C and Fortran API. Instead, you should
use the C++ documentation, employing the mappings described here to figure out
the equivalent syntax in C or Fortran.

C API
*****

Before you start writing your own C program that calls OpenMM, be sure you can
build and run the two C examples that are supplied with OpenMM (see Chapter :numref:`openmm-tutorials`\ ).
These can be built from the supplied :code:`Makefile` on Linux and Mac, or
supplied :code:`NMakefile` and Visual Studio solution files on Windows.

The example programs are :code:`HelloArgonInC` and
:code:`HelloSodiumChlorideInC`\ . The argon example serves as a quick check that
your installation is set up properly and you know how to build a C program that
is linked with OpenMM. It will also tell you whether OpenMM is executing on the
GPU or is running (slowly) on the Reference platform. However, the argon example
is not a good template to follow for your own programs. The sodium chloride
example, though necessarily simplified, is structured roughly in the way we
recommended you set up your own programs to call OpenMM. Please be sure you have
both of these programs executing successfully on your machine before continuing.

Mechanics of using the C API
============================

The C API is generated automatically from the C++ API when OpenMM is built.
There are two resulting components: C bindings (functions to call), and C
declarations (in a header file). The C bindings are small :code:`extern`
(global) interface functions, one for every method of every OpenMM class, whose
signatures (name and arguments) are predictable from the class name and method
signatures. There are also “helper” types and functions provided for the few
cases in which the C++ behavior cannot be directly mapped into C. These
interface and helper functions are compiled in to the main OpenMM library so
there is nothing special you have to do to get access to them.

In the :code:`include` subdirectory of your OpenMM installation directory,
there is a machine-generated header file :code:`OpenMMCWrapper.h` that
should be #included in any C program that is to make calls to OpenMM functions.
That header contains declarations for all the OpenMM C interface functions and
related types. Note that if you follow our suggested structure, you will not
need to include this file in your :code:`main()` compilation unit but can
instead use it only in a local file that you write to provide a simple interface
to your existing code (see Chapter :numref:`openmm-tutorials`).

Mapping from the C++ API to the C API
=====================================

The automated generator of the C “wrappers” follows the translation strategy
shown in :autonumref:`Table,C API`\ . The idea is that if you see the construct on the left in
the C++ API documentation, you should interpret it as the corresponding
construct on the right in C. Please look at the supplied example programs to see
how this is done in practice.

==========================  =========================================  ===================================================
Construct                   C++ API declaration                        Equivalent in C API
==========================  =========================================  ===================================================
namespace                   OpenMM\::                                  OpenMM\_ (prefix)
class                       class OpenMM::ClassName                    typedef OpenMM_ClassName
constant                    OpenMM::RadiansPerDeg                      OpenMM_RadiansPerDeg (static constant)
class enum                  OpenMM::State::Positions                   OpenMM_State_Positions
constructor                 new OpenMM::ClassName()                    | OpenMM_ClassName* OpenMM_ClassName_create()
                                                                       | (additional constructors are _create_2(), etc.)
destructor                  | OpenMM::ClassName* thing;                | OpenMM_ClassName* thing;
                            | delete thing;                            | OpenMM_ClassName_destroy(thing);
class method                | OpenMM::ClassName* thing;                | OpenMM_ClassName* thing;
                            | thing->method(args);                     | OpenMM_ClassName_method(thing, args)
Boolean (type & constants)  | bool                                     | OpenMM_Boolean
                            | true, false                              | OpenMM_True(1), OpenMM_False(0)
string                      std::string                                char*
3-vector                    OpenMM::Vec3                               typedef OpenMM_Vec3
arrays                      | std::vector<std::string>                 | typedef OpenMM_StringArray
                            | std::vector<double>                      | typedef OpenMM_DoubleArray
                            | std::vector<Vec3>                        | typedef OpenMM_Vec3Array
                            | std::vector<std::pair<int,int>>          | typedef OpenMM_BondArray
                            | std::map<std::string,double>             | typedef OpenMM_ParameterArray
==========================  =========================================  ===================================================

:autonumber:`Table,C API`\ : Default mapping of objects from the C++ API to the C API
There are some exceptions to the generic translation rules shown in the table;
they are enumerated in the next section. And because there are no C++ API
equivalents to the array types, they are described in detail below.

Exceptions
==========

These two methods are handled somewhat differently in the C API than in the C++ API:

* **OpenMM::Context::getState()** The C version,
  :code:`OpenMM_Context_getState()`\ , returns a pointer to a heap allocated
  :code:`OpenMM_State` object. You must then explicitly destroy this
  :code:`State` object when you are done with it, by calling
  :code:`OpenMM_State_destroy()`\ .
* **OpenMM::Platform::loadPluginsFromDirectory()** The C version
  :code:`OpenMM_Platform_loadPluginsFromDirectory()` returns a heap-allocated
  :code:`OpenMM_StringArray` object containing a list of all the file names
  that were successfully loaded. You must then explicitly destroy this
  :code:`StringArray` object when you are done with it. Do not ignore the return
  value; if you do you’ll have a memory leak since the :code:`StringArray`
  will still be allocated.


(In the C++ API, the equivalent methods return references into existing memory
rather than new heap-allocated memory, so the returned objects do not need to be
destroyed.)

OpenMM_Vec3 helper type
=======================

Unlike the other OpenMM objects which are opaque and manipulated via pointers,
the C API provides an explicit definition for the C :code:`OpenMM_Vec3` type
that is compatible with the :code:`OpenMM::Vec3` type. The definition of
:code:`OpenMM_Vec3` is:

.. code-block:: c

    typedef struct {double x, y, z;} OpenMM_Vec3;

You can work directly with the individual fields of this type from your C
program if you want. For convenience, a scale() function is provided that
creates a new OpenMM_Vec3 from an old one and a scale factor:

.. code-block:: c

    OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale);

Array helper types
==================

C++ has built-in container types :code:`std::vector` and :code:`std::map`
which OpenMM uses to manipulate arrays of objects. These don’t have direct
equivalents in C, so we supply special array types for each kind of object for
which OpenMM creates containers. These are: string, double, Vec3, bond, and
parameter map. See :autonumref:`Table,C arrays` for the names of the C types for each of these
object arrays. Each of the array types provides these functions (prefixed by
:code:`OpenMM_` and the actual *Thing* name), with the syntax shown
conceptually since it differs slightly for each kind of object.

.. tabularcolumns:: |l|L|

=======================================================  =========================================================================================================================================================================================================
Function                                                 Operation
=======================================================  =========================================================================================================================================================================================================
*Thing*\ Array\* create(int size)                        Create a heap-allocated array of *Things*\ , with space pre-allocated to hold :code:`size` of them. You can start at :code:`size==0` if you want since these arrays are dynamically resizeable.
void destroy(\ *Thing*\ Array\*)                         Free the heap space that is currently in use for the passed-in array of *Things*\ .
int getSize(\ *Thing*\ Array\*)                          Return the current number of *Things* in this array. This means you can :code:`get()` and :code:`set()` elements up to :code:`getSize()-1`\ .
void resize(\ *Thing*\ Array\*, int size)                Change the size of this array to the indicated value which may be smaller or larger than the current size. Existing elements remain in their same locations as long as they still fit.
void append(\ *Thing*\ Array\*, *Thing*\ )               Add a *Thing* to the end of the array, increasing the array size by one. The precise syntax depends on the actual type of *Thing*\ ; see below.
void set(\ *Thing*\ Array\*, int index, *Thing*\ )       Store a copy of *Thing* in the indicated element of the array (indexed from 0). The array must be of length at least :code:`index+1`\ ; you can’t grow the array with this function.
*Thing* get(\ *Thing*\ Array\*, int index)               Retrieve a particular element from the array (indexed from 0). (For some Things the value is returned in arguments rather than as the function return.)
=======================================================  =========================================================================================================================================================================================================

:autonumber:`Table,C arrays`\ : Generic description of array helper types

Here are the exact declarations with deviations from the generic description
noted, for each of the array types.

OpenMM_DoubleArray
------------------

.. code-block:: c

    OpenMM_DoubleArray*
                OpenMM_DoubleArray_create(int size);
    void        OpenMM_DoubleArray_destroy(OpenMM_DoubleArray*);
    int         OpenMM_DoubleArray_getSize(const OpenMM_DoubleArray*);
    void        OpenMM_DoubleArray_resize(OpenMM_DoubleArray*, int size);
    void        OpenMM_DoubleArray_append(OpenMM_DoubleArray*, double value);
    void        OpenMM_DoubleArray_set(OpenMM_DoubleArray*, int index, double value);
    double      OpenMM_DoubleArray_get(const OpenMM_DoubleArray*, int index);

OpenMM_StringArray
------------------

.. code-block:: c

    OpenMM_StringArray*
                OpenMM_StringArray_create(int size);
    void        OpenMM_StringArray_destroy(OpenMM_StringArray*);
    int         OpenMM_StringArray_getSize(const OpenMM_StringArray*);
    void        OpenMM_StringArray_resize(OpenMM_StringArray*, int size);
    void        OpenMM_StringArray_append(OpenMM_StringArray*, const char* string);
    void        OpenMM_StringArray_set(OpenMM_StringArray*, int index, const char* string);
    const char* OpenMM_StringArray_get(const OpenMM_StringArray*, int index);

OpenMM_Vec3Array
----------------

.. code-block:: c

    OpenMM_Vec3Array*
                OpenMM_Vec3Array_create(int size);
    void        OpenMM_Vec3Array_destroy(OpenMM_Vec3Array*);
    int         OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array*);
    void        OpenMM_Vec3Array_resize(OpenMM_Vec3Array*, int size);
    void        OpenMM_Vec3Array_append(OpenMM_Vec3Array*, const OpenMM_Vec3 vec);
    void        OpenMM_Vec3Array_set(OpenMM_Vec3Array*, int index, const OpenMM_Vec3 vec);
    const OpenMM_Vec3*
                OpenMM_Vec3Array_get(const OpenMM_Vec3Array*, int index);

OpenMM_BondArray
----------------

Note that bonds are specified by pairs of integers (the atom indices). The
:code:`get()` method returns those in a pair of final arguments rather than as
its functional return.

.. code-block:: c

    OpenMM_BondArray*
                OpenMM_BondArray_create(int size);
    void        OpenMM_BondArray_destroy(OpenMM_BondArray*);
    int         OpenMM_BondArray_getSize(const OpenMM_BondArray*);
    void        OpenMM_BondArray_resize(OpenMM_BondArray*, int size);
    void        OpenMM_BondArray_append(OpenMM_BondArray*, int particle1, int particle2);
    void        OpenMM_BondArray_set(OpenMM_BondArray*, int index, int particle1, int particle2);
    void        OpenMM_BondArray_get(const OpenMM_BondArray*, int index,
                                     int* particle1, int* particle2);

OpenMM_ParameterArray
---------------------

OpenMM returns references to internal :code:`ParameterArrays` but does not
support user-created :code:`ParameterArrays`\ , so only the :code:`get()`
and :code:`getSize()` functions are available. Also, note that since this is
actually a map rather than an array, the “index” is the *name* of the
parameter rather than its ordinal.

.. code-block:: c

    int         OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray*);
    double      OpenMM_ParameterArray_get(const OpenMM_ParameterArray*, const char* name);


Fortran 95 API
*****************

Before you start writing your own Fortran program that calls OpenMM, be sure you
can build and run the two Fortran examples that are supplied with OpenMM (see
Chapter :numref:`openmm-tutorials`). These can be built from the supplied :code:`Makefile` on Linux
and Mac, or supplied :code:`NMakefile` and Visual Studio solution files on
Windows.

The example programs are :code:`HelloArgonInFortran` and
:code:`HelloSodiumChlorideInFortran`\ . The argon example serves as a quick
check that your installation is set up properly and you know how to build a
Fortran program that is linked with OpenMM. It will also tell you whether OpenMM
is executing on the GPU or is running (slowly) on the Reference platform.
However, the argon example is not a good template to follow for your own
programs. The sodium chloride example, though necessarily simplified, is
structured roughly in the way we recommended you set up your own programs to
call OpenMM. Please be sure you have both of these programs executing
successfully on your machine before continuing.

Mechanics of using the Fortran API
==================================

The Fortran API is generated automatically from the C++ API when OpenMM is
built. There are two resulting components: Fortran bindings (subroutines to
call), and Fortran declarations of types and subroutines (in the form of a
Fortran 95 module file). The Fortran bindings are small interface subroutines,
one for every method of every OpenMM class, whose signatures (name and
arguments) are predictable from the class name and method signatures. There are
also “helper” types and subroutines provided for the few cases in which the C++
behavior cannot be directly mapped into Fortran. These interface and helper
subroutines are compiled in to the main OpenMM library so there is nothing
special you have to do to get access to them.

Because Fortran is case-insensitive, calls to Fortran subroutines (however
capitalized) are mapped by the compiler into all-lowercase or all-uppercase
names, and different compilers use different conventions. The automatically-generated
OpenMM Fortran “wrapper” subroutines, which are generated in C and
thus case-sensitive, are provided in two forms for compatibility with the
majority of Fortran compilers, including Intel Fortran and gfortran. The two
forms are: (1) all-lowercase with a trailing underscore, and (2) all-uppercase
without a trailing underscore. So regardless of the Fortran compiler you are
using, it should find a suitable subroutine to call in the main OpenMM library.

In the :code:`include` subdirectory of your OpenMM installation directory,
there is a machine-generated module file :code:`OpenMMFortranModule.f90`
that must be compiled along with any Fortran program that is to make calls to
OpenMM functions. (You can look at the :code:`Makefile` or Visual Studio
solution file provided with the OpenMM examples to see how to build a program
that uses this module file.) This module file contains definitions for two
modules: :code:`MODULE OpenMM_Types` and :code:`MODULE OpenMM`\ ; however,
only the :code:`OpenMM` module will appear in user programs (it references
the other module internally). The modules contain declarations for all the
OpenMM Fortran interface subroutines, related types, and parameters (constants).
Note that if you follow our suggested structure, you will not need to
:code:`use` the :code:`OpenMM` module in your :code:`main()`
compilation unit but can instead use it only in a local file that you write to
provide a simple interface to your existing code (see Chapter :numref:`openmm-tutorials`).

Mapping from the C++ API to the Fortran API
===========================================

The automated generator of the Fortran “wrappers” follows the translation
strategy shown in :autonumref:`Table,Fortran API`\ . The idea is that if you see the construct on the
left in the C++ API documentation, you should interpret it as the corresponding
construct on the right in Fortran. Please look at the supplied example programs
to see how this is done in practice. Note that all subroutines and modules are
declared with “\ :code:`implicit none`\ ”, meaning that the type of every symbol
is declared explicitly and should not be inferred from the first letter of the
symbol name.

==========================  ===================================  ========================================================
Construct                   C++ API declaration                  Equivalent in Fortran API
==========================  ===================================  ========================================================
namespace                   OpenMM\::                            OpenMM\_ (prefix)
class                       class OpenMM::ClassName              type (OpenMM_ClassName)
constant                    OpenMM::RadiansPerDeg                parameter (OpenMM_RadiansPerDeg)
class enum                  OpenMM::State::Positions             parameter (OpenMM_State_Positions)
constructor                 new OpenMM::ClassName()              | type (OpenMM_ClassName) thing
                                                                 | call OpenMM_ClassName_create(thing)
                                                                 | (additional constructors are \_create_2(), etc.)
destructor                  | OpenMM::ClassName* thing;          | type (OpenMM_ClassName) thing
                            | delete thing;                      | call OpenMM_ClassName_destroy(thing)
class method                | OpenMM::ClassName* thing;          | type (OpenMM_ClassName) thing
                            | thing->method(args*)               | call OpenMM_ClassName_method(thing, args)
Boolean (type & constants)  | bool                               | integer*4
                            | true                               | parameter (OpenMM_True=1)
                            | false                              | parameter (OpenMM_False=0)
string                      std::string                          character(*)
3-vector                    OpenMM::Vec3                         real*8 vec(3)
arrays                      std::vector<std::string>             | type (OpenMM_StringArray)
                            std::vector<double>                  | type (OpenMM_DoubleArray)
                            std::vector<Vec3>                    | type (OpenMM_Vec3Array)
                            std::vector<std::pair<int,int>>      | type (OpenMM_BondArray)
                            std::map<std::string, double>        | type (OpenMM_ParameterArray)
==========================  ===================================  ========================================================

:autonumber:`Table,Fortran API`\ : Default mapping of objects from the C++ API to the Fortran API

Because there are no C++ API equivalents to the array types, they are described
in detail below.

OpenMM_Vec3 helper type
=======================

Unlike the other OpenMM objects which are opaque and manipulated via pointers,
the Fortran API uses an ordinary :code:`real*8(3)` array in
place of the :code:`OpenMM::Vec3` type.
You can work directly with the individual elements of this type from your
Fortran program if you want. For convenience, a :code:`scale()` function is
provided that creates a new Vec3 from an old one and a scale factor:

.. code-block:: fortran

    subroutine OpenMM_Vec3_scale(vec, scale, result)
    real*8 vec(3), scale, result(3)

No explicit :code:`type`\ :code:`(OpenMM_Vec3)` is provided in the Fortran
API since it is not needed.

Array helper types
==================

C++ has built-in container types :code:`std::vector` and :code:`std::map`
which OpenMM uses to manipulate arrays of objects. These don’t have direct
equivalents in Fortran, so we supply special array types for each kind of object
for which OpenMM creates containers. These are: string, double, Vec3, bond, and
parameter map. See :autonumref:`Table,Fortran arrays` for the names of the Fortran types for each of
these object arrays. Each of the array types provides these functions (prefixed
by :code:`OpenMM_` and the actual *Thing* name), with the syntax shown
conceptually since it differs slightly for each kind of object.

+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| Function                                  | Operation                                                                                              |
+===========================================+========================================================================================================+
| | subroutine create(array,size)           | Create a heap-allocated array of *Things*\ , with space pre-allocated to hold :code:`size` of them.    |
| | type (OpenMM\_\ *Thing*\ Array) array   | You can start at :code:`size`\ ==0 if you want since these arrays are dynamically resizeable.          |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine destroy(array)               | Free the heap space that is currently in use for the passed-in array of *Things*\ .                    |
| | type (OpenMM\_\ *Thing*\ Array) array   |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | function getSize(array)                 | Return the current number of *Things* in this array. This means you can :code:`get()` and              |
| | type (OpenMM\_\ *Thing*\ Array) array   | :code:`set()` elements up to :code:`getSize()`\ .                                                      |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine resize(array,size)           | Change the size of this array to the indicated value which may be smaller or larger than the           |
| | type (OpenMM\_\ *Thing*\ Array) array   | current size. Existing elements remain in their same locations as long as they still fit.              |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine append(array,elt)            | Add a *Thing* to the end of the array, increasing the array size by one. The precise syntax depends    |
| | type (OpenMM\_\ *Thing*\ Array) array   | on the actual type of *Thing*\ ; see below.                                                            |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine set(array,index,elt)         | Store a copy of :code:`elt` in the indicated element of the array (indexed from 1). The array must     |
| | type (OpenMM\_\ *Thing*\ Array) array   | be of length at least :code:`index`\ ; you can’t grow the array with this function.                    |
| | integer*4 size                          |                                                                                                        |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine get(array,index,elt)         | Retrieve a particular element from the array (indexed from 1).  Some *Things* require more than one    |
| | type (OpenMM\_\ *Thing*\ Array) array   | argument to return.                                                                                    |
| | integer*4 size                          |                                                                                                        |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+

:autonumber:`Table,Fortran arrays`\ : Generic description of array helper types

Here are the exact declarations with deviations from the generic description
noted, for each of the array types.

OpenMM_DoubleArray
------------------

.. code-block:: fortran

    subroutine OpenMM_DoubleArray_create(array, size)
        integer*4 size
        type (OpenMM_DoubleArray) array
    subroutine OpenMM_DoubleArray_destroy(array)
        type (OpenMM_DoubleArray) array
    function OpenMM_DoubleArray_getSize(array)
        type (OpenMM_DoubleArray) array
        integer*4 OpenMM_DoubleArray_getSize
    subroutine OpenMM_DoubleArray_resize(array, size)
        type (OpenMM_DoubleArray) array
        integer*4 size
    subroutine OpenMM_DoubleArray_append(array, value)
        type (OpenMM_DoubleArray) array
        real*8 value
    subroutine OpenMM_DoubleArray_set(array, index, value)
        type (OpenMM_DoubleArray) array
        integer*4 index
        real*8 value
    subroutine OpenMM_DoubleArray_get(array, index, value)
        type (OpenMM_DoubleArray) array
        integer*4 index
        real*8 value

OpenMM_StringArray
------------------

.. code-block:: fortran

    subroutine OpenMM_StringArray_create(array, size)
        integer*4 size
        type (OpenMM_StringArray) array
    subroutine OpenMM_StringArray_destroy(array)
        type (OpenMM_StringArray) array
    function OpenMM_StringArray_getSize(array)
        type (OpenMM_StringArray) array
        integer*4 OpenMM_StringArray_getSize
    subroutine OpenMM_StringArray_resize(array, size)
        type (OpenMM_StringArray) array
        integer*4 size
    subroutine OpenMM_StringArray_append(array, str)
        type (OpenMM_StringArray) array
        character(*) str
    subroutine OpenMM_StringArray_set(array, index, str)
        type (OpenMM_StringArray) array
        integer*4 index
        character(*) str
    subroutine OpenMM_StringArray_get(array, index, str)
        type (OpenMM_StringArray) array
        integer*4 index
        character(*)str

OpenMM_Vec3Array
----------------

.. code-block:: fortran

    subroutine OpenMM_Vec3Array_create(array, size)
        integer*4 size
        type (OpenMM_Vec3Array) array
    subroutine OpenMM_Vec3Array_destroy(array)
        type (OpenMM_Vec3Array) array
    function OpenMM_Vec3Array_getSize(array)
        type (OpenMM_Vec3Array) array
        integer*4 OpenMM_Vec3Array_getSize
    subroutine OpenMM_Vec3Array_resize(array, size)
        type (OpenMM_Vec3Array) array
        integer*4 size
    subroutine OpenMM_Vec3Array_append(array, vec)
        type (OpenMM_Vec3Array) array
        real*8 vec(3)
    subroutine OpenMM_Vec3Array_set(array, index, vec)
        type (OpenMM_Vec3Array) array
        integer*4 index
        real*8 vec(3)
    subroutine OpenMM_Vec3Array_get(array, index, vec)
        type (OpenMM_Vec3Array) array
        integer*4 index
        real*8 vec (3)

OpenMM_BondArray
----------------

Note that bonds are specified by pairs of integers (the atom indices). The
:code:`get()` method returns those in a pair of final arguments rather than as
its functional return.

.. code-block:: fortran

    subroutine OpenMM_BondArray_create(array, size)
        integer*4 size
        type (OpenMM_BondArray) array
    subroutine OpenMM_BondArray_destroy(array)
        type (OpenMM_BondArray) array
    function OpenMM_BondArray_getSize(array)
        type (OpenMM_BondArray) array
        integer*4 OpenMM_BondArray_getSize
    subroutine OpenMM_BondArray_resize(array, size)
        type (OpenMM_BondArray) array
        integer*4 size
    subroutine OpenMM_BondArray_append(array, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 particle1, particle2
    subroutine OpenMM_BondArray_set(array, index, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 index, particle1, particle2
    subroutine OpenMM_BondArray_get(array, index, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 index, particle1, particle2

OpenMM_ParameterArray
---------------------

OpenMM returns references to internal :code:`ParameterArrays` but does not
support user-created :code:`ParameterArrays`\ , so only the :code:`get()`
and :code:`getSize()` functions are available. Also, note that since this is
actually a map rather than an array, the “index” is the *name* of the
parameter rather than its ordinal.

.. code-block:: fortran

    function OpenMM_ParameterArray_getSize(array)
        type (OpenMM_ParameterArray) array
        integer*4 OpenMM_ParameterArray_getSize
    subroutine OpenMM_ParameterArray_get(array, name, param)
        type (OpenMM_ParameterArray) array
        character(*) name
        character(*) param


Python API
**********


Mapping from the C++ API to the Python API
==========================================

The Python API follows the C++ API as closely as possible. There are three
notable differences:

#. The :code:`getState()` method in the :code:`Context` class takes
   Pythonic-type arguments to indicate which state variables should be made
   available.  For example:
   ::

    myContext.getState(getEnergy=True, getForce=False, …)

#. Wherever the C++ API uses references to return multiple values from a method,
   the Python API returns a tuple.  For example, in C++ you would query a
   HarmonicBondForce for a bond’s parameters as follows:
   ::

    int particle1, particle2;
    double length, k;
    f.getBondParameters(i, particle1, particle2, length, k);

   In Python, the equivalent code is:
   ::

    [particle1, particle2, length, k] = f.getBondParameters(i)

#. Unlike C++, the Python API accepts and returns quantities with units attached
   to most values (see Section :numref:`units-and-dimensional-analysis` below for
   details).  In short, this means that while values in C++ have *implicit*
   units, the Python API returns objects that have values and *explicit* units.


Mechanics of using the Python API
=================================

When using the Python API, be sure to include the GPU support
libraries in your library path, just as you would for a C++ application.  This
is set with the :code:`LD_LIBRARY_PATH` environment variable on Linux,
:code:`DYLD_LIBRARY_PATH` on Mac, or :code:`PATH` on Windows.  See
Chapter :numref:`installing-openmm` for details.

The Python API is contained in the openmm package, while the units code is
contained in the openmm.units package.  (The application layer, described in the
Application Guide, is contained in the openmm.app package.)  A program
using it will therefore typically begin
::

    import openmm as mm
    import openmm.unit as unit

Creating and using OpenMM objects is then done exactly as in C++:
::

    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(1.2*unit.nanometer)
    system.addForce(nb)

Note that when setting the cutoff distance, we explicitly specify that it is in
nanometers.  We could just as easily specify it in different units:
::

    nb.setCutoffDistance(12*unit.angstrom)

The use of units in OpenMM is discussed in the next section.


.. _units-and-dimensional-analysis:

Units and dimensional analysis
==============================


Why does the Python API include units?
--------------------------------------

The C++ API for OpenMM uses an *implicit* set of units for physical
quantities such as lengths, masses, energies, etc.  These units are based on
daltons, nanometers, and picoseconds for the mass, length, and time dimensions,
respectively.  When using the C++ API, it is very important to ensure that
quantities being manipulated are always expressed in terms of these units.  For
example, if you read in a distance in Angstroms, you must multiply that distance
by a conversion factor to turn it into nanometers before using it in the C++
API.  Such conversions can be a source of tedium and errors.  This is true in
many areas of scientific programming.  Units confusion was blamed for the loss
of the Mars Climate Orbiter spacecraft in 1999, at a cost of more than $100
million.  Units were introduced in the Python API to minimize the chance of such
errors.

The Python API addresses the potential problem of conversion errors by using
quantities with explicit units.  If a particular distance is expressed in
Angstroms, the Python API will know that it is in Angstroms.  When the time
comes to call the C++ API, it will understand that the quantity must be
converted to nanometers.  You, the programmer, must declare upfront that the
quantity is in Angstrom units, and the API will take care of the details from
then on.  Using explicit units is a bit like brushing your teeth: it requires
some effort upfront, but it probably saves you trouble in the long run.

Quantities, units, and dimensions
---------------------------------

The explicit unit system is based on three concepts: Dimensions, Units, and
Quantities.

Dimensions are measurable physical concepts such as mass, length, time, and
energy.  Energy is actually a composite dimension based on mass, length, and
time.

A Unit defines a linear scale used to measure amounts of a particular physical
Dimension.  Examples of units include meters, seconds, joules, inches, and
grams.

A Quantity is a specific amount of a physical Dimension.  An example of a
quantity is “0.63 kilograms”.  A Quantity is expressed as a combination of a
value (e.g., 0.63), and a Unit (e.g., kilogram).  The same Quantity can be
expressed in different Units.

The set of BaseDimensions defined in the openmm.unit module includes:

* mass
* length
* time
* temperature
* amount
* charge
* luminous intensity


These are not precisely the same list of base dimensions used in the SI unit
system.  SI defines “current” (charge per time) as a base unit, while openmm.unit
uses “charge”.  And openmm.unit treats angle as a dimension, even though angle
quantities are often considered dimensionless.  In this case, we choose to err
on the side of explicitness, particularly because interconversion of degrees and
radians is a frequent source of unit headaches.

Units examples
--------------

Many common units are defined in the openmm.unit module.
::

    from openmm.unit import nanometer, angstrom, dalton

Sometimes you don’t want to type the full unit name every time, so you can
assign it a shorter name using the :code:`as` functionality:
::

    from openmm.unit import nanometer as nm

New quantities can be created from a value and a unit.  You can use either the
multiply operator (‘*’) or the explicit Quantity constructor:
::

    from simk.unit import nanometer, Quantity
    # construct a Quantity using the multiply operator
    bond_length = 1.53 * nanometer
    # equivalently using the explicit Quantity constructor
    bond_length = Quantity(1.53, nanometer)
    # or more verbosely
    bond_length = Quantity(value=1.53, unit=nanometer)

Arithmetic with units
---------------------

Addition and subtraction of quantities is only permitted between quantities that
share the same dimension.  It makes no sense to add a mass to a distance.  If
you attempt to add or subtract two quantities with different dimensions, an
exception will be raised.  This is a good thing; it helps you avoid errors.
::

    x = 5.0*dalton + 4.3*nanometer; # error

Addition or subtraction of quantities with the same dimension, but different
units, is fine, and results in a new quantity created using the correct
conversion factor between the units used.
::

    x = 1.3*nanometer + 5.6*angstrom; # OK, result in nanometers

Quantities can be added and subtracted.  Naked Units cannot.

Multiplying or dividing two quantities creates a new quantity with a composite
dimension.  For example, dividing a distance by a time results in a velocity.
::

    from openmm.unit import kilogram, meter, second
    a = 9.8 * meter / second**2; # acceleration
    m = 0.36 * kilogram; # mass
    F = m * a; # force in kg*m/s**2::


Multiplication or division of two Units results in a composite Unit.
::

    mps = meter / second

Unlike amount (moles), angle (radians) is arguably dimensionless.  But openmm.unit
treats angle as another dimension.   Use the trigonometric functions from the
openmm.unit module (not those from the Python math module!) when dealing with
Units and Quantities.
::

    from openmm.unit import sin, cos, acos
    x = sin(90.0*degrees)
    angle = acos(0.68); # returns an angle quantity (in radians)

The method :code:`pow()` is a built-in Python method that works with
Quantities and Units.
::

    area = pow(3.0*meter, 2)
    # or, equivalently
    area = (3.0*meter)**2
    # or
    area = 9.0*(meter**2)

The method :code:`sqrt()` is not as built-in as :code:`pow()`\ .  Do not
use the Python :code:`math.sqrt()` method with Units and Quantities.  Use
the :code:`openmm.unit.sqrt()` method instead:
::

    from openmm.unit import sqrt
    side_length = sqrt(4.0*meter**2)


Atomic scale mass and energy units are “per amount”
---------------------------------------------------

Mass and energy units at the atomic scale are specified “per amount” in the
openmm.unit module.  Amount (mole) is one of the seven fundamental dimensions in
the SI unit system.   The atomic scale mass unit, dalton, is defined as grams
per mole.  The dimension of dalton is therefore mass/amount, instead of simply
mass.  Similarly, the atomic scale energy unit, kilojoule_per_mole (and
kilocalorie_per_mole) has “per amount” in its dimension.  Be careful to always
use “per amount” mass and energy types at the atomic scale, and your dimensional
analysis should work out properly.

The energy unit kilocalories_per_mole does not have the same Dimension as the
macroscopic energy unit kilocalories.  Molecular scientists sometimes use the
word "kilocalories" when they mean "kilocalories per mole".  Use "kilocalories
per mole" or"kilojoules per mole" for molecular energies.  Use "kilocalories"
for the metabolic energy content of your lunch.  The energy unit
kilojoule_per_mole happens to go naturally with the units nanometer,
picoseconds, and dalton.  This is because 1 kilojoule/mole happens to be equal
to 1 gram-nanometer\ :sup:`2`\ /mole-picosecond\ :sup:`2`\ , and is therefore
consistent with the molecular dynamics unit system used in the C++ OpenMM API.

These "per mole" units are what you should be using for molecular calculations,
as long as you are using SI / cgs / calorie sorts of units.

SI prefixes
-----------

Many units with SI prefixes such as “milligram” (milli) and “kilometer” (kilo)
are provided in the openmm.unit module.  Others can be created by multiplying a
prefix symbol by a non-prefixed unit:
::

    from openmm.unit import mega, kelvin
    megakelvin = mega * kelvin
    t = 8.3 * megakelvin

Only grams and meters get all of the SI prefixes (from yotto-(10\ :sup:`-24`\ )
to yotta-(10\ :sup:`24`\ )) automatically.


Converting to different units
-----------------------------

Use the :code:`Quantity.in_units_of()` method to create a new Quantity with
different units.
::

    from openmm.unit import nanosecond, fortnight
    x = (175000*nanosecond).in_units_of(fortnight)

When you want a plain number out of a Quantity, use the :code:`value_in_unit()` method:
::

    from openmm.unit import femtosecond, picosecond
    t = 5.0*femtosecond
    t_just_a_number = t.value_in_unit(picoseconds)

Using :code:`value_in_unit()` puts the responsibility for unit analysis back
into your hands, and it should be avoided.  It is sometimes necessary, however,
when you are called upon to use a non-units-aware Python API.


Lists, tuples, vectors, numpy arrays, and Units
-----------------------------------------------

Units can be attached to containers of numbers to create a vector quantity.  The
openmm.unit module overloads the :code:`__setitem__` and
:code:`__getitem__` methods for these containers to ensure that Quantities go
in and out.
::

    >>> a = Vec3(1,2,3) * nanometers
    >>> print(a)
    (1, 2, 3) nm
    >>> print(a.in_units_of(angstroms))
    (10.0, 20.0, 30.0) A

    >>> s2 = [[1,2,3],[4,5,6]] * centimeter
    >>> print(s2)
    [[1, 2, 3], [4, 5, 6]] cm
    >>> print(s2/millimeter)
    [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]

    >>> import numpy
    >>> a = numpy.array([1,2,3]) * centimeter
    >>> print(a)
    [1 2 3] cm
    >>> print(a/millimeter)
    [ 10.  20.  30.]

Converting a whole list to different units at once is much faster than
converting each element individually.  For example, consider the following code
that prints out the position of every particle in a State, as measured in
Angstroms:
::

    for v in state.getPositions():
        print(v.value_in_unit(angstrom))

This can be rewritten as follows:
::

    for v in state.getPositions().value_in_unit(angstrom):
        print(v)

The two versions produce identical results, but the second one will run faster,
and therefore is preferred.
.. _amoeba-plugin:

AMOEBA Plugin
#############

OpenMM |version| provides a plugin that implements the AMOEBA polarizable atomic
multipole force field from Jay Ponder’s lab. The AMOEBA force field may be used
through OpenMM’s Python application layer. We have also created a modified
version of TINKER (referred to as TINKER-OpenMM here) that uses OpenMM to
accelerate AMOEBA simulations. TINKER-OpenMM can be created from a TINKER
package using three files made available through the OpenMM home page. OpenMM
AMOEBA Force and System objects containing AMOEBA forces can be serialized.

At present, AMOEBA is only supported on the CUDA and Reference platforms, not on
the OpenCL platform.

In the following sections, the individual forces and options available in the
plugin are listed, and the steps required to build and use the plugin and
TINKER-OpenMM are outlined. Validation results are also reported.  Benchmarks
can be found on the OpenMM wiki at http://wiki.simtk.org/openmm/Benchmarks.

OpenMM AMOEBA Supported Forces and Options
*******************************************


.. _supported-forces-and-options:

Supported Forces and Options
============================

The AMOEBA force terms implemented in OpenMM are listed in :autonumref:`Table,mapping from TINKER` along
with the supported and unsupported options. TINKER options that are not
supported for any OpenMM force include the grouping of atoms (e.g. protein
chains), the infinite polymer check, and no exclusion of particles from
energy/force calculations (‘active’/’inactive’ particles).  The virial is not
calculated for any force.

All rotation axis types are supported: ‘Z-then-X’, ‘Bisector’, ‘Z-Bisect’,
‘3-Fold’, ‘Z-Only’.


=================================  ==================================  ======================================================================================================================================================================================
TINKER Force                       OpenMM Force                        Option/Note
=================================  ==================================  ======================================================================================================================================================================================
ebond1 (bondterm)                  AmoebaBondForce                     bndtyp='HARMONIC' supported, 'MORSE' not implemented
Eangle71 (angleterm)               AmoebaAngleForce                    angtyp='HARMONIC' and 'IN-PLANE' supported; 'LINEAR' and 'FOURIER' not implemented
etors1a (torsionterm)              PeriodicTorsionForce                All options implemented; smoothing version(etors1b) not supported
etortor1 (tortorterm)              AmoebaTorsionTorsionForce           All options implemented
eopbend1 (opbendterm)              AmoebaOutOfPlaneBendForce           opbtyp = 'ALLINGER' implemented; 'W-D-C' not implemented
epitors1 (pitorsterm)              AmoebaPiTorsionForce                All options implemented
estrbnd1 (strbndterm)              AmoebaStretchBendForce              All options implemented
ehal1a (vdwterm)                   AmoebaVdwForce                      ehal1b(LIGHTS) not supported
empole1a (mpoleterm)               AmoebaMultipoleForce                poltyp = 'MUTUAL', 'DIRECT'  supported
empole1c (mpoleterm) PME           AmoebaMultipoleForce                poltyp = 'MUTUAL', 'DIRECT' supported; boundary= 'VACUUM' unsupported
esolv1 (solvateterm)               | AmoebaWcaDispersionForce,         Only born-radius=’grycuk’ and solvate=’GK’ supported; unsupported solvate settings:
                                   | AmoebaGeneralizedKirkwoodForce    ‘ASP’, ‘SASA’, ‘ONION’, ‘pb’, 'GB-HPMF’, 'Gk-HPMF’; SASA computation is based on ACE approximation
eurey1 (ureyterm)                  HarmonicBondForce                   All options implemented
=================================  ==================================  ======================================================================================================================================================================================

:autonumber:`Table,mapping from TINKER`\ :  Mapping between TINKER and OpenMM AMOEBA forces


Some specific details to be aware of are the following:

* Forces available in TINKER but not implemented in the OpenMM AMOEBA plugin
  include the following: angle-angle, out-of-plane distance, improper dihedral,
  improper torsion, stretch-torsion, charge-charge, atomwise charge-dipole,
  dipole-dipole, reaction field, ligand field, restraint, scf molecular orbital
  calculation; strictly speaking, these are not part of the AMOEBA force field.

* Implicit solvent in TINKER-OpenMM is implemented with key file entry ‘solvate
  GK’.  The entry ‘born-radius grycuk’ should also be included; only the ‘grycuk’
  option for calculating the Born radii is available in the plugin.

* In TINKER, the nonpolar cavity contribution to the solvation term is
  calculated using an algorithm that does not map well to GPUs.  Instead the
  OpenMM plugin uses the TINKER version of the ACE approximation to estimate the
  cavity contribution to the SASA.

* Calculations using the CUDA platform may be done in either single or double
  precision; for the Reference platform, double precision is used.  TINKER uses
  double precision.

* The TINKER parameter files for the AMOEBA force-field parameters are based on
  units of kilocalorie/Å, whereas OpenMM uses units of kilojoules/nanometer; both
  TINKER and OpenMM use picoseconds time units. Hence, in mapping the force-field
  parameters from TINKER files to OpenMM, many of the parameter values must be
  converted to the OpenMM units. The setup methods in the TINKER-OpenMM
  application perform the required conversions.


Supported Integrators
=====================

In addition to the limitations to the forces outlined above, TINKER-OpenMM can
only use either the ‘Verlet’ or ‘Stochastic’ integrators when the OpenMM plugin
is used; an equivalent to the TINKER ‘Beeman’ integrator is unavailable in
OpenMM.

OpenMM AMOEBA Validation
************************

OpenMM and TINKER 6.1.01 were each used to compute the atomic forces for
dihydrofolate reductase (DHFR) in implicit and explicit solvent.  Calculations
used the CUDA platform, and were repeated for both single and double precision.
For every atom, the relative difference between OpenMM and TINKER was computed
as 2·\|F\ :sub:`MM`\ –F\ :sub:`T`\ \|/(\|F\ :sub:`MM`\ \|+\|F\ :sub:`T`\ \|), where
F\ :sub:`MM` is the force computed by OpenMM and F\ :sub:`T` is the force
computed by TINKER.  The median over all atoms is shown in :autonumref:`Table,comparison to TINKER`\ .

Because OpenMM and TINKER use different approximations to compute the cavity
term, the differences in forces are much larger for implicit solvent than for
explicit solvent.  We therefore repeated the calculations, removing the cavity
term.  This yields much closer agreement between OpenMM and TINKER,
demonstrating that the difference comes entirely from that one term.

=========================  ==========================  ===================
Solvent Model              single                      double
=========================  ==========================  ===================
Implicit                   1.04·10\ :sup:`-2`          1.04·10\ :sup:`-2`
Implicit (no cavity term)  9.23·10\ :sup:`-6`          1.17·10\ :sup:`-6`
Explicit                   3.73·10\ :sup:`-5`          1.83·10\ :sup:`-7`
=========================  ==========================  ===================

:autonumber:`Table,comparison to TINKER`\ :  Median relative difference in forces between OpenMM and TINKER

.. _the-openmm-library-introduction:

Introduction
############


What Is the OpenMM Library?
***************************

OpenMM consists of two parts.  First, there is a set of libraries for performing
many types of computations needed for molecular simulations: force evaluation,
numerical integration, energy minimization, etc.  These libraries provide an
interface targeted at developers of simulation software, allowing them to easily
add simulation features to their programs.

Second, there is an “application layer”, a set of Python libraries providing a
high level interface for running simulations.  This layer is targeted at
computational biologists or other people who want to run simulations, and who
may or may not be programmers.

The first part of this guide focused on the application layer and described how to run
simulations with it.  We now turn to the lower level libraries.  We will assume
you are a programmer, that you are writing your own applications, and that you
want to add simulation features to those applications.  The following chapters
describe how to do that with OpenMM.

How to get started
==================

We have provided a number of files that make it easy to get started with OpenMM.
Pre-compiled binaries are provided for quickly getting OpenMM onto your computer
(See Chapter :numref:`installing-openmm` for set-up instructions).  We recommend that you then
compile and run some of the tutorial examples, described in Chapter :numref:`openmm-tutorials`.
These highlight key functions within OpenMM and teach you the basic programming concepts for using
OpenMM.  Once you are ready to begin integrating OpenMM into a specific software package, read
through Chapter :numref:`examples-of-openmm-integration` to see how other software developers have
done this.

License
========

Two different licenses are used for different parts of OpenMM.  The public API,
the low level API, the reference platform, the CPU platform, and the application
layer are all distributed under the MIT
license.  This is a very permissive license which allows them to be used in
almost any way, requiring only that you retain the copyright notice and
disclaimer when distributing them.

The CUDA and OpenCL platforms are distributed under the GNU Lesser General
Public License (LGPL).  This also allows you to use, modify, and distribute them
in any way you want, but it requires you to also distribute the source code for
your modifications.  This restriction applies only to modifications to OpenMM
itself; you need not distribute the source code to applications that use it.

OpenMM also uses several pieces of code that were written by other people and
are covered by other licenses.  All of these licenses are similar in their terms
to the MIT license, and do not significantly restrict how OpenMM can be used.

All of these licenses may be found in the “licenses” directory included with
OpenMM.


Design Principles
*****************

The design of the OpenMM API is guided by the following principles.

1. The API must support efficient implementations on a variety of architectures.

The most important consequence of this goal is that the API cannot provide
direct access to state information (particle positions, velocities, etc.) at all
times.  On some architectures, accessing this information is expensive.  With a
GPU, for example, it will be stored in video memory, and must be transferred to
main memory before outside code can access it.  On a distributed architecture,
it might not even be present on the local computer.  OpenMM therefore only
allows state information to be accessed in bulk, with the understanding that
doing so may be a slow operation.

2. The API should be easy to understand and easy to use.

This seems obvious, but it is worth stating as an explicit goal.  We are
creating OpenMM with the hope that many other people will use it.  To achieve
that goal, it should be possible for someone to learn it without an enormous
amount of effort.  An equally important aspect of being “easy to use” is being
easy to use *correctly*\ .  A well designed API should minimize the
opportunities for a programmer to make mistakes.  For both of these reasons,
clarity and simplicity are essential.

3. It should be modular and extensible.

We cannot hope to provide every feature any user will ever want.  For that
reason, it is important that OpenMM be easy to extend.  If a user wants to add a
new molecular force field, a new thermostat algorithm, or a new hardware
platform, the API should make that easy to do.

4. The API should be hardware independent.

Computer architectures are changing rapidly, and it is impossible to predict
what hardware platforms might be important to support in the future.  One of the
goals of OpenMM is to separate the API from the hardware.  The developers of a
simulation application should be able to write their code once, and have it
automatically take advantage of any architecture that OpenMM supports, even
architectures that do not yet exist when they write it.

Choice of Language
******************

Molecular modeling and simulation tools are written in a variety of languages:
C, C++, Fortran, Python, TCL, etc.  It is important that any of these tools be
able to use OpenMM.  There are two possible approaches to achieving this goal.

One option is to provide a separate version of the API for each language.  These
could be created by hand, or generated automatically with a wrapper generator
such as SWIG.  This would require the API to use only “lowest common
denominator” features that can be reasonably supported in all languages.  For
example, an object oriented API would not be an option, since it could not be
cleanly expressed in C or Fortran.

The other option is to provide a single version of the API written in a single
language.  This would permit a cleaner, simpler API, but also restrict the
languages it could be directly called from.  For example, a C++ API could not be
invoked directly from Fortran or Python.

We have chosen to use a hybrid of these two approaches.  OpenMM is based on an
object oriented C++ API.  This is the primary way to invoke OpenMM, and is the
only API that fully exposes all features of the library.  We believe this will
ultimately produce the best, easiest to use API and create the least work for
developers who use it.  It does require that any code which directly invokes
this API must itself be written in C++, but this should not be a significant
burden.  Regardless of what language we had chosen, developers would need to
write a thin layer for translating between their own application’s data model
and OpenMM.  That layer is the only part which needs to be written in C++.

In addition, we have created wrapper APIs that allow OpenMM to be invoked from
other languages.  The current release includes wrappers for C, Fortran, and
Python.  These wrappers support as many features as reasonably possible given
the constraints of the particular languages, but some features cannot be fully
supported.  In particular, writing plug-ins to extend the OpenMM API can only be
done in C++.

We are also aware that some features of C++ can easily lead to compatibility and
portability problems, and we have tried to avoid those features.  In particular,
we make minimal use of templates and avoid multiple inheritance altogether.  Our
goal is to support OpenMM on all major compilers and operating systems.

Architectural Overview
**********************

OpenMM is based on a layered architecture, as shown in the following diagram:


.. figure:: ../../images/ArchitectureLayers.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,OpenMM architecture`:  OpenMM architecture

At the highest level is the OpenMM public API.  This is the API developers
program against when using OpenMM within their own applications.  It is designed
to be simple, easy to understand, and completely platform independent.  This is
the only layer that many users will ever need to look at.

The public API is implemented by a layer of platform independent code.  It
serves as the interface to the lower level, platform specific code.  Most users
will never need to look at it.

The next level down is the OpenMM Low Level API (OLLA).  This acts as an
abstraction layer to hide the details of each hardware platform.  It consists of
a set of C++ interfaces that each platform must implement.  Users who want to
extend OpenMM will need to write classes at the OLLA level.  Note the different
roles played by the public API and the low level API: the public API defines an
interface for users to invoke in their own code, while OLLA defines an interface
that users must implement, and that is invoked by the OpenMM implementation
layer.

At the lowest level is hardware specific code that actually performs
computations.  This code may be written in any language and use any technologies
that are appropriate.  For example, code for GPUs will be written in stream
processing languages such as OpenCL or CUDA, code written to run on clusters
will use MPI or other distributed computing tools, code written for multicore
processors will use threading tools such as Pthreads or OpenMP, etc.  OpenMM
sets no restrictions on how these computational kernels are written.  As long as
they are wrapped in the appropriate OLLA interfaces, OpenMM can use them.

.. _the-openmm-public-api:

The OpenMM Public API
*********************

The public API is based on a small number of classes:

**System**\ : A System specifies generic properties of the system to be
simulated: the number of particles it contains, the mass of each one, the size
of the periodic box, etc.  The interactions between the particles are specified
through a set of Force objects (see below) that are added to the System.  Force
field specific parameters, such as particle charges, are not direct properties
of the System.  They are properties of the Force objects contained within the
System.

**Force**\ : The Force objects added to a System define the behavior of the
particles.  Force is an abstract class; subclasses implement specific behaviors.
The Force class is actually slightly more general than its name suggests.  A
Force can, indeed, apply forces to particles, but it can also directly modify
particle positions and velocities in arbitrary ways.  Some thermostats and
barostats, for example, can be implemented as Force classes.  Examples of Force
subclasses include HarmonicBondForce, NonbondedForce, and MonteCarloBarostat.

**Context**\ : This stores all of the state information for a simulation:
particle positions and velocities, as well as arbitrary parameters defined by
the Forces in the System.  It is possible to create multiple Contexts for a
single System, and thus have multiple simulations of that System in progress at
the same time.

**Integrator**\ : This implements an algorithm for advancing the simulation
through time.  It is an abstract class; subclasses implement specific
algorithms.  Examples of Integrator subclasses include LangevinIntegrator,
VerletIntegrator, and BrownianIntegrator.

**State**\ : A State stores a snapshot of the simulation at a particular point
in time.  It is created by calling a method on a Context.  As discussed earlier,
this is a potentially expensive operation.  This is the only way to query the
values of state variables, such as particle positions and velocities; Context
does not provide methods for accessing them directly.

Here is an example of what the source code to create a System and run a
simulation might look like:

.. code-block:: c

    System system;
    for (int i = 0; i < numParticles; ++i)
        system.addParticle(particle[i].mass);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    for (int i = 0; i < numBonds; ++i)
        bonds->addBond(bond[i].particle1, bond[i].particle2,
            bond[i].length, bond[i].k);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    system.addForce(angles);
    for (int i = 0; i < numAngles; ++i)
        angles->addAngle(angle[i].particle1, angle[i].particle2,
            angle[i].particle3, angle[i].angle, angle[i].k);
    // ...create and initialize other force field terms in the same way
    LangevinMiddleIntegrator integrator(temperature, friction, stepSize);
    Context context(system, integrator);
    context.setPositions(initialPositions);
    context.setVelocities(initialVelocities);
    integrator.step(10000);

We create a System, add various Forces to it, and set parameters on both the
System and the Forces.  We then create a LangevinMiddleIntegrator, initialize a
Context in which to run a simulation, and instruct the Integrator to advance the
simulation for 10,000 time steps.

The OpenMM Low Level API
************************

The OpenMM Low Level API (OLLA) defines a set of interfaces that users must
implement in their own code if they want to extend OpenMM, such as to create a
new Force subclass or support a new hardware platform.  It is based on the
concept of “kernels” that define particular computations to be performed.

More specifically, there is an abstract class called **KernelImpl**\ .
Instances of this class (or rather, of its subclasses) are created by
**KernelFactory** objects.  These classes provide the concrete implementations
of kernels for a particular platform.  For example, to perform calculations on a
GPU, one would create one or more KernelImpl subclasses that implemented the
computations with GPU kernels, and one or more KernelFactory subclasses to
instantiate the KernelImpl objects.

All of these objects are encapsulated in a single object that extends
**Platform**\ . KernelFactory objects are registered with the Platform to be
used for creating specific named kernels.  The choice of what implementation to
use (a GPU implementation, a multithreaded CPU implementation, an MPI-based
distributed implementation, etc.) consists entirely of choosing what Platform to
use.

As discussed so far, the low level API is not in any way specific to molecular
simulation; it is a fairly generic computational API.  In addition to defining
the generic classes, OpenMM also defines abstract subclasses of KernelImpl
corresponding to specific calculations.  For example, there is a class called
CalcHarmonicBondForceKernel to implement HarmonicBondForce and a class called
IntegrateLangevinStepKernel to implement LangevinIntegrator.  It is these
classes for which each Platform must provide a concrete subclass.

This architecture is designed to allow easy extensibility.  To support a new
hardware platform, for example, you create concrete subclasses of all the
abstract kernel classes, then create appropriate factories and a Platform
subclass to bind everything together.  Any program that uses OpenMM can then use
your implementation simply by specifying your Platform subclass as the platform
to use.

Alternatively, you might want to create a new Force subclass to implement a new
type of interaction.  To do this, define an abstract KernelImpl subclass
corresponding to the new force, then write the Force class to use it.  Any
Platform can support the new Force by providing a concrete implementation of
your KernelImpl subclass.  Furthermore, you can easily provide that
implementation yourself, even for existing Platforms created by other people.
Simply create a new KernelFactory subclass for your kernel and register it with
the Platform object.  The goal is to have a completely modular system.  Each
module, which might be distributed as an independent library, can either add new
features to existing platforms or support existing features on new platforms.

In fact, there is nothing “special” about the kernel classes defined by OpenMM.
They are simply KernelImpl subclasses that happen to be used by Forces and
Integrators that happen to be bundled with OpenMM.  They are treated exactly
like any other KernelImpl, including the ones you define yourself.

It is important to understand that OLLA defines an interface, not an
implementation.  It would be easy to assume a one-to-one correspondence between
KernelImpl objects and the pieces of code that actually perform calculations,
but that need not be the case.  For a GPU implementation, for example, a single
KernelImpl might invoke several GPU kernels.  Alternatively, a single GPU kernel
might perform the calculations of several KernelImpl subclasses.

.. _platforms:

Platforms
*********

This release of OpenMM contains the following Platform subclasses:

**ReferencePlatform**\ : This is designed to serve as reference code for
writing other platforms.  It is written with simplicity and clarity in mind, not
performance.

**CpuPlatform**\ : This platform provides high performance when running on
conventional CPUs.

**CudaPlatform**\ : This platform is implemented using the CUDA language, and
performs calculations on Nvidia GPUs.

**OpenCLPlatform**\ : This platform is implemented using the OpenCL language,
and performs calculations on a variety of types of GPUs and CPUs.

The choice of which platform to use for a simulation depends on various factors:

#. The Reference platform is much slower than the others, and therefore is
   rarely used for production simulations.
#. The CPU platform is usually the fastest choice when a fast GPU is not
   available.  However, it requires the CPU to support SSE 4.1.  That includes most
   CPUs made in the last several years, but this platform may not be available on
   some older computers.  Also, for simulations that use certain features
   (primarily the various “custom” force classes), it may be faster to use the
   OpenCL platform running on the CPU.
#. The CUDA platform can only be used with NVIDIA GPUs.  For using an AMD or
   Intel GPU, use the OpenCL platform.
#. The AMOEBA force field only works with the CUDA platform, not with the OpenCL
   platform.  It also works with the Reference and CPU platforms, but the performance
   is usually too slow to be useful on those platforms.
.. _platform-specific-properties:

Platform-Specific Properties
############################

When creating a Context, you can specify values for properties specific to a
particular Platform.  This is used to control how calculations are done in ways
that are outside the scope of the generic OpenMM API.

To do this, pass both the Platform object and a map of property values to the
Context constructor:

.. code-block:: c

    Platform& platform = Platform::getPlatformByName("OpenCL");
    map<string, string> properties;
    properties["DeviceIndex"] = "1";
    Context context(system, integrator, platform, properties);

After a Context is created, you can use the Platform’s \
:code:`getPropertyValue()` method to query the values of properties.

OpenCL Platform
***************

The OpenCL Platform recognizes the following Platform-specific properties:

* Precision: This selects what numeric precision to use for calculations.
  The allowed values are “single”, “mixed”, and “double”.  If it is set to
  “single”, nearly all calculations are done in single precision.  This is the
  fastest option but also the least accurate.  If it is set to “mixed”, forces are
  computed in single precision but integration is done in double precision.  This
  gives much better energy conservation with only a slight decrease in speed.
  If it is set to “double”, all calculations are done in double precision.  This
  is the most accurate option, but is usually much slower than the others.
* UseCpuPme: This selects whether to use the CPU-based PME
  implementation.  The allowed values are “true” or “false”.  Depending on your
  hardware, this might (or might not) improve performance.  To use this option,
  you must have FFTW (single precision, multithreaded) installed, and your CPU
  must support SSE 4.1.
* OpenCLPlatformIndex: When multiple OpenCL implementations are installed on
  your computer, this is used to select which one to use.  The value is the
  zero-based index of the platform (in the OpenCL sense, not the OpenMM sense) to use,
  in the order they are returned by the OpenCL platform API.  This is useful, for
  example, in selecting whether to use a GPU or CPU based OpenCL implementation.
* DeviceIndex: When multiple OpenCL devices are available on your
  computer, this is used to select which one to use.  The value is the zero-based
  index of the device to use, in the order they are returned by the OpenCL device
  API.


The OpenCL Platform also supports parallelizing a simulation across multiple
GPUs.  To do that, set the DeviceIndex property to a comma separated list
of values.  For example,

.. code-block:: c

    properties["DeviceIndex"] = "0,1";

This tells it to use both devices 0 and 1, splitting the work between them.

CUDA Platform
*************

The CUDA Platform recognizes the following Platform-specific properties:

* Precision: This selects what numeric precision to use for calculations.
  The allowed values are “single”, “mixed”, and “double”.  If it is set to
  “single”, nearly all calculations are done in single precision.  This is the
  fastest option but also the least accurate.  If it is set to “mixed”, forces are
  computed in single precision but integration is done in double precision.  This
  gives much better energy conservation with only a slight decrease in speed.
  If it is set to “double”, all calculations are done in double precision.  This
  is the most accurate option, but is usually much slower than the others.
* UseCpuPme: This selects whether to use the CPU-based PME implementation.
  The allowed values are “true” or “false”.  Depending on your hardware, this
  might (or might not) improve performance.  To use this option, you must have
  FFTW (single precision, multithreaded) installed, and your CPU must support SSE
  4.1.
* CudaCompiler: This specifies the path to the CUDA kernel compiler.  Versions
  of CUDA before 7.0 require a separate compiler executable.  If you do
  not specify this, OpenMM will try to locate the compiler itself.  Specify this
  only when you want to override the default location.  The logic used to pick the
  default location depends on the operating system:

  * Mac/Linux: It first looks for an environment variable called
    OPENMM_CUDA_COMPILER.  If that is set, its value is used.  Otherwise, the
    default location is set to /usr/local/cuda/bin/nvcc.
  * Windows: It looks for an environment variable called CUDA_BIN_PATH, then
    appends \nvcc.exe to it.  That environment variable is set by the CUDA
    installer, so it usually is present.

* TempDirectory: This specifies a directory where temporary files can be
  written while compiling kernels.  OpenMM usually can locate your operating
  system’s temp directory automatically (for example, by looking for the TEMP
  environment variable), so you rarely need to specify this.
* DeviceIndex: When multiple CUDA devices are available on your computer,
  this is used to select which one to use.  The value is the zero-based index of
  the device to use, in the order they are returned by the CUDA API.
* UseBlockingSync: This is used to control how the CUDA runtime
  synchronizes between the CPU and GPU.  If this is set to “true” (the default),
  CUDA will allow the calling thread to sleep while the GPU is performing a
  computation, allowing the CPU to do other work.  If it is set to “false”, CUDA
  will spin-lock while the GPU is working.  Setting it to "false" can improve performance slightly,
  but also prevents the CPU from doing anything else while the GPU is working.
* DeterministicForces: In some cases, the CUDA platform may compute forces
  in ways that are not fully deterministic (typically differing in what order a
  set of numbers get added together).  This means that if you compute the forces
  twice for the same particle positions, there may be tiny differences in the
  results.  In most cases this is not a problem, but certain algorithms depend
  on forces being exactly reproducible to the last bit.  If you set this
  property to "true", it will instead do these calculations in a way that
  produces fully deterministic results, at the cost of a small decrease in
  performance.

The CUDA Platform also supports parallelizing a simulation across multiple GPUs.
To do that, set the DeviceIndex property to a comma separated list of
values.  For example,

.. code-block:: c

    properties["DeviceIndex"] = "0,1";

This tells it to use both devices 0 and 1, splitting the work between them.

CPU Platform
************

The CPU Platform recognizes the following Platform-specific properties:

* Threads: This specifies the number of CPU threads to use.  If you do not
  specify this, OpenMM will select a default number of threads as follows:

  * If an environment variable called OPENMM_CPU_THREADS is set, its value is
    used as the number of threads.
  * Otherwise, the number of threads is set to the number of logical CPU cores
    in the computer it is running on.

  Usually the default value works well.  This is mainly useful when you are
  running something else on the computer at the same time, and you want to
  prevent OpenMM from monopolizing all available cores.

.. _platform-specific-properties-determinism:

Determinism
***********

Whether a simulation is deterministic will depend on what plaform you run on in
addition to what settings/methods you use. For instance, as of this writing,
using PME on the Reference, OpenCL, and double-precision CUDA will result in
deterministic simulations. Single-precision CUDA and CPU platforms are not
deterministic in this case. However, none of this behavior is guaranteed in
future versions. In many cases it will still result in an identical trajectory.
If determinism is a critical for your needs, you should carefully check to
ensure that your settings and platform allow for this.
.. default-domain:: py

.. _model-building-and-editing:

Model Building and Editing
##########################

Sometimes you have a PDB file that needs some work before you can simulate it.
Maybe it doesn’t contain hydrogen atoms (which is common for structures
determined by X-ray crystallography), so you need to add them.  Or perhaps you
want to simulate the system in explicit water, but the PDB file doesn’t contain
water molecules.  Or maybe it does contain water molecules, but they contain the
wrong number of interaction sites for the water model you want to use.  OpenMM’s
Modeller class can fix problems such as these.

To use it, create a :class:`Modeller` object, providing the initial :class:`Topology` and atom
positions.  You then can invoke various modelling functions on it.  Each one
modifies the system in some way, creating a new :class:`Topology` and list of positions.
When you are all done, you can retrieve them from the :class:`Modeller` and use them as
the starting point for your simulation:

.. samepage::
    ::

        ...
        pdb = PDBFile('input.pdb')
        modeller = Modeller(pdb.topology, pdb.positions)
        # ... Call some modelling functions here ...
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

    .. caption::

        :autonumber:`Example,Modeller outline`

Now let’s consider the particular functions you can call.

Adding Hydrogens
****************

Call the :meth:`addHydrogens` function to add missing hydrogen atoms:
::

    modeller.addHydrogens(forcefield)

The force field is needed to determine the positions for the hydrogen atoms.  If
the system already contains some hydrogens but is missing others, that is fine.
The Modeller will recognize the existing ones and figure out which ones need to
be added.

Some residues can exist in different protonation states depending on the pH and
on details of the local environment.  By default it assumes pH 7, but you can
specify a different value:
::

    modeller.addHydrogens(forcefield, pH=5.0)

For each residue, it selects the protonation state that is most common at the
specified pH.  In the case of Cysteine residues, it also checks whether the
residue participates in a disulfide bond when selecting the state to use.
Histidine has two different protonation states that are equally likely at
neutral pH.  It therefore selects which one to use based on which will form a
better hydrogen bond.

If you want more control, it is possible to specify exactly which protonation
state to use for particular residues.  For details, consult the API
documentation for the Modeller class.

By default, :class:`Modeller` loads information about hydrogens in standard
amino acids and nucleic acids.  You can call :meth:`loadHydrogenDefinitions` to
load definitions for other types of molecules.  In particular, if your system
contains carbohydrates that you plan to simulate with the GLYCAM force field, call
::

    Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')

All subsequent calls to :meth:`addHydrogens` will make use of the newly loaded
definitions.

Adding Solvent
**************

Call :meth:`addSolvent` to create a box of solvent (water and ions) around the model:
::

    modeller.addSolvent(forcefield)

This constructs a box of water around the solute, ensuring that no water
molecule comes closer to any solute atom than the sum of their van der Waals
radii.  It also determines the charge of the solute, and adds enough positive or
negative ions to make the system neutral.

When called as shown above, :meth:`addSolvent` expects that periodic box dimensions were
specified in the PDB file, and it uses them as the size for the water box.  If
your PDB file does not specify a box size, or if you want to use a different
size, you can specify one:
::

    modeller.addSolvent(forcefield, boxSize=Vec3(5.0, 3.5, 3.5)*nanometers)

This requests a 5 nm by 3.5 nm by 3.5 nm box.  For a non-rectangular box, you
can specify the three box vectors defining the unit cell:
::

    modeller.addSolvent(forcefield, boxVectors=(avec, bvec, cvec))

Another option is to specify a padding distance:
::

    modeller.addSolvent(forcefield, padding=1.0*nanometers)

This determines the largest size of the solute along any axis (x, y, or z).  It
then creates a cubic box of width (solute size)+2*(padding).  The above line
guarantees that no part of the solute comes closer than 1 nm to any edge of the
box.

Finally, you can specify the exact number of solvent molecules (including both
water and ions) to add.  This is useful when you want to solvate several different
conformations of the same molecule while guaranteeing they all have the same
amount of solvent:
::

    modeller.addSolvent(forcefield, numAdded=5000)

By default, :meth:`addSolvent` creates TIP3P water molecules, but it also supports other
water models:
::

    modeller.addSolvent(forcefield, model='tip5p')

Allowed values for the :code:`model` option are ``'tip3p'``, ``'tip3pfb'``, ``'spce'``,
``'tip4pew'``, ``'tip4pfb'``, and ``'tip5p'``.  Be sure to include the single quotes
around the value.

Another option is to add extra ion pairs to give a desired total ionic strength.
For example:
::

    modeller.addSolvent(forcefield, ionicStrength=0.1*molar)

This solvates the system with a salt solution whose ionic strength is 0.1 molar.
Note that when computing the ionic strength, it does *not* consider the ions
that were added to neutralize the solute.  It assumes those are bound to the
solute and do not contribute to the bulk ionic strength.

By default, Na\ :sup:`+` and Cl\ :sup:`-` ions are used, but you can specify
different ones using the :code:`positiveIon` and :code:`negativeIon`
options.  For example, this creates a potassium chloride solution:
::

    modeller.addSolvent(forcefield, ionicStrength=0.1*molar, positiveIon='K+')

Allowed values for :code:`positiveIon` are ``'Cs+'``, ``'K+'``, ``'Li+'``, ``'Na+'``, and
``'Rb+'``.  Allowed values for :code:`negativeIon` are ``'Cl-'``, ``'Br-'``, ``'F-'``, and
``'I-'``.  Be sure to include the single quotes around the value.  Also be aware
some force fields do not include parameters for all of these ion types, so you
need to use types that are supported by your chosen force field.

Adding a Membrane
*****************

If you want to simulate a membrane protein, you may need to create a membrane as
well.  You can do this by calling :meth:`addMembrane`.  Call it *instead* of
:meth:`addSolvent`, not in addition to it.  This one method adds the membrane,
solvent, and ions all at once, making sure the lipid head groups are properly
solvated.  For example, this creates a POPC membrane, ensuring at least 1 nm of
padding on all sides:
::

    modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1*nanometer)

The membrane is added in the XY plane, and the existing protein is assumed to already be oriented
and positioned correctly.  When possible, it is recommended to start with a model
from the `Orientations of Proteins in Membranes`_ (OPM) database.  Otherwise, it
is up to you to select the protein position yourself.

Because this method also adds solvent, it takes many of the same arguments as
:meth:`addSolvent`.  See the API documentation for details.

.. _`Orientations of Proteins in Membranes`: http://opm.phar.umich.edu

.. _adding-or-removing-extra-particles:

Adding or Removing Extra Particles
**********************************

“Extra particles” are particles that do not represent ordinary atoms.  This
includes the virtual interaction sites used in many water models, Drude
particles, etc.  If you are using a force field that involves extra particles,
you must add them to the :class:`Topology`.  To do this, call:
::

    modeller.addExtraParticles(forcefield)

This looks at the force field to determine what extra particles are needed, then
modifies each residue to include them.  This function can remove extra particles
as well as adding them.

Removing Water
**************

Call deleteWater to remove all water molecules from the system:
::

    modeller.deleteWater()

This is useful, for example, if you want to simulate it with implicit solvent.
Be aware, though, that this only removes water molecules, not ions or other
small molecules that might be considered “solvent”.

.. _saving-the-results:

Saving The Results
******************

Once you have finished editing your model, you can immediately use the resulting
:class:`Topology` object and atom positions as the input to a :class:`Simulation`.  If you plan to
simulate it many times, though, it is usually better to save the result to a new
PDB file, then use that as the input for the simulations.  This avoids the cost
of repeating the modelling operations at the start of every simulation, and also
ensures that all your simulations are really starting from exactly the same
structure.

The following example loads a PDB file, adds missing hydrogens, builds a solvent
box around it, performs an energy minimization, and saves the result to a new
PDB file.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *

        print('Loading...')
        pdb = PDBFile('input.pdb')
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        print('Adding hydrogens...')
        modeller.addHydrogens(forcefield)
        print('Adding solvent...')
        modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)
        print('Minimizing...')
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=100)
        print('Saving...')
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
        print('Done')

    .. caption::

        :autonumber:`Example,Modeller complete`

.. default-domain:: py

.. _the-openmm-application-layer-introduction:

Getting Started
###############

Introduction
************

The first thing to understand about the OpenMM “application layer” is that it is
not exactly an application in the traditional sense: there is no program called
“OpenMM” that you run.  Rather, it is a collection of libraries written in the
Python programming language.  Those libraries can easily be chained together to
create Python programs that run simulations.  But don’t worry!  You don’t need
to know anything about Python programming (or programming at all) to use it.
Nearly all molecular simulation applications ask you to write some sort of
“script” that specifies the details of the simulation to run.  With OpenMM, that
script happens to be written in Python.  But it is no harder to write than those
for most other applications, and this guide will teach you everything you need
to know.  There is even a graphical interface that can write the script for you
based on a simple set of options (see Section :numref:`the-script-builder-application`),
so you never need to type a single line of code!

On the other hand, if you don’t mind doing a little programming, this approach
gives you enormous power and flexibility.  Your script has complete access to
the entire OpenMM application programming interface (API), as well as the full
power of the Python language and libraries.  You have complete control over
every detail of the simulation, from defining the molecular system to analyzing
the results.


.. _installing-openmm:

Installing OpenMM
*****************

OpenMM is installed using the Conda package manager (https://docs.conda.io).
Conda is included as part of the Anaconda Python distribution, which you can
download from https://docs.continuum.io/anaconda/install.  This is a Python
distribution specifically designed for scientific applications, with many of the
most popular mathematical and scientific packages preinstalled.  Alternatively
you can use Miniconda (available from https://docs.conda.io/en/latest/miniconda.html),
which includes only Python itself, plus the Conda package manager.  That offers
a much smaller initial download, with the ability to then install only the
packages you want.

(A third option is to compile OpenMM from source.  This provides more flexibility,
but it is much more work, and there is rarely a need for anyone but advanced users
to compile from source.  Detailed instruction are in Chapter :numref:`compiling-openmm-from-source-code`.)

\1. Begin by installing the most recent 64 bit, Python 3.x version of either
Anaconda or Miniconda.

\2. (Optional) If you want to run OpenMM on a GPU, make sure you have installed
modern drivers from your vendor.

  * If you have an Nvidia GPU, download the latest drivers from
    https://www.nvidia.com/Download/index.aspx. CUDA itself will be provided by
    the :code:`cudatoolkit` package when you install :code:`openmm` in the next steps.
  * If you have an AMD GPU and are using Linux or Windows, download the latest
    version of the drivers from https://support.amd.com.  On OS X, OpenCL
    is included with the operating system and is supported on OS X 10.10.3 or
    later.

3. Open a command line terminal and type the following command
::

    conda install -c conda-forge openmm

With recent :code:`conda` versions (v4.8.4+), this will install a version of
OpenMM compiled with the latest version of CUDA supported by your drivers.
Alternatively you can request a version that is compiled for a specific CUDA
version with the command
::

    conda install -c conda-forge openmm cudatoolkit=10.0

where :code:`10.0` should be replaced with the particular CUDA version
you want to target.  We build packages for CUDA 9.2 and above on Linux,
and CUDA 10.0 and above on Windows.  Because different CUDA releases are
not binary compatible with each other, OpenMM can only work with the particular
CUDA version it was compiled with.

.. note::

    Prior to v7.5, conda packages for OpenMM where distributed through the
    :code:`omnia` channel (https://anaconda.org/omnia). Starting with v7.5,
    OpenMM will use the :code:`conda-forge` channel. Check the documentation
    for previous versions in case you want to install older packages.


4. Verify your installation by typing the following command:
::

    python -m openmm.testInstallation

This command confirms that OpenMM is installed, checks whether GPU acceleration
is available (via the OpenCL and/or CUDA platforms), and verifies that all
platforms produce consistent results.

.. default-domain:: py

.. _advanced-simulation-examples:

Advanced Simulation Examples
############################

In the previous chapter, we looked at some basic scripts for running simulations
and saw lots of ways to customize them.  If that is all you want to do—run
straightforward molecular simulations—you already know everything you need to
know.  Just use the example scripts and customize them in the ways described in
Section :numref:`simulation-parameters`.

OpenMM can do far more than that.  Your script has the full OpenMM API at its
disposal, along with all the power of the Python language and libraries.  In
this chapter, we will consider some examples that illustrate more advanced
techniques.  Remember that these are still only examples; it would be impossible
to give an exhaustive list of everything OpenMM can do.  Hopefully they will
give you a sense of what is possible, and inspire you to experiment further on
your own.

Starting in this section, we will assume some knowledge of programming, as well
as familiarity with the OpenMM API.  Consult this User's Guide and the OpenMM
API documentation if you are uncertain about how something works. You can also
usethe Python :code:`help` command.  For example,
::

    help(Simulation)

will print detailed documentation on the :class:`Simulation` class.

Simulated Annealing
*******************

Here is a very simple example of how to do simulated annealing.  The following
lines linearly reduce the temperature from 300 K to 0 K in 100 increments,
executing 1000 time steps at each temperature:

.. samepage::
    ::

        ...
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        for i in range(100):
            integrator.setTemperature(3*(100-i)*kelvin)
            simulation.step(1000)

    .. caption::

        :autonumber:`Example,simulated annealing`

This code needs very little explanation.  The loop is executed 100 times.  Each
time through, it adjusts the temperature of the :class:`LangevinMiddleIntegrator` and then
calls :code:`step(1000)` to take 1000 time steps.

Applying an External Force to Particles: a Spherical Container
**************************************************************

In this example, we will simulate a non-periodic system contained inside a
spherical container with radius 2 nm.  We implement the container by applying a
harmonic potential to every particle:

.. math::
    E(r) = \begin{cases}
           0          & r\le2\\
           100(r-2)^2 & r>2
           \end{cases}

where *r* is the distance of the particle from the origin, measured in nm.
We can easily do this using OpenMM’s :class:`CustomExternalForce` class.  This class
applies a force to some or all of the particles in the system, where the energy
is an arbitrary function of each particle’s (\ *x*\ , *y*\ , *z*\ )
coordinates.  Here is the code to do it:

.. samepage::
    ::

        ...
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                nonbondedCutoff=1*nanometer, constraints=None)
        force = CustomExternalForce('100*max(0, r-2)^2; r=sqrt(x*x+y*y+z*z)')
        system.addForce(force)
        for i in range(system.getNumParticles()):
            force.addParticle(i, [])
        integrator = LangevinMiddleIntegrator(300*kelvin, 91/picosecond, 0.004*picoseconds)
        ...

    .. caption::

        :autonumber:`Example,spherical container`

The first thing it does is create a :class:`CustomExternalForce` object and add it to the
:class:`System`.  The argument to :class:`CustomExternalForce` is a mathematical expression
specifying the potential energy of each particle.  This can be any function of *x*\ ,
*y*\ , and *z* you want.  It also can depend on global or per-particle
parameters.  A wide variety of restraints, steering forces, shearing forces,
etc. can be implemented with this method.

Next it must specify which particles to apply the force to.  In this case, we
want it to affect every particle in the system, so we loop over them and call
:meth:`addParticle` once for each one.  The two arguments are the index of
the particle to affect, and the list of per-particle parameter values (an empty
list in this case).  If we had per-particle parameters, such as to make the
force stronger for some particles than for others, this is where we would
specify them.

Notice that we do all of this immediately after creating the :class:`System`.  That is
not an arbitrary choice.

.. warning::

    If you add new forces to a :class:`System`, you must do so before creating the :class:`Simulation`.
    Once you create a :class:`Simulation`, modifying the :class:`System` will have no effect on that :class:`Simulation`.

Extracting and Reporting Forces (and other data)
************************************************

OpenMM provides reporters for three output formats: `PDB <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_, `PDBx/mmCIF <https://mmcif.wwpdb.org/pdbx-mmcif-home-page.html>`_ and `DCD <https://www.ks.uiuc.edu/Research/namd/2.6/ug/node13.html>`_.
All of those formats store only positions, not velocities, forces, or other data.  In this
section, we create a new reporter that outputs forces.  This illustrates two
important things: how to write a reporter, and how to query the simulation for
forces or other data.

Here is the definition of the :class:`ForceReporter` class:

.. samepage::
    ::

        class ForceReporter(object):
            def __init__(self, file, reportInterval):
                self._out = open(file, 'w')
                self._reportInterval = reportInterval

            def __del__(self):
                self._out.close()

            def describeNextReport(self, simulation):
                steps = self._reportInterval - simulation.currentStep%self._reportInterval
                return (steps, False, False, True, False, None)

            def report(self, simulation, state):
                forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
                for f in forces:
                    self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))

    .. caption::

        :autonumber:`Example,ForceReporter`

The constructor and destructor are straightforward.  The arguments to the
constructor are the output filename and the interval (in time steps) at which it
should generate reports.  It opens the output file for writing and records the
reporting interval.  The destructor closes the file.

We then have two methods that every reporter must implement:
:meth:`describeNextReport()` and :meth:`report()`.  A Simulation object
periodically calls :meth:`describeNextReport()` on each of its reporters to
find out when that reporter will next generate a report, and what information
will be needed to generate it.  The return value should be a six element tuple,
whose elements are as follows:

* The number of time steps until the next report.  We calculate this as
  *(report interval)*\ -\ *(current step)*\ %\ *(report interval)*\ .  For
  example, if we want a report every 100 steps and the simulation is currently on
  step 530, we will return 100-(530%100) = 70.
* Whether the next report will need particle positions.
* Whether the next report will need particle velocities.
* Whether the next report will need forces.
* Whether the next report will need energies.
* Whether the positions should be wrapped to the periodic box.  If None, it will
  automatically decide whether to wrap positions based on whether the System uses
  periodic boundary conditions.


When the time comes for the next scheduled report, the :class:`Simulation` calls
:meth:`report()` to generate the report.  The arguments are the :class:`Simulation`
object, and a :class:`State` that is guaranteed to contain all the information that was
requested by :meth:`describeNextReport()`\ .  A State object contains a
snapshot of information about the simulation, such as forces or particle
positions.  We call :meth:`getForces()` to retrieve the forces and convert
them to the units we want to output (kJ/mole/nm).  Then we loop over each value
and write it to the file.  To keep the example simple, we just print the values
in text format, one line per particle.  In a real program, you might choose a
different output format.

Now that we have defined this class, we can use it exactly like any other
reporter.  For example,
::

    simulation.reporters.append(ForceReporter('forces.txt', 100))

will output forces to a file called “forces.txt” every 100 time steps.

Computing Energies
******************

This example illustrates a different sort of analysis.  Instead of running a
simulation, assume we have already identified a set of structures we are
interested in.  These structures are saved in a set of PDB files.  We want to
loop over all the files in a directory, load them in one at a time, and compute
the potential energy of each one.  Assume we have already created our :class:`System` and
:class:`Simulation`.  The following lines perform the analysis:

.. samepage::
    ::

        import os
        for file in os.listdir('structures'):
            pdb = PDBFile(os.path.join('structures', file))
            simulation.context.setPositions(pdb.positions)
            state = simulation.context.getState(getEnergy=True)
            print(file, state.getPotentialEnergy())

    .. caption::

        :autonumber:`Example,computing energies`

We use Python’s :code:`listdir()` function to list all the files in the
directory.  We create a :class:`PDBFile` object for each one and call
:meth:`setPositions()` on the Context to specify the particle positions loaded
from the PDB file.  We then compute the energy by calling :meth:`getState()`
with the option :code:`getEnergy=True`\ , and print it to the console along
with the name of the file.

.. default-domain:: py

.. _creating-force-fields:

Creating Force Fields
#####################

OpenMM uses a simple XML file format to describe force fields.  It includes many
common force fields, but you can also create your own.  A force field can use
all the standard OpenMM force classes, as well as the very flexible custom force
classes.  You can even extend the ForceField class to add support for completely
new forces, such as ones defined in plugins.  This makes it a powerful tool for
force field development.

Basic Concepts
**************

Let’s start by considering how OpenMM defines a force field.  There are a small
number of basic concepts to understand.

Atom Types and Atom Classes
===========================

Force field parameters are assigned to atoms based on their “atom types”.  Atom
types should be the most specific identification of an atom that will ever be
needed.  Two atoms should have the same type only if the force field will always
treat them identically in every way.

Multiple atom types can be grouped together into “atom classes”.  In general,
two types should be in the same class if the force field usually (but not
necessarily always) treats them identically.  For example, the :math:`\alpha`\ -carbon of an
alanine residue will probably have a different atom type than the :math:`\alpha`\ -carbon of a
leucine residue, but both of them will probably have the same atom class.

All force field parameters can be specified either by atom type or atom class.
Classes exist as a convenience to make force field definitions more compact.  If
necessary, you could define everything in terms of atom types, but when many
types all share the same parameters, it is convenient to only have to specify
them once.

Residue Templates
=================

Types are assigned to atoms by matching residues to templates.  A template
specifies a list of atoms, the type of each one, and the bonds between them.
For each residue in the PDB file, the force field searches its list of templates
for one that has an identical set of atoms with identical bonds between them.
When matching templates, neither the order of the atoms nor their names matter;
it only cares about their elements and the set of bonds between them.  (The PDB
file reader does care about names, of course, since it needs to figure out which
atom each line of the file corresponds to.)

Forces
======

Once a force field has defined its atom types and residue templates, it must
define its force field parameters.  This generally involves one block of XML for
each Force object that will be added to the System.  The details are different
for each Force, but it generally consists of a set of rules for adding
interactions based on bonds and atom types or classes.  For example, when adding
a HarmonicBondForce, the force field will loop over every pair of bonded atoms,
check their types and classes, and see if they match any of its rules.  If so,
it will call :code:`addBond()` on the HarmonicBondForce.  If none of them
match, it simply ignores that pair and continues.

Writing the XML File
********************

The root element of the XML file must be a :code:`<ForceField>` tag:

.. code-block:: xml

    <ForceField>
    ...
    </ForceField>

The :code:`<ForceField>` tag contains the following children:

* An :code:`<AtomTypes>` tag containing the atom type definitions
* A :code:`<Residues>` tag containing the residue template definitions
* Optionally a :code:`<Patches>` tag containing patch definitions
* Zero or more tags defining specific forces


The order of these tags does not matter.  They are described in detail below.

<AtomTypes>
===========

The atom type definitions look like this:

.. code-block:: xml

    <AtomTypes>
     <Type name="0" class="N" element="N" mass="14.00672"/>
     <Type name="1" class="H" element="H" mass="1.007947"/>
     <Type name="2" class="CT" element="C" mass="12.01078"/>
     ...
    </AtomTypes>

There is one :code:`<Type>` tag for each atom type.  It specifies the name
of the type, the name of the class it belongs to, the symbol for its element,
and its mass in amu.  The names are arbitrary strings: they need not be numbers,
as in this example.  The only requirement is that all types have unique names.
The classes are also arbitrary strings, and in general will not be unique.  Two
types belong to the same class if they list the same value for the
:code:`class` attribute.

<Residues>
==========

The residue template definitions look like this:

.. code-block:: xml

    <Residues>
     <Residue name="ACE">
      <Atom name="HH31" type="710"/>
      <Atom name="CH3" type="711"/>
      <Atom name="HH32" type="710"/>
      <Atom name="HH33" type="710"/>
      <Atom name="C" type="712"/>
      <Atom name="O" type="713"/>
      <Bond atomName1="HH31" atomName2="CH3"/>
      <Bond atomName1="CH3" atomName2="HH32"/>
      <Bond atomName1="CH3" atomName2="HH33"/>
      <Bond atomName1="CH3" atomName2="C"/>
      <Bond atomName1="C" atomName2="O"/>
      <ExternalBond atomName="C"/>
     </Residue>
     <Residue name="ALA">
      ...
     </Residue>
     ...
    </Residues>

There is one :code:`<Residue>` tag for each residue template.  That in turn
contains the following tags:

* An :code:`<Atom>` tag for each atom in the residue.  This specifies the
  name of the atom and its atom type.
* A :code:`<Bond>` tag for each pair of atoms that are bonded to each
  other.  The :code:`atomName1` and :code:`atomName2` attributes are the names
  of the two bonded atoms.  (Some older force fields use the alternate tags
  :code:`to` and :code:`from` to specify the atoms by index instead of name.
  This is still supported for backward compatibility, but specifying atoms by
  name is recommended, since it makes the residue definition much easier to
  understand.)
* An :code:`<ExternalBond>` tag for each atom that will be bonded to an
  atom of a different residue.  :code:`atomName` is the name of the atom.
  (Alternatively, the deprecated :code:`from` tag may indicate the atom by
  index instead of name.)


The :code:`<Residue>` tag may also contain :code:`<VirtualSite>` tags,
as in the following example:

.. code-block:: xml

    <Residue name="HOH">
     <Atom name="O" type="tip4pew-O"/>
     <Atom name="H1" type="tip4pew-H"/>
     <Atom name="H2" type="tip4pew-H"/>
     <Atom name="M" type="tip4pew-M"/>
     <VirtualSite type="average3" siteName="M" atomName1="O" atomName2="H1" atomName3="H2"
         weight1="0.786646558" weight2="0.106676721" weight3="0.106676721"/>
     <Bond atomName1="O" atomName2="H1"/>
     <Bond atomName1="O" atomName2="H2"/>
    </Residue>

Each :code:`<VirtualSite>` tag indicates an atom in the residue that should
be represented with a virtual site.  The :code:`type` attribute may equal
:code:`"average2"`\ , :code:`"average3"`\ , :code:`"outOfPlane"`\ , or
:code:`"localCoords"`\ , which correspond to the TwoParticleAverageSite, ThreeParticleAverageSite,
OutOfPlaneSite, and LocalCoordinatesSite classes respectively.  The :code:`siteName`
attribute gives the name of the atom to represent with a virtual site.  The atoms
it is calculated based on are specified by :code:`atomName1`\ , :code:`atomName2`\ , etc.
(Some old force fields use the deprecated tags :code:`index`, :code:`atom1`,
:code:`atom2`, etc. to refer to them by index instead of name.)

The remaining attributes are specific to the virtual site class, and specify the
parameters for calculating the site position.  For a TwoParticleAverageSite,
they are :code:`weight1` and :code:`weight2`\ .  For a
ThreeParticleAverageSite, they are :code:`weight1`\ , :code:`weight2`\ , and
\ :code:`weight3`\ . For an OutOfPlaneSite, they are :code:`weight12`\ ,
:code:`weight13`\ , and :code:`weightCross`\ . For a LocalCoordinatesSite, they
are :code:`p1`\ , :code:`p2`\ , and :code:`p3` (giving the x, y, and z coordinates
of the site position in the local coordinate system), and :code:`wo1`\ ,
:code:`wx1`\ , :code:`wy1`\ , :code:`wo2`\ , :code:`wx2`\ , :code:`wy2`\ , ...
(giving the weights for computing the origin, x axis, and y axis).

<Patches>
=========

A "patch" is a set of rules for modifying a residue template (or possibly multiple
templates at once).  For example a terminal amino acid is slightly different from
one in the middle of a chain.  A force field could of course define multiple
templates for each amino acid (standard, N-terminal, C-terminal, and monomer),
but since the modifications are the same for nearly all amino acids, it is simpler
to include only the "standard" templates, along with a set of patches for
modifying terminal residues.

Here is an example of a patch definition:

.. code-block:: xml

    <Patches>
     <Patch name="NTER">
      <RemoveAtom name="H"/>
      <RemoveBond atomName1="N" atomName2="H"/>
      <AddAtom name="H1" type="H"/>
      <AddAtom name="H2" type="H"/>
      <AddAtom name="H3" type="H"/>
      <AddBond atomName1="N" atomName2="H1"/>
      <AddBond atomName1="N" atomName2="H2"/>
      <AddBond atomName1="N" atomName2="H3"/>
      <RemoveExternalBond atomName="N"/>
      <ChangeAtom name="N" type="N3"/>
     </Patch>
    </Patches>

There is one :code:`<Patch>` tag for each patch definition.  That in turn may
contain any of the following tags:

 * An :code:`<AddAtom>` tag indicates that an atom should be added to the
   template.  It specifies the name of the atom and its atom type.
 * A :code:`<ChangeAtom>` tag indicates that the type of an atom already present
   in the template should be altered.  It specifies the name of the atom and its
   new atom type.
 * A :code:`<RemoveAtom>` tag indicates that an atom should be removed from the
   template.  It specifies the name of the atom to remove.
 * An :code:`<AddBond>` tag indicates that a bond should be added to the
   template.  It specifies the names of the two bonded atoms.
 * A :code:`<RemoveBond>` tag indicates that a bond already present in the
   template should be removed.  It specifies the names of the two bonded atoms.
 * An :code:`<AddExternalBond>` tag indicates that a new external bond should be
   added to the template.  It specifies the name of the bonded atom.
 * A :code:`<RemoveExternalBond>` tag indicates that an external bond aleady
   present in the template should be removed.  It specifies the name of the
   bonded atom.

In addition to defining the patches, you also must identify which residue
templates each patch can be applied to.  This can be done in two ways.  The more
common one is to have each template identify the patches that can be applied to
it.  This is done with an :code:`<AllowPatch>` tag:

.. code-block:: xml

    <Residue name="ALA">
     <AllowPatch name="CTER"/>
     <AllowPatch name="NTER"/>
     ...
    </Residue>

Alternatively, the patch can indicate which residues it may be applied to.  This
is done with an :code:`<ApplyToResidue>` tag:

.. code-block:: xml

    <Patch name="NTER">
     <ApplyToResidue name="ALA"/>
     <ApplyToResidue name="ARG"/>
     ...
    </Patch>

A patch can alter multiple templates at once.  This is useful for creating bonds
between molecules, and allows the atom types in one residue to depend on the
identity of the other residue it is bonded to.  To create a multi-residue patch,
added a :code:`residues` attribute to the :code:`<Patch>` tag specifying how many
residues that patch covers.  Then whenever you refer to an atom, prefix its name
with the index of the residue it belongs to:

.. code-block:: xml

  <Patch name="Disulfide" residues="2">
    <RemoveAtom name="1:HG"/>
    <RemoveAtom name="2:HG"/>
    <AddBond atomName1="1:SG" atomName2="2:SG"/>
    <ApplyToResidue name="1:CYS"/>
    <ApplyToResidue name="2:CYS"/>
  </Patch>

In this example, the patch modifies two residues of the same type, but that need
not always be true.  Each :code:`<ApplyToResidue>` tag therefore indicates which
one of the residue templates it modifies may be of the specified type.  Similarly,
if a residue template includes an :code:`<AcceptPatch>` tag for a multi-residue
patch, it must specify the name of the patch, followed by the index of the residue
within that patch:

.. code-block:: xml

    <AllowPatch name="Disulfide:1"/>


Missing residue templates
=========================

.. CAUTION::
   These features are experimental, and their API is subject to change.

You can use the :meth:`getUnmatchedResidues()` method to get a list of residues
in the provided :code:`topology` object that do not currently have a matching
residue template defined in the :class:`ForceField`.
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

This is useful for identifying issues with prepared systems, debugging issues
with residue template definitions, or identifying which additional residues need
to be parameterized.

As a convenience for parameterizing new residues, you can also get a list of
residues and empty residue templates using :meth:`generateTemplatesForUnmatchedResidues`
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(topology)
    # Se the atom types
    for template in templates:
        for atom in template.atoms:
            atom.type = ... # set the atom types here
        # Register the template with the forcefield.
        forcefield.registerResidueTemplate(template)

If you find that templates seem to be incorrectly matched, another useful
function :meth:`getMatchingTemplates()` can help you identify which templates
are being matched:
::

    pdb = PDBFile('input.pdb')
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    templates = forcefield.getMatchingTemplates(topology)
    for (residue, template) in zip(pdb.topology.residues(), templates):
        print("Residue %d %s matched template %s" % (residue.id, residue.name, template.name))

<HarmonicBondForce>
===================

To add a HarmonicBondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <HarmonicBondForce>
     <Bond class1="C" class2="C" length="0.1525" k="259408.0"/>
     <Bond class1="C" class2="CA" length="0.1409" k="392459.2"/>
     <Bond class1="C" class2="CB" length="0.1419" k="374049.6"/>
     ...
    </HarmonicBondForce>

Every :code:`<Bond>` tag defines a rule for creating harmonic bond
interactions between atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1` and :code:`type2`\ ) or by class
(using the attributes :code:`class1` and :code:`class2`\ ).  For every
pair of bonded atoms, the force field searches for a rule whose atom types or
atom classes match the two atoms.  If it finds one, it calls
:code:`addBond()` on the HarmonicBondForce with the specified parameters.
Otherwise, it ignores that pair and continues.  :code:`length` is the
equilibrium bond length in nm, and :code:`k` is the spring constant in
kJ/mol/nm\ :sup:`2`\ .

<HarmonicAngleForce>
====================

To add a HarmonicAngleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <HarmonicAngleForce>
     <Angle class1="C" class2="C" class3="O" angle="2.094" k="669.44"/>
     <Angle class1="C" class2="C" class3="OH" angle="2.094" k="669.44"/>
     <Angle class1="CA" class2="C" class3="CA" angle="2.094" k="527.184"/>
     ...
    </HarmonicAngleForce>

Every :code:`<Angle>` tag defines a rule for creating harmonic angle
interactions between triplets of atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of three atoms in the system where the first is
bonded to the second, and the second to the third.  For each one, it searches
for a rule whose atom types or atom classes match the three atoms.  If it finds
one, it calls :code:`addAngle()` on the HarmonicAngleForce with the
specified parameters.  Otherwise, it ignores that set and continues.
:code:`angle` is the equilibrium angle in radians, and :code:`k` is the
spring constant in kJ/mol/radian\ :sup:`2`\ .

<PeriodicTorsionForce>
======================

To add a PeriodicTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <PeriodicTorsionForce>
     <Proper class1="HC" class2="CT" class3="CT" class4="CT" periodicity1="3" phase1="0.0"
         k1="0.66944"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="HC" periodicity1="3" phase1="0.0"
         k1="0.6276"/>
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" periodicity1="2"
         phase1="3.14159265359" k1="4.6024"/>
     <Improper class1="N" class2="C" class3="CT" class4="H" periodicity1="2"
         phase1="3.14159265359" k1="4.6024"/>
     ...
    </PeriodicTorsionForce>

Every child tag defines a rule for creating periodic torsion interactions
between sets of four atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by class
(using the attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
PeriodicTorsionForce with the specified parameters.  Otherwise, it ignores that
set and continues.  :code:`periodicity1` is the periodicity of the torsion,
\ :code:`phase1` is the phase offset in radians, and :code:`k1` is the
force constant in kJ/mol.

Each torsion definition can specify multiple periodic torsion terms to add to
its atoms.  To add a second one, just add three more attributes:
:code:`periodicity2`\ , :code:`phase2`\ , and :code:`k2`\ .  You can have as
many terms as you want.  Here is an example of a rule that adds three torsion
terms to its atoms:

.. code-block:: xml

    <Proper class1="CT" class2="CT" class3="CT" class4="CT"
        periodicity1="3" phase1="0.0" k1="0.75312"
        periodicity2="2" phase2="3.14159265359" k2="1.046"
        periodicity3="1" phase3="3.14159265359" k3="0.8368"/>

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" periodicity1="3" phase1="0.0" k1="1.046"/>

The :code:`<PeriodicTorsionForce>` tag also supports an optional
:code:`ordering` attribute to provide better compatibility with the way
impropers are assigned in different simulation packages:

 * :code:`ordering="default"` specifies the default behavior if the attribute
   is omitted.
 * :code:`ordering="amber"` produces behavior that replicates the behavior of
   AmberTools LEaP
 * :code:`ordering="charmm"` produces behavior more consistent with CHARMM
 * :code:`ordering="smirnoff"` allows multiple impropers to be added using
   exact matching to replicate the beheavior of `SMIRNOFF <https://open-forcefield-toolkit.readthedocs.io/en/latest/users/smirnoff.html>`_
   impropers

Different :code:`<PeriodicTorsionForce>` tags can specify different :code:`ordering`
values to be used for the sub-elements appearing within their tags.

<RBTorsionForce>
================

To add an RBTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <RBTorsionForce>
     <Proper class1="CT" class2="CT" class3="OS" class4="CT" c0="2.439272" c1="4.807416"
         c2="-0.8368" c3="-6.409888" c4="0" c5="0" />
     <Proper class1="C" class2="N" class3="CT" class4="C" c0="10.46" c1="-3.34720"
         c2="-7.1128" c3="0" c4="0" c5="0" />
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" c0="0.8368" c1="0"
         c2="-2.76144" c3="0" c4="3.3472" c5="0" />
     <Improper class1="N" class2="C" class3="CT" class4="H" c0="29.288" c1="-8.368"
         c2="-20.92" c3="0" c4="0" c5="0" />
     ...
    </RBTorsionForce>

Every child tag defines a rule for creating Ryckaert-Bellemans torsion
interactions between sets of four atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
RBTorsionForce with the specified parameters.  Otherwise, it ignores that set
and continues.  The attributes :code:`c0` through :code:`c5` are the
coefficients of the terms in the Ryckaert-Bellemans force expression.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" c0="2.439272" c1="4.807416"
        c2="-0.8368" c3="-6.409888" c4="0" c5="0" />

<CMAPTorsionForce>
==================

To add a CMAPTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CMAPTorsionForce>
     <Map>
      0.0 0.809 0.951 0.309
      -0.587 -1.0 -0.587 0.309
      0.951 0.809 0.0 -0.809
      -0.951 -0.309 0.587 1.0
     </Map>
     <Torsion map="0" class1="CT" class2="CT" class3="C" class4="N" class5="CT"/>
     <Torsion map="0" class1="N" class2="CT" class3="C" class4="N" class5="CT"/>
     ...
    </CMAPTorsionForce>

Each :code:`<Map>` tag defines an energy correction map.  Its content is the
list of energy values in kJ/mole, listed in the correct order for
CMAPTorsionForce’s :code:`addMap()` method and separated by white space.
See the API documentation for details.  The size of the map is determined from
the number of energy values.

Each :code:`<Torsion>` tag defines a rule for creating CMAP torsion
interactions between sets of five atoms.  The tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of five atoms that are bonded in sequence: 1 to
2, 2 to 3, 3 to 4, and 4 to 5.  For each one, it searches for a rule whose atom
types or atom classes match the five atoms.  If it finds one, it calls
:code:`addTorsion()` on the CMAPTorsionForce with the specified parameters.
Otherwise, it ignores that set and continues.  The first torsion is defined by
the sequence of atoms 1-2-3-4, and the second one by atoms 2-3-4-5.
:code:`map` is the index of the map to use, starting from 0, in the order they
are listed in the file.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of five atoms
where the middle three have classes CT, C, and N respectively:

.. code-block:: xml

    <Torsion map="0" class1="" class2="CT" class3="C" class4="N" class5=""/>

<NonbondedForce>
================

To add a NonbondedForce to the System, include a tag that looks like this:

.. code-block:: xml

    <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
     <Atom type="0" charge="-0.4157" sigma="0.32499" epsilon="0.71128"/>
     <Atom type="1" charge="0.2719" sigma="0.10690" epsilon="0.06568"/>
     <Atom type="2" charge="0.0337" sigma="0.33996" epsilon="0.45772"/>
     ...
    </NonbondedForce>

The :code:`<NonbondedForce>` tag has two attributes
:code:`coulomb14scale` and :code:`lj14scale` that specify the scale
factors between pairs of atoms separated by three bonds.  After setting the
nonbonded parameters for all atoms, the force field calls
:code:`createExceptionsFromBonds()` on the NonbondedForce, passing in these
scale factors as arguments.

Each :code:`<Atom>` tag specifies the nonbonded parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
:code:`charge` is measured in units of the proton charge, :code:`sigma`
is in nm, and :code:`epsilon` is in kJ/mole.

<GBSAOBCForce>
==============

To add a GBSAOBCForce to the System, include a tag that looks like this:

.. code-block:: xml

    <GBSAOBCForce>
     <Atom type="0" charge="-0.4157" radius="0.1706" scale="0.79"/>
     <Atom type="1" charge="0.2719" radius="0.115" scale="0.85"/>
     <Atom type="2" charge="0.0337" radius="0.19" scale="0.72"/>
     ...
    </GBSAOBCForce>

Each :code:`<Atom>` tag specifies the OBC parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
:code:`charge` is measured in units of the proton charge, :code:`radius`
is the GBSA radius in nm, and :code:`scale` is the OBC scaling factor.

<CustomBondForce>
=================

To add a CustomBondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomBondForce energy="scale*k*(r-r0)^2">
     <GlobalParameter name="scale" defaultValue="0.5"/>
     <PerBondParameter name="k"/>
     <PerBondParameter name="r0"/>
     <Bond class1="OW" class2="HW" r0="0.09572" k="462750.4"/>
     <Bond class1="HW" class2="HW" r0="0.15136" k="462750.4"/>
     <Bond class1="C" class2="C" r0="0.1525" k="259408.0"/>
     ...
    </CustomBondForce>

The energy expression for the CustomBondForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each bond as a function of its length *r*\ .  It also may depend on
an arbitrary list of global or per-bond parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerBondParameter>` tag to define a per-bond parameter.

Every :code:`<Bond>` tag defines a rule for creating custom bond
interactions between atoms.  Each tag may identify the atoms either by type
(using the attributes :code:`type1` and :code:`type2`\ ) or by class
(using the attributes :code:`class1` and :code:`class2`\ ).  For every
pair of bonded atoms, the force field searches for a rule whose atom types or
atom classes match the two atoms.  If it finds one, it calls
:code:`addBond()` on the CustomBondForce.  Otherwise, it ignores that pair and
continues.  The remaining attributes are the values to use for the per-bond
parameters.  All per-bond parameters must be specified for every
:code:`<Bond>` tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-bond parameter with the name “k”,
then every :code:`<Bond>` tag must include an attribute called :code:`k`\ .

<CustomAngleForce>
==================

To add a CustomAngleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomAngleForce energy="scale*k*(theta-theta0)^2">
     <GlobalParameter name="scale" defaultValue="0.5"/>
     <PerAngleParameter name="k"/>
     <PerAngleParameter name=" theta0"/>
     <Angle class1="HW" class2="OW" class3="HW" theta0="1.824218" k="836.8"/>
     <Angle class1="HW" class2="HW" class3="OW" theta0="2.229483" k="0.0"/>
     <Angle class1="C" class2="C" class3="O" theta0="2.094395" k="669.44"/>
     ...
    </CustomAngleForce>

The energy expression for the CustomAngleForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each angle as a function of the angle *theta*\ .  It also may depend
on an arbitrary list of global or per-angle parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerAngleParameter>` tag to define a per-angle parameter.

Every :code:`<Angle>` tag defines a rule for creating custom angle
interactions between triplets of atoms.  Each tag may identify the atoms either
by type (using the attributes :code:`type1`\ , :code:`type2`\ , ...) or by
class (using the attributes :code:`class1`\ , :code:`class2`\ , ...).  The
force field identifies every set of three atoms in the system where the first is
bonded to the second, and the second to the third.  For each one, it searches
for a rule whose atom types or atom classes match the three atoms.  If it finds
one, it calls :code:`addAngle()` on the CustomAngleForce.  Otherwise, it
ignores that set and continues. The remaining attributes are the values to use
for the per-angle parameters. All per-angle parameters must be specified for
every :code:`<Angle>` tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-angle parameter with the name “k”,
then every :code:`<Angle>` tag must include an attribute called :code:`k`\ .

<CustomTorsionForce>
====================

To add a CustomTorsionForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomTorsionForce energy="scale*k*(1+cos(per*theta-phase))">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerTorsionParameter name="k"/>
     <PerTorsionParameter name="per"/>
     <PerTorsionParameter name="phase"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="CT" per="3" phase="0.0" k="0.66944"/>
     <Proper class1="HC" class2="CT" class3="CT" class4="HC" per="3" phase="0.0" k="0.6276"/>
     ...
     <Improper class1="N" class2="C" class3="CT" class4="O" per="2" phase="3.14159265359"
         k="4.6024"/>
     <Improper class1="N" class2="C" class3="CT" class4="H" per="2" phase="3.14159265359"
         k="4.6024"/>
     ...
    </CustomTorsionForce>

The energy expression for the CustomTorsionForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each torsion as a function of the angle *theta*\ .  It also may
depend on an arbitrary list of global or per-torsion parameters.  Use a
:code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerTorsionParameter>` tag to define a per-torsion parameter.

Every child tag defines a rule for creating custom torsion interactions between
sets of four atoms.  Each tag may identify the atoms either by type (using the
attributes :code:`type1`\ , :code:`type2`\ , ...) or by class (using the
attributes :code:`class1`\ , :code:`class2`\ , ...).

The force field recognizes two different types of torsions: proper and improper.
A proper torsion involves four atoms that are bonded in sequence: 1 to 2, 2 to
3, and 3 to 4.  An improper torsion involves a central atom and three others
that are bonded to it: atoms 2, 3, and 4 are all bonded to atom 1.  The force
field begins by identifying every set of atoms in the system of each of these
types. For each one, it searches for a rule whose atom types or atom classes
match the four atoms.  If it finds one, it calls :code:`addTorsion()` on the
CustomTorsionForce with the specified parameters.  Otherwise, it ignores that
set and continues. The remaining attributes are the values to use for the per-
torsion parameters.  Every :code:`<Torsion>` tag must include one attribute
for every per-torsion parameter, and the attribute name must match the name of
the parameter.

You can also use wildcards when defining torsions.  To do this, simply leave the
type or class name for an atom empty.  That will cause it to match any atom.
For example, the following definition will match any sequence of atoms where the
second atom has class OS and the third has class P:

.. code-block:: xml

    <Proper class1="" class2="OS" class3="P" class4="" per="3" phase="0.0" k="0.66944"/>

<CustomNonbondedForce>
======================

To add a CustomNonbondedForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomNonbondedForce energy="scale*epsilon1*epsilon2*((sigma1+sigma2)/r)^12" bondCutoff="3">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerParticleParameter name="sigma"/>
     <PerParticleParameter name="epsilon"/>
     <Atom type="0" sigma="0.3249" epsilon="0.7112"/>
     <Atom type="1" sigma="0.1069" epsilon="0.0656"/>
     <Atom type="2" sigma="0.3399" epsilon="0.4577"/>
     ...
    </CustomNonbondedForce>

The energy expression for the CustomNonbondedForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each pairwise interaction as a function of the distance *r*\ .  It
also may depend on an arbitrary list of global or per-particle parameters.  Use
a :code:`<GlobalParameter>` tag to define a global parameter, and a
:code:`<PerParticleParameter>` tag to define a per-particle parameter.

The expression may also depend on computed values, each defined with a
:code:`<ComputedValue>` tag.  The tag should have two attributes, :code:`name`
with the name of the computed value, and :code:`expression` with the expression
used to compute it.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
After setting the nonbonded parameters for all atoms, the force field calls
:code:`createExclusionsFromBonds()` on the CustomNonbondedForce, passing in this
value as its argument.  To avoid creating exclusions, set :code:`bondCutoff` to 0.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomNonbondedForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomGBForce>
===============

To add a CustomGBForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomGBForce>
     <GlobalParameter name="solventDielectric" defaultValue="78.3"/>
     <GlobalParameter name="soluteDielectric" defaultValue="1"/>
     <PerParticleParameter name="charge"/>
     <PerParticleParameter name="radius"/>
     <PerParticleParameter name="scale"/>
     <ComputedValue name="I" type="ParticlePairNoExclusions">
        step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);
        U=r+sr2; C=2*(1/or1-1/L)*step(sr2-r-or1); L=max(or1, D); D=abs(r-sr2); sr2 =
        scale2*or2; or1 = radius1-0.009; or2 = radius2-0.009
     </ComputedValue>
     <ComputedValue name="B" type="SingleParticle">
      1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius); psi=I*or; or=radius-0.009
     </ComputedValue>
     <EnergyTerm type="SingleParticle">
      28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*
              (1/soluteDielectric-1/solventDielectric)*charge^2/B
     </EnergyTerm>
     <EnergyTerm type="ParticlePair">
      -138.935456*(1/soluteDielectric-1/solventDielectric)*charge1*charge2/f;
              f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))
     </EnergyTerm>
     <Atom type="0" charge="-0.4157" radius="0.1706" scale="0.79"/>
     <Atom type="1" charge="0.2719" radius="0.115" scale="0.85"/>
     <Atom type="2" charge="0.0337" radius="0.19" scale="0.72"/>
     ...
    </CustomGBForce>

The above (rather complicated) example defines a generalized Born model that is
equivalent to GBSAOBCForce.  The definition consists of a set of computed values
(defined by :code:`<ComputedValue>` tags) and energy terms (defined by
:code:`<EnergyTerm>` tags), each of which is evaluated according to a
mathematical expression.  See the API documentation for details.

The expressions may depend on an arbitrary list of global or per-atom
parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, and a :code:`<PerAtomParameter>` tag to define a per-atom
parameter.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomGBForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomHbondForce>
=========================

To add a CustomHbondForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomHbondForce particlesPerDonor="3" particlesPerAcceptor="2" bondCutoff="2"
        energy="scale*k*(distance(a1,d1)-r0)^2*(angle(a1,d1,d2)-theta0)^2">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerDonorParameter name="theta0"/>
     <PerAcceptorParameter name="k"/>
     <PerAcceptorParameter name="r0"/>
     <Donor class1="H" class2="N" class3="C" theta0="2.1"/>
     <Acceptor class1="O" class2="C" k="115.0" r0="0.2"/>
     ...
    </CustomHbondForce>

The energy expression for the CustomHbondForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each donor-acceptor interaction as a function of various particle coordinates,
distances, and angles.  See the API documentation for details.  :code:`particlesPerDonor`
specifies the number of particles that make up a donor group, and :code:`particlesPerAcceptor`
specifies the number of particles that make up an acceptor group.

The expression may depend on an arbitrary list of global, per-donor, or
per-acceptor parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, a :code:`<PerDonorParameter>` tag to define a per-donor parameter,
and a :code:`<PerAcceptorParameter>` tag to define a per-acceptor parameter.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
If any atom of a donor is within the specified distance (measured in bonds) of
any atom of an acceptor, an exclusion is added to prevent them from interacting
with each other.  If a donor and an acceptor share any atom in common, that is a
bond distance of 0, so they are always excluded.

Every :code:`<Donor>` or :code:`<Acceptor>` tag defines a rule for creating donor
or acceptor groups.  The number of atoms specified in each one must match the
value of :code:`particlesPerDonor` or :code:`particlesPerAcceptor` specified in the
parent tag. Each tag may identify the atoms either by type (using the attributes
:code:`type1`\ , :code:`type2`\ , ...) or by class (using the attributes
:code:`class1`\ , :code:`class2`\ , ...).  The force field considers every atom
in the system (if the number of atoms is 1), every pair of bonded atoms (if the number
of atoms is 2), or every set of three atoms where the first is bonded to the second
and the second to the third (if the number of atoms is 3).  For each one, it searches
for a rule whose atom types or atom classes match the atoms.  If it finds one,
it calls :code:`addDonor()` or :code:`addAcceptor()` on the CustomHbondForce.
Otherwise, it ignores that set and continues. The remaining attributes are the
values to use for the per-donor and per-acceptor parameters. All parameters must
be specified for every tag, and the attribute name must match the name of the
parameter.  For instance, if there is a per-donor parameter with the name “k”,
then every :code:`<Donor>` tag must include an attribute called :code:`k`\ .

CustomHbondForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

<CustomManyParticleForce>
=========================

To add a CustomManyParticleForce to the System, include a tag that looks like this:

.. code-block:: xml

    <CustomManyParticleForce particlesPerSet="3" permutationMode="UniqueCentralParticle"
        bondCutoff="3" energy="scale*(distance(p1,p2)-r1)*(distance(p1,p3)-r1)">
     <GlobalParameter name="scale" defaultValue="1"/>
     <PerParticleParameter name="r"/>
     <TypeFilter index="0" types="1,2"/>
     <Atom type="0" r="0.31" filterType="0"/>
     <Atom type="1" r="0.25" filterType="0"/>
     <Atom type="2" r="0.33" filterType="1"/>
     ...
    </CustomManyParticleForce>

The energy expression for the CustomManyParticleForce is specified by the
:code:`energy` attribute.  This is a mathematical expression that gives the
energy of each interaction as a function of various particle coordinates,
distances, and angles.  See the API documentation for details.  :code:`particlesPerSet`
specifies the number of particles involved in the interaction and
:code:`permutationMode` specifies the permutation mode.

The expression may depend on an arbitrary list of global or per-atom
parameters.  Use a :code:`<GlobalParameter>` tag to define a global
parameter, and a :code:`<PerAtomParameter>` tag to define a per-atom
parameter.

Exclusions are created automatically based on the :code:`bondCutoff` attribute.
After setting the nonbonded parameters for all atoms, the force field calls
:code:`createExclusionsFromBonds()` on the CustomManyParticleForce, passing in this
value as its argument.  To avoid creating exclusions, set :code:`bondCutoff` to 0.

Type filters may be specified with a :code:`<TypeFilter>` tag.  The :code:`index`
attribute specifies the index of the particle to apply the filter to, and
:code:`types` is a comma separated list of allowed types.

Each :code:`<Atom>` tag specifies the parameters for one atom type
(specified with the :code:`type` attribute) or atom class (specified with
the :code:`class` attribute).  It is fine to mix these two methods, having
some tags specify a type and others specify a class.  However you do it, you
must make sure that a unique set of parameters is defined for every atom type.
In addition, each :code:`<Atom>` tag must include the :code:`filterType`
attribute, which specifies the atom type for use in type filters.
The remaining attributes are the values to use for the per-atom parameters. All
per-atom parameters must be specified for every :code:`<Atom>` tag, and the
attribute name must match the name of the parameter.  For instance, if there is
a per-atom parameter with the name “radius”, then every :code:`<Atom>` tag
must include an attribute called :code:`radius`\ .

CustomManyParticleForce also allows you to define tabulated functions.  See Section
:numref:`tabulated-functions` for details.

Writing Custom Expressions
==========================

The custom forces described in this chapter involve user defined algebraic
expressions.  These expressions are specified as character strings, and may
involve a variety of standard operators and mathematical functions.

The following operators are supported: + (add), - (subtract), * (multiply), /
(divide), and ^ (power).  Parentheses “(“and “)” may be used for grouping.

The following standard functions are supported: sqrt, exp, log, sin, cos, sec,
csc, tan, cot, asin, acos, atan, sinh, cosh, tanh, erf, erfc, min, max, abs,
floor, ceil, step, delta, select. step(x) = 0 if x < 0, 1 otherwise.
delta(x) = 1 if x is 0, 0 otherwise.  select(x,y,z) = z if x = 0, y otherwise.
Some custom forces allow additional functions to be defined from tabulated values.

Numbers may be given in either decimal or exponential form.  All of the
following are valid numbers: 5, -3.1, 1e6, and 3.12e-2.

The variables that may appear in expressions are specified in the API
documentation for each force class.  In addition, an expression may be followed
by definitions for intermediate values that appear in the expression.  A
semicolon “;” is used as a delimiter between value definitions.  For example,
the expression
::

    a^2+a*b+b^2; a=a1+a2; b=b1+b2

is exactly equivalent to
::

    (a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2

The definition of an intermediate value may itself involve other intermediate
values.  All uses of a value must appear *before* that value’s definition.

.. _tabulated-functions:

Tabulated Functions
===================

Some forces, such as CustomNonbondedForce and CustomGBForce, allow you to define
tabulated functions.  To define a function, include a :code:`<Function>` tag inside the
:code:`<CustomNonbondedForce>` or :code:`<CustomGBForce>` tag:

.. code-block:: xml

    <Function name="myfn" type="Continuous1D" min="-5" max="5">
    0.983674857694 -0.980096396266 -0.975743130031 -0.970451936613 -0.964027580076
    -0.956237458128 -0.946806012846 -0.935409070603 -0.921668554406 -0.905148253645
    -0.885351648202 -0.861723159313 -0.833654607012 -0.800499021761 -0.761594155956
    -0.716297870199 -0.664036770268 -0.604367777117 -0.537049566998 -0.46211715726
    -0.379948962255 -0.291312612452 -0.197375320225 -0.099667994625 0.0
    0.099667994625 0.197375320225 0.291312612452 0.379948962255 0.46211715726
    0.537049566998 0.604367777117 0.664036770268 0.716297870199 0.761594155956
    0.800499021761 0.833654607012 0.861723159313 0.885351648202 0.905148253645
    0.921668554406 0.935409070603 0.946806012846 0.956237458128 0.964027580076
    0.970451936613 0.975743130031 0.980096396266 0.983674857694 0.986614298151
    0.989027402201
    </Function>

The tag’s attributes define the name of the function, the type of function, and
the range of values for which it is defined.  The required set of attributed
depends on the function type:

.. tabularcolumns:: |l|L|

============  =======================================================
Type          Required Attributes
============  =======================================================
Continuous1D  min, max
Continuous2D  xmin, ymin, xmax, ymax, xsize, ysize
Continuous3D  xmin, ymin, zmin, xmax, ymax, zmax, xsize, ysize, zsize
Discrete1D
Discrete2D    xsize, ysize
Discrete3D    xsize, ysize, zsize
============  =======================================================


The "min" and "max" attributes define the range of the independent variables for
a continuous function.  The "size" attributes define the size of the table along
each axis.  The tabulated values are listed inside the body of the tag, with
successive values separated by white space.  See the API documentation for more
details.


Residue Template Parameters
===========================

In forces that use an :code:`<Atom>` tag to define parameters for atom types or
classes, there is an alternate mechanism you can also use: defining those
parameter values in the residue template.  This is useful for situations that
come up in certain force fields.  For example, :code:`NonbondedForce` and
:code:`GBSAOBCForce` each have a :code:`charge` attribute.  If you only have to
define the charge of each atom type once, that is more convenient and avoids
potential bugs.  Also, many force fields have a different charge for each atom
type, but Lennard-Jones parameters that are the same for all types in a class.
It would be preferable not to have to repeat those parameter values many times
over.

When writing a residue template, you can add arbitrary additional attributes
to each :code:`<Atom>` tag.  For example, you might include a :code:`charge`
attribute as follows:

.. code-block:: xml

   <Atom name="CA" type="53" charge="0.0381"/>

When writing the tag for a force, you can then include a
:code:`<UseAttributeFromResidue>` tag inside it.  This indicates that a
specified attribute should be taken from the residue template.  Finally, you
simply omit that attribute in the force's own :code:`<Atom>` tags.  For example:

.. code-block:: xml

    <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
     <UseAttributeFromResidue name="charge"/>
     <Atom class="CX" sigma="0.339966950842" epsilon="0.4577296"/>
     ...
    </NonbondedForce>

Notice that the :code:`charge` attribute is missing, and that the parameters
are specified by class, not by type.  This means that sigma and epsilon only
need to be specified once for each class.  The atom charges, which are different
for each type, are taken from the residue template instead.


Including Other Files
=====================

Sometimes it is useful to split a force field definition into multiple files,
but still be able to use the force field by specifying only a single file.  You
can accomplish this with the :code:`<Include>` tag.  For example:

.. code-block:: xml

    <ForceField>
     <Include file="definitions.xml"/>
     ...
    </ForceField>

The :code:`file` attribute gives the path of the file to include.  It may be
relative either to the directory containing the parent XML file (the one with
the :code:`<Include>` tag) or the OpenMM data directory (the one containing
built in force fields).


Using Multiple Files
********************

If multiple XML files are specified when a ForceField is created, their
definitions are combined as follows.

* A file may refer to atom types and classes that it defines, as well as those
  defined in previous files.  It may not refer to ones defined in later files.
  This means that the order in which files are listed when calling the ForceField
  constructor is potentially significant.
* Forces that involve per-atom parameters (such as NonbondedForce or
  GBSAOBCForce) require parameter values to be defined for every atom type.  It
  does not matter which file those types are defined in.  For example, files that
  define explicit water models generally define a small number of atom types, as
  well as nonbonded parameters for those types.  In contrast, files that define
  implicit solvent models do not define any new atom types, but provide parameters
  for all the atom types that were defined in the main force field file.
* For other forces, the files are effectively independent.  For example, if two
  files each include a :code:`<HarmonicBondForce>` tag, bonds will be created
  based on the rules in the first file, and then more bonds will be created based
  on the rules in the second file.  This means you could potentially end up with
  multiple bonds between a single pair of atoms.


Extending ForceField
********************

The ForceField class is designed to be modular and extensible.  This means you
can add support for entirely new force types, such as ones implemented with
plugins.

Adding new force types
======================

For every force class, there is a “generator” class that parses the
corresponding XML tag, then creates Force objects and adds them to the System.
ForceField maintains a map of tag names to generator classes.  When a ForceField
is created, it scans through the XML files, looks up the generator class for
each tag, and asks that class to create a generator object based on it.  Then,
when you call :code:`createSystem()`\ ,  it loops over each of its generators
and asks each one to create its Force object.  Adding a new Force type therefore
is simply a matter of creating a new generator class and adding it to
ForceField’s map.

The generator class must define two methods.  First, it needs a static method
with the following signature to parse the XML tag and create the generator:
::

    @staticmethod
    def parseElement(element, forcefield):

:code:`element` is the XML tag (an xml.etree.ElementTree.Element object) and
:code:`forcefield` is the ForceField being created.  This method should
create a generator and add it to the ForceField:
::

    generator = MyForceGenerator()
    forcefield._forces.append(generator)

It then should parse the information contained in the XML tag and configure the
generator based on it.

Second, it must define a method with the following signature:
::

    def createForce(self, system, data, nonbondedMethod, nonbondedCutoff, args):

When :code:`createSystem()` is called on the ForceField, it first creates
the System object, then loops over each of its generators and calls
:code:`createForce()` on each one.  This method should create the Force object
and add it to the System.  :code:`data` is a ForceField._SystemData object
containing information about the System being created (atom types, bonds,
angles, etc.), :code:`system` is the System object, and the remaining
arguments are values that were passed to :code:`createSystem()`\ .  To get a
better idea of how this works, look at the existing generator classes in
forcefield.py.

The generator class may optionally also define a method with the following
signature:
::

    def postprocessSystem(self, system, data, args):

If this method exists, it will be called after all Forces have been created.
This gives generators a chance to make additional changes to the System.

Finally, you need to register your class by adding it to ForceField’s map:
::

    forcefield.parsers['MyForce'] = MyForceGenerator.parseElement

The key is the XML tag name, and the value is the static method to use for
parsing it.

Now you can simply create a ForceField object as usual.  If an XML file contains
a :code:`<MyForce>` tag, it will be recognized and processed correctly.

Adding residue template generators
==================================

.. CAUTION::
   This feature is experimental, and its API is subject to change.

Typically, when :class:`ForceField` encounters a residue it does not have a template for,
it simply raises an :code:`Exception`, since it does not know how to assign atom types for
the unknown residue.

However, :class:`ForceField` has an API for registering *residue template generators* that are
called when a residue without an existing template is encountered.  These generators
may create new residue templates that match existing atom types and parameters, or can
even create new atom types and new parameters that are added to :class:`ForceField`. This
functionality can be useful for adding residue template generators that are able to
parameterize small molecules that are not represented in a protein or nucleic acid
forcefield, for example, or for creating new residue templates for post-translationally
modified residues, covalently-bound ligands, or unnatural amino acids or bases.

To register a new residue template generator named :code:`generator`, simply call the
:meth:`registerTemplateGenerator` method on an existing :class:`ForceField` object:
::

    forcefield.registerTemplateGenerator(generator)

This :code:`generator` function must conform to the following API:
::

    def generator(forcefield, residue):
        """
        Parameters
        ----------
        forcefield : openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not
            modify `forcefield`.

        The generator should either register a residue template directly with
        `forcefield.registerResidueTemplate(template)` or it should call `forcefield.loadFile(file)`
        to load residue definitions from an ffxml file.

        It can also use the `ForceField` programmatic API to add additional atom types (via
        `forcefield.registerAtomType(parameters)`) or additional parameters.

        """

The :class:`ForceField` object will be modified by the residue template generator as residues without previously
defined templates are encountered.  Because these templates are added to the :class:`ForceField` as new residue
types are encountered, subsequent residues will be parameterized using the same residue templates without
calling the :code:`generator` again.
.. default-domain:: py

.. _running-simulations:

Running Simulations
###################

.. _a-first-example:

A First Example
***************

Let’s begin with our first example of an OpenMM script. It loads a PDB file
called :file:`input.pdb` that defines a biomolecular system, parameterizes it using the Amber14 force field and TIP3P-FB water
model, energy minimizes it, simulates it for 10,000 steps with a Langevin
integrator, and saves a snapshot frame to a PDB file called :file:`output.pdb` every 1000 time
steps.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *
        from sys import stdout

        pdb = PDBFile('input.pdb')
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,PDB example`

You can find this script in the :file:`examples` folder of your OpenMM installation.
It is called :file:`simulatePdb.py`.  To execute it from a command line, go to your
terminal/console/command prompt window (see Section :numref:`installing-openmm`
on setting up the window to use OpenMM).  Navigate to the :file:`examples` folder by typing
::

    cd <examples_directory>

where the typical directory is :file:`/usr/local/openmm/examples` on Linux
and Mac machines and  :file:`C:\\Program Files\\OpenMM\\examples` on Windows
machines.

Then type
::

    python simulatePdb.py

You can name your own scripts whatever you want.  Let’s go through the script line
by line and see how it works.
::

    from openmm.app import *
    from openmm import *
    from openmm.unit import *
    from sys import stdout

These lines are just telling the Python interpreter about some libraries we will
be using.  Don’t worry about exactly what they mean.  Just include them at the
start of your scripts.
::

    pdb = PDBFile('input.pdb')

This line loads the PDB file from disk.  (The :file:`input.pdb` file in the :file:`examples`
directory contains the villin headpiece in explicit solvent.)  More precisely,
it creates a :class:`PDBFile` object, passes the file name :file:`input.pdb` to it as an
argument, and assigns the object to a variable called :code:`pdb`\ .  The
:class:`PDBFile` object contains the information that was read from the file: the
molecular topology and atom positions.  Your file need not be called
:file:`input.pdb`.  Feel free to change this line to specify any file you want,
though it must contain all of the atoms needed by the force field.
(More information on how to add missing atoms and residues using OpenMM tools can be found in Chapter :numref:`model-building-and-editing`.)
Make sure you include the single quotes around the file name.  OpenMM also can load
files in the newer PDBx/mmCIF format: just change :class:`PDBFile` to :class:`PDBxFile`.
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

This line specifies the force field to use for the simulation.  Force fields are
defined by XML files.  OpenMM includes XML files defining lots of standard force fields (see Section :numref:`force-fields`).
If you find you need to extend the repertoire of force fields available,
you can find more information on how to create these XML files in Chapter :numref:`creating-force-fields`.
In this case we load two of those files: :file:`amber14-all.xml`, which contains the
Amber14 force field, and :file:`amber14/tip3pfb.xml`, which contains the TIP3P-FB water model.  The
:class:`ForceField` object is assigned to a variable called :code:`forcefield`\ .
::

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)

This line combines the force field with the molecular topology loaded from the
PDB file to create a complete mathematical description of the system we want to
simulate.  (More precisely, we invoke the :class:`ForceField` object’s :meth:`.createSystem`
function.  It creates a :class:`System` object, which we assign to the variable
:code:`system`\ .)  It specifies some additional options about how to do that:
use particle mesh Ewald for the long range electrostatic interactions
(:code:`nonbondedMethod=PME`\ ), use a 1 nm cutoff for the direct space
interactions (\ :code:`nonbondedCutoff=1*nanometer`\ ), and constrain the length
of all bonds that involve a hydrogen atom (\ :code:`constraints=HBonds`\ ).
Note the way we specified the cutoff distance 1 nm using :code:`1*nanometer`:
This is an example of the powerful units tracking and automatic conversion facility
built into the OpenMM Python API that makes specifying unit-bearing quantities
convenient and less error-prone.  We could have equivalently specified
:code:`10*angstrom` instead of :code:`1*nanometer` and achieved the same result.
The units system will be described in more detail later, in Section :numref:`units-and-dimensional-analysis`.
::

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

This line creates the integrator to use for advancing the equations of motion.
It specifies a :class:`LangevinMiddleIntegrator`, which performs Langevin dynamics,
and assigns it to a variable called :code:`integrator`\ .  It also specifies
the values of three parameters that are specific to Langevin dynamics: the
simulation temperature (300 K), the friction coefficient (1 ps\ :sup:`-1`\ ), and
the step size (0.004 ps).  Lots of other integration methods are also available.
For example, if you wanted to simulate the system at constant energy rather than
constant temperature you would use a :code:`VerletIntegrator`\ .  The available
integration methods are listed in Section :numref:`integrators`.
::

    simulation = Simulation(pdb.topology, system, integrator)

This line combines the molecular topology, system, and integrator to begin a new
simulation.  It creates a :class:`Simulation` object and assigns it to a variable called
\ :code:`simulation`\ .  A :class:`Simulation` object manages all the processes
involved in running a simulation, such as advancing time and writing output.
::

    simulation.context.setPositions(pdb.positions)

This line specifies the initial atom positions for the simulation: in this case,
the positions that were loaded from the PDB file.
::

    simulation.minimizeEnergy()

This line tells OpenMM to perform a local energy minimization.  It is usually a
good idea to do this at the start of a simulation, since the coordinates in the
PDB file might produce very large forces.
::

    simulation.reporters.append(PDBReporter('output.pdb', 1000))

This line creates a “reporter” to generate output during the simulation, and
adds it to the :class:`Simulation` object’s list of reporters.  A :class:`PDBReporter` writes
structures to a PDB file.  We specify that the output file should be called
:file:`output.pdb`, and that a structure should be written every 1000 time steps.
::

    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True))

It can be useful to get regular status reports as a simulation runs so you can
monitor its progress.  This line adds another reporter to print out some basic
information every 1000 time steps: the current step index, the potential energy
of the system, and the temperature.  We specify :code:`stdout` (not in
quotes) as the output file, which means to write the results to the console.  We
also could have given a file name (in quotes), just as we did for the
:class:`PDBReporter`, to write the information to a file.
::

    simulation.step(10000)

Finally, we run the simulation, integrating the equations of motion for 10,000
time steps.  Once it is finished, you can load the PDB file into any program you
want for analysis and visualization (VMD_, PyMol_, AmberTools_, etc.).

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMol: http://www.pymol.org
.. _AmberTools: http://ambermd.org

.. _using_amber_files:

Using AMBER Files
*****************

OpenMM can build a system in several different ways.  One option, as shown
above, is to start with a PDB file and then select a force field with which to
model it.  Alternatively, you can use AmberTools_ to model your system.  In that
case, you provide a :class:`prmtop` file and an :class:`inpcrd` file.  OpenMM loads the files and
creates a :class:`System` from them.  This is illustrated in the following script.  It can be
found in OpenMM’s :file:`examples` folder with the name :file:`simulateAmber.py`.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *
        from sys import stdout

        prmtop = AmberPrmtopFile('input.prmtop')
        inpcrd = AmberInpcrdFile('input.inpcrd')
        system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,AMBER example`

This script is very similar to the previous one.  There are just a few
significant differences:
::

    prmtop = AmberPrmtopFile('input.prmtop')
    inpcrd = AmberInpcrdFile('input.inpcrd')

In these lines, we load the prmtop file and inpcrd file.  More precisely, we
create :class:`AmberPrmtopFile` and :class:`AmberInpcrdFile` objects and assign them to the
variables :code:`prmtop` and :code:`inpcrd`\ , respectively.  As before,
you can change these lines to specify any files you want.  Be sure to include
the single quotes around the file names.

.. note::

    The :class:`AmberPrmtopFile` reader provided by OpenMM only supports "new-style"
    :file:`prmtop` files introduced in AMBER 7. The AMBER distribution still contains a number of
    example files that are in the "old-style" :file:`prmtop` format. These "old-style" files will
    not run in OpenMM.

Next, the :class:`System` object is created in a different way:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)

In the previous section, we loaded the topology
from a PDB file and then had the force field create a system based on it.  In
this case, we don’t need a force field; the :file:`prmtop` file already contains the
force field parameters, so it can create the system
directly.
::

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

Notice that we now get the topology from the :file:`prmtop` file and the atom positions
from the :file:`inpcrd` file.  In the previous section, both of these came from a PDB
file, but AMBER puts the topology and positions in separate files.  We also add the
following lines:
::

    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

For periodic systems, the :file:`prmtop` file specifies the periodic box vectors, just
as a PDB file does.  When we call :meth:`createSystem`, it sets those as the default
periodic box vectors, to be used automatically for all simulations.  However, the
:file:`inpcrd` may *also* specify periodic box vectors,
and if so we want to use those ones instead.  For example, if the system has been
equilibrated with a barostat, the box vectors may have changed during equilibration.
We therefore check to see if the :file:`inpcrd` file contained box vectors.  If so,
we call :meth:`setPeriodicBoxVectors` to tell it to use those ones, overriding the
default ones provided by the :class:`System`.

.. _using_gromacs_files:

Using Gromacs Files
*******************

A third option for creating your system is to use the Gromacs setup tools.  They
produce a :file:`gro` file containing the coordinates and a :file:`top` file containing the
topology.  OpenMM can load these exactly as it did the AMBER files.  This is
shown in the following script.  It can be found in OpenMM’s :file:`examples` folder
with the name :file:`simulateGromacs.py`.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *
        from sys import stdout

        gro = GromacsGroFile('input.gro')
        top = GromacsTopFile('input.top', periodicBoxVectors=gro.getPeriodicBoxVectors(),
                includeDir='/usr/local/gromacs/share/gromacs/top')
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(top.topology, system, integrator)
        simulation.context.setPositions(gro.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,Gromacs example`

This script is nearly identical to the previous one, just replacing
:class:`AmberInpcrdFile` and :class:`AmberPrmtopFile` with :class:`GromacsGroFile` and :class:`GromacsTopFile`.
Note that when we create the :class:`GromacsTopFile`, we specify values for two extra
options.  First, we specify
:code:`periodicBoxVectors=gro.getPeriodicBoxVectors()`\ .  Unlike OpenMM and
AMBER, which can store periodic unit cell information with the topology, Gromacs
only stores it with the coordinates.  To let :class:`GromacsTopFile` create a :class:`Topology`
object, we therefore need to tell it the periodic box vectors that were loaded
from the :file:`gro` file.  You only need to do this if you are simulating a periodic
system.  For implicit solvent simulations, it usually can be omitted.

Second, we specify :code:`includeDir='/usr/local/gromacs/share/gromacs/top'`\ .  Unlike AMBER,
which stores all the force field parameters directly in a :file:`prmtop` file, Gromacs just stores
references to force field definition files that are installed with the Gromacs
application.  OpenMM needs to know where to find these files, so the
:code:`includeDir` parameter specifies the directory containing them.  If you
omit this parameter, OpenMM will assume the default location :file:`/usr/local/gromacs/share/gromacs/top`,
which is often where they are installed on
Unix-like operating systems.  So in :autonumref:`Example,Gromacs example` we actually could have omitted
this parameter, but if the Gromacs files were installed in any other location,
we would need to include it.

.. _using-charmm-files:

Using CHARMM Files
******************

Yet another option is to load files created by the CHARMM setup tools, or other compatible
tools such as VMD.  Those include a :file:`psf` file containing topology information, and an
ordinary PDB file for the atomic coordinates.  (Coordinates can also be loaded from CHARMM
coordinate or restart files using the :class:`CharmmCrdFile` and :class:`CharmmRstFile` classes).  In addition,
you must provide a set of files containing the force
field definition to use.  This can involve several different files with varying formats and
filename extensions such as :file:`par`, :file:`prm`, :file:`top`, :file:`rtf`, :file:`inp`,
and :file:`str`.  To do this, load all the definition files into a :class:`CharmmParameterSet`
object, then include that object as the first parameter when you call :meth:`createSystem`
on the :class:`CharmmPsfFile`.

.. samepage::
    ::

        from openmm.app import *
        from openmm import *
        from openmm.unit import *
        from sys import stdout, exit, stderr

        psf = CharmmPsfFile('input.psf')
        pdb = PDBFile('input.pdb')
        params = CharmmParameterSet('charmm22.rtf', 'charmm22.prm')
        system = psf.createSystem(params, nonbondedMethod=NoCutoff,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(psf.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        simulation.step(10000)

    .. caption::

        :autonumber:`Example,CHARMM example`

Note that both the CHARMM and XPLOR versions of the :file:`psf` file format are supported.

.. _the-script-builder-application:

The OpenMM-Setup Application
****************************

One way to create your own scripts is to start with one of the examples given
above and customize it to suit your needs, but there's an even easier option.
OpenMM-Setup is a graphical application that walks you through the whole process
of loading your input files and setting options.  It then generates a complete
script, and can even run it for you.

.. figure:: ../../images/OpenMMSetup.png
   :align: center
   :width: 100%

   :autonumber:`Figure,openmm setup`:  The OpenMM-Setup application

To install OpenMM-Setup, open a command line terminal and type the following command
::

    conda install -c conda-forge openmm-setup

You can then launch it by typing the command
::

    openmm-setup

It will automatically open a window in your web browser displaying the user interface.

OpenMM-Setup is far more than just a script generator.  It can fix problems in
your input files, add missing atoms, build membranes and water boxes, and much
more.  It is a very easy way to quickly do all necessary preparation and setup.
We highly recommend it to all users of OpenMM, from novices to experts.

.. _simulation-parameters:

Simulation Parameters
*********************

Now let’s consider lots of ways you might want to customize your script.

Platforms
=========

When creating a :class:`Simulation`, you can optionally tell it what :class:`Platform` to use.
OpenMM includes four platforms: :class:`Reference`, :class:`CPU`, :class:`CUDA`, and :class:`OpenCL`.  For a
description of the differences between them, see Section :numref:`platforms`.  There are three ways in which
the :class:`Platform` can be chosen:

1. By default, OpenMM will try to select the fastest available :class:`Platform`.  Usually its choice will
be reasonable, but sometimes you may want to change it.

2. Alternatively, you can set the :envvar:`OPENMM_DEFAULT_PLATFORM` environment variable to the name
of the :class:`Platform` to use.  This overrides the default logic.

3. Finally, you can explicitly specify a :class:`Platform` object in your script when you create the
:class:`Simulation`.  The following lines specify to use the :class:`CUDA` platform:
::

    platform = Platform.getPlatformByName('CUDA')
    simulation = Simulation(prmtop.topology, system, integrator, platform)

The platform name should be one of :code:`OpenCL`, :code:`CUDA`, :code:`CPU`, or
:code:`Reference`.

You also can specify platform-specific properties that customize how
calculations should be done.  See Chapter :numref:`platform-specific-properties` for details of the
properties that each Platform supports.  For example, the following lines specify to parallelize
work across two different GPUs (CUDA devices 0 and 1), doing all computations in
double precision:
::

    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': '0,1', 'Precision': 'double'}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

.. _force-fields:

Force Fields
============

When you create a force field, you specify one or more XML files from which to
load the force field definition.  Most often, there will be one file to define
the main force field, and possibly a second file to define the water model
(either implicit or explicit).  For example:
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

In some cases, one XML file may load several others.  For example, :file:`amber14-all.xml`
is really just a shortcut for loading several different files that together make up
the AMBER14 force field.  If you need finer grained control over which parameters
are loaded, you can instead specify the component files individually.

Be aware that some force fields and water models include "extra particles", such
as lone pairs or Drude particles.  Examples include the CHARMM polarizable force
field and all of the 4 and 5 site water models.  To use these force fields, you
must first add the extra particles to the :class:`Topology`.  See section
:numref:`adding-or-removing-extra-particles` for details.

The force fields described below are the ones that are bundled with OpenMM.
Additional force fields are available online at https://github.com/choderalab/openmm-forcefields.

Amber14
-------

The Amber14\ :cite:`Maier2015` force field is made up of various files that define
parameters for proteins, DNA, RNA, lipids, water, and ions.

.. tabularcolumns:: |l|L|

===================================  ============================================
File                                 Parameters
===================================  ============================================
:file:`amber14/protein.ff14SB.xml`   Protein (recommended)
:file:`amber14/protein.ff15ipq.xml`  Protein (alternative)
:file:`amber14/DNA.OL15.xml`         DNA (recommended)
:file:`amber14/DNA.bsc1.xml`         DNA (alternative)
:file:`amber14/RNA.OL3.xml`          RNA
:file:`amber14/lipid17.xml`          Lipid
:file:`amber14/GLYCAM_06j-1.xml`     Carbohydrates and glycosylated proteins\ :cite:`Kirschner2007`
:file:`amber14/tip3p.xml`            TIP3P water model\ :cite:`Jorgensen1983` and ions
:file:`amber14/tip3pfb.xml`          TIP3P-FB water model\ :cite:`Wang2014` and ions
:file:`amber14/tip4pew.xml`          TIP4P-Ew water model\ :cite:`Horn2004` and ions
:file:`amber14/tip4pfb.xml`          TIP4P-FB water model\ :cite:`Wang2014` and ions
:file:`amber14/spce.xml`             SPC/E water model\ :cite:`Berendsen1987` and ions
===================================  ============================================

As a convenience, the file :file:`amber14-all.xml` can be used as a shortcut to include
:file:`amber14/protein.ff14SB.xml`, :file:`amber14/DNA.OL15.xml`, :file:`amber14/RNA.OL3.xml`,
and :file:`amber14/lipid17.xml`.  In most cases, you can simply include that file,
plus one of the water models, such as :file:`amber14/tip3pfb.xml` for the TIP3P-FB
water model and ions\ :cite:`Wang2014`:
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

GLYCAM is not included by default, since it is quite large.  If your system contains
carbohydrates, include that file as well:
::

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml', 'amber14/GLYCAM_06j-1.xml')

Be aware that GLYCAM works somewhat differently from most force fields.  It uses
its own nonstandard `naming convention <http://legacy.glycam.org/docs/forcefield/glycam-naming-2/index.html>`_
for carbohydrates, and requires your input file to follow it.  If any residues have
names different from what it expects, GLYCAM will be unable to assign parameters
to them.

.. tip:: The solvent model XML files included under the :file:`amber14/` directory
         include both water *and* ions compatible with that water model, so if you
         mistakenly specify :file:`tip3p.xml` instead of :file:`amber14/tip3p.xml`,
         you run the risk of having :class:`ForceField` throw an exception since
         :file:`tip3p.xml` will be missing parameters for ions in your system.

The converted parameter sets come from the `AmberTools 17 release <http://ambermd.org/AmberTools.php>`_
and were converted using the `openmm-forcefields <https://github.com/choderalab/openmm-forcefields>`_ package and `ParmEd <https://github.com/parmed/parmed>`_.

CHARMM36
--------

The CHARMM36\ :cite:`Best2012` force field provides parameters for proteins, DNA,
RNA, lipids, carbohydrates, water, ions, and various small molecules (see `here <http://mackerell.umaryland.edu/charmm_ff.shtml#refs>`_
for full references).

.. tabularcolumns:: |l|L|

=================================  ============================================
File                               Parameters
=================================  ============================================
:file:`charmm36.xml`               Protein, DNA, RNA, lipids, carbohydrates, and small molecules
:file:`charmm36/water.xml`         Default CHARMM water model (a modified version of TIP3P\ :cite:`Jorgensen1983`) and ions
:file:`charmm36/spce.xml`          SPC/E water model\ :cite:`Berendsen1987` and ions
:file:`charmm36/tip3p-pme-b.xml`   TIP3P-PME-B water model\ :cite:`Price2004` and ions
:file:`charmm36/tip3p-pme-f.xml`   TIP3P-PME-F water model\ :cite:`Price2004` and ions
:file:`charmm36/tip4pew.xml`       TIP4P-Ew water model\ :cite:`Horn2004` and ions
:file:`charmm36/tip4p2005.xml`     TIP4P-2005 water model\ :cite:`Abascal2005` and ions
:file:`charmm36/tip5p.xml`         TIP5P water model\ :cite:`Mahoney2000` and ions
:file:`charmm36/tip5pew.xml`       TIP5P-Ew water model\ :cite:`Rick2004` and ions
=================================  ============================================

The file :file:`charmm36.xml` bundles everything but the water and ions into a single
file.  In most cases, you can simply include that file, plus one of the water models,
such as :file:`charmm36/water.xml`, which specifies the default CHARMM water model
(a modified version of TIP3P\ :cite:`Jorgensen1983`) and ions:
::

    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')

.. warning:: Drude polarizable sites and lone pairs are not yet supported
             by `ParmEd <https://github.com/parmed/parmed>`_ and the CHARMM36 forcefields
             that depend on these features are not included in this port.
             To use the CHARMM 2019 polarizable force field\ :cite:`Lopes2013`,
             include the single file :file:`charmm_polar_2019.xml`.

.. tip:: The solvent model XML files included under the :file:`charmm36/` directory
         include both water *and* ions compatible with that water model, so if you
         mistakenly specify :file:`tip3p.xml` instead of :file:`charmm36/water.xml`,
         you run the risk of having :class:`ForceField` raise an exception due to
         missing parameters for ions in your system.

.. tip:: CHARMM makes extensive use of patches, which are automatically combined with
         residue templates to create an expanded library of patched residue templates
         by :class:`ForceField`. That means that patched residues, such as ``ACE`` and
         ``NME`` patched termini, must occur as a single residue in order for :class:`ForceField`
         to correctly match the residue template and apply parameters. Since these
         patched residues are not standard PDB residues, :class:`Modeller` does not know
         how to add hydrogens to these nonstandard residues, and your input topologies
         must already contain appropriate hydrogens. This can often cause problems when
         trying to read in PDB files from sources such as `CHARMM-GUI <http://charmm-gui.org/>`_
         that do not generate PDB files that comply with the `PDB standard <http://www.wwpdb.org/documentation/file-format>`_.
         If you're using files from `CHARMM-GUI <http://charmm-gui.org/>`_, it's easiest to load
         the PSF file directly, as discussed in Section :numref:`using-charmm-files`.

.. tip:: Trying to read in PDB files from sources such as `CHARMM-GUI <http://charmm-gui.org/>`_
         that do not generate PDB files that comply with the `PDB standard <http://www.wwpdb.org/documentation/file-format>`_
         and omit ``CONECT`` records specifying bonds between residues (such as cysteines)
         or include ``CONECT`` records specifying non-chemical ``H-H`` bonds in waters
         can cause issues with the detection and parameter assignment for disulfide bonds.
         Make sure the files you read in comply with the appropriate standards regarding
         additional bonds and nonstandard residue definitions. If you're using files from
         `CHARMM-GUI <http://charmm-gui.org/>`_, it's easiest to load
         the PSF file directly, as discussed in Section :numref:`using-charmm-files`.

The converted parameter sets come from the `CHARMM36 July 2017 update <http://mackerell.umaryland.edu/charmm_ff.shtml>`_
and were converted using the `openmm-forcefields <https://github.com/choderalab/openmm-forcefields>`_ package and `parmed <https://github.com/parmed/parmed>`_.

Implicit Solvent
----------------

The Amber and CHARMM force fields described above can be used with any of the Generalized
Born implicit solvent models from AMBER.  To use them, include an extra file when
creating the ForceField.  For example,
::

    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

.. tabularcolumns:: |l|L|

==========================  ==================================================================================================================================
File                        Implicit Solvent Model
==========================  ==================================================================================================================================
:file:`implicit/hct.xml`    Hawkins-Cramer-Truhlar GBSA model\ :cite:`Hawkins1995` (corresponds to igb=1 in AMBER)
:file:`implicit/obc1.xml`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ I parameters (corresponds to igb=2 in AMBER).
:file:`implicit/obc2.xml`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ II parameters (corresponds to igb=5 in AMBER).
:file:`implicit/gbn.xml`    GBn solvation model\ :cite:`Mongan2007` (corresponds to igb=7 in AMBER).
:file:`implicit/gbn2.xml`   GBn2 solvation model\ :cite:`Nguyen2013` (corresponds to igb=8 in AMBER).
==========================  ==================================================================================================================================

The only nonbonded methods that are supported with implicit solvent are :code:`NoCutoff` (the default),
:code:`CutoffNonPeriodic`, and :code:`CutoffPeriodic.`  If you choose to use a nonbonded cutoff with
implicit solvent, it is usually best to set the cutoff distance larger than is typical with explicit solvent.
A cutoff of 2 nm gives good results in most cases.  Periodic boundary conditions are not usually used
with implicit solvent.  In fact, the lack of need for periodicity and the artifacts it creates is one
of the advantages of implicit solvent.  The option is still offered, since it could be useful in some
unusual situations.

You can further control the solvation model in a few ways.  First, you can
specify the dielectric constants to use for the solute and solvent:
::

    system = forcefield.createSystem(topology, soluteDielectric=1.0, solventDielectric=80.0)

If they are not specified, the solute and solvent dielectric constants default to 1.0 and
78.5, respectively.

You also can model the effect of a non-zero salt concentration by specifying the
Debye-Huckel screening parameter\ :cite:`Srinivasan1999`:
::

    system = forcefield.createSystem(topology, implicitSolventKappa=1.0/nanometer)

The screening parameter can be calculated as

.. math::
  \kappa = 367.434915 \sqrt{\frac{I}{\epsilon T}}

where :math:`I` is the ionic strength in moles/liter, :math:`\epsilon` is the solvent
dielectric constant, and :math:`T` is the temperature in Kelvin.

AMOEBA
------

The AMOEBA polarizable force field provides parameters for proteins, nucleic acids, water, and ions.

.. tabularcolumns:: |l|L|

=============================  ================================================================================
File                           Parameters
=============================  ================================================================================
:file:`amoeba2018.xml`         AMOEBA 2018\ :cite:`Shi2013`
:file:`amoeba2018_gk.xml`      Generalized Kirkwood solvation model\ :cite:`Schnieders2007` for use with AMOEBA 2018 force field
:file:`amoeba2013.xml`         AMOEBA 2013.  This force field is deprecated.  It is
                               recommended to use AMOEBA 2018 instead.
:file:`amoeba2013_gk.xml`      Generalized Kirkwood solvation model for use with AMOEBA 2013 force field
:file:`amoeba2009.xml`         AMOEBA 2009\ :cite:`Ren2002`.  This force field is deprecated.  It is
                               recommended to use AMOEBA 2018 instead.
:file:`amoeba2009_gk.xml`      Generalized Kirkwood solvation model for use with AMOEBA 2009 force field
=============================  ================================================================================

For explicit solvent simulations, just include the single file :file:`amoeba2018.xml`.
AMOEBA also supports implicit solvent using a Generalized Kirkwood model.  To enable
it, also include :file:`amoeba2018_gk.xml`.

The older AMOEBA 2009 and 2013 force fields are provided only for backward compatibility, and are not
recommended for most simulations.

CHARMM Polarizable Force Field
------------------------------

To use the CHARMM 2019 polarizable force field\ :cite:`Lopes2013`, include the
single file :file:`charmm_polar_2019.xml`.  It includes parameters for proteins, lipids,
water, and ions.  When using this force field, remember to add extra particles to
the :class:`Topology` as described in section :numref:`adding-or-removing-extra-particles`.
This force field also requires that you use one of the special integrators that
supports Drude particles.  The options are DrudeLangevinIntegrator, DrudeNoseHooverIntegrator,
and DrudeSCFIntegrator.

Older Force Fields
------------------

OpenMM includes several older force fields as well.  For most simulations, the
newer force fields described above are preferred over any of these, but they are
still useful for reproducing older results.

.. tabularcolumns:: |l|L|

=============================  ================================================================================
File                           Force Field
=============================  ================================================================================
:code:`amber96.xml`            Amber96\ :cite:`Kollman1997`
:code:`amber99sb.xml`          Amber99\ :cite:`Wang2000` with modified backbone torsions\ :cite:`Hornak2006`
:code:`amber99sbildn.xml`      Amber99SB plus improved side chain torsions\ :cite:`Lindorff-Larsen2010`
:code:`amber99sbnmr.xml`       Amber99SB with modifications to fit NMR data\ :cite:`Li2010`
:code:`amber03.xml`            Amber03\ :cite:`Duan2003`
:code:`amber10.xml`            Amber10 (documented in the AmberTools_ manual as `ff10`)
:code:`charmm_polar_2013.xml`  2013 version of the CHARMM polarizable force field\ :cite:`Lopes2013`
=============================  ================================================================================

Several of these force fields support implicit solvent.  To enable it, also
include the corresponding OBC file.

.. tabularcolumns:: |l|L|

=========================  =================================================================================================
File                       Implicit Solvation Model
=========================  =================================================================================================
:code:`amber96_obc.xml`    GBSA-OBC solvation model\ :cite:`Onufriev2004` for use with Amber96 force field
:code:`amber99_obc.xml`    GBSA-OBC solvation model for use with Amber99 force field and its variants
:code:`amber03_obc.xml`    GBSA-OBC solvation model for use with Amber03 force field
:code:`amber10_obc.xml`    GBSA-OBC solvation model for use with Amber10 force field
=========================  =================================================================================================

Note that the GBSA-OBC parameters in these files are those used in TINKER.\ :cite:`Tinker`
They are designed for use with Amber force fields, but they are different from
the parameters found in the AMBER application.

Water Models
------------

The following files define popular water models.  They can be used with force fields
that do not provide their own water models.  When using Amber14 or CHARMM36, use
the water files included with those force fields instead, since they also include
ion parameters.

.. tabularcolumns:: |l|L|

===================  ============================================
File                 Water Model
===================  ============================================
:code:`tip3p.xml`    TIP3P water model\ :cite:`Jorgensen1983`
:code:`tip3pfb.xml`  TIP3P-FB water model\ :cite:`Wang2014`
:code:`tip4pew.xml`  TIP4P-Ew water model\ :cite:`Horn2004`
:code:`tip4pfb.xml`  TIP4P-FB water model\ :cite:`Wang2014`
:code:`tip5p.xml`    TIP5P water model\ :cite:`Mahoney2000`
:code:`spce.xml`     SPC/E water model\ :cite:`Berendsen1987`
:code:`swm4ndp.xml`  SWM4-NDP water model\ :cite:`Lamoureux2006`
===================  ============================================

Small molecule parameters
=========================

The OpenMM force fields above include pregenerated templates for biopolymers
and solvents. If your system instead contain small molecules, it is often
necessary to generate these parameters on the fly.


There are two options for doing this within the OpenMM ``app`` ecosystem:

Small molecule residue template generators
------------------------------------------

One approach is to use residue template generators for small molecules from the
openmmforcefields_  conda package.
You can install this via conda with:

.. code-block:: bash

    $ conda install -c conda-forge openmmforcefields

You can then add a small molecule residue template generator using the Open Force
Field Initiative small molecule force fields using the following example:

::

    # Create an OpenFF Molecule object for benzene from SMILES
    from openff.toolkit.topology import Molecule
    molecule = Molecule.from_smiles('c1ccccc1')
    # Create the SMIRNOFF template generator with the default installed force field (openff-1.0.0)
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
    # Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
    from openmm.app import ForceField
    forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
    # Register the SMIRNOFF template generator
    forcefield.registerTemplateGenerator(smirnoff.generator)

Alternatively, you can use the older `AMBER GAFF small molecule force field <http://ambermd.org/antechamber/gaff.html>`_:

::

    # Create an OpenFF Molecule object for benzene from SMILES
    from openff.toolkit.topology import Molecule
    molecule = Molecule.from_smiles('c1ccccc1')
    # Create the GAFF template generator
    from openmmforcefields.generators import GAFFTemplateGenerator
    gaff = GAFFTemplateGenerator(molecules=molecule)
    # Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
    from openmm.app import ForceField
    forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
    # Register the GAFF template generator
    forcefield.registerTemplateGenerator(gaff.generator)
    # You can now parameterize an OpenMM Topology object that contains the specified molecule.
    # forcefield will load the appropriate GAFF parameters when needed, and antechamber
    # will be used to generate small molecule parameters on the fly.
    from openmm.app import PDBFile
    pdbfile = PDBFile('t4-lysozyme-L99A-with-benzene.pdb')
    system = forcefield.createSystem(pdbfile.topology)

More documentation can be found on the openmmforcefields_ page.

Managing force fields with ``SystemGenerator``
----------------------------------------------

As an alternative to explicitly registering template generators, the openmmforcefields_
package provides a ``SystemGenerator`` facility to simplify biopolymer and
small molecule force field management. To use this, you can simply specify the
small molecule force field you want to use:

::

    # Define the keyword arguments to feed to ForceField
    from openmm import unit, app
    forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
    # Initialize a SystemGenerator using GAFF
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'], small_molecule_forcefield='gaff-2.11', forcefield_kwargs=forcefield_kwargs, cache='db.json')
    # Create an OpenMM System from an OpenMM Topology object
    system = system_generator.create_system(openmm_topology)
    # Alternatively, create an OpenMM System from an OpenMM Topology object and a list of OpenFF Molecule objects
    molecules = Molecule.from_file('molecules.sdf', file_format='sdf')
    system = system_generator.create_system(openmm_topology, molecules=molecules)

The ``SystemGenerator`` will match any instances of the molecules found in ``molecules.sdf`` to those that appear in ``topology``.
Note that the protonation and tautomeric states must match exactly between the ``molecules`` read and those appearing in the Topology.
See the openmmforcefields_ documentation for more details.

.. _openmmforcefields: http://github.com/openmm/openmmforcefields

AMBER Implicit Solvent
======================


When creating a system from a prmtop file you do not specify force field files,
so you need a different way to tell it to use implicit solvent.  This is done
with the :code:`implicitSolvent` parameter:
::

    system = prmtop.createSystem(implicitSolvent=OBC2)

OpenMM supports all of the Generalized Born models used by AMBER.  Here are the
allowed values for :code:`implicitSolvent`\ :

.. tabularcolumns:: |l|L|

=============  ==================================================================================================================================
Value          Meaning
=============  ==================================================================================================================================
:code:`None`   No implicit solvent is used.
:code:`HCT`    Hawkins-Cramer-Truhlar GBSA model\ :cite:`Hawkins1995` (corresponds to igb=1 in AMBER)
:code:`OBC1`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ I parameters (corresponds to igb=2 in AMBER).
:code:`OBC2`   Onufriev-Bashford-Case GBSA model\ :cite:`Onufriev2004` using the GB\ :sup:`OBC`\ II parameters (corresponds to igb=5 in AMBER).
               This is the same model used by the GBSA-OBC files described in Section :numref:`force-fields`.
:code:`GBn`    GBn solvation model\ :cite:`Mongan2007` (corresponds to igb=7 in AMBER).
:code:`GBn2`   GBn2 solvation model\ :cite:`Nguyen2013` (corresponds to igb=8 in AMBER).
=============  ==================================================================================================================================


You can further control the solvation model in a few ways.  First, you can
specify the dielectric constants to use for the solute and solvent:
::

    system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0,
            solventDielectric=80.0)

If they are not specified, the solute and solvent dielectrics default to 1.0 and
78.5, respectively.  These values were chosen for consistency with AMBER, and
are slightly different from those used elsewhere in OpenMM: when building a
system from a force field, the solvent dielectric defaults to 78.3.

You also can model the effect of a non-zero salt concentration by specifying the
Debye-Huckel screening parameter\ :cite:`Srinivasan1999`:
::

    system = prmtop.createSystem(implicitSolvent=OBC2, implicitSolventKappa=1.0/nanometer)


Nonbonded Interactions
======================


When creating the system (either from a force field or a prmtop file), you can
specify options about how nonbonded interactions should be treated:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer)

The :code:`nonbondedMethod` parameter can have any of the following values:

.. tabularcolumns:: |l|L|

=========================  ===========================================================================================================================================================================================================================================
Value                      Meaning
=========================  ===========================================================================================================================================================================================================================================
:code:`NoCutoff`           No cutoff is applied.
:code:`CutoffNonPeriodic`  The reaction field method is used to eliminate all interactions beyond a cutoff distance.  Not valid for AMOEBA.
:code:`CutoffPeriodic`     The reaction field method is used to eliminate all interactions beyond a cutoff distance.  Periodic boundary conditions are applied, so each atom interacts only with the nearest periodic copy of every other atom.  Not valid for AMOEBA.
:code:`Ewald`              Periodic boundary conditions are applied.  Ewald summation is used to compute long range Coulomb interactions.  (This option is rarely used, since PME is much faster for all but the smallest systems.)  Not valid for AMOEBA.
:code:`PME`                Periodic boundary conditions are applied.  The Particle Mesh Ewald method is used to compute long range Coulomb interactions.
:code:`LJPME`              Periodic boundary conditions are applied.  The Particle Mesh Ewald method is used to compute long range interactions for both Coulomb and Lennard-Jones.
=========================  ===========================================================================================================================================================================================================================================


When using any method other than :code:`NoCutoff`\ , you should also specify a
cutoff distance.  Be sure to specify units, as shown in the examples above. For
example, :code:`nonbondedCutoff=1.5*nanometers` or
:code:`nonbondedCutoff=12*angstroms` are legal values.

When using :code:`Ewald`, :code:`PME`, or :code:`LJPME`\ , you can optionally specify an
error tolerance for the force computation.  For example:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            ewaldErrorTolerance=0.00001)

The error tolerance is roughly equal to the fractional error in the forces due
to truncating the Ewald summation.  If you do not specify it, a default value of
0.0005 is used.

Another optional parameter when using a cutoff is :code:`switchDistance`.  This
causes Lennard-Jones interactions to smoothly go to zero over some finite range,
rather than being sharply truncated at the cutoff distance.  This can improve
energy conservation.  To use it, specify a distance at which the interactions
should start being reduced.  For example:
::

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            switchDistance=0.9*nanometer)


Nonbonded Forces for AMOEBA
---------------------------

For the AMOEBA force field, the valid values for the :code:`nonbondedMethod`
are :code:`NoCutoff` and :code:`PME`\ .  The other nonbonded methods,
:code:`CutoffNonPeriodic`\ , :code:`CutoffPeriodic`\ , and :code:`Ewald`
are unavailable for this force field.

For implicit solvent runs using AMOEBA, only the :code:`nonbondedMethod`
option :code:`NoCutoff` is available.

Lennard-Jones Interaction Cutoff Value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition, for the AMOEBA force field a cutoff for the Lennard-Jones
interaction independent of the value used for the electrostatic interactions may
be specified using the keyword :code:`vdwCutoff`\ .
::

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            ewaldErrorTolerance=0.00001, vdwCutoff=1.2*nanometer)

If :code:`vdwCutoff` is not specified, then the value of
:code:`nonbondedCutoff` is used for the Lennard-Jones interactions.

Specifying the Polarization Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the AMOEBA force field, OpenMM allows the induced dipoles to be
calculated in any of three different ways.  The slowest but potentially most
accurate method is to iterate the calculation until the dipoles converge to a
specified tolerance.  To select this, specify :code:`polarization='mutual'`.
Use the :code:`mutualInducedTargetEpsilon` option to select the tolerance; for
most situations, a value of 0.00001 works well.  Alternatively you can specify
:code:`polarization='extrapolated'`.  This uses an analytic approximation
:cite:`Simmonett2015` to estimate what the fully converged dipoles will be without
actually continuing the calculation to convergence.  In many cases this can be
significantly faster with only a small loss in accuracy.  Finally, you can
specify :code:`polarization='direct'` to use the direct polarization
approximation, in which induced dipoles depend only on the fixed multipoles, not
on other induced dipoles.  This is even faster, but it produces very different
forces from mutual polarization, so it should only be used with force fields
that have been specifically parameterized for use with this approximation.

Here are examples of using each method:
::

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='mutual', mutualInducedTargetEpsilon=0.00001)

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='extrapolated')

    system = forcefield.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        vdwCutoff=1.2*nanometer, polarization='direct')


Implicit Solvent and Solute Dielectrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For implicit solvent simulations using the AMOEBA force field, the
:file:`amoeba2013_gk.xml` file should be included in the initialization of the force
field:
::

    forcefield = ForceField('amoeba2009.xml', 'amoeba2009_gk.xml')

Only the :code:`nonbondedMethod` option :code:`NoCutoff` is available
for implicit solvent runs using AMOEBA.  In addition, the solvent and solute
dielectric values can be specified for implicit solvent simulations:
::

    system=forcefield.createSystem(nonbondedMethod=NoCutoff, soluteDielectric=2.0,
            solventDielectric=80.0)

The default values are 1.0 for the solute dielectric and 78.3 for the solvent
dielectric.

Constraints
===========


When creating the system (either from a force field or an AMBER :file:`prmtop` file), you can
optionally tell OpenMM to constrain certain bond lengths and angles.  For
example,
::

    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds)

The :code:`constraints` parameter can have any of the following values:

.. tabularcolumns:: |l|L|

================  =============================================================================================================================================
Value             Meaning
================  =============================================================================================================================================
:code:`None`      No constraints are applied.  This is the default value.
:code:`HBonds`    The lengths of all bonds that involve a hydrogen atom are constrained.
:code:`AllBonds`  The lengths of all bonds are constrained.
:code:`HAngles`   The lengths of all bonds are constrained.  In addition, all angles of the form H-X-H or H-O-X (where X is an arbitrary atom) are constrained.
================  =============================================================================================================================================


The main reason to use constraints is that it allows one to use a larger
integration time step.  With no constraints, one is typically limited to a time
step of about 1 fs for typical biomolecular force fields like AMBER or CHARMM.
With :code:`HBonds` constraints, this can be increased to about 2 fs for Verlet
dynamics, or about 4 fs for Langevin dynamics.  With :code:`HAngles`\ , it can
sometimes be increased even further.

Regardless of the value of this parameter, OpenMM makes water molecules
completely rigid, constraining both their bond lengths and angles.  You can
disable this behavior with the :code:`rigidWater` parameter:
::

    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)

Be aware that flexible water may require you to further reduce the integration
step size, typically to about 0.5 fs.

.. note::

   The AMOEBA forcefield is designed to be used without constraints, so by
   default OpenMM makes AMOEBA water flexible.  You can still force it to be
   rigid by specifying :code:`rigidWater=True`.

Heavy Hydrogens
===============


When creating the system (either from a force field or an AMBER :file:`prmtop` file), you can
optionally tell OpenMM to increase the mass of hydrogen atoms.  For example,
::

    system = prmtop.createSystem(hydrogenMass=4*amu)

This applies only to hydrogens that are bonded to heavy atoms, and any mass
added to the hydrogen is subtracted from the heavy atom.  This keeps their total
mass constant while slowing down the fast motions of hydrogens.  When combined
with constraints (typically :code:`constraints=AllBonds`\ ), this often allows a
further increase in integration step size.

.. _integrators:

Integrators
===========


OpenMM offers a choice of several different integration methods.  You select
which one to use by creating an integrator object of the appropriate type.
Detailed descriptions of all these integrators can be found in Chapter
:numref:`integrators-theory`.  In addition to these built in integrators, lots of
others are available as part of the `OpenMMTools <https://openmmtools.readthedocs.io>`_ package.

Langevin Middle Integrator
--------------------------

In the examples of the previous sections, we used Langevin integration:
::

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

The three parameter values in this line are the simulation temperature (300 K),
the friction coefficient (1 ps\ :sup:`-1`\ ), and the step size (0.004 ps).  You
are free to change these to whatever values you want.  Be sure to specify units
on all values.  For example, the step size could be written either as
:code:`0.004*picoseconds` or :code:`4*femtoseconds`\ .  They are exactly
equivalent.  Note that :code:`LangevinMiddleIntegrator` is a leapfrog
integrator, so the velocities are offset by half a time step from the positions.

Langevin Integrator
-------------------

:code:`LangevinIntegrator` is very similar to :code:`LangevinMiddleIntegrator`,
but it uses a different discretization of the Langevin equation.
:code:`LangevinMiddleIntegrator` tends to produce more accurate configurational
sampling, and therefore is preferred for most applications.  Also note that
:code:`LangevinIntegrator`\ , like :code:`LangevinMiddleIntegrator`\ , is a leapfrog
integrator, so the velocities are offset by half a time step from the positions.

Nosé-Hoover Integrator
----------------------

The :code:`NoseHooverIntegrator` uses the same "middle" leapfrog propagation
algorithm as :code:`LangevinMiddleIntegrator`, but replaces the stochastic
temperature control with a velocity scaling algorithm that produces more
accurate transport properties :cite:`Basconi2013`.  This velocity scaling
results from propagating a chain of extra variables, which slightly reduces the
computational efficiency with respect to :code:`LangevinMiddleIntegrator`.  The
thermostated integrator is minimally created with syntax analogous to the
:code:`LangevinMiddleIntegrator` example above::

    NoseHooverIntegrator integrator(300*kelvin, 1/picosecond,
                                    0.004*picoseconds);

The first argument specifies the target temperature.  The second specifies the
frequency of interaction with the heat bath: a lower value interacts minimally,
yielding the microcanonical ensemble in the limit of a zero frequency, while a
larger frequency will perturb the system greater, keeping it closer to the
target temperature.  The third argument is the integration timestep that, like
the other arguments, must be specified with units.  For initial equilibration
to the target temperature, a larger interaction frequency is recommended,
*e.g.* 25 ps\ :sup:`-1`.

This integrator supports lots of other options, including the ability to couple
different parts of the system to thermostats at different temperatures. See the
API documentation for details.

Leapfrog Verlet Integrator
--------------------------

A leapfrog Verlet integrator can be used for running constant energy dynamics.
The command for this is:
::

    integrator = VerletIntegrator(0.002*picoseconds)

The only option is the step size.

Brownian Integrator
-------------------

Brownian (diffusive) dynamics can be used by specifying the following:
::

    integrator = BrownianIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

The parameters are the same as for Langevin dynamics: temperature (300 K),
friction coefficient (1 ps\ :sup:`-1`\ ), and step size (0.002 ps).

Variable Time Step Langevin Integrator
--------------------------------------

A variable time step Langevin integrator continuously adjusts its step size to
keep the integration error below a specified tolerance.  In some cases, this can
allow you to use a larger average step size than would be possible with a fixed
step size integrator.  It also is very useful in cases where you do not know in
advance what step size will be stable, such as when first equilibrating a
system.  You create this integrator with the following command:
::

    integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.001)

In place of a step size, you specify an integration error tolerance (0.001 in
this example).  It is best not to think of this value as having any absolute
meaning.  Just think of it as an adjustable parameter that affects the step size
and integration accuracy.  Smaller values will produce a smaller average step
size.  You should try different values to find the largest one that produces a
trajectory sufficiently accurate for your purposes.

Variable Time Step Leapfrog Verlet Integrator
---------------------------------------------

A variable time step leapfrog Verlet integrator works similarly to the variable
time step Langevin integrator in that it continuously adjusts its step size to
keep the integration error below a specified tolerance.  The command for this
integrator is:
::

    integrator = VariableVerletIntegrator(0.001)

The parameter is the integration error tolerance (0.001), whose meaning is the
same as for the Langevin integrator.

Multiple Time Step Integrator
-----------------------------

The :class:`MTSIntegrator` class implements the rRESPA multiple time step
algorithm\ :cite:`Tuckerman1992`.  This allows some forces in the system to be evaluated more
frequently than others.  For details on how to use it, consult the API
documentation.

Multiple Time Step Langevin Integrator
--------------------------------------

:class:`MTSLangevinIntegrator` is similar to :class:`MTSIntegrator`, but it uses
the Langevin method to perform constant temperature dynamics.  For details on
how to use it, consult the API documentation.

Compound Integrator
-------------------

The :class:`CompoundIntegrator` class is useful for cases where you want to use
multiple integration algorithms within a single simulation.  It allows you to
create multiple integrators, then switch back and forth between them.  For
details on how to use it, consult the API documentation.

Temperature Coupling
====================


If you want to run a simulation at constant temperature, using a Langevin
integrator (as shown in the examples above) is usually the best way to do it.
OpenMM does provide an alternative, however: you can use a Verlet integrator,
then add an Andersen thermostat to your system to provide temperature coupling.

To do this, we can add an :class:`AndersenThermostat` object to the :class:`System` as shown below.
::

    ...
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)
    system.addForce(AndersenThermostat(300*kelvin, 1/picosecond))
    integrator = VerletIntegrator(0.002*picoseconds)
    ...

The two parameters of the Andersen thermostat are the temperature (300 K) and
collision frequency (1 ps\ :sup:`-1`\ ).

Pressure Coupling
=================


All the examples so far have been constant volume simulations.  If you want to
run at constant pressure instead, add a Monte Carlo barostat to your system.
You do this exactly the same way you added the Andersen thermostat in the
previous section:
::

    ...
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    ...

The parameters of the Monte Carlo barostat are the pressure (1 bar) and
temperature (300 K).  The barostat assumes the simulation is being run at
constant temperature, but it does not itself do anything to regulate the
temperature.

.. warning::

    It is therefore critical that you always use it along with a Langevin integrator or
    Andersen thermostat, and that you specify the same temperature for both the barostat
    and the integrator or thermostat.  Otherwise, you will get incorrect results.

There also is an anisotropic barostat that scales each axis of the periodic box
independently, allowing it to change shape.  When using the anisotropic
barostat, you can specify a different pressure for each axis.  The following
line applies a pressure of 1 bar along the X and Y axes, but a pressure of 2 bar
along the Z axis:
::

    system.addForce(MonteCarloAnisotropicBarostat((1, 1, 2)*bar, 300*kelvin))

Another feature of the anisotropic barostat is that it can be applied to only
certain axes of the periodic box, keeping the size of the other axes fixed.
This is done by passing three additional parameters that specify whether the
barostat should be applied to each axis.  The following line specifies that the
X and Z axes of the periodic box should not be scaled, so only the Y axis can
change size.
::

    system.addForce(MonteCarloAnisotropicBarostat((1, 1, 1)*bar, 300*kelvin,
            False, True, False))

There is a third barostat designed specifically for simulations of membranes.
It assumes the membrane lies in the XY plane, and treats the X and Y axes of the
box differently from the Z axis.  It also applies a uniform surface tension in
the plane of the membrane.  The following line adds a membrane barostat that
applies a pressure of 1 bar and a surface tension of 200 bar*nm.  It specifies
that the X and Y axes are treated isotropically while the Z axis is free to
change independently.
::

    system.addForce(MonteCarloMembraneBarostat(1*bar, 200*bar*nanometer,
        MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree, 300*kelvin))

See the API documentation for details about the allowed parameter values and
their meanings.


Energy Minimization
===================


As seen in the examples, performing a local energy minimization takes a single
line in the script:
::

    simulation.minimizeEnergy()

In most cases, that is all you need.  There are two optional parameters you can
specify if you want further control over the minimization.  First, you can
specify a tolerance for when the energy should be considered to have converged:
::

    simulation.minimizeEnergy(tolerance=5*kilojoule/mole)

If you do not specify this parameter, a default tolerance of 10 kJ/mole is used.

Second, you can specify a maximum number of iterations:
::

    simulation.minimizeEnergy(maxIterations=100)

The minimizer will exit once the specified number of iterations is reached, even
if the energy has not yet converged.  If you do not specify this parameter, the
minimizer will continue until convergence is reached, no matter how many
iterations it takes.

These options are independent.  You can specify both if you want:
::

    simulation.minimizeEnergy(tolerance=0.1*kilojoule/mole, maxIterations=500)

Removing Center of Mass Motion
==============================


By default, :class:`System` objects created with the OpenMM application tools add
a :class:`CMMotionRemover` that removes all center of mass motion at every time step so the
system as a whole does not drift with time.  This is almost always what you
want.  In rare situations, you may want to allow the system to drift with time.
You can do this by specifying the :code:`removeCMMotion` parameter when you
create the System:
::

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
            removeCMMotion=False)

Writing Trajectories
====================


OpenMM can save simulation trajectories to disk in three formats: PDB_,
`PDBx/mmCIF`_, and DCD_.  All of these are widely supported formats, so you
should be able to read them into most analysis and visualization programs.

.. _PDB: http://www.wwpdb.org/documentation/format33/v3.3.html
.. _PDBx/mmCIF: http://mmcif.wwpdb.org
.. _DCD: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

To save a trajectory, just add a “reporter” to the simulation, as shown in the
example scripts above:
::

    simulation.reporters.append(PDBReporter('output.pdb', 1000))

The two parameters of the :class:`PDBReporter` are the output filename and how often (in
number of time steps) output structures should be written.  To use PDBx/mmCIF or
DCD format, just replace :class:`PDBReporter` with :class:`PDBxReporter` or
:class:`DCDReporter`.  The parameters represent the same values:
::

    simulation.reporters.append(DCDReporter('output.dcd', 1000))

Recording Other Data
====================


In addition to saving a trajectory, you may want to record other information
over the course of a simulation, such as the potential energy or temperature.
OpenMM provides a reporter for this purpose also.  Create a :class:`StateDataReporter`
and add it to the simulation:
::

    simulation.reporters.append(StateDataReporter('data.csv', 1000, time=True,
            kineticEnergy=True, potentialEnergy=True))

The first two parameters are the output filename and how often (in number of
time steps) values should be written.  The remaining arguments specify what
values should be written at each report.  The available options are
:code:`step` (the index of the current time step), :code:`time`\ ,
:code:`progress` (what percentage of the simulation has completed),
:code:`remainingTime` (an estimate of how long it will take the simulation to
complete), :code:`potentialEnergy`\ , :code:`kineticEnergy`\ ,
:code:`totalEnergy`\ , :code:`temperature`\ , :code:`volume` (the volume
of the periodic box), :code:`density` (the total system mass divided by the
volume of the periodic box), and :code:`speed` (an estimate of how quickly
the simulation is running).  If you include either the :code:`progress` or
:code:`remainingTime` option, you must also include the :code:`totalSteps`
parameter to specify the total number of time steps that will be included in the
simulation.  One line is written to the file for each report containing the
requested values.  By default the values are written in comma-separated-value
(CSV) format.  You can use the :code:`separator` parameter to choose a
different separator.  For example, the following line will cause values to be
separated by spaces instead of commas:
::

    simulation.reporters.append(StateDataReporter('data.txt', 1000, progress=True,
            temperature=True, totalSteps=10000, separator=' '))


Saving Simulation Progress and Results
==========================================

There are three built-in ways to save the results of your simulation in OpenMM
(additional methods can be written yourself or imported through other packages
like mdtraj or parmed). If you are simply interested in saving the structure,
you can write it out as a PDB file using :code:`PDBFile.writeFile()`.  You can
see an example of this in the modeller section :numref:`saving-the-results`.

If you are hoping to save more information than just positions, you can use
:code:`simulation.saveState()`. This will save the entire state of the
simulation, including positions, velocities, box dimensions and much more in an
XML file. This same file can be loaded back into OpenMM and used to continue
the simulation. Importantly, because this file is a text file, it can be
transfered between different platforms and different versions of OpenMM. A
potential downside to this approach is that state files are often quite large,
and may not fit all use cases. Here's an example of how to use it:
::

    simulation.saveState('output.xml')

To load the simulation back in:
::

    simulation.loadState('output.xml')

There is a third way to save your simulation, known as a checkpoint file, which
will save the entire simulation as a binary file. It will allow you to exactly
continue a simulation if the need arises (though whether the simulation is
deterministic depends on platform and methods, see
:numref:`platform-specific-properties-determinism`). There are important caveats
to this approach, however. This binary can only be used to restart simulations
on machines with the same hardware and the same OpenMM version as the one that
saved it. Therefore, it should only be used when it's clear that won't be an
issue.
::

    simulation.saveCheckpoint('state.chk')

And can be loaded back in like this:
::

    simulation.loadCheckpoint('state.chk')

Finally, OpenMM comes with a built-in reporter for saving checkpoints, the
:class:`CheckpointReporter`, which can be helpful in restarting simulations
that failed unexpectedly or due to outside reasons (e.g. server crash). To save
a checkpoint file every 5,000 steps, for example:
::

    simulation.reporters.append(CheckpointReporter('checkpnt.chk', 5000))

Note that the checkpoint reporter will overwrite the last checkpoint file.


Enhanced Sampling Methods
=========================

In many situations, the goal of a simulation is to sample the range of configurations
accessible to a system.  It does not matter whether the simulation represents a
single, physically realistic trajectory, only whether it produces a correct distribution
of states.  In this case, a variety of methods can be used to sample configuration
space much more quickly and efficiently than a single physical trajectory would.
These are known as enhanced sampling methods.  OpenMM offers several that you
can choose from.  They are briefly described here.  Consult the API documentation
for more detailed descriptions and example code.

Simulated Tempering
-------------------

Simulated tempering\ :cite:`Marinari1992` involves making frequent changes to the
temperature of a simulation.  At high temperatures, it can quickly cross energy barriers
and explore a wide range of configurations.  At lower temperatures, it more thoroughly
explores each local region of configuration space.  This is a powerful method to
speed up sampling when you do not know in advance what motions you want to sample.
Simply specify the range of temperatures to simulate and the algorithm handles
everything for you mostly automatically.

Metadynamics
------------

Metadynamics\ :cite:`Barducci2008` is used when you do know in advance what
motions you want to sample.  You specify one or more collective variables, and the
algorithm adds a biasing potential to make the simulation explore a wide range of
values for those variables.  It does this by periodically adding "bumps" to the biasing
potential at the current values of the collective variables.  This encourages the simulation
to move away from regions it has already explored and sample a wide range of values.
At the end of the simulation, the biasing potential can be used to calculate the
free energy of the system as a function of the collective variables.

Accelerated Molecular Dynamics (aMD)
------------------------------------

aMD\ :cite:`Hamelberg2007` is another method that can be used when you do not know in
advance what motions you want to accelerate.  It alters the potential energy surface
by adding a "boost" potential whenever the potential energy is below a threshold.
This makes local minima shallower and allows more frequent transitions between them.
The boost can be applied to the total potential energy, to just a subset of interactions
(typically the dihedral torsions), or both.  There are separate integrator classes
for each of these options: :class:`AMDIntegrator`, :class:`AMDForceGroupIntegrator`,
and :class:`DualAMDIntegrator`.
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-core-library:

The Core Library
################

OpenMM is based on a layered architecture, as shown in the following diagram:

.. figure:: ../images/ArchitectureLayers.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,Architecture Layers`\ : OpenMM architecture

The public API layer consists of the classes you access when using OpenMM in an
application: System; Force and its subclasses; Integrator and its subclasses;
and Context.  These classes define a public interface but do no computation.

The next layer down consists of “implementation” classes that mirror the public
API classes: ContextImpl, ForceImpl, and a subclass of ForceImpl for each
subclass of Force (HarmonicBondForceImpl, NonbondedForceImpl, etc.).  These
objects are created automatically when you create a Context.  They store
information related to a particular simulation, and define methods for
performing calculations.

Note that, whereas a Force is logically “part of” a System, a ForceImpl is
logically “part of” a Context.  (See :autonumref:`Figure,API Relationships`\ .)  If you create many Contexts
for simulating the same System, there is still only one System and only one copy
of each Force in it.  But there will be separate ForceImpls for each Context,
and those ForceImpls store information related to their particular Contexts.


.. figure:: ../images/SystemContextRelationships.jpg
   :align: center

   :autonumber:`Figure,API Relationships`\ : Relationships between public API and implementation layer objects

Also note that there is no “IntegratorImpl” class, because it is not needed.
Integrator is already specific to one Context.  Many Contexts can all simulate
the same System, but each of them must have its own Integrator, so information
specific to one simulation can be stored directly in the Integrator.

The next layer down is the OpenMM Low Level API (OLLA).  The important classes
in this layer are: Platform; Kernel; KernelImpl and its subclasses; and
KernelFactory.  A Kernel is just a reference counted pointer to a KernelImpl;
the real work is done by KernelImpl objects (or more precisely, by instances of
its subclasses).  A KernelFactory creates KernelImpl objects, and a Platform
ties together a set of KernelFactories, as well as defining information that
applies generally to performing computations with that Platform.

All of these classes (except Kernel) are abstract.  A particular Platform
provides concrete subclasses of all of them.  For example, the reference
platform defines a Platform subclass called ReferencePlatform, a KernelFactory
subclass called ReferenceKernelFactory, and a concrete subclass of each abstract
KernelImpl type: ReferenceCalcNonbondedForceKernel extends
CalcNonbondedForceKernel (which in turn extends KernelImpl),
ReferenceIntegrateVerletStepKernel extends IntegrateVerletStepKernel, and so on.

We can understand this better by walking through the entire sequence of events
that takes place when you create a Context.  As an example, suppose you create a
System; add a NonbondedForce to it; create a VerletIntegrator; and then create a
Context for them using the reference Platform.  Here is what happens.

#. The Context constructor creates a ContextImpl.
#. The ContextImpl calls :code:`createImpl()` on each Force in the System,
   which creates an instance of the appropriate ForceImpl subclass.
#. The ContextImpl calls :code:`contextCreated()` on the Platform(), which
   in turn calls :code:`setPlatformData()` on the ContextImpl.  This allows
   Platform-specific information to be stored in a ContextImpl.  Every Platform has
   its own mechanism for storing particle masses, constraint definitions, particle
   positions, and so on.  ContextImpl therefore allows the Platform to create an
   arbitrary block of data and store it where it can be accessed by that Platform’s
   kernels.
#. The ContextImpl  calls :code:`createKernel()` on the Platform several
   times to get instances of various kernels that it needs:
   CalcKineticEnergyKernel, ApplyConstraintsKernel, etc.

   #. For each kernel, the Platform looks up which KernelFactory has been
      registered for that particular kernel.  In this case, it will be a
      ReferenceKernelFactory.
   #. It calls :code:`createKernelImpl()` on the KernelFactory, which
      creates and returns an instance of an appropriate KernelImpl subclass:
      ReferenceCalcKineticEnergyKernel, ReferenceApplyConstraintsKernel, etc.

#. The ContextImpl loops over all of its ForceImpls and calls
   :code:`initialize()` on each one.

   #. Each ForceImpl asks the Platform to create whatever kernels it needs.  In
      this example, NonbondedForceImpl will request a CalcNonbondedForceKernel, and
      get back a ReferenceCalcNonbondedForceKernel.

#. The ContextImpl calls :code:`initialize()` on the Integrator which, like
   the other objects, requests kernels from the Platform.  In this example,
   VerletIntegrator requests an IntegrateVerletStepKernel and gets back a
   ReferenceIntegrateVerletStepKernel.


At this point, the Context is fully initialized and ready for doing computation.
Reference implementations of various KernelImpls have been created, but they are
always referenced through abstract superclasses.  Similarly, data structures
specific to the reference Platform have been created and stored in the
ContextImpl, but the format and content of these structures is opaque to the
ContextImpl.  Whenever it needs to access them (for example, to get or set
particle positions), it does so through a kernel (UpdateStateDataKernel in this
case).

Now suppose that you call :code:`step()` on the VerletIntegrator.  Here is
what happens to execute each time step.

#. The VerletIntegrator calls :code:`updateContextState()` on the
   ContextImpl.  This gives each Force an opportunity to modify the state of the
   Context at the start of each time step.

   #. The ContextImpl loops over its ForceImpls and calls
      :code:`updateContextState()` on each one.  In this case, our only ForceImpl is
      a NonbondedForceImpl, which returns without doing anything.  On the other hand,
      if we had an AndersenThermostat in our System, its ForceImpl would invoke a
      kernel to modify particle velocities.

#. The VerletIntegrator calls :code:`calcForcesAndEnergy()` on the
   ContextImpl to request that the forces be computed.

   #. The ContextImpl calls :code:`beginComputation()` on its
      CalcForcesAndEnergyKernel.  This initializes all the forces to zero and does any
      other initialization the Platform requires before forces can be computed.  For
      example, some Platforms construct their nonbonded neighbor lists at this point.
   #. The ContextImpl loops over its ForceImpls and calls
      :code:`calcForcesAndEnergy()` on each one.  In this case, we have a
      NonbondedForceImpl which invokes its CalcNonbondedForceKernel to compute forces.
   #. Finally, the ContextImpl calls :code:`finishComputation()` on its
      CalcForcesAndEnergyKernel.  This does any additional work needed to determine
      the final forces, such as summing the values from intermediate buffers.

#. Finally, the VerletIntegrator invokes its IntegrateVerletStepKernel.  This
   takes the forces, positions, and velocities that are stored in a Platform-
   specific format in the ContextImpl, uses them to compute new positions and
   velocities, and stores them in the ContextImpl.

.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _writing-plugins:

Writing Plugins
###############

A plugin is a dynamic library that adds new features to OpenMM.  It is typically
stored in the :code:`lib/plugins` directory inside your OpenMM installation,
and gets loaded along with all other plugins when the user calls
::

    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

It is also possible to load plugins from a different directory, or to load them
individually by calling :code:`Platform::loadPluginLibrary()`\ .

Every plugin must implement two functions that are declared in the
PluginInitializer.h header file:
::

    extern "C" void registerPlatforms();
    extern "C" void registerKernelFactories();

When a plugin is loaded, these two functions are invoked to register any
Platforms and KernelFactories defined by the plugin.  When many plugins are
loaded at once by calling :code:`Platform::loadPluginsFromDirectory()`\ ,
:code:`registerPlatforms()` is first called on all of them, then
:code:`registerKernelFactories()` is called on all of them.  This allows one
plugin to define a Platform, and a different plugin to add KernelFactories to
it; the Platform is guaranteed to be registered by the first plugin before the
second plugin tries to add its KernelFactories, regardless of what order the
plugins happen to be loaded in.

Creating New Platforms
**********************

One common type of plugin defines a new Platform.  There are four such plugins
that come with OpenMM: one for the Reference platform, one for the CPU Platform,
one for the CUDA Platform, and one for the OpenCL Platform.

To define a new Platform, you must create subclasses of the various abstract
classes in the OpenMM Low Level API: a subclass of Platform, one or more
subclasses of KernelFactory, and a subclass of each KernelImpl.  That is easy to
say, but a huge amount of work to actually do.  There are many different
algorithms involved in computing forces, enforcing constraints, performing
integration, and so on, all of which together make up a Platform.  Of course,
there is no requirement that every Platform must implement every possible
feature.  If you do not provide an implementation of a particular kernel, it
simply means your Platform cannot be used for any simulation that requires that
kernel; if a user tries to do so, an exception will be thrown.

Your plugin’s :code:`registerPlatforms()` function should create an instance
of your Platform subclass, then register it by calling
:code:`Platform::registerPlatform()`\ .  You also must register the
KernelFactory for each kernel your Platform supports.  This can be done in the
:code:`registerKernelFactories()` function, or more simply, directly in the
Platform’s constructor.  You can use as many different KernelFactories as you
want for different kernels, but usually it is simplest to use a single
KernelFactory for all of them.  The support for multiple KernelFactories exists
primarily to let plugins add new features to existing Platforms, as described in
the next section.

Creating New Forces
*******************

Another common type of plugin defines new Forces and provides implementations of
them for existing Platforms.  (Defining new Integrators is not specifically
discussed here, but the process is very similar.)  There are two such plugins
that come with OpenMM.  They implement the AMOEBA force field and Drude
oscillators, respectively.

As an example, suppose you want to create a new Force subclass called
StringForce that uses the equations of String Theory to compute the interactions
between particles.  You want to provide implementations of it for all four
standard platforms: Reference, CPU, CUDA, and OpenCL.

The first thing to realize is that this *cannot* be done with only a plugin
library.  Plugins are loaded dynamically at runtime, and they relate to the low
level API; but you must also provide a public API.  Users of your class need to
create StringForce objects and call methods on them.  That means providing a
header file with the class declaration, and a (non-plugin) library with the
class definition to link their code against.  The implementations for particular
Platforms can be in plugins, but the public API class itself cannot.  Or to put
it differently, the full “plugin” (from the user’s perspective) consists of
three parts: the library OpenMM loads at runtime (which is what OpenMM considers
to be the “plugin”), a second library for users to link their code against, and
a header file for them to include in their source code.

To define the API, you will need to create the following classes:

#. StringForce.  This is the public API for your force, and users will directly
   link against the library containing it.
#. StringForceImpl.  This is the ForceImpl subclass corresponding to
   StringForce.  It should be defined in the same library as StringForce, and
   StringForce’s :code:`createImpl()` method should create an instance of it.
#. CalcStringForceKernel.  This is an abstract class that extends KernelImpl,
   and defines the API by which StringForceImpl invokes its kernel.  You only need
   to provide a header file for it, not an implementation; those will be provided
   by Platforms.


Now suppose you are writing the OpenCL implementation of StringForce.  Here are
the classes you need to write:

#. OpenCLCalcStringForceKernel.  This extends CalcStringForceKernel and provides
   implementations of its virtual methods.  The code for this class will probably
   be very complicated (and if it actually works, worth a Nobel Prize).  It may
   execute many different GPU kernels and create its own internal data structures.
   But those details are entirely internal to your own code.  As long as this class
   implements the virtual methods of CalcStringForceKernel, you can do anything you
   want inside it.
#. OpenCLStringForceKernelFactory.  This is a KernelFactory subclass that knows
   how to create instances of OpenCLCalcStringForceKernel.


Both of these classes should be packaged into a dynamic library (.so on Linux,
.dylib on Mac, .dll on Windows) that can be loaded as a plugin.  This library
must also implement the two functions from PluginInitializer.h.
:code:`registerPlatforms()` will do nothing, since this plugin does not
implement any new Platforms.  :code:`registerKernelFactories()` should call
\ :code:`Platform::getPlatformByName("OpenCL")` to get the OpenCL Platform,
then create a new OpenCLStringForceKernelFactory and call
:code:`registerKernelFactory()` on the Platform to register it.  If the OpenCL
Platform is not available, you should catch the exception then return without
doing anything.  Most likely this means there is no OpenCL runtime on the
computer your code is running on.
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++


.. _the-reference-platform:

The Reference Platform
######################

The reference Platform is written with simplicity and clarity in mind, not
performance.  (It is still not always as simple or clear as one might hope, but
that is the goal.)  When implementing a new feature, it is recommended to create
the reference implementation first, then use that as a model for the versions in
other Platforms.

When using the reference Platform, the “platform-specific data” stored in
ContextImpl is of type ReferencePlatform::PlatformData, which is declared in
ReferencePlatform.h.  It has fields for storing positions, velocities, box
vectors, and other types of data.

The PlatformData’s vector of forces contains one element for each particle.  At
the start of each force evaluation, all elements of it are set to zero.  Each
Force adds its own contributions to the vector, so that at the end, it contains
the total force acting on each particle.

There are a few additional classes that contain useful static methods.
SimTKOpenMMUtilities has various utility functions, of which the most important
is a random number generator.  ReferenceForce provides methods for calculating
the displacement between two positions, optionally taking periodic boundary
conditions into account.

.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++


.. _the-cpu-platform:

The CPU Plaform
###############

CpuPlatform is a subclass of ReferencePlatform.  It provides optimized versions
of a small number of kernels, while using the reference implementations for all
the others.  Any kernel implementation written for the reference Platform will
work equally well with the CPU platform.  Of course, if that kernel happens to
be a performance bottleneck, you will probably want to write an optimized
version of it.  But many kernels have negligible effect on performance, and for
these you can just use the same implementation for both platforms.

If you choose to do that, you can easily support both platforms with a single
plugin library.  Just implement :code:`registerKernelFactories()` like this:
::

    extern "C" void registerKernelFactories() {
        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
            Platform& platform = Platform::getPlatform(i);
            if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
                // Create and register your KernelFactory.
            }
        }
    }

The loop identifies every ReferencePlatform, either an instance of the base
class or of a subclass, and registers a KernelFactory for every one.
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-cuda-platform:

The CUDA Platform
#################

The CUDA platform is very similar to the OpenCL platform, and most of the
previous chapter applies equally well to it, just changing “OpenCL” to “Cuda” in
class names.  There are a few differences worth noting.

Compiling Kernels
*****************

Like the OpenCL platform, the CUDA platform compiles all its kernels at runtime.
Unlike OpenCL, CUDA does not have built in support for runtime compilation.
OpenMM therefore needs to implement this itself by writing the source code out
to disk, invoking the nvcc compiler as a separate process, and then loading the
compiled kernel in from disk.

For the most part, you can ignore all of this.  Just call
:code:`createModule()` on the CudaContext, passing it the CUDA source code.
It takes care of the details of compilation and loading, returning a CUmodule
object when it is done.  You can then call :code:`getKernel()` to look up
individual kernels in the module (represented as CUfunction objects) and
:code:`executeKernel()` to execute them.

The CUDA platform does need two things to make this work: a directory on disk
where it can write out temporary files, and the path to the nvcc compiler.
These are specified by the “CudaTempDirectory” and “CudaCompiler” properties
when you create a new Context.  It often can figure out suitable values for them
on its own, but sometimes it needs help.  See the “Platform-Specific Properties”
chapter of the User's Manual for details.

Accumulating Forces
*******************

The OpenCL platform, as described in Section :numref:`computing-forces`\ , uses two types of buffers for
accumulating forces: a set of floating point buffers, and a single fixed point
buffer.  In contrast, the CUDA platform uses *only* the fixed point buffer
(represented by the CUDA type :code:`long` :code:`long`\ ).  This means
the CUDA platform only works on devices that support 64 bit atomic operations
(compute capability 1.2 or higher).
########################
OpenMM Developer's Guide
########################

.. only:: latex

   .. include:: license.rst

.. toctree::
   :maxdepth: 3
   :numbered:
   
   01_introduction.rst
   02_core_library.rst
   03_writing_plugins.rst
   04_reference_platform.rst
   05_cpu_platform.rst
   06_opencl_platform.rst
   07_cuda_platform.rst
   08_common_compute.rst

.. only:: html

   .. include:: license.rst
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-opencl-platform:

The OpenCL Platform
###################

The OpenCL Platform is much more complicated than the reference Platform.  It
also provides many more tools to simplify your work, but those tools themselves
can be complicated to use correctly.  This chapter will attempt to explain how
to use some of the most important ones.  It will *not* teach you how to
program with OpenCL.  There are many tutorials on that subject available
elsewhere, and this guide assumes you already understand it.

Overview
********

When using the OpenCL Platform, the “platform-specific data” stored in
ContextImpl is of type OpenCLPlatform::PlatformData, which is declared in
OpenCLPlatform.h.  The most important field of this class is :code:`contexts`
, which is a vector of OpenCLContexts.  (There is one OpenCLContext for each
device you are using.  The most common case is that you are running everything
on a single device, in which case there will be only one OpenCLContext.
Parallelizing computations across multiple devices is not discussed here.)  The
OpenCLContext stores most of the important information about a simulation:
positions, velocities, forces, an OpenCL CommandQueue used for executing
kernels, workspace buffers of various sorts, etc.  It provides many useful
methods for compiling and executing kernels, clearing and reducing buffers, and
so on.  It also provides access to three other important objects: the
OpenCLIntegrationUtilities, OpenCLNonbondedUtilities, and OpenCLBondedUtilities.
These are discussed below.

Allocation of device memory is generally done through the OpenCLArray class.  It
takes care of much of the work of memory management, and provides a simple
interface for transferring data between host and device memory.

Every kernel is specific to a particular OpenCLContext, which in turn is
specific to a particular OpenMM::Context.  This means that kernel source code
can be customized for a particular simulation.  For example, values such as the
number of particles can be turned into compile-time constants, and specific
versions of kernels can be selected based on the device being used or on
particular aspects of the system being simulated.
:code:`OpenCLContext::createProgram()` makes it easy to specify a list of
preprocessor definitions to use when compiling a kernel.

The normal way to execute a kernel is by calling :code:`executeKernel()` on
the OpenCLContext.  It allows you to specify the total number of work-items to
execute, and optionally the size of each work-group.  (If you do not specify a
work-group size, it uses 64 as a default.)  The number of work-groups to launch
is selected automatically based on the work-group size, the total number of
work-items, and the number of compute units in the device it will execute on.

Numerical Precision
*******************

The OpenCL platform supports three precision modes:

#. **Single**\ : All values are stored in single precision, and nearly all
   calculations are done in single precision.  The arrays of positions, velocities,
   forces, and energies (returned by the OpenCLContext’s :code:`getPosq()`\ ,
   :code:`getVelm()`\ , :code:`getForce()`\ , :code:`getForceBuffers()`\ , and
   :code:`getEnergyBuffer()` methods) are all of type :code:`float4` (or
   :code:`float` in the case of :code:`getEnergyBuffer()`\ ).
#. **Mixed**\ : Forces are computed and stored in single precision, but
   integration is done in double precision.  The velocities have type
   :code:`double4`\ .  The positions are still stored in single precision to avoid
   adding overhead to the force calculations, but a second array of type
   :code:`float4` is created to store “corrections” to the positions (returned by
   the OpenCLContext’s getPosqCorrection() method).  Adding the position and the
   correction together gives the full double precision position.
#. **Double**\ : Positions, velocities, forces, and energies are all stored in
   double precision, and nearly all calculations are done in double precision.


You can call :code:`getUseMixedPrecision()` and
:code:`getUseDoublePrecision()` on the OpenCLContext to determine which mode
is being used.  In addition, when you compile a kernel by calling
:code:`createKernel()`\ , it automatically defines two types for you to make it
easier to write kernels that work in any mode:

#. :code:`real` is defined as :code:`float` in single or mixed precision
   mode, :code:`double` in double precision mode.
#. :code:`mixed` is defined as :code:`float` in single precision mode,
   :code:`double` in mixed or double precision mode.


It also defines vector versions of these types (\ :code:`real2`\ ,
:code:`real4`\ , etc.).

.. _computing-forces:

Computing Forces
****************

When forces are computed, they can be stored in either of two places.  There is
an array of :code:`long` values storing them as 64 bit fixed point values, and
a collection of buffers of :code:`real4` values storing them in floating point
format.  Most GPUs support atomic operations on 64 bit integers, which allows
many threads to simultaneously record forces without a danger of conflicts.
Some low end GPUs do not support this, however, especially the embedded GPUs
found in many laptops.  These devices write to the floating point buffers, with
careful coordination to make sure two threads will never write to the same
memory location at the same time.

At the start of a force calculation, all forces in all buffers are set to zero.
Each Force is then free to add its contributions to any or all of the buffers.
Finally, the buffers are summed to produce the total force on each particle.
The total is recorded in both the floating point and fixed point arrays.

The size of each floating point buffer is equal to the number of particles, rounded up to the
next multiple of 32.  Call :code:`getPaddedNumAtoms()` on the OpenCLContext
to get that number.  The actual force buffers are obtained by calling
:code:`getForceBuffers()`\ .  The first *n* entries (where *n* is the
padded number of atoms) represent the first force buffer, the next *n*
represent the second force buffer, and so on.  More generally, the *i*\ ’th
force buffer’s contribution to the force on particle *j* is stored in
element :code:`i*context.getPaddedNumAtoms()+j`\ .

The fixed point buffer is ordered differently.  For atom *i*\ , the x component
of its force is stored in element :code:`i`\ , the y component in element
:code:`i+context.getPaddedNumAtoms()`\ , and the z component in element
:code:`i+2*context.getPaddedNumAtoms()`\ .  To convert a value from floating
point to fixed point, multiply it by 0x100000000 (2\ :sup:`32`\ ),
then cast it to a :code:`long`\ .  Call :code:`getLongForceBuffer()` to get the
array of fixed point values.

The potential energy is also accumulated in a set of buffers, but this one is
simply a list of floating point values.  All of them are set to zero at the
start of a computation, and they are summed at the end of the computation to
yield the total energy.

The OpenCL implementation of each Force object should define a subclass of
ComputeForceInfo, and register an instance of it by calling :code:`addForce()` on
the OpenCLContext.  It implements methods for determining whether particular
particles or groups of particles are identical.  This is important when
reordering particles, and is discussed below.


Nonbonded Forces
****************

Computing nonbonded interactions efficiently is a complicated business in the
best of cases.  It is even more complicated on a GPU.  Furthermore, the
algorithms must vary based on the type of processor being used, whether there is
a distance cutoff, and whether periodic boundary conditions are being applied.

The OpenCLNonbondedUtilities class tries to simplify all of this.  To use it you
need provide only a piece of code to compute the interaction between two
particles.  It then takes responsibility for generating a neighbor list, looping
over interacting particles, loading particle parameters from global memory, and
writing the forces and energies to the appropriate buffers.  All of these things
are done using an algorithm appropriate to the processor you are running on and
high level aspects of the interaction, such as whether it uses a cutoff and
whether particular particle pairs need to be excluded.

Of course, this system relies on certain assumptions, the most important of
which is that the Force can be represented as a sum of independent pairwise
interactions.  If that is not the case, things become much more complicated.
You may still be able to use features of OpenCLNonbondedUtilities, but you
cannot use the simple mechanism outlined above.  That is beyond the scope of
this guide.

To define a nonbonded interaction, call :code:`addInteraction()` on the
OpenCLNonbondedUtilities, providing a block of OpenCL source code for computing
the interaction.  This block of source code will be inserted into the middle of
an appropriate kernel.  At the point where it is inserted, various variables
will have been defined describing the interaction to compute:

#. :code:`atom1` and :code:`atom2` are the indices of the two
   interacting particles.
#. :code:`r`\ , :code:`r2`\ , and :code:`invR` are the distance *r*
   between the two particles, *r*\ :sup:`2`\ , and 1/\ *r* respectively.
#. :code:`isExcluded` is a :code:`bool` specifying whether this pair of
   particles is marked as an excluded interaction.  (Excluded pairs are not skipped
   automatically, because in some cases they still need to be processed, just
   differently from other pairs.)
#. :code:`posq1` and :code:`posq2` are :code:`real4`\ s containing the
   positions (in the xyz fields) and charges (in the w fields) of the two
   particles.
#. Other per-particle parameters may be specified, as described below.


The following preprocessor macros will also have been defined:

#. :code:`NUM_ATOMS` is the total number of particles in the system.
#. :code:`PADDED_NUM_ATOMS` is the padded number of particles in the system.
#. :code:`USE_CUTOFF` is defined if and only if a cutoff is being used
#. :code:`USE_PERIODIC` is defined if and only if periodic boundary
   conditions are being used.
#. :code:`CUTOFF` and :code:`CUTOFF_SQUARED` are the cutoff distance and
   its square respectively (but only defined if a cutoff is being used).


Finally, two output variables will have been defined:

#. You should add the energy of the interaction to :code:`tempEnergy`\ .
#. You should add the derivative of the energy with respect to the inter-particle
   distance to :code:`dEdR`\ .


You can also define arbitrary per-particle parameters by calling
:code:`addParameter()` on the OpenCLNonbondedUtilities.  You provide an array
in device memory containing the set of values, and the values for the two
interacting particles will be loaded and stored into variables called
:code:`<name>1` and :code:`<name>2`\ , where <name> is the name you specify
for the parameter.  Note that nonbonded interactions are not computed until
after :code:`calcForcesAndEnergy()` has been called on every ForceImpl, so
it is possible to make the parameter values change with time by modifying them
inside :code:`calcForcesAndEnergy()`\ .  Also note that the length of the
array containing the parameter values must equal the *padded* number of
particles in the system.

Finally, you can specify arbitrary other memory objects that should be passed as
arguments to the interaction kernel by calling :code:`addArgument()`\ .  The
rest of the kernel ignores these arguments, but you can make use of them in your
interaction code.

Consider a simple example.  Suppose we want to implement a nonbonded interaction
of the form *E*\ =\ *k*\ :sub:`1`\ *k*\ :sub:`2`\ *r*\ :sup:`2`\ ,
where *k* is a per-particle parameter.  First we create a parameter as
follows
::

    nb.addParameter(ComputeParameterInfo(kparam, "kparam", "float", 1));

where :code:`nb` is the OpenCLNonbondedUtilities for the context.  Now we
call :code:`addInteraction()` to define an interaction with the following
source code:
::

    #ifdef USE_CUTOFF
    if (!isExcluded && r2 < CUTOFF_SQUARED) {
    #else
    if (!isExcluded) {
    #endif
        tempEnergy += kparam1*kparam2*r2;
        dEdR += 2*kparam1*kparam2*r;
    }

An important point is that this code is executed for every pair of particles in
the *padded* list of atoms.  This means that some interactions involve
padding atoms, and should not actually be included.  You might think, then, that
the above code is incorrect and we need another check to filter out the extra
interactions:
::

    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)

This is not necessary in our case, because the :code:`isExcluded` flag is
always set for interactions that involve a padding atom.  If our force did not
use excluded interactions (and so did not check :code:`isExcluded`\ ), then we
would need to add this extra check.  Self interactions are a similar case: we do
not check for :code:`(atom1 == atom2)` because the exclusion flag prevents
them from being processed, but for some forces that check is necessary.

Bonded Forces
*************

Just as OpenCLNonbondedUtilities simplifies the task of creating nonbonded
interactions, OpenCLBondedUtilities simplifies the process for many types of
bonded interactions.  A “bonded interaction” means one that is applied to small,
fixed groups of particles.  This includes bonds, angles, torsions, etc.  The
important point is that the list of particles forming a “bond” is known in
advance and does not change with time.

Using OpenCLBondedUtilities is very similar to the process described above.  You
provide a block of OpenCL code for evaluating a single interaction.  This block
of code will be inserted into the middle of a kernel that loops over all
interactions and evaluates each one.  At the point where it is inserted, the
following variables will have been defined describing the interaction to
compute:

#. :code:`index` is the index of the interaction being evaluated.
#. :code:`atom1`\ , :code:`atom2`\ , ... are the indices of the interacting
   particles.
#. :code:`pos1`\ , :code:`pos2`\ , ... are :code:`real4`\ s containing the
   positions (in the xyz fields) of the interacting particles.


A variable called :code:`energy` will have been defined for accumulating the
total energy of all interactions.  Your code should add the energy of the
interaction to it.  You also should define :code:`real4` variables called
:code:`force1`\ , :code:`force2`\ , ... and store the force on each atom into
them.

As a simple example, the following source code implements a pairwise interaction
of the form *E*\ =\ *r*\ :sup:`2`\ :
::

    real4 delta = pos2-pos1;
    energy += delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    real4 force1 = 2.0f*delta;
    real4 force2 = -2.0f*delta;

To use it, call :code:`addInteraction()` on the Context’s
OpenCLBondedUtilities object.  You also provide a list of the particles involved
in every bonded interaction.

Exactly as with nonbonded interactions, you can call :code:`addArgument()`
to specify arbitrary memory objects that should be passed as arguments to the
interaction kernel.  These might contain per-bond parameters (use
:code:`index` to look up the appropriate element) or any other information you
want.

Reordering of Particles
***********************

Nonbonded calculations are done a bit differently in the OpenCL Platform than in
most CPU based codes.  In particular, interactions are computed on blocks of 32
particles at a time (which is why the number of particles needs to be padded to
bring it up to a multiple of 32), and the neighbor list actually lists pairs of
\ *blocks*\ , not pairs of individual particles, that are close enough to
interact with each other.

This only works well if sequential particles tend to be close together so that
blocks are spatially compact.  This is generally true of particles in a
macromolecule, but it is not true for solvent molecules.  Each water molecule,
for example, can move independently of other water molecules, so particles that
happen to be sequential in whatever order the molecules were defined in need not
be spatially close together.

The OpenCL Platform addresses this by periodically reordering particles so that
sequential particles are close together.  This means that what the OpenCL
Platform calls particle *i* need not be the same as what the System calls
particle *i*\ .

This reordering is done frequently, so it must be very fast.  If all the data
structures describing the structure of the System and the Forces acting on it
needed to be updated, that would make it prohibitively slow.  The OpenCL
Platform therefore only reorders particles in ways that do not alter any part of
the System definition.  In practice, this means exchanging entire molecules; as
long as two molecules are truly identical, their positions and velocities can be
exchanged without affecting the System in any way.

Every Force can contribute to defining the boundaries of molecules, and to
determining whether two molecules are identical.  This is done through the
ComputeForceInfo it adds to the OpenCLContext.  It can specify two types of
information:

#. Given a pair of particles, it can say whether those two particles are
   identical (as far as that Force is concerned).  For example, a Force object
   implementing a Coulomb force would check whether the two particles had equal
   charges.
#. It can define *particle groups*\ .  The OpenCL Platform will ensure that
   all the particles in a group are part of the same molecule.  It also can specify
   whether two groups are identical to each other.  For example, in a Force
   implementing harmonic bonds, each group would consist of the two particles
   connected by a bond, and two groups would be identical if they had the same
   spring constants and equilibrium lengths.


Integration Utilities
*********************

The OpenCLContext’s OpenCLIntegrationUtilities provides features that are used
by many integrators.  The two most important are random number generation and
constraint enforcement.

If you plan to use random numbers, you should call
:code:`initRandomNumberGenerator()` during initialization, specifying the
random number seed to use.  Be aware that there is only one random number
generator, even if multiple classes make use of it.  If two classes each call
:code:`initRandomNumberGenerator()` and request different seeds, an exception
will be thrown.  If they each request the same seed, the second call will simply
be ignored.

For efficiency, random numbers are generated in bulk and stored in an array in
device memory, which you can access by calling :code:`getRandom()`\ .  Each
time you need to use a block of random numbers, call
:code:`prepareRandomNumbers()`\ , specifying how many values you need.  It will
register that many values as having been used, and return the index in the array
at which you should start reading values.  If not enough unused values remain in
the array, it will generate a new batch of random values before returning.

To apply constraints, simply call :code:`applyConstraints()`\ .  For numerical
accuracy, the constraint algorithms do not work on particle positions directly,
but rather on the *displacements* taken by the most recent integration step.
These displacements must be stored in an array which you can get by calling
:code:`getPosDelta()`\ .  That is, the constraint algorithms assume the actual
(unconstrained) position of each particle equals the position stored in the
OpenCLContext plus the delta stored in the OpenCLIntegrationUtilities.  It then
modifies the deltas so that all distance constraints are satisfied.  The
integrator must then finish the time step by adding the deltas to the positions
and storing them into the main position array.
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _common-compute:

Common Compute
##############

Common Compute is not a platform, but it shares many elements of one.  It exists
to reduce code duplication between the OpenCL and CUDA platforms.  It allows a
single implementation to be written for most kernels that can be used by both
platforms.

OpenCL and CUDA are very similar to each other.  Their computational models are
nearly identical.  For example, each is based around launching kernels that are
executed in parallel by many threads.  Each of them groups threads into blocks,
with more communication and synchronization permitted between the threads
in a block than between ones in different blocks.  They have very similar memory
hierarchies: high latency global memory, low latency local/shared memory that
can be used for communication between the threads of a block, and local variables
that are visible only to a single thread.

Even their languages for writing kernels are very similar.  Here is an OpenCL
kernel that adds two arrays together, storing the result in a third array.
::

    __kernel void addArrays(__global const float* restrict a,
                            __global const float* restrict b,
                            __global float* restrict c
                            int length) {
        for (int i = get_global_id(0); i < length; i += get_global_size(0))
            c[i] = a[i]+b[i];
    }

Here is the corresponding CUDA kernel.
::

    __extern "C" __global__ void addArrays(const float* __restrict__ a,
                                           const float* __restrict__ b,
                                           _float* __restrict__ c
                                           int length) {
        for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < length; i += blockDim.x*gridDim.x)
            c[i] = a[i]+b[i];
    }

The difference between them is largely just a mechanical find-and-replace.
After many years of writing and maintaining nearly identical kernels by hand,
it finally occurred to us that the translation could be done automatically by
the compiler.  Simply by defining a few preprocessor macros, the following
kernel can be compiled equally well either as OpenCL or as CUDA.
::

    KERNEL void addArrays(GLOBAL const float* RESTRICT a,
                          GLOBAL const float* RESTRICT b,
                          GLOBAL float* RESTRICT c
                          int length) {
        for (int i = GLOBAL_ID; i < length; i += GLOBAL_SIZE)
            c[i] = a[i]+b[i];
    }

Writing Device Code
*******************

When compiling kernels with the Common Compute API, the following macros are
defined.

+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|Macro                          |OpenCL Definition                                           |CUDA Definition                             |
+===============================+============================================================+============================================+
|:code:`KERNEL`                 |:code:`__kernel`                                            |:code:`extern "C" __global__`               |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`DEVICE`                 |                                                            |:code:`__device__`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL`                  |:code:`__local`                                             |:code:`__shared__`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_ARG`              |:code:`__local`                                             |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL`                 |:code:`__global`                                            |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`RESTRICT`               |:code:`restrict`                                            |:code:`__restrict__`                        |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_ID`               |:code:`get_local_id(0)`                                     |:code:`threadIdx.x`                         |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_SIZE`             |:code:`get_local_size(0)`                                   |:code:`blockDim.x`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL_ID`              |:code:`get_global_id(0)`                                    |:code:`(blockIdx.x*blockDim.x+threadIdx.x)` |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL_SIZE`            |:code:`get_global_size(0)`                                  |:code:`(blockDim.x*gridDim.x)`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GROUP_ID`               |:code:`get_group_id(0)`                                     |:code:`blockIdx.x`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`NUM_GROUPS`             |:code:`get_num_groups(0)`                                   |:code:`gridDim.x`                           |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`SYNC_THREADS`           |:code:`barrier(CLK_LOCAL_MEM_FENCE+CLK_GLOBAL_MEM_FENCE);`  |:code:`__syncthreads();`                    |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`SYNC_WARPS`             | | if SIMT width >= 32:                                     | | if compute capability >= 7.0:            |
|                               | | :code:`mem_fence(CLK_LOCAL_MEM_FENCE)`                   | | :code:`__syncwarp();`                    |
|                               | | otherwise:                                               | | otherwise empty                          |
|                               | | :code:`barrier(CLK_LOCAL_MEM_FENCE)`                     |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`MEM_FENCE`              |:code:`mem_fence(CLK_LOCAL_MEM_FENCE+CLK_GLOBAL_MEM_FENCE);`|:code:`__threadfence_block();`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`ATOMIC_ADD(dest, value)`|:code:`atom_add(dest, value)`                               |:code:`atomicAdd(dest, value)`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+

A few other symbols may or may not be defined based on the device you are running on:
:code:`SUPPORTS_DOUBLE_PRECISION` and :code:`SUPPORTS_64_BIT_ATOMICS`\ .  You
can use :code:`#ifdef` blocks with these symbols to conditionally compile code
based on the features supported by the device.  In addition, the CUDA compiler
defines the symbol :code:`__CUDA_ARCH__`\ , so you can check for this symbol if
you want to have different code blocks for CUDA and OpenCL.

Both OpenCL and CUDA define vector types like :code:`int2` and :code:`float4`\ .
The types they support are different but overlapping.  When writing common code,
use only the vector types that are supported by both OpenCL and CUDA: 2, 3, and 4
element vectors of type :code:`short`\ , :code:`int`\ , :code:`float`\ , and
:code:`double`\ .

CUDA uses functions to construct vector values, such as :code:`make_float2(x, y)`\ .
OpenCL instead uses a typecast like syntax: :code:`(float2) (x, y)`\ .  In common
code, use the CUDA style :code:`make_` functions.  OpenMM provides definitions
of these functions when compiling as OpenCL.

In CUDA, vector types are simply data structures.  You can access their elements,
but not do much more with them.  In contrast, OpenCL's vectors are mathematical
types.  All standard math operators are defined for them, as well as geometrical
functions like :code:`dot()` and :code:`cross()`\ .  When compiling kernels as
CUDA, OpenMM provides definitions of these operators and functions.

OpenCL also supports "swizzle" notation for vectors.  For example, if :code:`f`
is a :code:`float4` you can construct a vector of its first three elements
by writing :code:`f.xyz`\ , or you can swap its first two elements by writing
:code:`f.xy = f.yx`\ .  Unfortunately, there is no practical way to support this
in CUDA, so swizzle notation cannot be used in common code.  Because stripping
the final element from a four component vector is such a common operation, OpenMM
provides a special function for doing it: :code:`trimTo3(f)` is a vector of its
first three elements.

64 bit integers are another data type that needs special handling.  Both OpenCL
and CUDA support them, but they use different names for them: :code:`long` in OpenCL,
:code:`long long` in CUDA.  To work around this inconsistency, OpenMM provides
the typedefs :code:`mm_long` and :code:`mm_ulong` for signed and unsigned 64 bit
integers in device code.

Writing Host Code
*****************

Host code for Common Compute is very similar to host code for OpenCL or CUDA.
In fact, most of the classes provided by the OpenCL and CUDA platforms are
subclasses of Common Compute classes.  For example, OpenCLContext and
CudaContext are both subclasses of ComputeContext.  When writing common code,
each KernelImpl should expect a ComputeContext to be passed to its constructor.
By using the common API provided by that abstract class, it can be used for
either OpenCL or CUDA just based on the particular context passed to it at
runtime.  Similarly, OpenCLNonbondedUtilities and CudaNonbondedUtilities are
subclasses of the abstract NonbondedUtilities class, and so on.

ArrayInterface is an abstract class defining the interface for arrays stored on
the device.  OpenCLArray and CudaArray are both subclasses of it.  To simplify
code that creates and uses arrays, there is also a third subclass called
ComputeArray.  It acts as a wrapper around an OpenCLArray or CudaArray,
automatically creating an array of the appropriate type for the current
platform.  In practice, just follow these rules:

  1. Whenever you need to create an array, make it a ComputeArray.

  2. Whenever you write a function that expects an array to be passed to it,
     declare the type to be ArrayInterface.

If you do these two things, all differences between platforms will be handled
automatically.

OpenCL and CUDA have quite different APIs for compiling and invoking kernels.
To hide these differences, OpenMM provides a set of abstract classes.  To compile
device code, pass the source code to :code:`compileProgram()` on the ComputeContext.
This returns a ComputeProgram.  You can then call its :code:`createKernel()`
method to get a ComputeKernel object, which has methods for setting arguments
and invoking the kernel.

Sometimes you need to refer to vector types in host code, such as to set the
value for a kernel argument or to access the elements of an array.  OpenCL and
CUDA both define types for them, but they have different names, and in any case
you want to avoid using OpenCL-specific or CUDA-specific types in common code.
OpenMM therefore defines types for vectors in host code.  They have the same
names as the corresponding types in device code, only with the prefix :code:`mm_`\ ,
for example :code:`mm_int2` and :code:`mm_float4`\ .

Three component vectors need special care in this context, because the platforms
define them differently.  In OpenCL, a three component vector is essentially a
four component vector whose last component is ignored.  For example,
:code:`sizeof(float3)` is 12 in CUDA but 16 in OpenCL.  Within a kernel this
distinction can usually be ignored, but when communicating between host and
device it becomes vitally important.  It is generally best to avoid storing
three component vectors in arrays or passing them as arguments.  There are no
:code:`mm_` host types defined for three component vectors, because CUDA and
OpenCL would require them to be defined in different ways.
.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

Introduction
############

This guide describes the internal architecture of the OpenMM library.  It is
targeted at developers who want to add features to OpenMM, either by modifying
the core library directly or by writing plugins.  If you just want to write
applications that use OpenMM, you do not need to read this guide; the User's
Manual tells you everything you need to know.  This guide is intended for
people who want to contribute to OpenMM itself.

It is organized as follows:

* Chapter :numref:`the-core-library` describes the architecture of the core OpenMM library.  It
  discusses how the high level and low level APIs relate to each other, and the
  flow of execution between them.
* Chapter :numref:`writing-plugins` describes in detail how to write a plugin.  It focuses on the two
  most common types of plugins: those which define new Forces, and those which
  implement new Platforms.
* Chapter :numref:`the-reference-platform` discusses the architecture of the reference Platform, providing
  information relevant to writing reference implementations of new features.
* Chapter :numref:`the-cpu-platform` discusses the architecture of the CPU Platform, providing
  information relevant to writing CPU implementations of new features.
* Chapter :numref:`the-opencl-platform` discusses the architecture of the OpenCL Platform, providing
  information relevant to writing OpenCL implementations of new features.
* Chapter :numref:`the-cuda-platform` discusses the architecture of the CUDA Platform, providing
  information relevant to writing CUDA implementations of new features.
* Chapter :numref:`common-compute` describes the Common Compute framework, which lets you
  write a single implementation of a feature that can be used for both OpenCL and CUDA.


This guide assumes you are already familiar with the public API and how to use
OpenMM in applications.  If that is not the case, you should first read the
User's Manual and work through some of the example programs.  Pay especially
close attention to the “Introduction to the OpenMM Library” chapter, since it
introduces concepts that are important in understanding this guide.


.. currentmodule:: openmm.openmm

OpenMM Python API
=================

The Python API provides information about the classes and methods available in OpenMM for Python developers.

OpenMM consists of two parts. First, there is a set of :ref:`libraries <library>` for performing many types of computations needed for molecular simulations: force evaluation, numerical integration, energy minimization, etc.

Second, there is an :ref:`application layer <app>`, a set of Python libraries providing a high level interface for running simulations. This layer is targeted at computational biologists or other people who want to run simulations, and who may or may not be programmers.

See the user guide for more details.


.. toctree::
   :maxdepth: 2

   app
   library
{{ objname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
