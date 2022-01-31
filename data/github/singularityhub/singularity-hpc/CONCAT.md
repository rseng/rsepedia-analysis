# CHANGELOG

This is a manually generated log to track changes to the repository for each release.
Each section should include general headers such as **Implemented enhancements**
and **Merged pull requests**. Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes (recipe file format? image file format?)
 - migration guidance (how to convert images?)
 - changed behaviour (recipe sections work differently)

The versions coincide with releases on pip. Only major versions will be released as tags on Github.

## [0.0.x](https://github.scom/singularityhub/singularity-hpc/tree/master) (0.0.x)
 - Central config should be read first (and updated with user config) (0.0.4)
 - Add support for oras pull for singularity (0.0.39)
 - Bug with setting a nested value (0.0.38)
 - Adding quotes around tcl descriptions (0.0.37)
 - fixing bug with container install (does not honor module directory) (0.0.36)
 - `.version` file should be written in top level of module folders [#450](https://github.com/singularityhub/singularity-hpc/issues/450) (0.0.35)
  - tcl module functions need `$@` to handle additional arguments
 - Tcl modules use shell functions for bash, to export to child shells (0.0.34)
   - fixed missing singularity -B flag for custom home feature
   - fixed singularity.tcl to always replace $ with \$ for custom home feature
   - fixed missing features.home occurrence in docker.tcl
 - New features for X11 and custom home (0.0.33)
 - Adding singularity and docker specific options (0.0.32)
   - Adding more documentation on aliases
   - Bugfix to output error message if path does not exist for inspect
 - Fixing bug with using variables in recipe options (0.0.31)
 - Removing un-used defaults and lmod_base (0.0.30)
 - Allow environment variables in settings (0.0.29)
   - User settings file creation and use with shpc config inituser
   - registry is now a list to support multiple registry locations
   - config supports add/remove to append/delete from list
 - Add test for docker and podman (0.0.28)
   - namespace as format string for command named renamed to repository
   - shpc test/uninstall should be run for all tests
   - bug with uninstall
   - adding user and group ids to docker and podman commands
   - remove conflict with entire module name, should only be for aliases
   - fix documentation of commands in Lua module template
 - Cleanup of module files and docker requirement bugfix (0.0.27)
 - Docker support (0.0.26)
   - added an enable_tty setting to allow disabling add -t in recipes
   - enforced not adding extra settings to setting.yml
   - settings are validated when an update is attempted
 - Podman support (0.0.25)
   - container_tech is now a settings.yml variable
   - allow for custom test interpreter "test_shell"
   - adding config edit for interactive edit
   - adding show --filter to filter list of modules shown
   - split documentation into user and developer guide
   - added automated spell checking with crate-ci/typos
   - support for container features
   - shpc get -e to show path to an envrironment file
   - added support for env section to render to bound environment file
   - added support for container features like gpu
   - allowing for a `container:tag` convention to be used for commands.
   - list defaults to showing modules/tags one per line, unless --short used
 - adding namespaces to make module install, show, inspect easier (0.0.24)
 - container url does not render in docs (0.0.23)
 - addition of tcl modules, removal of un-needed database (0.0.22)
 - first update of all containers, and bugfix to pull (0.0.21)
 - allowing for a custom container base, container_base (0.0.2)
   - adding more nvidia containers, rapidsai and paraview
 - more packages and shell for a container (0.0.19)
 - adding description to docgen, and descriptions for all containers (0.0.18)
   - descriptions are now required
 - adding docgen command (0.0.17)
 - additional recipes and bug fix for uninstall (0.0.16)
 - updating testing to be with ubuntu (on native) (0.0.15)
 - start to addition of shpc test (0.0.14)
 - Bugfix and better documentation for show/pull/list (0.0.13)
 - Adding support for pull from a GitHub release! (0.0.11)
 - Initial creation of project (0.0.1)

# Singularity Registry HPC (shpc)

[![GitHub actions status](https://github.com/singularityhub/singularity-hpc/workflows/singularity-hpc/badge.svg?branch=main)](https://github.com/singularityhub/singularity-hpc/actions?query=branch%3Amain+workflow%3Asingularity-hpc)
[![DOI](https://zenodo.org/badge/354130612.svg)](https://zenodo.org/badge/latestdoi/354130612)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03311/status.svg)](https://doi.org/10.21105/joss.03311)


![https://raw.githubusercontent.com/singularityhub/singularity-hpc/main/docs/assets/img/shpc.png](https://raw.githubusercontent.com/singularityhub/singularity-hpc/main/docs/assets/img/shpc.png)

Singularity HPC is optimized for managing containers in an HPC environment. Currently, this includes
module technologies:

 - [Lmod](https://lmod.readthedocs.io/en/latest/)
 - [Environment Modules](http://modules.sourceforge.net/)

And container technologies:

 - [Singularity](https://github.com/sylabs/singularity)
 - [Podman](https://podman.io)
 - [Docker](https://docker.io)


You can use shpc if you are:

1. a linux administrator wanting to manage containers as modules for your cluster
2. a cluster user that wants to maintain your own folder of custom modules
3. a cluster user that simply wants to pull Singularity images as GitHub packages.

A module technology is required in all cases.

📖️ Read the [documentation](https://singularity-hpc.readthedocs.io/en/latest/) 📖️
⭐️ Browse the [container module collection](https://singularityhub.github.io/singularity-hpc/) ⭐️
 
## 🎨️ Previous Art 🎨️

There are other tools that you might be interested in!

 - [VA Research Computing](https://www.rc.virginia.edu/userinfo/rivanna/software/containers/) has a similar system, but I couldn't find any code.
 - [Community Collections](https://github.com/community-collections/community-collections)
 - [Spack](https://spack.readthedocs.io/en/latest/module_file_support.html) installs modules for software built from source (not containers).
 

## License

This code is licensed under the MPL 2.0 [LICENSE](LICENSE).
---
title: 'Collaborative Container Modules with Singularity Registry HPC'
tags:
  - containers
  - singularity
  - linux
  - registry
authors:
 - name: Vanessa Sochat
   orcid: 0000-0002-4387-3819
   affiliation: 1
 - name: Alec Scott
   orcid: 0000-0001-6255-1308
   affiliation: 2
affiliations:
 - name: Lawrence Livermore National Lab, Livermore, CA, USA
   index: 1
 - name: University of Arizona Research Computing, Tuscon, AZ, USA
   index: 2
date: 17 April 2021
bibliography: paper.bib
---

# Summary

Portability and reproducibility of complex software stacks is essential for researchers to perform their work. High Performance Computing (HPC) environments add another level of complexity, where possibly conflicting dependencies must co-exist. Although container technologies like Singularity [@Kurtzer2017-xj] make it possible to "bring your own environment," without any form of central strategy to manage containers, researchers who seek reproducibility via using containers are tasked with managing their own container collection, often not taking care to ensure that a particular digest or version is used. The reproducibility of the work is at risk, as they cannot easily install and use containers, nor can they share their software with others.

Singularity Registry HPC (shpc) is the first of its kind to provide an easy means for a researcher to add their research software for sharing and collaboration with other researchers to an existing collection of over 200 popular scientific libraries [@da2017biocontainers; @noauthor_undated-kp, @gorgolewski2017bids; @gamblin2015spack; @autamus]. The software installs containers as environment modules [@McLay2011-wu] that are easy to use and read documentation for, and exposes aliases for commands in the container that the researcher can add to their pipeline without thinking about complex interactions with a container. The simple addition of an entry to the registry maintained by shpc comes down to adding a yaml file, and after doing this, another researcher can easily install the same software, down to the digest, to reproduce the original work.


## Statement of Need

Using environment modules [@McLay2011-wu] on HPC clusters is common.
Although writing the recipes can be complex, it's a fairly common practice for cluster administrators to provide
a set of natively installed recipes for their users [@noauthor_undated-bt], or for researchers to develop and deploy their own software via containers. Even well-known package managers like Spack [@noauthor_undated-ae] and EasyBuild [@noauthor_undated-dj] expose software as modules. However, these package manager approaches don't always ensure reproducibility, or ease of development for the researcher. They typically require relying on some subset of system software, the underlying operating system, or even making changes to the system, which is not under the researcher's control. Although using containers in this context has been discussed previously [@noauthor_undated-rj; @noauthor_undated-rc], the majority of these approaches and tools do not make the process of developing and installing container modules easy. The single researcher must either convince a cluster administrator to install dependencies needed for their software, or build a container and manually move and interact with it on the cluster. All of these small challenges come together to make it harder for a researcher to develop and manage their own software, and subsequently to share their approach to reproduce the work. Using Singularity, Podman, or other container technologies installed via Singularity Registry HPC offers a solution to this challenge. The only requirement is the container technology software, and writing a simple configuration file for the registry. By clearly defining commands, and pinning exact versions of scientific software, researchers on high performance computing
clusters can have more confidence in the reproducibility of their work [@Santana-Perez2015-wo; @Boettiger2014-cz; @Wandell2015-yt].

## Usage

Installing shpc is as easy as cloning the repository and installing in place:

```bash
$ git clone https://github.com/singularityhub/singularity-hpc
$ cd singularity-hpc
$ pip install -e .
```

While the defaults are suitable for most, the researcher can customize the location
of registry metadata files, the module directory to which modules are installed, and the directory to which containers are installed.
The user can then use `shpc show` to see readily available recipes, or browse the [library](https://singularityhub.github.io/singularity-hpc/) for an easily searchable interface. Installation comes down to installing a chosen module, loading it, and using it:

```bash
$ shpc install biocontainers/samtools
$ module load biocontainers/samtools
$ samtools
```

The samtools command would be an alias for a much more complicated command that the researcher
would not need to remember, including global options, command options, the full path
to the container, and the path to the executable inside the container. The container
pulled would also be a specific digest, making it easier to reproduce the researcher's workflow.
Finally, although containers typically only provide one entrypoint, there is no limit
to the number of aliases that can be exposed for easy usage.


## Collaborative Registry for Reproducible Science

Creating a registry entry for a scientific container comes down to writing 
a simple `container.yaml` file with basic metadata and description,
the definition of any and all important entrypoints, and the digests to pull.
As soon as a researcher puts their container in an online registry and adds the
entry, new versions of the container are automatically discovered by shpc,
and can be installed by the researcher when they choose.
The user does not need to look in advance for a version if they want the latest provided
by the registry. Software is easy to search for, and with a simple command, the user can quickly see complete
documentation and commands available:

```bash
$ module spider samtools
```

Every module includes a help section with a description, 
a complete list of commands that map to interactions with the container,
and a list of environment variables also available. 
The module system also provides command line completion, so the user can
use tabs to autocomplete the module names. Along with the expected commands to
use the container primary software (e.g., samtools) commands to shell, 
the software automatically generates alias to inspect, run, shell, or inspect 
container metadata. E.g., here is the shell command:

```bash
$ samtools-shell
```

Another compelling example is using a notebook provided by Jupyter Stacks [@cook2017opinionated]. 
Running notebooks that can be exposed via networking ports tends to be very complicated.
A full interaction might look like the following:

```bash
# pull the latest container, a moving target
$ singularity pull docker://jupyter/minimal-notebook

$ singularity exec --home ${HOME} --bind ${HOME}/.local:/home/joyvan/.local \
   minimal-notebook_latest.sif \
   jupyter notebook --no-browser --port=12345 --ip 0.0.0.0
```

With Singularity Registry HPC, the interaction to run the notebook can be figured
out and written down once, and provided as a community recipe. In this case, it's
exposed as the command "run-notebook":

```bash
$ run-notebook
```

which automatically selects a random port, binds the expected directories, and 
starts the notebook. The registry recipes are collaborative in nature because anyone
can open a pull request with a new recipe, or request a container be added by opening
an issue. Automation also ensures that adding and testing new containers, or working on the
code base is easy. Once a container is added, no further work is needed to update
versions for it. By way of a GitHub bot [@noauthor_undated-eh], both the latest version and newly available tags are 
updated automatically, following any filters that the recipe creator has provided for which tags should be added. Finally, on merge to the main branch, the documentation and library are also automatically updated.

## Conclusion

Singularity Registry HPC is the first local container registry that supports
reproducibility, easy of use, and portability of research software.
You can read more about it on the GitHub repository (https://github.com/singularityhub/singularity-hpc) or
the main documentation site (https://singularity-hpc.readthedocs.io).

# References
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
reported by contacting the project team leader @vsoch. All
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
# Contributor's Agreement

This code is licensed under the MPL 2.0 [LICENSE](LICENSE).

# Contributing

When contributing to Singularity Registry Global Client, it is 
important to properly communicate the gist of the contribution. 
If it is a simple code or editorial fix, simply explaining this 
within the GitHub Pull Request (PR) will suffice. But if this is a larger 
fix or Enhancement, it should be first discussed with the project
leader or developers.

Please note we have a code of conduct, described below. Please follow it in
all your interactions with the project members and users.

## Pull Request Process

1. Bug fix PRs should be sent to both the master and development branches.
   Feature enhancements should only be submitted against the development
   branch.
2. Follow the existing code style precedent. This does not need to be strictly
   defined as there are many thousands of lines of examples. Note the lack
   of tabs anywhere in the project, parentheses and spacing, documentation
   style, source code layout, variable scoping, and follow the project's
   standards.
3. Test your PR locally, and provide the steps necessary to test for the
   reviewers.
4. The project's default copyright and header have been included in any new
   source files.
5. All (major) changes to Singularity Registry must be documented in
   [docs](docs). If your PR changes a core functionality, please 
   include clear description of the changes in your PR so that the docs 
   can be updated, or better, submit another PR to update the docs directly.
6. If necessary, update the README.md.
7. The pull request will be reviewed by others, and the final merge must be
   done by the Singularity project lead, @vsoch (or approved by her).


# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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

### Our Responsibilities

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
reported by contacting the project leader @vsoch. All
complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. The project
team is obligated to maintain confidentiality with regard to the reporter of
an incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers, contributors and users who do not follow or enforce the
Code of Conduct in good faith may face temporary or permanent repercussions 
with their involvement in the project as determined by the project's leader(s).

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
---
name: Bug report
about: Report a bug with Singularity Registry HPC Client

---

**Describe the bug**

**To Reproduce**
Steps to reproduce the behavior:

**Expected behavior**
A clear and concise description of what you expected to happen.

**Version of Singularity and Singularity Registry HPC Client**

Anything else?
---
name: Documentation or Tutorial Request
about: What can we explain better, or what issue did you find?

---


---
name: Feature request
about: Request a feature for Singularity Registry HPC Client

---

**Is your feature request related to a problem? Please describe.**

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.
---
name: Question
about: What's on your mind?

---


---
layout: container
name:  "{{ name }}"
maintainer: "@vsoch"
github: "{{ github_url }}"
updated_at: "{{ creation_date }}"
container_url: "{{ container_url }}"
{% if aliases %}aliases:{% for alias in aliases %}
 - "{{ alias.name }}"
{% endfor %}{% endif %}
versions:{% for version in versions %}
 - "{{ version }}"{% endfor %}
{% if description %}description: "{{ description }}"{% endif %}
---

This module is a singularity container wrapper for {{ name }}.
{% if description %}{{ description }}{% endif %}
After [installing shpc](#install) you will want to install this container module:

```bash
$ shpc install {{ name }}
```

Or a specific version:

```bash
$ shpc install {{ name }}:{{ versions.0 }}
```

And then you can tell lmod about your modules folder:

```bash
$ module use ./modules
```

And load the module, and ask for help, or similar.

```bash
$ module load {{ name }}/{{ versions.0 }}
$ module help {{ name }}/{{ versions.0 }}
```

You can use tab for auto-completion of module names or commands that are provided.

<br>

### Commands

When you install this module, you'll be able to load it to make the following commands accessible.
Examples for both Singularity, Podman, and Docker (container technologies supported) are included.

#### {|module_name|}-run:

```bash
$ singularity run {% if bindpaths %}-B {{ bindpaths }} {% endif %}<container>
$ podman run --rm {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container>
$ docker run --rm {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container>
```

#### {|module_name|}-shell:

```bash
$ singularity shell -s {{ shell }} {% if bindpaths %}-B {{ bindpaths }} {% endif %}<container>
$ podman run --it --rm --entrypoint {{ shell }} {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container>
$ docker run --it --rm --entrypoint {{ shell }} {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container>
```

#### {|module_name|}-exec:

```bash
$ singularity exec -s {{ shell }} {% if bindpaths %}-B {{ bindpaths }} {% endif %}<container> "$@"
$ podman run --it --rm --entrypoint "" {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container> "$@"
$ docker run --it --rm --entrypoint "" {% if bindpaths %}-v {{ bindpaths }} {% endif %} -v ${PWD} -w ${PWD} <container> "$@"
```

#### {|module_name|}-inspect:

Podman and Docker only have one inspect type.

```bash
$ podman inspect <container>
$ docker inspect <container>
```

#### {|module_name|}-inspect-runscript:

```bash
$ singularity inspect -r <container>
```

#### {|module_name|}-inspect-deffile:

```bash
$ singularity inspect -d <container>
```

{% if aliases %}{% for alias in aliases %}
#### {{ alias.name }}
       
```bash
$ singularity exec {% if bindpaths %}-B {{ bindpaths }} {% endif %}{% if alias.options %}{{ alias.options }} {% endif %}<container> {{ alias.command }}
$ podman run --it --rm --entrypoint {{ alias.entrypoint }} {% if bindpaths %}-v {{ bindpaths }} {% endif %} {% if alias.options %}{{ alias.options }} {% endif %} -v ${PWD} -w ${PWD} <container> -c "{{ alias.args }} $@"
$ docker run --it --rm --entrypoint {{ alias.entrypoint }} {% if bindpaths %}-v {{ bindpaths }} {% endif %} {% if alias.options %}{{ alias.options }} {% endif %} -v ${PWD} -w ${PWD} <container> -c "{{ alias.args }} $@"
```

{% endfor %}{% else %}

#### {|module_name|}

```bash
$ singularity run {% if bindpaths %}-B {{ bindpaths }}{% endif %}<container>
$ podman run --rm {% if bindpaths %}-v {{ bindpaths }}{% endif %} -v ${PWD} -w ${PWD} <container>
$ docker run --rm {% if bindpaths %}-v {{ bindpaths }}{% endif %} -v ${PWD} -w ${PWD} <container>
```
{% endif %}

In the above, the `<container>` directive will reference an actual container provided
by the module, for the version you have chosen to load. An environment file in the
module folder will also be bound. Note that although a container
might provide custom commands, every container exposes unique exec, shell, run, and
inspect aliases. For anycommands above, you can export:

 - SINGULARITY_OPTS: to define custom options for singularity (e.g., --debug)
 - SINGULARITY_COMMAND_OPTS: to define custom options for the command (e.g., -b)
 - PODMAN_OPTS: to define custom options for podman or docker
 - PODMAN_COMMAND_OPTS: to define custom options for the command

<br>
  
### Install

You can install shpc locally (for yourself or your user base) as follows:

```bash
$ git clone https://github.com/singularityhub/singularity-hpc
$ cd singularity-hpc
$ pip install -e .
```

Have any questions, or want to request a new module or version? [ask for help!](https://github.com/singularityhub/singularity-hpc/issues)
.. _manual-main:

==========================
Singularity Registry (HPC)
==========================

.. image:: https://img.shields.io/github/stars/singularityhub/singularity-hpc?style=social
    :alt: GitHub stars
    :target: https://github.com/singularityhub/singularity-hpc/stargazers



Singularity Registry HPC (shpc) allows you to install containers as modules.
Currently, this includes:

 - `Lmod <https://lmod.readthedocs.io/en/latest/>`_
 - `Environment Modules <http://modules.sourceforge.net/>`_

And container technologies:

 - `Singularity <https://github.com/sylabs/singularity>`_
 - `Podman <https://podman.io>`_
 - `Docker <https://www.docker.com/>`_

And coming soon:

 - `Shifter <https://github.com/NERSC/shifter>`_
 - `Sarus <https://github.com/eth-cscs/sarus>`_


You can use shpc if you are:

1. a linux administrator wanting to manage containers as modules for your cluster
2. a cluster user that wants to maintain your own folder of custom modules
3. a cluster user that simply wants to pull Singularity images as GitHub packages.

The library contains a collection of module recipes that will install containers,
so you can easily use them or write your own. 
To see the code, head over to the `repository <https://github.com/singularityhub/singularity-hpc/>`_.
To browse modules available as containers, see `the library <https://singularityhub.github.io/singularity-hpc/>`_.


.. _main-getting-started:

-----------------------------------------------
Getting started with Singularity Registry (HPC)
-----------------------------------------------

Singularity Registry HPC (shpc) can be installed from pypi or directly from the repository. See :ref:`getting_started-installation` for
installation, and then the :ref:`getting-started` section for using the client.
You can browse modules available at the `Singularity HPC Library <https://singularityhub.github.io/singularity-hpc/>`_.

.. _main-support:

-------
Support
-------

* For **bugs and feature requests**, please use the `issue tracker <https://github.com/singularityhub/singularity-hpc/issues>`_.
* For **contributions**, visit Caliper on `Github <https://github.com/singularityhub/singularity-hpc>`_.

---------
Resources
---------

`GitHub Repository <https://github.com/singularityhub/singularity-hpc>`_
    The code for shpc on GitHub.

`Singularity HPC Library <https://singularityhub.github.io/singularity-hpc/>`_
    Shows modules available to install as containers.

`Autamus Registry <https://autamus.io>`_
    Provides many of the shpc container modules, built directly from spack.


.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 2

   getting_started/index
   getting_started/user-guide
   getting_started/developer-guide

.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1

    api_reference/shpc
    api_reference/internal/modules
.. _api_reference_shpc:

Singularity Registry HPC
========================

These sections detail the internal functions for shpc.

.. automodule:: shpc
shpc package
============

Submodules
----------

shpc.client module
--------------------------

.. automodule:: shpc.client
    :members:
    :undoc-members:
    :show-inheritance:


shpc.logger module
------------------

.. automodule:: shpc.logger
    :members:
    :undoc-members:
    :show-inheritance:


shpc.main module
----------------

.. automodule:: shpc.main
    :members:
    :undoc-members:
    :show-inheritance:


shpc.main.container module
--------------------------

.. automodule:: shpc.main.container
    :members:
    :undoc-members:
    :show-inheritance:


shpc.main.lmod
--------------

.. automodule:: shpc.main.lmod
    :members:
    :undoc-members:
    :show-inheritance:



Internal API
============

These pages document the entire internal API of SHPC.

.. toctree::
   :maxdepth: 4

   shpc
.. _getting_started-use-cases:

=========
Use Cases
=========

Linux Administrator
===================

If you are a linux administrator, you likely want to clone the repository
directly (or use a release when they are available). Then you can install modules
for your users from the local ``registry`` folder, create your own module files
(and contribute them to the repository if they are useful!) and update the
``module_base`` to be where you install modules.

.. code-block:: console

    # an absolute path
    $ shpc config module_base:/opt/lmod/shpc


If you pull or otherwise update the install of shpc, the module files will update
as well. For example, if you start first by seeing what modules are available
to install:

.. code-block:: console

    $ shpc show


And then install a module to your shpc modules directory:

.. code-block:: console

    $ shpc install tensorflow/tensorflow
    Module tensorflow/tensorflow:2.2.2 was created.


Make sure that lmod knows about the folder

.. code-block:: console

    $ module use /opt/lmod/shpc
     
(And likely if you administer an Lmod install you have your preferred way of
doing this). And then you can use your modules just as you would that are provided on
your cluster.

.. code-block:: console

    $ module load tensorflow/tensorflow/2.2.2


You should then be able to use any of the commands that the tensorflow container
provides, e.g., python and python-shell:

.. code-block:: console

    $ python
    Python 3.6.9 (default, Oct 8 2020, 12:12:24) 
    [GCC 8.4.0] on linux
    Type “help”, “copyright”, “credits” or “license” for more information.
    >>> quit()

    $ tensorflow-tensorflow-shell
    ________                _______________         
    ___ __/__________________________________ ____/__ /________   __
    __ / _ _ \_ __ \_ ___/ __ \_ ___/_ /_  __ /_ __ \_ | /| / /
    _ /  / __/ / / /(__ )/ /_/ / /  _ __/  _ / / /_/ /_ |/ |/ / 
    /_/  \___//_/ /_//____/ \____//_/  /_/   /_/ \____/____/|__/
    You are running this container as user with ID 34633 and group 34633,
    which should map to the ID and group for your user on the Docker host. Great!
    Singularity> quit()


If you want to inspect aliases available or singularity commands to debug:

.. code-block:: console

    $ module spider tensorflow/tensorflow/2.2.2/module
    ----------------------------------------------------------------------------------------------------------------------------
     tensorflow/tensorflow/2.2.2: tensorflow/tensorflow/2.2.2/module
    ----------------------------------------------------------------------------------------------------------------------------
      This module can be loaded directly: module load tensorflow/tensorflow/2.2.2/module
      Help:
       This module is a singularity container wrapper for tensorflow/tensorflow v2.2.2
       Commands include:
        - tensorflow-tensorflow-shell:
               singularity shell -s /bin/bash /usr/WS2/sochat1/singularity-hpc/modules/tensorflow/tensorflow/2.2.2/tensorflow-tensorflow-2.2.2-sha256:e2cde2bb70055511521d995cba58a28561089dfc443895fd5c66e65bbf33bfc0.sif
        - python:
               singularity exec --nv /usr/WS2/sochat1/singularity-hpc/modules/tensorflow/tensorflow/2.2.2/tensorflow-tensorflow-2.2.2-sha256:e2cde2bb70055511521d995cba58a28561089dfc443895fd5c66e65bbf33bfc0.sif /usr/local/bin/python”)



Cluster User
============

If you are a cluster user, you can easily install shpc to your own space
(e.g., in ``$HOME`` or ``$SCRATCH`` where you keep software) and then
use the defaults for the lmod base (the modules folder that is created alongside
the install) and the registry. You can also pull the repository to get updated
registry entries.  If you haven't yet, clone the repository:


.. code-block:: console

    $ git clone git@github.com:singularityhub/singularity-hpc.git
    $ cd singularity-hpc
    
You can then see modules available for install:


.. code-block:: console

    $ shpc show


And install a module to your local modules folder.


.. code-block:: console

    $ shpc install python
    Module python/3.9.2-slim was created.


Finally, you can add the module folder to those that lmod knows about:

.. code-block:: console

    $ module use $HOME/singularity-hpc/modules
     
And then you can use your modules just as you would that are provided on
your cluster.

.. code-block:: console

    $ module load python/3.9.2-slim


An error will typically be printed if there is a conflict with another module name, and it's
up to you to unload the conflicting module(s) and try again. For this module, 
since we didn't use a prefix the container python will be exposed
as "python" - an easier one to see is "python-shell" - each container exposes
a shell command so you can quickly get an interactive shell. 
Every installed entry will have it's named suffixed with "shell" if you quickly want
an interactive session. For example:

.. code-block:: console

    $ python-shell
    Singularity>


And of course running just "python" gives you the Python interpreter. If you
don't know the command that you need, or want to see help for the module you
loaded, just do:

.. code-block:: console

    $ module spider python/3.9.2-slim/module
    ----------------------------------------------------------------------------------------------------------------------------
    python/3.9.2-slim: python/3.9.2-slim/module
    ----------------------------------------------------------------------------------------------------------------------------
      This module can be loaded directly: module load python/3.9.2-slim/module
      Help:
       This module is a singularity container wrapper for python v3.9.2-slim
       Commands include:
        - python-shell:
           singularity shell -s /bin/bash /usr/WS2/sochat1/singularity-hpc/modules/python/3.9.2-slim/python-3.9.2-slim-    sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef.sif
        - python:
           singularity exec /usr/WS2/sochat1/singularity-hpc/modules/python/3.9.2-slim/python-3.9.2-slim-sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef.sif /usr/local/bin/python”)
    
    
The above not only shows you the description, but also the commands if you
need to debug. If you want to see metadata about the container (e.g., labels,
singularity recipe) then you can do:

.. code-block:: console

    $ module whatis python/3.9.2-slim
    python/3.9.2-slim/module             : Name    : python/3.9.2-slim
    python/3.9.2-slim/module             : Version   : module
    python/3.9.2-slim/module             : URL     : https://hub.docker.com/_/python
    python/3.9.2-slim/module             : Singularity Recipe  : bootstrap: docker 
    from: python@sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef
    python/3.9.2-slim/module             : org.label-schema.build-arch  : amd64
    python/3.9.2-slim/module             : org.label-schema.build-date  : Sunday_4_April_2021_19:56:56_PDT
    python/3.9.2-slim/module             : org.label-schema.schema-version  : 1.0
    python/3.9.2-slim/module             : org.label-schema.usage.singularity.deffile.bootstrap  : docker
    python/3.9.2-slim/module             : org.label-schema.usage.singularity.deffile.from  : python@sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef
    python/3.9.2-slim/module             : org.label-schema.usage.singularity.version  : 3.7.1-1.el7


Adding Options
--------------

By default, some of the commands will come with singularity options. For example,
a container intended for gpu is always going to give you the ``--nv`` flag. However,
it could be the case that you want to define custom options at the time of use.
In this case, you can export the following custom environment variables to add them:

**SINGULARITY_OPTS**: will provide additional options to the base Singularity command, such as ``--debug``
**SINGULARITY_COMMAND_OPTS**: will provide additional options to the command (e.g., exec), such as ``--cleanenv``.


Custom Images that are Added
============================

If you add a custom image, the interaction is similar, whether you are a cluster
user or administrator. First, let's say we pull a container:

.. code-block:: console

    $ singularity pull docker://vanessa/salad
    
And we add it to our unique namespace in the modules folder:

.. code-block:: console

    $ shpc add salad_latest.sif vanessa/salad:latest
    
    
We can again load the custom module:

.. code-block:: console

    $ module load vanessa/salad/latest


Since we didn't define any aliases via a registry entry, the defaults provided
are to run the container (the squashed unique resource identifier, ``vanessa-salad-latest``
or the same shell, ``vanessa-salad-latest-shell``. Of course you can check this if you don't know:

.. code-block:: console

    $ module spider vanessa/salad/latest/module 
    --------------------------------------------------------------------------------------------------------------------------------------------------------
     vanessa/salad/latest: vanessa/salad/latest/module
    --------------------------------------------------------------------------------------------------------------------------------------------------------
      This module can be loaded directly: module load vanessa/salad/latest/module
      Help:
       This module is a singularity container wrapper for vanessa-salad-latest vNone
       Commands include:
        - vanessa-salad-latest-shell:
           singularity shell -s /bin/bash /usr/WS2/sochat1/singularity-hpc/modules/vanessa/salad/latest/vanessa-salad-latest-sha256:71d1f3e42c1ceee9c02295577c9c6dfba4f011d9b8bce82ebdbb6c187b784b35.sif
        - vanessa-salad-latest: singularity run /usr/WS2/sochat1/singularity-hpc/modules/vanessa/salad/latest/vanessa-salad-latest-sha256:71d1f3e42c1ceee9c02295577c9c6dfba4f011d9b8bce82ebdbb6c187b784b35.sif


And then use them! For example, the command without ``-shell`` just runs the container:

.. code-block:: console

    $ vanessa-salad-latest
     You think you have problems? I’m a fork.  
                /\
               //\\
               // \\
             ^  \\ //  ^
            / \  ) (  / \ 
            ) (  ) (  ) (
            \ \_/ /\ \_/ /
             \__ _)(_ __/
              \ \ / /
               ) \/ (
               | /\ |
               | )( |
               | )( |
               | \/ |
               )____(
              /   \
              \______/ 


And the command with shell does exactly that.

.. code-block:: console

    $ vanessa-salad-latest-shell
    Singularity> exit

If you need more robust commands than that, it's recommended to define your
own registry entry. If you think it might be useful to others, please contribute it
to the repository!


Pull Singularity Images
=======================

Singularity Registry HPC tries to support researchers that cannot afford to
pay for a special Singularity registry, and perhaps don't want to pull
from a Docker URI. For this purpose, you can use the `Singularity Deploy <https://github.com/singularityhub/singularity-deploy>`_
template to create containers as releases associated with the same GitHub
repository, and then pull them down directly with the shpc client with
the ``gh://`` unique resource identifier as follows:

.. code-block:: console

    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:latest
    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:salad
    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:pokemon


In the example above, our repository is called ``singularityhub/singularity-deploy``,
and in the root we have three recipes:

 - Singularity (builds to latest)
 - Singularity.salad
 - Singularity.pokemon

And in the ``VERSION`` file in the root, we have ``0.0.1`` which corresponds with
the GitHub release. This will pull to a container.  For example:

.. code-block:: console

    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:latest
    singularity pull --name /home/vanessa/Desktop/Code/singularity-hpc/singularityhub-singularity-deploy.latest.sif https://github.com/singularityhub/singularity-deploy/releases/download/0.0.1/singularityhub-singularity-deploy.latest.sif
    /home/vanessa/Desktop/Code/singularity-hpc/singularityhub-singularity-deploy.latest.sif

And then you are ready to go!

.. code-block:: console

    $ singularity shell singularityhub-singularity-deploy.latest.sif 
    Singularity> 


See the `Singularity Deploy <https://github.com/singularityhub/singularity-deploy>`_ repository
for complete details for how to set up your container! Note that this uri (``gh://``)
can also be used in a registry entry.
.. _getting_started-user-guide:

==========
User Guide
==========

Singularity Registry HPC (shpc) will allow you to install Singularity containers as
modules. This means that you can install them as a cluster admin, or as a cluster user.
This getting started guide will walk you through setting up a local registry,
either for yourself or your user base. If you haven't read :ref:`getting_started-installation`
you should do that first.

Why shpc?
=========

Singularity Registry HPC is created to be modular, meaning that we support a distinct
set of container technologies and module systems. The name of the library "Singularity
Registry HPC" does not refer specifically to the container technology "Singularity,"
but more generally implies the same spirit -- a single entity that is "one library to rule them all!"


What is a registry?
===================

A registry consists of a database of local containers configuration files, ``container.yaml``
files organized in the root of the shpc install in one of the ``registry`` folders. The namespace
is organized by Docker unique resources identifiers. When you install an identifier
as we saw above, the container binaries and customized module files are added to 
the ``module_dir`` defined in your settings, which defaults to ``modules`` in the
root of the install. You should see the :ref:`getting_started-developer-guide`
for more information about contributing containers to this registry.


Really Quick Start
==================

Once you have shpc installed, make sure you tell shpc what your module software is
(note that you only need to run this command if you aren't using Lmod, which is the
default).

.. code-block:: console

    $ shpc config set module_sys:tcl
    $ shpc config set module_sys:lmod  # default


You can then easily install, load, and use modules:

.. code-block:: console

    $ shpc install biocontainers/samtools
    $ module load biocontainers/samtools
    $ samtools


The above assumes that you've installed the software, and have already
added the modules folder to be seen by your module software. If your module
software doesn't see the module, remember that you need to have done:

.. code-block:: console

    $ module use $(pwd)/modules


We walk through these steps in more detail in the next section.


Quick Start
===========

After  :ref:`getting_started-installation`, and let's say shpc is installed 
at ``~/singularity-hpc`` you can edit your settings in ``settings.yaml``.
Importantly, make sure your shpc install is configured to use the right module
software, which is typicall lmod or tcl. Here is how to change from the default 
"lmod" to "tcl" and then back:

.. code-block:: console

    $ shpc config set module_sys:tcl
    $ shpc config set module_sys:lmod # this is the default, which we change back to!


Once you have the correct module software indicated, try installing a container:

.. code-block:: console

    $ shpc install python
    
Make sure that the local ./modules folder can be seen by your module software
(you can run this in a bash profile or manually, and note that if you want to 
use Environment Modules, you need to add ``--module-sys tcl``).

.. code-block:: console

    $ module use ~/singularity-hpc/modules


And then load the module!

.. code-block:: console

    $ module load python/3.9.2-slim

If the module executable has a conflict with something already loaded, it
will tell you, and it's up to you to unload the conflicting modules before you
try loading again. If you want to quickly see commands that are supported, use module
help:

.. code-block:: console

    $ module help python/3.9.2-slim

If you want to add the modules folder to your modules path more permanently,
you can add it to ``MODULEPATH`` in your bash profile.

.. code-block:: console

    export MODULEPATH=$HOME/singularity-hpc/modules:$MODULEPATH


For more detailed tutorials, you should continue reading,
and see :ref:`getting_started-use-cases`. Also see the :ref:`getting_started-commands-config` for how to update configuration values with ``shpc config``.


Setup
=====

Setup includes, after installation, editing any configuration values to
customize your install. The configuration file will default to ``shpc/settings.yml``
in the installed module, however you can create your own user settings file to
take preference over this one as follows:

.. code-block:: console

    $ shpc config userinit


When you create a user settings file (or provide a custom settings file one off to
the client) the shpc default settings will be read first, and then updated by your file.
We do this so that if the default file updates and your user settings is missing a variable,
we still use the default. The defaults in either file are likely suitable for most. For any configuration value 
that you might set, the following variables are available to you:

 - ``$install_dir``: the shpc folder
 - ``$root_dir``: the parent directory of shpc (where this README.md is located)


Additionally, the variables ``module_base``, ``container_base``, and ``registry``
can be set with environment variables that will be expanded at runtime. You cannot
use the protected set of substitution variables (``$install_dir`` and ``$install_root``)
as environment variables, as they will be subbed in by shpc before environment
variable replacement. A summary table of variables is included below, and then further discussed in detail.


.. list-table:: Title
   :widths: 25 65 10
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - module_sys
     - Set a default module system. Currently lmod and tcl are supported
     - lmod
   * - registry
     - A list of full paths to one or more registry folders (with subfolders with container.yaml recipes)
     - [$root_dir/registry]
   * - module_base
     - The install directory for modules
     - $root_dir/modules
   * - container_base
     - Where to install containers. If not defined, they are installed alongside modules.
     - null
   * - container_tech
     - The container technology to use (singularity or podman)
     - singularity
   * - updated_at
     - a timestamp to keep track of when you last saved
     - never
   * - default_version
     - A boolean to indicate generating a .version file (LMOD or lua modules only)
     - true
   * - singularity_module
     - if defined, add to module script to load this Singularity module first
     - null
   * - module_name
     - Format string for module commands exec,shell,run (not aliases) can include ``{{ registry }}``, ``{{ repository }}``, ``{{ tool }}`` and ``{{ version }}``
     - ``'{{ tool }}'``
   * - bindpaths
     - string with comma separated list of paths to binds. If set, expored to SINGULARITY_BINDPATH
     - null
   * - singularity_shell
     - exported to SINGULARITY_SHELL
     - /bin/sh
   * - podman_shell
     - The shell used for podman
     - /bin/sh
   * - docker_shell
     - The shell used for docker
     - /bin/sh
   * - test_shell
     - The shell used for the test.sh file
     - /bin/bash
   * - namespace
     - Set a default module namespace that you want to install from.
     - null
   * - environment_file
     - The name of the environment file to generate and bind to the container.
     - 99-shpc.sh
   * - enable_tty
     - For container technologies that require -t for tty, enable (add) or disable (do not add)
     - true
   * - config_editor
     - The editor to use for your config editing
     - vim
   * - features
     - A key, value paired set of features to add to the container (see table below)
     - All features default to null


These settings will be discussed in more detail in the following sections.

Features
--------

Features are key value pairs that you can set to a determined set of values
to influence how your module files are written. For example, if you set the
gpu feature to "nvidia" in your settings file:

.. code-block:: yaml

    container_features:
      gpu: "nvidia"


and a container.yaml recipe has a gpu:true container feature to say "this container
supports gpu":

.. code-block:: yaml

    features:
      gpu: true
     
Given that you are installing a module for a Singularity container, the ``--nv``
option will be added. Currently, the following features are supported:


.. list-table:: Title
   :widths: 10 40 25 25
   :header-rows: 1

   * - Name
     - Description
     - Default
     - Options
   * - gpu
     - If the container technology supports it, add flags to indicate using gpu.
     - null
     - nvidia, amd, null
   * - x11
     - Bind mount ~/.Xauthority or a custom path
     - null
     - true (uses default path ~/.Xauthority), false/null (do not enable) or a custom path to an x11 file
   * - home
     - Specify and bind mount a custom home path
     - null
     - custom path for the home, or false/null


Modules Folder
--------------

The first thing you want to do is configure your module location, if you want it different
from the default. The path can be absolute or relative to ``$install_dir`` (the shpc
directory) or ``$root_dir`` (one above that) in your
configuration file at ``shpc/settings.yml``. If you are happy
with module files being stored in a ``modules`` folder in the present working
directory, you don't need to do any configuration. Otherwise, you can customize
your install:

.. code-block:: console

    # an absolute path
    $ shpc config set module_base:/opt/lmod/modules

    # or a path relative to a variable location remember to escape the "$"
    $ shpc config set module_base:\$root_dir/modules


This directory will be the base where lua files are added, and container are stored.
For example, if you were to add a container with unique resource identifier `python/3.8`
you would see:

.. code-block:: console

    $install_dir/modules/
    └── python
        └── 3.9.2
            ├── module.lua
            └── python-3.9.2.sif

Although your module path might have multiple locations, Singularity Registry HPC 
assumes this one location to install container modules to in order to ensure
a unique namespace. 


Container Images Folder
-----------------------

If you don't want your container images (sif files) to live alongside your
module files, then you should define the ``container_base`` to be something
non-null (a path that exists). For example:

.. code-block:: console

    $ mkdir -p /tmp/containers
    $ shpc config set container_base:/tmp/containers


The same hierarchy will be preserved as to not put all containers in the same
directory.


Registry
--------

The registry parameter is a list of one or more registry locations (filesystem
directories) where shpc will search for ``container.yaml`` files. The default
registry shipped with shpc is the folder in the root of the repository, but 
you can add or remove entries via the config variable ``registry``


.. code-block:: console

    # change to your own registry of container yaml configs
    $ shpc config add registry:/opt/lmod/registry


# Note that "add" is used for lists of things (e.g., the registry config variable is a list)
and "set" is used to set a key value pair.


Module Names
------------

The setting ``module_name`` is a format string in `Jinja2 <https://jinja.palletsprojects.com/en/3.0.x/>`_ 
that is used to generate your module command names. For each module, in addition
to aliases that are custom to the module, a set of commands for run, inspect, exec,
and shell are generated. These commands will use the ``module_name`` format string
to determine their names. For example, for a python container with the default ``module_name``
of "{{ tool }}" we will derive the following aliases for a Singularity module:

.. code-block:: console

    python-shell
    python-run
    python-exec
    python-inspect-deffile
    python-inspect-runscript

A container identifier is parsed as follows:

.. code-block:: console

    # quay.io   /biocontainers/samtools:latest
    # <registry>/ <repository>/  <tool>:<version>


So by default, we use tool because it's likely closest to the command that is wanted.
But let's say you had two versions of samtools - the namespaces would conflict! You
would want to change your format string to ``{{ repository }}-{{ tool }}`` to be
perhaps "biocontainers-samtools-exec" and "another-samtools-exec." 
If you change the format string to ``{{ tool }}-{{ version }}`` you would see:

.. code-block:: console

    python-3.9.5-alpine-shell
    python-3.9.5-alpine-run
    python-3.9.5-alpine-exec
    python-3.9.5-alpine-deffile
    python-3.9.5-alpine-runscript


And of course you are free to add any string that you wish, e.g., ``plab-{{ tool }}``

.. code-block:: console

    plab-python-shell

These prefixes are currently only provided to the automatically generated
commands. Aliases that are custom to the container are not modified.


Module Software
---------------

The default module software is currently Lmod, and there is also support for environment
modules that only use tcl (tcl). If you
are interested in adding another module type, please `open an issue <https://github.com/singularityhub/singularity-hpc>`_ and
provide description and links to what you have in mind. You can either specify the
module software on the command line:


.. code-block:: console

    $ shpc install --module-sys tcl python


or you can set the global variable to what you want to use (it defaults to lmod):

.. code-block:: console

    $ shpc config set module_sys:tcl


The command line argument, if provided, always over-rides the default.


Container Technology
--------------------

The default container technology to pull and then provide to users is Singularity,
and we have also recently added Podman and Docker, and will add support for Shifter and Sarus soon.
Akin to module software, you can specify the container technology to use on a global
setting, or via a one-off command:


.. code-block:: console

    $ shpc install --container-tech podman python


or for a global setting:

.. code-block:: console

    $ shpc config set container_tech:podman


If you would like support for a different container technology that has not been
mentioned, please also `open an issue <https://github.com/singularityhub/singularity-hpc>`_ and
provide description and links to what you have in mind.

.. _getting_started-commands:


Commands
========

The following commands are available! For any command, the default module system
is lmod, and you can change this to tcl by way of adding the ``--module-sys`` argument
after your command of interest.

.. code-block:: console

    $ shpc <command> --module-sys tcl <args>


.. _getting_started-commands-config:


Config
------

If you want to edit a configuration value, you can either edit the ``shpc/settings.yml``
file directly, or you can use ``shpc config``, which will accept:

 - set to set a parameter and value
 - get to get a parameter by name
 - add to add a value to a parameter that is a list (e.g., registry)
 - remove to remove a value from a parameter that is a list

The following example shows changing the default module_base path from the install directory modules folder.

.. code-block:: console

    # an absolute path
    $ shpc config set module_base:/opt/lmod/modules

    # or a path relative to the install directory, remember to escape the "$"
    $ shpc config set module_base:\$install_dir/modules


And then to get values:

.. code-block:: console

    $ shpc config get module_base


And to add and remove a value to a list:

.. code-block:: console

    $ shpc config add registry:/tmp/registry
    $ shpc config remove registry:/tmp/registry


You can also open the config in the editor defined in settings at ``config_editor``

.. code-block:: console

    $ shpc config edit
    

which defaults to vim.

Show and Install
----------------

The most basic thing you might want to do is install an already existing
recipe in the registry. You might first want to show the known registry entries
first. To show all entries, you can run:

.. code-block:: console

    $ shpc show
    tensorflow/tensorflow
    python
    singularityhub/singularity-deploy

The default will not show versions available. To flatten out this list and include versions for each, you can do:

.. code-block:: console

    $ shpc show --versions
    tensorflow/tensorflow:2.2.2
    python:3.9.2-slim
    python:3.9.2-alpine
    singularityhub/singularity-deploy:salad


To filter down the result set, use ``--filter``:


.. code-block:: console

    $ shpc show --filter bio
    biocontainers/bcftools
    biocontainers/vcftools
    biocontainers/bedtools
    biocontainers/tpp


To get details about a package, you would then add it's name to show:

.. code-block:: console

    $ shpc show python


And then you can install a version that you like (or don't specify to default to
the latest, which in this case is 3.9.2-slim). You will see the container pulled, 
and then a message to indicate that the module was created. 


.. code-block:: console
    
    $ shpc install python
    ...
    Module python/3.9.2 is created.


.. code-block:: console

    $ tree modules/
    modules/
    └── python
        └── 3.9.2
            ├── module.lua
            └── python-3.9.2.sif

    2 directories, 2 files
    

You can also install a specific tag (as shown in list).
    
.. code-block:: console

    $ shpc install python:3.9.2-alpine
    

Note that Lmod is the default for the module system, and Singularity for
the container technology.
If you don't have any module software on your system, you can now test interacting
with the module via the :ref:`getting_started-development` instructions.


Install Private Images
----------------------

What about private containers on Docker Hub? If you have a private image, you can
simply use `Singularity remote login <https://github.com/sylabs/singularity-userdocs/blob/master/singularity_and_docker.rst#singularity-cli-remote-command>`_ before attempting the install and everything should work.

Namespace
---------

Let's say that you are exclusively using continers in the namespace ghcr.io/autamus.

.. code-block:: console

    registry/ghcr.io/
    └── autamus
        ├── abi-dumper
        ├── abyss
        ├── accumulo
        ├── addrwatch
        ...
        ├── xrootd
        ├── xz
        └── zlib


It can become arduous to type the entire namespace every time! For this purpose,
you can set a namespace:

.. code-block:: console

    $ shpc namespace use ghcr.io/autamus

And then instead of asking to install clingo as follows:

.. code-block:: console

    $ shpc install ghcr.io/autamus/clingo
    

You can simply ask for:


.. code-block:: console

    $ shpc install clingo
    
    
And when you are done, unset the namespace.


.. code-block:: console

    $ shpc namespace unset


Note that you can also set the namespace as any other setting:

.. code-block:: console

    $ shpc config set namespace:ghcr.io/autamus

Namespaces currently work with:

 - install
 - uninstall
 - show
 - add
 - check

List
----

Once a module is installed, you can use list to show installed modules (and versions).
The default list will flatten out module names and tags into a single list
to make it easy to copy paste:

.. code-block:: console

    $ shpc list
        biocontainers/samtools:v1.9-4-deb_cv1
                        python:3.9.2-alpine
                        python:3.9.5-alpine
                        python:3.9.2-slim
                      dinosaur:fork
                 vanessa/salad:latest
                         salad:latest
      ghcr.io/autamus/prodigal:latest
      ghcr.io/autamus/samtools:latest
        ghcr.io/autamus/clingo:5.5.0


However, if you want a shorter version that shows multiple tags alongside
each unique module name, just add ``--short``:

.. code-block:: console

    $ shpc list --short

        biocontainers/samtools: v1.9-4-deb_cv1
                        python: 3.9.5-alpine, 3.9.2-alpine, 3.9.2-slim
                      dinosaur: fork
                 vanessa/salad: latest
                         salad: latest
      ghcr.io/autamus/prodigal: latest
      ghcr.io/autamus/samtools: latest
        ghcr.io/autamus/clingo: 5.5.0


Inspect
-------

Once you install a module, you might want to inspect the associated container! You
can do that as follows:

.. code-block:: console

    $ shpc inspect python:3.9.2-slim
    👉️ ENVIRONMENT 👈️
    /.singularity.d/env/10-docker2singularity.sh : #!/bin/sh
    export PATH="/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    export LANG="${LANG:-"C.UTF-8"}"
    export GPG_KEY="${GPG_KEY:-"E3FF2839C048B25C084DEBE9B26995E310250568"}"
    export PYTHON_VERSION="${PYTHON_VERSION:-"3.9.2"}"
    export PYTHON_PIP_VERSION="${PYTHON_PIP_VERSION:-"21.0.1"}"
    export PYTHON_GET_PIP_URL="${PYTHON_GET_PIP_URL:-"https://github.com/pypa/get-pip/raw/b60e2320d9e8d02348525bd74e871e466afdf77c/get-pip.py"}"
    export PYTHON_GET_PIP_SHA256="${PYTHON_GET_PIP_SHA256:-"c3b81e5d06371e135fb3156dc7d8fd6270735088428c4a9a5ec1f342e2024565"}"
    /.singularity.d/env/90-environment.sh : #!/bin/sh
    # Custom environment shell code should follow

    👉️ LABELS 👈️
    org.label-schema.build-arch : amd64
    org.label-schema.build-date : Sunday_4_April_2021_20:51:45_MDT
    org.label-schema.schema-version : 1.0
    org.label-schema.usage.singularity.deffile.bootstrap : docker
    org.label-schema.usage.singularity.deffile.from : python@sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef
    org.label-schema.usage.singularity.version : 3.6.0-rc.4+501-g42a030f8f

    👉️ DEFFILE 👈️
    bootstrap: docker
    from: python@sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef


We currently don't show the runscript, as they can be very large. However, if you want
to see it:

    $ shpc inspect --runscript python:3.9.2-slim


Or to get the entire metadata entry dumped as json to the terminal:

.. code-block:: console

    $ shpc inspect --json python:3.9.2-slim


.. _getting_started-commands-test:



Test
----

Singularity HPC makes it easy to test the full flow of installing and interacting
with modules. This functionality requires a module system (e.g., Lmod) to be installed,
and the assumption is that the test is being run in a shell environment where any
supporting modules (e.g., loading Singularity or Podman) would be found if needed.
This is done by way of extending the exported ``$MODULEPATH``. To run a test, you
can do:

.. code-block:: console

    shpc test python


If you don't have it, you can run tests in the provided docker container. 

.. code-block:: console

    docker build -t singularity-hpc .
    docker run --rm -it singularity-hpc shpc test python


Note that the ``Dockerfile.tcl`` builds an equivalent container with tcl modules.

.. code-block:: console

    $ docker build -f Dockerfile.tcl -t singularity-hpc .


If you want to stage a module install (e.g., install to a temporary directory and not remove it) do:


.. code-block:: console

    shpc test --stage python


To do this with Docker you would do:

.. code-block:: console

    $ docker run --rm -it singularity-hpc bash
    [root@1dfd9fe90443 code]# shpc test --stage python
    ...
    /tmp/shpc-test.fr1ehcrg


And then the last line printed is the directory where the stage exists,
which is normally cleaned up. You can also choose to skip testing the module
(e.g., lmod):


.. code-block:: console

    shpc test --skip-module python


Along with testing the container itself (the commands are defined in the ``tests``
section of a ``container.yaml``.


.. code-block:: console

    shpc test --skip-module --commands python


Uninstall
---------

To uninstall a module, since we are targeting a module folder, instead of
providing a container unique resource identifier like `python:3.9.2-alpine`,
we provide the module path relative to your module directory. E.g.,

.. code-block:: console

    $ shpc uninstall python:3.9.2-alpine


You can also uninstall an entire family  of modules:

.. code-block:: console

    $ shpc uninstall python

The uninstall will go up to the top level module folder but not remove it
in the case that you've added it to your ``MODULEPATH``.

Pull
----

Singularity Registry HPC tries to support researchers that cannot afford to
pay for a special Singularity registry, and perhaps don't want to pull
from a Docker URI. For this purpose, you can use the `Singularity Deploy <https://github.com/singularityhub/singularity-deploy>`_
template to create containers as releases associated with the same GitHub
repository, and then pull them down directly with the shpc client with
the ``gh://`` unique resource identifier as follows:

.. code-block:: console

    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:latest
    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:salad
    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:pokemon


In the example above, our repository is called ``singularityhub/singularity-deploy``,
and in the root we have three recipes:

 - Singularity (builds to latest)
 - Singularity.salad
 - Singularity.pokemon

And in the ``VERSION`` file in the root, we have ``0.0.1`` which corresponds with
the GitHub release. This will pull to a container.  For example:

.. code-block:: console

    $ shpc pull gh://singularityhub/singularity-deploy/0.0.1:latest
    singularity pull --name /home/vanessa/Desktop/Code/singularity-hpc/singularityhub-singularity-deploy.latest.sif https://github.com/singularityhub/singularity-deploy/releases/download/0.0.1/singularityhub-singularity-deploy.latest.sif
    /home/vanessa/Desktop/Code/singularity-hpc/singularityhub-singularity-deploy.latest.sif

And then you are ready to go!

.. code-block:: console

    $ singularity shell singularityhub-singularity-deploy.latest.sif 
    Singularity> 


See the `Singularity Deploy <https://github.com/singularityhub/singularity-deploy>`_ repository
for complete details for how to set up your container! Note that this uri (``gh://``)
can also be used in a registry entry.


Shell
-----

If you want a quick way to shell into an installed module's container
(perhaps to look around or debug without the module software being available) you can use
``shell``. For example:

.. code-block:: console

    shpc shell vanessa/salad:latest
    Singularity> /code/salad fork

     My life purpose: I cut butter.  
    
                       ________  .====
                      [________>< :===
                                 '==== 



If you want to interact with the shpc Python client directly, you can
do shell without a module identifier. This will give you a python terminal,
which defaults to ipython, and then python and
bypython (per what is available on your system). To start a shell:

.. code-block:: console

    $ shpc shell


or with a specific interpreter:

.. code-block:: console

    $ shpc shell -i python


And then you can interact with the client, which will be loaded.

.. code-block:: python

    client
    [shpc-client]

    client.list()
    python

    client.install('python')



Show
----

As shown above, show is a general command to show the metadata file for a registry entry:

.. code-block:: console

    $ shpc show python
    docker: python
    latest:
      3.9.2-slim: sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef
    tags:
      3.9.2-slim: sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef
      3.9.2-alpine: sha256:23e717dcd01e31caa4a8c6a6f2d5a222210f63085d87a903e024dd92cb9312fd
    filter:
    - 3.9.*
    maintainer: '@vsoch'
    url: https://hub.docker.com/_/python
    aliases:
      python: /usr/local/bin/python

Or without any arguments, it will show a list of all registry entries available:

.. code-block:: console

    $ shpc show
    python


Check
-----

How do you know if there is a newer version of a package to install? In
the future, if you pull updates from the main repository, we will have a bot
running that updates container versions (digests) as well as tags. Here
is how to check if a module (the tag) is up to date.

.. code-block:: console

    $ shpc check tensorflow/tensorflow
    ⭐️ latest tag 2.2.2 is up to date. ⭐️


And if you want to check a specific digest for tag (e.g., if you use "latest" it
is subject to change!)

.. code-block:: console

    $ shpc check tensorflow/tensorflow:2.2.2
    ⭐️ tag 2.2.2 is up to date. ⭐️

As a trick, you can loop through registry entries with ``shpc show``. The return
value will be 0 is there are no updates, and 1 otherwise. This is a trick
we use to check for new recipes to test.

.. code-block:: console



Add
---

It might be the case that you have a container locally, and you want to
make it available as a module (without pulling it from a registry). Although
this is discouraged because it means you will need to manually maintain
versions, shpc does support the "add" command to do this. You can simply provide
the container path and the unique resource identifier:

.. code-block:: console

    $ shpc add salad_latest.sif vanessa/salad:latest

If the unique resource identifier corresponds with a registry entry, you
will not be allowed to create it, as this would create a namespace conflict.
Since we don't have a configuration file to define custom aliases, the container
will just be exposed as it's command to run it.

Get
---

If you want to quickly get the path to a container binary, you can use get.

.. code-block:: console

    $ shpc get vanessa/salad:latest
    /home/vanessa/Desktop/Code/singularity-hpc/modules/vanessa/salad/latest/vanessa-salad-latest-sha256:8794086402ff9ff9f16c6facb93213bf0b01f1e61adf26fa394b78587be5e5a8.sif

    $ shpc get tensorflow/tensorflow:2.2.2
    /home/vanessa/Desktop/Code/singularity-hpc/modules/tensorflow/tensorflow/2.2.2/tensorflow-tensorflow-2.2.2-sha256:e2cde2bb70055511521d995cba58a28561089dfc443895fd5c66e65bbf33bfc0.sif

If you select a higher level module directory or there is no sif, you'll see:

.. code-block:: console

    $ shpc get tensorflow/tensorflow
    tensorflow/tensorflow is not a module tag folder, or does not have a sif binary.


You can add ``-e`` to get the environment file:


.. code-block:: console

    $ shpc get -e tensorflow/tensorflow


We could update this command to allow for listing all sif files within a top level
module folder (for different versions). Please open an issue if this would be useful for
you.
.. _getting_started-installation:

============
Installation
============

Singularity Registry HPC (shpc) can be installed from pypi, or from source. 
In all cases, a module technology is required such as `lmod (install intstructions) <https://lmod.readthedocs.io/en/latest/030_installing.html>`_
or `environment modules (install instructions) <https://modules.readthedocs.io/en/latest/INSTALL.html>`_.
Having module software installed means that the ``module`` command should be on your path.
Once you are ready to install shpc along your module software, it's recommended that you create a virtual environment, if you have not already
done so.


Virtual Environment
===================

The recommended approach is to install from the repository directly, whether
you use pip or another setup approach, and to install a `known release <https://github.com/singularityhub/singularity-hpc/releases/>`_. Here is how to clone the repository and do a local install.

.. code:: console

    # Install release 0.0.24
    $ git clone -b 0.0.24 git@github.com:singularityhub/singularity-hpc
    $ cd singularity-hpc
    $ pip install -e .[all]

or (an example with python setuptools and installing from the main branch)

.. code:: console

    $ git clone git@github.com:singularityhub/singularity-hpc
    $ cd singularity-hpc
    $ python setup.py develop


if you install to a system python, meaning either of these commands:

.. code:: console

    $ python setup.py install
    $ pip install .

You will need to put the registry files elsewhere (update the ``registry`` config argument to the path), as they will not be installed
alongside the package. The same is the case for modules - if you install to system
python it's recommended to define ``module_base`` as something else, unless you
can write to your install location. Installing locally ensures that you
can easily store your module files along with the install (the default until you
change it). Installation of singularity-hpc adds an executable, `shpc` to your path.

.. code:: console

    $ which shpc
    /opt/conda/bin/shpc


This executable should be accessible by an administrator, or anyone that you want
to be able to manage containers. Your user base will be interacting with your
containers via Lmod, so they do not need access to `shpc`. 
If you are a user creating your own folder of modules, you can add them
to your module path.

Once it's installed, you should be able to inspect the client!


.. code-block:: console

    $ shpc --help


You'll next want to configure and create your registry, discussed next in :ref:`getting-started`.
Generally, remember that your modules will be installed in
the ``modules`` folder, and container recipes are nested in ``registry``. If you don't
want your container images (sif files) installed alongside your module recipes,
then you can define ``container_base`` to be somewhere else. You
can change these easily with ``shpc config``, as they are defined via these
variables in the config:
 

.. code-block:: console
 

    $ shpc config set registry:/<DIR>
    $ shpc config set module_base:/<DIR> 
    $ shpc config set container_base:/<DIR> 


Also importantly, if you are using environment modules (Tcl) and not LMOD, you need
to tell shpc about this (as it defaults to LMOD):

.. code-block:: console

    $ shpc config set module_sys:tcl

You can also easily (manually) update any settings in the ``shpc/settings.yaml`` file:

.. code-block:: console

    $ shpc config edit

Take a look at this file for other configuration settings, and see the :ref:`getting-started` 
pages for next steps for setup and configuration, and interacting with your modules.

.. warning::

    You must have your container technology of choice installed and on your $PATH
    to install container modules.
     

Environment Modules
-------------------

If you are using `Environment Modules (tcl) <http://modules.sourceforge.net/>`_
and you find that your aliases do not expand, you can use `shopt <https://www.gnu.org/software/bash/manual/html_node/The-Shopt-Builtin.html>`_ to fix this issue:

.. code-block:: console

    $ shopt expand_aliases || true
    $ shopt -s expand_aliases


Pypi
====

The module is available in pypi as `singularity-hpc <https://pypi.org/project/singularity-hpc/>`_,
and this is primarily to have a consistent means for release, and an interface to show the package. Since the registry
files will not install and you would need to change the registry path
and module base (making it hard to update from the git remote) we do not
encourage you to install from pip unless you know exactly what you are doing.
.. _getting_started-developer-guide:

===============
Developer Guide
===============

This developer guide includes more complex interactions like contributing
registry entries and building containers. If you haven't read :ref:`getting_started-installation`
you should do that first.


Creating a Registry
===================

A registry consists of a database of local containers files, which are added
to the module system as executables for your user base. This typically means that you are a
linux administrator of your cluster, and shpc should be installed for you to use
(but your users will not be interacting with it).

The Registry Folder
-------------------

Although you likely will add custom containers, it's very likely that you
want to provide a set of core containers that are fairly standard, like Python
and other scientific packages. For this reason, Singularity Registry HPC
comes with a registry folder, or a folder with different containers and versions
that you can easily install. For example, here is a recipe for a Python 3.9.2 container
that would be installed to your modules as we showed above:

.. code-block:: yaml

    docker: python
    latest:
      3.9.2: sha256:7d241b7a6c97ffc47c72664165de7c5892c99930fb59b362dd7d0c441addc5ed
    tags:
      3.9.2: sha256:7d241b7a6c97ffc47c72664165de7c5892c99930fb59b362dd7d0c441addc5ed
      3.9.2-alpine: sha256:23e717dcd01e31caa4a8c6a6f2d5a222210f63085d87a903e024dd92cb9312fd
    filter:
    - 3.9.*
    maintainer: '@vsoch'
    url: https://hub.docker.com/_/python
    aliases:
      python: python

And then you would install the module file and container as follows:

.. code-block:: console

    $ shpc install python:3.9.2

But since latest is already 3.9.2, you could leave out the tag:

.. code-block:: console

    $ shpc install python


The module folder will be generated, with the structure discussed in the USer Guide. 
Currently, any new install will re-pull the container only if the hash is different, and only re-create the module otherwise.

Contributing Registry Recipes
-----------------------------

If you want to add a new registry file, you are encouraged to contribute it here
for others to use. You should:

1. Add the recipe to the ``registry`` folder in its logical namespace, either a docker or GitHub uri
2. The name of the recipe should be ``container.yaml``. You can use another recipe as a template, or see details in :ref:`getting_started-writing-registry-entries`
3. You are encouraged to add tests and then test with ``shpc test``. See :ref:`getting_started-commands-test` for details about testing.
4. You should generally choose smaller images (if possible) and define aliases (entrypoints) for the commands that you think would be useful.

A shell entrypoint for the container will be generated automatically.
When you open a pull request, a maintainer can apply
the ``container-recipe`` label and it will test your new or updated recipes accordingly.
Once your recipe is added to the repository, the versions will be automatically
updated with a nightly run. This means that you can pull the repository to get
updated recipes, and then check for updates (the bot to do this is not developed yet):


.. code-block:: console

    $ shpc check python
    ==> You have python 3.7 installed, but the latest is 3.8. Would you like to install?
    yes/no : yes


It's reasonable that you can store your recipes alongside these files, in the ``registry``
folder. If you see a conflict and want to request allowing for a custom install path
for recipes, please open an issue.


.. _getting_started-writing-registry-entries:


Writing Registry Entries
========================

An entry in the registry is a container.yaml file that lives in the ``registry``
folder. You should create subfolders based on a package name. Multiple versions
will be represented in the same file, and will install to the admin user's module
folder with version subfolders. E.g., two registry entries, one for python
(a single level name) and for tensorflow (a more nested name) would look like
this:

.. code-block:: console

    registry/
    ├── python
    │       └── container.yaml
    └── tensorflow
        └── tensorflow
            └── container.yaml


And this is what gets installed to the modules folder, where each is kept in
a separate directory based on version.

.. code-block:: console

    $ tree modules/
    modules/
    └── python
        └── 3.9.2
            ├── module.lua
            └── python-3.9.2.sif

    2 directories, 2 files

So different versions could exist alongside one another.

Registry Yaml Files
===================

Docker Hub
----------

The typical registry yaml file will reference a container from a registry,
one or more versions, and a maintainer GitHub alias that can be pinged
for any issues:


.. code-block:: yaml

    docker: python
    latest:
      3.9.2-slim: "sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef"
    tags:
      3.9.2-slim: "sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef"
      3.9.2-alpine: "sha256:23e717dcd01e31caa4a8c6a6f2d5a222210f63085d87a903e024dd92cb9312fd"
    filter:
      - "3.9.*"
    maintainer: "@vsoch"
    url: https://hub.docker.com/_/python
    aliases:
      python: /usr/local/bin/python


The above shows the simplest form of representing an alias, where each is
a key (python) and value (/usr/local/bin/python) set.


Aliases
-------

Each recipe has an optional section for defining aliases in the modulefile; there are two ways of defining them. In the python sample recipe above the simple form is used, using key value pairs:

.. code-block:: yaml

    aliases:
      python: /usr/local/bin/python

This format is container technology agnostic, because the command (``python``) and executable it targets (``/usr/local/bin/python``) would be consistent between
Podman and Singularity, for example. A second form is allowed, using dicts, in those cases where the command requires to specify custom options for the container runtime. For instance, suppose the python interpreter above requires an isolated shell environment (``--cleanenv`` in Singularity):

.. code-block:: yaml

    aliases:
    - name: python
      command: /usr/local/bin/python
      singularity_options: --cleanenv


Or perhaps the container required the docker options ``-it`` because it was an interactive, terminal session:

.. code-block:: yaml

    aliases:
    - name: python
      command: /usr/local/bin/python
      docker_options: -it


For each of the above, depending on the prefix of options that you choose, it will write them into the module files for Singularity and Docker, respectively.
This means that if you design a new registry recipe, you should consider how to run it for both kinds of technology. Also note that ``docker_options`` are
those that will also be used for Podman.


Environment Variables
---------------------

Finally, each recipe has an optional section for environment variables. For
example, the container ``vanessa/salad`` shows definition of one environment
variable:

.. code-block:: yaml

    docker: vanessa/salad
    url: https://hub.docker.com/r/vanessa/salad
    maintainer: '@vsoch'
    description: A container all about fork and spoon puns.
    latest:
      latest: sha256:e8302da47e3200915c1d3a9406d9446f04da7244e4995b7135afd2b79d4f63db
    tags:
      latest: sha256:e8302da47e3200915c1d3a9406d9446f04da7244e4995b7135afd2b79d4f63db
    aliases:
      salad: /code/salad
    env:
      maintainer: vsoch

And then during build, this variable is written to a ``99-shpc.sh`` file that
is mounted into the countainer. For the above, the following will be written:

.. code-block:: console

    export maintainer=vsoch

If a recipe does not have environment variables in the container.yaml, you have
two options for adding a variable after install. For a more permanent solution,
you can update the container.yaml file and install again. The container won't
be re-pulled, but the environment file will be re-generated. If you want to 
manually add them to the container, each module folder will have an environment
file added regardless of having this section or not, so you can export them there.
When you shell, exec, or run the container (all but inspect) you should be able
to see your environment variables:

.. code-block:: console

    $ echo $maintainer
    vsoch


Oras
----

As of version 0.0.39 Singularity Registry HPC has support for oras, meaning
we can use the Singularity client to pull an oras endpoint. Instead of using
``docker:`` in the recipe, the container.yaml might look like this:

.. code-block:: yaml

    oras: ghcr.io/singularityhub/github-ci
    url: https://github.com/singularityhub/github-ci/pkgs/container/github-ci
    maintainer: '@vsoch'
    description: An example SIF on GitHub packages to pull with oras
    latest:
      latest: sha256:227a917e9ce3a6e1a3727522361865ca92f3147fd202fa1b2e6a7a8220d510b7
    tags:
      latest: sha256:227a917e9ce3a6e1a3727522361865ca92f3147fd202fa1b2e6a7a8220d510b7


And then given the ``container.yaml`` file located in ``registry/ghcr.io/singularityhub/github-ci/`` 
you would install with shpc and the Singularity container backend as follows:

.. code-block:: console

    $ shpc install ghcr.io/singularityhub/github-ci


**Important**: You should retrieve the image sha from the container registry and 
not from the container on your computer, as the two will often be different depending
on metadata added.

Singularity Deploy
------------------

Using `Singularity Deploy <https://github.com/singularityhub/singularity-deploy>`_
you can easily deploy a container as a GitHub release! See the repository for
details. The registry entry should look like:

.. code-block:: yaml

    gh: singularityhub/singularity-deploy
    latest:
      salad: "0.0.1"
    tags:
      salad: "0.0.1"
    maintainer: "@vsoch"
    url: https://github.com/singularityhub/singularity-deploy
    aliases:
      salad: /code/salad

Where ``gh`` corresponds to the GitHub repository, the tags are the
extensions of your Singularity recipes in the root, and the "versions"
(e.g., 0.0.1) are the release numbers. There are examples in the registry
(as shown above) for details.


Choosing Containers to Contribute
---------------------------------

How should you choose container bases to contribute? You might consider using
smaller images, when possible (take advantage of multi-stage builds) and
for aliases, make sure (if possible) that you use full paths. If there is a
directive that you need for creating the module file that isn't there, please
open an issue so it can be added. Finally, if you don't have time to contribute directly, suggesting an idea via an issue or Slack to a maintainer (@vsoch).


Registry Yaml Fields
====================

Fields include:

.. list-table:: Title
   :widths: 25 65 10
   :header-rows: 1

   * - Name
     - Description
     - Required
   * - docker
     - A Docker uri, which should include the registry but not tag
     - true
   * - tags
     - A list of available tags
     - true
   * - latest
     - The latest tag, along with the digest that will be updated by a bot in the repository (e.g., tag: digest)
     - true
   * - maintainer
     - The GitHub alias of a maintainer to ping in case of trouble
     - true
   * - filter
     - A list of patterns to use for adding new tags. If not defined, all are added 
     - false
   * - aliases
     - Named entrypoints for container (dict)
     - false
   * - url
     - Documentation or other url for the container uri
     - false
   * - description
     - Additional information for the registry entry
     - false
   * - env
     - A list of environment variables to be defined in the container (key value pairs, e.g. var: value)
     - false
   * - features
     - Optional key, value paired set of features to enable for the container. Currently allowed keys: *gpu* *home* and *x11*.
     - varies


A complete table of features is shown here. The

Fields include:

.. list-table:: Title
   :widths: 20 20 20 10 10 10
   :header-rows: 1

   * - Name
     - Description
     - Container.yaml Values
     - Settings.yaml Values
     - Default
     - Supported
   * - gpu
     - Add flags to the container to enable GPU support (typically amd or nvidia)
     - true or false
     - null, amd, or nvidia
     - null
     - Singularity
   * - x11
     - Indicate to bind an Xauthority file to allow x11
     - true or false
     - null, true (uses default ~/.Xauthority) or bind path
     - null
     - Singularity
   * - home
     - Indicate a custom home to bind
     - true or false
     - null, or path to a custom home
     - null
     - Singularity, Docker


For bind paths (e.g., home and x11) you can do a single path to indicate the same
source and destination (e.g., /my/path) or a double for customization of that (e,g., /src:/dest).
Other supported (but not yet developed) fields could include different unique
resource identifiers to pull/obtain other kinds of containers. For this
current version, since we are assuming HPC and Singularity, we will typically
pull a Docker unique resource identifier with singularity, e.g.,:


.. code-block:: console

    $ singularity pull docker://python:3.9.2


Updating Registry Yaml Files
============================

We will be developing a GitHub action that automatically parses new versions
for a container, and then updates the registry packages. The algorithm we will
use is the following:

 - If docker, retrieve all tags for the image
 - Update tags:
   - if one or more filters ("filter") are defined, add new tags that match
   - otherwise, add all new tags
 - If latest is defined and a version string can be parsed, update latest
 - For each of latest and tags, add new version information


.. _getting_started-development:

Development or Testing
======================

If you first want to test singularity-hpc (shpc) with an Lmod installed in 
a container, a ``Dockerfile`` is provided for Lmod, and ``Dockerfile.tcl``
for tcl modules. The assumption is that
you have a module system installed on your cluster or in the container. If not, you
can find instructions `here for lmod <https://lmod.readthedocs.io/en/latest/030_installing.html>`_
or `here for tcl <https://modules.readthedocs.io/en/latest/INSTALL.html>`_.


.. code-block:: console
    
    $ docker build -t singularity-hpc .

If you are developing the library and need the module software, you can easily bind your
code as follows:


.. code-block:: console

    $ docker run -it --rm -v $PWD/:/code singularity-hpc

Once you are in the container, you can direct the module software to use your module files:

.. code-block:: console

    $ module use /code/modules

Then you can use spider to see the modules:

.. code-block:: console

    # module spider python

    --------------------------------------------------------------------------------------------------------------------------------------------------------------
      python/3.9.2: python/3.9.2/module
    --------------------------------------------------------------------------------------------------------------------------------------------------------------

        This module can be loaded directly: module load python/3.9.2/module
    ```


or ask for help directly!

.. code-block:: console

    # module help python/3.9.2-slim

    ----------------------------------------------------- Module Specific Help for "python/3.9.2-slim/module" ------------------------------------------------------
    This module is a singularity container wrapper for python v3.9.2-slim


    Container:

     - /home/vanessa/Desktop/Code/singularity-hpc/modules/python/3.9.2-slim/python-3.9.2-slim-sha256:85ed629e6ff79d0bf796339ea188c863048e9aedbf7f946171266671ee5c04ef.sif

    Commands include:

     - python-run:
           singularity run <container>
     - python-shell:
           singularity shell -s /bin/bash <container>
     - python-exec:
           singularity exec -s /bin/bash <container> "$@"
     - python-inspect-runscript:
           singularity inspect -r <container>
     - python-inspect-deffile:
           singularity inspect -d <container>

     - python:
           singularity exec <container> /usr/local/bin/python"


    For each of the above, you can export:

     - SINGULARITY_OPTS: to define custom options for singularity (e.g., --debug)
     - SINGULARITY_COMMAND_OPTS: to define custom options for the command (e.g., -b)


Note that you typically can't run or execute containers within another container, but 
you can interact with the module system. Also notice that for every container, we expose easy
commands to shell, run, exec, and inspect. The custom commands (e.g., Python) are then provided below that.

Make sure to write to files outside of the container so you don't muck with permissions.
Since we are using module use, this means that you can create module files as a user
or an admin - it all comes down to who has permission to write to the modules
folder, and of course use it. Note that I have not tested this on an HPC system
but plan to shortly.
.. _getting-started:

===============
Getting Started
===============

Singularity Registry (HPC) is a tool that makes it easy to install containers as
Lmod modules. You can create your own registry entries (e.g., a specification
to pull a particular container and expose some number of entrypoints) or
the library also provides you with a community set.

If you have any questions or issues, please `let us know <https://github.com/singularityhub/singularity-hpc/issues>`_

.. toctree::
   :maxdepth: 2

   installation
   user-guide
   developer-guide
   use-cases
