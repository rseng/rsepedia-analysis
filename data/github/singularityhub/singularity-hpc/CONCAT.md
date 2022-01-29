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
