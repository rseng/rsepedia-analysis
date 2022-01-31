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

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at info@nextflow.io. All
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

# GIT README 

Shortcuts to recurrent Git commands 

## Sync remote fork 

```bash
git fetch public
git checkout master
git merge public/master
```

List branch remote

    git branch -vv 

Checkout from a remote branch 

    git co -b <local name> upstream/master

Push to a remote upstream branch

    git push <remote> <local branch>:<remote branch>

eg:

    git push upstream foo:master

Read more [here](https://help.github.com/articles/syncing-a-fork/).

## Pull options

  git config pull.rebase false  # merge (the default strategy)
  git config pull.rebase true   # rebase
  git config pull.ff only       # fast-forward only


## Subtree  

The `tests` directory is a Git subtree created with the 
following commands: 

    git remote add tests git@github.com:nextflow-io/tests.git
    git subtree add --squash --prefix=tests/ tests integration


To pull changes from the [tests repo](https://github.com/nextflow-io/tests) use this command: 

    git subtree pull --squash --prefix=tests/ tests integration

To push changes to the [tests repo](https://github.com/nextflow-io/tests) use this command: 

    git subtree push --prefix=tests/ tests integration


Read more [here](https://andrey.nering.com.br/2016/git-submodules-vs-subtrees/).

## Stash shortcuts

    git stash list
    git stash pop
    git stash pop stash@{1}
    git showtool stash@{0}
    git stash drop
    git stash drop stash@{1}
    git stash clear
    git diff stash
    git diff stash@{1} [other]

## Misc 

Find a commit in any branch introducing a change

    git log -S <whatever> --source --all

Reset last merge pushed 

    git reset --hard HEAD@{1}

    Read more https://stackoverflow.com/a/11722640/395921
    
## GPG keys 

To sign Git commits with a GPG key on Mac use [GPG Suite](https://gpgtools.org/), import your key, then: 

    git config --global gpg.program /usr/local/MacGPG2/bin/gpg2
    git config --global user.signingkey <your key> 
    git config --global commit.gpgsign true 
    git config --global format.signoff true ## TO AVOID TO SPECIFY -S option each time

Read more: 
https://gist.github.com/danieleggert/b029d44d4a54b328c0bac65d46ba4c65


## Change Git history root 

https://stackoverflow.com/questions/4515580/how-do-i-remove-the-old-history-from-a-git-repository

# Pull master while in a another branch

    git fetch origin master:master

# Alternative to `git pull` 

Use `git fetch && git rebase`

![Nextflow logo](https://github.com/nextflow-io/trademark/blob/master/nextflow2014_no-bg.png)

*"Dataflow variables are spectacularly expressive in concurrent programming"*
<br>[Henri E. Bal , Jennifer G. Steiner , Andrew S. Tanenbaum](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.145.7873)


![Nextflow CI](https://github.com/nextflow-io/nextflow/workflows/Nextflow%20CI/badge.svg)
[![Nextflow version](https://img.shields.io/github/release/nextflow-io/nextflow.svg?colorB=26af64&style=popout)](https://github.com/nextflow-io/nextflow/releases/latest)
[![Chat on Gitter](https://img.shields.io/gitter/room/nextflow-io/nextflow.svg?colorB=26af64&style=popout)](https://gitter.im/nextflow-io/nextflow)
[![Nextflow Twitter](https://img.shields.io/twitter/url/https/nextflowio.svg?colorB=26af64&&label=%40nextflow&style=popout)](https://twitter.com/nextflowio)
[![Nextflow Publication](https://img.shields.io/badge/Published-Nature%20Biotechnology-26af64.svg?colorB=26af64&style=popout)](https://www.nature.com/articles/nbt.3820)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?colorB=26af64&style=popout)](http://bioconda.github.io/recipes/nextflow/README.html)
[![Nextflow license](https://img.shields.io/github/license/nextflow-io/nextflow.svg?colorB=26af64&style=popout)](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

Quick overview
==============
Nextflow is a bioinformatics workflow manager that enables the development of portable and reproducible workflows.
It supports deploying workflows on a variety of execution platforms including local, HPC schedulers, AWS Batch,
Google Cloud Life Sciences, and Kubernetes. Additionally, it provides support for manage your workflow dependencies
through built-in support for Conda, Docker, Singularity, and Modules.

## Contents
- [Rationale](#rationale)
- [Quick start](#quick-start)
- [Documentation](#documentation)
- [Tool Management](#tool-management)
  - [Conda environments](#conda-environments)
  - [Docker and Singularity](#containers)
  - [Environment Modules](#environment-modules)
- [HPC Schedulers](#hpc-schedulers)
  - [SGE](#hpc-schedulers)
  - [Univa Grid Engine](#hpc-schedulers)
  - [LSF](#hpc-schedulers)
  - [SLURM](#hpc-schedulers)
  - [PBS/Torque](#hpc-schedulers)
  - [HTCondor (experimental)](#hpc-schedulers)
- [Cloud Support](#cloud-support)
  - [AWS Batch](#cloud-support)
  - [AWS EC2](#cloud-support)
  - [Google Cloud](#cloud-support)
  - [Google Genomics Pipelines](#cloud-support)
  - [Kubernetes](#cloud-support)
- [Community](#community)
- [Build from source](#build-from-source)
- [Contributing](#contributing)
- [License](#license)
- [Citations](#citations)
- [Credits](#credits)


Rationale
=========

With the rise of big data, techniques to analyse and run experiments on large datasets are increasingly necessary.

Parallelization and distributed computing are the best ways to tackle this problem, but the tools commonly available to the bioinformatics community often lack good support for these techniques, or provide a model that fits badly with the specific requirements in the bioinformatics domain and, most of the time, require the knowledge of complex tools or low-level APIs.

Nextflow framework is based on the dataflow programming model, which greatly simplifies writing parallel and distributed pipelines without adding unnecessary complexity and letting you concentrate on the flow of data, i.e. the functional logic of the application/algorithm.

It doesn't aim to be another pipeline scripting language yet, but it is built around the idea that the Linux platform is the *lingua franca* of data science, since it provides many simple command line and scripting tools, which by themselves are powerful, but when chained together facilitate complex data manipulations.

In practice, this means that a Nextflow script is defined by composing many different processes. Each process can execute a given bioinformatics tool or scripting language, to which is added the ability to coordinate and synchronize the processes execution by simply specifying their inputs and outputs.



Quick start
============

Download the package
--------------------

Nextflow does not require any installation procedure, just download the distribution package by copying and pasting
this command in your terminal:

```
curl -fsSL https://get.nextflow.io | bash
```

It creates the ``nextflow`` executable file in the current directory. You may want to move it to a folder accessible from your ``$PATH``.

Download from Conda
-------------------

Nextflow can also be installed from Bioconda

```
conda install -c bioconda nextflow 
```

Documentation
=============

Nextflow documentation is available at this link http://docs.nextflow.io


HPC Schedulers
==============

*Nextflow* supports common HPC schedulers, abstracting the submission of jobs from the user. 

Currently the following clusters are supported:

  + [SGE](https://www.nextflow.io/docs/latest/executor.html#sge)
  + [Univa Grid Engine](https://www.nextflow.io/docs/latest/executor.html#sge)
  + [LSF](https://www.nextflow.io/docs/latest/executor.html#lsf)
  + [SLURM](https://www.nextflow.io/docs/latest/executor.html#slurm)
  + [PBS/Torque](https://www.nextflow.io/docs/latest/executor.html#pbs-torque)
  + [HTCondor (beta)](https://www.nextflow.io/docs/latest/executor.html#htcondor)
  + [Moab (beta)](https://www.nextflow.io/docs/latest/executor.html#moab)

For example to submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

```nextflow
process {
  executor='sge'
  queue='<your execution queue>'
}
```

In doing that, processes will be executed by Nextflow as SGE jobs using the `qsub` command. Your 
pipeline will behave like any other SGE job script, with the benefit that *Nextflow* will 
automatically and transparently manage the processes synchronisation, file(s) staging/un-staging, etc.  


Cloud support
=============
*Nextflow* also supports running workflows across various clouds and cloud technologies. Managed solutions from major 
cloud providers are also supported through AWS Batch, Azure Batch and Google Cloud compute services. 
Additionally, *Nextflow* can run workflows on either on-prem or managed cloud Kubernetes clusters. 

Currently supported cloud platforms:
  + [AWS Batch](https://www.nextflow.io/docs/latest/awscloud.html#aws-batch)
  + [Azure Batch](https://azure.microsoft.com/en-us/services/batch/)
  + [Google Cloud Life Sciences](https://cloud.google.com/life-sciences)
  + [Kubernetes](https://www.nextflow.io/docs/latest/kubernetes.html)



Tool management
================

Containers
----------------

*Nextflow* has first class support for containerization. It supports both [Docker](https://www.nextflow.io/docs/latest/docker.html) and [Singularity](https://www.nextflow.io/docs/latest/singularity.html) container engines. Additionally, *Nextflow* can easily switch between container engines enabling workflow portability. 

```nextflow
process samtools {
  container 'biocontainers/samtools:1.3.1'

  """
  samtools --version 
  """

}
```

Conda environments
------------------

[Conda environments](https://www.nextflow.io/docs/latest/conda.html) provide another option for managing software packages in your workflow. 


Environment Modules
-------

[Environment modules](https://www.nextflow.io/docs/latest/process.html#module) commonly found in HPC environments can also be used to manage the tools used in a *Nextflow* workflow. 


Community
=========

You can post questions, or report problems by using the Nextflow [discussion forum](https://groups.google.com/forum/#!forum/nextflow)
or the [Nextflow channel on Gitter](https://gitter.im/nextflow-io/nextflow).

*Nextflow* also hosts a yearly workshop showcasing researcher's workflows and advancements in the langauge. Talks from the past workshops are available on the [Nextflow YouTube Channel](https://www.youtube.com/channel/UCB-5LCKLdTKVn2F4V4KlPbQ)

The [nf-core](https://nf-co.re/) project is a community effort aggregating high quality *Nextflow* workflows which can be used by the community. 


Build from source
=================

Required dependencies
---------------------

* Compiler Java 8 or later
* Runtime Java 8 or later

Build from source
-----------------

*Nextflow* is written in [Groovy](http://groovy-lang.org) (a scripting language for the JVM). A pre-compiled,
ready-to-run, package is available at the [Github releases page](https://github.com/nextflow-io/nextflow/releases),
thus it is not necessary to compile it in order to use it.

If you are interested in modifying the source code, or contributing to the project, it worth knowing that
the build process is based on the [Gradle](http://www.gradle.org/) build automation system.

You can compile *Nextflow* by typing the following command in the project home directory on your computer:

```bash
make compile
```

The very first time you run it, it will automatically download all the libraries required by the build process.
It may take some minutes to complete.

When complete, execute the program by using the `launch.sh` script in the project directory.

The self-contained runnable Nextflow packages can be created by using the following command:

```bash
make pack
```            

Once compiled use the script `./launch.sh` as a replacement for the usual `nextflow` command.

The compiled packages can be locally installed using the following command:

```bash
make install
```

A self-contained distribution can be created with the command: `make pack`.  To include support of GA4GH and its dependencies in the binary, use `make packGA4GH` instead.


IntelliJ IDEA
---------------

Nextflow development with [IntelliJ IDEA](https://www.jetbrains.com/idea/) requires the latest version of the IDE (2019.1.2 or later).

If you have it installed in your computer, follow the steps below in order to use it with Nextflow:

1. Clone the Nextflow repository to a directory in your computer.
2. Open IntelliJ IDEA and choose "Import project" in the "File" menu bar.
3. Select the Nextflow project root directory in your computer and click "OK".
4. Then, choose the "Gradle" item in the "external module" list and click on "Next" button.
5. Confirm the default import options and click on "Finish" to finalize the project configuration.
6. When the import process complete, select the "Project structure" command in the "File" menu bar.
7. In the showed dialog click on the "Project" item in the list of the left, and make sure that
   the "Project SDK" choice on the right contains Java 8.
8. Set the code formatting options with setting provided [here](https://github.com/nextflow-io/nextflow/blob/master/CONTRIBUTING.md#ide-settings).



Contributing
============

Project contribution are more than welcome. See the [CONTRIBUTING](CONTRIBUTING.md) file for details.


Build servers
=============

  * [Travis-CI](https://travis-ci.org/nextflow-io/nextflow)
  * [GitHub Actions](https://github.com/nextflow-io/nextflow/actions)

License
=======

The *Nextflow* framework is released under the Apache 2.0 license.

Citations
=========

If you use Nextflow in your research, please cite:

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:[10.1038/nbt.3820](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3820.html)

Credits
=======

Nextflow is built on two great pieces of open source software, namely <a href='http://groovy-lang.org' target='_blank'>Groovy</a>
and <a href='http://www.gpars.org/' target='_blank'>Gpars</a>.

YourKit is kindly supporting this open source project with its full-featured Java Profiler.
Read more http://www.yourkit.com
# CONTRIBUTING TO NEXTFLOW

This guide documents the best way to make various types of contributions to Nextflow,
including what is required before submitting a code change.

Contributing to Nextflow doesn't just mean writing code. Helping new users on the mailing list,
testing releases and bug fixes, and improving documentation are all essential and valuable contributions. In fact, proposing
significant code changes usually first requires gaining experience and credibility within the
community by helping in other ways. This is also a guide to becoming an effective contributor.


## Contributing by Helping Other Users

A great way to contribute to Nextflow is to help answer user questions on the [discussion forum](https://groups.google.com/forum/#!forum/nextflow)
or the [Gitter channel](https://gitter.im/nextflow-io/nextflow). There are always many new Nextflow users;
taking a few minutes to help answer a question is a very valuable community service.

Contributors should ideally subscribe to these channels and follow them in order to keep up to date
on what's happening in Nextflow. Answering questions is an excellent and visible way to help the
community and also demonstrates your expertise.


## Contributing Documentation Changes

To propose a change to release documentation (that is, the docs that appear under http://docs.nextflow.io),
edit the documentation source files in Nextflow's [docs/](https://github.com/nextflow-io/nextflow/tree/master/docs)
directory, whose README file shows how to build the documentation locally to test your changes.

Then open a pull request with the proposed changes.


## Contributing Bug Reports

Filling a bug report is likely the simplest and most useful way to contribute to the project.
It helps us to identify issues and provide patches and therefore to make Nextflow more stable
and useful.

Report a bug using the "New issue" button in the
[issues page](https://github.com/nextflow-io/nextflow/issues) of this project.

A good bug report should include a minimal executable test case able to replicate the reported bug.

Follow the instructions in the bug [report template](https://github.com/nextflow-io/nextflow/blob/master/.github/issue_template.md) that is shown when filling the bug report out.

## Contributing Bug Fixes

Contributing bug fixes is the best way to gain experience and credibility within the community
and also to become an effective project contributor.

If you are a novice with the Nextflow code base, start by looking at issues marked
with the [help wanted](https://github.com/nextflow-io/nextflow/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22)
label.

If you have doubts on how to fix an issue, ask for help from senior contributors commenting
in the issue page.

## Contributing New Features

Before contributing a new feature, submit a new feature proposal in the
[issues page](https://github.com/nextflow-io/nextflow/issues) of the project and discuss it
with the community.

This is important to identify possible overlaps with other planned features and avoid misunderstandings and conflicts.

## Contributing Code Changes

When submitting a contribution, you will be required to sign a [Developer Certificate of Origin (DCO)](https://developercertificate.org/) to certify that you are the author of the source code or otherwise you have the right to submit it to the project. 

Contributor signatures are provided by adding a `Signed-off-by` line to the commit message 
as shown below, or by using the `-s` option for the [git commit command](https://help.github.com/articles/signing-commits/).

```
This is my commit message

Signed-off-by: Random J Developer <random@developer.example.org>
```

The process is automatically managed by the [Probot](https://probot.github.io/apps/dco/) app for GitHub.


## IDE settings

The suggested development environment is [IntelliJ IDEA](https://www.jetbrains.com/idea/download/). See the [README](https://github.com/nextflow-io/nextflow/#intellij-idea) for a short primer on how to import
and configure Nextflow to work with it.

Nextflow does not impose a strict code formatting style, however the following setting should be applied:

* Use spaces for indentation
* Tab size: 4
* Indent: 4
* Use single class import
* Class count to use import with `*`: 99
* Names count to use static import with `*`: 99
* Imports layout:
    * \<blank line>
    * `import org.junit.*`
    * `import spock.lang.*`
    * \<blank line>
    * `import java.*`
    * `import javax.*`
    * \<blank line>
    * *all other imports*
    * *all other static imports*

New files must include the appropriate license header boilerplate and the author name(s) and contact email(s) ([see for example](https://github.com/nextflow-io/nextflow/blob/e8945e8b6fc355d3f2eec793d8f288515db2f409/modules/nextflow/src/main/groovy/nextflow/Const.groovy#L1-L15)).
# Nextflow tests

This repository contains a collection of scripts used to validate Nextflow builds 

[![Build Status](https://travis-ci.org/nextflow-io/tests.svg?branch=master)](https://travis-ci.org/nextflow-io/tests)
Hi! Thanks for contributing to Nextflow project.

When submitting a Pull Request please make sure to not include
in the changeset any modification in these files:

* `nextflow`
* `docs/conf.py`
* `modules/nf-commons/src/main/nextflow/Const.groovy`

Also, please sign-off the DCO [1] to certify you are the author of the contribution
and you adhere to Nextflow open source license [2] adding a `Signed-off-by` line to
the contribution commit message. For more details check [3].

1. https://developercertificate.org/
2. https://github.com/nextflow-io/nextflow/blob/master/COPYING
3. https://github.com/apps/dco

# Action 

## Syntax

https://help.github.com/en/articles/workflow-syntax-for-github-actions
https://help.github.com/en/articles/contexts-and-expression-syntax-for-github-actions
https://help.github.com/en/articles/virtual-environments-for-github-actions#environment-variables
https://help.github.com/en/articles/configuring-docker-for-use-with-github-package-registry
https://help.github.com/en/articles/virtual-environments-for-github-actions#creating-and-using-secrets-encrypted-variables

## Java 

Java VMs has to match the ones at this link https://static.azul.com/zulu/bin 

Check the name *-jdk(x.y.z)
---
name: General question 
about: Need for help on Nextflow language and usage
---

Hi! Thanks for using Nextflow. 

If you need help about Nextflow scripting language, 
configuration options and general Nextflow usage the better 
channels to post this kind of questions are: 

* GitHub discussions: https://github.com/nextflow-io/nextflow/discussions
* Google group: https://groups.google.com/forum/#!forum/nextflow
* Gitter channel: https://gitter.im/nextflow-io/nextflow


Also you may also want to have a look at the patterns page 
for common solutions to recurrent implementation problems: 
http://nextflow-io.github.io/patterns/index.html

---
name: New feature
about: Submit a new feature proposal
---

## New feature

Hi! Thanks for using Nextflow and submitting the proposal 
for a new feature or the enhancement of an existing functionality. 

Please replace this text providing a short description of your 
proposal.

## Usage scenario 

(What's the main usage case and the deployment scenario addressed by this proposal)

## Suggest implementation 

(Highlight the main building blocks of a possible implementation and/or related components)


---
name: Bug report
about: Create a report to help us improve
---

## Bug report 

(Please follow this template replacing the text between parentheses with the requested information)

### Expected behavior and actual behavior

(Give an brief description of the expected behavior 
and actual behavior)

### Steps to reproduce the problem

(Provide a test case that reproduce the problem either with a self-contained script or GitHub repository)

### Program output 

(Copy and paste here output produced by the failing execution. Please highlight it as a code block. Whenever possible upload the `.nextflow.log` file.)

### Environment 

* Nextflow version: [?] 
* Java version: [?]
* Operating system: [macOS, Linux, etc]
* Bash version: (use the command `$SHELL --version`)

### Additional context

(Add any other context about the problem here)
# Nextflow plugins

This directory should contain plugin subprojects for Nextflow. 

The `build.gradle` defines the main actions for each plugin.

## Plugin subproject structure 
 
### Plugin structure 

The plugin subproject defines its own `build.gradle` and setup the required dependencies. 

Minimal dependencies shown below: 

``` 
dependencies {
    compileOnly project(':nextflow')
    compileOnly 'org.slf4j:slf4j-api:1.7.10'
    compileOnly 'org.pf4j:pf4j:3.4.1'

    testImplementation project(':nextflow')
    testImplementation "org.codehaus.groovy:groovy:3.0.5"
    testImplementation "org.codehaus.groovy:groovy-nio:3.0.5"
}
``` 

Each plugin subproject directory name has to begin with the prefix `nf-` and must include 
a file named `src/resources/META-INF/MANIFEST.MF` which contains the plugin metadata. 
The manifest content looks like the following:

```
Manifest-Version: 1.0
Plugin-Class: the.plugin.ClassName
Plugin-Id: the-plugin-id
Plugin-Provider: Some Provider Name
Plugin-Version: 0.0.0
```   
  
## Environment variables 

* `NXF_PLUGINS_MODE`: Define the plugin system execution mode, either *prod* for production or *dev* for development
  (see below for details).  
* `NXF_PLUGINS_DIR`: the path where the plugins archives are stored/loaded. Default is `$NXF_HOME/plugins` for 
  production mode and `$PWD/plugins` for dev mode. 
* `NXF_PLUGINS_DEFAULT`: Whenever use the default plugins when no plugin is specified in the config file.   
* `NXF_PLUGINS_DEV: Comma separate separated list of development plugins root directories

## Development environment

When running in development the plugin system uses the `DevPluginClasspath` to load plugins classes 
from each plugin project build path e.g. `$PWD/plugins/nf-amazon/build/classes` and 
`$PWD/plugins/nf-amazon/build/target/libs` (for deps libraries).    


## The plugins repository 

The plugins meta-info are published via a GitHub repository at [https://github.com/nextflow-io/plugins]()
and accessible through the URL [https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json]().

The repository index has the following structure: 

```
[
  {
    "id": "nf-amazon",
    "releases": [
      {
        "version": "0.2.0",
        "url": "https://github.com/nextflow-io/nf-amazon/releases/download/0.2.0/nf-amazon-0.2.0.zip",
        "date": "2020-10-12T10:05:44.28+02:00",
        "sha512sum": "9e9e33695c1a7c051271..."
      }
    ]
  },
  :
]
```     


## Plugins 

A plugin is a ZIP file holding either the plugin classes and the required dependencies JAR file.  

Nextflow core plugins are stored in the corresponding GitHub project release page. However, this is not a
strict requirement, it has been chosen to simplify the build deployment process and provide a more consistent 
download experience keeping all of them with the GitHub [nextflow-io](https://github.com/nextflow-io) organization.   

### The installation process 

Plugins need to be declared in the `nextflow.config` file using the plugins scope, eg. 

```
plugins {
    id 'nf-amazon@0.2.0'
}
```     

If the plugins is not locally available Nextflow check in the repository index for the download URL, 
download the ZIP in a temporary file, unzip and store the final plugin in the directory specified by the 
variable `NXF_PLUGINS_DIR` (default: `$NXF_HOME/plugins`). 

Finally, since each Nextflow run can have a different set of plugins (and version) requirement, each Nextflow 
instance keep local plugins directory root in directory `$PWD/.nextflow/plr/<unique id>` symlinking the exact list 
of plugins directory required for the current Nextflow instance.

If no plugins are specified in the nextflow.config file, Nextflow default plugins are automatically added. 
The default plugins list is defined in the Nextflow resources file included in the distribution runtime 
`./modules/nextflow/src/main/resources/META-INF/plugins-info.txt`. 

To disable the use of defualt plugins set the following variable `NXF_PLUGINS_DEFAULT=false`.

## Gradle Tasks 

### makeZip
    
Creates the plugin the zip file and the json meta file in the
subproject `build/libs` directory.

```
» ls -l1 $PWD/plugins/nf-tower/build/libs/
nf-tower-0.1.0.jar
nf-tower-0.1.0.json
nf-tower-0.1.0.zip
```               

### copyPluginLibs

Copies plugin dependencies jar files in the plugin build directory ie. `$PWD/plugins/nf-amazon/build/target/libs`. 
This is only needed when launching the plugin in *development* mode. 

### copyPluginZip

Copies the plugin zip file to the root project build dir ie. `$PWD/build/plugins/`.

### uploadPlugin

Uploads the plugin ZIP and meta (JSON) files to the corresponding GitHub repository. Options: 

* `release`: the plugin version e.g. `1.0.1`
* `repo`: the GitHub repository name e.g. `nf-amazon`
* `owner`: the GitHub owning organization e.g. `nextflow-io`
* `skipExisting`: do not upload a file already existing (only if the checksum is the same, default: `true`). 
* `dryRun`: execute the tasks without uploading file (default: `false`).
* `overwrite`: prevent to overwrite a remote file already existing (default: `false`).
* `userName`: the user name for authenticate GitHub API requests
* `authToken`: the personal token to authenticate GitHub API requests  

### upload

Uploads the plugin both zip and meat files. 

### publishIndex

Upload the plugins index to the repository hosted at [https://github.com/nextflow-io/plugins](https://github.com/nextflow-io/plugins), which makes 
them accessible through the URL [https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json). 


## Links 

* https://pf4j.org/
* https://proandroiddev.com/understanding-gradle-the-build-lifecycle-5118c1da613f
# Azure plugin for Nextflow 

This plugin implements the support for Azure Blob storage as fie system 
provider (via JSR203 interface) and Azure Batch executor  for Nextflow 

## Compile & configuration 

Compile nextflow as usual 

```
make compile
``` 

Create the `nextflow.config` file with the following content 

```
plugins { 
  id 'nf-azure'
}

azure {
  storage {
    accountKey = "<YOUR STORAGE ACCOUNT KEY>"
    accountName = "<YOUR STORAGE ACCOUNT KEY>"
  }

  batch {
    endpoint = 'https://<YOUR BATCH ACCOUNT NAME>.westeurope.batch.azure.com' 
    accountName = '<YOUR BATCH ACCOUNT NAME>' 
    accountKey = '<YOUR BATCH ACCOUNT KEY>'
  }
}

process.executor = 'azurebatch'
workDir = 'az://<YOUR DATA CONTAINER>/work'
```

Then run the a pipeline as shown below

```
./launch.sh run rnaseq-nf 
```


## Todo 

* Currently, the Blob storage service uses NettyHttpClient and Batch service 
uses OkHttp client, duplicating the number of required libraries. In principle 
the Blob service can use OkHttp, adding the following deps, however using that
Nextflow hangs during the shutdown, apparently because the connection pool used 
by the blob service is not closed timely. 

        compile('com.azure:azure-storage-blob:12.9.0') {
            exclude group: 'org.slf4j', module: 'slf4j-api'
            exclude group: 'com.azure', module: 'azure-core-http-netty'
        }
        compile('com.azure:azure-core-http-okhttp:1.3.3') {
            exclude group: 'org.slf4j', module: 'slf4j-api'
        }

* Remove invalid directory from .command.run PATH for project having `bin/` folder  
* Add the configuration for the region
* Make the backend endpoint optional 



### Links
* https://github.com/Azure/azure-sdk-for-java/wiki
* https://github.com/Azure/azure-sdk-for-java/tree/master/sdk/storage/azure-storage-blob-nio
* https://github.com/Azure/azure-sdk-for-java/blob/master/sdk/storage/azure-storage-blob-nio/src/samples/java/com/azure/storage/blob/nio/ReadmeSamples.java

# SQL DB plugin for Nextflow

This plugin provides an extension to implement built-in support for SQL DB access and manipulation in Nextflow scripts. 

It provides the ability to create a Nextflow channel from SQL queries and to populate database tables. The current version 
provides out-of-the-box support for the following databases: 

* [H2](https://www.h2database.com)
* [MySQL](https://www.mysql.com/) 
* [MariaDB](https://mariadb.org/)
* [PostgreSQL](https://www.postgresql.org/)
* [Sqlite](https://www.sqlite.org/index.html)
* [DuckDB](https://duckdb.org/)
                    
NOTE: THIS IS A PREVIEW TECHNOLOGY, FEATURES AND CONFIGURATION SETTINGS CAN CHANGE IN FUTURE RELEASES.

This repository only holds plugin artefacts. Source code is available at this [link](https://github.com/nextflow-io/nextflow/tree/master/plugins/nf-sqldb).

## Get started 
  
Make sure to have Nextflow 21.08.0 or later. Add the following snippet to your `nextflow.config` file. 

```
plugins {
  id 'nf-sqldb@0.2.0'
}
```
                                                              
The above declaration allows the use of the SQL plugin functionalities in your Nextflow pipelines. See the section 
below to configure the connection properties with a database instance. 

## Configuration

The target database connection coordinates are specified in the `nextflow.config` file using the
`sql.db` scope. The following are available

| Config option 	                    | Description 	                |
|---	                                |---	                        |
| `sql.db.'<DB-NAME>'.url`      | The database connection URL based on Java [JDBC standard](https://docs.oracle.com/javase/tutorial/jdbc/basics/connecting.html#db_connection_url). 
| `sql.db.'<DB-NAME>'.driver`   | The database driver class name (optional).
| `sql.db.'<DB-NAME>'.user`     | The database connection user name.
| `sql.db.'<DB-NAME>'.password` | The database connection password.

For example:

```
sql {
    db {
        foo {
              url = 'jdbc:mysql://localhost:3306/demo'
              user = 'my-user'
              password = 'my-password'
            }
    }
}

```

The above snippet defines SQL DB named *foo* that connects to a MySQL server running locally at port 3306 and
using `demo` schema, with `my-name` and `my-password` as credentials.

## Available operations

This plugin adds to the Nextflow DSL the following extensions that allows performing of queries and populating database tables.

### fromQuery

The `fromQuery` factory method allows for performing a query against a SQL database and creating a Nextflow channel emitting
a tuple for each row in the corresponding result set. For example:

```
ch = channel.sql.fromQuery('select alpha, delta, omega from SAMPLE', db: 'foo')
```

### sqlInsert

The `sqlInsert` operator provided by this plugin allows populating a database table with the data emitted
by a Nextflow channels and therefore produced as result by a pipeline process or an upstream operator. For example:

```
channel
    .of('Hello','world!')
    .map( it -> tuple(it, it.length) )
    .sqlInsert( into: 'SAMPLE', columns: 'NAME, LEN', db: 'foo' )

```

The above example creates and performs the following two SQL statements into the database with name `foo` as defined
in the `nextflow.config` file.

```
INSERT INTO SAMPLE (NAME, LEN) VALUES ('HELLO', 5);
INSERT INTO SAMPLE (NAME, LEN) VALUES ('WORLD!', 6);
```

NOTE: the target table (e.g. `SAMPLE` in the above example) must be created ahead.

The following options are available:

| Operator option 	| Description 	                |
|---	            |---	                        |
| `db`              | The database handle. It must must a `sql.db` name defined in the `nextflow.config` file.
| `into`            | The database table name into with the data needs to be stored.
| `columns`         | The database table column names to be filled with the channel data. The column names order and cardinality must match the tuple values emitted by the channel. The columns can be specified as a `List` object or a comma-separated value string.
| `statement`       | The SQL `insert` statement to be performed to insert values in the database using `?` as placeholder for the actual values, for example: `insert into SAMPLE(X,Y) values (?,?)`. When provided the `into` and `columsn` parameters are ignored.
| `batch`           | The number of insert statements that are grouped together before performing the SQL operations (default: `10`). 
| `setup`           | A SQL statement that's executed before the first insert operation. This is useful to create the target DB table. NOTE: the underlying DB should support the *create table if not exist* idiom (i.e. the plugin will execute this time every time the script is run).

## Query CSV files

The SQL plugin includes the [H2](https://www.h2database.com/html/main.html) database engine that allows the query of CSV files
as DB tables using SQL statements.

For example, create CSV file using the snippet below:

```
cat <<EOF > test.csv
foo,bar
1,hello
2,ciao
3,hola
4,bonjour
EOF
```

To query this file in a Nextflow script use the following snippet:

```nextflow
    channel
          .sql
          .fromQuery("SELECT * FROM CSVREAD('test.csv') where foo>=2;")
          .view()
```


The `CSVREAD` function provided by the H2 database engine allows the access of a CSV file in your computer file system,
you can replace `test.csv` with a CSV file path of your choice. The `foo>=2` condition shows how to define a filtering
clause using the conventional SQL WHERE constrains. 

## Important 

This plugin is not expected to be used to store and access a pipeline status in a synchronous manner during the pipeline 
execution. 

This means that if your script has a `sqlInsert` operation followed by a successive `fromQuery` operation, the query 
may *not* contain previously inserted data due to the asynchronous nature of Nextflow operators.

The SQL support provided by this plugin is meant to be used to fetch DB data from a previous run or to populate DB tables
for storing or archival purpose.
# Gradle plugins for Nextflow build 

The plugin is uploaded to the Gradle plugins portal. 

To check the status use [this link](https://plugins.gradle.org/u/nextflowio). 

Use the task `../gradlew publishPlugins` in the `buildSrc` directory to publish the plugin.# Nextflow Documentation 

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/) which 
uses the [reStructuredText](https://en.wikipedia.org/wiki/ReStructuredText) file format.

A quick intro to reStructuredText is available at [this link](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

To edit and contribute to the documentation you only need a text editor to change the
appropriate `.rst` files in this directory.

Once you have edited the documentation files verify the your changes are correctly applied
using the command below to generate the HTML files:

```
make html
```


Then open the `_build/html/index.html` file with your browser and navigate the documentation
you have modified.


### Dependencies

Sphinx can be installed either with

```
pip install -U Sphinx
```

or

```
conda install sphinx
```

### Theme 

Docs uses the [sphinx_rtd_theme](https://github.com/readthedocs/sphinx_rtd_theme) theme. 

To update it, clone the bove repo, then copy `sphinx_rtd_theme/sphinx_rtd_theme` directory 
into `docs/_themes`.  

### License

Nextflow documentation is distributed under 
[Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) license](https://creativecommons.org/licenses/by-sa/4.0/).
.. _awscloud-page:

************
Amazon Cloud
************

AWS security credentials
=========================

Nextflow uses the `AWS security credentials <https://docs.aws.amazon.com/general/latest/gr/aws-sec-cred-types.html>`_
to make programmatic calls to AWS services.

You can provide your AWS access keys using the standard AWS variables shown below:

    * ``AWS_ACCESS_KEY_ID``
    * ``AWS_SECRET_ACCESS_KEY``
    * ``AWS_DEFAULT_REGION``

If ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY`` are not defined in the environment, Nextflow will attempt to
retrieve credentials from your ``~/.aws/credentials`` or ``~/.aws/config`` files. The ``default`` profile can be
overridden via the environmental variable ``AWS_PROFILE`` (or ``AWS_DEFAULT_PROFILE``).

Alternatively AWS credentials can be specified in the Nextflow configuration file.

See :ref:`AWS configuration<config-aws>` for more details.

.. note:: Credentials can also be provided by using an IAM Instance Role. The benefit of this approach is that
  it spares you from managing/distributing AWS keys explicitly.
  Read the `IAM Roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_ documentation
  and `this blog post <https://aws.amazon.com/blogs/security/granting-permission-to-launch-ec2-instances-with-iam-roles-passrole-permission/>`_ for more details.

AWS IAM policies
=================

`IAM policies <https://docs.aws.amazon.com/IAM/latest/UserGuide/access_policies.html>`_ are the mechanism used by AWS to
defines permissions for IAM identities. In order to access certain AWS services, the proper policies must be
attached to the identity associated to the AWS credentials.

Minimal permissions policies to be attached to the AWS account used by Nextflow are:

- To interface AWS Batch::

  "batch:DescribeJobQueues"
  "batch:CancelJob"
  "batch:SubmitJob"
  "batch:ListJobs"
  "batch:DescribeComputeEnvironments"
  "batch:TerminateJob"
  "batch:DescribeJobs"
  "batch:RegisterJobDefinition"
  "batch:DescribeJobDefinitions"

- To be able to see the `EC2 <https://aws.amazon.com/ec2/>`_ instances::

  "ecs:DescribeTasks"
  "ec2:DescribeInstances"
  "ec2:DescribeInstanceTypes"
  "ec2:DescribeInstanceAttribute"
  "ecs:DescribeContainerInstances"
  "ec2:DescribeInstanceStatus"

- To pull container images stored in the `ECR <https://aws.amazon.com/ecr/>`_ repositories::

  "ecr:GetAuthorizationToken"
  "ecr:BatchCheckLayerAvailability"
  "ecr:GetDownloadUrlForLayer"
  "ecr:GetRepositoryPolicy"
  "ecr:DescribeRepositories"
  "ecr:ListImages"
  "ecr:DescribeImages"
  "ecr:BatchGetImage"
  "ecr:GetLifecyclePolicy"
  "ecr:GetLifecyclePolicyPreview"
  "ecr:ListTagsForResource"
  "ecr:DescribeImageScanFindings"

S3 policies
------------
Nextflow requires policies also to access `S3 buckets <https://aws.amazon.com/s3/>`_ in order to:

1. use the workdir
2. pull input data
3. publish results

Depending on the pipeline configuration, the above actions can be done all in a single bucket but, more likely, spread across multiple
buckets. Once the list of buckets used by the pipeline is identified, there are two alternative ways to give Nextflow access to these buckets:

1. grant access to all buckets by attaching the policy "s3:*" to the AIM identity. This works only if buckets do not set their own access policies (see point 2);
2. for a more fine grained control, assign to each bucket the following policy (replace the placeholders with the actual values)::

	{
    "Version": "2012-10-17",
    "Id": "<my policy id>",
    "Statement": [
        {
            "Sid": "<my statement id>",
            "Effect": "Allow",
            "Principal": {
                "AWS": "<ARN of the nextflow identity>"
            },
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:DeleteObject"
            ],
            "Resource": "arn:aws:s3:::<bucket name>/*"
        },
        {
            "Sid": "AllowSSLRequestsOnly",
            "Effect": "Deny",
            "Principal": "*",
            "Action": "s3:*",
            "Resource": [
                "arn:aws:s3:::<bucket name>",
                "arn:aws:s3:::<bucket name>/*"
            ],
            "Condition": {
                "Bool": {
                    "aws:SecureTransport": "false"
                }
            }
        }
    ]
	}

See the `bucket policy documentation <https://docs.aws.amazon.com/config/latest/developerguide/s3-bucket-policy.html>`_
for additional details.

.. _awscloud-batch:

AWS Batch
=========

.. note::
    Requires Nextflow version `0.26.0` or later.

`AWS Batch <https://aws.amazon.com/batch/>`_ is a managed computing service that allows the execution of containerised
workloads in the Amazon cloud infrastructure. It dynamically provisions the optimal quantity and type of compute
resources (e.g., CPU or memory optimized compute resources) based on the volume and specific resource requirements
of the jobs submitted.

Nextflow provides a built-in support for AWS Batch which allows the seamless deployment of a Nextflow pipeline
in the cloud offloading the process executions as Batch jobs.

.. _awscloud-batch-config:

AWS CLI
--------

Nextflow requires to access the `AWS command line tool <https://aws.amazon.com/cli/>`_ (``aws``) from the container in
which the job runs in order to stage the required input files and to copy back the resulting output files in the S3 storage.

The ``aws`` tool can be made available in the container in two ways:

1 - installed in the Docker image(s) used during the pipeline execution

2 - installed in a custom `AMI (Amazon Machine Image) <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html>`_ to
use in place of the default AMI when configuring AWS Batch (see next section).

The latter approach is preferred because it allows the use of existing Docker images without the need to add
the AWS CLI tool to them.

See the sections below to learn how to create a custom AMI and install the AWS CLI tool to it.

Get started
-------------

1 - In the AWS Console, create a `Compute environment <http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html>`_ (CE) in your AWS Batch Service.
    1.1 - if are using a custom AMI (see following sections), the AMI ID must be specified in the CE configuration
    1.2 - make sure to select an AMI (either custom or existing) with Docker installed (see following sections)
    1.3 - make sure the policy ``AmazonS3FullAccess`` (granting access to S3 buckets) is attached to the instance role configured for the CE
    1.4 - if you plan to use Docker images from Amazon ECS container, make sure the ``AmazonEC2ContainerServiceforEC2Role`` policy is also attached to the instance role

2 - In the AWS Console, create (at least) one `Job Queue <https://docs.aws.amazon.com/batch/latest/userguide/job_queues.html>`_
and bind it to the Compute environment

3 - In the AWS Console, create an S3 storage's bucket for the bucket-dir (see below) and others for the input data and
results, if/as needed

4 - Make sure your pipeline processes specifies one or more Docker containers by using the :ref:`process-container` directive.

5 - Container images need to be published in a Docker registry such as `Docker Hub <https://hub.docker.com/>`_,
`Quay <https://quay.io/>`_ or `ECS Container Registry <https://aws.amazon.com/ecr/>`_ that can be reached by ECS Batch.

Configuration
-------------

When configuring your pipeline:

1 - import the `nf-amazon` plugin
2 - specify the AWS Batch :ref:`executor<awsbatch-executor>`
3 - specify one or more AWS Batch queues for the execution by using the :ref:`process-queue` directive
4 - specify the AWS job container properties by using the :ref:`process-containerOptions` directive.

An example ``nextflow.config`` file is shown below::

    plugins {
        id 'nf-amazon'
    }

    process {
        executor = 'awsbatch'
        queue = 'my-batch-queue'
        container = 'quay.io/biocontainers/salmon'
        containerOptions = '--shm-size 16000000 --ulimit nofile=1280:2560 --ulimit nproc=16:32'
    }

    aws {
        batch {
            // NOTE: this setting is only required if the AWS CLI tool is installed in a custom AMI
            cliPath = '/home/ec2-user/miniconda/bin/aws'
        }
        region = 'us-east-1'
    }

Different queues bound to the same or different Compute environments can be configured according to each process' requirements.

Container Options
=================

As of version ``21.12.1-edge``, the use of the Nextflow :ref:`process-containerOptions` directive is supported to fine control
the properties of the container execution associated with each Batch job.

Not all the standard container options are supported by AWS Batch. These are the options accepted ::


    -e, --env string
        Set environment variables (format: <name> or <name>=<value>)
    --init
        Run an init inside the container that forwards signals and reaps processes
    --memory-swap int
        The total amount of swap memory (in MiB) the container can use: '-1' to enable unlimited swap
    --memory-swappiness int
        Tune container memory swappiness (0 to 100) (default -1)
    --privileged
        Give extended privileges to the container
    --read-only
        Mount the container's root filesystem as read only
    --shm-size int
        Size (in MiB) of /dev/shm
    --tmpfs string
        Mount a tmpfs directory (format: <path>:<options>,size=<int>), size is in MiB
    -u, --user string
        Username or UID (format: <name|uid>[:<group|gid>])
    --ulimit string
        Ulimit options (format: <type>=<soft limit>[:<hard limit>])

Container options must be passed in their long from for "--option value" or short form "-o value", if available.

Few examples::

  containerOptions '--tmpfs /run:rw,noexec,nosuid,size=128 --tmpfs /app:ro,size=64'

  containerOptions '-e MYVAR1 --env MYVAR2=foo2 --env MYVAR3=foo3 --memory-swap 3240000 --memory-swappiness 20 --shm-size 16000000'

  containerOptions '--ulimit nofile=1280:2560 --ulimit nproc=16:32 --privileged'


Check the `AWS doc <https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html>`_ for further details.

Custom AMI
==========
There are several reasons why you might need to create your own `AMI (Amazon Machine Image) <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html>`_
to use in your Compute environments. Typically:

1 - you do not want to modify your existing Docker images and prefer to install the CLI tool on the hosting environment

2 - the existing AMI (selected from the marketplace) does not have Docker installed

3 - you need to attach a larger storage to your EC2 instance (the default ECS instance AMI has only a 30G storage
volume which may not be enough for most data analysis pipelines)

4 - you need to install additional software, not available in the Docker image used to execute the job

Create your custom AMI
----------------------
In the EC2 Dashboard, click the `Launch Instance` button, then choose `AWS Marketplace` in the left pane and enter
`ECS` in the search box. In result list select `Amazon ECS-Optimized Amazon Linux 2 AMI`, then continue as usual to
configure and launch the instance.

.. note:: The selected instance has a bootstrap volume of 8GB and a second EBS volume 30G for computation which is
  hardly enough for real world genomic workloads. Make sure to specify an amount of storage in the second volume
  large enough for the needs of your pipeline execution.

When the instance is running, SSH into it (or connect with the Session Manager service), install the AWS CLI tool
or any other tool that may be required (see next sections).

Once done that, create a new AMI by using the *Create Image* option in the EC2 Dashboard or the AWS command line tool.

The new AMI ID needs to be specified when creating the Batch Compute Environment.

.. warning:: Any installation must be completed on the EC2 instance BEFORE creating the AMI.

.. _aws-cli:

AWS CLI installation
--------------------

.. warning:: The `AWS CLI tool <https://aws.amazon.com/cli>`_ must to be installed in your custom AMI
  by using a self-contained package manager such as `Conda <https://conda.io>`_.

The reason is that when the AWS CLI tool executes using Conda it will use the version of python supplied by Conda.
If you don't use Conda and install the AWS CLI using something like `pip <https://pypi.org/project/pip/>`_ the ``aws``
command will attempt to run using the version of python found in the running container which won't be able to find
the necessary dependencies.

The following snippet shows how to install AWS CLI with `Miniconda <https://conda.io/miniconda.html>`_ in the home folder::

    cd $HOME
    sudo yum install -y bzip2 wget
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
    $HOME/miniconda/bin/conda install -c conda-forge -y awscli
    rm Miniconda3-latest-Linux-x86_64.sh

When complete, verify that the AWS CLI package works correctly::

    $ ./miniconda/bin/aws --version
    aws-cli/1.19.79 Python/3.8.5 Linux/4.14.231-173.361.amzn2.x86_64 botocore/1.20.79


.. note:: The ``aws`` tool will be placed in a directory named ``bin`` in the main installation folder.
  Modifying this directory structure, after the installation, will cause the tool to not work properly.

To configure Nextflow to use this installation, specify the ``cliPath`` parameter in the :ref:`AWS Batch<config-aws-batch>`
configuration as shown below::

    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'

Replace the path above with the one matching the location where ``aws`` tool is installed in your AMI.

.. warning:: The grandparent directory of the ``aws`` tool will be mounted into the container at the same path as the host,
  e.g. ``/home/ec2-user/miniconda``, which will shadow existing files in the container.
  Ensure you use a path that is not already present in the container.

.. note:: Using a version of Nextflow prior 19.07.x the config setting `executor.awscli` should be used
  instead of `aws.batch.cliPath`.

Docker installation
---------------------------------------
Docker is required by Nextflow to execute tasks on AWS Batch. `Amazon ECS-Optimized Amazon Linux 2 AMI` has Docker installed,
however if you create your AMI starting from a different AMI that does not have Docker installed, you need to do it manually.

The following snippet shows how to install Docker on an Amazon EC2 instance::

    sudo yum update -y
    sudo amazon-linux-extras install docker
    sudo yum install docker
    sudo service docker start

Then, add the ``ec2-user`` to the docker group so you can execute Docker commands without using ``sudo``::

    sudo usermod -a -G docker ec2-user

You may have to reboot your instance to provide permissions for the ``ec2-user`` to access the Docker daemon. This has
to be done BEFORE creating the AMI from the current EC2 instance.

Amazon ECS container agent installation
---------------------------------------
The `ECS container agent <https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ECS_agent.html>`_ is a component
of Amazon Elastic Container Service (Amazon ECS) and is responsible for managing containers on behalf of Amazon ECS.
AWS Batch uses Amazon ECS to execute containerized jobs and therefore requires the agent to be installed on compute
resources within your Compute environments.

The ECS container agent is included in the `Amazon ECS-Optimized Amazon Linux 2 AMI`, but if you select a different AMI
you can also install it on any EC2 instance that supports the Amazon ECS specification.

To install the agent, follow these steps::

    sudo amazon-linux-extras disable docker
    sudo amazon-linux-extras install -y ecs
    sudo systemctl enable --now ecs

To test the installation::

    curl -s http://localhost:51678/v1/metadata | python -mjson.tool (test)

.. note:: The ``AmazonEC2ContainerServiceforEC2Role`` policy must be attached to the instance role in order to be able to
    connect the EC2 instance created by the Compute Environment to the ECS container.

Jobs & Execution
================

Custom job definition
---------------------

Nextflow automatically creates the Batch `Job definitions <http://docs.aws.amazon.com/batch/latest/userguide/job_definitions.html>`_
needed to execute your pipeline processes. Therefore it's not required to define them before running your workflow.

However you may still need to specify a custom `Job Definition` to fine control the configuration settings
of a specific job e.g. to define custom mount paths or other Batch Job special settings.

To do that first create a *Job Definition* in the AWS Console (or with other means). Note the name of the *Job Definition*
you created. You can then associate a process execution with this *Job definition* by using the :ref:`process-container`
directive and specifing, in place of the container image name, the Job definition name prefixed by the
``job-definition://`` string, as shown below::

  process.container = 'job-definition://your-job-definition-name'


Pipeline execution
------------------

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads.

Pipeline input data can be stored either locally or in a `S3 <https://aws.amazon.com/s3/>`_ bucket.
The pipeline execution must specifies a AWS Storage bucket where jobs intermediate results are stored with the
``-bucket-dir`` command line options. For example::

  nextflow run my-pipeline -bucket-dir s3://my-bucket/some/path


.. warning::
  The bucket path should include at least a top level directory name e.g. use ``s3://my-bucket/work``
  not just ``s3://my-bucket``. 

Hybrid workloads
----------------

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment
of hybrid workloads in which some jobs are execute in the local computer or local computing cluster and
some jobs are offloaded to AWS Batch service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the AWS Batch :ref:`configuration <awscloud-batch-config>` only to a subset of processes in your workflow.
For example::


  aws {
      region = 'eu-west-1'
      batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
      }
  }

  process {
      withLabel: bigTask {
        executor = 'awsbatch'
        queue = 'my-batch-queue'
        container = 'my/image:tag'
    }
  }


The above configuration snippet will deploy the execution with AWS Batch only for processes annotated
with the :ref:`process-label` ``bigTask``, the remaining process with run in the local computer.

Volume mounts
-------------

User provided container volume mounts can be provided as shown below::

  aws {
    region = 'eu-west-1'
    batch {
        volumes = '/tmp'
    }
  }

Multiple volumes can be specified using a comma separated paths. The usual Docker volume mount syntax
can be used to specify complex volumes for which the container paths is different from the host paths
or to specify *read-only* option. For example::

  aws {
    region = 'eu-west-1'
    batch {
        volumes = ['/tmp', '/host/path:/mnt/path:ro']
    }
  }


The above snippet defines two volume mounts the jobs executed in your pipeline. The first mounting the
host path ``/tmp`` in the same path in the container and using *read-write* access mode. The second
mounts the path ``/host/path`` in the host environment to the ``/mnt/path`` in the container using the
*read-only* access mode.

.. note:: This feature requires Nextflow version 19.07.x or later.

Troubleshooting
---------------

**Problem**: The Pipeline execution terminates with an AWS error message similar to the one shown below::

    JobQueue <your queue> not found


Make sure you have defined a AWS region in the Nextflow configuration file and it matches the region
in which your Batch environment has been created.

**Problem**: A process execution fails reporting the following error message::

  Process <your task> terminated for an unknown reason -- Likely it has been terminated by the external system

This may happen when Batch is unable to execute the process script. A common cause of this problem is that the
Docker container image you have specified uses a non standard `entrypoint <https://docs.docker.com/engine/reference/builder/#entrypoint>`_
which does not allow the execution of the Bash launcher script required by Nextflow to run the job.

This may also happen if the AWS CLI doesn't run correctly.

Other places to check for error information:

- The ``.nextflow.log`` file.
- The Job execution log in the AWS Batch dashboard.
- The `CloudWatch <https://aws.amazon.com/cloudwatch/>`_ logs found in the ``/aws/batch/job`` log group.

**Problem**: A process execution is stalled in the ``RUNNABLE`` status and the pipeline output is similar to the one below::

    executor >  awsbatch (1)
    process > <your process> (1) [  0%] 0 of ....

It may happen that the pipeline execution hangs indefinitely because one of the jobs is held in the queue and never gets
executed. In AWS Console, the queue reports the job as ``RUNNABLE`` but it never moves from there.

There are multiple reasons why this can happen. They are mainly related to the Compute Environment workload/configuration,
the docker service or container configuration, network status, etc.

This `AWS page <https://aws.amazon.com/premiumsupport/knowledge-center/batch-job-stuck-runnable-status/>`_ provides several
resolutions and tips to investigate and work around the issue.

Advanced configuration
======================

Read :ref:`AWS Batch configuration<config-aws-batch>` section to learn more about advanced Batch configuration options.
.. _script-page:

******************
Nextflow scripting
******************


The Nextflow scripting language is an extension of the Groovy programming language.
Groovy is a powerful programming language for the Java virtual machine. The Nextflow
syntax has been specialized to ease the writing of computational pipelines in a declarative manner.

Nextflow can execute any piece of Groovy code or use any library for the JVM platform.

For a detailed description of the Groovy programming language, reference these links:

* `Groovy User Guide <http://groovy-lang.org/documentation.html>`_
* `Groovy Cheat sheet <http://www.cheat-sheets.org/saved-copy/rc015-groovy_online.pdf>`_
* `Groovy in Action <http://www.manning.com/koenig2/>`_


Below you can find a crash course in the most important language constructs used in the Nextflow scripting language.

.. warning:: Nextflow uses ``UTF-8`` as the default file character encoding for source and application files. Make sure
  to use the ``UTF-8`` encoding when editing Nextflow scripts with your favourite text editor.

Language basics
==================


Hello world
------------

To print something is as easy as using one of the ``print`` or ``println`` methods.
::

    println "Hello, World!"

The only difference between the two is that the ``println`` method implicitly appends a `new line` character
to the printed string.


Variables
----------

To define a variable, simply assign a value to it::

    x = 1
    println x

    x = new java.util.Date()
    println x

    x = -3.1499392
    println x

    x = false
    println x

    x = "Hi"
    println x


Lists
------

A List object can be defined by placing the list items in square brackets::

    myList = [1776, -1, 33, 99, 0, 928734928763]

You can access a given item in the list with square-bracket notation (indexes start at 0)::

    println myList[0]

In order to get the length of the list use the ``size`` method::

    println myList.size()


Learn more about lists:

* `Groovy Lists tutorial <http://groovy-lang.org/groovy-dev-kit.html#Collections-Lists>`_
* `Groovy List SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html>`_
* `Java List SDK <http://docs.oracle.com/javase/7/docs/api/java/util/List.html>`_


Maps
-----

Maps are used to store `associative arrays` or `dictionaries`. They are unordered collections of heterogeneous, named data::

    scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]


Note that each of the values stored in the map can be of a different type. ``Brett`` is an integer, ``Pete`` is a string,
and ``Andrew`` is a floating-point number.

We can access the values in a map in two main ways::

    println scores["Pete"]
    println scores.Pete


To add data to or modify a map, the syntax is similar to adding values to list::

    scores["Pete"] = 3
    scores["Cedric"] = 120


Learn more about maps:

* `Groovy Maps tutorial <http://groovy-lang.org/groovy-dev-kit.html#Collections-Maps>`_
* `Groovy Map SDK <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html>`_
* `Java Map SDK <http://docs.oracle.com/javase/7/docs/api/java/util/Map.html>`_


.. _script-multiple-assignment:

Multiple assignment
----------------------

An array or a list object can used to assign to multiple variables at once::

    (a, b, c) = [10, 20, 'foo']
    assert a == 10 && b == 20 && c == 'foo'

The three variables on the left of the assignment operator are initialized by the corresponding item in the list.

Read more about `Multiple assignment <http://www.groovy-lang.org/semantics.html#_multiple_assignment>`_ in the Groovy documentation.


Conditional Execution
----------------------

One of the most important features of any programming language is the ability to execute different code under
different conditions. The simplest way to do this is to use the ``if`` construct::

    x = Math.random()
    if( x < 0.5 ) {
        println "You lost."
    }
    else {
        println "You won!"
    }



Strings
-------

Strings can be defined by enclosing text in single or double quotes (``'`` or ``"`` characters)::

    println "he said 'cheese' once"
    println 'he said "cheese!" again'


Strings can be concatenated with ``+``::

    a = "world"
    print "hello " + a + "\n"


.. _string-interpolation:

String interpolation
--------------------

There is an important difference between single-quoted and double-quoted strings:
Double-quoted strings support variable interpolations, while single-quoted strings do not.

In practice, double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the ``$`` character,
or the value of any expression by using the ``${expression}`` syntax, similar to Bash/shell scripts::

    foxtype = 'quick'
    foxcolor = ['b', 'r', 'o', 'w', 'n']
    println "The $foxtype ${foxcolor.join()} fox"

    x = 'Hello'
    println '$x + $y'

This code prints::

    The quick brown fox
    $x + $y


Multi-line strings
-------------------

A block of text that span multiple lines can be defined by delimiting it with triple single or double quotes::

    text = """
        hello there James
        how are you today?
        """

.. note:: Like before, multi-line strings inside double quotes support variable interpolation, while
   single-quoted multi-line strings do not.


As in Bash/shell scripts, terminating a line in a multi-line string with a ``\`` character prevents a
a `new line` character from separating that line from the one that follows::

    myLongCmdline = """ blastp \
                    -in $input_query \
                    -out $output_file \
                    -db $blast_database \
                    -html
                    """

    result = myLongCmdline.execute().text

In the preceding example, ``blastp`` and its ``-in``, ``-out``, ``-db`` and ``-html`` switches and
their arguments are effectively a single line.

.. _implicit-variables:

Implicit variables
==================

Script implicit variables
-------------------------

The following variables are implicitly defined in the script global execution scope:

=============== ========================
Name            Description
=============== ========================
``baseDir``     The directory where the main workflow script is located (deprecated in favour of ``projectDir`` since ``20.04.0``).
``launchDir``   The directory where the workflow is run (requires version ``20.04.0`` or later).
``moduleDir``   The directory where a module script is located for DSL2 modules or the same as ``projectDir`` for a non-module script (requires version ``20.04.0`` or later).
``nextflow``    Dictionary like object representing nextflow runtime information (see :ref:`metadata-nextflow`).
``params``      Dictionary like object holding workflow parameters specifing in the config file or as command line options.
``projectDir``  The directory where the main script is located (requires version ``20.04.0`` or later).
``workDir``     The directory where tasks temporary files are created.
``workflow``    Dictionary like object representing workflow runtime information (see :ref:`metadata-workflow`).
=============== ========================


Configuration implicit variables
--------------------------------

The following variables are implicitly defined in the Nextflow configuration file:

=============== ========================
Name            Description
=============== ========================
``baseDir``     The directory where the main workflow script is located (deprecated in favour of ``projectDir`` since ``20.04.0``).
``launchDir``   The directory where the workflow is run (requires version ``20.04.0`` or later).
``projectDir``  The directory where the main script is located (requires version ``20.04.0`` or later).
=============== ========================


Process implicit variables
--------------------------

In the process definition scope it's available the ``task`` implicit variable which allow accessing
the current task configuration directives. For examples::

    process foo {
      script:
      """
      some_tool --cpus $task.cpus --mem $task.memory
      """
    }


In the above snippet the ``task.cpus`` report the value for the :ref:`cpus directive<process-cpus>` and
the ``task.memory`` the current value for :ref:`memory directive<process-memory>` depending on the actual
setting given in the workflow configuration file.

See :ref:`Process directives <process-directives>` for details.


.. _script-closure:

Closures
=========

Briefly, a closure is a block of code that can be passed as an argument to a function.
Thus, you can define a chunk of code and then pass it around as if it were a string or an integer.

More formally, you can create functions that are defined as `first class objects`.

::

    square = { it * it }


The curly brackets around the expression ``it * it`` tells the script interpreter to treat this expression as code.
The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled the function object is assigned to the variable ``square`` as any other variable assignments shown previously.
Now we can do something like this::

    println square(9)

and get the value 81.


This is not very interesting until we find that we can pass the function ``square`` as an argument to other functions or methods.
Some built-in functions take a function like this as an argument. One example is the ``collect`` method on lists::

    [ 1, 2, 3, 4 ].collect(square)


This expression says: Create an array with the values 1, 2, 3 and 4, then call its ``collect`` method, passing in the
closure we defined above. The ``collect`` method runs through each item in the array, calls the closure on the item,
then puts the result in a new array, resulting in::

    [ 1, 4, 9, 16 ]


For more methods that you can call with closures as arguments, see the `Groovy GDK documentation <http://docs.groovy-lang.org/latest/html/groovy-jdk/>`_.


By default, closures take a single parameter called ``it``, but you can also create closures with multiple, custom-named parameters.
For example, the method ``Map.each()`` can take a closure with two arguments, to which it binds the `key` and the associated `value`
for each key-value pair in the ``Map``. Here, we use the obvious variable names ``key`` and ``value`` in our closure::


    printMapClosure = { key, value ->
        println "$key = $value"
    }

    [ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ].each(printMapClosure)


Prints::


    Yue = Wu
    Mark = Williams
    Sudha = Kumari

A closure has two other important features. First, it can access variables in the scope where it is defined,
so that it can `interact` with them.

Second, a closure can be defined in an `anonymous` manner, meaning that it is not given a name,
and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment::

    myMap = ["China": 1 , "India" : 2, "USA" : 3]

    result = 0
    myMap.keySet().each( { result+= myMap[it] } )

    println result


Learn more about closures in the `Groovy documentation <http://groovy-lang.org/closures.html>`_

.. _script-regexp:

Regular expressions
====================

Regular expressions are the Swiss Army knife of text processing. They provide the programmer with the ability to match
and extract patterns from strings.

Regular expressions are available via the ``~/pattern/`` syntax and the ``=~`` and ``==~`` operators.

Use ``=~`` to check whether a given pattern occurs anywhere in a string::

    assert 'foo' =~ /foo/       // return TRUE
    assert 'foobar' =~ /foo/    // return TRUE


Use ``==~`` to check whether a string matches a given regular expression pattern exactly.
::

    assert 'foo' ==~ /foo/       // return TRUE
    assert 'foobar' ==~ /foo/    // return FALSE


It is worth noting that the ``~`` operator creates a Java ``Pattern`` object from the given string,
while the ``=~`` operator creates a Java ``Matcher`` object.
::

    x = ~/abc/
    println x.class
    // prints java.util.regex.Pattern

    y = 'some string' =~ /abc/
    println y.class
    // prints java.util.regex.Matcher


Regular expression support is imported from Java. Java's regular expression language and API is documented in the
`Pattern Java documentation <http://download.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html>`_.

You may also be interested in this post: `Groovy: Don't Fear the RegExp <https://web.archive.org/web/20170621185113/http://www.naleid.com/blog/2008/05/19/dont-fear-the-regexp>`_.


String replacement
--------------------

To replace pattern occurrences in a given string, use the ``replaceFirst`` and ``replaceAll`` methods::

     x = "colour".replaceFirst(/ou/, "o")
     println x
     // prints: color

     y = "cheesecheese".replaceAll(/cheese/, "nice")
     println y
     // prints: nicenice



Capturing groups
----------------

You can match a pattern that includes groups.  First create a matcher object with the ``=~`` operator.
Then, you can index the matcher object to find the matches: ``matcher[0]`` returns a list representing the first match
of the regular expression in the string. The first list element is the string that matches the entire regular expression, and
the remaining elements are the strings that match each group.

Here's how it works::

    programVersion = '2.7.3-beta'
    m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/

    assert m[0] ==  ['2.7.3-beta', '2', '7', '3', 'beta']
    assert m[0][1] == '2'
    assert m[0][2] == '7'
    assert m[0][3] == '3'
    assert m[0][4] == 'beta'



Applying some syntactic sugar, you can do the same in just one line of code::

    programVersion = '2.7.3-beta'
    (full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

    println full    // 2.7.3-beta
    println major   // 2
    println minor   // 7
    println patch   // 3
    println flavor  // beta


Removing part of a string
-------------------------

You can remove part of a ``String`` value using a regular expression pattern. The first match found is
replaced with an empty String::

    // define the regexp pattern
    wordStartsWithGr = ~/(?i)\s+Gr\w+/

    // apply and verify the result
    ('Hello Groovy world!' - wordStartsWithGr) == 'Hello world!'
    ('Hi Grails users' - wordStartsWithGr) == 'Hi users'


Remove the first 5-character word from a string::

    assert ('Remove first match of 5 letter word' - ~/\b\w{5}\b/) == 'Remove  match of 5 letter word'


Remove the first number with its trailing whitespace from a string::

    assert ('Line contains 20 characters' - ~/\d+\s+/) == 'Line contains characters'



.. _script-file-io:

Files and I/O
==============

To access and work with files, use the ``file`` method, which returns a file system object
given a file path string::

  myFile = file('some/path/to/my_file.file')


The ``file`` method can reference either `files` or `directories`, depending on what the string path refers to in the
file system.

When using the wildcard characters ``*``, ``?``, ``[]`` and ``{}``, the argument is interpreted as a `glob`_ path matcher
and the ``file`` method returns a list object holding the paths of files whose names match the specified pattern, or an
empty list if no match is found::

  listOfFiles = file('some/path/*.fa')

.. note:: Two asterisks (``**``) in a glob pattern works like ``*`` but matches any number of directory components in a
          file system path.

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden
files in the result list, add the optional parameter ``hidden``::

  listWithHidden = file('some/path/*.fa', hidden: true)

Here are ``file``'s available options:

=============== ===================
Name            Description
=============== ===================
glob            When ``true`` interprets characters ``*``, ``?``, ``[]`` and ``{}`` as glob wildcards, otherwise handles them as normal characters (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` follows symbolic links during directory tree traversal, otherwise treats them as files (default: ``true``)
checkIfExists   When ``true`` throws an exception of the specified path do not exist in the file system (default: ``false``)
=============== ===================


.. tip:: If you are a Java geek you will be interested to know that the ``file`` method returns a
  `Path <http://docs.oracle.com/javase/8/docs/api/java/nio/file/Path.html>`_ object, which allows
  you to use the usual methods you would in a Java program.

See also: :ref:`Channel.fromPath <channel-path>`.

.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

Basic read/write
------------------

Given a file variable, declared using the ``file`` method as shown in the previous example, reading a file
is as easy as getting the value of the file's ``text`` property, which returns the file content
as a string value::

  print myFile.text


Similarly, you can save a string value to a file by simply assigning it to the file's ``text`` property::

  myFile.text = 'Hello world!'


.. note:: Existing file content is overwritten by the assignment operation, which also implicitly creates
          files that do not exist.

In order to append a string value to a file without erasing existing content, you can use the ``append`` method::

  myFile.append('Add this line\n')

Or use the `left shift` operator, a more idiomatic way to append text content to a file::

  myFile << 'Add a line more\n'


Binary data can managed in the same way, just using the file property ``bytes`` instead of ``text``. Thus, the following
example reads the file and returns its content as a byte array::

  binaryContent = myFile.bytes

Or you can save a byte array data buffer to a file, by simply writing::

  myFile.bytes = binaryBuffer


.. warning:: The above methods read and write ALL the file content at once, in a single variable or buffer. For this
  reason they are not suggested when dealing with big files, which require a more memory efficient approach, for example
  reading a file line by line or by using a fixed size buffer.


Read a file line by line
--------------------------

In order to read a text file line by line you can use the method ``readLines()`` provided by the file object, which
returns the file content as a list of strings::

    myFile = file('some/my_file.txt')
    allLines  = myFile.readLines()
    for( line : allLines ) {
        println line
    }


This can also be written in a more idiomatic syntax::

    file('some/my_file.txt')
        .readLines()
        .each { println it }


.. note:: The method ``readLines()`` reads all the file content at once and returns a list containing all the lines. For
  this reason, do not use it to read big files.


To process a big file, use the method ``eachLine``, which reads only a single line at a time into memory::

    count = 0
    myFile.eachLine {  str ->
            println "line ${count++}: $str"
        }



Advanced file reading operations
-----------------------------------

The classes ``Reader`` and ``InputStream`` provide fine control for reading text and binary files, respectively._


The method ``newReader`` creates a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
for the given file that allows you to read the content as single characters, lines or arrays of characters::

    myReader = myFile.newReader()
    String line
    while( line = myReader.readLine() ) {
        println line
    }
    myReader.close()


The method ``withReader`` works similarly, but automatically calls the ``close`` method for you when you have finished
processing the file. So, the previous example can be written more simply as::

    myFile.withReader {
        String line
        while( line = it.readLine() ) {
            println line
        }
    }

The methods ``newInputStream`` and ``withInputStream`` work similarly. The main difference is that they create an
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object useful for writing binary
data.

Here are the most important methods for reading from files:

=============== ==============
Name            Description
=============== ==============
getText         Returns the file content as a string value
getBytes        Returns the file content as byte array
readLines       Reads the file line by line and returns the content as a list of strings
eachLine        Iterates over the file line by line, applying the specified :ref:`closure <script-closure>`
eachByte        Iterates over the file byte by byte, applying the specified :ref:`closure <script-closure>`
withReader      Opens a file for reading and lets you access it with a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object
withInputStream Opens a file for reading and lets you access it with an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object
newReader       Returns a `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ object to read a text file
newInputStream  Returns an `InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ object to read a binary file
=============== ==============


Read the Java documentation for `Reader <http://docs.oracle.com/javase/7/docs/api/java/io/Reader.html>`_ and
`InputStream <http://docs.oracle.com/javase/7/docs/api/java/io/InputStream.html>`_ classes to learn more about
methods available for reading data from files.


Advanced file writing operations
----------------------------------

The ``Writer`` and ``OutputStream`` classes provide fine control for writing text and binary files,
respectively, including low-level operations for single characters or bytes, and support for big files.

For example, given two file objects ``sourceFile`` and ``targetFile``, the following code copies the
first file's content into the second file, replacing all ``U`` characters with ``X``::

    sourceFile.withReader { source ->
        targetFile.withWriter { target ->
            String line
            while( line=source.readLine() ) {
                target << line.replaceAll('U','X')
            }
        }
    }


Here are the most important methods for writing to files:

=================== ==============
Name                Description
=================== ==============
setText             Writes a string value to a file
setBytes            Writes a byte array to a file
write               Writes a string to a file, replacing any existing content
append              Appends a string value to a file without replacing existing content
newWriter           Creates a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object that allows you to save text data to a file
newPrintWriter      Creates a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object that allows you to write formatted text to a file
newOutputStream     Creates an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object that allows you to write binary data to a file
withWriter          Applies the specified closure to a `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_ object, closing it when finished
withPrintWriter     Applies the specified closure to a `PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ object, closing it when finished
withOutputStream    Applies the specified closure to an `OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ object, closing it when finished
=================== ==============

Read the Java documentation for the `Writer <http://docs.oracle.com/javase/7/docs/api/java/io/Writer.html>`_,
`PrintWriter <http://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html>`_ and
`OutputStream <http://docs.oracle.com/javase/7/docs/api/java/io/OutputStream.html>`_ classes to learn more about
methods available for writing data to files.


List directory content
----------------------

Let's assume that you need to walk through a directory of your choice. You can define the ``myDir`` variable
that points to it::

    myDir = file('any/path')

The simplest way to get a directory list is by using the methods ``list`` or ``listFiles``,
which return a collection of first-level elements (files and directories) of a directory::

    allFiles = myDir.list()
    for( def file : allFiles ) {
        println file
    }

.. note:: The only difference between ``list`` and ``listFiles`` is that the former returns a list of strings, and the latter a
   list of file objects that allow you to access file metadata, e.g. size, last modified time, etc.


The ``eachFile`` method allows you to iterate through the first-level elements only
(just like ``listFiles``). As with other `each-` methods, ``eachFiles`` takes a closure as a parameter::

    myDir.eachFile { item ->
        if( item.isFile() ) {
            println "${item.getName()} - size: ${item.size()}"
        }
        else if( item.isDirectory() ) {
            println "${item.getName()} - DIR"
        }
    }


Several variants of the above method are available. See the table below for a complete list.

=================== ==================
Name                Description
=================== ==================
eachFile            Iterates through first-level elements (files and directories). `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFile(groovy.io.FileType,%20groovy.lang.Closure)>`_
eachDir             Iterates through first-level directories only. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDir(groovy.lang.Closure)>`_
eachFileMatch       Iterates through files and dirs whose names match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachDirMatch        Iterates through directories whose names match the given filter. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirMatch(java.lang.Object,%20groovy.lang.Closure)>`_
eachFileRecurse     Iterates through directory elements depth-first. `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachFileRecurse(groovy.lang.Closure)>`_
eachDirRecurse      Iterates through directories depth-first (regular files are ignored). `Read more <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/io/File.html#eachDirRecurse(groovy.lang.Closure)>`_
=================== ==================


See also: Channel :ref:`channel-path` method.


Create directories
------------------

Given a file variable representing a nonexistent directory, like the following::

  myDir = file('any/path')

the method ``mkdir`` creates a directory at the given path, returning ``true`` if the directory is created
successfully, and ``false`` otherwise::

   result = myDir.mkdir()
   println result ? "OK" : "Cannot create directory: $myDir"

.. note:: If the parent directories do not exist, the above method will fail and return ``false``.

The method ``mkdirs`` creates the directory named by the file object, including any nonexistent parent directories::

    myDir.mkdirs()


Create links
------------

Given a file, the method ``mklink`` creates a *file system link* for that file using the path specified as a parameter::

  myFile = file('/some/path/file.txt')
  myFile.mklink('/user/name/link-to-file.txt')


Table of optional parameters:

==================  ================
Name                Description
==================  ================
hard                When ``true`` creates a *hard* link, otherwise creates a *soft* (aka *symbolic*) link. (default: ``false``)
overwrite           When ``true`` overwrites any existing file with the same name, otherwise throws a `FileAlreadyExistsException <http://docs.oracle.com/javase/8/docs/api/java/nio/file/FileAlreadyExistsException.html>`_ (default: ``false``)
==================  ================


Copy files
----------

The method ``copyTo`` copies a file into a new file or into a directory, or copies a directory to a new
directory::

  myFile.copyTo('new_name.txt')


.. note:: If the target file already exists, it will be replaced by the new one. Note also that, if the target is
  a directory, the source file will be copied into that directory, maintaining the file's original name.


When the source file is a directory, all its content is copied to the target directory::

  myDir = file('/some/path')
  myDir.copyTo('/some/new/path')


  If the target path does not exist, it will be created automatically.

.. tip:: The ``copyTo`` method mimics the semantics of the Linux command ``cp -r <source> <target>``, with the
         following caveat: While Linux tools often treat paths ending with a slash (e.g. ``/some/path/name/``)
         as directories, and those not (e.g. ``/some/path/name``) as regular files, Nextflow (due to its use of
         the Java files API) views both these paths as the same file system object. If the path exists, it is
         handled according to its actual type (i.e. as a regular file or as a directory). If the path does not
         exist, it is treated as a regular file, with any missing parent directories created automatically.



Move files
----------

You can move a file by using the method ``moveTo``::

  myFile = file('/some/path/file.txt')
  myFile.moveTo('/another/path/new_file.txt')


.. note:: When a file with the same name as the target already exists, it will be replaced by the source. Note
          also that, when the target is a directory, the file will be moved to (or within) that directory,
          maintaining the file's original name.

When the source is a directory, all the directory content is moved to the target directory::

  myDir = file('/any/dir_a')
  myDir.moveTo('/any/dir_b')


Please note that the result of the above example depends on the existence of the target directory. If the target
directory exists, the source is moved into the target directory, resulting in the path::

  /any/dir_b/dir_a

If the target directory does not exist, the source is just renamed to the target name, resulting in the path::

  /any/dir_b


.. tip:: The ``moveTo`` method mimics the semantics of the Linux command ``mv <source> <target>``, with the
         same caveat as that given for ``copyTo``, above.


Rename files
------------

You can rename a file or directory by simply using the ``renameTo`` file method::

  myFile = file('my_file.txt')
  myFile.renameTo('new_file_name.txt')


Delete files
------------

The file method ``delete`` deletes the file or directory at the given path, returning ``true`` if the
operation succeeds, and ``false`` otherwise::

  myFile = file('some/file.txt')
  result = myFile.delete()
  println result ? "OK" : "Cannot delete: $myFile"


.. note:: This method deletes a directory ONLY if it does not contain any files or sub-directories. To
          delete a directory and ALL its content (i.e. removing all the files and sub-directories it may
          contain), use the method ``deleteDir``.


Check file attributes
---------------------

The following methods can be used on a file variable created by using the ``file`` method:

==================  ================
Name                Description
==================  ================
getName             Gets the file name e.g. ``/some/path/file.txt`` -> ``file.txt``
getBaseName         Gets the file name without its extension e.g. ``/some/path/file.tar.gz`` -> ``file.tar``
getSimpleName       Gets the file name without any extension e.g. ``/some/path/file.tar.gz`` -> ``file``
getExtension        Gets the file extension e.g. ``/some/path/file.txt`` -> ``txt``
getParent           Gets the file parent path e.g. ``/some/path/file.txt`` -> ``/some/path``
size                Gets the file size in bytes
exists              Returns ``true`` if the file exists, or ``false`` otherwise
isEmpty             Returns ``true`` if the file is zero length or does not exist, ``false`` otherwise
isFile              Returns ``true`` if it is a regular file e.g. not a directory
isDirectory         Returns ``true`` if the file is a directory
isHidden            Returns ``true`` if the file is hidden
lastModified        Returns the file last modified timestamp i.e. a long as Linux epoch time
==================  ================


For example, the following line prints a file name and size::

  println "File ${myFile.getName() size: ${myFile.size()}"


.. tip:: The invocation of any method name starting with the ``get`` prefix can be shortcut
    omitting the `get` prefix and ending ``()`` parentheses. Therefore writing ``myFile.getName()``
    is exactly the same of ``myFile.name`` and ``myFile.getBaseName()`` is the same of ``myFile.baseName``
    and so on.


Get and modify file permissions
-------------------------------

Given a file variable representing a file (or directory), the method ``getPermissions`` returns a
9-character string representing the file's permissions using the
`Linux symbolic notation <http://en.wikipedia.org/wiki/File_system_permissions#Symbolic_notation>`_
e.g. ``rw-rw-r--``::

    permissions = myFile.getPermissions()


Similarly, the method ``setPermissions`` sets the file's permissions using the same notation::

    myFile.setPermissions('rwxr-xr-x')


A second version of the ``setPermissions`` method sets a file's permissions given three digits representing,
respectively, the `owner`, `group` and `other` permissions::

    myFile.setPermissions(7,5,5)


Learn more about `File permissions numeric notation <http://en.wikipedia.org/wiki/File_system_permissions#Numeric_notation>`_.

HTTP/FTP files
--------------

Nextflow provides transparent integration of HTTP/S and FTP protocols for handling remote resources
as local file system objects. Simply specify the resource URL as the argument of the `file` object::

    pdb = file('http://files.rcsb.org/header/5FID.pdb')

Then, you can access it as a local file as described in the previous sections::

    println pdb.text

The above one-liner prints the content of the remote PDB file. Previous sections provide code examples
showing how to stream or copy the content of files.

.. note:: Write and list operations are not supported for HTTP/S and FTP files.


Counting records
----------------

countLines
^^^^^^^^^^

The ``countLines`` methods counts the lines in a text files.
::

    def sample = file('/data/sample.txt')
    println sample.countLines()


Files whose name ends with the ``.gz`` suffix are expected to be GZIP compressed and
automatically uncompressed.

countFasta
^^^^^^^^^^

The ``countFasta`` method counts the number of records in `FASTA <https://en.wikipedia.org/wiki/FASTA_format>`_
formatted file.
::

    def sample = file('/data/sample.fasta')
    println sample.countFasta()

Files whose name ends with the ``.gz`` suffix are expected to be GZIP compressed and
automatically uncompressed.

countFastq
^^^^^^^^^^

The ``countFastq`` method counts the number of records in a `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_
formatted file.
::

    def sample = file('/data/sample.fastq')
    println sample.countFastq()

Files whose name ends with the ``.gz`` suffix are expected to be GZIP compressed and
automatically uncompressed.


.. _perfanalysis-page:

***********************
Tracing & visualisation
***********************


.. _execution-log:

Execution log
=============

The ``nextflow log`` command shows information about executed pipelines in the current folder::

  nextflow log <run name> [options]

.. note:: Both the :ref:`execution report <execution-report>` and the :ref:`trace report <trace-report>` must be specified when the pipeline is first called. By contrast, the ``log`` option is useful after a pipeline has already run and is available for every executed pipeline.

By default, ``log`` prints the list of executed pipelines::

  $ nextflow log
  TIMESTAMP            RUN NAME         SESSION ID                            COMMAND
  2016-08-01 11:44:51  grave_poincare   18cbe2d3-d1b7-4030-8df4-ae6c42abaa9c  nextflow run hello
  2016-08-01 11:44:55  small_goldstine  18cbe2d3-d1b7-4030-8df4-ae6c42abaa9c  nextflow run hello -resume
  2016-08-01 11:45:09  goofy_kilby      0a1f1589-bd0e-4cfc-b688-34a03810735e  nextflow run rnatoy -with-docker

Specifying a run name or session id prints tasks executed by that pipeline run::

  $ nextflow log goofy_kilby
  /Users/../work/0b/be0d1c4b6fd6c778d509caa3565b64
  /Users/../work/ec/3100e79e21c28a12ec2204304c1081
  /Users/../work/7d/eb4d4471d04cec3c69523aab599fd4
  /Users/../work/8f/d5a26b17b40374d37338ccfe967a30
  /Users/../work/94/dfdfb63d5816c9c65889ae34511b32


Customizing fields
------------------

By default, only the task execution paths are printed. A custom list of fields to print can be provided via the ``-f`` (``-fields``) option. For example::

  $ nextflow log goofy_kilby -f hash,name,exit,status
  0b/be0d1c  buildIndex (ggal_1_48850000_49020000.Ggal71.500bpflank)  0  COMPLETED
  ec/3100e7  mapping (ggal_gut)                                       0  COMPLETED
  7d/eb4d44  mapping (ggal_liver)                                     0  COMPLETED
  8f/d5a26b  makeTranscript (ggal_liver)                              0  COMPLETED
  94/dfdfb6  makeTranscript (ggal_gut)                                0  COMPLETED

The fields accepted by the ``-f`` options are the ones in the :ref:`trace report<trace-fields>`, as well as: script, stdout, stderr, env. List available fields using the ``-l`` (``-list-fields``) option.

The ``script`` field is useful for examining script commands run in each task::

  $ nextflow log goofy_kilby -f name,status,script
  align_genome      COMPLETED
     bowtie --index /data/genome input.fastq > output
  ...

Templates
---------

The ``-t`` option allows a template (string or file) to be specified. This makes it possible to create complex custom report in any text based format.  For example you could save this markdown snippet to a file::

  ## $name

  script:

      $script

  exist status: $exit
  task status: $status
  task folder: $folder

Then, the following command will output a markdown file containing the script, exit status and folder of all executed tasks::

  nextflow log goofy_kilby -t my-template.md > execution-report.md


Filtering
---------

The ``filter`` option makes it possible to select which entries to include in the log report. Any valid groovy boolean expression on the log fields can be used to define the filter condition. For example::

  nextflow log goofy_kilby -filter 'name =~ /foo.*/ && status == "FAILED"'

.. warning:: The ``log`` command replaces the deprecated ``history`` command.


.. _execution-report:

Execution report
================

Nextflow can create an HTML execution report: a single document which includes many useful metrics
about a workflow execution. The report is organised in the three main sections: `Summary`, `Resources` and `Tasks`
(see below for details).

To enable the creation of this report add the ``-with-report`` command line option when launching the pipeline
execution. For example::

  nextflow run <pipeline name> -with-report [file name]

The report file name can be specified as an optional parameter following the report option.


Summary
-------

The `Summary` section reports the execution status, the launch command, overall execution time and some
other workflow metadata. You can see an example below:

.. image:: images/report-summary-min.png


Resource Usage
---------------

The `Resources` sections plots the distributions of resource usages for each workflow process
using the interactive `plotly.js  <https://plot.ly/javascript/>`_ plotting library.

Plots are shown for CPU, memory, job duration and disk I/O. They have two (or three) tabs with the raw values and a percentage representation showing what proportion of the requested resources
were used. These plots are very helpful to check that job pipeline requests are efficient.

.. image:: images/report-resource-cpu.png

Learn more about how resource usage are computed in the :ref:`Metrics documentation <metrics-page>`.

Tasks
-----

Finally the `Tasks` section lists all executed tasks reporting for each of them, the status, the actual command script
and many other runtime metrics. You can see an example below:

.. image:: images/report-tasks-min.png


.. note:: Nextflow collect these metrics running a background process for each job in the target environment.
  Make sure the following tools are available ``awk``, ``date``, ``grep``, ``egrep``, ``ps``, ``sed``, ``tail``, ``tee`` in the
  system where the jobs are executed. Moreover some of these metrics are not reported when using a Mac OSX system. See the note
  message about that in the `Trace report`_ below.

.. warning:: A common problem when using a third party container image is that it does not ship one or more of the
  above utilities resulting in an empty execution report.

Please read :ref:`Report scope <config-report>` section to learn more about the execution report configuration details.

.. _trace-report:

Trace report
============

Nextflow creates an execution tracing file that contains some useful information about each process executed in your pipeline
script, including: submission time, start time, completion time, cpu and memory used.

In order to create the execution trace file add the ``-with-trace`` command line option when launching the pipeline execution.
For example::

  nextflow run <pipeline name> -with-trace

It will create a file named ``trace.txt`` in the current directory. The content looks like the above example:

======= ========= ========= =============== =========== ======== ======================= =========== =========== ======= =========== =========== =========== ===========
task_id hash      native_id   name          status      exit     submit                  duration    walltime    %cpu    rss         vmem        rchar       wchar
======= ========= ========= =============== =========== ======== ======================= =========== =========== ======= =========== =========== =========== ===========
19      45/ab752a 2032      blast (1)       COMPLETED   0        2014-10-23 16:33:16.288 1m          5s          0.0%    29.8 MB     354 MB      33.3 MB     0
20      72/db873d 2033      blast (2)       COMPLETED   0        2014-10-23 16:34:17.211 30s         10s         35.7%   152.8 MB    428.1 MB    192.7 MB    1 MB
21      53/d13188 2034      blast (3)       COMPLETED   0        2014-10-23 16:34:17.518 29s         20s         4.5%    289.5 MB    381.6 MB    33.3 MB     0
22      26/f65116 2035      blast (4)       COMPLETED   0        2014-10-23 16:34:18.459 30s         9s          6.0%    122.8 MB    353.4 MB    33.3 MB     0
23      88/bc00e4 2036      blast (5)       COMPLETED   0        2014-10-23 16:34:18.507 30s         19s         5.0%    195 MB      395.8 MB    65.3 MB     121 KB
24      74/2556e9 2037      blast (6)       COMPLETED   0        2014-10-23 16:34:18.553 30s         12s         43.6%   140.7 MB    432.2 MB    192.7 MB    182.7 MB
28      b4/0f9613 2041      exonerate (1)   COMPLETED   0        2014-10-23 16:38:19.657 1m 30s      1m 11s      94.3%   611.6 MB    693.8 MB    961.2 GB    6.1 GB
32      af/7f2f57 2044      exonerate (4)   COMPLETED   0        2014-10-23 16:46:50.902 1m 1s       38s         36.6%   115.8 MB    167.8 MB    364 GB      5.1 GB
33      37/ab1fcc 2045      exonerate (5)   COMPLETED   0        2014-10-23 16:47:51.625 30s         12s         59.6%   696 MB      734.6 MB    354.3 GB    420.4 MB
31      d7/eabe51 2042      exonerate (3)   COMPLETED   0        2014-10-23 16:45:50.846 3m 1s       2m 6s       130.1%  703.3 MB    760.9 MB    1.1 TB      28.6 GB
36      c4/d6cc15 2048      exonerate (6)   COMPLETED   0        2014-10-23 16:48:48.718 3m 1s       2m 43s      116.6%  682.1 MB    743.6 MB    868.5 GB    42 GB
30      4f/1ad1f0 2043      exonerate (2)   COMPLETED   0        2014-10-23 16:45:50.961 10m 2s      9m 16s      95.5%   706.2 MB    764 MB      1.6 TB      172.4 GB
52      72/41d0c6 2055      similarity (1)  COMPLETED   0        2014-10-23 17:13:23.543 30s         352ms       0.0%    35.6 MB     58.3 MB     199.3 MB    7.9 MB
57      9b/111b5e 2058      similarity (6)  COMPLETED   0        2014-10-23 17:13:23.655 30s         488ms       0.0%    108.2 MB    158 MB      317.1 MB    9.8 MB
53      3e/bca30f 2061      similarity (2)  COMPLETED   0        2014-10-23 17:13:23.770 30s         238ms       0.0%    6.7 MB      29.6 MB     190 MB      91.2 MB
54      8b/d45b47 2062      similarity (3)  COMPLETED   0        2014-10-23 17:13:23.808 30s         442ms       0.0%    108.1 MB    158 MB      832 MB      565.6 MB
55      51/ac19c6 2064      similarity (4)  COMPLETED   0        2014-10-23 17:13:23.873 30s         6s          0.0%    112.7 MB    162.8 MB    4.9 GB      3.9 GB
56      c3/ec5f4a 2066      similarity (5)  COMPLETED   0        2014-10-23 17:13:23.948 30s         616ms       0.0%    10.4 MB     34.6 MB     238 MB      8.4 MB
98      de/d6c0a6 2099      matrix (1)      COMPLETED   0        2014-10-23 17:14:27.139 30s         1s          0.0%    4.8 MB      42 MB       240.6 MB    79 KB
======= ========= ========= =============== =========== ======== ======================= =========== =========== ======= =========== =========== =========== ===========

.. _trace-fields:

The following table shows the fields that can be included in the execution report:

======================= ===============
Name                    Description
======================= ===============
task_id                 Task ID.
hash                    Task hash code.
native_id               Task ID given by the underlying execution system e.g. POSIX process PID when executed locally, job ID when executed by a grid engine, etc.
process                 Nextflow process name.
tag                     User provided identifier associated this task.
name                    Task name.
status                  Task status.
exit                    POSIX process exit status.
module                  Environment module used to run the task.
container               Docker image name used to execute the task.
cpus                    The cpus number request for the task execution.
time                    The time request for the task execution
disk                    The disk space request for the task execution.
memory                  The memory request for the task execution.
attempt                 Attempt at which the task completed.
submit                  Timestamp when the task has been submitted.
start                   Timestamp when the task execution has started.
complete                Timestamp when task execution has completed.
duration                Time elapsed to complete since the submission.
realtime                Task execution time i.e. delta between completion and start timestamp.
queue                   The queue that the executor attempted to run the process on.
%cpu                    Percentage of CPU used by the process.
%mem                    Percentage of memory used by the process.
rss                     Real memory (resident set) size of the process. Equivalent to ``ps -o rss`` .
vmem                    Virtual memory size of the process. Equivalent to ``ps -o vsize`` .
peak_rss                Peak of real memory. This data is read from field ``VmHWM`` in ``/proc/$pid/status`` file.
peak_vmem               Peak of virtual memory. This data is read from field ``VmPeak`` in ``/proc/$pid/status`` file.
rchar                   Number of bytes the process read, using any read-like system call from files, pipes, tty, etc. This data is read from file ``/proc/$pid/io``.
wchar                   Number of bytes the process wrote, using any write-like system call. This data is read from file ``/proc/$pid/io``.
syscr                   Number of read-like system call invocations that the process performed. This data is read from file ``/proc/$pid/io``.
syscw                   Number of write-like system call invocations that the process performed. This data is read from file ``/proc/$pid/io``.
read_bytes              Number of bytes the process directly read from disk. This data is read from file ``/proc/$pid/io``.
write_bytes             Number of bytes the process originally dirtied in the page-cache (assuming they will go to disk later). This data is read from file ``/proc/$pid/io``.
vol_ctxt                Number of voluntary context switches.
inv_ctxt                Number of involuntary context switches.
env                     The variables defined in task execution environment.
workdir                 The directory path where the task was executed.
script                  The task command script.
scratch                 The value of the process ``scratch`` directive.
error_action            The action applied on errof task failure.
======================= ===============

.. note:: These numbers provide an estimation of the resources used by running tasks. They should not be intended as an alternative
  to low level performance analysis provided by other tools and they may not be fully accurate, in particular for very short-lived tasks
  (running for less than one second).

Trace report layout and other configuration settings can be specified by using the ``nextflow.config`` configuration file.

Please read :ref:`Trace scope <config-trace>` section to learn more about it.

.. _timeline-report:

Timeline report
===============

Nextflow can render an HTML timeline for all processes executed in your pipeline. An example of the timeline
report is shown below:

.. image:: images/timeline-min.png


Each bar represents a process run in the pipeline execution. The bar length represents the task duration time (wall-time).
The colored area in each bar represents the real execution time. The grey area to the *left* of the colored area represents
the task scheduling wait time. The grey area to the *right* of the colored area represents the task termination time
(clean-up and file un-staging). The numbers on the x-axis represent the time in absolute units eg. minutes, hours, etc.

Each bar displays two numbers: the task duration time and the virtual memory size peak.

As each process can spawn many tasks, colors are used to identify those tasks belonging to the same process.


To enable the creation of the timeline report add the ``-with-timeline`` command line option when launching the pipeline
execution. For example::

  nextflow run <pipeline name> -with-timeline [file name]

The report file name can be specified as an optional parameter following the timeline option.

.. _dag-visualisation:

DAG visualisation
=================

A Nextflow pipeline is implicitly modelled by a direct acyclic graph (DAG). The vertices in the graph represent
the pipeline's processes and operators, while the edges represent the data connections (i.e. channels) between them.

The pipeline execution DAG can be outputted by adding the ``-with-dag`` option to the run command line.
It creates a file named ``dag.dot`` containing a textual representation of the pipeline execution graph
in the `DOT format <http://www.graphviz.org/content/dot-language>`_.

The execution DAG can be rendered in a different format by specifying an output file name which has an extension
corresponding to the required format. For example::

    nextflow run <script-name> -with-dag flowchart.png


List of supported file formats:

============ ====================
Extension     File format
============ ====================
dot           Graphviz DOT file
html          HTML file
pdf           PDF file (*)
png           PNG file (*)
svg           SVG file (*)
gexf          Graph Exchange XML file (Gephi)
============ ====================

.. warning:: The file formats marked with a `*` require the `Graphviz <http://www.graphviz.org>`_ tool installed
  in your computer.

The DAG produced by Nextflow for the `Shootstrap <https://github.com/cbcrg/shootstrap/>`_ pipeline:

.. image:: images/dag.png

.. _weblog-service:

Weblog via HTTP
===============

Nextflow is able to send detailed workflow execution metadata and runtime statistics to a HTTP endpoint.
To enable this feature use  the ``-with-weblog`` as shown below::

  nextflow run <pipeline name> -with-weblog [url]

Workflow events are sent as HTTP POST requests to the given URL. The message is formatted using the
following JSON structure::

   {
        "runName": <run name>,
        "runId": <uuid>,
        "event": <started|process_submitted|process_started|process_completed|error|completed>,
        "utcTime": <UTC timestamp>,
        "trace": { ... },
        "metadata": { ... }
   }

The JSON object contains the following attributes:

================== ================
Attribute          Description
================== ================
runName            The workflow execution run name.
runId              The workflow execution unique ID.
event              The workflow execution event. One of ``started``, ``process_submitted``, ``process_started``, ``process_completed``, ``error``, ``completed``.
utcTime            The UTC timestamp in ISO 8601 format.
trace              A process runtime information as described in the :ref:`trace fields<trace-fields>` section. This attribute is only provided for the following events: ``process_submitted``, ``process_started``, ``process_completed``, ``error``.
metadata           The workflow metadata including the :ref:`config manifest<config-manifest>`. For a list of all fields, have a look at the bottom message examples. This attribute is only provided for the following events: ``started``, ``completed``.
================== ================

.. warning::
  The content of the ``trace`` attribute depends on the settings for the `Trace report <trace-report>`_ defined in the
  ``nextflow.config`` file. See the :ref:`Trace configuration<config-trace>` section to learn more.


Weblog Started example message
------------------------------

When a workflow execution is started, a message like the following is posted to the specified end-point. Be aware that the
properties in the parameter scope will look different for your workflow. This is an example output from the ``nf-core/hlatyping``
pipeline with the weblog feature enabled::


  {
    "runName": "friendly_pesquet",
    "runId": "170aa09c-105f-49d0-99b4-8eb6a146e4a7",
    "event": "started",
    "utcTime": "2018-10-07T11:42:08Z",
    "metadata": {
            "params": {
                "container": "nfcore/hlatyping:1.1.4",
                "help": false,
                "outdir": "results",
                "bam": true,
                "singleEnd": false,
                "single-end": false,
                "reads": "data/test*{1,2}.fq.gz",
                "seqtype": "dna",
                "solver": "glpk",
                "igenomes_base": "./iGenomes",
                "multiqc_config": "/Users/sven1103/.nextflow/assets/nf-core/hlatyping/conf/multiqc_config.yaml",
                "clusterOptions": false,
                "cluster-options": false,
                "enumerations": 1,
                "beta": 0.009,
                "prefix": "hla_run",
                "base_index": "/Users/sven1103/.nextflow/assets/nf-core/hlatyping/data/indices/yara/hla_reference_",
                "index": "/Users/sven1103/.nextflow/assets/nf-core/hlatyping/data/indices/yara/hla_reference_dna",
                "custom_config_version": "master",
                "custom_config_base": "https://raw.githubusercontent.com/nf-core/configs/master"
            },
            "workflow": {
                "start": "2019-03-25T12:09:52Z",
                "projectDir": "/Users/sven1103/.nextflow/assets/nf-core/hlatyping",
                "manifest": {
                    "nextflowVersion": ">=18.10.1",
                    "defaultBranch": "master",
                    "version": "1.1.4",
                    "homePage": "https://github.com/nf-core/hlatyping",
                    "gitmodules": null,
                    "description": "Precision HLA typing from next-generation sequencing data.",
                    "name": "nf-core/hlatyping",
                    "mainScript": "main.nf",
                    "author": null
                },
                "complete": null,
                "profile": "docker,test",
                "homeDir": "/Users/sven1103",
                "workDir": "/Users/sven1103/git/nextflow/work",
                "container": "nfcore/hlatyping:1.1.4",
                "commitId": "4bcced898ee23600bd8c249ff085f8f88db90e7c",
                "errorMessage": null,
                "repository": "https://github.com/nf-core/hlatyping.git",
                "containerEngine": "docker",
                "scriptFile": "/Users/sven1103/.nextflow/assets/nf-core/hlatyping/main.nf",
                "userName": "sven1103",
                "launchDir": "/Users/sven1103/git/nextflow",
                "runName": "shrivelled_cantor",
                "configFiles": [
                    "/Users/sven1103/.nextflow/assets/nf-core/hlatyping/nextflow.config"
                ],
                "sessionId": "7f344978-999c-480d-8439-741bc7520f6a",
                "errorReport": null,
                "scriptId": "2902f5aa7f297f2dccd6baebac7730a2",
                "revision": "master",
                "exitStatus": null,
                "commandLine": "./launch.sh run nf-core/hlatyping -profile docker,test -with-weblog 'http://localhost:4567'",
                "nextflow": {
                              "version": "19.03.0-edge",
                              "build": 5137,
                              "timestamp": "2019-03-28T14:46:55Z"
                            },
                },
                "stats": {
                    "computeTimeFmt": "(a few seconds)",
                    "cachedCount": 0,
                    "cachedDuration": 0,
                    "failedDuration": 0,
                    "succeedDuration": 0,
                    "failedCount": 0,
                    "cachedPct": 0.0,
                    "cachedCountFmt": "0",
                    "succeedCountFmt": "0",
                    "failedPct": 0.0,
                    "failedCountFmt": "0",
                    "ignoredCountFmt": "0",
                    "ignoredCount": 0,
                    "succeedPct": 0.0,
                    "succeedCount": 0,
                    "ignoredPct": 0.0
                },
                "resume": false,
                "success": false,
                "scriptName": "main.nf",
                "duration": null
            }
        }
  }


Weblog Completed example message
--------------------------------

Once a process is completed, a message like the following is posted to the specified end-point::

  {
    "runName": "friendly_pesquet",
    "runId": "170aa09c-105f-49d0-99b4-8eb6a146e4a7",
    "event": "process_completed",
    "utcTime": "2018-10-07T11:45:30Z",
    "trace": {
        "task_id": 2,
        "status": "COMPLETED",
        "hash": "a1/0024fd",
        "name": "make_ot_config",
        "exit": 0,
        "submit": 1538912529498,
        "start": 1538912529629,
        "process": "make_ot_config",
        "tag": null,
        "module": [

        ],
        "container": "nfcore/hlatyping:1.1.1",
        "attempt": 1,
        "script": "\n    configbuilder --max-cpus 2 --solver glpk > config.ini\n    ",
        "scratch": null,
        "workdir": "/home/sven1103/git/hlatyping-workflow/work/a1/0024fd028375e2b601aaed44d112e3",
        "queue": null,
        "cpus": 1,
        "memory": 7516192768,
        "disk": null,
        "time": 7200000,
        "env": "PATH=/home/sven1103/git/hlatyping-workflow/bin:$PATH\n",
        "error_action": null,
        "complete": 1538912730599,
        "duration": 201101,
        "realtime": 69,
        "%cpu": 0.0,
        "%mem": 0.1,
        "vmem": 54259712,
        "rss": 10469376,
        "peak_vmem": 20185088,
        "peak_rss": 574972928,
        "rchar": 7597,
        "wchar": 162,
        "syscr": 16,
        "syscw": 4083712,
        "read_bytes": 4096,
        "write_bytes": 0,
        "native_id": 27185
    }
  }
.. _dsl2-page:

******
DSL 2
******

Nextflow provides a syntax extension that allows the definition of module libraries and
simplifies the writing of complex data analysis pipelines.

To enable this feature you need to define the following directive at the beginning of
your workflow script::

    nextflow.enable.dsl=2


Function
========

Nextflow allows the definition of custom functions in the workflow script using the following syntax::

    def <function name> ( arg1, arg, .. ) {
        <function body>
    }

For example::

    def foo() {
        'Hello world'
    }

    def bar(alpha, omega) {
        alpha + omega
    }


The above snippet defines two simple functions, that can be invoked in the workflow script as ``foo()`` which
returns the ``Hello world`` string and ``bar(10,20)`` which returns the sum of two parameters (``30`` in this case).

.. tip:: Functions implicitly return the result of the last evaluated statement.

The keyword ``return`` can be used to explicitly exit from a function and return the specified value.
For example::

    def fib( x ) {
        if( x <= 1 )
            return x
        else
            fib(x-1) + fib(x-2)
    }

Process
=======

Process definition
------------------

The new DSL separates the definition of a process from its invocation. The process definition follows the usual 
syntax as described in the :ref:`process documentation <process-page>`. The only difference is that the
``from`` and ``into`` channel declarations have to be omitted.

Then a process can be invoked as a function in the ``workflow`` scope, passing the expected
input channels as parameters as if it were a custom function. For example::

    nextflow.enable.dsl=2

    process foo {
        output:
          path 'foo.txt'
        script:
          """
          your_command > foo.txt
          """
    }

     process bar {
        input:
          path x
        output:
          path 'bar.txt'
        script:
          """
          another_command $x > bar.txt
          """
    }

    workflow {
        data = channel.fromPath('/some/path/*.txt')
        foo()
        bar(data)
    }


.. warning::
  A process component can be invoked only once in the same workflow context.


Process composition
-------------------

Processes having matching *input-output* declaration can be composed so that the output
of the first process is passed as input to the next process. Taking in consideration
the previous example, it's possible to write the following::

    workflow {
        bar(foo())
    }


Process output
---------------

A process output can also be accessed using the ``out`` attribute on the corresponding
process object. For example::

    workflow {
        foo()
        bar(foo.out)
        bar.out.view()
    }


When a process defines two or more output channels, each of them can be accessed
using the array element operator e.g. ``out[0]``, ``out[1]``, etc. or using
*named outputs* (see below).

Process named output
--------------------

The ``emit`` option can be added to the process output definition to assign a name identifier. This name
can be used to reference the channel within the caller scope. For example::

    process foo {
      output:
        path '*.bam', emit: samples_bam

      '''
      your_command --here
      '''
    }
    
    workflow {
        foo()
        foo.out.samples_bam.view()
    }
    
Process named stdout
--------------------

The ``emit`` option can be used also to name the stdout::

    process sayHello {
        input:
            val cheers
        output:
            stdout emit: verbiage
        script:
        """
        echo -n $cheers
        """
    }

    workflow {
        things = channel.of('Hello world!', 'Yo, dude!', 'Duck!')
        sayHello(things)
        sayHello.out.verbiage.view()
    }

Workflow
========

Workflow definition
--------------------

The ``workflow`` keyword allows the definition of sub-workflow components that enclose the
invocation of one or more processes and operators::

    workflow my_pipeline {
        foo()
        bar( foo.out.collect() )
    }


For example, the above snippet defines a workflow component, named ``my_pipeline``, that can be invoked from
another workflow component definition as any other function or process with ``my_pipeline()``.


Workflow parameters
---------------------

A workflow component can access any variable and parameter defined in the outer scope::

        params.data = '/some/data/file'

        workflow my_pipeline {
            if( params.data )
                bar(params.data)
            else
                bar(foo())
        }


Workflow input
---------------

A workflow component can declare one or more input channels using the ``take`` keyword. For example::

        workflow my_pipeline {
            take: data
            main:
              foo(data)
              bar(foo.out)
        }

.. warning::
  When the ``take`` keyword is used, the beginning of the workflow body must be identified with the
  ``main`` keyword.

Then, the input can be specified as an argument in the workflow invocation statement::

    workflow {
        my_pipeline( channel.from('/some/data') )
    }

.. note::
  Workflow inputs are by definition *channel* data structures. If a basic data type is provided
  instead, ie. number, string, list, etc. it's implicitly converted to a :ref:`channel value <channel-type-value>`
  (ie. non-consumable).


Workflow output
----------------

A workflow component can declare one or more output channels using the ``emit`` keyword. For example::

    workflow my_pipeline {
        main:
          foo(data)
          bar(foo.out)
        emit:
          bar.out
    }

Then, the result of the ``my_pipeline`` execution can be accessed using the ``out`` property, i.e.
``my_pipeline.out``. When multiple output channels are declared, use the array bracket notation
to access each output channel as described for the `Process output`_ definition.

Workflow named output
---------------------
If the output channel is assigned to an identifier in the ``emit`` declaration, such identifier can be used
to reference the channel within the caller scope. For example::

     workflow my_pipeline {
        main:
          foo(data)
          bar(foo.out)
        emit:
          my_data = bar.out
     }

Then, the result of the above snippet can accessed using ``my_pipeline.out.my_data``.


Entry point of execution
------------------------

A workflow definition which does not declare any name (also known as *implicit workflow*) is
the entry point of execution for the workflow application.

.. note::
  Implicit workflow definition is ignored when a script is included as module. This
  allows the writing of a workflow script that can be used either as a library module and as
  application script. 

.. tip::
  An alternative named workflow entry can be specified using the ``-entry`` command line option.


Workflow composition
--------------------

Workflows defined in your script or imported with `Module inclusion`_ can be invoked and composed
as any other process in your application.

::

    workflow flow1 {
        take: data
        main:
            foo(data)
            bar(foo.out)
        emit:
            bar.out
    }

    workflow flow2 {
        take: data
        main:
            foo(data)
            baz(foo.out)
        emit:
            baz.out
    }

    workflow {
        take: data
        main:
          flow1(data)
          flow2(flow1.out)
    }


.. note::
    Nested workflow execution determines an implicit scope. Therefore the same process can be
    invoked in two different workflow scopes, like for example ``foo`` in the above snippet that
    is used either in ``flow1`` and ``flow2``. The workflow execution path along with the
    process names defines the process *fully qualified name* that is used to distinguish the
    two different process invocations, i.e. ``flow1:foo`` and ``flow2:foo`` in the above example.

.. tip::
    The process fully qualified name can be used as a valid :ref:`process selector <config-process-selectors>` in the
    ``nextflow.config`` file and it has priority over the process simple name.


Modules
=======

The new DSL allows the definition of *module scripts* that
can be included and shared across workflow applications.

A module script (or simply, module) can contain the definition of functions, processes and workflows
as described in the previous sections.

.. note::
    Functions, processes and workflows are globally referred as *components*.

Module inclusion
----------------

A component defined in a module script can be imported into another Nextflow script using the ``include`` keyword.

For example::

    include { foo } from './some/module'

    workflow {
        data = channel.fromPath('/some/data/*.txt')
        foo(data)
    }

The above snippet includes a process with name ``foo`` defined in the module script in the main
execution context. This way, `foo`` can be invoked in the ``workflow`` scope.

Nextflow implicitly looks for the script file ``./some/module.nf`` resolving the path
against the *including* script location.

.. note:: Relative paths must begin with the ``./`` prefix. Also, the ``include`` statement must be defined **outside** of the workflow definition.

Multiple inclusions
-------------------

A Nextflow script allows the inclusion of an arbitrary number of modules and components. When multiple
components need to be included from the same module script, the component names can be
specified in the same inclusion using the curly brackets notation as shown below::

    include { foo; bar } from './some/module'

    workflow {
        data = channel.fromPath('/some/data/*.txt')
        foo(data)
        bar(data)
    }


Module aliases
--------------

When including a module component, it's possible to specify an *alias* with the ``as`` keyword.
This allows the inclusion and the invocation of components with the same name
in your script using different names. For example::

    include { foo } from './some/module'
    include { foo as bar } from './other/module'

    workflow {
        foo(some_data)
        bar(other_data)
    }

The same is possible when including the same component multiple times from the same module script as shown below::

    include { foo; foo as bar } from './some/module'

    workflow {
        foo(some_data)
        bar(other_data)
    }


Module parameters
-----------------

A module script can define one or more parameters using the same syntax of a Nextflow workflow script::

    params.foo = 'Hello'
    params.bar = 'world!'

    def sayHello() {
        println "$params.foo $params.bar"
    }


Then, parameters are inherited from the including context. For example::

    params.foo = 'Hola'
    params.bar = 'Mundo'

    include {sayHello} from './some/module'

    workflow {
        sayHello()
    }

The above snippet prints::

    Hola Mundo


.. note::
  The module inherits the parameters defined *before* the ``include`` statement, therefore any further
  parameter set later is ignored.

.. tip::
  Define all pipeline parameters at the beginning of the script *before*
  any ``include`` declaration.

The option ``addParams`` can be used to extend the module parameters without affecting the external
scope. For example::


    include {sayHello} from './some/module' addParams(foo: 'Ciao')

    workflow {
        sayHello()
    }


The above snippet prints::

    Ciao world!


Finally, the include option ``params`` allows the specification of one or more parameters without
inheriting any value from the external environment. 

.. _module-templates:

Module templates
-----------------
The module script can be defined in an external :ref:`template <process-template>` file. With DSL2 the template file
can be placed under the ``templates`` directory where the module script is located.

For example, let's suppose to have a project L with a module script defining 2 processes (P1 and P2) and both use templates.
The template files can be made available under the local ``templates`` directory::

	Project L
		|-myModules.nf
		|-templates
			|-P1-template.sh
			|-P2-template.sh

Then, we have a second project A with a workflow that includes P1 and P2::

	Pipeline A
		|-main.nf

Finally, we have a third project B with a workflow that includes again P1 and P2::

	Pipeline B
		|-main.nf

With the possibility to keep the template files inside the project L, A and B can use the modules defined in L without any changes.
A future prject C would do the same, just cloning L (if not available on the system) and including its module script.

Beside promoting sharing modules across pipelines, there are several advantages in keeping the module template under the script path::
1 - module components are *self-contained*
2 - module components can be tested independently from the pipeline(s) importing them
3 - it is possible to create libraries of module components

Ultimately, having multiple template locations allows a more structured organization within the same project. If a project
has several module components, and all them use templates, the project could group module scripts and their templates as needed. For example::

	baseDir
		|-main.nf
		|-Phase0-Modules
			|-mymodules1.nf
			|-mymodules2.nf
			|-templates
				|-P1-template.sh
				|-P2-template.sh
		|-Phase1-Modules
			|-mymodules3.nf
			|-mymodules4.nf
			|-templates
				|-P3-template.sh
				|-P4-template.sh
		|-Phase2-Modules
			|-mymodules5.nf
			|-mymodules6.nf
			|-templates
				|-P5-template.sh
				|-P6-template.sh
				|-P7-template.sh


Channel forking
===============

Using the new DSL, Nextflow channels are automatically forked when connecting two or more consumers.

For example::

    channel
        .from('Hello','Hola','Ciao')
        .set{ cheers }

    cheers
        .map{ it.toUpperCase() }
        .view()

    cheers
        .map{ it.reverse() }
        .view()


The same is valid for the result (channel) of a process execution. Therefore a process output can be consumed by
two or more processes without the need to fork it using the :ref:`operator-into` operator, making the
writing of workflow scripts more fluent and readable.


Pipes
=====

The *pipe* operator
-------------------

Nextflow processes and operators can be composed using the ``|`` *pipe* operator. For example::

    process foo {
        input: val data
        output: val result
        exec:
          result = "$data world"
    }

    workflow {
       channel.from('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
    }



The above snippet defines a process named ``foo`` and invokes it passing the content of the
``data`` channel. The result is then piped to the :ref:`operator-map` operator which converts each string
to uppercase and finally, the last :ref:`operator-view` operator prints it.


The *and* operator
------------------

The ``&`` *and* operator allows feeding of two or more processes with the content of the same
channel(s). For example::

    process foo {
      input: val data
      output: val result
      exec:
        result = "$data world"
    }

    process bar {
        input: val data
        output: val result
        exec:
          result = data.toUpperCase()
    }

    workflow {
       channel.from('Hello') | map { it.reverse() } | (foo & bar) | mix | view
    }


In the above snippet the channel emitting the ``Hello`` string is piped with the :ref:`operator-map`
which reverses the string value. Then, the result is passed to both ``foo`` and ``bar``
processes which are executed in parallel. Each process outputs a channel, and the two channels are merged
into a single channel using the :ref:`operator-mix` operator. Finally the result is printed
using the :ref:`operator-view` operator.

.. tip::
  The break-line operator ``\`` can be used to split long pipe concatenations
  over multiple lines.


The above snippet can be written as shown below::

    workflow {
       channel.from('Hello') \
         | map { it.reverse() } \
         | (foo & bar) \
         | mix \
         | view
    }



DSL2 migration notes
=====================

* DSL2 final version is activated using the declaration ``nextflow.enable.dsl=2`` in place of ``nextflow.preview.dsl=2``.
* Process inputs of type ``set`` have to be replaced with :ref:`tuple <process-input-tuple>`.
* Process outputs of type ``set`` have to be replaced with :ref:`tuple <process-out-tuple>`.
* Process output option ``mode flatten`` is no longer available. Replace it using the :ref:`operator-flatten` operator on the corresponding output channel.
* Anonymous and unwrapped includes are not supported anymore. Replace them with an explicit module inclusion. For example::

        include './some/library'
        include bar from './other/library'

        workflow {
          foo()
          bar()
        }

  Should be replaced with::

        include { foo } from './some/library'
        include { bar } from './other/library'

        workflow {
          foo()
          bar()
        }
        
* The use of unqualified value and file elements into input tuples is not allowed anymore. Replace them with a corresponding
  ``val`` or ``path`` qualifier::

        process foo {
        input:
          tuple X, 'some-file.bam'
         script:
           '''
           your_command --in $X some-file.bam
           '''
        }

  Use::

        process foo {
        input:
          tuple val(X), path('some-file.bam')
         script:
           '''
           your_command --in $X some-file.bam
           '''
        }


* The use of unqualified value and file elements into output tuples is not allowed anymore. Replace them with a corresponding
  ``val`` or ``path`` qualifier::


        process foo {
        output:
          tuple X, 'some-file.bam'

        script:
           X = 'some value'
           '''
           your_command > some-file.bam
           '''
        }

  Use::

        process foo {
        output:
          tuple val(X), path('some-file.bam')

        script:
           X = 'some value'
           '''
           your_command > some-file.bam
           '''
        }


* Operator :ref:`channel-bind1` has been deprecated by DSL2 syntax
* Operator :ref:`channel-bind2` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-choice` has been deprecated by DSL2 syntax. Use :ref:`operator-branch` instead.
* Operator :ref:`operator-close` has been deprecated by DSL2 syntax.
* Operator :ref:`channel-create` has been deprecated by DSL2 syntax.
* Operator ``countBy`` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-into` has been deprecated by DSL2 syntax since it's not needed anymore.
* Operator ``fork`` has been renamed to :ref:`operator-multimap`.
* Operator ``groupBy`` has been deprecated by DSL2 syntax. Replace it with :ref:`operator-grouptuple`
* Operator ``print`` and ``println`` have been deprecated by DSL2 syntax. Use :ref:`operator-view` instead.
* Operator :ref:`operator-merge` has been deprecated by DSL2 syntax. Use :ref:`operator-join` instead.
* Operator :ref:`operator-separate` has been deprecated by DSL2 syntax.
* Operator :ref:`operator-spread` has been deprecated with DSL2 syntax. Replace it with :ref:`operator-combine`.
* Operator route has been deprecated by DSL2 syntax.

.. _config-page:

*************
Configuration
*************

Configuration file
==================

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations.
Since each configuration file can contain conflicting settings, the sources are ranked to decide which
settings to are applied. All possible configuration sources are reported below, listed in order
of priority:

1. Parameters specified on the command line (``--something value``)
2. Parameters provided using the ``-params-file`` option
3. Config file specified using the ``-c my_config`` option
4. The config file named ``nextflow.config`` in the current directory
5. The config file named ``nextflow.config`` in the workflow project directory
6. The config file ``$HOME/.nextflow/config``
7. Values defined within the pipeline script itself (e.g. ``main.nf``)

When more than one of these ways of specifying configurations are used, they are merged, so that the settings in the
first override the same ones that may appear in the second one, and so on.

.. tip:: If you want to ignore any default configuration files and use only the custom one use the command line option
  ``-C <config file>``.

Config syntax
--------------

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax::

  name = value

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (``true``, ``false``) do not.
Also note that values are typed, meaning for example that, ``1`` is different from ``'1'``, since the first is interpreted
as the number one, while the latter is interpreted as a string value.


Config Variables
----------------

Configuration properties can be used as variables in the configuration file itself, by using the usual
``$propertyName`` or ``${expression}`` syntax.


For example::

     propertyOne = 'world'
     anotherProp = "Hello $propertyOne"
     customPath = "$PATH:/my/app/folder"

Please note, the usual rules for :ref:`string-interpolation` are applied, thus a string containing a variable
reference must be wrapped in double-quote chars instead of single-quote chars.

The same mechanism allows you to access environment variables defined in the hosting system. Any variable whose name is
not defined in the Nextflow configuration file(s) is supposed to be a reference to an environment variable with that name.
So, in the above example the property ``customPath`` is defined as the current system ``PATH`` to which
the string ``/my/app/folder`` is appended.


Config comments
------------------

Configuration files use the same conventions for comments used by the Groovy or Java programming languages. Thus, use ``//`` to comment
a single line or ``/*`` .. ``*/`` to comment a block on multiple lines.


Config include
--------------

A configuration file can include one or more configuration files using the keyword ``includeConfig``. For example::

    process.executor = 'sge'
    process.queue = 'long'
    process.memory = '10G'

    includeConfig 'path/foo.config'

When a relative path is used, it is resolved against the actual location of the including file.


Config scopes
=============

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope
identifier or grouping the properties in the same scope using the curly brackets notation. This is shown in the
following example::

   alpha.x  = 1
   alpha.y  = 'string value..'

   beta {
        p = 2
        q = 'another string ..'
   }


.. _config-env:

Scope `env`
-----------

The ``env`` scope allows the definition one or more variable that will be exported in the environment where the
workflow tasks will be executed.

Simply prefix your variable names with the ``env`` scope or surround them by curly brackets, as shown below::

   env.ALPHA = 'some value'
   env.BETA = "$HOME/some/path"

   env {
        DELTA = 'one more'
        GAMMA = "/my/path:$PATH"
   }

.. tip:: In the above example, variables like `$HOME` and `$PATH` are evaluated when the workflow is launched. If
  you want these variables to be evaluated during task execution, escape them with `\$`. This difference is important
  for variables like `$PATH`, which may be very different in the workflow environment versus the task environment.


Scope `params`
--------------

The ``params`` scope allows you to define parameters that will be accessible in the pipeline script. Simply prefix the
parameter names with the ``params`` scope or surround them by curly brackets, as shown below::

     params.custom_param = 123
     params.another_param = 'string value .. '

     params {

        alpha_1 = true
        beta_2 = 'another string ..'

     }



.. _config-process:

Scope `process`
---------------

The ``process`` configuration scope allows you to provide the default configuration for the processes in your pipeline.

You can specify here any property described in the :ref:`process directive<process-directives>` and the executor sections.
For examples::

  process {
    executor='sge'
    queue='long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
  }


By using this configuration all processes in your pipeline will be executed through the SGE cluster, with the specified
settings.

.. _config-process-selectors:

Process selectors
^^^^^^^^^^^^^^^^^

The ``withLabel`` selectors allow the configuration of all processes annotated with a :ref:`process-label` directive as
shown below::

    process {
        withLabel: big_mem {
            cpus = 16
            memory = 64.GB
            queue = 'long'
        }
    }

The above configuration example assigns 16 cpus, 64 Gb of memory and the ``long`` queue to all processes annotated
with the ``big_mem`` label.


In the same manner, the ``withName`` selector allows the configuration of a specific process in your pipeline by its name.
For example::

    process {
        withName: hello {
            cpus = 4
            memory = 8.GB
            queue = 'short'
        }
    }

.. tip:: Either label and process names do not need to be enclosed with quote characters, provided the name
  does include special characters (e.g. ``-``, ``!``, etc) or it's not a keyword or a built-in type identifier.
  In case of doubt, you can enclose the label names or the process names with single or double quote characters.

.. _config-selector-expressions:

Selector expressions
^^^^^^^^^^^^^^^^^^^^

Both label and process name selectors allow the use of a regular expression in order to apply the same configuration
to all processes matching the specified pattern condition. For example::

    process {
        withLabel: 'foo|bar' {
            cpus = 2
            memory = 4.GB
        }
    }

The above configuration snippet sets 2 cpus and 4 GB of memory to the processes annotated with with a label ``foo``
and ``bar``.

A process selector can be negated prefixing it with the special character ``!``. For example::

    process {
        withLabel: 'foo' { cpus = 2 }
        withLabel: '!foo' { cpus = 4 }
        withName: '!align.*' { queue = 'long' }
    }

The above configuration snippet sets 2 cpus for the processes annotated with the ``foo`` label and 4 cpus to all processes
*not* annotated with that label. Finally it sets the use of ``long`` queue to all process whose name does *not* start
with ``align``.

.. _config-selectors-priority:

Selectors priority
^^^^^^^^^^^^^^^^^^

When mixing generic process configuration and selectors the following priority rules are applied (from lower to higher):

1. Process generic configuration.
2. Process specific directive defined in the workflow script.
3. ``withLabel`` selector definition.
4. ``withName`` selector definition.

For example::

    process {
        cpus = 4
        withLabel: foo { cpus = 8 }
        withName: bar { cpus = 32 }
    }

Using the above configuration snippet, all workflow processes use 4 cpus if not otherwise specified in the workflow
script. Moreover processes annotated with the ``foo`` label use 8 cpus. Finally the process named ``bar``
uses 32 cpus.


.. _config-executor:

Scope `executor`
----------------

The ``executor`` configuration scope allows you to set the optional executor settings, listed in the following table.

===================== =====================
Name                  Description
===================== =====================
name                  The name of the executor to be used e.g. ``local``, ``sge``, etc.
queueSize             The number of tasks the executor will handle in a parallel manner (default: ``100``).
pollInterval          Determines how often a poll occurs to check for a process termination.
dumpInterval          Determines how often the executor status is written in the application log file (default: ``5min``).
queueStatInterval     Determines how often the queue status is fetched from the cluster system. This setting is used only by grid executors (default: ``1min``).
exitReadTimeout       Determines how long the executor waits before return an error status when a process is terminated but the `exit` file does not exist or it is empty. This setting is used only by grid executors (default: ``270 sec``).
killBatchSize         Determines the number of jobs that can be `killed` in a single command execution (default: ``100``).
submitRateLimit       Determines the max rate of job submission per time unit, for example ``'10sec'`` eg. max 10 jobs per second or ``'50/2min'`` i.e. 50 job submissions every 2 minutes (default: `unlimited`).
perJobMemLimit        Specifies Platform LSF *per-job* memory limit mode. See :ref:`lsf-executor`.
jobName               Determines the name of jobs submitted to the underlying cluster executor e.g. ``executor.jobName = { "$task.name - $task.hash" }`` Note: when using this option you need to make sure the resulting job name matches the validation constraints of the underlying batch scheduler.
cpus                  The maximum number of CPUs made available by the underlying system (only used by the ``local`` executor).
memory                The maximum amount of memory made available by the underlying system (only used by the ``local`` executor).
===================== =====================



The executor settings can be defined as shown below::

    executor {
        name = 'sge'
        queueSize = 200
        pollInterval = '30 sec'
    }


When using two (or more) different executors in your pipeline, you can specify their settings separately by prefixing
the executor name with the symbol ``$`` and using it as special scope identifier. For example::

  executor {
    $sge {
        queueSize = 100
        pollInterval = '30sec'
    }

    $local {
        cpus = 8
        memory = '32 GB'
    }
  }

The above configuration example can be rewritten using the dot notation as shown below::

    executor.$sge.queueSize = 100
    executor.$sge.pollInterval = '30sec'
    executor.$local.cpus = 8
    executor.$local.memory = '32 GB'

.. _config-docker:

Scope `docker`
--------------

The ``docker`` configuration scope controls how `Docker <https://www.docker.com>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Docker execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
legacy              Uses command line options removed since version 1.10.x (default: ``false``).
sudo                Executes Docker run command as ``sudo`` (default: ``false``).
tty                 Allocates a pseudo-tty (default: ``false``).
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
remove              Clean-up the container after the execution (default: ``true``). For details see: https://docs.docker.com/engine/reference/run/#clean-up---rm .
runOptions          This attribute can be used to provide any extra command line options supported by the ``docker run`` command. For details see: https://docs.docker.com/engine/reference/run/ .
registry            The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. ``http://``.
fixOwnership        Fixes ownership of files created by the docker container.
engineOptions       This attribute can be used to provide any option supported by the Docker engine i.e. ``docker [OPTIONS]``.
mountFlags          Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`
================== ================

The above options can be used by prefixing them with the ``docker`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    docker {
        enabled = true
        temp = 'auto'
    }



Read :ref:`docker-page` page to learn more how use Docker containers with Nextflow.


.. _config-singularity:

Scope `singularity`
-------------------

The ``singularity`` configuration scope controls how `Singularity <https://sylabs.io/singularity/>`_ containers are executed
by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Singularity execution (default: ``false``).
engineOptions       This attribute can be used to provide any option supported by the Singularity engine i.e. ``singularity [OPTIONS]``.
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
runOptions          This attribute can be used to provide any extra command line options supported by the ``singularity exec``.
noHttps             Turn this flag to ``true`` to pull the Singularity image with http protocol (default: ``false``).
autoMounts          When ``true`` Nextflow automatically mounts host paths in the executed container. It requires the `user bind control` feature enabled in your Singularity installation (default: ``false``).
cacheDir            The directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all computing nodes.
pullTimeout         The amount of time the Singularity pull can last, exceeding which the process is terminated (default: ``20 min``).
================== ================


Read :ref:`singularity-page` page to learn more how use Singularity containers with Nextflow.

.. _config-podman:

Scope `podman`
--------------

The ``podman`` configuration scope controls how `Podman <https://podman.io/>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Podman execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
remove              Clean-up the container after the execution (default: ``true``).
runOptions          This attribute can be used to provide any extra command line options supported by the ``podman run`` command.
registry            The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. ``http://``.
engineOptions       This attribute can be used to provide any option supported by the Docker engine i.e. ``podman [OPTIONS]``.
mountFlags          Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`
================== ================

The above options can be used by prefixing them with the ``podman`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    podman {
        enabled = true
        temp = 'auto'
    }



Read :ref:`podman-page` page to learn more how use Podman containers with Nextflow.

.. _config-charliecloud:

Scope `charliecloud`
--------------------

The ``charliecloud`` configuration scope controls how `Charliecloud <https://hpc.github.io/charliecloud/>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Charliecloud execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
runOptions          This attribute can be used to provide any extra command line options supported by the ``ch-run`` command.
cacheDir            The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all computing nodes.
pullTimeout         The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: ``20 min``).
================== ================

The above options can be used by prefixing them with the ``charliecloud`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    charliecloud {
        enabled = true
    }



Read :ref:`charliecloud-page` page to learn more how use Charliecloud containers with Nextflow.

.. _config-manifest:

Scope `manifest`
----------------

The ``manifest`` configuration scope allows you to define some meta-data information needed when publishing your pipeline project on GitHub, BitBucket or GitLab, or when running your pipeline.

The following settings are available:

================== ================
Name                Description
================== ================
author              Project author name (use a comma to separate multiple names).
defaultBranch       Git repository default branch (default: ``master``).
recurseSubmodules   Turn this flag to ``true`` to pull submodules recursively from the Git repository
description         Free text describing the workflow project.
doi                 Project related publication DOI identifier.
homePage            Project home page URL.
mainScript          Project main script (default: ``main.nf``).
name                Project short name.
nextflowVersion     Minimum required Nextflow version.
version             Project version number.
================== ================

The above options can be used by prefixing them with the ``manifest`` scope or surrounding them by curly
brackets. For example::

    manifest {
        homePage = 'http://foo.com'
        description = 'Pipeline does this and that'
        mainScript = 'foo.nf'
        version = '1.0.0'
    }


To learn how to publish your pipeline on GitHub, BitBucket or GitLab code repositories read :ref:`sharing-page`
documentation page.

Nextflow version
^^^^^^^^^^^^^^^^

The ``nextflowVersion`` setting allows you to specify a minimum required version to run the pipeline.
This may be useful to ensure that a specific version is used::

    nextflowVersion = '1.2.3'        // exact match
    nextflowVersion = '1.2+'         // 1.2 or later (excluding 2 and later)
    nextflowVersion = '>=1.2'        // 1.2 or later
    nextflowVersion = '>=1.2, <=1.5' // any version in the 1.2 .. 1.5 range
    nextflowVersion = '!>=1.2'       // with ! prefix, stop execution if current version
                                        does not match required version.

.. _config-trace:

Scope `trace`
-------------

The ``trace`` scope allows you to control the layout of the execution trace file generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the execution trace report file (default: ``false``).
fields              Comma separated list of fields to be included in the report. The available fields are listed at :ref:`this page <trace-fields>`
file                Trace file name (default: ``trace.txt``).
sep                 Character used to separate values in each row (default: ``\t``).
raw                 When ``true`` turns on raw number report generation i.e. date and time are reported as milliseconds and memory as number of bytes
overwrite           When ``true`` overwrites an existing trace file instead of rolling it.
================== ================

The above options can be used by prefixing them with the ``trace`` scope or surrounding them by curly
brackets. For example::

    trace {
        enabled = true
        file = 'pipeline_trace.txt'
        fields = 'task_id,name,status,exit,realtime,%cpu,rss'
    }


To learn more about the execution report that can be generated by Nextflow read :ref:`trace-report` documentation page.

.. _config-dag:

Scope `dag`
-------------

The ``dag`` scope allows you to control the layout of the execution graph file generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the execution graph report file (default: ``false``).
file                Graph file name (default: ``dag.dot``).
================== ================

The above options can be used by prefixing them with the ``dag`` scope or surrounding them by curly
brackets. For example::

    dag {
        enabled = true
        file = 'pipeline_dag.html'
    }


To learn more about the execution graph that can be generated by Nextflow read :ref:`dag-visualisation` documentation page.

.. _config-aws:

Scope `aws`
-----------

The ``aws`` scope allows you to configure the access to Amazon S3 storage. Use the attributes ``accessKey`` and ``secretKey``
to specify your bucket credentials. For example::


    aws {
        accessKey = '<YOUR S3 ACCESS KEY>'
        secretKey = '<YOUR S3 SECRET KEY>'
        region = '<REGION IDENTIFIER>'
    }

Click the following link to learn more about `AWS Security Credentials <http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html>`_.

Advanced client configuration options can be set by using the ``client`` attribute. The following properties can be used:

=========================== ================
Name                        Description
=========================== ================
s3Acl                       Allow the setting of a predefined bucket permissions also known as *canned ACL*. Permitted values are ``Private``, ``PublicRead``, ``PublicReadWrite``, ``AuthenticatedRead``, ``LogDeliveryWrite``, ``BucketOwnerRead``, ``BucketOwnerFullControl`` and ``AwsExecRead``. See `Amazon docs <https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl>`_ for details.
connectionTimeout           The amount of time to wait (in milliseconds) when initially establishing a connection before giving up and timing out.
endpoint                    The AWS S3 API entry point e.g. `s3-us-west-1.amazonaws.com`.
maxConnections              The maximum number of allowed open HTTP connections.
maxErrorRetry               The maximum number of retry attempts for failed retryable requests.
protocol                    The protocol (i.e. HTTP or HTTPS) to use when connecting to AWS.
proxyHost                   The proxy host to connect through.
proxyPort                   The port on the proxy host to connect through.
proxyUsername               The user name to use when connecting through a proxy.
proxyPassword               The password to use when connecting through a proxy.
signerOverride              The name of the signature algorithm to use for signing requests made by the client.
socketSendBufferSizeHint    The Size hint (in bytes) for the low level TCP send buffer.
socketRecvBufferSizeHint    The Size hint (in bytes) for the low level TCP receive buffer.
socketTimeout               The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out.
storageEncryption           The S3 server side encryption to be used when saving objects on S3 (currently only AES256 is supported)
userAgent                   The HTTP user agent header passed with all HTTP requests.
uploadMaxThreads            The maximum number of threads used for multipart upload.
uploadChunkSize             The size of a single part in a multipart upload (default: `10 MB`).
uploadStorageClass          The S3 storage class applied to stored objects, one of [`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`] (default: `STANDARD`).
uploadMaxAttempts           The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).
uploadRetrySleep            The time to wait after a failed upload attempt to retry the part upload (default: `100ms`).
=========================== ================

For example::

    aws {
        client {
            maxConnections = 20
            connectionTimeout = 10000
            uploadStorageClass = 'INTELLIGENT_TIERING'
            storageEncryption = 'AES256'
        }
    }


.. _config-aws-batch:

Advanced Batch configuration options can be set by using the ``batch`` attribute. The following properties can be used (required version `19.07.0` or later):

=========================== ================
Name                        Description
=========================== ================
cliPath                     The path where the AWS command line tool is installed in the host AMI.
jobRole                     The AWS Job Role ARN that needs to be used to execute the Batch Job.
volumes                     One or more container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. ``/host/path:/mount/path[:ro|rw]``. Multiple mounts can be specifid separating them with a comma or using a list object.
delayBetweenAttempts        Delay between download attempts from S3 (default `10 sec`).
maxParallelTransfers        Max parallel upload/download transfer operations *per job* (default: ``4``).
maxTransferAttempts         Max number of downloads attempts from S3 (default: `1`).
=========================== ================

.. _config-cloud:

Scope `cloud`
-------------

.. note::
    The ``cloud`` configuration scope has been retired.


.. _config-conda:

Scope `conda`
-------------

The ``conda`` scope allows for the definition of the configuration settings that control the creation of a Conda environment
by the Conda package manager.

The following settings are available:

================== ================
Name                Description
================== ================
cacheDir            Defines the path where Conda environments are stored. When using a compute cluster make sure to provide a shared file system path accessible from all computing nodes.
createOptions       Defines any extra command line options supported by the ``conda create`` command. For details see: https://docs.conda.io/projects/conda/en/latest/commands/create.html.
createTimeout       Defines the amount of time the Conda environment creation can last. The creation process is terminated when the timeout is exceeded (default: ``20 min``).
useMamba            Uses the ``mamba`` binary instead of ``conda`` to create the conda environments. For details see: https://github.com/mamba-org/mamba.
================== ================


.. _config-k8s:

Scope `k8s`
-----------

The ``k8s`` scope allows the definition of the configuration settings that control the deployment and execution of
workflow applications in a Kubernetes cluster.

The following settings are available:

================== ================
Name                Description
================== ================
autoMountHostPaths  Automatically mounts host paths in the job pods. Only for development purpose when using a single node cluster (default: ``false``).
context             Defines the Kubernetes `configuration context name <https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/>`_ to use.
namespace           Defines the Kubernetes namespace to use (default: ``default``).
serviceAccount      Defines the Kubernetes `service account name <https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/>`_ to use.
launchDir           Defines the path where the workflow is launched and the user data is stored. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/<user-name>``.
workDir             Defines the path where the workflow temporary data is stored. This must be a path in a shared K8s persistent volume (default:``<user-dir>/work``).
projectDir          Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/projects``).
pod                 Allows the definition of one or more pod configuration options such as environment variables, config maps, secrets, etc. It allows the same settings as the :ref:`process-pod` process directive.
pullPolicy          Defines the strategy to be used to pull the container image e.g. ``pullPolicy: 'Always'``.
runAsUser           Defines the user ID to be used to run the containers. Shortcut for the ``securityContext`` option.
securityContext     Defines the `security context <https://kubernetes.io/docs/tasks/configure-pod-container/security-context/>`_ for all pods.
storageClaimName    The name of the persistent volume claim where store workflow result data.
storageMountPath    The path location used to mount the persistent volume claim (default: ``/workspace``).
storageSubPath      The path in the persistent volume to be mounted (default: root).
volumeClaims        (deprecated)
================== ================

See the :ref:`k8s-page` documentation for more details.

.. _config-timeline:

Scope `timeline`
----------------

The ``timeline`` scope allows you to enable/disable the processes execution timeline report generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the timeline report file (default: ``false``).
file                Timeline file name (default: ``timeline.html``).
overwrite           When ``true`` overwrites an existing timeline file instead of rolling it.
================== ================

.. _config-mail:

Scope `mail`
------------

The ``mail`` scope allows you to define the mail server configuration settings needed to send email messages.

================== ================
Name                Description
================== ================
from                Default email sender address.
smtp.host           Host name of the mail server.
smtp.port           Port number of the mail server.
smtp.user           User name to connect to  the mail server.
smtp.password       User password to connect to the mail server.
smtp.proxy.host     Host name of an HTTP web proxy server that will be used for connections to the mail server.
smtp.proxy.port     Port number for the HTTP web proxy server.
smtp.*              Any SMTP configuration property supported by the Java Mail API (see link below).
debug               When ``true`` enables Java Mail logging for debugging purpose.
================== ================

.. note:: Nextflow relies on the `Java Mail API <https://javaee.github.io/javamail/>`_ to send email messages.
  Advanced mail configuration can be provided by using any SMTP configuration property supported by the Java Mail API.
  See the `table of available properties at this link <https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties>`_.

For example, the following snippet shows how to configure Nextflow to send emails through the
`AWS Simple Email Service <https://aws.amazon.com/ses/>`_::

    mail {
        smtp.host = 'email-smtp.us-east-1.amazonaws.com'
        smtp.port = 587
        smtp.user = '<Your AWS SES access key>'
        smtp.password = '<Your AWS SES secret key>'
        smtp.auth = true
        smtp.starttls.enable = true
        smtp.starttls.required = true
    }

.. _config-notification:

Scope `notification`
--------------------

The ``notification`` scope allows you to define the automatic sending of a notification email message
when the workflow execution terminates.

================== ================
Name                Description
================== ================
enabled             Enables the sending of a notification message when the workflow execution completes.
to                  Recipient address for the notification email. Multiple addresses can be specified separating them with a comma.
from                Sender address for the notification email message.
template            Path of a template file which provides the content of the notification message.
binding             An associative array modelling the variables in the template file.
================== ================

The notification message is sent my using the STMP server defined in the configuration :ref:`mail scope<config-mail>`.

If no mail configuration is provided, it tries to send the notification message by using the external mail command
eventually provided by the underlying system (eg. ``sendmail`` or ``mail``).

.. _config-report:

Scope `report`
--------------

The ``report`` scope allows you to define configuration setting of the workflow :ref:`execution-report`.

================== ================
Name                Description
================== ================
enabled             If ``true`` it create the workflow execution report.
file                The path of the created execution report file (default: ``report.html``).
overwrite           When ``true`` overwrites existing report file instead of rolling it.
================== ================


.. _config-tower:

Scope `tower`
--------------

The ``tower`` configuration scope controls the settings for the `Nextflow Tower <https://tower.nf>`_ monitoring and tracing service.

The following settings are available:

================== ================
Name                Description
================== ================
enabled            When ``true`` Nextflow sends the workflow tracing and execution metrics to the Nextflow Tower service (default: ``false``).
accessToken        The unique access token specific to your account on an instance of Tower.
endpoint           The endpoint of your Tower deployment (default: ``https://tower.nf``).
workspaceId        The ID of the Tower workspace where the run should be added (default: the launching user personal workspace).
================== ================

The above options can be used by prefixing them with the ``tower`` scope or surrounding them by curly
brackets, as shown below::

    tower {
      enabled = true
      accessToken = '<YOUR TOKEN>'
      workspaceId = '<YOUR WORKSPACE ID>'
    }

.. tip::
  Your ``accessToken`` can be obtained from your Tower instance in the `Tokens page <https://tower.nf/tokens>`.

.. tip:: 
  The Tower workspace ID can also the specified using the environment variable ``TOWER_WORKSPACE_ID`` (config file has priority over the environment variable). 

.. _config-weblog:

Scope `weblog`
--------------

The ``weblog`` scope allows you to send detailed :ref:`trace scope<trace-fields>` information as HTTP POST request to a webserver, shipped as a JSON object.

Detailed information about the JSON fields can be found in the :ref:`weblog description<weblog-service>`.

================== ================
Name                Description
================== ================
enabled             If ``true`` it will send HTTP POST requests to a given url.
url                The url where to send HTTP POST requests (default: ``http:localhost``).
================== ================

.. _config-miscellaneous:

Miscellaneous
-------------

There are additional variables that can be defined within a configuration file that do not have a dedicated scope.

These are defined alongside other scopes, but the option is assigned as typically variable.

================== ================
Name                Description
================== ================
cleanup             If ``true``, on a successful completion of a run all files in *work* directory are automatically deleted.
================== ================

.. warning:: 
    The use of the above ``cleanup`` option will prevent the use of the *resume* feature on subsequent executions of that pipeline run. 
    Also, be aware that deleting all scratch files can take a lot of time especially when using shared file system or remote cloud storage.

.. _config-profiles:

Config profiles
===============

Configuration files can contain the definition of one or more *profiles*. A profile is a set of configuration attributes
that can be activated/chosen when launching a pipeline execution by using the ``-profile`` command line option.

Configuration profiles are defined by using the special scope ``profiles`` which group the attributes that belong
to the same profile using a common prefix. For example::

    profiles {

        standard {
            process.executor = 'local'
        }

        cluster {
            process.executor = 'sge'
            process.queue = 'long'
            process.memory = '10GB'
        }

        cloud {
            process.executor = 'cirrus'
            process.container = 'cbcrg/imagex'
            docker.enabled = true
        }

    }


This configuration defines three different profiles: ``standard``, ``cluster`` and ``cloud`` that set different process
configuration strategies depending on the target runtime platform. By convention the ``standard`` profile is implicitly used
when no other profile is specified by the user.

.. tip:: Two or more configuration profiles can be specified by separating the profile names
    with a comma character, for example::

        nextflow run <your script> -profile standard,cloud



.. danger:: When using the *profiles* feature in your config file do NOT set attributes in the same scope both
  inside and outside a ``profiles`` context. For example::

        process.cpus = 1

        profiles {
          foo {
            process.memory = '2 GB'
          }

          bar {
            process.memory = '4 GB'
          }
        }

  In the above example the ``process.cpus`` attribute is not correctly applied because the ``process`` scope is also
  used in the ``foo`` and ``bar`` profile contexts.

The above feature requires version 0.28.x or higher.

.. _config-env-vars:

Environment variables
=====================

The following environment variables control the configuration of the Nextflow runtime and
the Java virtual machine used by it.

=============================== ================
Name                            Description
=============================== ================
NXF_HOME                        Nextflow home directory (default: ``$HOME/.nextflow``).
NXF_VER                         Defines what version of Nextflow to use.
NXF_ORG                         Default `organization` prefix when looking for a hosted repository (default: ``nextflow-io``).
NXF_GRAB                        Provides extra runtime dependencies downloaded from a Maven repository service.
NXF_OPTS                        Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of ``-Dkey[=value]`` properties.
NXF_JVM_ARGS                    Allows the setting Java VM options. This is similar to ``NXF_OPTS`` however it's only applied the JVM running Nextflow and not to any java pre-launching commands (requires ``21.12.1-edge`` or later).
NXF_CLASSPATH                   Allows the extension of the Java runtime classpath with extra JAR files or class folders.
NXF_ASSETS                      Defines the directory where downloaded pipeline repositories are stored (default: ``$NXF_HOME/assets``)
NXF_PID_FILE                    Name of the file where the process PID is saved when Nextflow is launched in background.
NXF_WORK                        Directory where working files are stored (usually your *scratch* directory)
NXF_TEMP                        Directory where temporary files are stored
NXF_DEBUG                       Defines scripts debugging level: ``1`` dump task environment variables in the task log file; ``2`` enables command script execution tracing; ``3`` enables command wrapper execution tracing.
NXF_EXECUTOR                    Defines the default process executor e.g. `sge`
NXF_CONDA_CACHEDIR              Directory where Conda environments are store. When using a computing cluster it must be a shared folder accessible from all computing nodes.
NXF_SINGULARITY_CACHEDIR        Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all computing nodes.
NXF_SINGULARITY_LIBRARYDIR      Directory where remote Singularity images are retrieved. It should be a directory accessible to all computing nodes (requires: ``21.09.0-edge`` or later).
NXF_CHARLIECLOUD_CACHEDIR       Directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible from all computing nodes.
NXF_JAVA_HOME                   Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the ``JAVA_HOME`` variable if defined.
NXF_OFFLINE                     When ``true`` disables the project automatic download and update from remote repositories (default: ``false``).
NXF_CLOUD_DRIVER                Defines the default cloud driver to be used if not specified in the config file or as command line option, either ``aws`` or ``google``.
NXF_ANSI_LOG                    Enables/disables ANSI console output (default ``true`` when ANSI terminal is detected).
NXF_ANSI_SUMMARY                Enables/disables ANSI completion summary: `true|false` (default: print summary if execution last more than 1 minute).
NXF_SCM_FILE                    Defines the path location of the SCM config file (requires version ``20.10.0`` or later).
NXF_PARAMS_FILE                 Defines the path location of the pipeline parameters file (requires version ``20.10.0`` or later).
NXF_DISABLE_JOBS_CANCELLATION   Disables the cancellation of child jobs on workflow execution termination (requires ``21.12.0-edge`` or later).
JAVA_HOME                       Defines the path location of the Java VM installation used to run Nextflow.
JAVA_CMD                        Defines the path location of the Java binary command used to launch Nextflow.
HTTP_PROXY                      Defines the HTTP proxy server. As of version ``21.06.0-edge``, proxy authentication is supported providing the credentials in the proxy URL e.g. ``http://user:password@proxy-host.com:port``.
HTTPS_PROXY                     Defines the HTTPS proxy server. As of version ``21.06.0-edge``, proxy authentication is supported providing the credentials in the proxy URL e.g. ``https://user:password@proxy-host.com:port``.
FTP_PROXY                       Defines the FTP proxy server. Proxy authentication is supported providing the credentials in the proxy URL e.g. ``ftp://user:password@proxy-host.com:port``. FTP proxy support requires version ``21.06.0-edge`` or later.
NO_PROXY                        Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.
=============================== ================
.. _process-page:

*********
Processes
*********

In Nextflow a `process` is the basic processing `primitive` to execute a user script.

The process definition starts with the keyword ``process``, followed by process name and finally the process `body`
delimited by curly brackets. The process body must contain a string which represents the command or, more generally,
a script that is executed by it. A basic process looks like the following example::

  process sayHello {

      """
      echo 'Hello world!' > file
      """

  }


A process may contain five definition blocks, respectively: directives,
inputs, outputs, when clause and finally the process script. The syntax is defined as follows:

::

  process < name > {

     [ directives ]

     input:
      < process inputs >

     output:
      < process outputs >

     when:
      < condition >

     [script|shell|exec]:
     < user script to be executed >

  }


.. _process-script:

Script
======

The `script` block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains one and only one script block, and it must be the last statement when the process contains
input and output declarations.

The entered string is executed as a `Bash <http://en.wikipedia.org/wiki/Bash_(Unix_shell)>`_ script in the
`host` system. It can be any command, script or combination of them, that you would normally use in terminal shell
or in a common Bash script.

The only limitation to the commands that can be used in the script statement is given by the availability of those
programs in the target execution system.


The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts
composed by multiple commands spanning over multiple lines. For example::

    process doMoreThings {

      """
      blastp -db $db -query query.fa -outfmt 6 > blast_result
      cat blast_result | head -n 10 | cut -f 2 > top_hits
      blastdbcmd -db $db -entry_batch top_hits > sequences
      """

    }

As explained in the script tutorial section, strings can be defined by using a single-quote
or a double-quote, and multi-line strings are defined by three single-quote or three double-quote characters.

There is a subtle but important difference between them. Like in Bash, strings delimited by a ``"`` character support
variable substitutions, while strings delimited by ``'`` do not.

In the above code fragment the ``$db`` variable is replaced by the actual value defined somewhere in the
pipeline script.

.. warning:: Since Nextflow uses the same Bash syntax for variable substitutions in strings, you need to manage them
  carefully depending on if you want to evaluate a variable in the Nextflow context - or - in the Bash environment execution.

When you need to access a system environment variable  in your script you have two options. The first choice is as
easy as defining your script block by using a single-quote string. For example::

    process printPath {

       '''
       echo The path is: $PATH
       '''

    }

The drawback of this solution is that you will not able to access variables defined in the pipeline script context,
in your script block.

To fix this, define your script by using a double-quote string and `escape` the system environment variables by
prefixing them with a back-slash ``\`` character, as shown in the following example::


    process doOtherThings {

      """
      blastp -db \$DB -query query.fa -outfmt 6 > blast_result
      cat blast_result | head -n $MAX | cut -f 2 > top_hits
      blastdbcmd -db \$DB -entry_batch top_hits > sequences
      """

    }

In this example the ``$MAX`` variable has to be defined somewhere before, in the pipeline script.
`Nextflow` replaces it with the actual value before executing the script. Instead, the ``$DB`` variable
must exist in the script execution environment and the Bash interpreter will replace it with the actual value.

.. tip::
  Alternatively you can use the :ref:`process-shell` block definition which allows a script to contain both
  Bash and Nextflow variables without having to escape the first.

Scripts `à la carte`
--------------------

The process script is interpreted by Nextflow as a Bash script by default, but you are not limited to it.

You can use your favourite scripting language (e.g. Perl, Python, Ruby, R, etc), or even mix them in the same pipeline.

A pipeline may be composed by processes that execute very different tasks. Using `Nextflow` you can choose the scripting
language that better fits the task carried out by a specified process. For example for some processes `R` could be
more useful than `Perl`, in other you may need to use `Python` because it provides better access to a library or an API, etc.

To use a scripting other than Bash, simply start your process script with the corresponding
`shebang <http://en.wikipedia.org/wiki/Shebang_(Unix)>`_ declaration. For example::

    process perlStuff {

        """
        #!/usr/bin/perl

        print 'Hi there!' . '\n';
        """

    }

    process pyStuff {

        """
        #!/usr/bin/python

        x = 'Hello'
        y = 'world!'
        print "%s - %s" % (x,y)
        """

    }


.. tip:: Since the actual location of the interpreter binary file can change across platforms, to make your scripts
   more portable it is wise to use the ``env`` shell command followed by the interpreter's name, instead of the absolute
   path of it. Thus, the `shebang` declaration for a Perl script, for example,
   would look like: ``#!/usr/bin/env perl`` instead of the one in the above pipeline fragment.


Conditional scripts
-------------------

Complex process scripts may need to evaluate conditions on the input parameters or use traditional flow control
statements (i.e. ``if``, ``switch``, etc) in order to execute specific script commands, depending on the current
inputs configuration.

Process scripts can contain conditional statements by simply prefixing the script block with the keyword ``script:``.
By doing that the interpreter will evaluate all the following statements as a code block that must return the
script string to be executed. It's much easier to use than to explain, for example::


    seq_to_align = ...
    mode = 'tcoffee'

    process align {
        input:
        file seq_to_aln from sequences

        script:
        if( mode == 'tcoffee' )
            """
            t_coffee -in $seq_to_aln > out_file
            """

        else if( mode == 'mafft' )
            """
            mafft --anysymbol --parttree --quiet $seq_to_aln > out_file
            """

        else if( mode == 'clustalo' )
            """
            clustalo -i $seq_to_aln -o out_file
            """

        else
            error "Invalid alignment mode: ${mode}"

    }


In the above example the process will execute the script fragment depending on the value of the ``mode`` parameter.
By default it will execute the ``tcoffee`` command, changing the ``mode`` variable to ``mafft`` or ``clustalo`` value,
the other branches will be executed.

.. _process-template:

Template
--------

Process script can be externalised by using *template* files which can be reused across different processes and tested
independently from the overall pipeline execution.

A template is simply a shell script file that Nextflow is able to execute by using the ``template`` function
as shown below::

    process template_example {

        input:
        val STR from 'this', 'that'

        script:
        template 'my_script.sh'

    }


Nextflow looks for the ``my_script.sh`` template file in the directory ``templates`` that must exist in the same folder
where the Nextflow script file is located (any other location can be provided by using an absolute template path).

.. note::
  When using :ref:`DSL2 <dsl2-page>` Nextflow looks for the specified file name also in the ``templates`` directory
  located in the same folder where the module script is placed. See :ref:`module templates <module-templates>`.


The template script can contain any piece of code that can be executed by the underlying system. For example::

  #!/bin/bash
  echo "process started at `date`"
  echo $STR
  :
  echo "process completed"



.. tip::
  Note that the dollar character (``$``) is interpreted as a Nextflow variable placeholder, when the script is run as a
  Nextflow template, while it is evaluated as a Bash variable when it is run alone. This can be very useful to test
  your script autonomously, i.e. independently from Nextflow execution. You only need to provide a Bash environment
  variable for each the Nextflow variable existing in your script. For example, it would be possible to execute the above
  script entering the following command in the shell terminal: ``STR='foo' bash templates/my_script.sh``


.. _process-shell:

Shell
-----

The ``shell`` block is a string statement that defines the *shell* command executed by the process to carry out its task.
It is an alternative to the :ref:`process-script` definition with an important difference, it uses
the exclamation mark ``!`` character as the variable placeholder for Nextflow variables in place of the usual dollar character.

In this way it is possible to use both Nextflow and Bash variables in the same piece of code without having to escape
the latter and making process scripts more readable and easy to maintain. For example::

    process myTask {

        input:
        val str from 'Hello', 'Hola', 'Bonjour'

        shell:
        '''
        echo User $USER says !{str}
        '''

    }



In the above trivial example the ``$USER`` variable is managed by the Bash interpreter, while ``!{str}`` is handled
as a process input variable managed by Nextflow.

.. note::

    - Shell script definition requires the use of single-quote ``'`` delimited strings. When using double-quote ``"``
      delimited strings, dollar variables are interpreted as Nextflow variables as usual. See :ref:`string-interpolation`.

    - Exclamation mark prefixed variables always need to be enclosed in curly brackets i.e. ``!{str}`` is a valid 
      variable while ``!str`` is ignored.

    - Shell script supports the use of the file :ref:`process-template` mechanism. The same rules are applied to the variables
      defined in the script template.

.. _process-native:

Native execution
----------------

Nextflow processes can execute native code other than system scripts as shown in the previous paragraphs.

This means that instead of specifying the process command to be executed as a string script, you can
define it by providing one or more language statements, as you would do in the rest of the pipeline script.
Simply starting the script definition block with the ``exec:`` keyword, for example::

    x = Channel.from( 'a', 'b', 'c')

    process simpleSum {
        input:
        val x

        exec:
        println "Hello Mr. $x"
    }

Will display::

    Hello Mr. b
    Hello Mr. a
    Hello Mr. c


.. _process-stub:

Stub
====

.. warning::
    This is an incubating feature. It may change in future versions.

As of version 20.11.0-edge it's possible to define a command *stub* that replaces the actual process command, when
the `-stub-run` or `-stub` command line option. ::

    process INDEX {
        input:
          path transcriptome
        output:
          path 'index'

        script:
          """
          salmon index --threads $task.cpus -t $transcriptome -i index
          """

        stub:
          """
          mkdir index
          touch index/seq.bin
          touch index/info.json
          touch index/refseq.bin
          """
    }

This feature is meant to allow the fast prototyping and test of the workflow logic without using the real
commands. The developer can use it to provide a dummy command which is expected to mimic the execution
of the real one in a quicker manner. This can also be used as an alternative for the *dry-run* feature.

.. tip::
    The ``stub`` block can be defined before or after the process ``script`` definition.
    When the execution is run with the option `-stub-run` and a process is not implementing the ``stub`` command the
    real is executed.


.. _process-input:

Inputs
======

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

The `input` block defines from which channels the process expects to receive data. You can only define one
input block at a time and it must contain one or more input declarations.

The input block follows the syntax shown below::

    input:
      <input qualifier> <input name> [from <source channel>] [attributes]


An input definition starts with an input `qualifier` and the input `name`, followed by the keyword ``from`` and
the actual channel over which inputs are received. Finally some input optional attributes can be specified.

.. note:: When the input name is the same as the channel name, the ``from`` part of the declaration can be omitted.

The input qualifier declares the `type` of data to be received. This information is used by Nextflow to apply the
semantic rules associated to each qualifier and handle it properly depending on the target execution platform
(grid, cloud, etc).

The qualifiers available are the ones listed in the following table:

=========== =============
Qualifier   Semantic
=========== =============
val         Lets you access the received input value by its name in the process script.
env         Lets you use the received value to set an environment variable named
            as the specified input name.
file        Lets you handle the received value as a file, staging it properly in the execution context.
path        Lets you handle the received value as a path, staging the file properly in the execution context.
stdin       Lets you forward the received value to the process `stdin` special file.
tuple       Lets you handle a group of input values having one of the above qualifiers.
each        Lets you execute the process for each entry in the input collection.
=========== =============


Input of generic values
-----------------------

The ``val`` qualifier allows you to receive data of any type as input. It can be accessed in the process script
by using the specified input name, as shown in the following example::

    num = Channel.from( 1, 2, 3 )

    process basicExample {
      input:
      val x from num

      "echo process job $x"

    }


In the above example the process is executed three times, each time a value is received from the channel ``num``
and used to process the script. Thus, it results in an output similar to the one shown below::

    process job 3
    process job 1
    process job 2

.. note:: The `channel` guarantees that items are delivered in the same order as they have been sent - but -
  since the process is executed in a parallel manner, there is no guarantee that they are processed in the
  same order as they are received. In fact, in the above example, value ``3`` is processed before the others.


When the ``val`` has the same name as the channel from where the data is received, the ``from`` part can be omitted.
Thus the above example can be written as shown below::

    num = Channel.from( 1, 2, 3 )

    process basicExample {
      input:
      val num

      "echo process job $num"

    }


Input of files
--------------

The ``file`` qualifier allows the handling of file values in the process execution context. This means that
Nextflow will stage it in the process execution directory, and it can be access in the script by using the name
specified in the input declaration. For example::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file query_file from proteins

      "blastp -query ${query_file} -db nr"

    }

In the above example all the files ending with the suffix ``.fa`` are sent over the channel ``proteins``.
Then, these files are received by the process which will execute a `BLAST` query on each of them.

When the file input name is the same as the channel name, the ``from`` part of the input declaration can be omitted.
Thus, the above example could be written as shown below::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file proteins

      "blastp -query $proteins -db nr"

    }


It's worth noting that in the above examples, the name of the file in the file-system is not touched, you can
access the file even without knowing its name because you can reference it in the process script by using the
variable whose name is specified in the input file parameter declaration.

There may be cases where your task needs to use a file whose name is fixed, it does not have to change along
with the actual provided file. In this case you can specify its name by specifying the ``name`` attribute in the
input file parameter declaration, as shown in the following example::

    input:
        file query_file name 'query.fa' from proteins


Or alternatively using a shorter syntax::

    input:
        file 'query.fa' from proteins


Using this, the previous example can be re-written as shown below::

    proteins = Channel.fromPath( '/some/path/*.fa' )

    process blastThemAll {
      input:
      file 'query.fa' from proteins

      "blastp -query query.fa -db nr"

    }


What happens in this example is that each file, that the process receives, is staged with the name ``query.fa``
in a different execution context (i.e. the folder where the job is executed) and an independent process
execution is launched.

.. tip:: This allows you to execute the process command various time without worrying the files names changing.
  In other words, `Nextflow` helps you write pipeline tasks that are self-contained and decoupled by the execution
  environment. This is also the reason why you should avoid whenever possible to use absolute or relative paths
  referencing files in your pipeline processes.


.. TODO describe that file can handle channels containing any data type not only file


Multiple input files
--------------------

A process can declare as input file a channel that emits a collection of values, instead of a simple value.

In this case, the script variable defined by the input file parameter will hold a list of files. You can
use it as shown before, referring to all the files in the list, or by accessing a specific entry using the
usual square brackets notation.

When a target file name is defined in the input parameter and a collection of files is received by the process,
the file name will be appended by a numerical suffix representing its ordinal position in the list. For example::

    fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size:3)

    process blastThemAll {
        input:
        file 'seq' from fasta

        "echo seq*"

    }

Will output::

    seq1 seq2 seq3
    seq1 seq2 seq3
    ...

The target input file name can contain the ``*`` and ``?`` wildcards, that can be used
to control the name of staged files. The following table shows how the wildcards are
replaced depending on the cardinality of the received input collection.

============ ============== ==================================================
Cardinality   Name pattern     Staged file names
============ ============== ==================================================
 any         ``*``           named as the source file
 1           ``file*.ext``   ``file.ext``
 1           ``file?.ext``   ``file1.ext``
 1           ``file??.ext``  ``file01.ext``
 many        ``file*.ext``   ``file1.ext``, ``file2.ext``, ``file3.ext``, ..
 many        ``file?.ext``   ``file1.ext``, ``file2.ext``, ``file3.ext``, ..
 many        ``file??.ext``  ``file01.ext``, ``file02.ext``, ``file03.ext``, ..
 many        ``dir/*``       named as the source file, created in ``dir`` subdirectory
 many        ``dir??/*``     named as the source file, created in a progressively indexed subdirectory e.g. ``dir01/``, ``dir02/``, etc.
 many        ``dir*/*``      (as above)
============ ============== ==================================================

The following fragment shows how a wildcard can be used in the input file declaration::


    fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size:3)

    process blastThemAll {
        input:
        file 'seq?.fa' from fasta

        "cat seq1.fa seq2.fa seq3.fa"

    }


.. note:: Rewriting input file names according to a named pattern is an extra feature and not at all obligatory.
  The normal file input constructs introduced in the `Input of files`_ section are valid for collections of
  multiple files as well. To handle multiple input files preserving the original file names, use the ``*`` wildcard as
  name pattern or a variable identifier.

Dynamic input file names
------------------------

When the input file name is specified by using the ``name`` file clause or the short `string` notation, you
are allowed to use other input values as variables in the file name string. For example::


  process simpleCount {
    input:
    val x from species
    file "${x}.fa" from genomes

    """
    cat ${x}.fa | grep '>'
    """
  }


In the above example, the input file name is set by using the current value of the ``x`` input value.

This allows the input files to be staged in the script working directory with a name that is coherent
with the current execution context.

.. tip:: In most cases, you won't need to use dynamic file names, because each process is executed in its 
  own private temporary directory, and input files are automatically staged to this directory by Nextflow. 
  This guarantees that input files with the same name won't overwrite each other.


Input of type 'path'
--------------------

The ``path`` input qualifier was introduced by Nextflow version 19.10.0 and it's a drop-in replacement
for the ``file`` qualifier, therefore it's backward compatible with the syntax
and the semantic for the input ``file`` described above.

The important difference between ``file`` and ``path`` qualifier is that the first expects the
values received as input to be *file* objects. When inputs is a different type, it automatically
coverts to a string and saves it to a temporary files. This can be useful in some uses cases,
but it turned out to be tricky in most common cases.

The ``path`` qualifier instead interprets string values as the path location of the input file
and automatically converts to a file object.

::

    process foo {
      input:
        path x from '/some/data/file.txt'
      """
        your_command --in $x
      """
    }


.. note::
    Provided input value should represent an absolute path location i.e. the string value
    **must** be prefixed with a `/` character or with a supported URI protocol i.e. ``file://``,
    ``http://``, ``s3://``, etc. and it cannot contains special characters (e.g. ``\n``, etc.).



The option ``stageAs`` allow you to control how the file should be named in the task work
directory, providing a specific name or a name pattern as described in the `Multiple input files`_
section::


    process foo {
      input:
        path x, stageAs: 'data.txt' from '/some/data/file.txt'
      """
        your_command --in data.txt
      """
    }


.. tip::
    The ``path`` qualifier should be preferred over ``file`` to handle process input files
    when using Nextflow 19.10.0 or later.


Input of type 'stdin'
---------------------

The ``stdin`` input qualifier allows you the forwarding of the value received from a channel to the
`standard input <http://en.wikipedia.org/wiki/Standard_streams#Standard_input_.28stdin.29>`_
of the command executed by the process. For example::

    str = Channel.from('hello', 'hola', 'bonjour', 'ciao').map { it+'\n' }

    process printAll {
       input:
       stdin str

       """
       cat -
       """

    }

It will output::

    hola
    bonjour
    ciao
    hello




Input of type 'env'
-------------------

The ``env`` qualifier allows you to define an environment variable in the process execution context based
on the value received from the channel. For example::

    str = Channel.from('hello', 'hola', 'bonjour', 'ciao')

    process printEnv {

        input:
        env HELLO from str

        '''
        echo $HELLO world!
        '''

    }

::

    hello world!
    ciao world!
    bonjour world!
    hola world!


.. _process-input-set:

Input of type 'set'
-------------------

.. warning:: The `set` input type has been deprecated. See `tuple` instead.


.. _process-input-tuple:

Input of type 'tuple'
---------------------


The ``tuple`` qualifier allows you to group multiple parameters in a single parameter definition. It can be useful
when a process receives, in input, tuples of values that need to be handled separately. Each element in the tuple
is associated to a corresponding element with the ``tuple`` definition. For example::

     values = Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

     process tupleExample {
         input:
         tuple val(x), file('latin.txt') from values

         """
         echo Processing $x
         cat - latin.txt > copy
         """

     }


In the above example the ``tuple`` parameter is used to define the value ``x`` and the file ``latin.txt``,
which will receive a value from the same channel.

In the ``tuple`` declaration items can be defined by using the following qualifiers: ``val``, ``env``, ``file`` and ``stdin``.

A shorter notation can be used by applying the following substitution rules:

============== =======
long            short
============== =======
val(x)          x
file(x)         (not supported)
file('name')    'name'
file(x:'name')  x:'name'
stdin           '-'
env(x)          (not supported)
============== =======

Thus the previous example could be rewritten as follows::

      values = Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

      process tupleExample {
          input:
          tuple x, 'latin.txt' from values

          """
          echo Processing $x
          cat - latin.txt > copy
          """
      }

File names can be defined in *dynamic* manner as explained in the `Dynamic input file names`_ section.


Input repeaters
---------------

The ``each`` qualifier allows you to repeat the execution of a process for each item in a collection,
every time a new data is received. For example::

  sequences = Channel.fromPath('*.fa')
  methods = ['regular', 'expresso', 'psicoffee']

  process alignSequences {
    input:
    file seq from sequences
    each mode from methods

    """
    t_coffee -in $seq -mode $mode > result
    """
  }


In the above example every time a file of sequences is received as input by the process,
it executes *three* tasks running a T-coffee alignment with a different value for the ``mode`` parameter.
This is useful when you need to `repeat` the same task for a given set of parameters.

Since version 0.25+ input repeaters can be applied to files as well. For example::

    sequences = Channel.fromPath('*.fa')
    methods = ['regular', 'expresso']
    libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]

    process alignSequences {
      input:
      file seq from sequences
      each mode from methods
      each file(lib) from libraries

      """
      t_coffee -in $seq -mode $mode -lib $lib > result
      """
    }


.. note:: When multiple repeaters are declared, the process is executed for each *combination* of them.

In the latter example for any sequence input file emitted by the ``sequences`` channel are executed 6 alignments,
3 using the ``regular`` method against each library files, and other 3 by using the ``expresso`` method always
against the same library files.


.. hint:: If you need to repeat the execution of a process over n-tuple of elements instead a simple values or files,
  create a channel combining the input values as needed to trigger the process execution multiple times.
  In this regard, see the :ref:`operator-combine`, :ref:`operator-cross` and :ref:`operator-phase` operators.

.. _process-understand-how-multiple-input-channels-work:

Understand how multiple input channels work
-------------------------------------------

A key feature of processes is the ability to handle inputs from multiple channels.

When two or more channels are declared as process inputs, the process stops until
there's a complete input configuration ie. it receives an input value from all the channels declared
as input.

When this condition is verified, it consumes the input values coming from the respective channels,
and spawns a task execution, then repeat the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel
cause the process execution to stop even if there are other values in other channels.

For example::

  process foo {
    echo true
    input:
    val x from Channel.from(1,2)
    val y from Channel.from('a','b','c')
    script:
     """
     echo $x and $y
     """
  }


The process ``foo`` is executed two times because the first input channel only provides two values and therefore
the ``c`` element is discarded. It prints::

    1 and a
    2 and b


.. warning:: A different semantic is applied when using *Value channel* a.k.a. *Singleton channel*.

This kind of channel is created by the :ref:`Channel.value <channel-value>` factory method or implicitly
when a process input specifies a simple value in the ``from`` clause.

By definition, a *Value channel* is bound to a single value and it can be read unlimited times without
consuming its content.

These properties make that when mixing a *value channel* with one or more (queue) channels,
it does not affect the process termination which only depends by the other channels and its
content is applied repeatedly.

To better understand this behavior compare the previous example with the following one::

  process bar {
    echo true
    input:
    val x from Channel.value(1)
    val y from Channel.from('a','b','c')
    script:
     """
     echo $x and $y
     """
  }

The above snippet executes the ``bar`` process three times because the first input is a *value channel*, therefore
its content can be read as many times as needed. The process termination is determined by the content of the second
channel. It prints::


  1 and a
  1 and b
  1 and c

See also: :ref:`channel-types`.

Outputs
=======

The `output` declaration block allows you to define the channels used by the process to send out the results produced.
You can only define one output block at a time and it must contain one or more output declarations.

The output block follows the syntax shown below::
    output:
      <output qualifier> <output name> [into <target channel>[,channel,..]] [attribute [,..]]

Output definitions start by an output `qualifier` and the output `name`, followed by the keyword ``into`` and
one or more channels over which outputs are sent. Finally some optional attributes can be specified.

.. note:: When the output name is the same as the channel name, the ``into`` part of the declaration can be omitted.


.. TODO the channel is implicitly created if does not exist

The qualifiers that can be used in the output declaration block are the ones listed in the following table:

=========== =============
Qualifier   Semantic
=========== =============
val         Sends variables with the name specified over the output channel.
file        Sends a file produced by the process with the name specified over the output channel.
path        Sends a file produced by the process with the name specified over the output channel (replaces ``file``).
env         Sends the variable defined in the process environment with the name specified over the output channel.
stdout      Sends the executed process `stdout` over the output channel.
tuple       Sends multiple values over the same output channel.
=========== =============


Output values
-------------

The ``val`` qualifier allows you to output a `value` defined in the script context. In a common usage scenario,
this is a value which has been defined in the `input` declaration block, as shown in the following example::

   methods = ['prot','dna', 'rna']

   process foo {
     input:
     val x from methods

     output:
     val x into receiver

     """
     echo $x > file
     """

   }

   receiver.view { "Received: $it" }


Valid output values are value literals, input value identifiers, variables accessible in the process scope and
value expressions. For example::

    process foo {
      input:
      file fasta from 'dummy'

      output:
      val x into var_channel
      val 'BB11' into str_channel
      val "${fasta.baseName}.out" into exp_channel

      script:
      x = fasta.name
      """
      cat $x > file
      """
    }




Output files
------------

The ``file`` qualifier allows you to output one or more files, produced by the process, over the specified channel.
For example::


    process randomNum {

       output:
       file 'result.txt' into numbers

       '''
       echo $RANDOM > result.txt
       '''

    }

    numbers.subscribe { println "Received: " + it.text }


In the above example the process, when executed, creates a file named ``result.txt`` containing a random number.
Since a file parameter using the same name is declared between the outputs, when the task is completed that
file is sent over the ``numbers`` channel. A downstream `process` declaring the same channel as `input` will
be able to receive it.

.. note:: If the channel specified as output has not been previously declared in the pipeline script, it
  will implicitly be created by the output declaration itself.


.. TODO explain Path object

Multiple output files
---------------------

When an output file name contains a ``*`` or ``?`` wildcard character it is interpreted as a `glob`_ path matcher.
This allows you to *capture* multiple files into a list object and output them as a sole emission. For example::

    process splitLetters {

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }

    letters
        .flatMap()
        .subscribe { println "File: ${it.name} => ${it.text}" }

It prints::

    File: chunk_aa => H
    File: chunk_ab => o
    File: chunk_ac => l
    File: chunk_ad => a

.. note:: In the above example the operator :ref:`operator-flatmap` is used to transform the list of files emitted by
  the ``letters`` channel into a channel that emits each file object independently.

Some caveats on glob pattern behavior:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directory paths.
* When a two stars pattern ``**`` is used to recourse across directories, only file paths are matched
  i.e. directories are not included in the result list.

.. warning:: Although the input files matching a glob output declaration are not included in the
   resulting output channel, these files may still be transferred from the task scratch directory
   to the target task work directory. Therefore, to avoid unnecessary file copies it is recommended
   to avoid the usage of loose wildcards when defining output files e.g. ``file '*'`` .
   Instead, use a prefix or a postfix naming notation to restrict the set of matching files to
   only the expected ones e.g. ``file 'prefix_*.sorted.bam'``. 

By default all the files matching the specified glob pattern are emitted by the channel as a sole (list) item.
It is also possible to emit each file as a sole item by adding the ``mode flatten`` attribute in the output file
declaration.

By using the ``mode`` attribute the previous example can be re-written as shown below::

    process splitLetters {

        output:
        file 'chunk_*' into letters mode flatten

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }

    letters .subscribe { println "File: ${it.name} => ${it.text}" }


.. warning::
    The option ``mode`` is deprecated as of version 19.10.0. Use the operator :ref:`operator-collect`
    in the downstream process instead.

Read more about glob syntax at the following link `What is a glob?`_

.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
.. _What is a glob?: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob

.. _process-dynoutname:

Dynamic output file names
-------------------------

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated
string which references values defined in the input declaration block or in the script global context.
For example::


  process align {
    input:
    val x from species
    file seq from sequences

    output:
    file "${x}.aln" into genomes

    """
    t_coffee -in $seq > ${x}.aln
    """
  }

In the above example, each time the process is executed an alignment file is produced whose name depends
on the actual value of the ``x`` input.

.. tip:: The management of output files is a very common misunderstanding when using Nextflow.
  With other tools it is generally necessary to organize the output files into some kind of directory 
  structure or to guarantee a unique file name scheme, so that result files won't overwrite each other 
  and that they can be referenced univocally by downstream tasks.

  With Nextflow, in most cases, you don't need to take care of naming output files, because each task is executed 
  in its own unique temporary directory, so files produced by different tasks can never override each other.
  Also meta-data can be associated with outputs by using the :ref:`tuple output <process-out-tuple>` qualifier, instead of
  including them in the output file name.

  To sum up, the use of output files with static names over dynamic ones is preferable whenever possible, 
  because it will result in a simpler and more portable code.

.. _process-out-path:

Output path
-----------

The ``path`` output qualifier was introduced by Nextflow version 19.10.0 and it's a drop-in replacement
for the ``file`` output qualifier, therefore it's backward compatible with the syntax
and the semantic for the input ``file`` described above.

The main advantage of ``path`` over the ``file`` qualifier is that it allows the specification
of a number of outputs to fine-control the output files.

============== =====================
Name            Description
============== =====================
glob            When ``true`` the specified name is interpreted as a glob pattern (default: ``true``)
hidden          When ``true`` hidden files are included in the matching output files (default: ``false``)
followLinks     When ``true`` target files are return in place of any matching symlink (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``any``, or ``file`` if the specified file name pattern contains a `**` - double star - symbol)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
includeInputs   When ``true`` any input files matching an output file glob pattern are included.
============== =====================


.. warning::
    Breaking change: the ``file`` qualifier interprets ``:`` as path separator, therefore ``file 'foo:bar'``
    captures both files ``foo`` and ``bar``. The ``path`` qualifier interprets it as just a plain file name character,
    and therefore the output definition ``path 'foo:bar'`` captures the output file with name ``foo:bar``.


.. tip::
    The ``path`` qualifier should be preferred over ``file`` to handle process output files
    when using Nextflow 19.10.0 or later.

.. _process-stdout:

Output 'stdout' special file
----------------------------

The ``stdout`` qualifier allows you to `capture` the `stdout` output of the executed process and send it over
the channel specified in the output parameter declaration. For example::

    process sayHello {
        output:
        stdout ch

        """
        echo Hello world!
        """
    }

    ch.view { print "I say..  $it" }

In the above example ``ch`` represents an arbitrary channel variable that holds the process outputs.

.. _process-env:

Output 'env'
------------

The ``env`` qualifier allows you to capture a variable defined in the process execution environment
and send it over the channel specified in the output parameter declaration::

    process myTask {
        output:
        env FOO into target
        script:
        '''
        FOO=$(ls -la)
        '''
    }

    target.view { "directory content: $it" }


.. _process-set:

Output 'set' of values
----------------------

.. warning:: The `set` output type has been deprecated. See `tuple` instead.


.. _process-out-tuple:

Output 'tuple' of values
------------------------

The ``tuple`` qualifier allows you to send multiple values into a single channel. This feature is useful
when you need to `group together` the results of multiple executions of the same process, as shown in the following
example::

    query_ch = Channel.fromPath '*.fa'
    species_ch = Channel.from 'human', 'cow', 'horse'

    process blast {

    input:
      val species from query_ch
      file query from species_ch

    output:
      tuple val(species), file('result') into blastOuts

    script:
      """
      blast -db nr -query $query > result
      """
    }


In the above example a `BLAST` task is executed for each pair of ``species`` and ``query`` that are received.
When the task completes a new tuple containing the value for ``species`` and the file ``result`` is sent to the ``blastOuts`` channel.


A `tuple` declaration can contain any combination of the following qualifiers, previously described: ``val``, ``file`` and ``stdout``.

.. tip:: Variable identifiers are interpreted as `values` while strings literals are interpreted as `files` by default,
  thus the above output `tuple` can be rewritten using a short notation as shown below.

::

    output:
        tuple species, 'result' into blastOuts



File names can be defined in a dynamic manner as explained in the :ref:`process-dynoutname` section.

Optional Output
---------------

In most cases a process is expected to generate output that is added to the output channel.  However, there are situations where it is valid for a process to `not` generate output. In these cases ``optional true`` may be added to the output declaration, which tells Nextflow not to fail the process if the declared output is not created.

::

    output:
        file("output.txt") optional true into outChannel

In this example, the process is normally expected to generate an ``output.txt`` file, but in the cases where the file is legitimately missing, the process does not fail. ``outChannel`` is only populated by those processes that do generate ``output.txt``. 


When
====

The ``when`` declaration allows you to define a condition that must be verified in order to execute the process.
This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters. For example::


    process find {
      input:
      file proteins
      val type from dbtype

      when:
      proteins.name =~ /^BB11.*/ && type == 'nr'

      script:
      """
      blastp -query $proteins -db nr
      """

    }


.. _process-directives:

Directives
==========

Using the `directive` declarations block you can provide optional settings that will affect the execution of the current
process.

They must be entered at the top of the process `body`, before any other declaration blocks (i.e. ``input``, ``output``, etc) 
and have the following syntax::

    name value [, value2 [,..]]

Some directives are generally available to all processes, some others depends on the `executor` currently defined.

The directives are:

* `accelerator`_
* `afterScript`_
* `beforeScript`_
* `cache`_
* `cpus`_
* `conda`_
* `container`_
* `containerOptions`_
* `clusterOptions`_
* `disk`_
* `echo`_
* `errorStrategy`_
* `executor`_
* `ext`_
* `label`_
* `machineType`_
* `maxErrors`_
* `maxForks`_
* `maxRetries`_
* `memory`_
* `module`_
* `penv`_
* `pod`_
* `publishDir`_
* `queue`_
* `scratch`_
* `stageInMode`_
* `stageOutMode`_
* `storeDir`_
* `tag`_
* `time`_


accelerator
-----------

The ``accelerator`` directive allows you to specify the hardware accelerator requirement for the task execution
e.g. *GPU* processor. For example::

    process foo {
        accelerator 4, type: 'nvidia-tesla-k80'

        script:
        """
        your_gpu_enabled --command --line
        """
    }


The above examples will request 4 GPUs of type `nvidia-tesla-k80`.


.. note:: This directive is only supported by :ref:`awsbatch-executor`, :ref:`google-lifesciences-executor` and :ref:`k8s-executor` executors.

.. tip:: The accelerator ``type`` option value depends by the target execution platform. Refer to the target
  platform documentation for details on the available accelerators. `AWS <https://aws.amazon.com/batch/faqs/?#GPU_Scheduling_>`_
  `Google <https://cloud.google.com/compute/docs/gpus/>`_
  `Kubernetes <https://kubernetes.io/docs/tasks/manage-gpus/scheduling-gpus/#clusters-containing-different-types-of-gpus>`_.


afterScript
-----------

The ``afterScript`` directive allows you to execute a custom (Bash) snippet immediately *after* the main process has run.
This may be useful to clean up your staging area.

beforeScript
------------

The ``beforeScript`` directive allows you to execute a custom (Bash) snippet *before* the main process script is run.
This may be useful to initialise the underlying cluster environment or for other custom initialisation.

For example::

    process foo {

      beforeScript 'source /cluster/bin/setup'

      """
      echo bar
      """

    }

.. _process-cache:

cache
-----

The ``cache`` directive allows you to store the process results to a local cache. When the cache is enabled *and*
the pipeline is launched with the :ref:`resume <getstart-resume>` option, any following attempt to execute the process,
along with the same inputs, will cause the process execution to be skipped, producing the stored data as
the actual results.

The caching feature generates a unique `key` by indexing the process script and inputs. This key is used
to identify univocally the outputs produced by the process execution.


The cache is enabled by default, you can disable it for a specific process by setting the ``cache``
directive to ``false``. For example:: 

  process noCacheThis {
    cache false

    script:
    <your command string here>
  }

The ``cache`` directive possible values are shown in the following table:

===================== =================
Value                 Description
===================== =================
``false``             Disable cache feature.
``true`` (default)    Enable caching. Cache keys are created indexing input files meta-data information (name, size and last update timestamp attributes).
``'deep'``            Enable caching. Cache keys are created indexing input files content.
``'lenient'``         Enable caching. Cache keys are created indexing input files path and size attributes (this policy provides a workaround for incorrect caching invalidation observed on shared file systems due to inconsistent files timestamps; requires version 0.32.x or later).
===================== =================


.. _process-conda:

conda
-----

The ``conda`` directive allows for the definition of the process dependencies using the `Conda <https://conda.io>`_
package manager.

Nextflow automatically sets up an environment for the given package names listed by in the ``conda`` directive.
For example::

  process foo {
    conda 'bwa=0.7.15'

    '''
    your_command --here
    '''
  }


Multiple packages can be specified separating them with a blank space eg. ``bwa=0.7.15 fastqc=0.11.5``.
The name of the channel from where a specific package needs to be downloaded can be specified using the usual
Conda notation i.e. prefixing the package with the channel name as shown here ``bioconda::bwa=0.7.15``.

The ``conda`` directory also allows the specification of a Conda environment file
path or the path of an existing environment directory. See the :ref:`conda-page` page for further details.

.. _process-container:

container
---------

The ``container`` directive allows you to execute the process script in a `Docker <http://docker.io>`_ container.

It requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the
*local* executor or the cluster nodes when the pipeline is deployed through a *grid* executor.

For example::


    process runThisInDocker {

      container 'dockerbox:tag'

      """
      <your holy script here>
      """

    }


Simply replace in the above script ``dockerbox:tag`` with the Docker image name you want to use.

.. tip:: This can be very useful to execute your scripts into a replicable self-contained environment or to deploy your pipeline in the cloud.

.. note:: This directive is ignored for processes :ref:`executed natively <process-native>`.


.. _process-containerOptions:

containerOptions
----------------

The ``containerOptions`` directive allows you to specify any container execution option supported by the underlying
container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings
only for a specific process e.g. mount a custom path::


  process runThisWithDocker {

      container 'busybox:latest'
      containerOptions '--volume /data/db:/db'

      output: file 'output.txt'

      '''
      your_command --data /db > output.txt
      '''
  }


.. warning:: This feature is not supported by :ref:`k8s-executor` and :ref:`azurebatch-executor` executors.

.. _process-cpus:

cpus
----

The ``cpus`` directive allows you to define the number of (logical) CPU required by the process' task.
For example::

    process big_job {

      cpus 8
      executor 'sge'

      """
      blastp -query input_sequence -num_threads ${task.cpus}
      """
    }


This directive is required for tasks that execute multi-process or multi-threaded commands/tools and it is meant
to reserve enough CPUs when a pipeline task is executed through a cluster resource manager.

See also: `penv`_, `memory`_, `time`_, `queue`_, `maxForks`_

.. _process-clusterOptions:

clusterOptions
--------------

The ``clusterOptions`` directive allows the usage of any `native` configuration option accepted by your cluster submit command.
You can use it to request non-standard resources or use settings that are specific to your cluster and not supported
out of the box by Nextflow.

.. note:: This directive is taken in account only when using a grid based executor:
  :ref:`sge-executor`, :ref:`lsf-executor`, :ref:`slurm-executor`, :ref:`pbs-executor`, :ref:`pbspro-executor`,
  :ref:`moab-executor` and :ref:`condor-executor` executors.

.. _process-disk:

disk
----

The ``disk`` directive allows you to define how much local disk storage the process is allowed to use. For example::

    process big_job {

        disk '2 GB'
        executor 'cirrus'

        """
        your task script here
        """
    }

The following memory unit suffix can be used when specifying the disk value:

======= =============
Unit    Description
======= =============
B       Bytes
KB      Kilobytes
MB      Megabytes
GB      Gigabytes
TB      Terabytes
======= =============

.. note:: This directive currently is taken in account only by the :ref:`ignite-executor`
  and the :ref:`condor-executor` executors.

See also: `cpus`_, `memory`_ `time`_, `queue`_ and `Dynamic computing resources`_.

.. _process-echo:

echo
----

By default the `stdout` produced by the commands executed in all processes is ignored.
Setting the ``echo`` directive to ``true`` you can forward the process `stdout` to the current top
running process `stdout` file, showing it in the shell terminal.

For example::

    process sayHello {
      echo true

      script:
      "echo Hello"
    }

::

    Hello

Without specifying ``echo true`` you won't see the ``Hello`` string printed out when executing the above example.


.. _process-page-error-strategy:

errorStrategy
-------------

The ``errorStrategy`` directive allows you to define how an error condition is managed by the process. By default when
an error status is returned by the executed script, the process stops immediately. This in turn forces the entire pipeline
to terminate.

Table of available error strategies:

============== ==================
Name            Executor
============== ==================
``terminate``   Terminates the execution as soon as an error condition is reported. Pending jobs are killed (default)
``finish``      Initiates an orderly pipeline shutdown when an error condition is raised, waiting the completion of any submitted job.
``ignore``      Ignores processes execution errors.
``retry``       Re-submit for execution a process returning an error condition.
============== ==================


When setting the ``errorStrategy`` directive to ``ignore`` the process doesn't stop on an error condition,
it just reports a message notifying you of the error event.

For example::

    process ignoreAnyError {
       errorStrategy 'ignore'

       script:
       <your command string here>
    }

.. tip:: By definition a command script fails when it ends with a non-zero exit status.

The ``retry`` `error strategy`, allows you to re-submit for execution a process
returning an error condition. For example::

    process retryIfFail {
       errorStrategy 'retry'

       script:
       <your command string here>
    }


The number of times a failing process is re-executed is defined by the `maxRetries`_ and `maxErrors`_ directives.

.. note:: More complex strategies depending on the task exit status
  or other parametric values can be defined using a dynamic ``errorStrategy``
  directive. See the `Dynamic directives`_ section for details.

See also: `maxErrors`_, `maxRetries`_ and `Dynamic computing resources`_.

.. _process-executor:

executor
--------

The `executor` defines the underlying system where processes are executed. By default a process uses the executor
defined globally in the ``nextflow.config`` file.

The ``executor`` directive allows you to configure what executor has to be used by the process, overriding the default
configuration. The following values can be used:

=====================  ==================
Name                   Executor
=====================  ==================
``local``              The process is executed in the computer where `Nextflow` is launched.
``sge``                The process is executed using the Sun Grid Engine / `Open Grid Engine <http://gridscheduler.sourceforge.net/>`_.
``uge``                The process is executed using the `Univa Grid Engine <https://en.wikipedia.org/wiki/Univa_Grid_Engine/>`_ job scheduler.
``lsf``                The process is executed using the `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ job scheduler.
``slurm``              The process is executed using the SLURM job scheduler.
``pbs``                The process is executed using the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ job scheduler.
``pbspro``             The process is executed using the `PBS Pro <https://www.pbsworks.com/>`_ job scheduler.
``moab``               The process is executed using the `Moab <http://www.adaptivecomputing.com/moab-hpc-basic-edition/>`_ job scheduler.
``condor``             The process is executed using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ job scheduler.
``nqsii``              The process is executed using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ job scheduler.
``ignite``             The process is executed using the `Apache Ignite <https://ignite.apache.org/>`_ cluster.
``k8s``                The process is executed using the `Kubernetes <https://kubernetes.io/>`_ cluster.
``awsbatch``           The process is executed using the `AWS Batch <https://aws.amazon.com/batch/>`_ service.
``google-pipelines``   The process is executed using the `Google Genomics Pipelines <https://cloud.google.com/genomics/>`_ service.
=====================  ==================

The following example shows how to set the process's executor::


   process doSomething {

      executor 'sge'

      script:
      <your script here>

   }


.. note:: Each executor provides its own set of configuration options that can set be in the `directive` declarations block.
   See :ref:`executor-page` section to read about specific executor directives.

.. _process-ext:

ext
---

The ``ext`` is a special directive used as *namespace* for user custom process directives. This can be useful for
advanced configuration options. For example::

    process mapping {
      container "biocontainers/star:${task.ext.version}"

      input:
      file genome from genome_file
      set sampleId, file(reads) from reads_ch

      """
      STAR --genomeDir $genome --readFilesIn $reads
      """
    }

In the above example, the process uses a container whose version is controlled by the ``ext.version`` property.
This can be defined in the ``nextflow.config`` file as shown below::

    process.ext.version = '2.5.3'


.. _process-machineType:

machineType
-----------

The ``machineType`` can be used to specify a predefined Google Compute Platform `machine type <https://cloud.google.com/compute/docs/machine-types>`_
when running using the :ref:`Google Pipeline <google-pipelines>` executor.

This directive is optional and if specified overrides the cpus and memory directives::

    process foo {
      machineType 'n1-highmem-8'

      """
      <your script here>
      """
    }

.. note:: This feature requires Nextflow 19.07.0 or later.
    
See also: `cpus`_ and `memory`_.

.. _process-maxErrors:

maxErrors
---------

The ``maxErrors`` directive allows you to specify the maximum number of times a process can fail when using the ``retry`` `error strategy`.
By default this directive is disabled, you can set it as shown in the example below::

    process retryIfFail {
      errorStrategy 'retry'
      maxErrors 5

      """
      echo 'do this as that .. '
      """
    }
    
.. note:: This setting considers the **total** errors accumulated for a given process, across all instances. If you want
  to control the number of times a process **instance** (aka task) can fail, use ``maxRetries``.

See also: `errorStrategy`_ and `maxRetries`_.

.. _process-maxForks:

maxForks
--------

The ``maxForks`` directive allows you to define the maximum number of process instances that can be executed in parallel.
By default this value is equals to the number of CPU cores available minus 1.

If you want to execute a process in a sequential manner, set this directive to one. For example::

    process doNotParallelizeIt {

       maxForks 1

       '''
       <your script here>
       '''

    }

.. _process-maxRetries:

maxRetries
----------

The ``maxRetries`` directive allows you to define the maximum number of times a process instance can be
re-submitted in case of failure. This value is applied only when using the ``retry`` `error strategy`. By default
only one retry is allowed, you can increase this value as shown below::

    process retryIfFail {
        errorStrategy 'retry'
        maxRetries 3

        """
        echo 'do this as that .. '
        """
    }


.. note:: There is a subtle but important difference between ``maxRetries`` and the ``maxErrors`` directive.
    The latter defines the total number of errors that are allowed during the process execution (the same process can
    launch different execution instances), while the ``maxRetries`` defines the maximum number of times the same process
    execution can be retried in case of an error.

See also: `errorStrategy`_ and `maxErrors`_.


.. _process-memory:

memory
------

The ``memory`` directive allows you to define how much memory the process is allowed to use. For example::

    process big_job {

        memory '2 GB'
        executor 'sge'

        """
        your task script here
        """
    }


The following memory unit suffix can be used when specifying the memory value:

======= =============
Unit    Description
======= =============
B       Bytes
KB      Kilobytes
MB      Megabytes
GB      Gigabytes
TB      Terabytes
======= =============

.. This setting is equivalent to set the ``qsub -l virtual_free=<mem>`` command line option.

See also: `cpus`_, `time`_, `queue`_ and `Dynamic computing resources`_.


.. _process-module:

module
------

`Environment Modules <http://modules.sourceforge.net/>`_ is a package manager that allows you to dynamically configure
your execution environment and easily switch between multiple versions of the same software tool.

If it is available in your system you can use it with Nextflow in order to configure the processes execution
environment in your pipeline.

In a process definition you can use the ``module`` directive to load a specific module version to be used in the
process execution environment. For example::

  process basicExample {

    module 'ncbi-blast/2.2.27'

    """
    blastp -query <etc..>
    """
  }

You can repeat the ``module`` directive for each module you need to load. Alternatively multiple modules
can be specified in a single ``module`` directive by separating all the module names by using a ``:``
(colon) character as shown below::

   process manyModules {

     module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

     """
     blastp -query <etc..>
     """
  }


.. _process-penv:

penv
----

The ``penv`` directive  allows you to define the `parallel environment` to be used when submitting a parallel task to the
:ref:`SGE <sge-executor>` resource manager. For example::

    process big_job {

      cpus 4
      penv 'smp'
      executor 'sge'

      """
      blastp -query input_sequence -num_threads ${task.cpus}
      """
    }

This configuration depends on the parallel environment provided by your grid engine installation. Refer to your
cluster documentation or contact your admin to learn more about this.

.. note:: This setting is available when using the :ref:`sge-executor` executor.

See also: `cpus`_, `memory`_, `time`_

.. _process-pod:

pod
---

The ``pod`` directive allows the definition of pods specific settings, such as environment variables, secrets
and config maps when using the :ref:`k8s-executor` executor.

For example::

  process your_task {
    pod env: 'FOO', value: 'bar'

    '''
    echo $FOO
    '''
  }

The above snippet defines an environment variable named ``FOO`` which value is ``bar``.

The ``pod`` directive allows the definition of the following options:

================================================= =================================================
``label: <K>, value: <V>``                        Defines a pod label with key ``K`` and value ``V``.
``annotation: <K>, value: <V>``                   Defines a pod annotation with key ``K`` and value ``V``.
``env: <E>, value: <V>``                          Defines an environment variable with name ``E`` and whose value is given by the ``V`` string.
``env: <E>, fieldPath: <V>``                      Defines an environment variable with name ``E`` and whose value is given by the ``V`` `field path <https://kubernetes.io/docs/tasks/inject-data-application/environment-variable-expose-pod-information/>`_.
``env: <E>, config: <C/K>``                       Defines an environment variable with name ``E`` and whose value is given by the entry associated to the key with name ``K`` in the `ConfigMap <https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/>`_ with name ``C``.
``env: <E>, secret: <S/K>``                       Defines an environment variable with name ``E`` and whose value is given by the entry associated to the key with name ``K`` in the `Secret <https://kubernetes.io/docs/concepts/configuration/secret/>`_ with name ``S``.
``config: <C/K>, mountPath: </absolute/path>``    The content of the `ConfigMap <https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/>`_ with name ``C`` with key ``K`` is made available to the path ``/absolute/path``. When the key component is omitted the path is interpreted as a directory and all the `ConfigMap` entries are exposed in that path.
``secret: <S/K>, mountPath: </absolute/path>``    The content of the `Secret <https://kubernetes.io/docs/concepts/configuration/secret/>`_ with name ``S`` with key ``K`` is made available to the path ``/absolute/path``. When the key component is omitted the path is interpreted as a directory and all the `Secret` entries are exposed in that path.
``volumeClaim: <V>, mountPath: </absolute/path>`` Mounts a `Persistent volume claim <https://kubernetes.io/docs/concepts/storage/persistent-volumes/>`_ with name ``V`` to the specified path location. Use the optional `subPath` parameter to mount a directory inside the referenced volume instead of its root. The volume may be mounted with `readOnly: true`, but is read/write by default.
``imagePullPolicy: <V>``                          Specifies the strategy to be used to pull the container image e.g. ``imagePullPolicy: 'Always'``.
``imagePullSecret: <V>``                          Specifies the secret name to access a private container image registry. See `Kubernetes documentation <https://kubernetes.io/docs/concepts/containers/images/#specifying-imagepullsecrets-on-a-pod>`_ for details.
``runAsUser: <UID>``                              Specifies the user ID to be used to run the container. Shortcut for the ``securityContext`` option.
``securityContext: <V>``                          Specifies the pod security context. See `Kubernetes security context <https://kubernetes.io/docs/tasks/configure-pod-container/security-context/>`_ for details.
``nodeSelector: <V>``                             Specifies which node the process will run on. See `Kubernetes nodeSelector <https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#nodeselector>`_ for details.
``affinity: <V>``                                 Specifies affinity for which nodes the process should run on. See `Kubernetes affinity <https://kubernetes.io/docs/concepts/scheduling-eviction/assign-pod-node/#affinity-and-anti-affinity>`_ for details.
``automountServiceAccountToken: <V>``             Specifies whether to `automount service account token <https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/>`_ into process pods. If ``V`` is true, service account token is automounted into task pods (default).
``priorityClassName: <V>``                        Specifies the `priority class name <https://kubernetes.io/docs/concepts/scheduling-eviction/pod-priority-preemption/>`_ for pods.
================================================= =================================================

When defined in the Nextflow configuration file, a pod setting can be defined using the canonical
associative array syntax. For example::

  process {
    pod = [env: 'FOO', value: 'bar']
  }

When more than one setting needs to be provides they must be enclosed in a list definition as shown below::

  process {
    pod = [ [env: 'FOO', value: 'bar'], [secret: 'my-secret/key1', mountPath: '/etc/file.txt'] ]
  }


.. _process-publishDir:

publishDir
----------

The ``publishDir`` directive allows you to publish the process output files to a specified folder. For example::


    process foo {

        publishDir '/data/chunks'

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }


The above example splits the string ``Hola`` into file chunks of a single byte. When complete the ``chunk_*`` output files
are published into the ``/data/chunks`` folder.

.. tip:: The ``publishDir`` directive can be specified more than one time in to publish the output files
  to different target directories. This feature requires version 0.29.0 or higher.

By default files are published to the target folder creating a *symbolic link* for each process output that links
the file produced into the process working directory. This behavior can be modified using the ``mode`` parameter.

Table of optional parameters that can be used with the ``publishDir`` directive:

=============== =================
Name            Description
=============== =================
mode            The file publishing method. See the following table for possible values.
overwrite       When ``true`` any existing file in the specified folder will be overridden (default: ``true`` during normal
                pipeline execution and ``false`` when pipeline execution is `resumed`).
pattern         Specifies a `glob`_ file pattern that selects which files to publish from the overall set of output files.
path            Specifies the directory where files need to be published. **Note**: the syntax ``publishDir '/some/dir'`` is a shortcut for ``publishDir path: '/some/dir'``.
saveAs          A closure which, given the name of the file being published, returns the actual file name or a full path where the file is required to be stored.
                This can be used to rename or change the destination directory of the published files dynamically by using
                a custom strategy.
                Return the value ``null`` from the closure to *not* publish a file.
                This is useful when the process has multiple output files, but you want to publish only some of them.
enabled         Enable or disable the publish rule depending on the boolean value specified (default: ``true``).
tags            Allow to associate tags with the target file e.g. ``tag: [FOO: 'Hello world']`` (EXPERIMENTAL, currently only supported by files stored on AWS S3, requires version ``21.12.0-edge`` or later).
=============== =================

Table of publish modes:

=============== =================
 Mode           Description
=============== =================
symlink         Creates an absolute `symbolic link` in the published directory for each process output file (default).
rellink         Creates a relative `symbolic link` in the published directory for each process output file.
link            Creates a `hard link` in the published directory for each process output file.
copy            Copies the output files into the published directory.
copyNoFollow    Copies the output files into the published directory without following symlinks ie. copies the links themselves. 
move            Moves the output files into the published directory. **Note**: this is only supposed to be used for a `terminating` process i.e. a process whose output is not consumed by any other downstream process.
=============== =================

.. note:: The `mode` value needs to be specified as a string literal i.e. enclosed by quote characters. Multiple parameters
  need to be separated by a colon character. For example:

::

    process foo {

        publishDir '/data/chunks', mode: 'copy', overwrite: false

        output:
        file 'chunk_*' into letters

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }


.. warning:: Files are copied into the specified directory in an *asynchronous* manner, thus they may not be immediately
  available in the published directory at the end of the process execution. For this reason files published by a process
  must not be accessed by other downstream processes.


.. _process-queue:

queue
-----

The ``queue`` directory allows you to set the `queue` where jobs are scheduled when using a grid based executor
in your pipeline. For example::

    process grid_job {

        queue 'long'
        executor 'sge'

        """
        your task script here
        """
    }


Multiple queues can be specified by separating their names with a comma for example::

    process grid_job {

        queue 'short,long,cn-el6'
        executor 'sge'

        """
        your task script here
        """
    }


.. note:: This directive is taken in account only by the following executors: :ref:`sge-executor`, :ref:`lsf-executor`,
  :ref:`slurm-executor` and :ref:`pbs-executor` executors.


.. _process-label:

label
-----

The ``label`` directive allows the annotation of processes with mnemonic identifier of your choice.
For example::

  process bigTask {

    label 'big_mem'

    '''
    <task script>
    '''
  }


The same label can be applied to more than a process and multiple labels can be applied to the same
process using the ``label`` directive more than one time.

.. note:: A label must consist of alphanumeric characters or ``_``, must start with an alphabetic character
  and must end with an alphanumeric character.

Labels are useful to organise workflow processes in separate groups which can be referenced
in the configuration file to select and configure subset of processes having similar computing requirements.

See the :ref:`config-process-selectors` documentation for details.


.. _process-scratch:

scratch
-------

The ``scratch`` directive allows you to execute the process in a temporary folder that is local to the execution node.

This is useful when your pipeline is launched by using a `grid` executor, because it allows you to decrease the NFS
overhead by running the pipeline processes in a temporary directory in the local disk of the actual execution node.
Only the files declared as output in the process definition will be copied in the pipeline working area.

In its basic form simply specify ``true`` at the directive value, as shown below::

  process simpleTask {

    scratch true

    output:
    file 'data_out'

    '''
    <task script>
    '''
  }


By doing this, it tries to execute the script in the directory defined by the variable ``$TMPDIR`` in the execution node.
If this variable does not exist, it will create a new temporary directory by using the Linux command ``mktemp``.

A custom environment variable, other than ``$TMPDIR``, can be specified by simply using it as the scratch value, for
example::

  scratch '$MY_GRID_TMP'

Note, it must be wrapped by single quotation characters, otherwise the variable will be evaluated in the
pipeline script context.

You can also provide a specific folder path as scratch value, for example::

  scratch '/tmp/my/path'

By doing this, a new temporary directory will be created in the specified path each time a process is executed.

Finally, when the ``ram-disk`` string is provided as ``scratch`` value, the process will be execute in the node
RAM virtual disk.

Summary of allowed values:

=========== ==================
scratch     Description
=========== ==================
false       Do not use the scratch folder.
true        Creates a scratch folder in the directory defined by the ``$TMPDIR`` variable; fallback to ``mktemp /tmp`` if that variable do not exists.
$YOUR_VAR   Creates a scratch folder in the directory defined by the ``$YOUR_VAR`` environment variable; fallback to ``mktemp /tmp`` if that variable do not exists.
/my/tmp     Creates a scratch folder in the specified directory.
ram-disk    Creates a scratch folder in the RAM disk ``/dev/shm/`` (experimental).
=========== ==================

.. _process-storeDir:

storeDir
--------

The ``storeDir`` directive allows you to define a directory that is used as `permanent` cache for your process results.

In more detail, it affects the process execution in two main ways:

#. The process is executed only if the files declared in the `output` clause do not exist in the directory specified by
   the ``storeDir`` directive. When the files exist the process execution is skipped and these files are used as
   the actual process result.

#. Whenever a process is successfully completed the files listed in the `output` declaration block are moved into the directory
   specified by the ``storeDir`` directive.

The following example shows how to use the ``storeDir`` directive to create a directory containing a BLAST database
for each species specified by an input parameter::

  genomes = Channel.fromPath(params.genomes)

  process formatBlastDatabases {

    storeDir '/db/genomes'

    input:
    file species from genomes

    output:
    file "${dbName}.*" into blastDb

    script:
    dbName = species.baseName
    """
    makeblastdb -dbtype nucl -in ${species} -out ${dbName}
    """

  }


.. warning:: The ``storeDir`` directive is meant for long term process caching and should not be used to
    output the files produced by a process to a specific folder or organise result data in `semantic` directory structure.
    In these cases you may use the `publishDir`_ directive instead.

.. note:: The use of AWS S3 path is supported however it requires the installation of the `AWS CLI tool <https://aws.amazon.com/cli/>`_
  (i.e. ``aws``) in the target computing node.

.. _process-stageInMode:

stageInMode
-----------

The ``stageInMode`` directive defines how input files are staged-in to the process work directory. The following values
are allowed:

======= ==================
Value   Description
======= ==================
copy    Input files are staged in the process work directory by creating a copy.
link    Input files are staged in the process work directory by creating an (hard) link for each of them.
symlink Input files are staged in the process work directory by creating a symbolic link with an absolute path for each of them (default).
rellink Input files are staged in the process work directory by creating a symbolic link with a relative path for each of them.
======= ==================


.. _process-stageOutMode:

stageOutMode
------------

The ``stageOutMode`` directive defines how output files are staged-out from the scratch directory to the process work
directory. The following values are allowed:

======= ==================
Value   Description
======= ==================
copy    Output files are copied from the scratch directory to the work directory.
move    Output files are moved from the scratch directory to the work directory.
rsync   Output files are copied from the scratch directory to the work directory by using the ``rsync`` utility.
======= ==================

See also: `scratch`_.


.. _process-tag:

tag
---

The ``tag`` directive allows you to associate each process execution with a custom label, so that it will be easier
to identify them in the log file or in the trace execution report. For example::

    process foo {
      tag "$code"

      input:
      val code from 'alpha', 'gamma', 'omega'

      """
      echo $code
      """
    }

The above snippet will print a log similar to the following one, where process names contain the tag value::

    [6e/28919b] Submitted process > foo (alpha)
    [d2/1c6175] Submitted process > foo (gamma)
    [1c/3ef220] Submitted process > foo (omega)


See also :ref:`Trace execution report <trace-report>`


.. _process-time:

time
----

The ``time`` directive allows you to define how long a process is allowed to run. For example::

    process big_job {

        time '1h'

        """
        your task script here
        """
    }



The following time unit suffixes can be used when specifying the duration value:

+---------------------------------+--------------+
| Unit                            | Description  |
+=================================+==============+
| `ms`, `milli`, `millis`         | Milliseconds |
+---------------------------------+--------------+
| `s`, `sec`, `second`, `seconds` | Seconds      |
+---------------------------------+--------------+
| `m`, `min`, `minute`, `minutes` | Minutes      |
+---------------------------------+--------------+
| `h`, `hour`, `hours`            | Hours        |
+---------------------------------+--------------+
| `d`, `day`, `days`              | Days         |
+---------------------------------+--------------+

Multiple units can be used in a single declaration, for example: ``'1day 6hours 3minutes 30seconds'``

.. note:: This directive is taken in account only when using one of the following grid based executors:
  :ref:`sge-executor`, :ref:`lsf-executor`, :ref:`slurm-executor`, :ref:`pbs-executor`,
  :ref:`condor-executor` and :ref:`awsbatch-executor` executors.

See also: `cpus`_, `memory`_, `queue`_ and `Dynamic computing resources`_.


Dynamic directives
------------------

A directive can be assigned *dynamically*, during the process execution, so that its actual value can be evaluated
depending on the value of one, or more, process' input values.

In order to be defined in a dynamic manner the directive's value needs to be expressed by using a :ref:`closure <script-closure>`
statement, as in the following example::

    process foo {

      executor 'sge'
      queue { entries > 100 ? 'long' : 'short' }

      input:
      set entries, file(x) from data

      script:
      """
      < your job here >
      """
    }

In the above example the `queue`_ directive is evaluated dynamically, depending on the input value ``entries``. When it is
bigger than 100, jobs will be submitted to the queue ``long``, otherwise the ``short`` one will be used.

All directives can be assigned to a dynamic value except the following:

* `executor`_
* `maxForks`_


.. tip::
  Directives taking a string value containing one or more variables are always resolved in a dynamic manner, and therefore
  it's semantically equivalent to the above above syntax. Therefore the above directive can also be written as::

    queue "${ entries > 100 ? 'long' : 'short' }"

  Note however the latter syntax can be used both for directive main argument (like in the ``queue`` example) and for directive
  optional named attributes. Instead the closure based syntax is only resolved dynamically for the directive main argument.

.. note:: You can retrieve the current value of a dynamic directive in the process script by using the implicit variable ``task``
  which holds the directive values defined in the current process instance.

For example::


   process foo {

      queue { entries > 100 ? 'long' : 'short' }

      input:
      set entries, file(x) from data

      script:
      """
      echo Current queue: ${task.queue}
      """
    }


Dynamic computing resources
---------------------------

It's a very common scenario that different instances of the same process may have very different needs in terms of computing resources. 
In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. 
Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your jobs.

The `Dynamic directives`_ evaluation feature can be used to modify the amount of computing resources requested in case
of a process failure and try to re-execute it using a higher limit. For example::


    process foo {

        memory { 2.GB * task.attempt }
        time { 1.hour * task.attempt }

        errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
        maxRetries 3

        script:
        <your job here>

    }


In the above example the `memory`_ and execution `time`_ limits are defined dynamically. The first time the process
is executed the ``task.attempt`` is set to ``1``, thus it will request a two GB of memory and one hour of maximum execution
time.

If the task execution fail reporting an exit status in the range between 137 and 140, the task is re-submitted (otherwise terminates immediately).
This time the value of ``task.attempt`` is ``2``, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.

The directive `maxRetries`_ set the maximum number of time the same task can be re-executed.

Dynamic Retry with backoff
--------------------------

There are cases in which the required execution resources may be temporary unavailable e.g.
network congestion. In these cases immediately re-executing the task will likely result in
the identical error. A retry with an exponential backoff delay can better recover these error
conditions::

    process foo {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
      maxRetries 5
      script:
      '''
      your_command --here
      '''
    }
.. _sharing-page:

****************
Pipeline sharing
****************

Nextflow seamlessly integrates with `BitBucket <http://bitbucket.org/>`_ [#]_, `GitHub <http://github.com>`_,
and `GitLab <http://gitlab.com>`_ hosted code repositories and sharing platforms. This feature allows you to manage your
project code in a more consistent manner or use other people's Nextflow pipelines, published through BitBucket/GitHub/GitLab,
in a quick and transparent way.

How it works
============

When you launch a script execution with Nextflow, it will look for a file with the pipeline name you've specified.
If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified).
If it is found, the repository is automatically downloaded to your computer and executed. This repository is
stored in the Nextflow home directory, that is by default the ``$HOME/.nextflow`` path, and thus will be reused for any further
executions.

Running a pipeline
==================

To launch the execution of a pipeline project, hosted in a remote code repository, you simply need to specify its `qualified` name
or the repository URL after the ``run`` command. The qualified name is formed by two parts: the `owner` name and the
`repository` name separated by a ``/`` character.

In other words if a Nextflow project is hosted, for example, in a GitHub repository at the address
``http://github.com/foo/bar``, it can be executed by entering the following command in your shell terminal::

    nextflow run foo/bar

or using the project URL::

    nextflow run http://github.com/foo/bar


.. note:: In the first case, if your project is hosted on a service other than GitHub, you will need to specify this hosting
    service in the command line by using the ``-hub`` option. For example ``-hub bitbucket`` or ``-hub gitlab``.
    In the second case, i.e. when using the project URL as name, the ``-hub`` option is not needed.

You can try this feature out by simply entering the following command in your shell terminal::

    nextflow run nextflow-io/hello

It will download a trivial `Hello` example from the repository published at the following address
http://github.com/nextflow-io/hello and execute it in your computer.

If the `owner` part in the pipeline name is omitted, Nextflow will look for a pipeline between the ones you have
already executed having a name that matches the name specified. If none is found it will try to download
it using the `organisation` name defined by the environment variable ``NXF_ORG`` (which by default is ``nextflow-io``).


.. tip:: To access a private repository, specify the access credentials by using the ``-user`` command
    line option, then the program will ask you to enter the password interactively.
    Private repositories access credentials can also be defined in the `SCM configuration file`_.


Handling revisions
==================

Any Git branch, tag or commit ID defined in a project repository, can be used to specify the revision that you want to execute
when launching a pipeline by adding the ``-r`` option to the run command line. So for example you could enter::

    nextflow run nextflow-io/hello -r mybranch

or ::

    nextflow run nextflow-io/hello -r v1.1


It will execute two different project revisions corresponding to the Git tag/branch having that names.

Commands to manage projects
===========================

The following commands allows you to perform some basic operations that can be used to manage your projects.

.. note:: Nextflow is not meant to replace functionalities provided by the `Git <https://git-scm.com/>`_ tool. You may still need it to create new
  repositories or commit changes, etc.

Listing available projects
--------------------------

The ``list`` command allows you to list all the projects you have downloaded in your computer. For example::

    nextflow list

This prints a list similar to the following one::

    cbcrg/ampa-nf
    cbcrg/piper-nf
    nextflow-io/hello
    nextflow-io/examples


Showing project information
---------------------------

By using the ``info`` command you can show information from a downloaded project. For example::

     project name: nextflow-io/hello
     repository  : http://github.com/nextflow-io/hello
     local path  : $HOME/.nextflow/assets/nextflow-io/hello
     main script : main.nf
     revisions   :
     * master (default)
       mybranch
       v1.1 [t]
       v1.2 [t]

Starting from the top it shows: 1) the project name; 2) the Git repository URL; 3) the local folder where the
project has been downloaded; 4) the script that is executed when launched; 5) the list of available
revisions i.e. branches and tags. Tags are marked with a ``[t]`` on the right, the current checked-out revision is
marked with a ``*`` on the left.

Pulling or updating a project
-----------------------------

The ``pull`` command allows you to download a project from a GitHub repository or to update it if
that repository has already been downloaded. For example::

    nextflow pull nextflow-io/examples

Altenatively, you can use the repository URL as the name of the project to pull::

    nextflow pull https://github.com/nextflow-io/examples


Downloaded pipeline projects are stored in the folder ``$HOME/.nextflow/assets`` in your computer.


Viewing the project code
-------------------------

The ``view`` command allows you to quickly show the content of the pipeline script you have downloaded. For example::

    nextflow view nextflow-io/hello

By adding the ``-l`` option to the example above it will list the content of the repository.


Cloning a project into a folder
-------------------------------

The ``clone`` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example::

    nextflow clone nextflow-io/hello target-dir

If the destination directory is omitted the specified project is cloned to a directory with the same name as the
pipeline base name (e.g. `hello`) in the current folder.

The clone command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push
back your changes by using the usual Git/GitHub workflow.

Deleting a downloaded project
-----------------------------

Downloaded pipelines can be deleted by using the ``drop`` command, as shown below::

    nextflow drop nextflow-io/hello


.. _sharing-scm-file:

SCM configuration file
=======================

The file ``$HOME/.nextflow/scm`` allows you to centralise the security credentials required to access private project
repositories on Bitbucket, GitHub and GitLab source code management (`SCM`) platforms or to manage the configuration properties
of private server installations (of the same platforms).

The configuration properties for each SCM platform are defined inside the ``providers`` section,
properties for the same provider are grouped together with a common name and delimited with curly brackets as in this example::

    providers {
        <provider-name> {
            property = value
            :
        }
    }


In the above template replace `<provider-name>` with one of the "default" servers (i.e. ``bitbucket``, ``github`` or ``gitlab``)
or a custom identifier representing a private SCM server installation.

The following configuration properties are supported for each provider configuration:

=================== ==============
Name                Description
=================== ==============
user                User name required to access private repositories on the SCM server.
password            User password required to access private repositories on the SCM server.
token               Private API access token (used only when the specified platform is ``gitlab``).
:sup:`*` platform   SCM platform name, either: ``github``, ``gitlab`` or ``bitbucket``.
:sup:`*` server     SCM server name including the protocol prefix e.g. ``https://github.com``.
:sup:`*` endpoint   SCM API `endpoint` URL e.g. ``https://api.github.com`` (default: the same value specified for ``server``).
=================== ==============

The attributes marked with a * are only required when defining the configuration of a private SCM server.

.. tip::
  A custom location for the SCM file can be specified using the ``NXF_SCM_FILE`` environment variable (requires
 version ``20.10.0`` or later).

BitBucket credentials
---------------------

Create a ``bitbucket`` entry in the `SCM configuration file`_ specifying your user name and app password, as shown below::

    providers {

        bitbucket {
            user = 'me'
            password = 'my-secret'
        }

    }


.. note::
   App passwords are substitute passwords for a user account which you can use for scripts and integrating
   tools to avoid putting your real password into configuration files.
   Learn more at `this link <https://support.atlassian.com/bitbucket-cloud/docs/app-passwords/>`_.

BitBucket Server credentials
-----------------------------

`BitBucket Server <https://confluence.atlassian.com/bitbucketserver>`_ is a self-hosted Git repository and management
platform.

.. note::
    BitBucket Server uses different API from the `BitBucket <https://bitbucket.org/>`_ cloud service. Make sure to
    use the right configuration whether you are using the cloud service or a self-hosted installation.

To access your local BitBucket Server create an entry in the `SCM configuration file`_ specifying as shown below::

        providers {

            mybitbucket {
                platform = 'bitbucketserver'
                server = 'https://your.bitbucket.host.com'
                endpoint = 'https://your.bitbucket.host.com'
                user = 'your-user'
                password = 'your-password or your-token'
            }

        }


GitHub credentials
------------------

Create a ``github`` entry in the `SCM configuration file`_ specifying your user name and access token as shown below::

    providers {

        github {
            user = 'your-user-name'
            password = 'your-personal-access-token;'
        }

    }

.. tip:: GitHub requires the use of the personal access token (PAT) in place of password field when accessing APIs.
  Learn more about PAT and how to create it at `this link <https://docs.github.com/en/github/authenticating-to-github/keeping-your-account-and-data-secure/creating-a-personal-access-token>`_.


GitLab credentials
-------------------

Create a ``gitlab`` entry in the `SCM configuration file`_ specifying the user name, password and your API access token
that can be found in your GitLab `account page <https://gitlab.com/profile/account>`_ (sign in required). For example::

    providers {

        gitlab {
            user = 'me'
            password = 'my-secret'
            token = 'YgpR8m7viH_ZYnC8YSe8'
        }

    }

.. tip::
    The GitLab *token* string can be used as the ``password`` value in the above setting.
    When doing that the ``token`` field can be omitted.


Gitea credentials
-----------------

`Gitea <https://gitea.io>`_ is a Git repository server with GitHub-like GUI access. Since Gitea installation is quite
easy, it is suitable for building a private development environment in your network. To access your Gitea server, you
have to provide all the credential information below::

    providers {

        mygitea {
            server = 'http://your-domain.org/gitea'
            endpoint = 'http://your-domain.org/gitea/api/v1'
            platform = 'gitea'
            user = 'your-user'
            password = 'your-password'
            token = 'your-api-token'
        }

    }


See `Gitea documentation <https://docs.gitea.io/en-us/api-usage/>`_ about how to enable API access on your
server and how to issue a token.

Azure Repos credentials
-----------------------

Nextflow has a builtin support for `Azure Repos <https://azure.microsoft.com/en-us/services/devops/repos/>`_, a Git source
code management service hosted in the Azure cloud. To access your Azure Repos with Nextflow provide the repository credentials
using the configuration snippet shown below:

    providers {

        azurerepos {
            user = 'your-user-name'
            password = 'your-personal-access-token'
        }

    }

.. tip::
  The Personal access token can be generated in the repository `Clone Repository` dialog.

Private server configuration
============================

Nextflow is able to access repositories hosted on private BitBucket, GitHub, GitLab and Gitea server installations.

In order to use a private SCM installation you will need to set the server name and access credentials
in your `SCM configuration file`_ .

If, for example, the host name of your private GitLab server is ``gitlab.acme.org``, you will need to have in the
``$HOME/.nextflow/scm`` file a configuration like the following::

    providers {

        mygit {
            server = 'http://gitlab.acme.org'
            platform = 'gitlab'
            user = 'your-user'
            password = 'your-password'
            token = 'your-api-token'
        }

    }


Then you will be able to run/pull a project with Nextflow using the following command line::

    $ nextflow run foo/bar -hub mygit

Or, in alternative, using the Git clone URL::

    $ nextflow run http://gitlab.acme.org/foo/bar.git

.. note::
    You must also specify the server API endpoint URL if it differs from the server
    base URL. For example, for GitHub Enterprise V3, add
    ``endpoint = 'https://git.your-domain.com/api/v3'``.

.. warning:: When accessing a private SCM installation over ``https`` and that server uses a custom SSL certificate
  you may need to import such certificate into your local Java keystore. Read more
  `here <https://docs.oracle.com/javase/tutorial/security/toolsign/rstep2.html>`_.

Local repository configuration
==============================

Nextflow is also able to handle repositories stored in a local or shared file system. The repository
must be created as a `bare repository <https://mijingo.com/blog/what-is-a-bare-git-repository>`_.


Having, for example. a bare repository store at path ``/shared/projects/foo.git``, Nextflow is able
to run it using the following syntax::

  $ nextflow run file:/shared/projects/foo.git

See `Git documentation <https://git-scm.com/book/en/v2/Git-on-the-Server-Getting-Git-on-a-Server>`_ for
more details about how create and manage bare repositories.

Publishing your pipeline
========================

In order to publish your Nextflow pipeline to GitHub (or any other supported platform) and allow other people to use it,
you only need to create a GitHub repository containing all your project script and data files. If you don't know how to
do it, follow this simple tutorial that explains how `create a GitHub repository <https://help.github.com/articles/create-a-repo>`_.

Nextflow only requires that the main script in your pipeline project is called ``main.nf``. A different name can be
used by specifying the ``manifest.mainScript`` attribute in the ``nextflow.config`` file that must be
included in your project. For example::

  manifest.mainScript = 'my_very_long_script_name.nf'

To learn more about this and other project metadata information, that can be defined in the Nextflow configuration file,
read the :ref:`Manifest <config-manifest>` section on the Nextflow configuration page.

Once you have uploaded your pipeline project to GitHub other people can execute it simply
using the project name or the repository URL.

For if your GitHub account name is ``foo`` and you have uploaded a project into a repository named ``bar`` the
repository URL will be ``http://github.com/foo/bar`` and people will able to download and run it by using either
the command::

    nextflow run foo/bar

or

::

    nextflow run http://github.com/foo/bar

See the `Running a pipeline`_ section for more details on how to run Nextflow projects.

Manage dependencies
=====================

Computational pipelines are rarely composed by a single script. In real world applications they depend on dozens of other components.
These can be other scripts, databases, or applications compiled for a platform native binary format.

External dependencies are the most common source of problems when sharing a piece of software, because the
users need to have an identical set of tools and the same configuration to be able to use it. In many cases this has proven to be
a painful and error prone process, that can severely limit the ability to reproduce computational results on a system other than
the one on which it was originally developed.

Nextflow tackles this problem by integrating GitHub, BitBucket and GitLab sharing platforms and
`Docker <http://www.docker.com>`_ containers technology.

The use of a code management system is important to keep together all the dependencies of your
pipeline project and allows you to track the changes of the source code in a consistent manner.

Moreover to guarantee that a pipeline is reproducible it should be self-contained i.e. it should have ideally no
dependencies on the hosting environment. By using Nextflow you can achieve this goal following these methods:

Third party scripts
--------------------

Any third party script that does not need to be compiled (Bash, Python, Perl, etc) can be included in the pipeline
project repository, so that they are distributed with it.

Grant the execute permission to these files and copy them into a folder named ``bin/`` in the root directory of your
project repository. Nextflow will automatically add this folder to the ``PATH`` environment variable, and the scripts
will automatically be accessible in your pipeline without the need to specify an absolute path to invoke them.

System environment
--------------------

Any environment variable that may be required by the tools in your pipeline can be defined in the ``nextflow.config`` file
by using the ``env`` scope and including it in the root directory of your project. For example::

  env {
    DELTA = 'foo'
    GAMMA = 'bar'
  }


See the :ref:`config-page` page to learn more about the Nextflow configuration file.

Resource manager
--------------------

When using Nextflow you don't need to write the code to parallelize your pipeline for a specific grid engine/resource
manager because the parallelization is defined implicitly and managed by the Nextflow runtime. The target execution
environment is parametrized and defined in the configuration file, thus your code is free from this kind of dependency.

Bootstrap data
--------------------

Whenever your pipeline requires some files or dataset to carry out any initialization step, you
can include this data in the pipeline repository itself and distribute them together.

To reference this data in your pipeline script in a portable manner (i.e. without the need to use a static absolute path)
use the implicit variable ``baseDir`` which locates the base directory of your pipeline project.

For example, you can create a folder named ``dataset/`` in your repository root directory and copy there the
required data file(s) you may need, then you can access this data in your script by writing::

   sequences = file("$baseDir/dataset/sequences.fa")
   sequences.splitFasta {
        println it
    }

User inputs
-------------

Nextflow scripts can be easily parametrised to allow users to provide their own input data. Simply declare on the
top of your script all the parameters it may require as shown below::

  params.my_input = 'default input file'
  params.my_output = 'default output path'
  params.my_flag = false
  ..

The actual parameter values can be provided when launching the script execution on the command line
by prefixed the parameter name with a double minus character i.e. ``--``, for example::

  nextflow run <your pipeline> --my_input /path/to/input/file --my_output /other/path --my_flag true




Binary applications
--------------------

Docker allows you to ship any binary dependencies that you may have in your pipeline to a portable image
that is downloaded on-demand and can be executed on any platform where a Docker engine is installed.

In order to use it with Nextflow, create a Docker image containing the tools needed by your pipeline and make it available
in the `Docker registry <https://registry.hub.docker.com>`_.

Then declare in the ``nextflow.config`` file, that you will include in your project, the name of the Docker image you
have created. For example::

  process.container = 'my-docker-image'
  docker.enabled = true

In this way when you launch the pipeline execution, the Docker image will be automatically downloaded and used to run
your tasks.

Read the :ref:`docker-page` page to learn more on how to use Docker containers with Nextflow.


This mix of technologies makes it possible to write self-contained and truly reproducible pipelines which require
zero configuration and can be reproduced in any system having a Java VM and a Docker engine installed.


.. [#] BitBucket provides two types of version control system: `Git` and `Mercurial`. Nextflow supports only `Git` based repositories.
.. _google-page:

************
Google Cloud
************

Requirements
============

Nextflow
--------
The support for Google Cloud requires Nextflow version ``20.01.0`` or later. To install it define the following variables
in your system environment::

    export NXF_VER=20.01.0
    export NXF_MODE=google

.. note:: As of version ``21.04.0`` or later the above variables are not required anymore and therefore should not be used.

Credentials
-----------

Credentials for submitting requests to the Google LifeSciences API are picked up from your
environment using `Application Default Credentials <https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http>`_.
Application Default Credentials are designed to use the credentials most natural to the
environment in which a tool runs.

The most common case will be to pick up your end-user Google credentials from your
workstation. You can create these by running the command::

    gcloud auth application-default login 

and running through the authentication flow. This will write a credential file to your gcloud
configuration directory that will be used for any tool you run on your workstation that
picks up default credentials.

The next most common case would be when running on a Compute Engine VM. In this case,
Application Default Credentials will pick up the Compute Engine Service Account
credentials for that VM.

See the `Application Default Credentials <https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http>`_ documentation for how to enable other use cases.


Finally, the ``GOOGLE_APPLICATION_CREDENTIALS`` environment variable can be used to specify location
of the Google credentials file.

If you don't have it, the credentials file can be download from the Google Cloud Console following these steps:

* Open the `Google Cloud Console <https://console.cloud.google.com>`_
* Go to APIs & Services → Credentials
* Click on the *Create credentials* (blue) drop-down and choose *Service account key*, in the following page
* Select an existing *Service account* or create a new one if needed
* Select JSON as *Key type*
* Click the *Create* button and download the JSON file giving a name of your choice e.g. ``creds.json``.

Then, define the following variable replacing the path in the example with the one of your
credentials file just downloaded::

    export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json

.. _google-lifesciences:

Cloud Life Sciences
===================

`Cloud Life Sciences <https://cloud.google.com/life-sciences/>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Cloud Life Sciences API which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions through the Google Cloud service.

.. note::
  This features requires Nextflow ``20.01.0-edge`` or later.

.. warning::
  This API works well for coarse-grained workloads i.e. long running jobs. It's not suggested the use
  this feature for pipelines spawning many short lived tasks.

.. _google-lifesciences-config:

Configuration
-------------

Make sure to have defined in your environment the ``GOOGLE_APPLICATION_CREDENTIALS`` variable.
See the section `Requirements`_ for details.

.. tip:: Make sure to have enabled Cloud Life Sciences API to use this feature. To learn how to enable it
  follow `this link <https://cloud.google.com/life-sciences/docs/quickstart>`_.

Create a ``nextflow.config`` file in the project root directory. The config must specify the following parameters:

* Google Life Sciences as Nextflow executor i.e. ``process.executor = 'google-lifesciences'``.
* The Docker container images to be used to run pipeline tasks e.g. ``process.container = 'biocontainers/salmon:0.8.2--1'``.
* The Google Cloud `project` ID to run in e.g. ``google.project = 'rare-lattice-222412'``.
* The Google Cloud `region` or `zone`. This is where the Compute Engine VMs will be started.
  You need to specify either one, **not** both. Multiple regions or zones can be specified by
  separating them with a comma e.g. ``google.zone = 'us-central1-f,us-central-1-b'``.

Example::

    process {
        executor = 'google-lifesciences'
        container = 'your/container:latest'
    }

    google {
        project = 'your-project-id'
        zone = 'europe-west1-b'
    }


.. warning:: Make sure to specify in the above setting the project ID not the project name.

.. Note:: A container image must be specified to deploy the process execution. You can use a different Docker image for
  each process using one or more :ref:`config-process-selectors`.

The following configuration options are available:

============================================== =================
Name                                           Description
============================================== =================
google.project                                 The Google Project Id to use for the pipeline execution.
google.region                                  The Google *region* where the computation is executed in Compute Engine VMs. Multiple regions can be provided separating them by a comma. Do not specify if a zone is provided. See  `available Compute Engine regions and zones <https://cloud.google.com/compute/docs/regions-zones/>`_
google.zone                                    The Google *zone* where the computation is executed in Compute Engine VMs. Multiple zones can be provided separating them by a comma. Do not specify if a region is provided. See  `available Compute Engine regions and zones <https://cloud.google.com/compute/docs/regions-zones/>`_
google.location                                The Google *location* where the job executions are deployed to Cloud Life Sciences API. See  `available Cloud Life Sciences API locations <https://cloud.google.com/life-sciences/docs/concepts/locations>`_ (default: the same as the region or the zone specified).
google.enableRequesterPaysBuckets              When ``true`` uses the configured Google project id as the billing project for storage access. This is required when accessing data from *requester pays enabled* buckets. See `Requester Pays on Google Cloud Storage documentation  <https://cloud.google.com/storage/docs/requester-pays>`_ (default: ``false``)
google.lifeSciences.cpuPlatform                Set the minimum CPU Platform e.g. `'Intel Skylake'`. See `Specifying a minimum CPU Platform for VM instances <https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications>`_ (default: none).
google.lifeSciences.bootDiskSize               Set the size of the virtual machine boot disk e.g `50.GB` (default: none).
google.lifeSciences.copyImage                  The container image run to copy input and output files. It must include the ``gsutil`` tool (default: ``google/cloud-sdk:alpine``).
google.lifeSciences.debug                      When ``true`` copies the `/google` debug directory in that task bucket directory (default: ``false``)
google.lifeSciences.preemptible                When ``true`` enables the usage of *preemptible* virtual machines or ``false`` otherwise (default: ``true``)
google.lifeSciences.usePrivateAddress          When ``true`` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: ``false``). Requires version ``20.03.0-edge`` or later.
google.lifeSciences.network                    Set network name to attach the VM's network interface to. The value will be prefixed with global/networks/ unless it contains a /, in which case it is assumed to be a fully specified network resource URL. If unspecified, the global default network is used. Requires version ``21.03.0-edge`` or later.
google.lifeSciences.serviceAccountEmail        Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used. Requires version ``20.05.0-edge`` or later.
google.lifeSciences.subnetwork                 Define the name of the subnetwork to attach the instance to must be specified here, when the specified network is configured for custom subnet creation. The value is prefixed with `regions/subnetworks/` unless it contains a `/`, in which case it is assumed to be a fully specified subnetwork resource URL. Requires version ``21.03.0-edge`` or later.
google.lifeSciences.sshDaemon                  When ``true`` runs SSH daemon in the VM carrying out the job to which it's possible to connect for debugging purposes (default: ``false``).
google.lifeSciences.sshImage                   The container image used to run the SSH daemon (default: ``gcr.io/cloud-genomics-pipelines/tools``).
google.lifeSciences.keepAliveOnFailure         When ``true`` and a task complete with an unexpected exit status the associated computing node is kept up for 1 hour. This options implies ``sshDaemon=true`` (default: ``false``, requires Nextflow version ``21.06.0-edge`` or later).
google.storage.delayBetweenAttempts            Delay between download attempts from Google Storage (default `10 sec`, requires version ``21.06.0-edge`` or later).
google.storage.maxParallelTransfers            Max parallel upload/download transfer operations *per job* (default: ``4``, requires version ``21.06.0-edge`` or later).
google.storage.maxTransferAttempts             Max number of downloads attempts from Google Storage (default: `1`, requires version ``21.06.0-edge`` or later).
google.storage.parallelThreadCount             Defines the value for the option ``GSUtil:parallel_thread_count`` used by ``gsutil`` for transfer input and output data (default: ``1``, requires version ``21.06.0-edge`` or later).
google.storage.downloadMaxComponents           Defines the value for the option ``GSUtil:sliced_object_download_max_components`` used by ``gsutil`` for transfer input and output data (default: ``8``, requires version ``21.06.0-edge`` or later).
============================================== =================


Process definition
------------------
Processes can be defined as usual and by default the ``cpus`` and ``memory`` directives are used to instantiate a custom
machine type with the specified compute resources.  If ``memory`` is not specified, 1GB of memory is allocated per cpu.
A persistent disk will be created with size corresponding to the ``disk`` directive.  If ``disk`` is not specified, the
instance default is chosen to ensure reasonable I/O performance.

The process ``machineType`` directive may optionally be used to specify a predefined Google Compute Platform `machine type <https://cloud.google.com/compute/docs/machine-types>`_
If specified, this value overrides the ``cpus`` and ``memory`` directives.
If the ``cpus`` and ``memory`` directives are used, the values must comply with the allowed custom machine type `specifications <https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#specifications>`_ .  Extended memory is not directly supported, however high memory or cpu predefined
instances may be utilized using the ``machineType`` directive

Examples::

    process custom_resources_task {
        cpus 8
        memory '40 GB'
        disk '200 GB'

        """
        <Your script here>
        """
    }

    process predefined_resources_task {
        machineType 'n1-highmem-8'

        """
        <Your script here>
        """
    }

.. note:: This feature requires Nextflow 19.07.0 or later.

Pipeline execution
------------------

The pipeline can be launched either in a local computer or a cloud instance. Pipeline input data can be stored either
locally or in a Google Storage bucket.

The pipeline execution must specify a Google Storage bucket where the workflow's intermediate results are stored using
the ``-work-dir`` command line options. For example::

    nextflow run <script or project name> -work-dir gs://my-bucket/some/path


.. tip:: Any input data **not** stored in a Google Storage bucket will automatically be transferred to the
  pipeline work bucket. Use this feature with caution being careful to avoid unnecessary data transfers.

Preemptible instances
---------------------

Preemptible instances are supported adding the following setting in the Nextflow config file::

    google {
        lifeSciences.preemptible = true
    }

Since this type of virtual machines can be retired by the provider before the job completion, it is advisable
to add the following retry strategy to your config file to instruct Nextflow to automatically re-execute a job
if the virtual machine was terminated preemptively::

    process {
      errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
      maxRetries = 5
    }

.. note:: Preemptible instances have a `runtime limit <https://cloud.google.com/compute/docs/instances/preemptible>`_ of 24 hours.

.. tip:: For an exhaustive list of all possible error codes, please refer to the official Google LifeSciences `documentation <https://cloud.google.com/life-sciences/docs/troubleshooting#error_codes>`_.

Hybrid execution
----------------

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment
of hybrid workloads in which some jobs are executed in the local computer or local computing cluster and
some other jobs are offloaded to Google Pipelines service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the Google Pipelines *executor* only to a subset of processes in your workflow.
For example::


    process {
        withLabel: bigTask {
            executor = 'google-lifesciences'
            container = 'my/image:tag'
        }
    }

    google {
        project = 'your-project-id'
        zone = 'europe-west1-b'
    }


Then deploy the workflow execution using the ``-bucket-dir`` to specify a Google Storage path
for the jobs computed by the Google Pipeline service and, optionally, the ``-work-dir`` to
specify the local storage for the jobs computed locally::

    nextflow run <script or project name> -bucket-dir gs://my-bucket/some/path

.. warning:: The Google Storage path needs to contain at least sub-directory. Don't use only the
  bucket name e.g. ``gs://my-bucket``.

Quotas
------

Compute resources in Google Cloud are subject to `resource quotas <https://cloud.google.com/compute/quotas>`_ which may affect your ability to run pipelines at scale. You can request quota increases, and your quotas may automatically increase over time as you use the platform. In particular, GPU quotas are initially set to 0, so you must explicitly request a quota increase in order to use GPUs. Initially you can request an increase to 1 GPU at a time, and after one billing cycle you may be able to increase it further.

Limitations
-----------

* Currently it's not possible to specify a disk type different from the default one assigned
  by the service depending on the chosen instance type.



Troubleshooting
---------------

* Make sure to have enabled Compute Engine API, Life Sciences API and Cloud Storage Service in the
  `APIs & Services Dashboard <https://console.cloud.google.com/apis/dashboard>`_ page.

* Make sure to have enough compute resources to run your pipeline in your project
  `Quotas <https://console.cloud.google.com/iam-admin/quotas>`_ (i.e. Compute Engine CPUs,
  Compute Engine Persistent Disk, Compute Engine In-use IP addresses, etc).

* Make sure your security credentials allows you to access any Google Storage bucket
  where input data and temporary files are stored.

* Check the directory ``google/`` created in the task work directory (in the bucket storage) created
  when on job failure and containing useful information of the job execution. The creation
  can be enabled as default setting the option ``google.lifeSciences.debug = true`` in the
  Nextflow config file

* Enable the optional SSH daemon in the job VM using the option ``google.lifeSciences.sshDaemon = true``

* Make sure you are choosing a `location` where  `Cloud Life Sciences API is available <https://cloud.google.com/life-sciences/docs/concepts/locations>`_,
  and a `region` or `zone` where `Compute Engine is available <https://cloud.google.com/compute/docs/regions-zones/>`_.

.. _secrets-page:

*******
Secrets
*******


As of version ``21.09.0-edge``, Nextflow adds the built-in support for pipeline secrets to allow users to handle
and manage sensitive information for pipeline execution in a safe manner.

.. warning::
    This is a preview feature, therefore options and syntax may change in future release.

How it works
============

This feature allows decoupling the use secrets in your pipelines from the pipeline code and configuration files.
Secrets are instead managed by Nextflow and store separately into a local store only accessible to the secrets
owner.

When the pipeline execution is launched Nextflow inject the secrets in pipeline jobs without leaking them
into temporary execution files. The secrets are accessible into the job command via environment variables.

.. note::
  This feature needs to be enabled by settings the following environment variable in the launching environment::

        export NXF_ENABLE_SECRETS=true


Command line
============

When enabling this feature Nextflow provides a new command named ``secrets``. This command allow four simple
operations:

===================== =====================
Operation               Description
===================== =====================
``list``                List secrets available in the current store e.g. ``nextflow secrets list``.
``get``                 Allows retrieving a secret value e.g. ``nextflow secrets get -n FOO``.
``put``                 Allows creating creating a new secret or overriding an existing one e.g. ``nextflow secrets put -n FOO -v "Hello world"``
``delete``              Allows deleting an existing secret e.g. ``nextflow secrets delete -n FOO``.
===================== =====================

Configuration file
==================

Once create the secrets can be used in the pipeline configuration file as implicit variables using the ``secrets`` scope::

    aws {
      accessKey = secrets.MY_ACCESS_KEY
      secretKey = secrets.MY_SECRET_KEY
    }

The above above snippet access the secrets ``MY_ACCESS_KEY`` and ``MY_SECRET_KEY`` previously and assign them to
the corresponding AWS credentials settings.

.. note::
    Secrets **cannot** be assigned to pipeline parameters. 


Process secrets
===============

Secrets can be access by pipeline processes by using the `secret` directive. For example::

    process someJob {
        secret 'MY_ACCESS_KEY'
        secret 'MY_SECRET_KEY'

        """
        your_command --access $MY_ACCESS_KEY --secret $MY_SECRET_KEY
        """
    }

The above snippet runs a command in with the variables ``MY_ACCESS_KEY`` and ``MY_SECRET_KEY`` are injected in the
process execution environment holding the values defines in the secret store.

.. warning::
    This feature is only available when using the local or batch scheduler executions e.g. Slurm, Grid Engine, etc.
.. _cli-page:

*****************************
Command line interface (CLI)
*****************************

`Nextflow` provides a robust command line interface for the management and 
execution pipelines. The top-level interface consists of two aspects, 
*options* and *commands*.

Here's what you'll see at the top-level upon invoking the Nextflow CLI. ::


    $ nextflow
    Usage: nextflow [options] COMMAND [arg...]


.. _cli-options:

Options
============

The top-level options are meant to be invoked in relation to the core 
Nextflow application and are applied to all commands. For options 
specific to any command, refer the CLI Commands section.

.. note::
  Nextflow options only use single dash prefix e.g. ``-foo``. Do not confuse
  double dash notation e.g. ``--foo`` that is instead used for
  :ref:`Pipeline parameters <cli-params>`.

An overview of the top-level options. ::


    $ nextflow
    Usage: nextflow [options] COMMAND [arg...]

    Options:
    -C
        Use the specified configuration file(s) overriding any defaults
    -D
        Set JVM properties
    -bg
        Execute nextflow in background
    -c, -config
        Add the specified file to configuration set
    -d, -dockerize
        Launch nextflow via Docker (experimental)
    -h
        Print this help
    -log
        Set nextflow log file path
    -q, -quiet
        Do not print information messages
    -syslog
        Send logs to syslog server (eg. localhost:514)
    -v, -version
        Print the program version

    Commands...

---------------------------
Hard configuration override
---------------------------

Use the specified configuration file(s) overriding any defaults.

**Usage** ::

   $ nextflow -C my.config COMMAND [arg...]


**Description**

The ``-C`` option is used to override *all* settings specified in the default config file. 
For soft override, please refer the ``-c`` option.


**Examples**


- Override **any** default configuration with a custom configuration file. ::
    
  $ nextflow -C my.config run nextflow-io/hello


--------------------
JVM properties
--------------------

Set JVM properties.

**Usage**
::

   $ nextflow -Dkey=value COMMAND [arg...]

**Description**

This options allows the definition of custom Java system properties that can be used to 
properly configure or fine tuning the JVM instance used by the Nextflow runtime.
 
For specifying other JVM level options, please refer to the :ref:`config-env-vars` section.


**Examples**

Add `JVM properties` to the invoked pipeline. ::
    
    $ nextflow -Dfile.encoding=UTF-8 run nextflow-io/hello


-----------------------------
Execution as a background job
-----------------------------

Execute ``nextflow`` in the background.

**Usage**
::

   $ nextflow -bg COMMAND [arg...]

**Description**

The ``-bg`` option is used to invoke the nextflow execution in the background and allows 
the user to continue interacting with the terminal. This option is similar to ``nohup`` in 
behavior.


**Examples**

Invoke any execution as a background job. ::
    
    $ nextflow -bg run nextflow-io/hello 


---------------------------
Soft configuration override
---------------------------

Add the specified file to configuration set.

**Usage**

::

   $ nextflow -c nxf.config COMMAND [arg...]


**Description**

The ``-c`` option is used to append a new configuration to the default configuration. 
The ``-c`` option allows us to update the config in an additive manner. For **hard override**, 
refer the ``-C`` option.


**Examples**

Update **some** fields of the default config for any pipeline. ::

  $ nextflow -c nxf.config run nextflow-io/hello


-----------------------
Docker driven execution
-----------------------

Launch Nextflow via Docker.


.. warning::
    This is an experimental unsupported feature.

**Usage**
::

   $ nextflow -dockerize COMMAND [arg...]


**Description**

The ``-dockerize`` option is used to invoke the execution of **Nextflow** within a Docker container
itself without installing a Java VM in the hosting environment.

Note, this option is *not* needed to run containerised pipeline jobs. For invoking a pipeline with the ``docker`` profile or executor,
please to refer the ``-with-docker`` options the ``run`` command.

**Examples**

Invoke ``nextflow`` as a docker container to execute a pipeline. ::

   $ nextflow -dockerize run nextflow-io/hello


--------------------
Help
--------------------

Print the help message.

**Usage**
::

   $ nextflow -h

**Description**

The ``-h`` option prints out the overview of the CLI interface and enumerates the top-level *options* 
and *commands*.


--------------------
Execution logs
--------------------

Sets the path of the nextflow log file.

**Usage**
::

   $ nextflow -log custom.log COMMAND [arg...]


**Description**

The ``log`` option takes a path of the new log file which to be used instead of the 
default ``.nextflow.log`` or to save logs files to another directory.


**Examples**

Save all execution logs to the custom ``/var/log/nextflow.log`` file. ::

   $ nextflow -log /var/log/nextflow.log run nextflow-io/hello



--------------------
Quiet execution
--------------------

Disable the printing of information to the terminal.

**Usage**
::

    $ nextflow -q COMMAND [arg...]

**Description**

The ``-q`` option suppresses the banner, process related info and exits once the 
execution is completed. Please note that it does not affect any explicit print 
statement within a pipeline.


**Examples**

Invoke the pipeline execution without the banner and pipeline information. ::

   $ nextflow -q run nextflow-io/hello


---------------------------
Logging to a syslog server
---------------------------


Send logs to `Syslog <https://en.wikipedia.org/wiki/Syslog>`_ server endpoint.

**Usage** ::

    $ nextflow -syslog localhost:1234 COMMAND [arg...]


**Description**

The ``-syslog`` option is used to send logs to a `Syslog` logging server at the specified endpoint.


**Examples**

Send the logs to a `Syslog` server at specific endpoint. ::

    $ nextflow -syslog localhost:1234 run nextflow-io/hello


--------------------
Version
--------------------

Print the Nextflow version information.

**Usage**

::

    $ nextflow -v


**Description**

The ``-v`` option prints out information about *Nextflow* such as the version and build. 
The ``-version`` option in addition prints out the citation reference and official website.

**Examples**

- The short version. ::

      $ nextflow -v
      nextflow version 20.07.1.5412


- The full version info with citation and website link. ::

      $ nextflow -version
      N E X T F L O W
      version 20.07.1 build 5412
      created 24-07-2020 15:18 UTC (20:48 IDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


.. _cli-commands:

Commands
============

An overview of the Nextflow top-level commands. ::


    $ nextflow

    Usage: nextflow [options] COMMAND [arg...]
    
    Options...

    Commands:
    clean         Clean up project cache and work directories
    clone         Clone a project into a folder
    config        Print a project configuration
    console       Launch Nextflow interactive console
    drop          Delete the local copy of a project
    help          Print the usage help for a command
    info          Print project and system runtime information
    kuberun       Execute a workflow in a Kubernetes cluster (experimental)
    list          List all downloaded projects
    log           Print executions log and runtime info
    pull          Download or update a project
    run           Execute a pipeline project
    self-update   Update nextflow runtime to the latest available version
    view          View project script file(s)

--------------------
clean
--------------------

Clean up *cache* and *work* directories.

**Usage** ::

    $ nextflow clean [run_name|session_id] [options]


**Description**

Upon invocation within a directory, ``nextflow`` creates a project specific ``.nextflow.log`` 
file, ``.nextflow`` cache directory as well as a ``work`` directory. The ``clean`` command is 
designed to facilitate removal of these files from previous executions. 
A list of of run names and session ids can be generated by invoking ``nextflow log -q``.


**Options**

+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -after                    |            | Clean up runs executed *after* the specified one.                              |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -before                   |            | Clean up runs executed *before* the specified one.                             |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -but                      |            | Clean up all runs *except* the specified one.                                  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -dry-run, -n              |   false    | Print names of files to be removed without deleting them.                      | 
+---------------------------+------------+--------------------------------------------------------------------------------+
| -force, -f                |   false    | Force clean command.                                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |   false    | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -keep-logs, -k            |   false    | Removes only temporary files but retains execution log entries and metadata.   |                                           
+---------------------------+------------+--------------------------------------------------------------------------------+
| -quiet, -q                |   false    | Do not print names of files removed.                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

Dry run to remove work directories for the run name ``boring_euler``::

   $ nextflow clean boring_euler -n

   Would remove work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Would remove work/3f/70944c7a549b6221e1ccc7b4b21b62
   Would remove work/0e/2ebdba85f76f6068b21a1bcbf10cab

Remove work directories for the run name ``boring_euler``. ::

   $ nextflow clean boring_euler -f

   Removed work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
   Removed work/3f/70944c7a549b6221e1ccc7b4b21b62
   Removed work/0e/2ebdba85f76f6068b21a1bcbf10cab


Remove the execution entries *except* for a specific execution. ::

    $ nextflow clean -but tiny_leavitt -f

    Removed work/1f/f1ea9158fb23b53d5083953121d6b6
    Removed work/bf/334115deec60929dc18edf0010032a
    Removed work/a3/06521d75da296d4dd7f4f8caaddad8

Dry run to remove the execution data *before* a specific execution. ::

   $ nextflow clean -before tiny_leavitt -n

   Would remove work/5d/ad76f7b7ab3500cf616814ef644b61
   Would remove work/c4/69a82b080a477612ba8d8e4c27b579
   Would remove work/be/a4fa2aa38f76fd324958c81c2e4603
   Would remove work/54/39116773891c47a91e3c1733aad4de


Dry run to remove the execution data *after* a specific execution. ::

   $ nextflow clean -after focused_payne -n

   Would remove work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove work/bf/334115deec60929dc18edf0010032a
   Would remove work/a3/06521d75da296d4dd7f4f8caaddad8


Dry run to remove the temporary execution data for a specific execution, while keeping the log files. ::

   $ nextflow clean -keep-logs tiny_leavitt -n

   Would remove temp files from work/1f/f1ea9158fb23b53d5083953121d6b6
   Would remove temp files from work/bf/334115deec60929dc18edf0010032a
   Would remove temp files from work/a3/06521d75da296d4dd7f4f8caaddad8


--------------------
clone         
--------------------

Clone a remote project into a folder.


**Usage**


::

    $ nextflow clone [options] [project]


**Description**


The ``clone`` command downloads a pipeline from a Git-hosting platform into the *current directory* 
and modifies it accordingly. For downloading a pipeline into the global cache ``~/.nextflow/assets``,
please refer to the ``nextflow pull`` command.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -hub                      |  github    | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -r                        |  master    | Revision to clone - It can be a git ``branch``, ``tag`` or ``revision number`` |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -user                     |            | Private repository user name                                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**
Clone the latest revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello
    nextflow-io/hello cloned to: hello


Clone a specific revision of a pipeline. ::

    $ nextflow clone nextflow-io/hello -r v1.1
    nextflow-io/hello cloned to: hello


--------------------
config        
--------------------


Print the resolved pipeline configuration.

**Usage** ::

    $ nextflow config [options]


**Description**

The ``config`` command is used for printing the project's configuration i.e. the ``nextflow.config`` 
and is especially useful for understanding the resolved profiles and parameters that Nextflow will use 
run a pipeline. For in-depth information, please refer the :ref:`config-profiles` section.

**Options**

+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -flat                     |  false     | Print config using flat notation.                                              |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -profile                  |            | Choose a configuration profile.                                                |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -properties               |  false     | Print config using Java properties notation.                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -show-profiles, -a        |  false     | Show all configuration profiles.                                               |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -sort                     |  false     | Sort config attributes.                                                        |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

Print out the inferred config using a the default group key-value notation. ::

   $ nextflow config

   docker {
      enabled = true
   }

   process {
      executor = 'local'
   }

Print out the config using a flat notation. ::

   $ nextflow config -flat

   docker.enabled = true
   process.executor = 'local'


Print out the config using the Java properties notation. ::

   $ nextflow config -properties

   docker.enabled = true
   process.executor = local


Print out all profiles from the project's configuration. ::

   $ nextflow config -show-profiles

   docker {
      enabled = true
   }

   profiles {
      standard {
         process {
            executor = 'local'
         }
      }
      cloud {
         process {
            executor = 'cirrus'
            container = 'cbcrg/imagex'
         }
      }
   }

--------------------
console       
--------------------

Launch the *Nextflow* interactive console.

**Usage** ::

    $ nextflow console


**Description**

The ``console`` command is a wrapper over the Groovy *console* and provides a Graphic User 
Interface (GUI) and an interactive REPL (Read-Eval-Print-Loop) for quick experimentation.


**Options**

None available


**Examples**


Launch the ``console`` GUI. ::

  $ nextflow console


--------------------
drop          
--------------------

Delete the local copy of a project.


**Usage**

::

    $ nextflow drop [options] [project]


**Description**


The ``drop`` command is used to remove the projects which have been downloaded into the 
global cache. Please refer the ``list`` command for generating a list of downloaded pipelines.

**Options**

+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -f                        |  false     | Delete the repository without taking care of local changes.                    |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**


Drop the ``nextflow-io/hello`` project. ::

  $ nextflow drop nextflow-io/hello


Forcefully drop the ``nextflow-io/hello`` pipeline, ignoring any local changes. ::

  $ nextflow drop nextflow-io/hello -f


--------------------
help          
--------------------

Print the top-level help or specific help for a command.

**Usage**


::

    $ nextflow help [options] [command]


**Description**

The ``help`` command prints out the overview of the CLI interface and enumerates the top-level 
*options* and *commands*. Note that this command is equivalent to simply invoking ``nextflow`` 
at the command line.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

Invoke the ``help`` option for the ``drop`` command. ::

     $ nextflow help drop
 
     Delete the local copy of a project
     Usage: drop [options] name of the project to drop
        Options:
          -f
               Delete the repository without taking care of local changes
               Default: false
          -h, -help
               Print the command usage
               Default: false


--------------------
info          
--------------------


Print project or system runtime information.


**Usage**


::

    $ nextflow info [options] [project]


**Description**

The ``info`` command prints out the nextflow runtime information about the hardware as 
well as the software versions of the `Nextflow version and build`, `Operating System`
and `Groovy and Java runtime`. It can also be used to display information about a
specific project.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -check-updates, -u        |  false     | Check for remote updates.                                                      |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -d                        |  false     | Show detailed information.                                                     |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -o                        |  text      | Output format, either ``text``, ``json`` or ``yaml``.                          |
+---------------------------+------------+--------------------------------------------------------------------------------+



**Examples**

Display Nextflow runtime and system info::

    $ nextflow info

      Version: 20.07.1 build 5412
      Created: 24-07-2020 15:18 UTC (20:48 IDT)
      System: Mac OS X 10.15.6
      Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
      Encoding: UTF-8 (UTF-8)

Display information about a specific project::

    $ nextflow info nextflow-io/hello

      project name: nextflow-io/hello
      repository  : https://github.com/nextflow-io/hello
      local path  : /Users/evanfloden/.nextflow/assets/nextflow-io/hello
      main script : main.nf
      revisions   : 
      * master (default)
        mybranch
        testing
        v1.1 [t]
        v1.2 [t]


--------------------
kuberun       
--------------------

Deploy Nextflow into a Kubernetes cluster (experimental)

**Usage**

::

    $ nextflow kuberun [options] [project]


**Description**

The ``kuberun`` command builds upon the ``run`` command and offers a deep integration with 
the Kubernetes execution environment. This command deploys the Nextflow runtime as a Kubernetes 
pod and assumes that you've already installed the ``kubectl`` CLI. The ``kuberun`` command 
does not allow the execution of **local** Nextflow scripts. For more information please refer 
the :ref:`k8s-page` section.


**Options**


+---------------------------+-------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default     | Description                                                                    |
+===========================+=============+================================================================================+
| -E                        | false       | Exports all current system environment.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -ansi-log                 |             | Enable/disable ANSI console logging.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -bucket-dir               |             | Remote bucket where intermediate result files are stored.                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -cache                    |             | Enable/disable processes caching.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dsl2                     | false       | Execute the workflow using DSL2 syntax.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-channels            |             | Dump channels for debugging purpose.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-hashes              | false       | Dump task hash keys for debugging purpose.                                     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -e.                       | {}          | Add the specified variable to execution environment. Syntax: ``-e.key=value``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -entry                    |             | Entry workflow name to be executed.                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -h, -help                 | false       | Print the command usage.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -hub                      | github      | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -latest                   | false       | Pull latest changes before run.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -lib                      |             | Library extension path.                                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -name                     |             | Assign a mnemonic name to the a pipeline run.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -n, -namespace            |             | Specify the K8s namespace to use.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -offline                  | false       | Do not check for remote project updates.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -params-file              |             | Load script parameters from a JSON/YAML file.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -pod-image                |             | Specify the container image for the Nextflow pod.                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -process.                 | {}          | Set process options. Syntax ``-process.key=value``                             |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -profile                  |             | Choose a configuration profile.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -qs, -queue-size          |             | Max number of processes that can be executed in parallel by each executor.     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -resume                   |             | Execute the script using the cached results, useful to continue executions that|
|                           |             | was stopped by an error.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -r, -revision             |             | Revision of the project to run (either a git branch, tag or commit SHA number) |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -stub, -stub-run          |             | Execute the workflow replacing process scripts with command stubs.             |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -test                     |             | Test a script function with the name specified.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -user                     |             | Private repository user name.                                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -v, -volume-mount         |             | Volume claim mounts eg. ``my-pvc:/mnt/path``                                   |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-conda               |             | Use the specified Conda environment package or                                 |
|                           |             | file (must end with ``.yml|.yaml``)                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-dag                 | dag.dot     | Create pipeline DAG file.                                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-docker              |             | Enable process execution in a Docker container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -N, -with-notification    |             | Send a notification email on workflow completion to the specified recipients.  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-podman              |             | Enable process execution in a Podman container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-report              | report.html | Create processes execution html report.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-singularity         |             | Enable process execution in a Singularity container.                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-timeline            |timeline.html| Create processes execution timeline file.                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-tower               |             | Monitor workflow execution with Seqera Tower service.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-trace               | trace.txt   | Create processes execution tracing file.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-weblog              |             | Send workflow status messages via HTTP to target URL.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-docker           | false       | Disable process execution with Docker.                                         |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-podman           |             | Disable process execution in a Podman container.                               |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -w, -work-dir             | work        | Directory where intermediate result files are stored.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+



**Examples**

Execute a pipeline into a Kubernetes cluster. ::

     $ nextflow kuberun nextflow-io/hello 


--------------------
list          
--------------------

List all downloaded projects.

**Usage**

::

    $ nextflow list [options]



**Description**


The ``list`` commands prints a list of the projects which are already downloaded into the global cache ``~/.nextflow/assets``.


**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

List the downloaded pipelines. ::

    $ nextflow list

    nextflow-io/hello
    nextflow-hub/fastqc


--------------------
log           
--------------------

Print the execution history and log information.

**Usage** ::

    $ nextflow log [options] [run_name | session_id]


**Description**

The ``log`` command is used to query the execution metadata associated with pipelines executed 
by Nextflow. The list of executed pipelines can be generated by issuing ``nextflow log`` at the terminal. 
Instead of run name, it's also possible to use a session id. Moreover, this command contains multiple options 
to facilitate the queries and is especially useful while debugging a pipeline and while inspecting pipeline 
execution metadata.


**Options**

+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -after                    |            | Show log entries for runs executed *after* the specified one.                  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -before                   |            | Show log entries for runs executed *before* the specified one.                 |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -but                      |            | Show log entries for runs executed *but* the specified one.                    |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -filter, -F               |            | Filter log entires by a custom expression                                      |
|                           |            | e.g. ``process =~ /foo.*/ && status == 'COMPLETED'``                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -list-fields, -l          |  false     | Show all available fields.                                                     |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -quiet                    |  false     | Show only run names.                                                           |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -s                        |            | Character used to separate column values                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -template, -t             |            | Text template used to each record in the log.                                  |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

Listing the execution logs of previous invocations of all pipelines in a project. ::

    $ nextflow log

    TIMESTAMP          	DURATION	RUN NAME     	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2020-10-07 11:52:24	2.1s    	focused_payne	OK    	96eb04d6a4 	af6adaaa-ad4f-48a2-9f6a-b121e789adf5	nextflow run nextflow-io/hello -r master
    2020-10-07 11:53:00	3.1s    	tiny_leavitt 	OK    	e3b475a61b 	4d3b95c5-4385-42b6-b430-c865a70d56a4	nextflow run ./tutorial.nf
    2020-10-07 11:53:29	2.5s    	boring_euler 	OK    	e3b475a61b 	a6276975-7173-4208-ae09-ab9d6dce8737	nextflow run tutorial.nf


Listing only the *run names* of the execution logs of all pipelines invocations in a project. ::

    $ nextflow log -quiet

    focused_payne
    tiny_leavitt
    boring_euler

List the execution entries *only* a specific execution. ::

   $ nextflow log tiny_leavitt

   work/1f/f1ea9158fb23b53d5083953121d6b6
   work/bf/334115deec60929dc18edf0010032a
   work/a3/06521d75da296d4dd7f4f8caaddad8


List the execution entries *after* a specific execution. ::

    $ nextflow log -after tiny_leavitt

    work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
    work/3f/70944c7a549b6221e1ccc7b4b21b62
    work/0e/2ebdba85f76f6068b21a1bcbf10cab

List the execution entries *before* a specific execution. ::

    $ nextflow log -before tiny_leavitt

    work/5d/ad76f7b7ab3500cf616814ef644b61
    work/c4/69a82b080a477612ba8d8e4c27b579
    work/be/a4fa2aa38f76fd324958c81c2e4603
    work/54/39116773891c47a91e3c1733aad4de

List the execution entries *except* for a specific execution. ::

   $ nextflow log -but tiny_leavitt

    work/5d/ad76f7b7ab3500cf616814ef644b61
    work/c4/69a82b080a477612ba8d8e4c27b579
    work/be/a4fa2aa38f76fd324958c81c2e4603
    work/54/39116773891c47a91e3c1733aad4de

Filter specific fields from the execution log of a process. ::

    $ nextflow log tiny_leavitt -f 'process,exit,hash,duration'

    splitLetters	0	1f/f1ea91	112ms
    convertToUpper	0	bf/334115	144ms
    convertToUpper	0	a3/06521d	139ms

Filter fields from the execution log of a process based on a criteria. ::

    $ nextflow log tiny_leavitt -F 'process =~ /splitLetters/'

    work/1f/f1ea9158fb23b53d5083953121d6b6

--------------------
pull          
--------------------

Download or update a project.

**Usage** ::

    $ nextflow pull [options] [project]


**Description**


The ``pull`` command downloads a pipeline from a Git-hosting platform into the global cache ``~/.nextflow/assets`` 
and modifies it accordingly. For downloading a pipeline into a local directory, please refer to the ``nextflow clone`` command.


**Options**

+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -all                      |  false     | Update all downloaded projects.                                                |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -hub                      |  github    | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -r                        |            | Revision to run (either a git ``branch``, ``tag`` or commit ``SHA`` number).   |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -user                     |            | Private repository user name                                                   |
+---------------------------+------------+--------------------------------------------------------------------------------+


**Examples**

Download a new pipeline or pull the latest revision for a specific project. ::

    $ nextflow pull nextflow-io/hello

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]

Pull the latest revision for all downloaded projects. ::

    $ nextflow pull -all

    Checking nextflow-io/hello ...
    done - revision: 96eb04d6a4 [master]
    Checking nextflow-hub/fastqc ...
    done - revision: 087659b18e [master]

Download a specific revision of a new project or pull the latest revision for a specific project. ::

    $ nextflow pull nextflow-io/hello -r v1.1

    Checking nextflow-io/hello ...
    checkout-out at AnyObjectId[1c3e9e7404127514d69369cd87f8036830f5cf64] - revision: 1c3e9e7404 [v1.1]


--------------------
run           
--------------------

Execute a pipeline.

**Usage**

::

    $ nextflow run [options] [project]


**Description**


The ``run`` command is used to initiate the execution of the a pipeline script or
download a pipeline project. Along with serving the purpose of script execution, this command
facilitates rapid iterations, inspections of any pipeline as well as debugging.


**Options**

+---------------------------+-------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default     | Description                                                                    |
+===========================+=============+================================================================================+
| -E                        |  false      | Exports all current system environment.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -ansi-log                 |             | Enable/disable ANSI console logging.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -bucket-dir               |             | Remote bucket where intermediate result files are stored.                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -cache                    |             | Enable/disable processes caching.                                              |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dsl2                     | false       | Execute the workflow using DSL2 syntax.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-channels            |             | Dump channels for debugging purpose.                                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -dump-hashes              | false       | Dump task hash keys for debugging purpose.                                     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -e.                       | {}          | Add the specified variable to execution environment. Syntax: ``-e.key=value``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -entry                    |             | Entry workflow name to be executed.                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -h, -help                 | false       | Print the command usage.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -hub                      | github      | Service hub where the project is hosted. Options: ``gitlab`` or ``bitbucket``  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -latest                   | false       | Pull latest changes before run.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -lib                      |             | Library extension path.                                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -main-script              | main.nf     | The script file to be executed when launching a project directory or repository|
|                           |             | (requires version 20.09.1-edge or later).                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -name                     |             | Assign a mnemonic name to the a pipeline run.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -offline                  | false       | Do not check for remote project updates.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -params-file              |             | Load script parameters from a JSON/YAML file.                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -plugins                  |             | Comma separated list of plugin ids to be applied in the pipeline execution.    |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -process.                 | {}          | Set process options. Syntax ``-process.key=value``                             |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -profile                  |             | Choose a configuration profile.                                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -qs, -queue-size          |             | Max number of processes that can be executed in parallel by each executor.     |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -resume                   |             | Execute the script using the cached results, useful to continue executions that|
|                           |             | was stopped by an error.                                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -r, -revision             |             | Revision of the project to run                                                 |
|                           |             | (either a git ``branch``, ``tag`` or commit ``SHA`` number).                   |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -test                     |             | Test a script function with the name specified.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -user                     |             | Private repository user name.                                                  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-conda               |             | Use the specified Conda environment package or                                 |
|                           |             | file (must end with ``.yml|.yaml``)                                            |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-dag                 | dag.dot     | Create pipeline DAG file.                                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-docker              |             | Enable process execution in a Docker container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -N, -with-notification    |             | Send a notification email on workflow completion to the specified recipients.  |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-podman              |             | Enable process execution in a Podman container.                                |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-report              | report.html | Create processes execution html report.                                        |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-singularity         |             | Enable process execution in a Singularity container.                           |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-timeline            |timeline.html| Create processes execution timeline file.                                      |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-tower               |             | Monitor workflow execution with Seqera Tower service.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-trace               | trace.txt   | Create processes execution tracing file.                                       |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -with-weblog              |             | Send workflow status messages via HTTP to target URL.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-docker           | false       | Disable process execution with Docker.                                         |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -without-podman           |             | Disable process execution in a Podman container.                               |
+---------------------------+-------------+--------------------------------------------------------------------------------+
| -w, -work-dir             | work        | Directory where intermediate result files are stored.                          |
+---------------------------+-------------+--------------------------------------------------------------------------------+


**Examples**

- Run a specific revision of a downloaded pipeline. ::

    $ nextflow run nextflow-io/hello -r v1.1

    N E X T F L O W  ~  version 20.07.1
    Launching `nextflow-io/hello` [grave_cajal] - revision: 1c3e9e7404 [v1.1]


- Choose a `profile` for running the project. Assumes that a profile named ``docker`` has already been defined in the config file. ::

    $ nextflow run main.nf -profile docker


- Invoke the pipeline execution and generate the summary HTML report. For more information on the metrics, please refer the :ref:`perfanalysis-page` section::

    $ nextflow run main.nf -with-report


- Invoke the nextflow pipeline execution with a custom queue size. By default, the value of **queue-size** is the same as the number of available CPUs. ::

    $ nextflow run nextflow-io/hello -qs 4


- Execute the pipeline with DSL-2 syntax. ::

    $ nextflow run nextflow-io/hello -dsl2


- Invoke the pipeline with a specific workflow as the entry-point, this option is meant to be used with DSL-2.
For more information on DSL-2, please

refer the :ref:`dsl2-page`::

   $ nextflow run main.nf -entry workflow_A


- Invoke the nextflow pipeline execution with the integrated monitoring dashboard Tower. For more information, please refer to the `tower.nf <https://tower.nf>`_ website. ::

    $ nextflow run nextflow-io/hello -with-tower
 

--------------------
self-update   
--------------------




Update the nextflow runtime to the latest available version.


**Usage**

::

    $ nextflow self-update


**Description**

The ``self-update`` command directs the ``nextflow`` cli to update itself to the latest stable release.


**Examples**

Update Nextflow. ::

    $ nextflow self-update

          N E X T F L O W
          version 20.07.1 build 5412
          created 24-07-2020 15:18 UTC (20:48 IDT)
          cite doi:10.1038/nbt.3820
          http://nextflow.io


    Nextflow installation completed. Please note:
    - the executable file `nextflow` has been created in the folder: /usr/local/bin


--------------------
view          
--------------------

View a projects script file(s).

**Usage**

::

    $ nextflow view [options] [project]


**Description**


The ``view`` command is used to inspect the pipelines which are already stored in the global nextflow cache. 
For downloading a pipeline into the global cache ``~/.nextflow/assets``, please refer to the ``pull`` command.

**Options**


+---------------------------+------------+--------------------------------------------------------------------------------+
| Name, shorthand (if any)  | Default    | Description                                                                    | 
+===========================+============+================================================================================+
| -help, -h                 |  false     | Print the command usage.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -l                        |  false     | List repository content.                                                       |
+---------------------------+------------+--------------------------------------------------------------------------------+
| -q                        |  false     | Hide header line.                                                              |
+---------------------------+------------+--------------------------------------------------------------------------------+

**Examples**


Viewing the contents of a downloaded pipeline. ::

   $ nextflow view nextflow-io/hello

   == content of file: .nextflow/assets/nextflow-io/hello/main.nf
   #!/usr/bin/env nextflow

   cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

    process sayHello {
      echo true
      input:
        val x from cheers
      script:
        """
        echo '$x world!'
        """
    }


Listing the folder structure of the downloaded pipeline. ::

   $ nextflow view -l nextflow-io/hello

   == content of path: .nextflow/assets/nextflow-io/hello
   LICENSE
   README.md
   nextflow.config
   .gitignore
   circle.yml
   foo.nf
   .git
   .travis.yml
   main.nf


Viewing the contents of a downloaded pipeline without omitting the header. ::

   $ nextflow view -q nextflow-io/hello

   #!/usr/bin/env nextflow

   cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

    process sayHello {
      echo true
      input:
        val x from cheers
      script:
        """
        echo '$x world!'
        """
    }


.. _cli-params:

Pipeline parameters
====================

Pipeline script can use an arbitrary number of parameters that can be overridden either
using the command line or the Nextflow configuration file. Any script parameter can be specified
on the command line prefixing the parameter name with double dash characters e.g.::

    nextflow run <my script> --foo Hello

Then, the parameter can be accessed in the pipeline script using the ``params.foo`` identifier.

.. note::
  When the parameter name is formatted using the `camelCase` notation e.g. ``fooBar``, a second parameter
  is created with the same value using the `kebab-case` notation e.g. ``foo-bar``, and the other way around.

.. warning::
  When a command line parameters includes one or more glob characters i.e. wildcards like ``*`` or ``?``,
  the parameter value needs to be enclosed with double-quote character to prevent Bash expansion and preserve
  the glob characters. For example::

        nextflow run <my script> --files "*.fasta"
.. _metadata-page:

***********************
Workflow introspection
***********************

.. _metadata-workflow:

Runtime metadata
----------------

The implicit ``workflow`` object allows you to access some workflow and runtime metadata in your Nextflow scripts.
For example::

    println "Project : $workflow.projectDir"
    println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println "Cmd line: $workflow.commandLine"
    println "Manifest's pipeline version: $workflow.manifest.version"


.. tip:: To shortcut the access to multiple ``workflow`` properties you can use the Groovy
    `with <http://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Object.html#with(groovy.lang.Closure)>`_ method.


The following table lists the properties that can be accessed on the ``workflow`` object:

=========================== ===========================
Name                        Description
=========================== ===========================
scriptId                    Project main script unique hash ID.
scriptName                  Project main script file name.
scriptFile                  Project main script file path.
repository                  Project repository Git remote URL.
commitId                    Git commit ID of the executed workflow repository.
revision                    Git branch/tag of the executed workflow repository.
projectDir                  Directory where the workflow project is stored in the computer.
launchDir                   Directory where the workflow execution has been launched.
workDir                     Workflow working directory.
homeDir                     User system home directory.
userName                    User system account name.
configFiles                 Configuration files used for the workflow execution.
container                   Docker image used to run workflow tasks. When more than one image is used
                            it returns a map object containing `[process name, image name]` pair entries.
containerEngine             Returns the name of the container engine (e.g. docker or singularity) or null
                            if no container engine is enabled. 
commandLine                 Command line as entered by the user to launch the workflow execution.
profile                     Used configuration profile.
runName                     Mnemonic name assigned to this execution instance.
sessionId                   Unique identifier (UUID) associated to current execution.
resume                      Returns ``true`` whenever the current instance is resumed from a previous execution.
stubRun                     Returns ``true`` whenever the current instance is a stub-run execution .
start                       Timestamp of workflow at execution start.
manifest                    Entries of the workflow manifest.
:sup:`✝` complete           Timestamp of workflow when execution is completed.
:sup:`✝` duration           Time elapsed to complete workflow execution.
:sup:`*` success            Reports if the execution completed successfully.
:sup:`*` exitStatus         Exit status of the task that caused the workflow execution to fail.
:sup:`*` errorMessage       Error message of the task that caused the workflow execution to fail.
:sup:`*` errorReport        Detailed error of the task that caused the workflow execution to fail.
=========================== ===========================

| Properties marked with a `✝` are accessible only in the workflow completion handler.
| Properties marked with a `*` are accessible only in the workflow completion and error handlers. See the `Completion handler`_ section for details.
|

.. _metadata-nextflow:

Nextflow metadata
-----------------

The implicit ``nextflow`` object allows you to access the metadata information of the Nextflow runtime.

=========================== ===========================
Name                        Description
=========================== ===========================
nextflow.version            Nextflow runtime version number.
nextflow.build              Nextflow runtime build number.
nextflow.timestamp          Nextflow runtime compile timestamp.
=========================== ===========================

The method ``nextflow.version.matches`` allows you to check if the Nextflow runtime satisfies the version
requirement eventually needed by your workflow script. The required version string can be prefixed with the usual
comparison operators eg ``>``, ``>=``, ``=``, etc. or postfixed with the ``+`` operator to specify a minimal version
requirement. For example::

    if( !nextflow.version.matches('0.22+') ) {
        println "This workflow requires Nextflow version 0.22 or greater -- You are running version $nextflow.version"
        exit 1
    }


.. _metadata-completion-handler:

Completion handler
------------------

Due to the asynchronous nature of Nextflow the termination of a script does not correspond to the termination
of the running workflow. Thus some information, only available on execution completion, needs to be accessed by
using an asynchronous handler.

The ``onComplete`` event handler is invoked by the framework when the workflow execution is completed. It allows one
to access the workflow termination status and other useful information. For example::

    workflow.onComplete {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }

If you want an e-mail notification on completion, check :ref:`mail-page`.

.. _metadata-error-handler:

Error handler
-------------

The ``onError`` event handler is invoked by Nextflow when a runtime or process error caused the pipeline execution to stop.
For example::

    workflow.onError {
        println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    }

.. note:: Both the ``onError`` and ``onComplete`` handlers are invoked when an error condition is encountered.
    However the first is called as soon as the error is raised, while the second just before the pipeline execution
    is going to terminate. When using the ``finish`` :ref:`process-page-error-strategy`, between the two there could be
    a significant time gap depending on the time required to complete any pending job.


Decoupling metadata
-----------------------

The workflow event handlers can be defined also in the ``nextflow.config`` file. This is useful to
decouple the handling of pipeline events from the main script logic.

When the event handlers are included in a configuration file the only difference is that the ``onComplete`` and
the ``onError`` closures have to be defined by using the assignment operator as shown below::

    workflow.onComplete = {
        // any workflow property can be used here
        println "Pipeline complete"
        println "Command line: $workflow.commandLine"
    }


    workflow.onError = {
        println "Oops .. something when wrong"
    }


.. note:: It is possible to define a workflow event handlers both in the pipeline script **and** in the
  configuration file.

.. _example-page:

********
Examples
********

Basic pipeline
--------------

This example shows a pipeline that is made of two processes. The first process receives a
`FASTA formatted <http://en.wikipedia.org/wiki/FASTA_format>`_ file and splits it into file chunks whose names start with
the prefix ``seq_``.

The process that follows, receives these files and it simply `reverses` their content by using the ``rev`` command line tool.

.. raw:: html

    <style type='text/css'>
    .gist {font-size:13px;line-height:18px;margin-bottom:20px;width:100%}
    .gist pre {font-family:Menlo,Monaco,'Bitstream Vera Sans Mono','Courier New',monospace !important}
    .gist-meta {font-family:Helvetica,Arial,sans-serif;font-size:13px !important}
    .gist-meta a {color:#26a !important;text-decoration:none}
    .gist-meta a:hover{color:#0e4071 !important}
    .gist .file-data { background-color: white; }
    </style>
    <script src="https://gist.github.com/pditommaso/92fea4525cd66c286904.js"></script>


In more detail:

* line 1: The script starts with a `shebang <http://en.wikipedia.org/wiki/Shebang_(Unix)>`_ declaration. This allows you
  to launch your pipeline, as any other Bash script

* line 3: Declares a pipeline parameter named ``params.in`` that is initialized with the value ``$HOME/sample.fa``.This value
  can be overridden when launching the pipeline, by simply adding the option ``--in <value>`` to the script command line

* line 5: Defines a variable ``sequences`` holding a reference for the file whose name is specified by the ``params.in``
  parameter

* line 6: Defines a variable ``SPLIT`` whose value is ``gcsplit`` when the script is executed on a Mac OSX or ``csplit``
  when it runs on Linux. This is the name of the tool that is used to split the file.

* lines 8-20: The process that splits the provided file.

* line 10: Opens the `input` declaration block. The lines following this clause are interpreted as input definitions.

* line 11: Defines the process input file. This file is received from the variable ``sequences`` and will be named ``input.fa``.

* line 13: Opens the `output` declaration block. Lines following this clause are interpreted as output definitions.

* line 14: Defines that the process output files whose names match the pattern ``seq_*``. These files are sent over the
  channel ``records``.

* lines 16-18: The actual script executed by the process to split the provided file.

* lines 22-33: Defines the second process, that receives the splits produced by the previous process and reverses their
  content.

* line 24: Opens the `input` declaration block. Lines following this clause are interpreted as input definitions.

* line 25: Defines the process input file. This file is received through the channel ``records``.

* line 27: Opens the `output` declaration block. Lines following this clause are interpreted as output definitions.

* line 28: The `standard output` of the executed script is declared as the process output. This output is sent over the
  channel ``result``.

* lines 30-32: The actual script executed by the process to `reverse` the content of the received files.

* line 35: Prints a `result` each time a new item is received on the ``result`` channel.


.. tip:: The above example can manage only a single file at a time. If you want to execute it for two (or more) different files
   you will need to launch it several times.

   It is possible to modify it in such a way that it can handle any number of input files, as shown below.


In order to make the above script able to handle any number of files simply replace `line 3` with the following line::

  sequences = Channel.fromPath(params.in)


By doing this the ``sequences`` variable is assigned to the channel created by the :ref:`channel-path` method. This
channel emits all the files that match the pattern specified by the parameter ``params.in``.

Given that you saved the script to a file named ``example.nf`` and you have a list of `FASTA` files in a folder
named ``dataset/``, you can execute it by entering this command::

  nextflow example.nf --in 'dataset/*.fa'


.. warning:: Make sure you enclose the ``dataset/*.fa`` parameter value in single-quotation characters,
  otherwise the Bash environment will expand the ``*`` symbol to the actual file names and the example won't work.

More examples
-------------

You can find at `this link <https://github.com/nextflow-io/examples>`_ a collection of examples introducing Nextflow
scripting.

Check also `Awesome Nextflow <https://github.com/nextflow-io/awesome-nextflow/>`_ for a list
of pipelines developed by the Nextflow community... _metrics-page:

*******
Metrics
*******

This section details how the resource usage from the :ref:`Execution report <execution-report>` are computed.


CPU Usage
=========

The plot reports how much CPU resources were used by each process.

.. image:: images/report-resource-cpu.png

Let's illustrate how this plot behaves with several examples.

In the first example, let's consider the simple use case in which a process performs one task of pure computation using one CPU. Then, you expect the `Raw Usage` tab to report 100%. If the task is distributed over, 2, 3, 4, `etc.` CPUs, then the `Raw Usage` will be 200%, 300%, 400%, `etc.` respectively. The `% Allocated` tab just rescales the raw value usage with respect to the number of CPUs set with the ``cpus`` directive (if not set with the directive, the number of CPUs is set to 1, thus showing the same values as in the `Raw Usage` tab). Using the program `stress <https://people.seas.harvard.edu/~apw/stress/>`_  as follows would report 100% in the `Raw Usage` tab  and 50% in the `% Allocated` tab since the process asked twice the number of CPUs needed by the process::

  #!/usr/bin/env nextflow

  process CpuUsageEx1 {
    
    cpus 2

    """
    stress -c 1 -t 10 # compute square-root of random numbers during 10s using 1 CPU
    """
  }


In the second example, some time will be spent performing pure computation and some time just waiting. Using the program `stress <https://people.seas.harvard.edu/~apw/stress/>`_  and `sleep` as follows would report 75% in the `Raw Usage` tab::



  #!/usr/bin/env nextflow

  process CpuUsageEx2 {
    
    cpus 1

    """
    stress -c 1 -t 10 # compute square-root of random numbers during 10s using 1 CPU
    stress -c 1 -t 5 # compute square-root of random numbers during 5s using 1 CPU
    sleep 5 # use no CPU during 5s
    """
  }


Indeed, the percentage of the CPU that this process got is a weighted average taking into account the percentage of the CPU and duration of each individual program over the job duration (a.k.a. elapsed real time, real time or wall time ) as follows:


.. math::
  \frac{ 100\% \times 10s + 100\% \times 5s + 0\% \times 5s }{10s+5s+5s} = 75\%


The third example is similar to the second one except that the pure computation stage is performed in a single step forked on 2 CPUs::



  #!/usr/bin/env nextflow

  process CpuUsageEx3 {
    
    cpus 2

    """
    stress -c 2 -t 10 # compute square-root of random numbers during 10s using 2 CPUs
    sleep 10 # use no CPU during 10s
    """
  }


The `Raw Usage` tab would report 100% in the `Raw Usage` tab:

.. math::
  \frac{ 200\% \times 10s }{10s+10s} = 100\%

The `% Allocated` tab would report 50%, however, it would not be relevant to change the ``cpus`` directive from 2 to 1 as the process really uses 2 CPUs at it peak load.


.. hint:: The `stress <https://people.seas.harvard.edu/~apw/stress/>`_ program can be installed with ``sudo apt-get install stress`` or ``sudo yum install stress`` depending on your linux distribution.


Memory Usage
============

The plot has three tabs showing the usage of the physical memory (RAM), the virtual memory (vmem) and the percentage of RAM used by the process with respect to what was set in the ``memory`` directive. The peak usage during the execution of the process is reported for both physical and virtual memories.

.. hint::
  To better understand the memory usage plot, it is important to know that:

  - the total amount of memory used be a processs is the `virtual memory (vmem)`. The `vmem` contains all memory areas whether they are in the physical memory (RAM), in the Swap space, on the disk or shared with other processes,
  - the `resident set size (RSS)` is the amount of space of `physical memory (RAM)` held by a process,
  - the relationship is: vmem :math:`\geq` RSS + Swap,
  - the ``memory`` directive sets the RAM requested by the process.


Let's illustrate how this plot behaves with one example which relies on two C programs. 

The first program just allocates a variable of 1 GiB:

.. code-block:: c
   :linenos:
   :emphasize-lines: 31,43

    #include <stdio.h>
    #include <stdlib.h>
    #include <sys/resource.h>
    #include <sys/time.h>
    #include <sys/types.h>
    #include <unistd.h>
    #include <time.h>

    /* Get vmem and rss usage from /proc/<pid>/statm */
    static int mem_used(pid_t pid, unsigned long* vmem, unsigned long* rss) {
        FILE* file;
        char path[40];
        unsigned int page_size;

        snprintf(path, 40, "/proc/%ld/statm", (long) pid);
        file = fopen(path, "r");
        // vmem and rss are the first values in the file
        fscanf(file, "%lu %lu", vmem, rss);
        // values in statm are in pages so to get bytes we need to know page size
        page_size = (unsigned) getpagesize();
        *vmem = *vmem * page_size;
        *rss = *rss * page_size;

        fclose(file);
        return 0;
    }

    int main(int argc, char **argv) {
        unsigned char *address;
        char input;
        size_t size = 1024*1024*1024;  // 1 GiB
        unsigned long i;
        unsigned long vmem = 0;
        unsigned long rss = 0;
        pid_t pid;

        pid = getpid();
        printf("Pid: %ld\n", (long) pid);

        mem_used(pid, &vmem, &rss);
        printf("VMEM: %lu RSS: %lu\n", vmem, rss);

        address = malloc(size);
        printf("Allocated %d Bytes of memory\n", (int) size);

        mem_used(pid, &vmem, &rss);
        printf("VMEM: %lu RSS: %lu\n", vmem, rss);

        // Leave time for nextflow to get information
        sleep(15);
        
        free(address);
        return 0;
    }



The second program allocates a variable of 1 GiB and fills it with data:

.. code-block:: c
   :linenos:
   :emphasize-lines: 31,43,49-53

    #include <stdio.h>
    #include <stdlib.h>
    #include <sys/resource.h>
    #include <sys/time.h>
    #include <sys/types.h>
    #include <unistd.h>
    #include <time.h>

    /* Get vmem and rss usage from /proc/<pid>/statm */
    static int mem_used(pid_t pid, unsigned long* vmem, unsigned long* rss) {
        FILE* file;
        char path[40];
        unsigned int page_size;

        snprintf(path, 40, "/proc/%ld/statm", (long) pid);
        file = fopen(path, "r");
        // vmem and rss are the first values in the file
        fscanf(file, "%lu %lu", vmem, rss);
        // values in statm are in pages so to get bytes we need to know page size
        page_size = (unsigned) getpagesize();
        *vmem = *vmem * page_size;
        *rss = *rss * page_size;

        fclose(file);
        return 0;
    }

    int main(int argc, char **argv) {
        unsigned char *address;
        char input;
        size_t size = 1024*1024*1024;  // 1 GiB
        unsigned long i;
        unsigned long vmem = 0;
        unsigned long rss = 0;
        pid_t pid;

        pid = getpid();
        printf("Pid: %ld\n", (long) pid);

        mem_used(pid, &vmem, &rss);
        printf("VMEM: %lu RSS: %lu\n", vmem, rss);

        address = malloc(size);
        printf("Allocated %d Bytes of memory\n", (int) size);

        mem_used(pid, &vmem, &rss);
        printf("VMEM: %lu RSS: %lu\n", vmem, rss);

        printf("Filling memory with data...");
        fflush(stdout);  
        for (i = 0; i < size; i++) {
            *(address + i) = 123;
        }

        mem_used(pid, &vmem, &rss);
        printf("\nVMEM: %lu RSS: %lu\n", vmem, rss);

        // Leave time for nextflow to get information
        sleep(15);
        
        free(address);
        return 0;
    }


The first and second programs are executed in ``foo`` and ``bar`` processes respectively as follows::

  #!/usr/bin/env nextflow

  process foo {
  
      memory '1.5 GB'
  
      """
      memory_vmem_1GiB_ram_0Gib
      """

  }

  process bar {

      memory '1.5 GB'

      """
      memory_vmem_1GiB_ram_1Gib
      """

  }

The `Virtual (RAM + Disk swap)` tab shows that both ``foo`` and ``bar`` processes use the same amount of virtual memory (~1 GiB):

.. image:: images/report-resource-memory-vmem.png

However, the `Physical (RAM)` tab shows that only the ``bar`` process uses ~1 GiB of RAM while ``foo`` process uses  ~0 GiB:

.. image:: images/report-resource-memory-ram.png


As expected, the `% RAM Allocated` tab shows that 0% of the resource set in the ``memory`` directive was used for ``foo`` process while 67% (= 1 / 1.5) of the resource were used for ``bar`` process:

.. image:: images/report-resource-memory-pctram.png

.. warning:: Binary unit are used to report memory raw values. This means that 1KB = :math:`1024` bytes, 1 MB = :math:`1024^2` bytes, 1 GB = :math:`1024^3` bytes, `etc.`


Job Duration
============

The plot has two tabs the job duration (a.k.a. elapsed real time, real time or wall time ) in the `Raw Usage` tag and the percentage  of requested time used in the `% Allocated` tab with respect to the duration set in the ``time`` directive of the process.

.. image:: images/report-resource-job-duration.png


I/O Usage
=========

The plot has two tabs showing how many data were read and/or written each process. For example, the following processes read and write 1GB and 256MB of data respectively::

  #!/usr/bin/env nextflow

    process io_read_write_1G {
      """
      dd if=/dev/zero of=/dev/null bs=1G count=1
      """
    }


    process io_read_write_256M {
      """
      dd if=/dev/zero of=/dev/null bs=256M count=1
      """
    }

`Read` tab:

.. image:: images/report-resource-io-read.png


`Write` tab:

.. image:: images/report-resource-io-write.png


.. _plugins-page:

*********
Plugins
*********

Main concepts
=============

Nextflow is based on a plugins system that allows extending core functionalities via pluggable components
that are download and installed at runtime.

Currently the following functionalities are implemented as plugin components and they make part of the
Nextflow *default* plugins:

* ``nf-amazon``: Support for Amazon cloud.
* ``nf-azure``: Support for Azure cloud.
* ``nf-console``: Implement Nextflow `REPL console <https://www.nextflow.io/blog/2015/introducing-nextflow-console.html>`_.
* ``nf-ga4gh``: Support `GA4GH APIs <https://www.ga4gh.org/>`_.
* ``nf-google``: Support for Google cloud.
* ``nf-tower``: Support for `Nextflow Tower <https://tower.nf>`_ platform.


Configuration
==============

Nextflow *defaults* plugins do not require any configuration, they are automatically installed on-demand when
the corresponding feature is requested by a Nextflow pipeline.

To use **non-default** plugins in your pipeline execution it's required to declared them in the Nextflow configuration file
listing each plugin identifiers as shown below::

    plugins {
      id 'nf-hello@0.1.0'
    }


.. note::
  The plugin identifier is composed by the plugin name, followed by the ``@`` separator and finally the plugin version.

Alternatively, plugins can be required using the command line option ``-plugins`` e.g.::

    nextflow run <PIPELINE NAME> -plugins nf-hello@0.1.0


.. note::
  When using the ``-plugins`` CLI option any plugin declaration in the Nextflow config file is ignored.
  Multiple plugin Ids can be specified separating them with a comma character.


Index
======

Nextflow resolves plugins download location through the `Plugins index <https://github.com/nextflow-io/plugins/>`_.
The index stores for each plugin the available version, the creation date, checksum and the link from where the plugin
file is downloaded.

To add a new plugin to the Index, create a pull request including the request plugin metadata.

.. tip::
  The `nf-hello plugin <https://github.com/nextflow-io/nf-hello>`_ repository provides an bare minimal code example for
  the implementation of a Nextflow plugin.
.. _mail-page:

***********************
Mail & Notifications
***********************

Mail message
-------------

The built-in function ``sendMail`` allows you to send a mail message from a workflow script.

.. _mail-basic:

Basic mail
==========

The mail attributes are specified as named parameters or providing an equivalent associative array as argument.
For example::

        sendMail( to: 'you@gmail.com',
                  subject: 'Catch up',
                  body: 'Hi, how are you!',
                  attach: '/some/path/attachment/file.txt' )

therefore this is equivalent to write::

        mail = [ to: 'you@gmail.com',
                 subject: 'Catch up',
                 body: 'Hi, how are you!',
                 attach: '/some/path/attachment/file.txt' ]

        sendMail(mail)


The following parameters can be specified:

================== ================
Name                Description
================== ================
to :sup:`*`         The mail target recipients.
cc :sup:`*`         The mail CC recipients.
bcc :sup:`*`        The mail BCC recipients.
from :sup:`*`       The mail sender address.
subject             The mail subject.
charset             The mail content charset (default: ``UTF-8``).
text                The mail plain text content.
body                The mail body content. It can be either plain text or HTML content.
type                The mail body mime type. If not specified it's automatically detected.
attach              Single file or a list of files to be included as mail attachments.
================== ================

`*` Multiple email addresses can be specified separating them with a comma.

.. _mail-advanced:

Advanced mail
=============

An second version of the ``sendMail`` allows a more idiomatic syntax::

    sendMail {
        to 'you@gmail.com'
        from 'me@gmail.com'
        attach '/some/path/attachment/file.txt'
        attach '/other/path/image.png'
        subject 'Catch up'

        '''
        Hi there,
        Look! Multi-lines
        mail content!
        '''
    }

The same attributes listed in the table in the previous section are allowed.

.. note:: When it terminates with a string expression it's implicitly interpreted as the mail body content, therefore
  the ``body`` parameter can be omitted as shown above.

.. tip:: To send an `alternative` mail message that includes either text and HTML content use both the ``text`` and ``body`` attributes.
  The first must be used for the plain text content, while the second for the rich HTML message.

.. _mail-attachments:

Mail attachments
================

When using the curly brackets syntax, the ``attach`` parameter can be repeated two or more times to include
multiple attachments in the mail message.

Moreover for each attachment it's possible to specify one or more of the following optional attributes:

================== ================
Name                Description
================== ================
contentId           Defines the `Content-ID` header field for the attachment.
disposition         Defines the `Content-Disposition` header field for the attachment.
fileName            Defines the `filename` parameter of the "Content-Disposition" header field.
description         Defines the `Content-Description` header field for the attachment.
================== ================

For example::

    sendMail {
        to 'you@dot.com'
        attach '/some/file.txt', fileName: 'manuscript.txt'
        attach '/other/image.png', disposition: 'inline'
        subject 'Sending documents'
        '''
        the mail body
        '''
    }

.. _mail-config:

Mail configuration
==================

If no mail server configuration is provided, Nextflow tries to send the email by using the external mail command
eventually provided by the underlying system (eg. ``sendmail`` or ``mail``).

If your system does not provide access to none of the above you can configure a SMTP server in the ``nextflow.config`` file.
For example::

    mail {
        smtp.host = 'your.smtp-server.com'
        smtp.port = 475
        smtp.user = 'my-user'
    }

See the :ref:`mail scope <config-mail>` section to learn more the mail server configuration options.


Mail notification
-------------------

You can use the ``sendMail`` function with a :ref:`workflow completion handler <metadata-completion-handler>`
to notify the completion of a workflow completion. For example::

    workflow.onComplete {

        def msg = """\
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            """
            .stripIndent()

        sendMail(to: 'you@gmail.com', subject: 'My pipeline execution', body: msg)
    }


This is useful to send a custom notification message. Note however that Nextflow includes a built-in notification mechanism
which is the most convenient way to notify the completion of a workflow execution in most cases. Read the following
section to learn about it.

Workflow notification
---------------------

Nextflow includes a built-in workflow notification features that automatically sends a notification message
when a workflow execution terminates.

To enable simply specify the ``-N`` option when launching the pipeline execution. For example::

  nextflow run <pipeline name> -N <recipient address>

It will send a notification mail when the execution completes similar to the one shown below:

.. image:: images/workflow-notification-min.png


.. warning:: By default the notification message is sent by using the ``sendmail`` system tool which is assumed to be
    available in the computer where Nextflow is running. Make sure it's properly installed and configured.
    Alternatively provide the SMTP server configuration settings to use the Nextflow
    built-in mail support, which doesn't require any external system tool.

See the `Mail configuration`_ section to learn about the available mail delivery options and configuration settings.

Read :ref:`Notification scope <config-notification>` section to learn more about the workflow notification
configuration details.



.. _charliecloud-page:

************************
Charliecloud containers
************************

`Charliecloud <https://hpc.github.io/charliecloud>`_ is an alternative Container engine to Docker, that is better
suited for use in HPC environments. Its main advantage is that it can be used with unpriviledged permissions,
making use of user namespaces in the Linux kernel. Charliecloud is able to pull from Docker registries.

By using this feature any process in a Nextflow script can be transparently executed into a Charliecloud container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Charliecloud engine.

.. note::
    This feature requires Nextflow version ``21.03.0-edge`` or later and Charliecloud ``v0.22` or later.

.. warning::
    This is an incubating feature. The use in production environment is not recommended.

Prerequisites
==============

You will need Charliecloud version ``0.22`` or later installed on your execution environment e.g. your computer or a
distributed cluster, depending on where you want to run your pipeline.

How it works
=============

You won't need to modify your Nextflow script in order to run it with Charliecloud. Simply specify the docker image from
where the containers are started by using the ``-with-charliecloud`` command line option. For example::

  nextflow run <your script> -with-charliecloud [container]

Every time your script launches a process execution, Nextflow will run it into a charliecloud container created by using the
specified container image. In practice Nextflow will automatically wrap your processes and run them by executing the ``ch-run``
command with the container you have provided.

.. note:: A container image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Container image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/container'
    charliecloud.enabled = true

.. warning::
    If an absolute is provided, the container needs to be in the Charliecloud flat directory format.
    See section below for compatibility with Docker registries.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.

Charliecloud & Docker Hub
=========================

Nextflow is able to transparently pull remote container images stored in any Docker compatible registry and converts
them to the Charliecloud compatible format.

By default when a container name is specified, Nextflow checks if a container with that name exists in the local file
system. If it exists, it's used to execute the container. If a matching file does not exist,
Nextflow automatically tries to pull an image with the specified name from the Docker Hub.

To pull images from a third party Docker registry simply provide the URL to the image. If no URL is provided,
Docker Hub is assumed. For example this can be used to pull an image from quay.io and convert it automatically
to the Charliecloud container format::

    process.container = 'https://quay.io/biocontainers/multiqc:1.3--py35_2'
    charliecloud.enabled = true
 
Whereas this would pull from Docker Hub::

    process.container = 'nextflow/examples:latest'
    charliecloud.enabled = true


Multiple containers
===================

It is possible to specify a different container for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different containers
in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    charliecloud {
        enabled = true
    }

Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
=================

Charliecloud advanced configuration settings are described in :ref:`config-charliecloud` section in the Nextflow
configuration page.
.. _podman-page:

******************
Podman containers
******************

Nextflow integration with `Podman containers <http://www.podman.io>`_ technology allows you to write self-contained
and truly reproducible computational pipelines.

By using this feature any process in a Nextflow script can be transparently executed into a Podman container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Podman engine.

.. note::
    This feature requires Nextflow version ``20.01.0`` or later.

.. warning::
    This is an incubating feature. The use in production environment is not recommended.

Prerequisites
==============

You will need Podman installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline. Running in rootless mode requires appropriate OS configuration. Due to current
Podman limits using cpuset for cpus and memory such is only possible using sudo.


How it works
=============

You won't need to modify your Nextflow script in order to run it with Podman. Simply specify the Podman image from
where the containers are started by using the ``-with-podman`` command line option. For example::

  nextflow run <your script> -with-podman [OCI container image]

Every time your script launches a process execution, Nextflow will run it into a Podman container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``podman run``
command with the image you have provided.

.. note:: A OCI container image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Podman image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    podman.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Podman image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.


Multiple containers
=====================

It is possible to specify a different container image for each process definition in your pipeline script. Let's
suppose you have two processes named ``foo`` and ``bar``. You can specify two different container images for them
in the Nextflow script as shown below::

    process foo {
      container 'image_name_1'

      '''
      do this
      '''
    }

    process bar {
      container 'image_name_2'

      '''
      do that
      '''
    }


Alternatively, the same containers definitions can be provided by using the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    podman {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
==================

Podman advanced configuration settings are described in :ref:`config-podman` section in the Nextflow configuration page.


.. _singularity-page:

**********************
Singularity containers
**********************

`Singularity <http://singularity.lbl.gov/>`_ is a container engine alternative to Docker. The main advantages
of Singularity is that it can be used with unprivileged permissions and doesn't require a separate daemon process.
These, along other features, like for example the support for autofs mounts, makes Singularity a container engine
better suited the requirements of HPC workloads. Singularity is able to use existing Docker images, and pull from Docker
registries.

Nextflow provides built-in support for Singularity. This allows you to precisely control the execution environment
of the processes in your pipeline by running them in isolated containers along all their dependencies.

Moreover the support provided by Nextflow for different container technologies, allows the same pipeline to be
transparently executed both with Docker or Singularity containers, depending on the available engine in the target
execution platforms.


Prerequisites
=============

You will need Singularity installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

Images
======

Singularity makes use of a container image file, which physically contains the container. Refer to the `Singularity
documentation <https://www.sylabs.io/docs/>`_ to learn how create Singularity images.

Singularity allows paths that do not currently exist within the container to be created
and mounted dynamically by specifying them on the command line. However this feature is only supported on hosts
that support the `Overlay file system <https://en.wikipedia.org/wiki/OverlayFS>`_ and is not enabled by default.

.. note::
    Nextflow expects that data paths are defined system wide, and your Singularity images need to be created having the
    mount paths defined in the container file system.

If your Singularity installation support the `user bind control` feature,
enable the Nextflow support for that by defining the ``singularity.autoMounts = true`` setting in the Nextflow
configuration file.


How it works
============

The integration for Singularity follows the same execution model implemented for Docker. You won't need to modify your
Nextflow script in order to run it with Singularity. Simply specify the Singularity image
file from where the containers are started by using the ``-with-singularity`` command line option. For example::

  nextflow run <your script> -with-singularity [singularity image file]

Every time your script launches a process execution, Nextflow will run it into a Singularity container created by using the
specified image. In practice Nextflow will automatically wrap your processes and launch them by running the
``singularity exec`` command with the image you have provided.

.. note:: A Singularity image can contain any tool or piece of software you may need to carry out a process execution.
  Moreover the container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Singularity image as a command line parameter, you can define it in the Nextflow
configuration file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = '/path/to/singularity.img'
    singularity.enabled = true

In the above example replace ``/path/to/singularity.img`` with any Singularity image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. note::
   Unlike Docker, Nextflow does not mount automatically host paths in the container when using Singularity.
   It expects they are configure and mounted system wide by the Singularity runtime. If your Singularity installation
   allows `user defined bind points` read the :ref:`Singularity configuration <config-singularity>` section to learn
   how to enable Nextflow auto mounts.

.. warning::
    When a process input is a *symbolic link* file, make sure the linked file **is** stored in a host folder that is
    accessible from a bind path defined in your Singularity installation. Otherwise the process execution will fail
    because the launched container won't be able to access the linked file.


Multiple containers
===================

It is possible to specify a different Singularity image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Singularity images
specifing them in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    singularity {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Singularity & Docker Hub
========================

Nextflow is able to transparently pull remote container images stored in the `Singularity-Hub <https://singularity-hub.org/>`_,
`Singularity Library <https://cloud.sylabs.io/library/>`_, or any Docker compatible registry.

.. note:: This feature requires you have installed Singularity 2.3.x or higher

By default when a container name is specified, Nextflow checks if an image file with that name exists in the local file
system. If that image file exists, it's used to execute the container. If a matching file does not exist,
Nextflow automatically tries to pull an image with the specified name from the Docker Hub.

If you want Nextflow to check only for local file images, prefix the container name with the ``file://`` pseudo-protocol.
For example::

    process.container = 'file:///path/to/singularity.img'
    singularity.enabled = true

.. warning:: Note the use of triple ``/`` to specify an **absolute** file path, otherwise the path is interpreted as
 relative to the workflow launching directory.

To pull images from the Singularity Hub or a third party Docker registry simply prefix the image name
with the ``shub://``, ``docker://`` or ``docker-daemon://`` pseudo-protocol as required by the Singularity tool. For example::

    process.container = 'docker://quay.io/biocontainers/multiqc:1.3--py35_2'
    singularity.enabled = true
    
.. note:: As of Nextflow v0.27 you no longer need to specify `docker://` to pull from a Docker repository. Nextflow will automatically add it to your image name when Singularity is enabled. Additionally, the Docker engine will not work with containers specified as `docker://`. 

Nextflow version 18.10 introduced support for the `Singularity Library <https://cloud.sylabs.io/library/>`_ repository. This feature also requires Singularity 3.0::
   
   process.container = 'library://library/default/alpine:3.8'

The pseudo-protocol allows you to import Singularity using a local Docker installation instead of downloading
the container image from the Docker registry. It requires Nextflow 19.04.0 or later and Singularity 3.0.3 or later.

.. note:: This feature requires the availability of the ``singularity`` tool in the computer
  where the workflow execution is launched (other than the computing nodes).

Nextflow caches those images in the ``singularity`` directory in the pipeline work directory by default. However it is
suggest to provide a centralised caching directory by using either the ``NXF_SINGULARITY_CACHEDIR`` environment variable
or the ``singularity.cacheDir`` setting in the Nextflow config file.

As of version ``21.09.0-edge``, when looking for a Singularity image file, Nextflow before checks the *library* directory,
if the image file is not found, the *cache* directory is used as usual. The library directory can be defined either using
the ``NXF_SINGULARITY_LIBRARYDIR`` environment variable or the ``singularity.libraryDir`` configuration setting (the
latter overrides the former).

.. warning:: When using a computing cluster the Singularity cache directory must be a shared folder accessible
  to all computing nodes.

.. error::  When pulling Docker images Singularity may be unable to determine the container size if the image was
  stored by using an old Docker format, resulting in a pipeline execution error. See the Singularity documentation for details.

Advanced settings
=================

Singularity advanced configuration settings are described in :ref:`config-singularity` section in the Nextflow
configuration page.













.. _ignite-page:

*************
Apache Ignite
*************

.. danger::
  This feature has been phased out and it's not supported any more as of version 22.01.x.

Nextflow can be be deployed in a *cluster* mode by using `Apache Ignite <https://ignite.apache.org/>`_, an in-memory data-grid
and clustering platform.

Apache Ignite is packaged with Nextflow itself, so you won't need to install it separately or configure other third party
software.

.. _ignite-daemon:

Cluster daemon
--------------

In order to setup a cluster you will need to run a cluster daemon on each node within your cluster.
If you want to support the :ref:`Docker integration <docker-page>` feature provided by Nextflow, the Docker engine has
to be installed and must run in each node.

In its simplest form just launch the Nextflow daemon in each cluster node as shown below::

    nextflow node -bg

The command line option ``-bg`` launches the node daemon in the background. The output is stored in the log file ``.node-nextflow.log``. The daemon
process ``PID`` is saved in the file ``.nextflow.pid`` in the same folder.


Multicast discovery
===================

By default, the Ignite daemon tries to automatically discover all members in the cluster by using the network *multicast* discovery.
Note that NOT all networks support this feature (for example Amazon AWS does not).

.. tip::  To check if multicast is available in your network, use the `iperf <http://sourceforge.net/projects/iperf/>`_ tool.
  To test multicast, open a terminal on two machines within the network and run the following command in the first one
  ``iperf -s -u -B 228.1.2.4 -i 1`` and ``iperf -c 228.1.2.4 -u -T 32 -t 3 -i 1`` on the second.
  If data is being transferred then multicast is working.


Ignite uses the multicast group ``228.1.2.4`` and port ``47400`` by default. You can change these values, by using the
``cluster.join`` command line option, as shown below::

    nextflow node -cluster.join multicast:224.2.2.3:55701



In case multicast discovery is not available in your network, you can try one of the following alternative methods:

Shared file system
==================

Simply provide a path shared across the cluster by a network file system, as shown below::

    nextflow node -bg -cluster.join path:/net/shared/cluster


The cluster members will use that path to discover each other.


IP addresses
============

Provide a list of pre-configured IP addresses on the daemon launch command line, for example::

    nextflow node -cluster.join ip:10.0.2.1,10.0.2.2,10.0.2.4


Advanced options
=====================

The following cluster node configuration options can be used.

============================= ================
Name                          Description
============================= ================
join                          IP address(es) of one or more cluster nodes to which the daemon will join.
group                         The name of the cluster to which this node join. It allows the creation of separate cluster instances. Default: ``nextflow``
maxCpus                       Max number of CPUs that can be used by the daemon to run user tasks. By default it is equal to the number of CPU cores.
maxDisk                       Max amount of disk storage that can be used by user tasks eg. ``500 GB``.
maxMemory                     Max amount of memory that can be used by user tasks eg. ``64 GB``. Default total available memory.
interface                     Network interfaces that Ignite will use. It can be the interface IP address or name.
failureDetectionTimeout       Failure detection timeout is used to determine how long the communication or discovery SPIs should wait before considering a remote connection failed.
clientFailureDetectionTimeout Failure detection timeout is used to determine how long the communication or discovery SPIs should wait before considering a remote connection failed.
tcp.localAddress              Defines the local host IP address.
tcp.localPort                 Defines the local port to listen to.
tcp.localPortRange            Range for local ports.
tcp.reconnectCount            Number of times the node tries to (re)establish connection to another node.
tcp.networkTimeout            Defines the network timeout.
tcp.socketTimeout             Defines the socket operations timeout. This timeout is used to limit connection time and write-to-socket time. Note that when running Ignite on Amazon EC2, socket timeout must be set to a value significantly greater than the default (e.g. to 30000).
tcp.ackTimeout                Defines the timeout for receiving acknowledgement for sent message.
tcp.maxAckTimeout             Defines the maximum timeout for receiving acknowledgement for sent message.
tcp.joinTimeout               Defines the join timeout.
============================= ================



These options can be specified as command line parameters by adding the prefix ``-cluster.`` to them, as shown below::

    nextflow node -bg -cluster.maxCpus 4 -cluster.interface eth0

The same options can be entered into the Nextflow :ref:`configuration file<config-page>`, as shown below::

    cluster {
        join = 'ip:192.168.1.104'
        interface = 'eth0'
    }

Finally daemon options can be provided also as environment variables having the name in upper-case and by adding
the prefix ``NXF_CLUSTER_`` to them, for example::

    export NXF_CLUSTER_JOIN='ip:192.168.1.104'
    export NXF_CLUSTER_INTERFACE='eth0'


Pipeline execution
------------------

The pipeline execution needs to be launched in a `head` node i.e. a cluster node where the Nextflow node daemon
is **not** running. In order to execute your pipeline in the Ignite cluster you will need to use the Ignite executor,
as shown below::

   nextflow run <your pipeline> -process.executor ignite


If your network does no support multicast discovery, you will need to specify the `joining` strategy as you did for the
cluster daemons. For example, using a shared path::

    nextflow run <your pipeline> -process.executor ignite -cluster.join path:/net/shared/path



Execution with MPI
------------------

Nextflow is able to deploy and self-configure an Ignite cluster on demand, taking advantage of the Open `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_
standard that is commonly available in grid and supercomputer facilities.

In this scenario a Nextflow workflow needs to be executed as an MPI job. Under the hood Nextflow will launch a `driver`
process in the first of the nodes, allocated by your job request, and an Ignite daemon in the remaining nodes.

In practice you will need a launcher script to submit an MPI job request to your batch scheduler/resource manager.
The batch scheduler must reserve the computing nodes in an exclusive manner to avoid having multiple Ignite daemons
running on the same node. Nextflow must be launched using the ``mpirun`` utility, as if it were an MPI application,
specifying the ``--pernode`` option.

Platform LSF launcher
=====================

The following example shows a launcher script for the `Platform LSF <https://en.wikipedia.org/wiki/Platform_LSF/>`_ resource manager::

    #!/bin/bash
    #BSUB -oo output_%J.out
    #BSUB -eo output_%J.err
    #BSUB -J <job name>
    #BSUB -q <queue name>
    #BSUB -W 02:00
    #BSUB -x
    #BSUB -n 80
    #BSUB -M 10240
    #BSUB -R "span[ptile=16]"
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

It requests 5 nodes (80 processes, with 16 cpus per node). The ``-x`` directive allocates the node in an exclusive manner.
Nextflow needs to be executed using the ``-with-mpi`` command line option. It will automatically use ``ignite`` as the executor.

The variable ``NXF_CLUSTER_SEED`` must contain an integer value (in the range 0-16777216) that will unequivocally identify
your cluster instance. In the above example it is randomly generated by using the `shuf <http://linux.die.net/man/1/shuf>`_ Linux tool.

Univa Grid Engine launcher
==========================

The following example shows a launcher script for the `Univa Grid Engine <https://en.wikipedia.org/wiki/Univa_Grid_Engine>`_ (aka SGE)::

    #!/bin/bash
    #$ -cwd
    #$ -j y
    #$ -o <output file name>
    #$ -l virtual_free=10G
    #$ -q <queue name>
    #$ -N <job name>
    #$ -pe ompi 5
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode nextflow run <your-project-name> -with-mpi [pipeline parameters]

As in the previous script it allocates 5 processing nodes. UGE/SGE does not have an option to reserve a node in an exclusive
manner. A common workaround is to request the maximum amount of memory or cpus available in the nodes of your cluster.


Linux SLURM launcher
====================

When using Linux SLURM you will need to use ``srun`` instead ``mpirun`` in your launcher script. For example::

    #!/bin/bash
    #SBATCH --job-name=<job name>
    #SBATCH --output=<log file %j>
    #SBATCH --ntasks=5
    #SBATCH --cpus-per-task=16
    #SBATCH --tasks-per-node=1
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    srun nextflow run hello.nf -with-mpi

As before, this allocates 5 processing nodes (``--ntasks=5``) and each node will be able to use up to 16 cpus
(``--cpus-per-task=16``). When using SLURM it's not necessary to allocate computing nodes in an exclusive manner.
It's even possible to launch more than one Nextflow daemon instance per node, though not suggested.

To submit the pipeline execution create a file like the above, then use the following command::

    sbatch <launcher script name>

.. _azure-page:

************
Azure Cloud
************

Requirements
============

The support for Azure Cloud requires Nextflow version ``21.04.0`` or later. If you don't have it installed
use the following command to download it in your computer::

    curl get.nextflow.io | bash
    ./nextflow self-update


.. _azure-blobstorage:

Azure Blob Storage
===================

Nextflow has built-in support for `Azure Blob Storage <https://azure.microsoft.com/en-us/services/storage/blobs/>`_.
Files stored in a Azure blob container can be accessed transparently in your pipeline script like any other file
in the local file system.

The Blob storage account `name` and `key` needs to be provided in the Nextflow configuration file as shown below::

    azure {
      storage {
        accountName = "<YOUR BLOB ACCOUNT NAME>"
        accountKey = "<YOUR BLOB ACCOUNT KEY>"
      }
    }

As an alternative to the account key it can also used `Shared Access Token` using the setting ``sasToken`` in place
of the ``accountKey`` attribute.

.. tip::
    When creating the `Shared Access Token` make sure to allow the resource types `Container` and `Object` and allow
    the permissions: `Read`, `Write`, `Delete`, `List`, `Add`, `Create`.

.. tip::
    The value of ``sasToken`` is the token stripped by the character ``?`` from the beginning of the token.


Once the Blob Storage credentials are set you can access the files in the blob container like local files prepending
the file path with the ``az://`` prefix followed by the container name. For example, having a container named ``my-data``
container a file named ``foo.txt`` you can access it in your Nextflow script using the following fully qualified
file path ``az://my-data/foo.txt``.

Azure File Shares
==================

As of version nf-azure@0.11.0, Nextflow has built-in support also for `Azure Files <https://azure.microsoft.com/en-us/services/storage/files/>`_.
Files available in the serverless Azure File shares can be mounted concurrently on the nodes of a pool executing the pipeline.
These files become immediately available in the file system and can be referred as local files within the processes. This
is especially useful when a task needs to access large amount of data (such as genome indexes) during its execution. An
arbitrary number of File shares can be mounted on each pool noode.

The Azure File share must exist in the storage account configured for Blob Storage.
The name of the source Azure File share and mount path (the destination path where the files are mounted) must be provided.
Additional mount options (see the Azure Files documentation) can be set as well for further customisation of the mounting process.

For example::

	azure {
		storage {
			accountName = "<YOUR BLOB ACCOUNT NAME>"
			accountKey = "<YOUR BLOB ACCOUNT KEY>"
			fileShares {
				<YOUR SOURCE FILE SHARE NAME> {
					mountPath = "<YOUR MOUNT DESTINATION>"
					mountOptions = "<SAME AS MOUNT COMMAND>" //optional
				}
				<YOUR SOURCE FILE SHARE NAME> {
					mountPath = "<YOUR MOUNT DESTINATION>"
					mountOptions = "<SAME AS MOUNT COMMAND>" //optional
				}
			}
		}
	}

The files in the File share are available to the task in the directory:
`<YOUR MOUNT DESTINATION>/<YOUR SOURCE FILE SHARE NAME>`.

For instance, given the following configuration::

	azure {
		storage {
			...
			fileShares {
				dir1 {
					mountPath = "/mnt/mydata/"
				}
			}
		}
	}


The task can access to the File share in `/mnt/mydata/dir1`.

.. _azure-batch:

Azure Batch
============

`Azure Batch <https://docs.microsoft.com/en-us/azure/batch/>`_ is a managed computing service that allows the execution
of containerised workloads in the Azure cloud infrastructure.

Nextflow provides a built-in support for Azure Batch which allows the seamless deployment of a Nextflow pipeline in the cloud
offloading the process executions as Batch jobs.

Get started
-------------

1 - Create a Batch account in Azure portal. Take note of the account name and key.

2 - Make sure to adjust your quotas on the pipeline's needs. There are limits on certain resources associated with the
Batch account. Many of these limits are default quotas applied by Azure at the subscription or account level.
Quotas impact on the number of Pools, CPUs and Jobs you can create at any given time.

3 - Create a Storage account and, within, an Azure Blob Container in the same location where the Batch account was created.
Take note of the account name and key.

4 - If planning to use Azure files, create an Azure File share within the same Storage account and upload there
the data to mount on the pool nodes.

5 - Associate the Storage account with the Azure Batch account.

6 - Make sure your pipeline processes specify one or more Docker containers by using the :ref:`process-container` directive.

7 - The container images need to be published into Docker registry such as `Docker Hub <https://hub.docker.com/>`_,
`Quay <https://quay.io/>`_ or `Azure Container Registry <https://docs.microsoft.com/en-us/azure/container-registry/>`_
that can be reached by Azure Batch environment.


A minimal configuration looks like the following snippet::

    process {
      executor = 'azurebatch'
    }

    azure {
      storage {
        accountName = "<YOUR STORAGE ACCOUNT NAME>"
        accountKey = "<YOUR STORAGE ACCOUNT KEY>"
      }
      batch {
        location = '<YOUR LOCATION>'
        accountName = '<YOUR BATCH ACCOUNT NAME>'
        accountKey = '<YOUR BATCH ACCOUNT KEY>'
        autoPoolMode = true
      }
    }

In the above example, replace the location and the account placeholders with the value corresponding to your configuration and
save it to a file named ``nextflow.config``.

Given the previous configuration, launch the execution of the pipeline using the following command::

    nextflow run <PIPELINE NAME> -w az://YOUR-CONTAINER/work


Replacing ``<PIPELINE NAME>`` with a pipeline name e.g. ``nextflow-io/rnaseq-nf`` and ``YOUR-CONTAINER`` a blob
container in the storage account defined in the above configuration.

See the `Batch documentation <https://docs.microsoft.com/en-us/azure/batch/quick-create-portal>`_ for further
details about the configuration for the Azure Batch service.


Pools configuration
-------------------

When using the ``autoPoolMode`` setting Nextflow automatically creates a `pool` of computing nodes to execute the
jobs run by your pipeline. By default it only uses 1 compute node of ``Standard_D4_v3`` type.

The pool is not removed when the pipeline execution terminates, unless the configuration setting ``deletePoolsOnCompletion=true``
is added in your pipeline configuration file.

Pool specific settings, e.g. VM type and count, should be provided in the ``auto`` pool configuration scope, e.g. ::

    azure {
        batch {
            pools {
                auto {
                   vmType = 'Standard_D2_v2'
                   vmCount = 10
                }
            }
        }
    }



.. warning::
    Don't forget to clean up the Batch pools to avoid in extra charges in the Batch account or use the auto scaling feature.

.. warning::
   Make sure your Batch account has enough resources to satisfy the pipeline's requirements and the pool configuration.

.. warning::
   Nextflow uses the same pool Id across pipeline executions, if the pool features have not changed.
   Therefore, when using ``deletePoolsOnCompletion=true``, make sure the pool is completely removed from the Azure Batch account
   before re-running the pipeline. The following message is returned when the pool is still shutting down ::


    Error executing process > '<process name> (1)'
    Caused by:
        Azure Batch pool '<pool name>' not in active state

Named pools
-------------

If you want to have a more precise control on the computing nodes pools used in your pipeline using a different pool
depending on the task in your pipeline, you can use the Nextflow :ref:`process-queue` directive to specify the *ID* of a
Azure Batch compute pool that has to be used to run that process' tasks.

The pool is expected to be already available in the Batch environment, unless the setting ``allowPoolCreation=true`` is
provided in the ``batch`` setting in the pipeline configuration file. In the latter case Nextflow will create the pools on-demand.

The configuration details for each pool can be specified using a snippet as shown below::

    azure {
        batch {
            pools {
                foo {
                   vmType = 'Standard_D2_v2'
                   vmCount = 10
                }

                bar {
                    vmType = 'Standard_E2_v3'
                    vmCount = 5
                }
            }
        }
    }

The above example defines the configuration for two node pools. The first will provision 10 compute nodes of type ``Standard_D2_v2``,
the second 5 nodes of type ``Standard_E2_v3``. See the `Advanced settings`_ below for the complete list of available
configuration options.

Requirements on pre-existing named pools
----------------------------------------

When Nextflow is configured to use a pool already available in the Batch account, the target pool must satisfy the following
requirements:

1 - the pool must be declared as ``dockerCompatible`` (``Container Type`` property)

2 - the task slots per node must match with the number of cores for the selected VM. Nextflow would return an error like
"Azure Batch pool 'ID' slots per node does not match the VM num cores (slots: N, cores: Y)".

Pool autoscaling
----------------

Azure Batch can automatically scale pools based on parameters that you define, saving you time and money. With automatic scaling,
Batch dynamically adds nodes to a pool as task demands increase, and removes compute nodes as task demands decrease.

To enable this feature for pools created by Nextflow, add the option ``autoScale = true`` to the corresponding pool configuration scope.
For example, when using the ``autoPoolMode``, the setting looks like::

    azure {
        batch {
            pools {
                auto {
                   autoScale = true
                   vmType = 'Standard_D2_v2'
                   vmCount = 5
                   maxVmCount = 50
                }
            }
        }
    }

Nextflow uses the formula shown below to determine the number of VMs to be provisioned in the pool::

        // Get pool lifetime since creation.
        lifespan = time() - time("{{poolCreationTime}}");
        interval = TimeInterval_Minute * {{scaleInterval}};

        // Compute the target nodes based on pending tasks.
        // $PendingTasks == The sum of $ActiveTasks and $RunningTasks
        $samples = $PendingTasks.GetSamplePercent(interval);
        $tasks = $samples < 70 ? max(0, $PendingTasks.GetSample(1)) : max( $PendingTasks.GetSample(1), avg($PendingTasks.GetSample(interval)));
        $targetVMs = $tasks > 0 ? $tasks : max(0, $TargetDedicatedNodes/2);
        targetPoolSize = max(0, min($targetVMs, {{maxVmCount}}));

        // For first interval deploy 1 node, for other intervals scale up/down as per tasks.
        $TargetDedicatedNodes = lifespan < interval ? {{vmCount}} : targetPoolSize;
        $NodeDeallocationOption = taskcompletion;


The above formula initialises a pool with the number of VMs specified by the ``vmCount`` option, it scales up the pool on-demand,
based on the number of pending tasks up to ``maxVmCount`` nodes. If no jobs are submitted for execution, it scales down
to zero nodes automatically.

If you need a different strategy you can provide your own formula using the ``scaleFormula`` option.
See the `Azure Batch <https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling>`_ documentation for details.

Pool nodes
-----------
When Nextflow creates a pool of compute nodes, it selects:

* the virtual machine image reference to be installed on the node
* the Batch node agent SKU, a program that runs on each node and provides an interface between the node and the Batch service

Together, these settings determine the Operating System and version installed on each node.

By default, Nextflow creates CentOS 8-based pool nodes, but this behavior can be customised in the pool configuration.
Below the configurations for image reference/SKU combinations to select two popular systems.

* Ubuntu 20.04::

	sku = "batch.node.ubuntu 20.04"
	offer = "ubuntu-server-container"
	publisher = "microsoft-azure-batch"

* CentOS 8 (default)::

	sku = "batch.node.centos 8"
	offer = "centos-container"
	publisher = "microsoft-azure-batch"

See the `Advanced settings`_ below and `Azure Batch nodes <https://docs.microsoft.com/en-us/azure/batch/batch-linux-nodes>` documentation for more details.

Private container registry
--------------------------
As of version ``21.05.0-edge``, a private container registry from where to pull Docker images can be optionally specified as follows ::

    azure {
        registry {
            server =  '<YOUR REGISTRY SERVER>' // e.g.: docker.io, quay.io, <ACCOUNT>.azurecr.io, etc.
            userName =  '<YOUR REGISTRY USER NAME>'
            password =  '<YOUR REGISTRY PASSWORD>'
        }
    }


The private registry is not exclusive, rather it is an addition to the configuration.
Public images from other registries are still pulled (if requested by a Task) when a private registry is configured.

.. note::
  When using containers hosted into a private registry, the registry name must also be provided in the container name
  specified via the :ref:`container <process-container>` directive using the format: ``[server]/[your-organization]/[your-image]:[tag]``.
  Read more about image fully qualified image names in the `Docker documentation <https://docs.docker.com/engine/reference/commandline/pull/#pull-from-a-different-registry>`_.

Advanced settings
==================

The following configuration options are available:

============================================== =================
Name                                           Description
============================================== =================
azure.storage.accountName                       The blob storage account name
azure.storage.accountKey                        The blob storage account key
azure.storage.sasToken                          The blob storage shared access signature token. This can be provided as an alternative to the ``accountKey`` setting.
azure.storage.tokenDuration                     The duration of the shared access signature token created by Nextflow when the ``sasToken`` option is *not* specified (default: ``12h``).
azure.batch.accountName                         The batch service account name.
azure.batch.accountKey                          The batch service account key.
azure.batch.endpoint                            The batch service endpoint e.g. ``https://nfbatch1.westeurope.batch.azure.com``.
azure.batch.location                            The batch service location e.g. ``westeurope``. This is not needed when the endpoint is specified.
azure.batch.autoPoolMode                        Enable the automatic creation of batch pools depending on the pipeline resources demand (default: ``true``).
azure.batch.allowPoolCreation                   Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: ``false``).
azure.batch.deleteJobsOnCompletion              Enable the automatic deletion of jobs created by the pipeline execution (default: ``true``).
azure.batch.deletePoolsOnCompletion             Enable the automatic deletion of compute node pools upon pipeline completion (default: ``false``).
azure.batch.copyToolInstallMode                 Specify where the `azcopy` tool used by Nextflow. When ``node`` is specified it's copied once during the pool creation. When ``task`` is provider, it's installed for each task execution (default: ``node``).
azure.batch.pools.<name>.publisher              Specify the publisher of virtual machine type used by the pool identified with ``<name>`` (default: ``microsoft-azure-batch``, requires ``nf-azure@0.11.0``).
azure.batch.pools.<name>.offer                  Specify the offer type of the virtual machine type used by the pool identified with ``<name>`` (default: ``centos-container``, requires ``nf-azure@0.11.0``).
azure.batch.pools.<name>.sku                    Specify the ID of the Compute Node agent SKU which the pool identified with ``<name>`` supports (default: ``batch.node.centos 8``, requires ``nf-azure@0.11.0``).
azure.batch.pools.<name>.vmType                 Specify the virtual machine type used by the pool identified with ``<name>``.
azure.batch.pools.<name>.vmCount                Specify the number of virtual machines provisioned by the pool identified with ``<name>``.
azure.batch.pools.<name>.maxVmCount             Specify the max of virtual machine when using auto scale option.
azure.batch.pools.<name>.autoScale              Enable autoscaling feature for the pool identified with ``<name>``.
azure.batch.pools.<name>.fileShareRootPath      If mounting File Shares, this is the internal root mounting point. Must be ``/mnt/resource/batch/tasks/fsmounts`` for CentOS nodes or ``/mnt/batch/tasks/fsmounts`` for Ubuntu nodes (default is for CentOS, requires ``nf-azure@0.11.0``).
azure.batch.pools.<name>.scaleFormula           Specify the scale formula for the pool identified with ``<name>``. See Azure Batch `scaling documentation <https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling>`_ for details.
azure.batch.pools.<name>.scaleInterval          Specify the interval at which to automatically adjust the Pool size according to the autoscale formula. The minimum and maximum value are 5 minutes and 168 hours respectively (default: `10 mins`).
azure.batch.pools.<name>.schedulePolicy         Specify the scheduling policy for the pool identified with ``<name>``. It can be either ``spread`` or ``pack`` (default: ``spread``).
azure.batch.pools.<name>.privileged             Enable the task to run with elevated access. Ignored if `runAs` is set (default: ``false``).
azure.batch.pools.<name>.runAs                  Specify the username under which the task is run. The user must already exist on each node of the pool.
azure.registry.server                           Specify the container registry from which to pull the Docker images (default: ``docker.io``, requires ``nf-azure@0.9.8``).
azure.registry.userName                         Specify the username to connect to a private container registry (requires ``nf-azure@0.9.8``).
azure.registry.password                         Specify the password to connect to a private container registry (requires ``nf-azure@0.9.8``).
============================================== =================
.. _dnanexus-page:

****************
DNAnexus cloud
****************

Extend the Nextflow runtime for scalability and performance on the DNAnexus platform.
Check `Seqera Labs <https://seqera.io/dnanexus>`_ website for details.

.. raw:: html

    <script>
    window.location.replace("https://seqera.io/dnanexus");
    </script>


.. _channel-page:

********
Channels
********

Nextflow is based on the Dataflow programming model in which processes communicate through channels.

A channel has two major properties:

#. Sending a message is an `asynchronous` operation which completes immediately,
   without having to wait for the receiving process.

#. Receiving data is a blocking operation which stops the receiving process until the message has arrived.

.. _channel-types:

Channel types
=============

Nextflow distinguish two different kinds of channels: `queue channels` and `value channels`.

.. _channel-type-queue:

Queue channel
-------------

A `queue channel` is a non-blocking unidirectional FIFO queue which connects two processes or operators.

A queue channel is usually created using a factory method such as a `from`_, `fromPath`_, etc.
or chaining it with a channel operator such as :ref:`operator-map`, :ref:`operator-flatmap`, etc.

Queue channels are also created by process output declarations using the ``into`` clause.

.. note:: The definition implies that the same queue channel cannot be used more than one time as process
 output and more than one time as process input.

If you need to connect a process output channel to more than one process or operator use the
:ref:`operator-into` operator to create two (or more) copies of the same channel and use each
of them to connect a separate process.

.. _channel-type-value:

Value channel
-------------

A `value channel` a.k.a. *singleton channel* by definition is bound to a single value and it can be read
unlimited times without consuming its content.

.. tip:: For this reason a value channel can be used as input by more than one process.

A value channel is created using the `value`_ factory method or by operators returning
a single value, such us :ref:`operator-first`, :ref:`operator-last`, :ref:`operator-collect`,
:ref:`operator-count`, :ref:`operator-min`, :ref:`operator-max`, :ref:`operator-reduce`, :ref:`operator-sum`.


.. note:: A value channel is implicitly created by a process when an input specifies a simple value
  in the ``from`` clause.
  Moreover, a value channel is also implicitly created as output for a process whose
  inputs are only value channels.

For example::

    process foo {
      input:
      val x from 1
      output:
      file 'x.txt' into result

      """
      echo $x > x.txt
      """
    }

The process in the above snippet declares a single input which implicitly is a value channel.
Therefore also the ``result`` output is a value channel that can be read by more than one process.

See also: :ref:`process-understand-how-multiple-input-channels-work`.

.. _channel-factory:

Channel factory
===============

Channels may be created implicitly by the process output(s) declaration or explicitly using the following channel
factory methods.

The available factory methods are:

* `create`_
* `empty`_
* `from`_
* `fromPath`_
* `fromFilePairs`_
* `fromSRA`_
* `of`_
* `value`_
* `watchPath`_

.. tip::
  As of version 20.07.0 the prefix ``channel.`` has been introduced as an alias of ``Channel.``, therefore factory
  methods can be used either with the syntaxes ``channel.from()`` and ``Channel.from()``, and so on.

.. _channel-create:

create
------

.. warning::
    This method is deprecated and won't be available in DSL2 syntax.

Creates a new `channel` by using the ``create`` method, as shown below::

    channelObj = Channel.create()


.. _channel-of:

of
--

The ``of`` method allows you to create a channel emitting any sequence of values that are specified as the method argument,
for example::

    ch = Channel.of( 1, 3, 5, 7 )
    ch.view { "value: $it" }

The first line in this example creates a variable ``ch`` which holds a channel object. This channel emits the values
specified as a parameter in the ``of`` method. Thus the second line prints the following::

    value: 1
    value: 3
    value: 5
    value: 7


.. tip::
    Range of values are expanded accordingly.

::

    Channel
        .of(1..23, 'X', 'Y')
        .view()

Prints::

    1
    2
    3
    4
    :
    23
    X
    Y

.. note::
  This feature requires Nextflow version 19.10.0 of later.

See also: `fromList`_ factory method.

.. _channel-from:

from
----

.. warning::
  This method is deprecated and should only be used for backward compatibility in legacy code.
  Use `of`_ or `fromList`_ instead.

The ``from`` method allows you to create a channel emitting any sequence of values that are specified as the method argument,
for example::

    ch = Channel.from( 1, 3, 5, 7 )
    ch.subscribe { println "value: $it" }

The first line in this example creates a variable ``ch`` which holds a channel object. This channel emits the values
specified as a parameter in the ``from`` method. Thus the second line will print the following::

    value: 1
    value: 3
    value: 5
    value: 7


The following example shows how to create a channel from a `range` of numbers or strings::

    zeroToNine = Channel.from( 0..9 )
    strings = Channel.from( 'A'..'Z' )



.. note:: Note that when the ``from`` argument is an object implementing the (Java)
  `Collection <http://docs.oracle.com/javase/7/docs/api/java/util/Collection.html>`_ interface, the resulting channel
  emits the collection entries as individual emissions.

Thus the following two declarations produce an identical result even tough in the first case the items are specified
as multiple arguments while in the second case as a single list object argument::

    Channel.from( 1, 3, 5, 7, 9 )
    Channel.from( [1, 3, 5, 7, 9] )


But when more than one argument is provided, they are always managed as `single` emissions. Thus, the following example
creates a channel emitting three entries each of which is a list containing two elements::

    Channel.from( [1, 2], [5,6], [7,9] )



.. _channel-value:

value
-----

The `value` factory method is used to create a *value* channel. An optional not ``null`` argument
can be specified to bind the channel to a specific value. For example::


    expl1 = Channel.value()
    expl2 = Channel.value( 'Hello there' )
    expl3 = Channel.value( [1,2,3,4,5] )


The first line in the example creates an 'empty' variable. The second line creates a channel and binds a string to it.
Finally the last one creates a channel and binds a list object to it that will be emitted as a sole emission.


.. _channel-fromlist:

fromList
--------

The ``fromList`` method allows you to create a channel emitting the values provided as a list of elements,
for example::

    Channel
        .fromList( ['a', 'b', 'c', 'd'] )
        .view { "value: $it" }

Prints::

    value: a
    value: b
    value: c
    value: d


See also: `of`_ factory method.

.. note::
  This feature requires Nextflow version 19.10.0 of later.

.. _channel-path:

fromPath
--------

You can create a channel emitting one or more file paths by using the ``fromPath`` method and specifying a path string
as an argument. For example::

    myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )

The above line creates a channel and binds to it a `Path <http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html>`_
item referring the specified file.

.. note:: It does not check the file existence.

Whenever the ``fromPath`` argument contains a ``*`` or ``?`` wildcard character it is interpreted as a `glob`_ path matcher.
For example::

    myFileChannel = Channel.fromPath( '/data/big/*.txt' )


This example creates a channel and emits as many ``Path`` items as there are files with ``txt`` extension in the ``/data/big`` folder.

.. tip:: Two asterisks, i.e. ``**``, works like ``*`` but crosses directory boundaries.
  This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.

For example::

    files = Channel.fromPath( 'data/**.fa' )
    moreFiles = Channel.fromPath( 'data/**/*.fa' )
    pairFiles = Channel.fromPath( 'data/file_{1,2}.fq' )

The first line returns a channel emitting the files ending with the suffix ``.fa`` in the ``data`` folder `and` recursively
in all its sub-folders. While the second one only emits the files which have the same suffix in `any` sub-folder in the ``data`` path.
Finally the last example emits two files: ``data/file_1.fq`` and ``data/file_2.fq``.

.. note:: As in Linux Bash the ``*`` wildcard does not match against hidden files (i.e. files whose name start with a ``.`` character).

In order to include hidden files, you need to start your pattern with a period character or specify the ``hidden: true`` option. For example::

    expl1 = Channel.fromPath( '/path/.*' )
    expl2 = Channel.fromPath( '/path/.*.fa' )
    expl3 = Channel.fromPath( '/path/*', hidden: true )


The first example returns all hidden files in the specified path. The second one returns all hidden files
ending with the ``.fa`` suffix. Finally the last example returns all files (hidden and non-hidden) in that path.

By default a `glob`_ pattern only looks for `regular file` paths that match the specified criteria, i.e.
it won't return directory paths.

You may use the parameter ``type`` specifying the value ``file``, ``dir`` or ``any`` in order to define what kind of paths
you want. For example::

        myFileChannel = Channel.fromPath( '/path/*b', type: 'dir' )
        myFileChannel = Channel.fromPath( '/path/a*', type: 'any' )

The first example will return all `directory` paths ending with the ``b`` suffix, while the second will return any file
and directory starting with a ``a`` prefix.


=============== ===================
Name            Description
=============== ===================
glob            When ``true`` interprets characters ``*``, ``?``, ``[]`` and ``{}`` as glob wildcards, otherwise handles them as normal characters (default: ``true``)
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: ``true``)
relative        When ``true`` returned paths are relative to the top-most common directory (default: ``false``)
checkIfExists   When ``true`` throws an exception of the specified path do not exist in the file system (default: ``false``)
=============== ===================

.. note:: More than one path or glob pattern can be specified using a list as argument::

      Channel.fromPath( ['/some/path/*.fq', '/other/path/*.fastq'] )

  (requires version 0.31.x or later)

.. _channel-filepairs:

fromFilePairs
-------------

The ``fromFilePairs`` method creates a channel emitting the file pairs matching a `glob`_ pattern provided by the user.
The matching files are emitted as tuples in which the first element is the grouping key of the matching
pair and the second element is the list of files (sorted in lexicographical order). For example::

    Channel
        .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
        .view()

It will produce an output similar to the following::

    [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
    [SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
    [SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
    [SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
    [SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
    [SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]


.. note::
    The glob pattern must contain at least a star wildcard character.

Alternatively it is possible to implement a custom file pair grouping strategy providing a closure which,
given the current file as parameter, returns the grouping key.
For example::

    Channel
        .fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
        .view { ext, files -> "Files with the extension $ext are $files" }


Table of optional parameters available:

=============== ===================
Name            Description
=============== ===================
type            Type of paths returned, either ``file``, ``dir`` or ``any`` (default: ``file``)
hidden          When ``true`` includes hidden files in the resulting paths (default: ``false``)
maxDepth        Maximum number of directory levels to visit (default: `no limit`)
followLinks     When ``true`` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: ``true``)
size            Defines the number of files each emitted item is expected to hold (default: 2). Set to ``-1`` for any.
flat            When ``true`` the matching files are produced as sole elements in the emitted tuples (default: ``false``).
checkIfExists   When ``true`` throws an exception of the specified path do not exist in the file system (default: ``false``)
=============== ===================

.. note:: More than one glob pattern can be specified using a list as argument::

      Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )

  (requires version 0.31.x or later)


.. _channel-fromsra:

fromSRA
-------

The ``fromSRA`` method queries the `NCBI SRA <https://www.ncbi.nlm.nih.gov/sra>`_ database and returns a channel emitting
the FASTQ files matching the specified criteria i.e project or accession number(s). For example::

    Channel
        .fromSRA('SRP043510')
        .view()


It returns::

    [SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
    [SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
    [SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
    [SRR1448793, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/003/SRR1448793/SRR1448793.fastq.gz]
    [SRR1910483, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/003/SRR1910483/SRR1910483.fastq.gz]
    [SRR1910482, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/002/SRR1910482/SRR1910482.fastq.gz]
    (remaining omitted)

Multiple accession IDs can be specified using a list object::

    ids = ['ERR908507', 'ERR908506', 'ERR908505']
    Channel
        .fromSRA(ids)
        .view()

::

    [ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
    [ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
    [ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]


.. note:: Read pairs are implicitly managed are returned as a list of files.

.. tip:: Behind the scene it's uses the NCBI `ESearch <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch>`_
  API, therefore the ``fromSRA`` method allows the usage of any query term supported by this API.

Table of optional parameters available:

=============== ===================
Name            Description
=============== ===================
apiKey          NCBI user API key.
cache           Enable/disable the caching API requests (default: ``true``).
max             Maximum number of entries that can be retried (default: unlimited) .
protocol        Allow choosing the protocol for the resulting remote URLs. Available choices: ``ftp``, ``http``, ``https`` (default: ``ftp``).
=============== ===================


To access the NCBI search service the `NCBI API keys <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities>`_
should be provided either:

* Using the ``apiKey`` optional parameter e.g. ``Channel.fromSRA(ids, apiKey:'0123456789abcdef')``.
* Exporting the ``NCBI_API_KEY`` variable in your environment e.g. ``export NCBI_API_KEY=0123456789abcdef``.

.. note:: This feature requires Nextflow version 19.04.0 or later.

.. _channel-watch:

watchPath
---------

The ``watchPath`` method watches a folder for one or more files matching a specified pattern. As soon as
there is a file that meets the specified condition, it is emitted over the channel that is returned by the ``watchPath``
method. The condition on files to watch can be specified by using ``*`` or ``?`` wildcard characters i.e. by specifying
a `glob`_ path matching criteria.

For example::

     Channel
        .watchPath( '/path/*.fa' )
        .subscribe { println "Fasta file: $it" }


By default it watches only for new files created in the specified folder. Optionally, it is possible to provide a
second argument that specifies what event(s) to watch. The supported events are:

=========== ================
Name        Description
=========== ================
``create``  A new file is created (default)
``modify``  A file is modified
``delete``  A file is deleted
=========== ================

You can specified more than one of these events by using a comma separated string as shown below::

     Channel
        .watchPath( '/path/*.fa', 'create,modify' )
        .subscribe { println "File created or modified: $it" }


.. warning:: The ``watchPath`` factory waits endlessly for files that match the specified pattern and event(s).
  Thus, whenever you use it in your script, the resulting pipeline will never finish.

See also: `fromPath`_ factory method.


.. _channel-empty:

empty
-----

The ``empty`` factory method, by definition, creates a channel that doesn't emit any value.

See also: :ref:`operator-ifempty` and :ref:`operator-close` operators.


Binding values
==============

Since in `Nextflow` channels are implemented using `dataflow` variables or queues. Thus sending a message
is equivalent to `bind` a value to object representing the communication channel.


.. _channel-bind1:

bind
----

Channel objects provide a `bind( )` method which is the basic operation to send a message over the channel.
For example::

    myChannel = Channel.create()
    myChannel.bind( 'Hello world' )


.. _channel-bind2:

operator <<
-----------

The operator ``<<`` is just a syntax sugar for the `bind` method. Thus, the following example produce
an identical result as the previous one::

    myChannel = Channel.create()
    myChannel << 'Hello world'



Observing events
================


.. _channel-subscribe:

subscribe
---------

The ``subscribe`` method allows you to execute a user defined function each time a new value is emitted by the source channel.

The emitted value is passed implicitly to the specified function. For example::

    // define a channel emitting three values
    source = Channel.from ( 'alpha', 'beta', 'delta' )

    // subscribe a function to the channel printing the emitted values
    source.subscribe {  println "Got: $it"  }

::

    Got: alpha
    Got: beta
    Got: delta


.. note:: Formally the user defined function is a ``Closure`` as defined by the Groovy programming language on which
  the `Nextflow` scripts are based.

If needed the closure parameter can be defined explicitly, using a name other than ``it`` and, optionally,
specifying the expected value type, as shown in the following example::

    Channel
        .from( 'alpha', 'beta', 'lambda' )
        .subscribe { String str ->
            println "Got: ${str}; len: ${str.size()}"
         }

::

    Got: alpha; len: 5
    Got: beta; len: 4
    Got: lambda; len: 6

Read :ref:`script-closure` paragraph to learn more about `closure` feature.


onNext, onComplete, and onError
-------------------------------

The ``subscribe`` method may accept one or more of the following event handlers:

* ``onNext``: registers a function that is invoked whenever the channel emits a value.
  This is the same as using the ``subscribe`` with a `plain` closure as describe in the examples above.

* ``onComplete``: registers a function that is invoked after the `last` value is emitted by the channel.

* ``onError``: registers a function that it is invoked when an exception is raised while handling the
  ``onNext`` event. It will not make further calls to ``onNext`` or ``onComplete``.
  The ``onError`` method takes as its parameter the ``Throwable`` that caused the error.


For example::

    Channel
        .from( 1, 2, 3 )
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    2
    3
    Done



.. _glob: http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
.. _faq-page:

***
FAQ
***

How do I process multiple input files in parallel?
--------------------------------------------------

Q: *I have a collection of input files (e.g. carrots.fa, onions.fa, broccoli.fa). How can I specify that a process is performed on each input file in a parallel manner?*

A: The idea here is to create a ``channel`` that will trigger a process
execution for each of your files. First define a parameter that specifies where
the input files are:

::

    params.input = "data/*.fa"

Each of the files in the data directory can be made into a channel with:

::

    vegetable_datasets = Channel.fromPath(params.input)

From here, each time the variable ``vegetable_datasets`` is called as an
input to a process, the process will be performed on each of the files
in the vegetable datasets. For example, each input file may contain a
collection of unaligned sequences. We can specify a process to align
them as follows:

::

    process clustalw2_align {
        input:
        file vegetable_fasta from vegetable_datasets

        output:
        file "${vegetable_fasta.baseName}.aln" into vegetable_alns

        script:
        """
        clustalw2 -INFILE=${vegetable_fasta}
        """
    }

This would result in the alignment of the three vegetable fasta files
into ``carrots.aln``, ``onions.aln`` and ``broccoli.aln``.

These aligned files are now in the channel ``vegetable_alns`` and can be
used as input for a further process.

How do I get a unique ID based on the file name?
------------------------------------------------

*Q: How do I get a unique identifier based on a dataset file names (e.g. broccoli from broccoli.fa) and have the results going to a specific folder (e.g. results/broccoli/)?*

A: First we can specify a results directory as shown below:

::

    results_path = $PWD/results

The best way to manage this is to have the channel emit a tuple
containing both the file base name (``broccoli``) and the full file path
(``data/broccoli.fa``):

::

    datasets = Channel
                    .fromPath(params.input)
                    .map { file -> tuple(file.baseName, file) }

And in the process we can then reference these variables (``datasetID``
and ``datasetFile``):

::

    process clustalw2_align {
        publishDir "$results_path/$datasetID"

        input:
        set datasetID, file(datasetFile) from datasets

        output:
        set datasetID, file("${datasetID}.aln") into aligned_files

        script:
        """
        clustalw2 -INFILE=${datasetFile} -OUTFILE=${datasetID}.aln
        """
    }

In our example above would now have the folder ``broccoli`` in the results directory which would
contain the file ``broccoli.aln``.

If the input file has multiple extensions (e.g. ``brocolli.tar.gz``), you will want to use
``file.simpleName`` instead, to strip all of them (available since Nextflow 0.25+).


How do I use the same channel multiple times?
---------------------------------------------

*Q: Can a channel be used in two input statements? For example, I want carrots.fa to be aligned by both ClustalW and T-Coffee.*


A: A channel can be consumed only by one process or operator (except if channel only ever contains one item). You must
duplicate a channel before calling it as an input in different processes.
First we create the channel emitting the input files:

::

    vegetable_datasets = Channel.fromPath(params.input)

Next we can split it into two channels by using the :ref:`operator-into` operator:

::

    vegetable_datasets.into { datasets_clustalw; datasets_tcoffee }

Then we can define a process for aligning the datasets with *ClustalW*:

::

    process clustalw2_align {
        input:
        file vegetable_fasta from datasets_clustalw
        
        output:
        file "${vegetable_fasta.baseName}.aln" into clustalw_alns

        script:
        """
        clustalw2 -INFILE=${vegetable_fasta}
        """
    }

And a process for aligning the datasets with *T-Coffee*:

::

    process tcoffee_align {
        input:
        file vegetable_fasta from datasets_tcoffee
        
        output:
        file "${vegetable_fasta.baseName}.aln" into tcoffee_alns

        script:
        """
        t_coffee ${vegetable_fasta}
        """
    }

The upside of splitting the channels is that given our three unaligned
fasta files (``broccoli.fa``, ``onion.fa`` and ``carrots.fa``) six
alignment processes (three x ClustalW) + (three x T-Coffee) will be
executed as parallel processes.


How do I invoke custom scripts and tools?
-----------------------------------------

*Q: I have executables in my code, how should I call them in Nextflow?*

A: Nextflow will automatically add the directory ``bin`` into the ``PATH``
environmental variable. So therefore any executable in the ``bin``
folder of a Nextflow pipeline can be called without the need to
reference the full path.

For example, we may wish to reformat our *ClustalW* alignments from
Question 3 into *PHYLIP* format. We will use the handy tool
``esl-reformat`` for this task.

First we place copy (or create a symlink to) the ``esl-reformat``
executable to the project's bin folder. From above we see the *ClustalW*
alignments are in the channel ``clustalw_alns``:

::

    process phylip_reformat {
        input:
        file clustalw_alignment from clustalw_alns
        
        output:
        file "${clustalw_alignment.baseName}.phy" to clustalw_phylips

        script:
        """
        esl-reformat phylip ${clustalw_alignment} ${clustalw_alignment.baseName}.phy
        """
    }


    process generate_bootstrap_replicates {
        input:
        file clustalw_phylip from clustalw_phylips
        output:
        file "${clustalw_alignment.baseName}.phy" to clustalw_phylips

        script:
        """
        esl-reformat phylip ${clustalw_alignment} ${clustalw_alignment.baseName}.phy
        """
    }

How do I iterate over a process n times?
-----------------------------------------

To perform a process *n* times, we can specify the input to be
``each x from y..z``. For example:

::

    bootstrapReplicates=100

    process bootstrapReplicateTrees {
        publishDir "$results_path/$datasetID/bootstrapsReplicateTrees"

        input:
        each x from 1..bootstrapReplicates
        set val(datasetID), file(ClustalwPhylips)

        output:
        file "bootstrapTree_${x}.nwk" into bootstrapReplicateTrees

        script:
        // Generate Bootstrap Trees

        """
        raxmlHPC -m PROTGAMMAJTT -n tmpPhylip${x} -s tmpPhylip${x}
        mv "RAxML_bestTree.tmpPhylip${x}" bootstrapTree_${x}.nwk
        """
    }


How do I iterate over nth files from within a process?
------------------------------------------------------

*Q: For example, I have 100 files emitted by a channel. I wish to perform one process where I iterate over each file inside the process.*

A: The idea here to transform a channel emitting multiple items into a channel
that will collect all files into a list object and produce that list as a single emission. We do this using the ``collect()`` operator. The process script would then be able to iterate over
the files by using a simple for-loop.

This is also useful if all the items of a channel are required to be in the work directory.

::

    process concatenateBootstrapReplicates {
        publishDir "$results_path/$datasetID/concatenate"

        input:
        file bootstrapTreeList from bootstrapReplicateTrees.collect()

        output:
        file "concatenatedBootstrapTrees.nwk"

        // Concatenate Bootstrap Trees
        script:
        """
        for treeFile in ${bootstrapTreeList}
        do
            cat \$treeFile >> concatenatedBootstrapTrees.nwk
        done

        """
    }

How do I use a specific version of Nextflow?
------------------------------------------------------

*Q: I need to specify a version of Nextflow to use, or I need to pull a snapshot release.*

A: Sometimes it is necessary to use a different version of Nextflow for a specific feature or testing purposes. Nextflow is able to automatically pull versions when the ``NXF_VER`` environment variable is defined on the commandline. 

::

    NXF_VER=0.28.0 nextflow run main.nf
.. _k8s-page:

**********
Kubernetes
**********

`Kubernetes <https://kubernetes.io/>`_ is a cloud-native open-source system for deployment, scaling, and management of
containerized applications.

It provides clustering and file system abstractions that allows the execution of containerised workloads across
different cloud platforms and on-premises installations.

The built-in support for Kubernetes provided by Nextflow streamlines the execution of containerised workflows in
Kubernetes clusters.


Concepts
========

Kubernetes main abstraction is the `pod`. A `pod` defines the (desired) state of one or more containers i.e. required
computing resources, storage, network configuration.

Kubernetes abstracts also the storage provisioning through the definition of one more more persistent volumes that
allow containers to access to the underlying storage systems in a transparent and portable manner.

When using the ``k8s`` executor Nextflow deploys the workflow execution as a Kubernetes pod. This pod orchestrates
the workflow execution and submits a separate pod execution for each job that need to be carried out by the workflow
application.

.. image:: images/nextflow-k8s-min.png


Requirements
============

At least a `Persistent Volume <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#persistent-volumes>`_ with
``ReadWriteMany`` access mode has to be defined in the Kubernetes cluster (check the supported storage systems
at `this link <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#access-modes>`_).

Such volume needs to be accessible through a
`Persistent Volume Claim <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#persistentvolumeclaims>`_, which
will be used by Nextflow to run the application and store the scratch data and the pipeline final result.

The workflow application has to be containerised using the usual Nextflow :ref:`container<process-container>` directive.


Execution
=========

The workflow execution needs to be submitted from a computer able to connect to the Kubernetes cluster.

Nextflow uses the Kubernetes configuration file available at the path ``$HOME/.kube/config`` or the file specified
by the environment variable ``KUBECONFIG``.

You can verify such configuration with the command below::

    $ kubectl cluster-info
    Kubernetes master is running at https://your-host:6443
    KubeDNS is running at https://your-host:6443/api/v1/namespaces/kube-system/services/kube-dns:dns/proxy


To deploy and launch the workflow execution use the Nextflow command ``kuberun`` as shown below::

    nextflow kuberun <pipeline-name> -v vol-claim:/mount/path


This command will create and execute a pod running the nextflow orchestrator for the specified workflow.
In the above example replace ``<pipeline-name>`` with an existing nextflow project or the absolute path
of a workflow already deployed in the Kubernetes cluster.

The ``-v`` command line option is required to specify the volume claim name and mount path to use for the workflow
execution. In the above example replace ``vol-claim`` with the name of an existing persistent volume claim and
``/mount/path`` with the path where the volume is required to be mount in the container. Volume claims can also be
specified in the Nextflow configuration file, see the :ref:`Kubernetes configuration section<config-k8s>` for details.

Once the pod execution starts, the application in the foreground prints the console output produced by the running
workflow pod.

Interactive login
=================

For debugging purpose it's possible to execute a Nextflow pod and launch an interactive shell using the following command::

   nextflow kuberun login -v vol-claim:/mount/path

This command creates a pod, sets up the volume claim(s), configures the Nextflow environment and finally launches a Bash
login session.  

.. warning:: The pod is automatically destroyed once the shell session terminates. Do not use to start long running
  workflow executions in background.


Running in a pod
==================

The main convenience of the ``kuberun`` command is that it spares the user from manually creating a pod from
where the main Nextflow application is launched. In this scenario, the user environment is not containerised.

However there are scenarios in which Nextflow needs to be executed directly from a pod running in a
Kubernetes cluster. In these cases you will need to use the plain Nextflow ``run`` command and specify
the ``k8s`` executor and the required persistent volume claim in the ``nextflow.config`` file as shown below::

    process {
       executor = 'k8s'
    }

    k8s {
       storageClaimName = 'vol-claim'
       storageMountPath = '/mount/path'
       storageSubPath = '/my-data'
    }

In the above snippet replace ``vol-claim`` with the name of an existing persistent volume claim and replace
``/mount/path`` with the actual desired mount path (default: ``/workspace``) and ``storageSubPath``
with the directory in the volume to be mounted (default: ``/``).

.. warning:: The running pod must have been created with the same persistent volume claim name and mount as the
    one specified in your Nextflow configuration file.
    Note also that the ``run`` command does not support the ``-v`` option.

.. tip:: It is also possible to mount multiple volumes using the ``pod`` directive, setting such as ``k8s.pod = [ [volumeClaim: "other-pvc", mountPath: "/other" ]]``

Pod settings
============

The process :ref:`process-pod` directive allows the definition of pods specific settings, such as environment variables,
secrets and config maps when using the :ref:`k8s-executor` executor. See the :ref:`process-pod` directive for more details.

Limitation
==========

.. note::
  The ``kuberun`` command does not allow the execution of local Nextflow scripts and it's has been designed to
  provide a shortcut to simple pipeline deployment into a Kubernetes cluster.

  For stable pipeline deployment, Nextflow needs to be executed as a pod as mentioned in the `Running in a pod`_ section.
  In alternative take in consideration a managed provisioning service such as `Nextflow Tower <https://tower.nf>`_.

Advanced configuration
======================

Read :ref:`Kubernetes configuration<config-k8s>` and :ref:`executor <k8s-executor>` sections to learn more
about advanced configuration options.
.. _operator-page:

*******************
Operators
*******************

Nextflow `operators` are methods that allow you to connect channels to each other or to transform values
emitted by a channel applying some user provided rules.

Operators can be separated into seven groups:

* `Filtering operators`_
* `Transforming operators`_
* `Splitting operators`_
* `Combining operators`_
* `Forking operators`_
* `Maths operators`_
* `Other operators`_

.. note:: The operators :ref:`operator-set` and ``subscribe`` are *final* operators
  and therefore, if used, they must be the last operator in a chain of combined operators.


Filtering operators
===================

Given a channel, filtering operators allow you to select only the items that comply with a given rule.

The available filtering operators are:

* `distinct`_
* `filter`_
* `first`_
* `last`_
* `randomSample`_
* `take`_
* `unique`_
* `until`_

filter
---------

The ``filter`` operator allows you to get only the items emitted by a channel that satisfy a condition and discarding
all the others. The filtering condition can be specified by using either a :ref:`regular expression <script-regexp>`,
a literal value, a type `qualifier` (i.e. a Java class) or any boolean `predicate`.

The following example shows how to filter a channel by using a regular expression that returns only strings that
begin with ``a``::

    Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .filter( ~/^a.*/ )
        .view()

::

    a
    aa


The following example shows how to filter a channel by specifying the type qualifier ``Number`` so that only numbers
are returned::

    Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .filter( Number )
        .view()

::

    3
    4.5




Finally, a filtering condition can be defined by using any a boolean `predicate`. A predicate is expressed by
a :ref:`closure <script-closure>` returning a boolean value. For example the following fragment shows how filter
a channel emitting numbers so that the `odd` values are returned::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .filter { it % 2 == 1 }
        .view()

::

    1
    3
    5


.. tip:: In the above example the filter condition is wrapped in curly brackets,
  instead of round brackets, since it specifies a :ref:`closure <script-closure>` as the operator's argument.
  This is just a language syntax-sugar for ``filter({ it % 2 == 1 })``


unique
---------

The ``unique`` operator allows you to remove duplicate items from a channel and only emit single items with no repetition.

For example::

    Channel
        .from( 1,1,1,5,7,7,7,3,3 )
        .unique()
        .view()

::

    1
    5
    7
    3


You can also specify an optional :ref:`closure <script-closure>` that customizes the way it distinguishes between unique items.
For example::

    Channel
        .from(1,3,4,5)
        .unique { it % 2 }
        .view()

::

    1
    4


distinct
-----------

The ``distinct`` operator allows you to remove `consecutive` duplicated items from a channel, so that each emitted item
is different from the preceding one. For example::


    Channel
        .from( 1,1,2,2,2,3,1,1,2,2,3 )
        .distinct()
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    2
    3
    1
    2
    3
    Done



You can also specify an optional :ref:`closure <script-closure>` that customizes the way it distinguishes between distinct items.
For example::

    Channel
        .from( 1,1,2,2,2,3,1,1,2,4,6 )
        .distinct { it % 2 }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }


::

    1
    2
    3
    2
    Done


.. _operator-first:

first
--------

The ``first`` operator creates a channel that returns the first item emitted by the source channel, or eventually
the first item that matches an optional condition. The condition can be specified by using a :ref:`regular expression<script-regexp>`,
a Java `class` type or any boolean `predicate`. For example::


    // no condition is specified, emits the very first item: 1
    Channel
        .from( 1, 2, 3 )
        .first()
        .view()


    // emits the first String value: 'a'
    Channel
        .from( 1, 2, 'a', 'b', 3 )
        .first( String )
        .view()

    // emits the first item matching the regular expression: 'aa'
    Channel
        .from( 'a', 'aa', 'aaa' )
        .first( ~/aa.*/ )
        .view()

    // emits the first item for which the predicate evaluates to true: 4
    Channel
        .from( 1,2,3,4,5 )
        .first { it > 3 }
        .view()


randomSample
------------

The ``randomSample`` operator allows you to create a channel emitting the specified number of items randomly taken
from the channel to which is applied. For example::

  Channel
        .from( 1..100 )
        .randomSample( 10 )
        .view()

The above snippet will print 10 numbers in the range from 1 to 100.

The operator supports a second parameter that allows you to set the initial `seed` for the random number generator.
By setting it, the ``randomSample`` operator will always return the same pseudo-random sequence. For example::

  Channel
        .from( 1..100 )
        .randomSample( 10, 234 )
        .view()

The above example will print 10 random numbers in the range between 1 and 100. At each run of the script, the same 
sequence will be returned.

take
-------

The ``take`` operator allows you to filter only the first `n` items emitted by a channel. For example::

    Channel
        .from( 1,2,3,4,5,6 )
        .take( 3 )
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    2
    3
    Done

.. note:: By specifying the value ``-1`` the operator takes all values.

See also `until`_.

.. _operator-last:

last
-------

The ``last`` operator creates a channel that only returns the last item emitted by the source channel. For example::

    Channel
        .from( 1,2,3,4,5,6 )
        .last()
        .view()

::

    6


until
-----

The ``until`` operator creates a channel that returns the items emitted by the source channel and stop when
the condition specified is verified. For example::

  Channel
      .from( 3,2,1,5,1,5 )
      .until{ it==5 }
      .view()

::

  3
  2
  1

See also `take`_. 

Transforming operators
======================

Transforming operators are used to transform the items emitted by a channel to new values.

These operators are:

* `buffer`_
* `collate`_
* `collect`_
* `flatten`_
* `flatMap`_
* `groupBy`_
* `groupTuple`_
* `map`_
* `reduce`_
* `toList`_
* `toSortedList`_
* `transpose`_

.. _operator-map:

map
------

The ``map`` operator applies a function of your choosing to every item emitted by a channel, and 
returns the items so obtained as a new channel. The function applied is called the `mapping` function 
and is expressed with a :ref:`closure <script-closure>` as shown in the example below::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .map { it * it }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    1
    4
    9
    16
    25
    Done


.. _operator-flatmap:

flatMap
----------

The ``flatMap`` operator applies a function of your choosing to every item emitted by a channel, and
returns the items so obtained as a new channel. Whereas the `mapping` function returns a list of items,
this list is flattened so that each single item is emitted on its own.  

For example::

    // create a channel of numbers
    numbers = Channel.from( 1, 2, 3 )

    // map each number to a tuple (array), which items are emitted separately
    results = numbers.flatMap { n -> [ n*2, n*3 ] }

    // print the final results
    results.subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    2
    3
    4
    6
    6
    9
    Done


Associative arrays are handled in the same way, so that each array entry is emitted as a single `key-value` item. For example::

    Channel.from ( 1, 2, 3 )
           .flatMap { it -> [ number: it, square: it*it ] }
           .view { it.key + ': ' + it.value }

::

    number: 1
    square: 1
    number: 2
    square: 4
    number: 3
    square: 9


.. _operator-reduce:

reduce
---------

The ``reduce`` operator applies a function of your choosing to every item emitted by a channel.
Each time this function is invoked it takes two parameters: firstly the `i-th` emitted item
and secondly the result of the previous invocation of the function itself. The result is 
passed on to the next function call, along with the `i+1 th` item, until all the items are 
processed.

Finally, the ``reduce`` operator emits the result of the last invocation of your function 
as the sole output.

For example::

    Channel
        .from( 1, 2, 3, 4, 5 )
        .reduce { a, b -> println "a: $a b: $b"; return a+b }
        .view { "result = $it" }


It prints the following output::

	a: 1	b: 2
	a: 3	b: 3
	a: 6	b: 4
	a: 10	b: 5
	result = 15


.. note:: In a common usage scenario the first function parameter is used as an `accumulator` and
  the second parameter represents the `i-th` item to be processed.

Optionally you can specify a `seed` value in order to initialise the accumulator parameter
as shown below::

    myChannel.reduce( seedValue ) {  a, b -> ... }



groupBy
----------

.. warning::
    This operator is deprecated. Use the `groupTuple`_ operator instead.

The ``groupBy`` operator collects the values emitted by the source channel grouping them together using a `mapping`
function that associates each item with a key. When finished, it emits an associative
array that maps each key to the set of items identified by that key.  

For example::

    Channel
    	.from('hello','ciao','hola', 'hi', 'bonjour')
    	.groupBy { String str -> str[0] } 
    	.view()

:: 

    [ b:['bonjour'], c:['ciao'], h:['hello','hola','hi'] ]
    

The `mapping` function is an optional parameter. When omitted, the values are grouped
according to these rules:

* Any value of type ``Map`` is associated with the value of its first entry, or ``null`` when the map itself is empty.
* Any value of type ``Map.Entry`` is associated with the value of its ``key`` attribute.
* Any value of type ``Collection`` or ``Array`` is associated with its first entry.
* For any other value, the value itself is used as a key.


.. _operator-grouptuple:

groupTuple
----------

The ``groupTuple`` operator collects tuples (or lists) of values emitted by the source channel grouping together the
elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

In other words, the operator transforms a sequence of tuple like *(K, V, W, ..)* into a new channel emitting a sequence of
*(K, list(V), list(W), ..)*

For example::

   Channel
        .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
        .groupTuple()
        .view()

It prints::

    [1, [A, B, C]]
    [2, [C, A]]
    [3, [B, D]]

By default the first entry in the tuple is used as grouping key. A different key can be chosen by using the
``by`` parameter and specifying the index of the entry to be used as key (the index is zero-based). For example,
grouping by the second value in each tuple::

   Channel
        .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
        .groupTuple(by: 1)
        .view()

The result is::

    [[1, 2], A]
    [[1, 3], B]
    [[2, 1], C]
    [[3], D]


Available parameters:

=========== ============================
Field       Description
=========== ============================
by          The index (zero based) of the element to be used as grouping key.
            A key composed by multiple elements can be defined specifying a list of indices e.g. ``by: [0,2]``
sort        Defines the sorting criteria for the grouped items. See below for available sorting options.
size        The number of items the grouped list(s) has to contain. When the specified size is reached, the tuple is emitted.
remainder   When ``false`` incomplete tuples (i.e. with less than `size` grouped items)
            are discarded (default). When ``true`` incomplete tuples are emitted as the ending emission. Only valid when a ``size`` parameter
            is specified.
=========== ============================

Sorting options:

=============== ========================
Sort            Description
=============== ========================
false           No sorting is applied (default).
true            Order the grouped items by the item natural ordering i.e. numerical for number, lexicographic for string, etc. See http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html
hash            Order the grouped items by the hash number associated to each entry.
deep            Similar to the previous, but the hash number is created on actual entries content e.g. when the item is a file, the hash is created on the actual file content.
`custom`        A custom sorting criteria used to order the tuples element holding list of values. It can be specified by using either a :ref:`Closure <script-closure>` or a `Comparator <http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html>`_ object.
=============== ========================


.. tip:: You should always specify the number of expected elements in each tuple using the ``size`` attribute
  to allow the ``groupTuple`` operator to stream the collected values as soon as possible. However there
  are use cases in which each tuple has a different size depending on the grouping key. In this case use the
  built-in function ``groupKey`` that allows you to create a special grouping key object to which it's possible
  to associate the group size for a given key.


buffer
---------

The ``buffer`` operator gathers the items emitted by the source channel into subsets and emits these subsets separately.


There are a number of ways you can regulate how ``buffer`` gathers the items from
the source channel into subsets:

* ``buffer( closingCondition )``: starts to collect the items emitted by the channel into 
  a subset until the `closing condition` is verified. After that the subset is emitted 
  to the resulting channel and new items are gathered into a new subset. The process is repeated 
  until the last value in the source channel is sent. The ``closingCondition`` can be specified 
  either as a :ref:`regular expression <script-regexp>`, a Java class, a literal value, or a `boolean predicate`
  that has to be satisfied. For example::
  
    Channel
        .from( 1,2,3,1,2,3 ) 
        .buffer { it == 2 } 
        .view()

    // emitted values
    [1,2]
    [3,1,2]
  
  

* ``buffer( openingCondition, closingCondition )``: starts to gather the items emitted by the channel 
  as soon as one of the them verify the `opening condition` and it continues until there is one item
  which verify the `closing condition`. After that the subset is emitted and it continues applying the 
  described logic until the last channel item is emitted.
  Both conditions can be defined either as a :ref:`regular expression <script-regexp>`, a literal value,
  a Java class, or a `boolean predicate` that need to be satisfied. For example:: 
 
    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( 2, 4 ) 
        .view()

    // emits bundles starting with '2' and ending with'4'
    [2,3,4]
    [2,3,4]      
  

* ``buffer( size: n )``: transform the source channel in such a way that it emits tuples 
  made up of `n` elements. An incomplete tuple is discarded. For example::

    Channel
        .from( 1,2,3,1,2,3,1 ) 
        .buffer( size: 2 )
        .view()
        
    // emitted values 
    [1, 2]
    [3, 1]
    [2, 3]

If you want to emit the last items in a tuple containing less than `n` elements, simply 
add the parameter ``remainder`` specifying ``true``, for example::

    Channel
        .from( 1,2,3,1,2,3,1 )
        .buffer( size: 2, remainder: true )
        .view()

    // emitted values
    [1, 2]
    [3, 1]
    [2, 3]
    [1]



* ``buffer( size: n, skip: m )``: as in the previous example, it emits tuples containing `n` elements, 
  but skips ``m`` values before starting to collect the values for the next tuple (including the first emission). For example::

    Channel
        .from( 1,2,3,4,5,1,2,3,4,5,1,2 ) 
        .buffer( size:3, skip:2 )
        .view()
        
    // emitted values 
    [3, 4, 5]
    [3, 4, 5]

If you want to emit the remaining items in a tuple containing less than `n` elements, simply
add the parameter ``remainder`` specifying ``true``, as shown in the previous example.

See also: `collate`_ operator.


collate
---------

The ``collate`` operator transforms a channel in such a way that the emitted values are grouped in tuples containing `n` items. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3 )
        .view()

::

        [1, 2, 3]
        [1, 2, 3]
        [1]

As shown in the above example the last tuple may be incomplete e.g. contain fewer elements than the specified size.
If you want to avoid this, specify ``false`` as the second parameter. For example::

    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3, false )
        .view()

::

        [1, 2, 3]
        [1, 2, 3]


A second version of the ``collate`` operator allows you to specify, after the `size`, the `step` by which elements
are collected in tuples. For example::

    Channel
      .from(1,2,3,4)
      .collate( 3, 1 )
      .view()

::

    [1, 2, 3]
    [2, 3, 4]
    [3, 4]
    [4]

As before, if you don't want to emit the last items which do not complete a tuple, specify ``false`` as the third parameter.


See also: `buffer`_ operator.

.. _operator-collect:

collect
-------

The ``collect`` operator collects all the items emitted by a channel to a ``List`` and return
the resulting object as a sole emission. For example::

    Channel
        .from( 1, 2, 3, 4 )
        .collect()
        .view()

    # outputs
    [1,2,3,4]

An optional :ref:`closure <script-closure>` can be specified to transform each item before adding it to the resulting list.
For example::

    Channel
        .from( 'hello', 'ciao', 'bonjour' )
        .collect { it.length() }
        .view()

    # outputs
    [5,4,7]

.. Available parameters:
..
.. =========== ============================
.. Field       Description
.. =========== ============================
.. flat        When ``true`` nested list structures are normalised and their items are added to the resulting list object (default: ``true``).
.. sort        When ``true`` the items in the resulting list are sorted by their natural ordering. It is possible to provide a custom ordering criteria by using either a :ref:`closure <script-closure>` or a `Comparator <https://docs.oracle.com/javase/8/docs/api/java/util/Comparator.html>`_ object (default: ``false``).
.. =========== ============================

See also: `toList`_ and `toSortedList`_ operator.

.. _operator-flatten:

flatten
----------

The ``flatten`` operator transforms a channel in such a way that every item of type ``Collection`` or ``Array``
is flattened so that each single entry is emitted separately by the resulting channel. For example::

    Channel
    	.from( [1,[2,3]], 4, [5,[6]] )
    	.flatten()
    	.view()

:: 
    
    1
    2
    3
    4
    5
    6
    
    
See also: `flatMap`_ operator.



toList
---------

The ``toList`` operator collects all the items emitted by a channel to a ``List`` object
and emits the resulting collection as a single item. For example::

    Channel
    	.from( 1, 2, 3, 4 )
    	.toList() 
    	.subscribe onNext: { println it }, onComplete: { println 'Done' }
    	
::
 
    [1,2,3,4]
    Done

See also: `collect`_ operator.

toSortedList
---------------


The ``toSortedList`` operator collects all the items emitted by a channel to a ``List`` object where they are sorted
and emits the resulting collection as a single item. For example::

    Channel
    	.from( 3, 2, 1, 4 )
    	.toSortedList()
    	.subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    [1,2,3,4]
    Done

You may also pass a comparator closure as an argument to the ``toSortedList`` operator to customize the sorting criteria.  For example, to sort by the second element of a tuple in descending order::

    Channel
        .from( ["homer", 5], ["bart", 2], ["lisa", 10], ["marge", 3], ["maggie", 7])
        .toSortedList( { a, b -> b[1] <=> a[1] } )
        .view()

::

   [[lisa, 10], [maggie, 7], [homer, 5], [marge, 3], [bart, 2]]

See also: `collect`_ operator.

transpose
---------

The ``transpose`` operator transforms a channel in such a way that the emitted items are the result of a transposition
of all tuple elements in each item. For example::

    Channel.from([
       ['a', ['p', 'q'], ['u','v'] ],
       ['b', ['s', 't'], ['x','y'] ]
       ])
       .transpose()
       .view()

The above snippet prints::

    [a, p, u]
    [a, q, v]
    [b, s, x]
    [b, t, y]


Available parameters:

=========== ============================
Field       Description
=========== ============================
by          The index (zero based) of the element to be transposed.
            Multiple elements can be defined specifying as list of indices e.g. ``by: [0,2]``
remainder   When ``false`` incomplete tuples are discarded (default). When ``true`` incomplete tuples are emitted
            containing a `null` in place of a missing element.
=========== ============================


Splitting operators
====================

These operators are used to split items emitted by channels into chunks that can be processed by downstream
operators or processes.

The available splitting operators are:

* `splitCsv`_
* `splitFasta`_
* `splitFastq`_
* `splitText`_


splitCsv
---------

The ``splitCsv`` operator allows you to parse text items emitted by a channel, that are formatted using the
`CSV format <http://en.wikipedia.org/wiki/Comma-separated_values>`_, and split them into records or group them into
list of records with a specified length.

In the simplest case just apply the ``splitCsv`` operator to a channel emitting a CSV formatted text files or
text entries. For example::

    Channel
        .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
        .splitCsv()
        .view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }

The above example shows hows CSV text is parsed and is split into single rows. Values can be accessed
by its column index in the row object.

When the CSV begins with a header line defining the column names, you can specify the parameter ``header: true`` which
allows you to reference each value by its name, as shown in the following example::

    Channel
        .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
        .splitCsv(header: true)
        .view { row -> "${row.alpha} - ${row.beta} - ${row.gamma}" }

It will print ::

 10 - 20 - 30
 70 - 80 - 90

Alternatively you can provide custom header names by specifying a the list of strings in the ``header`` parameter
as shown below::


    Channel
        .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
        .splitCsv(header: ['col1', 'col2', 'col3'], skip: 1 )
        .view { row -> "${row.col1} - ${row.col2} - ${row.col3}" }


Available parameters:

=========== ============================
Field       Description
=========== ============================
by          The number of rows in each `chunk`
sep         The character used to separate the values (default: ``,``)
quote       Values may be quoted by single or double quote characters.
header      When ``true`` the first line is used as columns names. Alternatively it can be used to provide the list of columns names.
charset     Parse the content by using the specified charset e.g. ``UTF-8``
strip       Removes leading and trailing blanks from values (default: ``false``)
skip        Number of lines since the file beginning to ignore when parsing the CSV content.
limit       Limits the number of retrieved records for each file to the specified value.
decompress  When ``true`` decompress the content using the GZIP format before processing it (note: files whose name ends with ``.gz`` extension are decompressed automatically)
elem        The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)
=========== ============================


splitFasta
------------

The ``splitFasta`` operator allows you to split the entries emitted by a channel, that are formatted using the
`FASTA format <http://en.wikipedia.org/wiki/FASTA_format>`_. It returns a channel which emits text item
for each sequence in the received FASTA content.

The number of sequences in each text chunk produced by the ``splitFasta`` operator can be set by using
the ``by`` parameter. The following example shows how to read a FASTA file and split it into chunks containing 10 sequences
each::

   Channel
        .fromPath('misc/sample.fa')
        .splitFasta( by: 10 )
        .view()

.. warning:: By default chunks are kept in memory. When splitting big files specify the parameter ``file: true`` to save the
  chunks into files in order to not incur in a ``OutOfMemoryException``. See the available parameter table below for details.

A second version of the ``splitFasta`` operator allows you to split a FASTA content into record objects, instead
of text chunks. A record object contains a set of fields that let you access and manipulate the FASTA sequence
information with ease.


In order to split a FASTA content into record objects, simply use the ``record`` parameter specifying the map of
required the fields, as shown in the example below::

   Channel
        .fromPath('misc/sample.fa')
        .splitFasta( record: [id: true, seqString: true ])
        .filter { record -> record.id =~ /^ENST0.*/ }
        .view { record -> record.seqString }


.. note:: In this example, the file ``misc/sample.fa`` is split into records containing the ``id`` and the ``seqString`` fields
  (i.e. the sequence id and the sequence data). The following ``filter`` operator only keeps the sequences which ID
  starts with the ``ENST0`` prefix, finally the sequence content is printed by using the ``subscribe`` operator.

Available parameters:

=========== ============================
Field       Description
=========== ============================
by          Defines the number of sequences in each `chunk` (default: ``1``)
size        Defines the size in memory units of the expected chunks eg. `1.MB`.
limit       Limits the number of retrieved sequences for each file to the specified value.
record      Parse each entry in the FASTA file as record objects (see following table for accepted values)
charset     Parse the content by using the specified charset e.g. ``UTF-8``
compress    When ``true`` resulting file chunks are GZIP compressed. The ``.gz`` suffix is automatically added to chunk file names.
decompress  When ``true``, decompress the content using the GZIP format before processing it (note: files whose name ends with ``.gz`` extension are decompressed automatically)
file        When ``true`` saves each split to a file. Use a string instead of ``true`` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.
elem        The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)
=========== ============================


The following fields are available when using the ``record`` parameter:

=========== ============================
Field       Description
=========== ============================
id          The FASTA sequence identifier i.e. the word following the ``>`` symbol up to the first `blank` or `newline` character
header      The first line in a FASTA sequence without the ``>`` character
desc        The text in the FASTA header following the ID value
text        The complete FASTA sequence including the header
seqString   The sequence data as a single line string i.e. containing no `newline` characters
sequence    The sequence data as a multi-line string (always ending with a `newline` character)
width       Define the length of a single line when the ``sequence`` field is used, after that the sequence data continues on a new line.
=========== ============================



splitFastq
----------

The ``splitFastq`` operator allows you to split the entries emitted by a channel, that are formatted using the
`FASTQ format <http://en.wikipedia.org/wiki/FASTQ_format>`_. It returns a channel which emits a text chunk
for each sequence in the received item.

The number of sequences in each text chunk produced by the ``splitFastq`` operator is defined by the
parameter ``by``. The following example shows you how to read a FASTQ file and split it into chunks containing 10
sequences each::

   Channel
        .fromPath('misc/sample.fastq')
        .splitFastq( by: 10 )
        .view()


.. warning:: By default chunks are kept in memory. When splitting big files specify the parameter ``file: true`` to save the
  chunks into files in order to not incur in a ``OutOfMemoryException``. See the available parameter table below for details.


A second version of the ``splitFastq`` operator allows you to split a FASTQ formatted content into record objects,
instead of text chunks. A record object contains a set of fields that let you access and manipulate the FASTQ sequence
data with ease.

In order to split FASTQ sequences into record objects simply use the ``record`` parameter specifying the map of
the required fields, or just specify ``record: true`` as in the example shown below::

   Channel
        .fromPath('misc/sample.fastq')
        .splitFastq( record: true )
        .view { record -> record.readHeader }


Finally the ``splitFastq`` operator is able to split paired-end read pair FASTQ files. It must be applied to a channel
which emits tuples containing at least two elements that are the files to be splitted. For example::

    Channel
        .fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat:true)
        .splitFastq(by: 100_000, pe:true, file:true)
        .view()


.. note:: The ``fromFilePairs`` requires the ``flat:true`` option to have the file pairs as separate elements
  in the produced tuples.

.. warning:: This operator assumes that the order of the PE reads correspond with each other and both files contain
  the same number of reads.


Available parameters:

=========== ============================
Field       Description
=========== ============================
by          Defines the number of *reads* in each `chunk` (default: ``1``)
pe          When ``true`` splits paired-end read files, therefore items emitted by the source channel must be tuples in which at least two elements are the read-pair files to be splitted.
limit       Limits the number of retrieved *reads* for each file to the specified value.
record      Parse each entry in the FASTQ file as record objects (see following table for accepted values)
charset     Parse the content by using the specified charset e.g. ``UTF-8``
compress    When ``true`` resulting file chunks are GZIP compressed. The ``.gz`` suffix is automatically added to chunk file names.
decompress  When ``true`` decompress the content using the GZIP format before processing it (note: files whose name ends with ``.gz`` extension are decompressed automatically)
file        When ``true`` saves each split to a file. Use a string instead of ``true`` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.
elem        The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)
=========== ============================

The following fields are available when using the ``record`` parameter:

=============== ============================
Field           Description
=============== ============================
readHeader      Sequence header (without the ``@`` prefix)
readString      The raw sequence data
qualityHeader   Base quality header (it may be empty)
qualityString   Quality values for the sequence
=============== ============================

splitText
----------

The ``splitText`` operator allows you to split multi-line strings or text file items, emitted by a source channel
into chunks containing `n` lines, which will be emitted by the resulting channel.

For example::

   Channel
        .fromPath('/some/path/*.txt')
        .splitText()
        .view()


It splits the content of the files with suffix ``.txt``, and prints it line by line.

By default the ``splitText`` operator splits each item into chunks of one line. You can define the number of lines in each chunk by using
the parameter ``by``, as shown in the following example::


   Channel
        .fromPath('/some/path/*.txt')
        .splitText( by: 10 )
        .subscribe {
            print it;
            print "--- end of the chunk ---\n"
        }


An optional :ref:`closure <script-closure>` can be specified in order to `transform` the text chunks produced by the operator.
The following example shows how to split text files into chunks of 10 lines and transform them to capital letters::

     Channel
        .fromPath('/some/path/*.txt')
        .splitText( by: 10 ) { it.toUpperCase() }
        .view()


.. note:: Text chunks returned by the operator ``splitText`` are always terminated by a ``newline`` character.


Available parameters:

=========== ============================
Field       Description
=========== ============================
by          Defines the number of lines in each `chunk` (default: ``1``).
limit       Limits the number of retrieved lines for each file to the specified value.
charset     Parse the content by using the specified charset e.g. ``UTF-8``.
compress    When ``true`` resulting file chunks are GZIP compressed. The ``.gz`` suffix is automatically added to chunk file names.
decompress  When ``true``, decompress the content using the GZIP format before processing it (note: files whose name ends with ``.gz`` extension are decompressed automatically).
file        When ``true`` saves each split to a file. Use a string instead of ``true`` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in oder to save the split files into the specified folder.
elem        The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element).
keepHeader  Parses the first line as header and prepends it to each emitted chunk.
=========== ============================


Combining operators
=====================

The combining operators are:

* `cross`_
* `collectFile`_
* `combine`_
* `concat`_
* `join`_
* `merge`_
* `mix`_
* `phase`_
* `spread`_
* `tap`_


.. _operator-join:

join
-----

The ``join`` operator creates a channel that joins together the items emitted by two channels for which exists
a matching key. The key is defined, by default, as the first element in each item emitted.

For example::

  left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
  right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
  left.join(right).view()

The resulting channel emits::

  [Z, 3, 6]
  [Y, 2, 5]
  [X, 1, 4]

The `index` of a different matching element can be specified by using the ``by`` parameter.

The ``join`` operator can emit all the pairs that are incomplete, i.e. the items for which a matching element
is missing, by specifying the optional parameter ``remainder`` as shown below::

    left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
    right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
    left.join(right, remainder: true).view()

The above example prints::

    [Y, 2, 5]
    [Z, 3, 6]
    [X, 1, 4]
    [P, 7, null]


The following parameters can be used with the ``join`` operator:

=============== ========================
Name            Description
=============== ========================
by              The index (zero based) of the element to be used as grouping key.
                A key composed by multiple elements can be defined specifying a list of indices e.g. ``by: [0,2]``
remainder       When ``false`` incomplete tuples (i.e. with less than `size` grouped items)
                are discarded (default). When ``true`` incomplete tuples are emitted as the ending emission.
failOnDuplicate An error is reported when the same key is found more than once.
failOnMismatch  An error is reported when a channel emits a value for which there isn't a corresponding element in the joining channel. This option cannot be used with ``remainder``.
=============== ========================


.. _operator-merge:

merge
--------

.. warning::
    This operator is deprecated and it will be removed in upcoming release.

The ``merge`` operator lets you join items emitted by two (or more) channels into a new channel.

For example the following code merges two channels together, one which emits a series of odd integers
and the other which emits a series of even integers::

    odds  = Channel.from([1, 3, 5, 7, 9]);
    evens = Channel.from([2, 4, 6]);

    odds
        .merge( evens )
        .view()

::

    [1, 2]
    [3, 4]
    [5, 6]

An option closure can be provide to customise the items emitted by the resulting merged channel. For example::

    odds  = Channel.from([1, 3, 5, 7, 9]);
    evens = Channel.from([2, 4, 6]);

    odds
        .merge( evens ) { a, b -> tuple(b*b, a) }
        .view()

.. _operator-mix:

mix
------

The ``mix`` operator combines the items emitted by two (or more) channels into a single channel.


For example::

        c1 = Channel.from( 1,2,3 )
        c2 = Channel.from( 'a','b' )
        c3 = Channel.from( 'z' )

        c1 .mix(c2,c3)
           .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

        1
        2
        3
        'a'
        'b'
        'z'

.. note:: The items emitted by the resulting mixed channel may appear in any order,
  regardless of which source channel they came from. Thus, the following example
  it could be a possible result of the above example as well.

::

          'z'
          1
          'a'
          2
          'b'
          3


.. _operator-phase:

phase
--------

.. warning:: This operator is deprecated. Use the `join`_ operator instead.

The ``phase`` operator creates a channel that synchronizes the values emitted by two other channels,
in such a way that it emits pairs of items that have a matching key.

The key is defined, by default, as the first entry in an array, a list or map object,
or the value itself for any other data type.

For example::

        ch1 = Channel.from( 1,2,3 )
        ch2 = Channel.from( 1,0,0,2,7,8,9,3 )
        ch1 .phase(ch2) .view()

It prints::

    [1,1]
    [2,2]
    [3,3]


Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key,
as shown in the following example::


    ch1 = Channel.from( [sequence: 'aaaaaa', id: 1], [sequence: 'bbbbbb', id: 2] )
    ch2 = Channel.from( [val: 'zzzz', id: 3], [val: 'xxxxx', id: 1], [val: 'yyyyy', id: 2])
    ch1 .phase(ch2) { it -> it.id } .view()


It prints::

    [[sequence:aaaaaa, id:1], [val:xxxxx, id:1]]
    [[sequence:bbbbbb, id:2], [val:yyyyy, id:2]]


Finally, the ``phase`` operator can emit all the pairs that are incomplete, i.e. the items for which a matching element
is missing, by specifying the optional parameter ``remainder`` as shown below::

        ch1 = Channel.from( 1,0,0,2,5,3 )
        ch2 = Channel.from( 1,2,3,4 )
        ch1 .phase(ch2, remainder: true) .view()

It prints::

    [1, 1]
    [2, 2]
    [3, 3]
    [0, null]
    [0, null]
    [5, null]
    [null, 4]

See also `join`_ operator.

.. _operator-cross:

cross
-------

The ``cross`` operators allows you to combine the items of two channels in such a way that
the items of the source channel are emitted along with the items emitted by the target channel 
for which they have a matching key.  

The key is defined, by default, as the first entry in an array, a list or map object,
or the value itself for any other data type. For example:: 

	source = Channel.from( [1, 'alpha'], [2, 'beta'] )
	target = Channel.from( [1, 'x'], [1, 'y'], [1, 'z'], [2,'p'], [2,'q'], [2,'t'] )

	source.cross(target).view()

It will output:: 

	[ [1, alpha], [1, x] ]
	[ [1, alpha], [1, y] ]
	[ [1, alpha], [1, z] ]
	[ [2, beta],  [2, p] ]
	[ [2, beta],  [2, q] ]
	[ [2, beta],  [2, t] ]

The above example shows how the items emitted by the source channels are associated to the ones
emitted by the target channel (on the right) having the same key. 

There are two important caveats when using the ``cross`` operator:

	#. The operator is not `commutative`, i.e. the result of ``a.cross(b)`` is different from ``b.cross(a)`` 
	#. The source channel should emits items for which there's no key repetition i.e. the emitted 
	   items have an unique key identifier. 

Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key,
in a similar manner as shown for the `phase`_ operator.

collectFile
-----------

The ``collectFile`` operator allows you to gather the items emitted by a channel and save them to one or more files.
The operator returns a new channel that emits the collected file(s).

In the simplest case, just specify the name of a file where the entries have to be stored. For example::

    Channel
        .from('alpha', 'beta', 'gamma')
        .collectFile(name: 'sample.txt', newLine: true)
        .subscribe {
            println "Entries are saved to file: $it"
            println "File content is: ${it.text}"
        }



A second version of the ``collectFile`` operator allows you to gather the items emitted by a channel and group them together
into files whose name can be defined by a dynamic criteria. The grouping criteria is specified by a :ref:`closure <script-closure>`
that must return a pair in which the first element defines the file name for the group and the second element the actual
value to be appended to that file. For example::

     Channel
        .from('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
        .collectFile() { item ->
            [ "${item[0]}.txt", item + '\n' ]
        }
        .subscribe {
            println "File ${it.name} contains:"
            println it.text
        }

It will print::

    File 'B.txt' contains:
    Bonjour

    File 'C.txt' contains:
    Ciao

    File 'H.txt' contains:
    Halo
    Hola
    Hello


.. tip:: When the items emitted by the source channel are files, the grouping criteria can be omitted. In this case
  the items content will be grouped in file(s) having the same name as the source items.


The following parameters can be used with the ``collectFile`` operator:

=============== ========================
Name            Description
=============== ========================
``cache``       Controls the caching ability of the ``collectFile`` operator when using the *resume* feature. It follows the same semantic of the :ref:`process-cache` directive (default: ``true``).
``keepHeader``  Prepend the resulting file with the header fetched in the first collected file. The header size (ie. lines) can be specified by using the ``skip`` parameter (default: ``false``), to determine how many lines to remove from all collected files except for the first (where no lines will be removed).
``name``        Name of the file where all received values are stored.
``newLine``     Appends a ``newline`` character automatically after each entry (default: ``false``).
``seed``        A value or a map of values used to initialise the files content.
``skip``        Skip the first `n` lines eg. ``skip: 1``.
``sort``        Defines sorting criteria of content in resulting file(s). See below for sorting options.
``storeDir``    Folder where the resulting file(s) are be stored.
``tempDir``     Folder where temporary files, used by the collecting process, are stored.
=============== ========================

.. note:: The file content is sorted in such a way that it does not depend on the order on which
    entries have been added to it, this guarantees that it is consistent (i.e. do not change) across different executions
    with the same data.

The ordering of file's content can be defined by using the ``sort`` parameter. The following criteria
can be specified:

=============== ========================
Sort            Description
=============== ========================
``false``       Disable content sorting. Entries are appended as they are produced.
``true``        Order the content by the entries natural ordering i.e. numerical for number, lexicographic for string, etc. See http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html
``'index'``     Order the content by the incremental index number assigned to each entry while they are collected.
``'hash'``      Order the content by the hash number associated to each entry (default)
``'deep'``      Similar to the previous, but the hash number is created on actual entries content e.g. when the entry is a file the hash is created on the actual file content.
`custom`        A custom sorting criteria can be specified by using either a :ref:`Closure <script-closure>` or a `Comparator <http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html>`_ object.
=============== ========================

For example the following snippet shows how sort the content of the result file alphabetically::

     Channel
        .from('Z'..'A')
        .collectFile(name:'result', sort: true, newLine: true)
        .view { it.text }

It will print::

        A
        B
        C
        :
        Z


The following example shows how use a `closure` to collect and sort all sequences in a FASTA file from shortest to longest::

    Channel
         .fromPath('/data/sequences.fa')
         .splitFasta( record: [id: true, sequence: true] )
         .collectFile( name:'result.fa', sort: { it.size() } )  {
            it.sequence
          }
         .view { it.text }


.. warning:: The ``collectFile`` operator needs to store files in a temporary folder that is automatically deleted on 
  job completion. For performance reasons this folder is located in the machine's local storage,
 and it will require as much free space as are the data you are collecting. Optionally, an alternative temporary data
 folder can be specified by using the ``tempDir`` parameter.

.. _operator-combine:

combine
-------

The ``combine`` operator combines (cartesian product) the items emitted by two channels or by a channel and a ``Collection``
object (as right operand). For example::

    numbers = Channel.from(1,2,3)
    words = Channel.from('hello', 'ciao')
    numbers
        .combine(words)
        .view()

    # outputs
    [1, hello]
    [2, hello]
    [3, hello]
    [1, ciao]
    [2, ciao]
    [3, ciao]

A second version of the ``combine`` operator allows you to combine between them those items that share a common
matching key. The index of the key element is specified by using the ``by`` parameter (the index is zero-based,
multiple indexes can be specified with list a integers).
For example::

    left = Channel.from(['A',1], ['B',2], ['A',3])
    right = Channel.from(['B','x'], ['B','y'], ['A','z'], ['A', 'w'])

    left
        .combine(right, by: 0)
        .view()

    # outputs
    [A, 1, z]
    [A, 3, z]
    [A, 1, w]
    [A, 3, w]
    [B, 2, x]
    [B, 2, y]


See also `join`_, `cross`_, `spread`_ and `phase`_.

.. _operator-concat:

concat
--------

The ``concat`` operator allows you to `concatenate` the items emitted by two or more channels to a new channel, in such
a way that the items emitted by the resulting channel are in same order as they were when specified as operator arguments.

In other words it guarantees that given any `n` channels, the concatenation channel emits the items proceeding from the channel `i+1 th`
only after `all` the items proceeding from the channel `i th` were emitted.

For example::

    a = Channel.from('a','b','c')
    b = Channel.from(1,2,3)
    c = Channel.from('p','q')

    c.concat( b, a ).view()

It will output::

    p
    q
    1
    2
    3
    a
    b
    c

.. _operator-spread:

spread
---------

.. warning:: This operator is deprecated. See `combine`_ instead.

The ``spread`` operator combines the items emitted by the source channel with all the values in an array
or a ``Collection`` object specified as the operator argument. For example::

    Channel
        .from(1,2,3)
        .spread(['a','b'])
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

::

    [1, 'a']
    [1, 'b']
    [2, 'a']
    [2, 'b']
    [3, 'a']
    [3, 'b']
    Done




Forking operators
=================

The forking operators are:

* `branch`_
* `choice`_
* `multiMap`_
* `into`_
* `separate`_
* `tap`_

.. _operator-branch:

branch
------

.. note:: Requires Nextflow version ``19.08.0-edge`` or later.

The ``branch`` operator allows you to forward the items emitted by a source channel to one
or more output channels, `choosing` one out of them at a time.

The selection criteria is defined by specifying a :ref:`closure <script-closure>` that provides
one or more boolean expression, each of which is identified by a unique label. On the first expression 
that evaluates to a *true* value, the current item is bound to a named channel as the label identifier.
For example::

    Channel
        .from(1,2,3,40,50)
        .branch {
            small: it < 10
            large: it > 10
        }
        .set { result }

     result.small.view { "$it is small" }
     result.large.view { "$it is large" }

It shows::

    1 is small
    2 is small
    3 is small
    40 is large
    50 is large

.. note:: The above *small* and *large* strings maybe be printed interleaving each other
  due to the asynchronous execution of the ``view`` operator.

.. tip:: A default fallback condition can be specified using ``true`` as last branch condition. See the example below.

::

    Channel
        .from(1,2,3,40,50)
        .branch {
            small: it < 10
            large: it < 50
            other: true
        }


The value returned by each branch condition can be customised by specifying an optional expression statement(s)
just after the condition expression. For example::

       Channel
        .from(1,2,3,40,50)
        .branch {
            foo: it < 10
                return it+2

            bar: it < 50
                return it-2

            other: true
                return 0
        }


.. tip:: When the ``return`` keyword is omitted the value of the last expression statement is
  implicitly returned.

.. warning:: The branch evaluation closure must be specified inline, ie. it *cannot* be assigned to a
  variable and passed as argument to the operator, the way it can be done with other operators.

To create a branch criteria as variable that can be passed as an argument to more than one
``branch`` operator use the ``branchCriteria`` built-in method as shown below::

    def criteria = branchCriteria {
                    small: it < 10
                    large: it > 10
                    }

    Channel.from(1,2,30).branch(criteria).set { ch1 }
    Channel.from(10,20,1).branch(criteria).set { ch2 }


.. _operator-choice:

choice
------

.. warning:: The choice operator has been deprecated. Use `branch`_ instead.

The ``choice`` operator allows you to forward the items emitted by a source channel to two 
(or more) output channels, `choosing` one out of them at a time. 

The destination channel is selected by using a :ref:`closure <script-closure>` that must return the `index` number of the channel
where the item has to be sent. The first channel is identified by the index ``0``, the second as ``1`` and so on. 

The following example sends all string items beginning with ``Hello`` into ``queue1``, 
the others into ``queue2``  

::
  
    source = Channel.from 'Hello world', 'Hola', 'Hello John'
    queue1 = Channel.create()
    queue2 = Channel.create()

    source.choice( queue1, queue2 ) { a -> a =~ /^Hello.*/ ? 0 : 1 }

    queue1.view()

See also `branch`_ operator.

 .. _operator-multimap:

multiMap
--------

.. note:: Requires Nextflow version ``19.11.0-edge`` or later.

The multiMap operator allows you to forward the items emitted by a source channel to two
or more output channels mapping each input value as a separate element.

The mapping criteria is defined by specifying a :ref:`closure <script-closure>` that specify the
target channels labelled by a unique identifier followed by an expression statement that
evaluates the value to be assigned to such channel.

For example::

    Channel
        .from(1,2,3,4)
        .multiMap { it ->
            foo: it + 1
            bar: it * it
            }
        .set { result }

     result.foo.view { "foo $it" }
     result.bar.view { "bar $it" }

It prints::

    foo 2
    foo 3
    foo 4
    foo 5
    bar 1
    bar 4
    bar 9
    bar 16


.. tip:: The statement expression can be omitted when the value to be emitted is the same as
  the following one. If you need just need to forward the same value to multiple channels
  you can use the following the shorthand notation shown below.

::

   Channel
        .from(1,2,3)
        .multiMap { it -> foo: bar: it }
        .set { result }

As before this creates two channels but now both of them receive the same source items.


.. warning::
  The multi-map evaluation closure must be specified inline, ie. it *cannot* be assigned to a
  variable and passed as argument to the operator, the way it can be done with other operators.

To create a multi-map criteria as variable that can be passed as an argument to more than one
``multiMap`` operator use the ``multiMapCriteria`` built-in method as shown below::

    def criteria = multiMapCriteria {
                      small: it < 10
                      large: it > 10
                    }

    Channel.from(1,2,30).multiMap(criteria).set { ch1 }
    Channel.from(10,20,1).multiMap(criteria).set { ch2 }


.. _operator-into:

into
----

.. warning::
    The ``into`` operator is not available when using Nextflow DSL2 syntax.

The ``into`` operator connects a source channel to two or more target channels in such a way the values emitted by
the source channel are copied to the target channels. For example::

   Channel
        .from( 'a', 'b', 'c' )
        .into{ foo; bar }

    foo.view{ "Foo emit: " + it }
    bar.view{ "Bar emit: " + it }

::

    Foo emit: a
    Foo emit: b
    Foo emit: c
    Bar emit: a
    Bar emit: b
    Bar emit: c

.. note:: Note the use in this example of curly brackets and the ``;`` as channel names separator. This is needed
  because the actual parameter of ``into`` is a :ref:`closure <script-closure>` which defines the target channels
  to which the source one is connected.

A second version of the ``into`` operator takes an integer `n` as an argument and returns
a list of `n` channels, each of which emits a copy of the items that were emitted by the
source channel. For example::


    (foo, bar) = Channel.from( 'a','b','c').into(2)
    foo.view{ "Foo emit: " + it }
    bar.view{ "Bar emit: " + it }


.. note:: The above example takes advantage of the :ref:`multiple assignment <script-multiple-assignment>` syntax
  in order to assign two variables at once using the list of channels returned by the ``into`` operator.

See also `tap`_ and `separate`_ operators.


tap
---

The ``tap`` operator combines the functions of `into`_ and `separate`_ operators in such a way that
it connects two channels, copying the values from the source into the `tapped` channel. At the same
time it splits the source channel into a newly created channel that is returned by the operator itself.

The ``tap`` can be useful in certain scenarios where you may be required to concatenate multiple operations,
as in the following example::


    log1 = Channel.create()
    log2 = Channel.create()

    Channel
        .of ( 'a', 'b', 'c' )
        .tap ( log1 )
        .map { it * 2 }
        .tap ( log2 )
        .map { it.toUpperCase() }
        .view { "Result: $it" }

    log1.view { "Log 1: $it" }
    log2.view { "Log 2: $it" }

::

    Result: AA
    Result: BB
    Result: CC

    Log 1: a
    Log 1: b
    Log 1: c

    Log 2: aa
    Log 2: bb
    Log 2: cc


The ``tap`` operator also allows the target channel to be specified by using a closure. The advantage of this syntax
is that you won't need to previously create the target channel, because it is created implicitly by the operator itself.

Using the closure syntax the above example can be rewritten as shown below::

    Channel
        .of ( 'a', 'b', 'c' )
        .tap { log1 }
        .map { it * 2 }
        .tap { log2 }
        .map { it.toUpperCase() }
        .view { "Result: $it" }

    log1.view { "Log 1: $it" }
    log2.view { "Log 2: $it" }

See also `into`_ and `separate`_ operators.

.. _operator-separate:

separate
--------

.. warning:: The ``separate`` operator has been deprecated. Use `multiMap`_ instead.

The ``separate`` operator lets you copy the items emitted by the source channel into multiple 
channels, which each of these can receive a `separate` version of the same item. 

The operator applies a `mapping function` of your choosing to every item emitted by the source channel.
This function must return a list of as many values as there are output channels. Each entry in the result 
list will be assigned to the output channel with the corresponding position index. For example:: 

    queue1 = Channel.create()
    queue2 = Channel.create()

    Channel
        .from ( 2,4,8 ) 
        .separate( queue1, queue2 ) { a -> [a+1, a*a] }

    queue1.view { "Channel 1: $it" }
    queue2.view { "Channel 2: $it" }
	
::

	Channel 1: 3
	Channel 2: 4
	Channel 1: 5
	Channel 2: 16
	Channel 2: 64
	Channel 1: 9


When the `mapping function` is omitted, the source channel must emit tuples of values. In this case the operator ``separate``
splits the tuple in such a way that the value `i-th` in a tuple is assigned to the target channel with the corresponding position index.
For example::


     alpha = Channel.create()
     delta = Channel.create()

     Channel
        .from([1,2], ['a','b'], ['p','q'])
        .separate( alpha, delta )

     alpha.view { "first : $it" }
     delta.view { "second: $it" }

It will output::

        first : 1
        first : a
        first : p
        second: 2
        second: b
        second: q

A second version of the ``separate`` operator takes an integer `n` as an argument and returns a list of `n` channels,
each of which gets a value from the corresponding element in the list returned by the closure as explained above.
For example::	

    source = Channel.from(1,2,3)
    (queue1, queue2, queue3) = source.separate(3) { a -> [a, a+1, a*a] }

    queue1.view { "Queue 1 > $it" }
    queue2.view { "Queue 2 > $it" }
    queue3.view { "Queue 3 > $it" }

The output will look like the following fragment::

    Queue 1 > 1
    Queue 1 > 2
    Queue 1 > 3
    Queue 2 > 2
    Queue 2 > 3
    Queue 2 > 4
    Queue 3 > 1
    Queue 3 > 4
    Queue 3 > 9


.. warning:: In the above example, please note that since the ``subscribe`` operator is asynchronous,
  the output of ``channel2`` and ``channel3`` can be printed before the content of ``channel1``.

.. note:: The above example takes advantage of the :ref:`multiple assignment <script-multiple-assignment>` syntax
  in order to assign two variables at once using the list of channels returned by the ``separate`` operator.



See also: `multiMap`_, `into`_, `choice`_ and `map`_ operators.


Maths operators
================

This section talks about operators that performs maths operations on channels.

The maths operators are:

* `count`_
* `countBy`_
* `min`_
* `max`_
* `sum`_
* `toInteger`_

.. _operator-count:

count
--------

The ``count`` operator creates a channel that emits a single item: a number that represents the total number of
items emitted by the source channel. For example:: 

        Channel
            .from(9,1,7,5)
            .count()
            .view()
        // -> 4


An optional parameter can be provided in order to select which items are to be counted. 
The selection criteria can be specified either as a :ref:`regular expression <script-regexp>`, 
a literal value, a Java class, or a `boolean predicate` that needs to be satisfied. For example::


        Channel
            .from(4,1,7,1,1)
            .count(1)
            .view()
         // -> 3

        Channel
            .from('a','c','c','q','b')
            .count ( ~/c/ )
            .view()
        // -> 2
        
        Channel
            .from('a','c','c','q','b')
            .count { it <= 'c' }
            .view()
        // -> 4


.. _operator-countby:

countBy
----------

The ``countBy`` operator creates a channel which emits an associative array (i.e. ``Map`` object) 
that counts the occurrences of the emitted items in the source channel having the same key. 
For example::

    Channel
        .from( 'x', 'y', 'x', 'x', 'z', 'y' )
        .countBy()
        .view()

::

    [x:3, y:2, z:1]


An optional grouping criteria can be specified by using a :ref:`closure <script-closure>` 
that associates each item with the grouping key. For example::


    Channel
        .from( 'hola', 'hello', 'ciao', 'bonjour', 'halo' )
        .countBy { it[0] }
        .view()


::

    [h:3, c:1, b:1]


.. _operator-min:

min
------

The ``min`` operator waits until the source channel completes, and then emits the item that has the lowest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .min()
        .view { "Min value is $it" }

::

  Min value is 2

An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that returns the value to be compared. The example below shows how to find the string 
item that has the minimum length:: 

    Channel
    	.from("hello","hi","hey")
    	.min { it.size() } 
    	.view()

::

	 "hi"

Alternatively it is possible to specify a comparator function i.e. a :ref:`closure <script-closure>`
taking two parameters that represent two emitted items to be compared. For example:: 


    Channel
    	.from("hello","hi","hey")
    	.min { a,b -> a.size() <=> b.size() } 
    	.view()


.. _operator-max:

max
------

The ``max`` operator waits until the source channel completes, and then emits the item that has the greatest value.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .max()
        .view { "Max value is $it" }

::

  Max value is 8


An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that returns the value to be compared. The example below shows how to find the string 
item that has the maximum length:: 

    Channel
    	.from("hello","hi","hey")
    	.max { it.size() } 
    	.view()

::

	 "hello"

Alternatively it is possible to specify a comparator function i.e. a :ref:`closure <script-closure>`
taking two parameters that represent two emitted items to be compared. For example:: 


    Channel
    	.from("hello","hi","hey")
    	.max { a,b -> a.size() <=> b.size() } 
    	.view()


.. _operator-sum:

sum
------

The ``sum`` operator creates a channel that emits the sum of all the items emitted by the channel itself.
For example::

    Channel
        .from( 8, 6, 2, 5 )
        .sum()
        .view { "The sum is $it" }

::

    The sum is 21


An optional :ref:`closure <script-closure>` parameter can be specified in order to provide 
a function that, given an item, returns the value to be summed. For example:: 

	Channel
		.from( 4, 1, 7, 5 )
		.sum { it * it } 
		.view { "Square: $it" }

::

	Square: 91



toInteger
---------

The ``toInteger`` operator allows you to convert the string values emitted by a channel to ``Integer`` values. For
example::

    Channel
        .from( '1', '7', '12' )
        .toInteger()
        .sum()
        .view()



Other operators
========================

* `close`_
* `dump`_
* `ifEmpty`_
* `print`_
* `println`_
* `set`_
* `view`_

.. _operator-dump:

dump
----

The ``dump`` operator prints the items emitted by the channel to which is applied only when the option
``-dump-channels`` is specified on the ``run`` command line, otherwise it is ignored.

This is useful to enable the debugging of one or more channel content on-demand by using a command line option
instead of modifying your script code.

An optional ``tag`` parameter allows you to select which channel to dump. For example::

    Channel
        .from(1,2,3)
        .map { it+1 }
        .dump(tag:'foo')

    Channel
        .from(1,2,3)
        .map { it^2 }
        .dump(tag: 'bar')


Then you will be able to specify the tag ``foo`` or ``bar`` as an argument of the ``-dump-channels`` option to print
either the content of the first or the second channel. Multiple tag names can be specified separating them with a ``,``
character.


.. _operator-set:

set
----

The ``set`` operator assigns the channel to a variable whose name is specified as a closure parameter.
For example::

    Channel.from(10,20,30).set { my_channel }

This is semantically equivalent to the following assignment::

    my_channel = Channel.from(10,20,30)

However the ``set`` operator is more idiomatic in Nextflow scripting, since it can be used at the end
of a chain of operator transformations, thus resulting in a more fluent and readable operation.

.. _operator-ifempty:

ifEmpty
--------

The ``ifEmpty`` operator creates a channel which emits a default value, specified as the operator parameter, when the channel to which
is applied is *empty* i.e. doesn't emit any value. Otherwise it will emit the same sequence of entries as the original channel.

Thus, the following example prints::

    Channel .from(1,2,3) .ifEmpty('Hello') .view()

    1
    2
    3





Instead, this one prints::

    Channel.empty().ifEmpty('Hello') .view()

    Hello

The ``ifEmpty`` value parameter can be defined with a :ref:`closure <script-closure>`. In this case the result value of the closure evaluation
will be emitted when the empty condition is satisfied.

See also: :ref:`channel-empty` method.

.. _operator-print:

print
------

.. warning::
  The ``print`` operator is deprecated and not supported anymore when using DSL2 syntax. Use
  `view`_ instead.

The ``print`` operator prints the items emitted by a channel to the standard output.
An optional :ref:`closure <script-closure>` parameter can be specified to customise how items are printed.
For example::

  Channel
        .from('foo', 'bar', 'baz', 'qux')
        .print { it.toUpperCase() + ' ' }

It prints::

    FOO BAR BAZ QUX

See also: `println`_ and `view`_.

.. _operator-println:

println
--------

.. warning::
  The ``println`` operator is deprecated and not supported anymore when using DSL2 syntax. Use
  `view`_ instead.

The ``println`` operator prints the items emitted by a channel to the console standard output appending
a *new line* character to each of them. For example::

  Channel
        .from('foo', 'bar', 'baz', 'qux')
        .println()

It prints::

        foo
        bar
        baz
        qux


An optional closure parameter can be specified to customise how items are printed. For example::

  Channel
        .from('foo', 'bar', 'baz', 'qux')
        .view { "~ $it" }


It prints::

        ~ foo
        ~ bar
        ~ baz
        ~ qux

See also: `print`_ and `view`_.

.. _operator-view:

view
------

The ``view`` operator prints the items emitted by a channel to the console standard output. For example::

    Channel.from(1,2,3).view()

    1
    2
    3

Each item is printed on a separate line unless otherwise specified by using the ``newLine: false`` optional parameter.

How the channel items are printed can be controlled by using an optional closure parameter. The closure must return
the actual value of the item to be printed::

    Channel.from(1,2,3)
            .map { it -> [it, it*it] }
            .view { num, sqr -> "Square of: $num is $sqr" }

It prints::

    Square of: 1 is 1
    Square of: 2 is 4
    Square of: 3 is 9


.. note:: Both the *view* and `print`_ (or `println`_) operators consume the items emitted by the source channel to which they
    are applied. The main difference between them is that the former returns a newly created channel whose content
    is identical to the source channel while the latter does not. This allows the *view* operator to be chained like other operators.

.. _operator-close:

close
------

The ``close`` operator sends a termination signal over the channel, causing downstream processes or operators to stop.
In a common usage scenario channels are closed automatically by Nextflow, so you won't need to use this operator explicitly.

See also: :ref:`channel-empty` factory method.
.. Nextflow documentation master file, created by
   sphinx-quickstart on Sat May  5 16:57:28 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Nextflow's documentation!
=========================

Contents:

.. toctree::
   :maxdepth: 2

   getstarted
   basic
   script
   process
   channel
   operator
   executor
   config
   dsl2
   cli
   awscloud
   amazons3
   azure
   google
   conda
   docker
   shifter
   singularity
   podman
   charliecloud
   ignite
   kubernetes
   tracing
   sharing
   metadata
   mail
   plugins
   secrets
.. _getstart-page:

*******************
Get started
*******************

.. _getstart-requirement:

Requirements
============

`Nextflow` can be used on any POSIX compatible system (Linux, OS X, etc).
It requires Bash 3.2 (or later) and `Java 8 (or later, up to 17) <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_ to be installed.

For the execution in a cluster of computers the use a shared file system is required to allow
the sharing of tasks input/output files.

Windows system is supported through `WSL <https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_.

.. _getstart-install:

Installation
============

`Nextflow` is distributed as a self-installing package, which means that it does not require any special installation procedure.

It only needs two easy steps:

#.  Download the executable package by copying and pasting either one of the following commands in your terminal
    window: ``wget -qO- https://get.nextflow.io | bash``

    Or, if you prefer `curl`: ``curl -s https://get.nextflow.io | bash``

    This will create the ``nextflow`` main executable file in the current directory.

#.  Make the binary executable on your system by running ``chmod +x nextflow``.

#.  Optionally, move the ``nextflow`` file to a directory accessible by your ``$PATH`` variable
    (this is only required to avoid remembering and typing the full path to ``nextflow`` each time you need to run it).

.. tip:: Set ``export CAPSULE_LOG=none`` to make the dependency installation logs less verbose.

.. tip::
    If you don't have ``curl`` nor ``wget``, you can also download the Nextflow launcher script from the
    `project releases page <https://github.com/nextflow-io/nextflow/releases/latest>`_ on GitHub . Once downloaded
    continue the installation from the point 2 in the above guide.

.. note::
    To avoid downloading the dependencies, you could also use the ``nextflow-VERSION-all`` distribution available from Github for every `Nextflow` release version.
   #. Go to the `Github releases page <https://github.com/nextflow-io/nextflow/releases>`__ and unfold the `Assets` section for a release.
   #. Copy the URL of the ``nextflow-VERSION-all`` asset and issue the download command on your terminal. ``wget -qO- ASSET-URL``. It will create the completely self-contained ``nextflow-VERSION-all`` executable file in the current directory.

Updates
=======

Having Nextflow installed in your computer you can update to the latest version using the following command::

    nextflow self-update


.. tip::
  You can temporary switch or stick to a specific version of Nextflow just prefixing the ``nextflow`` command
  with the ``NXF_VER`` environment variable. For example::

    NXF_VER=20.04.0 nextflow run hello

Stable & Edge releases
======================

A *stable* version of Nextflow is released on a six-months basic schedule, in the 1st and 3rd quarter of every year.

Along with the stable release, an `edge` version is released on a monthly basis. This version is useful to test and
use most recent updates and experimental features.

To use the latest `edge` release run the following snippet in your shell terminal::

    export NXF_EDGE=1
    nextflow self-update


.. _getstart-first:

Your first script
==================

Copy the following example into your favourite text editor and save it to a file named ``tutorial.nf`` ::

    #!/usr/bin/env nextflow

    params.str = 'Hello world!'

    process splitLetters {

        output:
        file 'chunk_*' into letters

        """
        printf '${params.str}' | split -b 6 - chunk_
        """
    }


    process convertToUpper {

        input:
        file x from letters.flatten()

        output:
        stdout result

        """
        cat $x | tr '[a-z]' '[A-Z]'
        """
    }

    result.view { it.trim() }


This script defines two processes. The first splits a string into 6-character chunks, writing each one to a file with the prefix ``chunk_``,
and the second receives these files and transforms their contents to uppercase letters.
The resulting strings are emitted on the ``result`` channel and the final output is printed by the
``view`` operator.



Execute the script by entering the following command in your terminal::

   nextflow run tutorial.nf

It will output something similar to the text shown below::

    N E X T F L O W  ~  version 19.04.0
    executor >  local (3)
    [69/c8ea4a] process > splitLetters   [100%] 1 of 1 ✔
    [84/c8b7f1] process > convertToUpper [100%] 2 of 2 ✔
    HELLO
    WORLD!


You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It's worth noting that the process ``convertToUpper`` is executed in parallel, so there's no guarantee that the instance
processing the first split (the chunk `Hello`) will be executed before the one processing the second split (the chunk `world!`).

Thus, it is perfectly possible that you will get the final result printed out in a different order::

    WORLD!
    HELLO



.. tip:: The hexadecimal numbers, like ``22/7548fa``, identify the unique process execution. These numbers are
  also the prefix of the directories where each process is executed. You can inspect the files produced by them
  changing to the directory ``$PWD/work`` and using these numbers to find the process-specific execution path.

.. _getstart-resume:

Modify and resume
-----------------

`Nextflow` keeps track of all the processes executed in your pipeline. If you modify some parts of your script,
only the processes that are actually changed will be re-executed. The execution of the processes that are not changed
will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the ``convertToUpper`` process in the previous example, replacing the
process script with the string ``rev $x``, so that the process looks like this::

    process convertToUpper {

        input:
        file x from letters

        output:
        stdout result

        """
        rev $x
        """
    }

Then save the file with the same name, and execute it by adding the ``-resume`` option to the command line::

    nextflow run tutorial.nf -resume


It will print output similar to this::

    N E X T F L O W  ~  version 19.04.0
    executor >  local (2)
    [69/c8ea4a] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
    [d0/e94f07] process > convertToUpper [100%] 2 of 2 ✔
    olleH
    !dlrow


You will see that the execution of the process ``splitLetters`` is actually skipped (the process ID is the same), and
its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.


.. tip:: The pipeline results are cached by default in the directory ``$PWD/work``. Depending on your script, this folder
  can take of lot of disk space. If you are sure you won't resume your pipeline execution, clean this folder periodically.

.. _getstart-params:

Pipeline parameters
--------------------

Pipeline parameters are simply declared by prepending to a variable name the prefix ``params``, separated by dot character.
Their value can be specified on the command line by prefixing the parameter name with a double `dash` character, i.e. ``--paramName``

For the sake of this tutorial, you can try to execute the previous example specifying a different input
string parameter, as shown below::

  nextflow run tutorial.nf --str 'Bonjour le monde'


The string specified on the command line will override the default value of the parameter. The output
will look like this::

    N E X T F L O W  ~  version 19.04.0
    executor >  local (4)
    [8b/16e7d7] process > splitLetters   [100%] 1 of 1 ✔
    [eb/729772] process > convertToUpper [100%] 3 of 3 ✔
    m el r
    edno
    uojnoB


.. tip::
    As of version 20.11.0-edge any ``.`` (dot) character in a parameter name is interpreted as the delimiter
    or nested scope e.g. ``--foo.bar Hello`` will be accessible from the script as `params.foo.bar`.
    If you want to have a parameter name including a ``.`` (dot) character escape it using the back-slash character e.g.
    ``--foo\.bar Hello``
.. _docker-page:

******************
Docker containers
******************

Nextflow integration with `Docker containers <http://www.docker.io>`_ technology allows you to write self-contained
and truly reproducible computational pipelines.

By using this feature any process in a Nextflow script can be transparently executed into a Docker container. This may
be extremely useful to package the binary dependencies of a script into a standard and portable format that can be 
executed on any platform supporting the Docker engine.

Prerequisites
==============

You will need Docker installed on your execution environment e.g. your computer or a distributed cluster, depending
on where you want to run your pipeline.

If you are running Docker on Mac OSX make sure you are mounting your local ``/Users`` directory into the Docker VM as
explained in this excellent tutorial: `How to use Docker on OSX <http://viget.com/extend/how-to-use-docker-on-os-x-the-missing-guide>`_.


How it works
=============

You won't need to modify your Nextflow script in order to run it with Docker. Simply specify the Docker image from
where the containers are started by using the ``-with-docker`` command line option. For example::

  nextflow run <your script> -with-docker [docker image]

Every time your script launches a process execution, Nextflow will run it into a Docker container created by using the
specified image. In practice Nextflow will automatically wrap your processes and run them by executing the ``docker run``
command with the image you have provided.

.. note:: A Docker image can contain any tool or piece of software you may need to carry out a process execution. Moreover the
  container is run in such a way that the process result files are created in the hosting file system, thus
  it behaves in a completely transparent manner without requiring extra steps or affecting the flow in your pipeline.

If you want to avoid entering the Docker image as a command line parameter, you can define it in the Nextflow configuration
file. For example you can add the following lines in the ``nextflow.config`` file::

    process.container = 'nextflow/examples:latest'
    docker.enabled = true

In the above example replace ``nextflow/examples:latest`` with any Docker image of your choice.

Read the :ref:`config-page` page to learn more about the ``nextflow.config`` file and how to use it to configure
your pipeline execution.

.. warning::
    Nextflow automatically manages the file system mounts each time a container is launched depending on the process
    input files. Note, however, that when a process input is a *symbolic link* file, the linked file **must** be stored
    in the same folder where the symlink is located, or any its sub-folder. Otherwise the process execution will fail because the
    launched container won't be able to access the linked file.


Multiple containers
=====================

It is possible to specify a different Docker image for each process definition in your pipeline script. Let's
suppose you have two processes named ``foo`` and ``bar``. You can specify two different Docker images for them
in the Nextflow script as shown below::

    process foo {
      container 'image_name_1'

      '''
      do this
      '''
    }

    process bar {
      container 'image_name_2'

      '''
      do that
      '''
    }


Alternatively, the same containers definitions can be provided by using the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    docker {
        enabled = true
    }


Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.


Advanced settings 
==================

Docker advanced configuration settings are described in :ref:`config-docker` section in the Nextflow configuration page.













.. _amazons3-page:

*******************
Amazon S3 storage
*******************

Nextflow includes the support for Amazon S3 storage. Files stored in a S3 bucket can be accessed
transparently in your pipeline script like any other file in the local file system.

S3 path
---------
In order to access a S3 file you only need to prefix the file path with the ``s3`` schema and the `bucket` name
where it is stored.

For example if you need to access the file ``/data/sequences.fa`` stored in a bucket with name ``my-bucket``,
that file can be accessed using the following fully qualified path::

   s3://my-bucket/data/sequences.fa


The usual file operations can be applied on a path handle created using the above notation. For example the content
of a S3 file can be printed as shown below::

    println file('s3://my-bucket/data/sequences.fa').text


See section :ref:`script-file-io` to learn more about available file operations.




Security credentials
---------------------

Amazon access credentials can be provided in two ways:

#. Using AWS access and secret keys in your pipeline configuration.
#. Using IAM roles to grant access to S3 storage on Amazon EC2 instances.

AWS access and secret keys
===========================

The AWS access and secret keys can be specified by using the ``aws`` section in the ``nextflow.config`` configuration
file as shown below::

  aws {
    accessKey = '<Your AWS access key>'
    secretKey = '<Your AWS secret key>'
    region = '<AWS region identifier>'
  }


If the access credentials are not found in the above file, Nextflow looks for AWS credentials in a number of different
places, including environment variables and local AWS configuration files.


Nextflow looks for AWS credentials in the following order:

    #. the ``nextflow.config`` file in the pipeline execution directory
    #. the environment variables ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY``
    #. the environment variables ``AWS_ACCESS_KEY`` and ``AWS_SECRET_KEY``
    #. the `default` profile in the AWS credentials file located at ``~/.aws/credentials``
    #. the `default` profile in the AWS client configuration file located at ``~/.aws/config``
    #. the temporary AWS credentials provided by an IAM instance role. See `IAM Roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_ documentation for details.


More information regarding `AWS Security Credentials <http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html>`_
are available in Amazon documentation.

IAM roles Amazon EC2 instances
================================

When running your pipeline into a Ec2 instance, IAM roles can be used to grant access to AWS resources.

In this scenario, you only need to launch the Ec2 instance specifying a IAM role which includes a
`S3 full access` policy. Nextflow will detected and acquire automatically the access grant to the S3 storage,
without any further configuration.

Learn more about `Using IAM Roles to Delegate Permissions to Applications that Run on Amazon EC2 <http://docs.aws.amazon.com/IAM/latest/UserGuide/roles-usingrole-ec2instance.html>`_ on Amazon
documentation.

China regions
-------------

To use an AWS China region, please make sure to specify the corresponding AWS API S3 endpoint in the Nextflow configuration
file as shown below::

    aws.client.endpoint = "https://s3.cn-north-1.amazonaws.com.cn"

Read more about AWS API endpoints in the `AWS documentation <https://docs.aws.amazon.com/general/latest/gr/s3.html>`_

Advanced configuration
-----------------------

Read :ref:`AWS configuration<config-aws>` section to learn more about advanced S3 client configuration options.







.. _shifter-page:

******************
Shifter Containers
******************

`Shifter <https://docs.nersc.gov/programming/shifter/overview/>`__ is container engine alternative to
`Docker <https://www.docker.com>`__. Shifter works by converting Docker images to a common format that can then be
distributed and launched on HPC systems. The user interface to Shifter enables a user to select an image
from the `dockerhub registry <https://hub.docker.com/>`__ and then submit jobs which run entirely within the container.

Nextflow provides built-in support for Shifter. This allows you to control the execution environment of the processes
in your pipeline by running them in isolated containers along-side their dependencies.

Moreover the support provided by Nextflow for different container technologies, allows the same pipeline to be
transparently executed both with :ref:`Docker <_docker-page>`, :ref:`Singularity <_singularity-page>` or
:ref:`Shifter <_shifter-page>` containers, depending on the available engine in target execution platforms.

Prerequisites
=============

You need Shifter and Shifter image gateway installed in your execution environment, i.e: your personal computed or the
entry node of a distributed cluster. In the case of the distributed cluster case, you should have Shifter installed on
all of the compute nodes and the ``shifterimg`` command should also be available and Shifter properly setup to access the
Image gateway, for more information see the
`official documentation <https://github.com/NERSC/shifter/tree/master/doc>`__.

.. note:: This feature requires Shifter version 18.03 (or later) and Nextflow 19.10.0 (or later).

Images
======

Shifter converts a docker image to squashfs layers which are distributed and launched in the cluster. For more info on
how to Build Shifter images see the
`official documentation <https://docs.nersc.gov/programming/shifter/how-to-use/#building-shifter-images>`__.

How it works
============

The integration for Shifter, at this time, requires you to set up the following parameters in your config file::

  process.container = "dockerhub_user/image_name:image_tag"
  shifter.enabled = true

and it will always try to search the Docker Hub registry for the images.

.. note:: if you do not specify an image tag it will fetch the 'latest' tag by default.

Multiple Containers
===================

It is possible to specify a different Shifter image for each process definition in your pipeline script. For example,
let's suppose you have two processes named ``foo`` and ``bar``. You can specify two different Shifter images
specifying them in the ``nextflow.config`` file as shown below::

    process {
        withName:foo {
            container = 'image_name_1'
        }
        withName:bar {
            container = 'image_name_2'
        }
    }
    shifter {
        enabled = true
    }

Read the :ref:`Process scope <config-process>` section to learn more about processes configuration.
***************
Basic concepts
***************


`Nextflow` is a reactive workflow framework and a programming `DSL <http://en.wikipedia.org/wiki/Domain-specific_language>`_
that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many
simple but powerful command-line and scripting tools that, when chained together, facilitate complex
data manipulations.

`Nextflow` extends this approach, adding the ability to define complex program interactions and a high-level
parallel computational environment based on the `dataflow` programming model.


Processes and channels
----------------------

In practice a Nextflow pipeline script is made by joining together different processes.
Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state.
The only way they can communicate is via asynchronous FIFO queues, called `channels` in Nextflow.

Any process can define one or more channels as `input` and `output`. The interaction between these processes,
and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.

A Nextflow script looks like this::

    // Script parameters
    params.query = "/some/data/sample.fa"
    params.db = "/some/path/pdb"

    db = file(params.db)
    query_ch = Channel.fromPath(params.query)

    process blastSearch {
        input:
        file query from query_ch

        output:
        file "top_hits.txt" into top_hits_ch

        """
        blastp -db $db -query $query -outfmt 6 > blast_result
        cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
        """
    }

    process extractTopHits {
        input:
        file top_hits from top_hits_ch

        output:
        file "sequences.txt" into sequences_ch

        """
        blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
        """
    }



The above example defines two processes. Their execution order is not determined by the fact that the ``blastSearch``
process comes before ``extractTopHits`` in the script (it could also be written the other way around).

Instead, because the first process defines the channel ``top_hits_ch`` in its output declarations, and the
process ``extractTopHits`` defines the channel in its input declaration, a communication link is established.

This linking via the channels means that `extractTopHits` is waiting for the output of `blastSearch`, and then
runs `reactively` when the channel has contents.

.. TODO describe that both processes are launched at the same time

Read the :ref:`Channel <channel-page>` and :ref:`Process <process-page>` sections to learn more about these features.


Execution abstraction
---------------------

While a process defines `what` command or script has to be executed, the `executor` determines `how`
that script is actually run on the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline
development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, `Nextflow` provides an abstraction between the pipeline's functional logic and the underlying execution system.
Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud,
without modifying it, by simply defining the target execution platform in the configuration file.

The following batch schedulers are supported:

* `Open grid engine <http://gridscheduler.sourceforge.net/>`_
* `Univa grid engine <http://www.univa.com/>`_
* `Platform LSF <http://www.ibm.com/systems/technicalcomputing/platformcomputing/products/lsf/>`_
* `Linux SLURM <https://computing.llnl.gov/linux/slurm/>`_
* `PBS Works <http://www.pbsworks.com/gridengine/>`_
* `Torque <http://www.adaptivecomputing.com/products/open-source/torque/>`_
* `HTCondor <https://research.cs.wisc.edu/htcondor/>`_


The following cloud platforms are supported:

* `Amazon Web Services (AWS) <https://aws.amazon.com/>`_
* `Google Cloud Platform (GCP) <https://cloud.google.com/>`_
* `Kubernetes <https://kubernetes.io/>`_

Read the :ref:`executor-page` to learn more about the Nextflow executors.


Scripting language
------------------

`Nextflow` is designed to have a minimal learning curve, without having to pick up
a new programming language. In most cases, users can utilise their current skills to develop
Nextflow workflows. However, it also provides a powerful scripting DSL.

Nextflow scripting is an extension of the `Groovy programming language <http://en.wikipedia.org/wiki/Groovy_(programming_language)>`_,
which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java
in that it simplifies the writing of code and is more approachable.

Read the :ref:`script-page` section to learn about the Nextflow scripting language.


.. TODO Running pipeline


.. TODO Pipeline parameters


Configuration options
---------------------

Pipeline configuration properties are defined in a file named ``nextflow.config`` in the pipeline execution directory. 

This file can be used to define which executor to use, the process's environment variables, pipeline parameters etc. 

A basic configuration file might look like this::

	process { 
	  executor='sge'
	  queue = 'cn-el6' 
	}


Read the :ref:`config-page` section to learn more about the Nextflow configuration file and settings.



.. _conda-page:

******************
Conda environments
******************

`Conda <https://conda.io/>`_ is an open source package and environment management
system that simplifies the installation and the configuration of complex software packages
in a platform agnostic manner.

Nextflow has built-in support for Conda that allows the configuration of workflow dependencies
using Conda recipes and environment files.

This allows Nextflow applications to use popular tool collections
such as `Bioconda <https://bioconda.github.io>`_ whilst taking advantage of the configuration
flexibility provided by Nextflow.

Prerequisites
-------------

This feature requires Nextflow version 0.30.x or higher and the Conda or
`Miniconda <https://conda.io/miniconda.html>`_ package manager installed on your system.

How it works
------------

Nextflow  automatically creates and activates the Conda environment(s) given the dependencies
specified by each process.

Dependencies are specified by using the :ref:`process-conda` directive, providing either
the names of the required Conda packages, the path of a Conda environment yaml file or
the path of an existing Conda environment directory.

.. note:: Conda environments are stored on the file system. By default Nextflow instructs Conda to save
  the required environments in the pipeline work directory. Therefore the same environment can be created/saved
  multiple times across multiple executions when using a different work directory.

You can specify the directory where the Conda environments are stored using the ``conda.cacheDir``
configuration property (see the :ref:`configuration page <config-conda>` for details).
When using a computing cluster, make sure to use a shared file system path
accessible from all computing nodes.

.. warning:: The Conda environment feature is not supported by executors which use
  a remote object storage as a work directory eg. AWS Batch.

Use Mamba to resolve packages
=============================

It is also possible to use `mamba <https://github.com/mamba-org/mamba>`_ for speeding up creation of conda environments. For more information on how to enable this feature please refer :ref:`Conda <config-conda>`.

.. warning:: The use of ``mamba`` to create conda environments is an experimental feature and the functionality may change in future releases.


Use Conda package names
=======================

Conda package names can specified using the ``conda`` directive. Multiple package names can be specified
by separating them with a blank space.
For example::

  process foo {
    conda 'bwa samtools multiqc'

    '''
    your_command --here
    '''
  }


Using the above definition a Conda environment that includes BWA, Samtools and MultiQC tools is created and
activated when the process is executed.

The usual Conda package syntax and naming conventions can be used. The version of a package can be
specified after the package name as shown here ``bwa=0.7.15``.

The name of the channel where a package is located can be specified prefixing the package with
the channel name as shown here ``bioconda::bwa=0.7.15``.


Use Conda environment files
===========================

Conda environments can also be defined using one or more Conda environment files. This is a file that
lists the required packages and channels structured using the YAML format. For example::

    name: my-env
    channels:
      - conda-forge
      - bioconda
      - defaults
    dependencies:
      - star=2.5.4a
      - bwa=0.7.15

Read the Conda documentation for more details about how to create `environment files <https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-file-manually>`_.

The path of an environment file can be specified using the ``conda`` directive::

  process foo {
    conda '/some/path/my-env.yaml'

    '''
    your_command --here
    '''
  }

.. warning:: The environment file name **must** end with a ``.yml`` or ``.yaml`` suffix otherwise 
  it won't be properly recognised.


Alternatively it is also possible to provide the dependencies using a plain text file,
just listing each package name as a separate line. For example::

      bioconda::star=2.5.4a
      bioconda::bwa=0.7.15
      bioconda::multiqc=1.4


.. note:: Like before the extension matter, make sure such file ends with the ``.txt`` extension.


Use existing Conda environments
===============================

If you already have a local Conda environment, you can use it in your workflow specifying the
installation directory of such environment by using the ``conda`` directive::


  process foo {
    conda '/path/to/an/existing/env/directory'

    '''
    your_command --here
    '''
  }


Best practices
--------------

When a ``conda`` directive is used in any ``process`` definition within the workflow script, Conda tool is required for
the workflow execution.

Specifying the Conda environments in a separate configuration :ref:`profile <config-profiles>` is therefore
recommended to allow the execution via a command line option and to enhance the workflow portability. For example::
  
  profiles {
    conda {
      process.conda = 'samtools'
    }

    docker {
      process.container = 'biocontainers/samtools'
      docker.enabled = true
    }
  }

The above configuration snippet allows the execution either with Conda or Docker specifying ``-profile conda`` or
``-profile docker`` when running the workflow script.


Advanced settings
-----------------


Conda advanced configuration settings are described in the :ref:`Conda <config-conda>` section on the Nextflow configuration page.
.. _executor-page:

***********
Executors
***********

In the Nextflow framework architecture, the `executor` is the component that determines the system where a pipeline
process is run and supervises its execution.

The `executor` provides an abstraction between the pipeline processes and the underlying execution system. This
allows you to write the pipeline functional logic independently from the actual processing platform.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager
or the cloud by simply changing the executor definition in the Nextflow configuration file.

.. _local-executor:

Local
=====

The `local` executor is used by default. It runs the pipeline processes in the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads` and by taking advantage of multi-cores
architecture provided by the CPU.

In a common usage scenario, the `local` executor can be useful to develop and test your pipeline script in your computer,
switching to a cluster facility when you need to run it on production data.


.. _sge-executor:

SGE
===

The `SGE` executor allows you to run your pipeline script by using a `Sun Grid Engine <http://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_
cluster or a compatible platform (`Open Grid Engine <http://gridscheduler.sourceforge.net/>`_, `Univa Grid Engine <http://www.univa.com/products/grid-engine.php>`_, etc).

Nextflow manages each process as a separate grid job that is submitted to the cluster by using the ``qsub`` command.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SGE executor simply set to ``process.executor`` property to ``sge`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-memory`
* :ref:`process-penv`
* :ref:`process-time`
* :ref:`process-clusterOptions`

.. _lsf-executor:

LSF
===

The `LSF` executor allows you to run your pipeline script by using a `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``bsub`` command.

Being so, the pipeline must be launched from a node where the ``bsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the LSF executor simply set to ``process.executor`` property to ``lsf`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. note::

    LSF supports both *per-core* and *per-job* memory limit. Nextflow assumes that LSF works in the
    *per-core* memory limits mode, thus it divides the requested :ref:`process-memory` by the number of requested :ref:`process-cpus`.

    This is not required when LSF is configured to work in *per-job* memory limit mode. You will need to specified that
    adding the option ``perJobMemLimit`` in :ref:`config-executor` in the Nextflow configuration file.

    See also the `Platform LSF documentation <https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita>`_.


.. _slurm-executor:

SLURM
=====


The `SLURM` executor allows you to run your pipeline script by using the `SLURM <https://slurm.schedmd.com/documentation.html>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``sbatch`` command.

Being so, the pipeline must be launched from a node where the ``sbatch`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SLURM executor simply set to ``process.executor`` property to ``slurm`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. note:: SLURM `partitions` can be considered jobs queues. Nextflow allows you to set partitions by using the above ``queue``
    directive.

.. tip:: Nextflow does not provide a direct support for SLURM multi-clusters feature. If you need to
  submit workflow executions to a cluster that is not the current one, specify it setting the
  ``SLURM_CLUSTERS`` variable in the launching environment. 

.. _pbs-executor:

PBS/Torque
==========

The `PBS` executor allows you to run your pipeline script by using a resource manager belonging to the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ family of batch schedulers.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the PBS executor simply set the property ``process.executor = 'pbs'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _pbspro-executor:

PBS Pro
=======

The `PBS Pro` executor allows you to run your pipeline script by using the `PBS Pro <https://www.pbspro.org/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the PBS Pro executor simply set the property ``process.executor = 'pbspro'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _moab-executor:

Moab
====

The `Moab` executor allows you to run your pipeline script by using the
`Moab <https://en.wikipedia.org/wiki/Moab_Cluster_Suite>`_ resource manager by
`Adaptive Computing <http://www.adaptivecomputing.com/>`_.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``msub`` command provided
by the resource manager.

Being so, the pipeline must be launched from a node where the ``msub`` command is available, that is, in a common usage
scenario, the compute cluster `login` node.

To enable the `Moab` executor simply set the property ``process.executor = 'moab'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _nqsii-executor:

NQSII
=====

The `NQSII` executor allows you to run your pipeline script by using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the NQSII executor simply set the property ``process.executor = 'nqsii'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _condor-executor:

HTCondor
========

The `HTCondor` executor allows you to run your pipeline script by using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ resource manager.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``condor_submit`` command.

Being so, the pipeline must be launched from a node where the ``condor_submit`` command is available, that is, in a
common usage scenario, the cluster `head` node.

.. note::
  The HTCondor executor for Nextflow does not support at this time the HTCondor ability to transfer input/output data to
  the corresponding job computing node. Therefore the data needs to be made accessible to the computing nodes using
  a shared file system directory from where the Nextflow workflow has to be executed (or specified via the ``-w`` option).

To enable the HTCondor executor simply set to ``process.executor`` property to ``condor`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-disk`
* :ref:`process-clusterOptions`


.. _ignite-executor:

Ignite
======

The `Ignite` executor allows you to run a pipeline by using the `Apache Ignite <https://ignite.apache.org/>`_ clustering
technology that is embedded with the Nextflow runtime.

To enable this executor set the property ``process.executor = 'ignite'`` in the ``nextflow.config`` file.

The amount of resources requested by each task submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`

Read the :ref:`ignite-page` section in this documentation to learn how to configure Nextflow to deploy and run an
Ignite cluster in your infrastructure.

.. _k8s-executor:

Kubernetes
==========

Nextflow provides built-in support for `Kubernetes <http://kubernetes.io/>`_ clustering technology. It allows
you to deploy and transparently run a Nextflow pipeline in a Kubernetes cluster.

The following directives can be used to define the amount of computing resources needed and the container(s) to use:

* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-container`

See the :ref:`Kubernetes documentation <k8s-page>` to learn how to deploy a workflow execution in a Kubernetes cluster.

.. _awsbatch-executor:

AWS Batch
==========

Nextflow supports `AWS Batch <https://aws.amazon.com/batch/>`_ service which allows submitting jobs in the cloud
without having to spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks,
which makes deploying pipelines much simpler.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'awsbatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads. Moreover a S3 bucket must be used as pipeline work directory.

See the :ref:`AWS Batch<awscloud-batch>` page for further configuration details.

.. _azurebatch-executor:

Azure Batch
============

Nextflow supports `Azure Batch <https://azure.microsoft.com/en-us/services/batch/>`_ service which allows submitting jobs in the cloud
without having to spin out and manage a cluster of virtual machines. Azure Batch uses Docker containers to run tasks,
which makes deploying pipelines much simpler.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'azurebatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a cloud virtual machine. The latter is suggested for heavy or long
running workloads. Moreover a Azure Blob storage container must be used as pipeline work directory.

See the :ref:`Azure Batch <azure-batch>` page for further configuration details.

.. _google-lifesciences-executor:

Google Life Sciences
====================

`Google Cloud Life Sciences <https://cloud.google.com/life-sciences>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Life Sciences API which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as pipelines (it requires Nextflow 20.01.0 or later).

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor set the property ``process.executor = 'google-lifesciences'`` in the ``nextflow.config`` file.

See the :ref:`Google Life Sciences <google-lifesciences>` page for further configuration details.

.. _ga4ghtes-executor:

GA4GH TES
=========

.. warning:: This is an experimental feature and it may change in a future release. It requires Nextflow
  version 0.31.0 or later.

The `Task Execution Schema <https://github.com/ga4gh/task-execution-schemas>`_ (TES) project
by the `GA4GH <https://www.ga4gh.org>`_ standardisation initiative is an effort to define a
standardized schema and API for describing batch execution tasks in portable manner.

Nextflow includes an experimental support for the TES API providing a ``tes`` executor which allows
the submission of workflow tasks to a remote execution back-end exposing a TES API endpoint.

To use this feature define the following variables in the workflow launching environment::

    export NXF_MODE=ga4gh
    export NXF_EXECUTOR=tes
    export NXF_EXECUTOR_TES_ENDPOINT='http://back.end.com'
    

It is important that the endpoint is specified without the trailing slash; otherwise, the resulting URLs will be not
normalized and the requests to TES will fail.

Then you will be able to run your workflow over TES using the usual Nextflow command line. Be sure to specify the Docker
image to use, i.e.::

    nextflow run rnaseq-nf -with-docker alpine

.. note:: If the variable ``NXF_EXECUTOR_TES_ENDPOINT`` is omitted the default endpoint is ``http://localhost:8000``.

.. tip:: You can use a local `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ server using the following launch
  command line::

  ./funnel server --Server.HTTPPort 8000 --LocalStorage.AllowedDirs $HOME run

  (tested with version 0.8.0 on macOS)

.. warning:: Make sure the TES back-end can access the workflow work directory when
  data is exchanged using a local or shared file system.


**Known limitation**

* Automatic deployment of workflow scripts in the `bin` folder is not supported.
* Process output directories are not supported. For details see `#76 <https://github.com/ga4gh/task-execution-schemas/issues/76>`_.
* Glob patterns in process output declarations are not supported. For details see `#77 <https://github.com/ga4gh/task-execution-schemas/issues/77>`_.

.. _oar-executor:

OAR
===

The `OAR` executor allows you to run your pipeline script by using the `OAR <https://oar.imag.fr>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``oarsub`` command.

Being so, the pipeline must be launched from a node where the ``oarsub`` command is available, that is, in a common usage scenario, the cluster `head` node.

To enable the `OAR` executor simply set to ``process.executor`` property to ``oar`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

**Known limitation**

* `clusterOptions` should be given, if more than one, semicolon separated. It ensures the `OAR` batch script to be accurately formatted.

```
clusterOptions = '-t besteffort;--project myproject'
```
