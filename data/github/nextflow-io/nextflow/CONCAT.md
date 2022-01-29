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
