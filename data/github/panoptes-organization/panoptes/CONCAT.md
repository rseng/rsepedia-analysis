# ![alt text](panoptes/static/src/img/brand/panoptes.png "panoptes")


Bioinformaticians and data scientists, rely on computational frameworks (e.g. [snakemake](https://snakemake.readthedocs.io/en/stable/), [nextflow](https://www.nextflow.io/), [CWL](https://www.commonwl.org/), [WDL](https://software.broadinstitute.org/wdl/)) to process, analyze and integrate data of various types. Such frameworks allow scientists to combine software and custom tools of different origin in a unified way, which lets them reproduce the results of others, or reuse the same pipeline on different datasets. One of the fundamental issues is that the majority of the users execute multiple pipelines at the same time, or execute a multistep pipeline for a big number of datasets, or both, making it hard to track the execution of the individual steps or monitor which of the processed datasets are complete. panoptes is a tool that monitors the execution of such workflows.

panoptes is a service that can be used by:
- Data scientists, bioinformaticians, etc. that want to have a general overview of the progress of their pipelines and the status of their jobs
- Administrations that want to monitor their servers
- Web developers that want to integrate the service in bigger web applications

**Note:** panoptes is in early development stage and the first proof of concept server will support only workflows written in [snakemake](https://snakemake.readthedocs.io/en/stable/).

# Installation

## Basic installation process

### Requirements

- Python>=3.6
- virtualenv
- [sqlite3](https://www.sqlite.org/download.html)

### Option 1: Install via pypi and run server

Create virtual environment
```bash
virtualenv -p `which python3` venv
```

Activate virtual environment
```bash
source venv/bin/activate
```

Install via pypi
```bash
pip install panoptes-ui
```
Run server
```bash
panoptes
```
Server should run on: 127.0.0.1:5000

By default it should generate an sqlite database: .panoptes.db

### Option 2: Install via conda and run server

Create conda environment
```bash
conda create --name panoptes
```

Activate conda environment
```bash
conda activate panoptes
```

Install via pypi OR conda
```bash
conda install -c panoptes-organization panoptes-ui
```
Run server
```bash
panoptes
```
Server should run on: 127.0.0.1:5000

By default it should generate an sqlite database: .panoptes.db

### Option 3: Install from source code and run server

Clone repo
```bash
git clone https://github.com/panoptes-organization/panoptes.git
```

Enter repo
```bash
cd panoptes
```

Create virtual environment
```bash
virtualenv -p `which python3` venv
```

Activate virtual environment
```bash
source venv/bin/activate
```

Install all requirements
```bash
pip install .
```

Run server
```bash
panoptes
```
Server should run on: 127.0.0.1:5000

By default it should generate an sqlite database: .panoptes.db 

## Docker installation

### Requirements

- docker
- docker-compose

### Build and run with docker-compose

Build
```bash
docker-compose build
```

Run
```bash
docker-compose up -d
```

Server should run on: http://127.0.0.1:8000

Stop
```bash
docker-compose down
```

### Run an example workflow

In order to run an example workflow please follow the instructions [here](https://github.com/panoptes-organization/snakemake_example_workflow)

### panoptes in action

[![Watch the video](https://img.youtube.com/vi/de-YSJmq_5s/hqdefault.jpg)](https://www.youtube.com/watch?v=de-YSJmq_5s)

### panoptes API

Panoptes provides the following API endpoints:

Endpoint | Method | Description 
-- | -- | --
`/api/service-info` | `GET` | Server status
`/api/workflows` | `GET` | Get all workflows
`/api/workflow/<workflow-id>` | `GET` | Get workflow status
`/api/workflow/<workflow-id>/jobs` | `GET` | Get all jobs of a workflow
`/api/workflow/<workflow-id>/job/<job-id>` | `GET` | Get job status
`/api/workflow/<workflow-id>` | `PUT` | Rename a workflow  <br>  Expects a dictionary with new name <br> (e.g. `{'name': 'my new workflow name'}`)
`/api/workflow/<workflow-id>` | `DELETE` | Delete a workflow
`/api/workflows/all` | `DELETE` | Clean up database

To communicate with panoptes the following endpoints are used by snakemake:

Endpoint | Method | Description 
-- | -- | --
`/api/service-info` | `GET` | Server status (same as above)
`/create_workflow` | `GET` | Get a unique id/name str(uuid.uuid4()) for each workflow
`/update_workflow_status` | `POST` | Panoptes receives a dictionary from snakemake that contains: <br> - A log message dictionary <br> - The current timestamp <br> - The unique id/name of the workflow. <br> (e.g. `{'msg': repr(msg), 'timestamp': time.asctime(), 'id': id}`)

# Contribute

Please see the [Contributing instructions](CONTRIBUTING.md).

## CI server

Changes in develop or master trigger a [Travis](https://travis-ci.com/panoptes-organization/panoptes) build (and runs tests)

# Contact

In case the [issues section](https://github.com/panoptes-organization/panoptes/issues) is not enough for you, you can also contact us via [discord](https://discord.gg/vMcZCVZ)
# List of contributors (Alphabetical order, last name)

- Dimitrios Afentoulis (@dafentoulis)
- Argyrios-Alexandros Gardelakos (@agardelakos)
- Foivos Gypas (@fgypas)
- Georgios Kostoulas (@gkostoulas)
- Johannes Köster (@johanneskoester)
- Georgios Ntalaperas (@gntalaperas)
- Dimitrios Rekoumis (@drekoumis)
- Vanessa Sochat (@vsoch)
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
reported by contacting the project team at foivos DOT gypas AT unibas.ch. All
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
# Contributing to panoptes

First off, thanks for taking the time to contribute!

The following is a set of guidelines for contributing to panoptes. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## How Can I Contribute?

- Report bugs
- Propose or implement features
- Submit code changes / fixes
- Discuss code

See here for a [short tutorial for GitHub's issue tracking
system](https://guides.github.com/features/issues/).

Please adhere to the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting bugs

Please use the project's
[issue tracker](https://github.com/panoptes-organization/panoptes/issues) to report bugs. If you have no experience in filing bug reports, see e.g.,
[these recommendations by the Mozilla Developer Network](https://developer.mozilla.org/en-US/docs/Mozilla/QA/Bug_writing_guidelines)
first. Briefly, it is important that bug reports contain enough detail,
background and, if applicable, _minimal_ reproducible sample code. Tell us
what you expect to happen, and what actually does happen.

## Implementing features and submitting fixes

Kindly use pull requests to submit changes to the code base. But please note
that this project is driven by a community that likes to act on consensus. So
in your own best interest, before just firing off a pull request after a lot of
work, please [open an issue](https://github.com/panoptes-organization/panoptes/issues)
to **discuss your proposed changes first**. Afterwards, please stick to the
following simple rules to make sure your pull request will indeed be merged:


1. Fork the repo and create a [_feature
   branch_](https://datasift.github.io/gitflow/IntroducingGitFlow.html) from
   branch `develop`
2. If you've added code that should be tested, add tests.
3. Ensure that all tests pass.
4. Document your code and update all relevant documentation.
5. Stick to the code and documentation style (see below).
6. Issue the pull request.

Don't forget to add your name and GitHub profile URL to the
[list of contributors](contributors.md).

Important: Note that all your contributions are understood to be covered by the
[same license](LICENSE.md) that covers the entire project.
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

**Screenshots or logs**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
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
