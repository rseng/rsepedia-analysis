Babuji, Y., Woodard, A., Li, Z., Katz, D. S., Clifford, B., Kumar, R., Lacinski, L., Chard, R., Wozniak, J., Foster, I., Wilde, M., and Chard, K., Parsl: Pervasive Parallel Programming in Python. 28th ACM International Symposium on High-Performance Parallel and Distributed Computing (HPDC). 2019.

or

```{tex}
@inproceedings{babuji19parsl,
  author       = {Babuji, Yadu and
                  Woodard, Anna and
                  Li, Zhuozhao and
                  Katz, Daniel S. and
                  Clifford, Ben and
                  Kumar, Rohan and
                  Lacinski, Lukasz and
                  Chard, Ryan and 
                  Wozniak, Justin and
                  Foster, Ian and 
                  Wilde, Mike and
                  Chard, Kyle},
  title        = {Parsl: Pervasive Parallel Programming in Python},
  booktitle    = {28th ACM International Symposium on High-Performance Parallel and Distributed Computing (HPDC)},
  doi          = {10.1145/3307681.3325400},
  year         = {2019},
  url          = {https://doi.org/10.1145/3307681.3325400}
}
```
# Parsl Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free and bullying-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

Discussions relating to pros/cons of various technologies, programming languages, and so on are welcome,
but these should be done with respect, taking proactive measure to ensure that all participants are heard
and feel confident that they can freely express their opinions.

We pledge to welcome questions and answer them respectfully, paying particular attention to those new to
the community. We pledge to provide respectful criticisms and feedback in forums, especially in discussion
threads resulting from code contributions.

We pledge to be conscientious of the perceptions of the wider community and to respond to criticism respectfully.
We will strive to model behaviors that encourage productive debate and disagreement, both within our community
and where we are criticized. We will treat those outside our community with the same respect as people within
our community.

We pledge to help the entire community follow the code of conduct, and to not remain silent when we see violations
of the code of conduct. We will take action when members of our community violate this code such as contacting
<a href="mailto:parsl-coc@googlegroups.com">parsl-coc@googlegroups.com</a> (all emails sent to this address will be treated with the strictest confidence)
or talking privately with the person.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members
* Respecting the work of others by recognizing acknowledgment/citation requests of original authors
* Being explicit about how we want our own work to be cited or acknowledged

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Sexist, racist, or otherwise exclusionary jokes
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
reported by contacting the project team at <a href="mailto:parsl-coc@googlegroups.com">parsl-coc@googlegroups.com</a>. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html,
and the [yt-project](https://yt-project.org)'s Code of Conduct.
Running the MPI Test
====================

This simple MPI test is designed to give you some very basic information:
    1. What are the ranks running on ?
    2. Args passed by a launch mechanism printed as argv[1], and argv[2]


Compile
=======

1. Load the appropriate MPI modules for your system
2. Compile the code with make:

    make clean; make

Running the app
===============

Make sure the right MPI modules is loaded. Run the app as a simple executable:

    ./mpi_hello

Or launch it with N ranks with an mpi launcher like mpirun:

    mpirun -n 8 mpi_hello# Description

Please include a summary of the change and (optionally) which issue is fixed. Please also include
relevant motivation and context.

Fixes # (issue)

## Type of change

Choose which options apply, and delete the ones which do not apply.

- Bug fix (non-breaking change that fixes an issue)
- New feature (non-breaking change that adds functionality)
- Breaking change (fix or feature that would cause existing functionality to not work as expected)
- Documentation update
- Code maintentance/cleanup
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior, for e.g:
1. Setup Parsl 0.8.0 with Python 3.6 on cluster
2. Run a test script
3. Wait 5 mins
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Environment**
 - OS: [e.g. ubuntu, centos, MacOS, windows]
 - Python version
 - Parsl version

**Distributed Environment**
- Where are you running the Parsl script from ? [e.g. Laptop/Workstation, Login node, Compute node]
- Where do you need the workers to run ? [e.g. Same as Parsl script, Compute nodes, Cloud nodes]
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. E.g. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context about the feature request here.
