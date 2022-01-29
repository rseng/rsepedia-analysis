[![Build status](https://github.com/kutaslab/fitgrid/actions/workflows/fitgrid-cid.yml/badge.svg)](https://github.com/kutaslab/fitgrid/actions)
[![Coverage](https://codecov.io/gh/kutaslab/fitgrid/branch/main/graph/badge.svg)](https://codecov.io/gh/kutaslab/fitgrid)
[![DOI](https://zenodo.org/badge/147436563.svg)](https://zenodo.org/badge/latestdoi/147436563)

# fitgrid

A Python library for regression modeling time-varying patterns of activity in sensor-array data streams on a 2-D grid.

We gratefully acknowledge the support of grant NICHD 5R01HD022614 for the development of these routines.

## Documentation

User guide, installation instructions, workflow and usage examples are available [here](https://kutaslab.github.io/fitgrid).

## Demo

Click this button to launch a demo notebook:

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/kutaslab/fitgrid/main?filepath=notebooks/Demo.ipynb)
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
reported by contacting the project team at turbach@ucsd.edu. All
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
Gold standard test data files have been moved to
`fitgrid/fitgrid/data`.

This directory is retained as landing site for temp fake data files
generated during testing.

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

Thank you for helping to improve fitgrid.

## Before you submit an issue please check your fitgrid installation.

Problems are unavoidable when fitgrid is installed with incompatible or missing
Python or R packages.

Please try to replicate the issue after first using mamba or conda to install the latest fitgrid [stable version](https://kutaslab.github.io/fitgrid/installation.html#fitgrid-stable-release) into a newly created conda virtual environment.

If the problem persists, please try to replicate the issue after installing the latest [pre-release version](https://kutaslab.github.io/fitgrid/installation.html#fitgrid-development-version) into a newly created conda virtual environment.

**Warning: `pip install fitgrid` is officially not supported.**


## Please provide the following information

### 1. Description
A clear and concise description of what the problem is, specifically:
- What you expected to happen and what actually happened.
- Anything you tried to solve the issue.

### 2. Minimal reproducible example
These are the shortest steps that reconstruct the problem.
- The **exact** character-for-character command(s) you ran to to install fitgrid, copy-paste is best.
- A Python code snippet or shell commands, the shorter the better, that runs and exhibits the issue.

### 3. Conda environment
Activate the conda environment that has fitgrid installed, run the following command in a terminal window, and upload the `fitgrid_issue.txt` file as an attachment with your issue.
```
conda list --explicit > fitgrid_issue.txt
```

### 4. System information
Please provide the specifics about your computer hardware architecture and
operating system version. For example:

- Linux, in a terminal window

```
	$ uname -mprsv
	Linux 3.10.0-957.21.3.el7.x86_64 #1 SMP Tue Jun 18 16:35:19 UTC 2019 x86_64 x86_64
	
	$ cat /etc/*-release
	CentOS Linux release 7.3.1611 (Core) 
	NAME="CentOS Linux"
	VERSION="7 (Core)"
	ID="centos"
	ID_LIKE="rhel fedora"
	VERSION_ID="7"
	PRETTY_NAME="CentOS Linux 7 (Core)"
	ANSI_COLOR="0;31"
	CPE_NAME="cpe:/o:centos:centos:7"
	HOME_URL="https://www.centos.org/"
	BUG_REPORT_URL="https://bugs.centos.org/"
	
	CENTOS_MANTISBT_PROJECT="CentOS-7"
	CENTOS_MANTISBT_PROJECT_VERSION="7"
	REDHAT_SUPPORT_PRODUCT="centos"
	REDHAT_SUPPORT_PRODUCT_VERSION="7"
	
	CentOS Linux release 7.3.1611 (Core) 
	CentOS Linux release 7.3.1611 (Core) 
```


- Mac OSX, in a terminal window

```
	$ uname -mprsv
	Darwin 19.6.0 Darwin Kernel Version 19.6.0: Thu Jun 18 20:49:00 PDT 2020; root:xnu-6153.141.1~1/RELEASE_X86_64 x86_64 i386

	$ sw_vers
	ProductName:	Mac OS X
	ProductVersion:	10.15.6
	BuildVersion:	19G2021
```

  
- [TODO: not officially supported] Windows, in a Windows Command Window

```
	C:\Users\some_user> systeminfo

	Host Name:                 DESKTOP-G57OVSM
	OS Name:                   Microsoft Windows 10 Home
	OS Version:                10.0.18362 N/A Build 18362
	OS Manufacturer:           Microsoft Corporation
	OS Configuration:          Standalone Workstation
	OS Build Type:             Multiprocessor Free
	Registered Owner:          some_user
	Registered Organization:
	Product ID:                00326-00840-79774-AAOEM
	Original Install Date:     7/29/2019, 6:13:59 AM
	System Boot Time:          9/9/2020, 9:07:46 AM
	System Manufacturer:       System manufacturer
	System Model:              System Product Name
	System Type:               x64-based PC
	Processor(s):              1 Processor(s) Installed.
    [01]: Intel64 Family 6 Model 158 Stepping 12 GenuineIntel ~3600 Mhz
	BIOS Version:              American Megatrends Inc. 0606, 8/31/2018
	Windows Directory:         C:\WINDOWS
	System Directory:          C:\WINDOWS\system32
	Boot Device:               \Device\HarddiskVolume2
	System Locale:             en-us;English (United States)
	Input Locale:              en-us;English (United States)
	Time Zone:                 (UTC-08:00) Pacific Time (US & Canada)
	Total Physical Memory:     16,305 MB
	Available Physical Memory: 13,598 MB
	Virtual Memory: Max Size:  18,737 MB
	Virtual Memory: Available: 14,542 MB
	Virtual Memory: In Use:    4,195 MB
```
