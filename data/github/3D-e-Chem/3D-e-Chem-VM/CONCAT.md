# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [1.4.4] - 2018-09-27

Based on [Chemical Analytics Platform](https://github.com/NLeSC/Chemical-Analytics-Platform) version 0.8.6

## [1.4.3] - 2018-02-09

### Changed

* 3D-e-Chem node collection v1.1.3
* Upgraded Virtualbox Guest Additions to v5.2.6
* Upgrade KNIME to v3.5.2
* Upgraded RDKit to 2017_09_3

### Fixed

- Made Python version 2 the default in KNIME

## [1.4.2] - 2017-12-06

### Changed

- adapted to new format of workflows repo

## [1.4.1] - 2017-11-12

Based on [Chemical Analytics Platform](https://github.com/NLeSC/Chemical-Analytics-Platform) version 0.8.5

### Added

- PLANTS example workflow
- Silicos-IT example workflows

### Changed

- Upgraded KripoDB to v2.3.1

## [1.4.0] - 2017-03-10

Based on [Chemical Analytics Platform](https://github.com/NLeSC/Chemical-Analytics-Platform) version 0.8.4

### Changed

- Upgraded KripoDB to v2.2.1

## Removed

- KripoDB PDB dataset, use web service to get latest dataset

## [1.3.0] - 2017-01-26

### Added

- KNIME molviewer example workflow (#27)

### Changed

- Switched to single 3D-e-Chem feature instead of every node separately (#28)
- Upgraded KNIME to v3.3.1
- Upgraded KripoDB to v2.1.0

## [1.2.0] - 2016-09-04

Based on [Chemical Analytics Platform](https://github.com/NLeSC/Chemical-Analytics-Platform) version 0.8.2

### Added

- Sygma Python library and KNIME node (#20)

### Changed

- Disabled Virtualbox 3D accleration, so webgl works
- 3D-e-Chem workflows from single GitHub repo (#21)
- Kripo fragments db source url (#22)
- Upgraded Knime to v3.2.1

### Fixed

- Kripodb sample data set renamed (#19)

## [1.1.1] - 2016-06-08

### Fixed

- Modified tanimoto node install fails (#17)

### Added

- Chemdb4VS workflow
- Klifs Knime nodes & example workflow (#14)
- Kripo fragments db for whole PDB
- Kripo bioisosteric replacement workflow (#16)
- Pymol sessions for Kripo and GPCRDB published (#18)

### Changed

- Updated KripoDB to 1.4.2
- Workflow examples are always downloaded to force latest version.

## [1.1.0] - 2016-05-07

### Changed

- Based on Chemical Analytics Platform 0.8.0
- Upgraded RDKit to 2016.03.1
- Upgraded Ubuntu to 16.04
- Upgraded Knime to v3.1.2
- Upgraded Virtualbox Guest Additions to v5.0.20

### Added

- ss-TEA Knime node

## [1.0.3] - 2016-05-03

### Added

- Kripodb GPCR sample dataset (#12)

### Changed

- Updated Knime to 3.1.2
- Updated Ubuntu to 14.04.4
- Updated KripoDB to 1.4.0

## [1.0.2] - 2016-01-23

### Fixed

- Allow vagrant user to update Knime
- self upgrade script gets overwritten by cap (#6)

## [1.0.1] - 2016-01-22

### Changed

- Based on Chemical Analytics Platform (https://github.com/NLeSC/Chemical-Analytics-Platform)

### Added

- 3D e Chem software, Knime nodes, Fingerprint cli

## [1.0.0] - 2016-01-16

Initial release with Knime, rdkit, postgresql, R and lxde.

### Added

- RDKit (#2)
- Postgresql (#1)
- Chembl postgresql script (#3)
# 3D e-Chem Virtual Machine

[![Build Status](https://travis-ci.org/3D-e-Chem/3D-e-Chem-VM.svg?branch=master)](https://travis-ci.org/3D-e-Chem/3D-e-Chem-VM)
[![DOI](https://zenodo.org/badge/19641/3D-e-Chem/3D-e-Chem-VM.svg)](https://zenodo.org/badge/latestdoi/19641/3D-e-Chem/3D-e-Chem-VM)

Scripts to create a Vagrant box using packer and ansible.

For available software/datasets/workflows inside Virtual machine see https://github.com/3D-e-Chem/3D-e-Chem-VM/wiki

# Usage

* VirtualBox, https://www.virtualbox.org
* Vagrant, https://www.vagrantup.com

Create a new directory and cd to it then start virtual machine with

```
vagrant init nlesc/3d-e-chem
vagrant up
```

Usage screencast on YouTube:

[![3D-e-Chem Virtual Machine screencast on YouTube](https://img.youtube.com/vi/zBv4rPfsLLc/0.jpg)](https://www.youtube.com/watch?v=zBv4rPfsLLc)

# Build

Requirements:

* Virtualbox, https://www.virtualbox.org/
* Packer, https://packer.io
* Enough disk space
  * Make sure temporary directory (/tmp by default on Linux) has enough space. Use TMPDIR environment variable to overwrite default location
* ova file `../Chemical-Analytics-Platform/output-virtualbox-iso/cap.ova` from build phase of https://github.com/NLeSC/Chemical-Analytics-Platform

```
packer build packer.json
```
# Test

Add box to Vagrant with

```
vagrant box remove --force --all nlesc/3d-e-chem
vagrant box add --name nlesc/3d-e-chem --force packer_virtualbox-ovf_virtualbox.box
```

Then use steps described at Usage chapter in a new directory.

# Push

Requirements:

* Vagrant cloud account, https://app.vagrantup.com/nlesc/boxes/3d-e-chem

Publish box on https://app.vagrantup.com/nlesc/boxes/3d-e-chem using the following steps:

1. Create a new version
2. Create a new provider
3. Choose `virtualbox` as provider
4. Choose Upload
5. Press `Continue to upload`
6. Upload the `packer_virtualbox-ivf_virtualbox.box` file generated by `vagrant package`
7. Edit version
8. Press `Release version`

# Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md).

# Cite

When using 3D-e-Chem-VM please cite one of the following:

* [Zenodo software release DOI](https://zenodo.org/badge/latestdoi/19641/3D-e-Chem/3D-e-Chem-VM)
* R. McGuire, S. Verhoeven, M. Vass, G. Vriend, I. J. P. De Esch, S. J. Lusher, R. Leurs, L. Ridder, A. J. Kooistra, T. Ritschel, C. de Graaf (2017). 3D-e-Chem-VM: Structural cheminformatics research infrastructure in a freely available Virtual Machine. Journal of Chemical Information and Modeling. doi: http://dx.doi.org/10.1021/acs.jcim.6b00686
# Contributor Covenant Code of Conduct

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
reported by contacting the project team at s.verhoeven(at)esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing

We love pull requests from everyone. By participating in this project, you
agree to abide by the [code of conduct](CODE_OF_CONDUCT.md).

Fork, then clone the repo:

    git clone git@github.com:your-username/3D-e-Chem-VM.git

Set up your machine:

    pip install ansible

And install packer see https://www.packer.io/downloads.html for instructions.

Make sure the syntax checks pass:

    ansible-playbook -i localhost, --syntax-check playbook.yml
    packer validate -syntax-only packer.json

Make your change. Make the syntax checks pass:

    ansible-playbook -i localhost, --syntax-check playbook.yml
    packer validate -syntax-only packer.json

Optionally test the change works by building a Vagrant box and testing it manually as explained in [README.md#Build](README.md#build) and [README.md#Test](README.md#test) chapters respectively.

Push to your fork and [submit a pull request](https://github.com/3D-e-Chem/3D-e-Chem-VM/compare/).

At this point you're waiting on us. We like to at least comment on pull requests
within three business days (and, typically, one business day). We may suggest
some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted:

* Write a [good commit message](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
