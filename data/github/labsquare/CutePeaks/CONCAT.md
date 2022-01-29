# CutePeaks

[![C/C++ CI](https://github.com/labsquare/CutePeaks/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/labsquare/CutePeaks/actions/workflows/c-cpp.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5148809.svg)](https://doi.org/10.5281/zenodo.5148809)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03457/status.svg)](https://doi.org/10.21105/joss.03457)


https://labsquare.github.io/CutePeaks/

A simple viewer for Sanger trace file made with Qt5.
Supports AB1 and SCF 3.0 file formats.
It has regular expression pattern finder and can export trace as SVG vector image.

![Preview](https://raw.githubusercontent.com/labsquare/CutePeaks/master/cutepeaks.gif)


# Statement of need
Despite the major use of Next Generation Sequencing, the Sanger method is still widely used in genetic labs as the gold standard to read target DNA sequences. Very few opensource software is available to explore Sanger trace data and most of labs staff still rely on proprietary software. Moreover, they are not always user-friendly and lack modern look and feel.

# State of field
[4peaks](https://nucleobytes.com/4peaks/) is a software widely used by biologists that benefits from a nice User interface. Sadly, it is only available on MacOS and source code is not opened to community enhancement. [Seqtrace](https://github.com/stuckyb/seqtrace) is the only standalone and opensource application we could find. However, it is written with GTK framework in Python 2, the latter being deprecated and slower than C++.

# Installation
## Windows 
[Download windows binary](https://github.com/labsquare/CutePeaks/releases/download/0.2.3/CutePeaks-win32.exe)

## MacOSX 
[Downlad MacOS binary](https://github.com/labsquare/CutePeaks/releases/download/0.2.3/cutepeaks-macos.dmg)

## Linux
Linux binary is available as [AppImage](http://appimage.org/).
Download the AppImage from [here](https://github.com/labsquare/CutePeaks/releases).
For ubuntu 21.04, Download this one [here](https://github.com/labsquare/CutePeaks/releases/download/0.2.3/cutepeaks-ubuntu_21-04-x86_64.AppImage)

Run it as follow:


    chmod +x cutepeaks-0.2.0-linux-x86_64.appimage
    ./cutepeaks-0.2.0-linux-x86_64.appimage


## Compilation
### Prerequisites
#### Install Qt ≥ 5.7

**From website**: Download Qt ≥ 5.7 from https://www.qt.io/.
Don't forget to check QtChart module during installation.

**From Ubuntu**: Qt 5.7 is not yet available with Ubuntu. But you can add a PPA to your software system.
For exemple for Xenial:

    sudo add-apt-repository ppa:beineri/opt-qt57-xenial
    sudo apt install qt57base qt57charts-no-lgpl
    source /opt/qt57/bin/qt57-env.sh

**From Fedora**: Qt 5.7 is available.

    sudo dnf install qt5-qtbase-devel qt5-qtcharts-devel

### Compile CutePeaks
Be sure you have the correct version of Qt (≥ 5.7) by using qmake --version. For exemple, if you have installed Qt from ppa:beineri, you will find it under /opt/qt57/bin/qmake. Then launch the compilation from CutePeaks folder as follow.

     /opt/qt57/bin/qmake --version
     /opt/qt57/bin/qmake
     make
     sudo make install

## Usage
CutePeaks supports following trace file formats:

- *.ab1
- *.scf

Example files are available here: 
https://github.com/labsquare/CutePeaks/tree/master/examples       
You can open those files from cutepeaks by clicking on *open* from the File menu.

## Features 
Once the file is open, cutepeaks allows you to : 
- Explore the trace from a scroll area. ( Finger geasture are supported with touch screen) 
- Scale the trace horizontally or vertically using 2 sliders at the bottom right.
- Select a subsequence with the mouse as with any text editor. Then you can cut or copy to the clipboard
- Make the reverse complement from the edit menu
- Display Sequence and metadata from the view menu
- Search for a regular expression pattern. Open the "Find Sequence..." from the edit menu
- Export trace or sequence to different format. ( e.g: Fasta, CSV, SVG or PNG image ) 

## Contributions / Bugs
See [how to contributing](https://github.com/labsquare/CutePeaks/edit/master/CONTRIBUTING.md)

## Licenses
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/gpl-3.0.txt.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
sacha@labsquare.org
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing
All assistance is welcome. You can contribute to the project from the following topic:

# Issues
Bug or feature requests can report from [github issue tracker](https://github.com/labsquare/CutePeaks).

# Pull request
Fell free to create pull request. Compilation test from [github](https://github.com/labsquare/CutePeaks/actions/workflows/c-cpp.yml) action are triggered from any pull request.

# Coding style
Please respect the Qt coding style described [here](https://wiki.qt.io/Qt_Coding_Style)

# Translation
Cutepeaks is only translated in French. You can translate to your language by using [Qt linguist](https://doc.qt.io/qt-5/linguist-manager.html)

# Chat
You can join us [on discord](https://discord.com/invite/7sSH4VSPKK). We are speaking french right now, but we can switch to english.
---
title: 'CutePeaks: A modern viewer for Sanger trace file'
tags:
  - Python
  - Sanger
  - genetics
  - Qt
  - GUI
authors:
  - name: Sacha Schutz^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-4563-7537
    affiliation: "1,2" # (Multiple affiliations must be quoted)
  - name: Charles Monod-Broca^[co-first author]
    orcid: 0000-0003-4095-8099
    affiliation: "2"
  - name: Anne-Sophie Denommé-Pichon
    orcid: 0000-0002-8986-8222
    affiliation: "3,4"

affiliations:
 - name: CHRU Brest, Hôpital Morvan, Laboratoire de Génétique Moléculaire, Brest, France
   index: 1
 - name: Univ Brest, Inserm, EFS, UMR 1078, GGB, 29200
   index: 2
 - name: Unité Fonctionnelle Innovation en Diagnostic génomique des maladies rares, FHU-TRANSLAD, CHU Dijon Bourgogne, Dijon, France
   index: 3
 - name: UMR1231 GAD, Inserm - Université Bourgogne-Franche Comté, Dijon, France
   index: 4

date: 15 June 2021
bibliography: paper.bib

---

# Summary
CutePeaks is a standalone Sanger trace viewer steered by a modern and user-friendly UI. Unlike other software, CutePeaks comes with two new features: searching for a regular expression and exporting the traces to SVG.    
CutePeaks is available for Linux, macOS and Windows at [https://labsquare.github.io/CutePeaks/](https://labsquare.github.io/CutePeaks/).

# Statement of need
Despite the major use of next generation sequencing, the Sanger method is still widely used in genetic labs as the gold standard to read target DNA sequences. Very few open source software is available to explore Sanger trace data and most of labs staff still rely on proprietary software. Moreover, they are not always user-friendly and lack modern look and feel. 

# State of fields
4peaks [@4Peaks] is software widely used by biologists that benefits from a nice user interface. Sadly, it is only available on macOS and source code is not open to community enhancement. Seqtrace [@seqtrace] is the only standalone and open source application we could find. However, it is written with the GTK framework in Python 2, the latter being deprecated and slower than C++. 

# Software overview
![CutePeaks screenshot with regular expression search bar.\label{fig:example}](figure.png)

CutePeaks is a cross-platform application implemented in C++ using the open source Qt5 framework. It can read FSA and ABIF file formats, and display the chromatogram with standard controllers.
The chromatogram is displayed in an interactive window allowing the user to move along the trace. It can also re-scale the plot dynamically using two slider controllers. Finger gestures are also available for scrolling upon using a touch screen.
Similarly to 4peaks software [@4Peaks], Phred quality scores are displayed behind the trace as a blue histogram. Base calling is displayed at the top of the viewing window, along with adjustable amino-acid translation.
The trace can be used as with a text editor. Navigating along the trace, copying the sequence to the clipboard or cutting it is done using standard keyboard shortcuts. Revert/complement is also possible.
An original feature of CutePeaks is the possibility to search for a sequence in the trace using a regular expression. This is especially useful to search for a sequence pattern. For example, the query A[CG]T will search for all instances of ACT or AGT. The query AC+T will select all instances of the form ACT, ACCT, ACCCCT, etc. Finally, the trace data can be exported to different formats, such as FASTA or SVG image, the latter being particularly useful for resolution-independent illustration.


# Installation

CutePeaks is hosted on the GitHub development platform. Continuous integration is provided by GitHub Actions.
For Linux, an AppImage is provided, that is, distribution agnostic and runs out of the box.
For Windows, a 32 bits binary compiled with mingw is provided and can be executed as a standalone application without administrator privileges. For macOS, a disk image is provided.

# Acknowledgements

We acknowledge contributions from Jérémie Roquet, and Francisco Pina-Martins.

# References
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
