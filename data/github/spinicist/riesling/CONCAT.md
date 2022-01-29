![Logo](riesling-logo.png)

[![Build](https://github.com/spinicist/riesling/workflows/Build/badge.svg)](https://github.com/spinicist/riesling/actions)
[![DOI](https://zenodo.org/badge/317237623.svg)](https://zenodo.org/badge/latestdoi/317237623)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03500/status.svg)](https://doi.org/10.21105/joss.03500)

## Radial Interstices Enable Speedy Low-Volume imagING

This is a reconstruction toolbox optimised for 3D ZTE MR images. There are many high quality MR recon toolboxes available, e.g. [BART](http://mrirecon.github.io/bart/), but these are mostly optimised for 2D sequences. 3D non-cartesian sequences present unique challenges for efficient reconstruction, so we wrote our own.

This toolbox was presented at ISMRM 2020.

## Authors

Tobias C Wood, Emil Ljungberg, Florian Wiesinger.

## Installation

Pre-compiled executables are provided for Linux and Mac OS X in a .tar.gz 
archive from http://github.com/spinicist/riesling/releases. Download the 
archive and extract it with `tar -xzf riesling-platform.tar.gz`. Then, move the 
resulting `riesling` executable to somewhere on your `$PATH`, for instance 
`/usr/local/bin`. That's it.

- MacOS Catalina or higher users should use `curl` to download the binary, i.e. 
  `curl -L https://github.com/spinicist/riesling/releases/download/v1.0/riesling-macos.tar.gz`. 
  This is because Safari now sets the quarantine attribute of all downloads, 
  which prevents them being run as the binary is unsigned. It is possible to 
  remove the quarantine flag with `xattr`, but downloading with `curl` is more 
  straightforward.
- The Linux executable is compiled on Ubuntu 16.04 with GLIBC version 2.3 and a 
  statically linked libc++. This means it will hopefully run on most modern 
  Linux distributions. Let us know if it doesn't.

## Compilation

If you wish to compile RIESLING yourself, compilation should hopefully be 
straightforward as long as you have access to a C++17 compiler (GCC 8 or higher,
Clang 7 or higher). RIESLING relies on `vcpkg` for dependency management. To 
download and compile RIESLING, follow these steps:

### 0. MacOS Dependencies
Install the [MacOS vcpkg dependencies](https://github.com/microsoft/vcpkg#installing-macos-developer-tools). This includes:
1. XCode from the AppStore
2. Run `$ xcode-select --install` in the terminal

You may also need to install `pkg-config` depending on your macOS version. This is easily installed with [Homebrew](https://brew.sh/) using
```
$ brew install pkg-config
```

**_(17/01/2022): The Apple Silicon Architecture (M1) is currently not supported. We are working on it though, stay tuned!)_**

### 0. Linux Dependencies
Install the [Linux vcpkg dependencies](https://github.com/microsoft/vcpkg#installing-linux-developer-tools). These includes:

These are `cmake`, `tar`, `curl`, `zip`, `unzip`, `pkg-config` \& `build-essentials`. You may be surprised by which distributions do not include these by default.

### 1. Clone repository
```
$ git clone https://github.com/spinicist/riesling
```

### 2. Compile
In the `riesling` folder execute
```
$ ./bootstraph.sh
```

## Usage

RIESLING comes as a single executable file with multiple commands, similar to 
`git` or `bart`. Type `riesling` to see a list of all the available commands. If you run a RIESLING command without any additional parameter RIESLING will output all available options for the given command.

RIESLING uses HDF5 (.h5) files but can also output NIFTI (.nii). To create an 
example digital phantom, use `riesling phantom`. RIESLING will append suffixes 
to input filenames when writing outputs to indicate which command was executed.

There are several reconstruction algorithms currently provided in RIESLING. 
Simple non-iterative recon is available with (`riesling recon`).

A separate examples repository https://github.com/spinicist/riesling-examples
contains Jupyter notebooks demonstrating most functionality. These can also be
run on the mybinder website.

## Documentation & Help

Documentation is available at https://riesling.readthedocs.io.

If you can't find an answer there or in the help strings, 
you can open an [issue](https://github.com/spinicist/riesling/issues), or find
the developers on Twitter ([@spinicist](https://twitter.com/spinicist)) or
e-mail tobias.wood@kcl.ac.uk.
# Contributing to RIESLING

Welcome to RIESLING. Any help is very much appreciated.

These guidelines exist to make it easy to get involved and ensure that RIESLING is welcoming to everyone.

## Asking Questions / Reporting Bugs

If you have general questions about how to use RIESLING programs, you can either open an [issue](https://github.com/spinicist/riesling/issues).

If you think you have found a bug, then opening a Github issue is the preferred avenue. Please check whether an identical or similar issue already exists first. When opening a new issue, please give as much information as you can, including
- The version of riesling you are using. Either give the version number if you downloaded the binaries, or a branch/git commit id if you compiled from source.
- The operating system you are running on.
- A description of the bug, including the full shell command you ran and any input files.

## Contributing Changes

If you want to edit the RIESLING code yourself to fix a problem or implement a new feature, you are very welcome to! Please follow this model for submitting a Pull Request:
- (Fork)[https://help.github.com/articles/fork-a-repo/] the riesling repo to your Github profile.
- Clone this repo, then create a new branch.
- Make your changes on this branch.
- Run clang-format. Riesling has a .clang-format file to avoid wasting time arguing over style.
- Squash/rebase your changes into a minimal number of commits.
- Make sure your copy of the `main` branch is up to date, and if necessary `rebase` your branch to the latest commit.
- Push your branch to Github.
- Open a Pull Request.

## Code of Conduct

- Anyone who participates in the development of RIESLING is expected to show respect and courtesy to all other community members at all times.
- Harrasment in any form towards any members will not be tolerated.
- All communication should be appropraite for a professional audience including people of different backgrounds.
- Be kind to others. Do not insult or put down other contributers.

## Thanks

Above all, thank you for using or contributing to RIESLING. It really is appreciated.
---
title: "Radial Interstices Enable Speedy Low-volume Imaging"
tags:
  - mri
  - reconstruction
  - ZTE
  - cpp
authors:
  - name: Tobias C Wood
    orcid: 0000-0001-7640-5520
    affiliation: 1
  - name: Emil Ljungberg
    orcid: 0000-0003-1456-7967
    affiliation: 1
  - name: Florian Wiesinger
    orcid: 0000-0002-5597-6057
    affiliation: "1,2"
affiliations:
  - name: Department of Neuroimaging, King's College London
    index: 1
  - name: GE Healthcare
    index: 2
date: 2021-06-25
bibliography: paper.bib
---

# Summary

- Standard MRI methods acquire Fourier-encoded data on a regularly spaced Cartesian grid. Noncartesian MRI methods abandon this grid and instead acquire data along arbitrary trajectories in k-space, which can lead to advantages such as motion robustness and more flexible pulse sequence design.
- Zero Echo Time (ZTE) imaging is a specialised form of MRI with an inherently noncartesian three-dimensional radial trajectory [@Ljungberg]. ZTE imaging has many interesting benefits, such as near silent operation, but also technical challenges such as dead-time gap artefacts [@Froidevaux].
- We have developed a toolbox, named Radial Interstices Enable Speedy Low-volume imagING (RIESLING), tuned for high performance reconstruction of MRI data acquired with 3D noncartesian k-space trajectories. While our group has a focus on ZTE imaging, RIESLING is suitable for all 3D noncartesian trajectories.

# Statement of Need

3D noncartesian trajectories can be challenging to reconstruct compared to other MRI trajectories. The major issues are:

- **Inseparability of the trajectory dimensions.** In Cartesian imaging it is often possible to consider one of the image dimensions separately to the others, leading to algorithmic complexity reduction and memory savings. 3D noncartesian must consider all dimensions simultaneously.
- **Oversampling requirements.** In order to obtain good image quality, it is necessary to oversample the reconstruction grid. Coupled with the inseparability of the 3D problem, this leads to large memory requirements.
- **Sample density compensation.** Noncartesian trajectories can lead to uneven k-space sampling, the correction of which is equivalent to preconditioning a linear system [@cgSENSE]. For 2D noncartesian trajectories, it is generally not necessary to correct for this in the context of an iterative reconstruction, as the problem will still converge. For 3D noncartesian trajectories, such preconditioning is essential to ensure reasonable convergence properties.

Existing MRI reconstruction toolboxes such as BART [@BART] or SigPy [@SigPy] provide high-quality implementations of many reconstruction algorithms but have prioritised applications other than 3D noncartesian reconstruction. We hence created a dedicated toolbox for this, including specific features such as:

- A configurable grid oversampling factor
- A thread-safe gridding implementation suitable for multi-core CPUs
- Integrated sample density compensation functions [@Zwart]

RIESLING (http://github.com/spinicist/riesling) is written with a modern C++ toolchain and utilizes Eigen [@Eigen] for all core operations, providing high performance. Data is stored in the HDF5 format allowing easy integration with other tools. We provide multiple reconstruction strategies, including classic conjugate-gradient SENSE [@cgSENSE] and Total-Generalized Variation [@TGV]. We are actively using RIESLING for our own studies and hope that it will be of interest to other groups using such sequences.

# References
# MATLAB Integration
Includes two main functions for reading and writing riesling k-space data

```
varargout = riesling_read(fname);
```

```
riesling_write(fname, data, traj, info, meta);
```

## Info struct

You can create an empty header info struct with

info = riesling_info();

## Examples
One examples is included:

1. `riesling_matlab_demo.m`: Demonstrating basic data input, output, manipulation and plotting.
