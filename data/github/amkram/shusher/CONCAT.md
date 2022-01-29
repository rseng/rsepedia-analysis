---
title: 'ShUShER: private browser-based placement of sensitive genome samples on phylogenetic trees'
tags:
  - phylogenetics
  - JavaScript
  - WebAssembly
  - React
  - SARS-CoV-2
authors:
  - name: Alexander Kramer
    orcid: 0000-0003-3630-3209
    affiliation: "1, 2"
  - name: Yatish Turakhia
    orcid: 0000-0001-5600-2900
    affiliation: "3"
  - name: Russell Corbett-Detig
    orcid: 0000-0001-6535-2478
    affiliation: "1, 2"
affiliations:
 - name: Department of Biomolecular Engineering, University of California Santa Cruz. Santa Cruz, CA 95064, USA
   index: 1
 - name: Genomics Institute, University of California Santa Cruz, Santa Cruz, CA 95064, USA
   index: 2
 - name: Department of Electrical and Computer Engineering, University of California, San Diego; San Diego, CA 92093, USA
   index: 3
date: 26 August 2021
bibliography: paper.bib
---

# Summary

[ShUShER](https://github.com/amkram/shusher) (Shh: private Ultrafast Sample placement on Existing tRees) is a browser-based application for local analysis of sensitive biological sequence data. It uses UShER, a previously developed algorithm [@Turakhia:2021], to place user-provided genome samples on an existing phylogenetic tree and return subtrees surrounding those samples. It displays visualizations of the results using Auspice, part of the Nextstrain project [@Hadfield:2018]. UCSC hosts a [web tool](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) with similar functionality, but it performs server-side computation on user samples, making it unsuitable for samples containing Protected Health Information (PHI) or other sensitive data. There are two main components of ShUShER: a port of the existing C++ UShER code to WebAssembly and a user interface built with React, which together perform computation on user-provided samples entirely client-side in a web browser. The web application is accessible at https://shusher.gi.ucsc.edu/.


# Statement of need

Phylogenetic trees are often used to help trace the origin, spread, and evolution of viruses. The continuously growing number of sequenced SARS-CoV-2 genomes has quickly overwhelmed the capabilities of many existing tree construction methods. UShER is a method that can efficiently place newly sequenced genomes on existing, large phylogenetic trees. Many researchers have been using the UShER web tool hosted at UCSC to place their samples on existing global trees comprising millions of SARS-CoV-2 sequences. However, some jurisdictions consider viral genomes PHI, in which case the viral sequence data cannot be transmitted over the Internet. Presently, researchers cannot use the UCSC web tool for such data and must instead install and run UShER locally with a command-line application. ShUShER is an alternative to this, providing a user-friendly platform for researchers to run UShER in their web browser using an existing tree of sequences while keeping their data private. Currently, ShUShER supports placement of user samples on a [global tree](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) of publicly available SARS-CoV-2 samples [@McBroome:2021], and it may be extended in the future to support private analysis of other pathogen sequences.

# Acknowledgements
ShUShER uses or adapts code from several open source projects. The authors thank the many developers of [Nextclade](https://github.com/nextstrain/nextclade), [Auspice](https://github.com/nextstrain/auspice), [auspice.us](https://github.com/nextstrain/auspice.us), and [UShER](https://github.com/yatisht/usher). We also thank Angie Hinrichs for developing and maintaining essential web infrastructure upon which ShUShER relies. This work was supported by CDC award BAA 200-2021-11554 and by a gift from the Eric and Wendy Schmidt Foundation. 

# References
<div align="center">
  <img src="web-app/public/img/logo.png" height=100/>
</div>
<div align="center">
  private (<strong>Sh</strong>h :shushing_face:) <strong>U</strong>ltrafast <strong>S</strong>ample placement on <strong>E</strong>xisting t<strong>R</strong>ees</strong>
</div>
<div align="center">
  
 <br />
  
  [![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](LICENSE)
  [![Integration Tests](https://github.com/amkram/shusher/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/amkram/shusher/actions/workflows/build_and_test.yml)<!---BEGIN_USHER_BADGE-->
<a target="_blank" href="https://github.com/yatisht/usher/tree/e33eec1198c708a69ff8aafcb86072a98b8fee7a"><img src="https://img.shields.io/badge/UShER%20Version-commit%20e33eec11-%235e0000"/></a> <a href="https://joss.theoj.org/papers/03edfaa561a1cfbc53be7d98c8461cf3"><img src="https://joss.theoj.org/papers/03edfaa561a1cfbc53be7d98c8461cf3/status.svg"></a>
<!---END_USHER_BADGE-->
 
  | :computer_mouse:	Access ShUShER <a target="_blank" href="https://shusher.gi.ucsc.edu">here</a>! |
| --- |
</div>
<div align="center">
    ShUShER is a browser tool for placing sensitive genome sequences on phylogenetic trees using <a target="_blank" href="https://github.com/yatisht/usher">UShER</a>.
  <h3>
    <a href="#usage">
      Usage
    </a>
    <span> | </span>
    <a href="#how-it-works">
      How it works
    </a>
    <span> | </span>
    <a href="#installation-for-developers">
      Installation
    </a>
  </h3>
</div>


## Contents
- [Usage](#usage)
  - [Loading samples](#loading-samples)
  - [Running UShER](#running-usher)
  - [Interpreting results](#interpreting-results)
  - [Visualizing subtrees](#visualizing-subtrees)
  - [Downloading data](#downloading-data)
- [How it works](#how-it-works)
- [Installation](#installation-for-developers)
  - [Running the web app locally](#running-the-web-app-locally)
  - [Compiling UShER to WebAssembly](#compiling-usher-to-webassembly)
- [Contributing](#contributing)
- [About this repository](#about-this-repository)
- [Acknowledgements](#acknowledgements)

## Usage
> :warning:	This tool is intended to be used <strong>only for sequences that cannot be shared publicly</strong>. If you do not have this requirement, please use the [UShER web tool](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) and submit your sequences to an INSDC member institution (NCBI, EMBL-EBI, or DDBJ) and GISAID

ShUShER is currently designed for use with SARS-CoV-2 genomes. The user supplies a set of samples in FASTA or VCF format, and the provided samples are placed on a continuously growing global tree ([read more](https://www.biorxiv.org/content/10.1101/2021.04.03.438321v1)). After placement, subtrees containing user samples can be visualized (using [Auspice](https://docs.nextstrain.org/projects/auspice/en/stable/)).

### Loading samples
Samples can be provided to ShUShER in either FASTA (`.fa`, `.fasta`, `.fna`) or VCF format (`.vcf`).

All samples must be in a single file. When you load your samples into ShUShER, they will not be uploaded to our servers and the data will remain on your computer.

### Running UShER

After loading your samples, two input fields will appear:

  | <img src="data/run_usher_input.png" height=120/> |
| --- |

The first field selects the existing tree to place your samples on. Currently, the only option is the global SARS-CoV-2 tree (maintained [here](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/)).

After UShER places your samples on the global tree, it will output subtrees containing your samples. The second field allows you to select how many closely-related samples from the global tree to include in each subtree.
### Interpreting results

After UShER has finished running, a table of information about your samples will be displayed.

Two numbers are reported for each sample:

>**Number of maximally parsimonious placements** is the number of potential placements in the tree with minimal parsimony score. A higher number indicates a less confident placement.

>**Parsimony score** is the number of mutations/changes that must be added to the tree when placing this sample. The higher the number, the more diverged the sample.

### Visualizing subtrees
Each sample in the table will have a button, e.g.
<img src="data/subtree_button.png" height=45/>
allowing you to open the subtree containing that sample in Auspice. The subtree visualization will open in a new browser tab (but data is not sent over the Internet).

### Downloading data
Newick files for each of the generated subtrees can be downloaded at the bottom of the Auspice visualization page.
<img src="data/download_data_full.png" height=400/>

## How it works

The ShUShER web app uses a ported version of UShER that can be run client-side in a web browser. The original C++ [code base](https://github.com/yatisht/usher) is compiled to WebAssembly with [Emscripten](https://emscripten.org/) and wrapped in a React frontend (read more about the port [here](https://github.com/amkram/shusher/tree/master/usher-port)). User-provided samples are not transmitted across the Internet, and computation is performed locally in the browser. We use a modified version of [Auspice](https://docs.nextstrain.org/projects/auspice/en/stable/) to display the subtrees computed by UShER. The visualization opens in a new browser tab, using [localStorage](https://developer.mozilla.org/en-US/docs/Web/API/Window/localStorage) to share data between tabs without transmitting any user data over the web.

FASTA to VCF conversion is performed by aligning each provided sample pairwise to the reference SARS-CoV-2 genome. The implementation of pairwise alignment is from [Nextclade](https://github.com/nextstrain/nextclade/blob/0ed4e6a1569dbd0b91e9d4861494e97861a11e7e/packages/web/src/algorithms/alignPairwise.ts).

## Installation (for developers)
>SHUShER currently only supports building on Linux systems, and has been tested on Ubuntu 20.04

If you would like to run ShUShER locally or modify the source, first download the source code, e.g.:
  
  `wget https://github.com/amkram/shusher/archive/refs/tags/latest.tar.gz`
  
  `tar xvzf latest.tar.gz`

The above command will download the latest tagged release of ShUShER. View all "Releases" in the right sidebar if you want to download a specific version. Alternately, cloning this repository will give you the latest, unreleased code, but may be unstable.

The downloaded source code contains code for building both the web app and the UShER port.

### Running the web app locally

Enter the `web-app` subdirectory and run

  `npm install`

To build the app, run

  `npm run build`
  
And to start the local server, run

  `npm start`
  
You should now be able to access ShUShER in your browser at `localhost:4000`

### Compiling UShER to WebAssembly

The directory `usher-port` contains the original C++ UShER code and a script that will compile it to WebAssembly. You only need to compile UShER yourself if you want to change the UShER source code. Otherwise, the web app will automatically use the most recent pre-compiled release from this repository.

#### 1. Install Dependencies

`sudo apt-get update`

`sudo apt-get install wget python3 build-essential cmake protobuf-compiler dh-autoreconf`

#### 2. Compile UShER 

`./installUbuntuWeb.sh`

This script will download the C++ library dependencies of UShER, make some modifications necessary for WebAssembly compilation, and then compile them using emscripten. Output in the `build` directory includes `usher.wasm`, `usher.js`, `usher.data`, and `usher.worker.js`, all of which are used by the ShUShER web app.

#### 3. Specify custom UShER code

By default, the web app grabs the latest tagged release of the WebAssembly UShER bundle from this repository. If you compiled UShER yourself using the above steps, you can tell ShUShER to use your compiled code instead.

In the `web-app` subdirectory, edit `package.json` and change the following line:

    config: {
      usherBundle: "latest"
    }

to

    config: {
      usherBundle: "[path to build output]"
    }

## Contributing
We welcome and encourage contributions to ShUShER from the community. If you would like to contribute, please read the contribution [guidelines](CONTRIBUTING.md) and [code of conduct](CODE_OF_CONDUCT.md).

## About this repository
`usher-port` contains the scripts and files needed to compile UShER to WebAssembly. See [here](https://github.com/amkram/shusher/tree/master/usher-port) for details on the process.


`web-app` contains the React application that uses the UShER port.

Twice a day, the UShER C++ source hosted in this repository is updated from the main [UShER repository](https://github.com/yatisht/usher).

Upon each push to the master branch, the integration test Github Action is run, which (1) compiles the latest source from the main UShER repo to a binary executable, (2) compiles UShER to WebAssembly with this repo's latest code, and (3) runs both on a sample file and compares the outputs, ensuring they are the same.

New releases are tagged periodically and pushed to the live web app.

## Acknowledgements
This project uses or adapts code from several open-source projects. We are grateful for their contributions.


Pairwise sequence alignment uses the implementation from [Nextclade](https://github.com/nextstrain/nextclade).


Visualization of subtrees is performed with [Auspice](https://github.com/nextstrain/auspice/blob/5132a1c1d063761eb02dc5434a8316c6d5be7085/docs/index.rst).


Scripts to modify the Auspice server are from [auspice.us](https://github.com/nextstrain/auspice.us).


Nextclade, Auspice, and auspice.us are part of the [Nextstrain](https://github.com/nextstrain) project.


The core functionality of this tool is a ported version of [UShER](https://github.com/yatisht/usher).
# Code of Conduct - ShUShER

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

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
advances
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they deem
inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at .
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://contributor-covenant.org/), version
[1.4](https://www.contributor-covenant.org/version/1/4/code-of-conduct/code_of_conduct.md) and
[2.0](https://www.contributor-covenant.org/version/2/0/code_of_conduct/code_of_conduct.md),
and was generated by [contributing-gen](https://github.com/bttger/contributing-gen).
<!-- omit in toc -->
# Contributing to ShUShER

We welcome and encourage contributions to ShUShER from the community. If you would like to contribute, please read the below guidelines. 

<!-- omit in toc -->
## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [I Have a Question](#i-have-a-question)
- [I Want To Contribute](#i-want-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)

## Code of Conduct

This project and everyone participating in it is governed by the
[ShUShER Code of Conduct](https://github.com/amkram/shusherblob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report unacceptable behavior
to [alex.kramer@ucsc.edu](mailto:alex.kramer@ucsc.edu).


## I Have a Question

> If you want to ask a question, we assume that you have read the available [Documentation](https://github.com/amkram/shusher/blob/master/README.md).

Before you ask a question, it is best to search for existing [Issues](https://github.com/amkram/shusherissues) that might help you. In case you have found a suitable issue and still need clarification, you can write your question in this issue. It is also advisable to search the internet for answers first.

If you then still feel the need to ask a question and need clarification, we recommend the following:

- Open an [Issue](https://github.com/amkram/shusher/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions (nodejs, npm, browser, etc.), depending on what seems relevant.

We will then review the issue as soon as possible.

## I Want To Contribute

> ### Legal Notice <!-- omit in toc -->
> When contributing to this project, you must agree that you have authored 100% of the content, that you have the necessary rights to the content and that the content you contribute may be provided under the project license.

### Reporting Bugs

<!-- omit in toc -->
#### Before Submitting a Bug Report

A good bug report shouldn't leave others needing to chase you up for more information. Therefore, we ask you to investigate carefully, collect information and describe the issue in detail in your report. Please complete the following steps in advance to help us fix any potential bug as fast as possible.

**If you are having issues with the web tool:**
- Please include in your bug report:
  - The browser and version you are using
  - If relevant, provide a copy of the FASTA or VCF file you are providing to ShUShER (only for samples that do not contain PHI). Otherwise, provide as much detail about the format and contents as you are able.
  - What went wrong? Include any errors or status messages you can see.
  - There may be additional information in the developer console. [Open the console](https://balsamiq.com/support/faqs/browserconsole/#mozilla-firefox) and copy the output into your bug report.

**If you are having issues running ShUShER locally:**
- Please include in your bug report:
  - OS, Platform and Version
  - Version of npm, Node.js, UShER, etc. depending on what seems relevant.
  - Modifications, if any, you have made to the source code.
  - Relevant errors and information messages.
 
<!-- omit in toc -->
#### How Do I Submit a Good Bug Report?

> You must never report security related issues, vulnerabilities or bugs to the issue tracker, or elsewhere in public. Instead sensitive bugs must be sent by email to [alex.kramer@ucsc.edu](mailto:alex.kramer@ucsc.edu).
<!-- You may add a PGP key to allow the messages to be sent encrypted as well. -->

We use GitHub issues to track bugs and errors. If you run into an issue with the project:

- Open an [Issue](https://github.com/amkram/shusher/issues/new). (Since we can't be sure at this point whether it is a bug or not, we ask you not to label the issue as a bug yet.)
- Explain the behavior you would expect and the actual behavior.
- Please provide as much context as possible and describe the *reproduction steps* that someone else can follow to recreate the issue on their own. This usually includes your code. For good bug reports you should isolate the problem and create a reduced test case.
- Provide the information you collected in the previous section.

Once it's filed:

- The project team will label the issue accordingly.
- We will try to reproduce the issue with your provided steps. If there are no reproduction steps or no obvious way to reproduce the issue, the team will ask you for those steps and mark the issue as `needs-repro`. Bugs with the `needs-repro` tag will not be addressed until they are reproduced.
- If the team is able to reproduce the issue, it will be marked `needs-fix`, as well as possibly other tags (such as `critical`).

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for ShUShER, **including completely new features and minor improvements to existing functionality**. Following these guidelines will help maintainers and the community to understand your suggestion and find related suggestions.

<!-- omit in toc -->
#### Before Submitting an Enhancement

- Read the [documentation](https://github.com/amkram/shusher/README.md) carefully and find out if the functionality is already covered.
- Perform a [search](https://github.com/amkram/shusher/issues) to see if the enhancement has already been suggested. If it has, add a comment to the existing issue instead of opening a new one.
- Find out whether your idea fits with the scope and aims of the project. Good enhancement suggestions are useful to the majority of our users and not just a small subset.

<!-- omit in toc -->
#### How Do I Submit a Good Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://github.com/amkram/shusher/issues).

- Use a **clear and descriptive title** for the issue to identify the suggestion.
- Provide a **step-by-step description of the suggested enhancement** in as many details as possible.
- **Describe the current behavior** and **explain which behavior you expected to see instead** and why. At this point you can also tell which alternatives do not work for you.
- You may want to **include screenshots and animated GIFs** which help you demonstrate the steps or point out the part which the suggestion is related to. You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast) or [this tool](https://github.com/GNOME/byzanz) on Linux. <!-- this should only be included if the project has a GUI -->
- **Explain why this enhancement would be useful** to most ShUShER users. You may also want to point out other projects that solved it better which could serve as inspiration.

<!-- omit in toc -->
## Attribution
This guide is based on the **contributing-gen**. [Make your own](https://github.com/bttger/contributing-gen)!
# UShER WebAssembly Port

This directory contains modifications and build scripts to compile UShER to WebAssembly. All of the below steps are performed in `installUbuntuWeb.sh`.

Most files in this directory are auto-updated upon changes to the main UShER code base. `installUbuntuWeb.sh` the steps necessary to compile UShER to WebAssembly.

# Source code changes
- line `path = boost::filesystem::canonical(outdir);` removed from `usher.cpp`

# Library compilation
UShER has three dependencies: [Protocol Buffers](https://developers.google.com/protocol-buffers), [Boost](https://www.boost.org/), and [oneTBB](https://github.com/oneapi-src/oneTBB). The libraries must be compiled to WebAssembly prior to linking with the UShER port. Details below.

## Boost
ShUShER uses Boost 1.76.0. Boost's build system `b2` ships with configuration support for Emscripten, yielding WebAssembly output. There is a bug in the build system that fails to recognize the correct generator for Emscripten. It is fixed by appending
    
    import generators ;
    generators.override emscripten.searched-lib-generator : searched-lib-generator ;
to the `tools/build/src/tools/emscripten.jam` file in the Boost source subdirectory.

ShUShER requires Boost's `filesystem`, `program_options`, `iostreams`, and `zlib` to compile. The following command compiles these libraries to WebAssembly:

    ./b2 -a toolset=emscripten link=static threading=multi --with-filesystem --with-program_options --with-iostreams -sZLIB_SOURCE=[zlib path] cflags="-DHAVE_UNISTD_H" 
 
 
The output of the above command is `bc` files for each library. These must be converted to `.a` files by running `emar rc` on each. `emar` is the Emscripten version of the GNU `ar` program. The resulting `.a`. files can be linked during compilation.

## Protocol Buffers
`installUbuntuWeb.sh` downloads the version of the `protoc` compiler matching that installed on the system (`protoc-compiler` must be installed before running the script).

Protocol Buffers doesn't require modifications to the source code before compiling. It is compiled to WebAssembly by calling:
       
    ./autogen.sh && emconfigure ./configure --disable-shared --enable-static --build=wasm32 --target=wasm32
    emmake make -j8
in the source directory.

## oneTBB
oneTBB requires a fair amount of source code modifications to successfully compile to WebAssembly. This modified version is hosted in a [separate repository
](https://github.com/amkram/oneTBB-2019-wasm) and pulled during the build process. At the present, multithreading does not work with ShUShER, but oneTBB is still included to avoid modifying the UShER source code.

# Compiling to WebAssembly
After the above steps are completed, the libraries are linked and UShER is compiled to WebAssembly with `emmake` and `make` (see `installUbuntuWeb.sh`).
# Transposed VCF
Transposed VCF supplement MAT with information about ambiguous bases at sample nodes to aid tree optimization. All variants of a sample is written contiguously (in rows instead of columns), so new samples can be concatenated.
## Converting VCF to Transposed VCF
```
transpose_vcf -v <input vcf> -o <output> -T <number of threads>
```
Output will be concatenated if exists. Currently, it can only utilize less than 5 threads, as it is throttled by decompression.
## Converting Transposed VCF to VCF
```
transposed_vcf_to_vcf -i <transposed vcf> -o <vcf output> -r <reference fasta file> -T <number of threads>
```
## Use transposed VCF for optimization
```
matOptimize -i <usher protobuf> -V <transposed VCF> -o <output usher protobuf> -r <radius (4-6) > -T <number of threads>
```
## File Format
<table>
  <tbody>
    <tr>
        <td>Length of compressed Block (8 bytes)</td>
    </tr>
    <tr>
        <td>
        Compressed Data (zlib)
        <table>
            <tbody>
                <tr>
                    <td> Sample name (null terminated) </td>
                </tr>
                <tr>
                    <td> Called mutations (null terminated)
                        <table>
                        <tbody>
                            <tr> 
                              <td> Position of mutation 1 (variant encoded) 
                            <tr> 
                              <td> Position of mutation 2 (variant encoded, omitted if not present) 
                            <tr> 
                              <td> {allele at mutation 2 (one hot , 4 bits),allele at mutation 1 (one hot , 4 bits)} 
                        </tbody>
                        </table>
                    </td>
                </tr>
                <tr>
                    <td>
                        Completely ambiguous mutations (Ns) (Null terminated)
                        <table>
                        <tbody>
                        <tr>
                            <td> End of range (inclusive, variant encoded) 
                        <tr>
                            <td> Strat of range (inclusive, variant encoded, omitted if the same as end of range)
                        </tbody>
                        </table>
                    </td>
                </tr>
            </tbody>
        </table>
    </tr>
  </tbody>
</table>
Contig name is not recorded, as coronavirus only have one contig. If we need to handle multiple contig, concatnate all contig as a long reference, and let position be the position in that concatnated long contig. No change to file format is necessary, just transpose_vcf will need to take faidx to map contig+position to posiition in the concatnated long reference.
# matUtils
matUtils is a toolkit for the manipulation and analysis of mutation-annotated trees (MAT). Please refer to our [wiki](https://usher-wiki.readthedocs.io/en/latest/matUtils.html) for more details. 
## Description

Please summarize any changes this pull request makes.

Fixes {ISSUE} (Please reference any relevant Issues this PR addresses using the reference issue button in the editor.)


## Type of change

- [ ] Bug fix
- [ ] Enhancement

## Tests

Please describe any tests you performed that verify the incoming changes.

## Label PR

Please label this PR with any relevant labels from the right-hand sidebar. (Delete the "Label PR" heading and section)
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
- Steps to reproduce the behavior.
- If relevant, please **provide a copy of the FASTA or VCF file** you are providing to ShUShER (only for samples that do not contain PHI). Otherwise, provide as much detail about the format and contents as you are able.

**Expected behavior**

A clear and concise description of what you expected to happen.

**Screenshots**

If applicable, add screenshots to help explain your problem.

**System Details (please complete the following information):**

 - OS: [e.g. Ubuntu 20.04]
 - Browser [e.g. Firefox 91.0.1]. Please include version if possible.
 - Version of npm, Node.js, UShER, etc. depending on what seems relevant.

**Additional context**

- If this bug report is related to usage of the web app, **please include the following additional information**:
  - Browser developer console output. Please follow the instructions [here](https://balsamiq.com/support/faqs/browserconsole/) to open the console and paste the output into the bug report.
- Modifications, if any, you have made to the source code.
- Relevant errors and information messages.
- Add any other context about the problem here.
This directory contains the React web app that runs the ported version of UShER (`../usher-port`). See the main [README](https://github.com/amkram/shusher/blob/master/README.md) for information on installation and usage.
