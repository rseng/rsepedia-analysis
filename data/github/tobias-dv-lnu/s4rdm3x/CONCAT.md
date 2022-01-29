![Core Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/Core%20Build/badge.svg)
![Core Tests](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/Core%20Tests/badge.svg)
[![](https://jitpack.io/v/tobias-dv-lnu/s4rdm3x.svg)](https://jitpack.io/#tobias-dv-lnu/s4rdm3x)

![v3xt Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/v3xt%20Build/badge.svg)
![CmdExRunner Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/CmdExRunner%20Build/badge.svg)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/tobias-dv-lnu/s4rdm3x/blob/master/LICENSE)

[![status](https://joss.theoj.org/papers/f4301adc3e9121a10354c355d91b5c1f/status.svg)](https://joss.theoj.org/papers/f4301adc3e9121a10354c355d91b5c1f)


# s4rdm3x
A tool suite to perform experiments in automatic mapping of source code to modular architecture definitions, also called the orphan adoption problem. It consists of a reusable base code (core) and two tools (v3xt & CMDExRunner).

## core
The base code provides Java bytecode analysis to extract a dependency graph (and naming information) as well as loading an architectural definition and source to module mapping. Furthermore it implements the HuGMe method and four attraction functions to map a source code file to an architectural module. The attraction functions are CountAttract, IRAttract, LSIAttract and NBAttract.

## v3xt
A tool that provides a GUI to define and run small scale experiments as well as visualize the results in real-time. This can be used to quickly try and assess new ideas and define larger experiments. Supports loading and saving of experiments definitions as experiments.

## CMDExRunner
A command line tool for executing experiments in parallel. It reads an experiment definition xml-file and distributes the experiments over a number of threads. Typically useful for running experiments in multicore computing clouds.

# Documentation
Documentation is available in [docs](docs) and published: https://tobias-dv-lnu.github.io/s4rdm3x/

# Licence
s4rdmex, v3xt, cmdexrunner
Copyright (c) 2020 Tobias Olsson

Released under

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
https://www.gnu.org/licenses/gpl-3.0.html

See LICENCE for further details



---
title: 's4rdm3x: A Tool Suite to Explore Code to Architecture Mapping Techniques'
tags:
  - Java
  - orphan adoption
  - clustering
  - reflexion modeling
  - architecture conformance checking
authors:
  - name: Tobias Olsson
    orcid: 0000-0003-1154-5308
    affiliation: 1
  - name: Morgan Ericsson
    orcid: 0000-0003-1173-5187
    affiliation: 1
  - name: Anna Wingkvist
    orcid: 0000-0002-0835-823X
    affiliation: 1
affiliations:
 - name: Department of Computer Science and Media Technology, Linnaeus University, Sweden
   index: 1
date: 13 April 2020
bibliography: paper.bib

---

# Summary

Architectural drift and erosion, where the implementation starts to deviate from the intended software architecture or the rules set by it, are common problems in long-lived software systems. This can be avoided by using techniques, such as Reflexion modeling [@murphy1995software], to validate that the implementation conforms to the indented architecture. Unfortunately, such techniques require a mapping from source code modules (e.g., classes) to elements of the architecture, something that is not always available or up to date. This is a known problem when, e.g., companies want to adopt static architecture conformance checking; the effort to manually create or bring this mapping up to date is just too time-consuming and error-prone [@Ali2017ArchitectureRequirements; @InfoRetrieval].

The ``s4rdm3x`` tool suite is designed for researchers and practitioners to study and evaluate algorithms that perform part of the mapping automatically, such as orphan-adoption clustering [@HuGMe] or information retrieval techniques [@InfoRetrieval]. It includes a graphical user interface to define and run experiments with mapping algorithms and their parameters, and visualize and explore the results of these experiments. The experiments can be executed locally or in a remote high-performance computing environment. The tool suite includes reference implementations of state of the art mapping algorithms and a set of Java systems with validated mappings between classes and architecture elements. The tool suite is extensible, so it is easy to add new mapping algorithms and visualizations to explore their performance.

# Statement of Need
To faciliate the further development and evaluiation of mapping techniques the software provides reference implementations of the current state-of-the-art mapping techniques and the means to implement new techniques and run experiments. It includes the HuGMe orphan adoption clustering method [@HuGMe], and four attraction functions to decide which architectural element a source code module should be mapped to: `CountAttract` [@HuGMe], `IRAttract`, `LSIAttract` [@InfoRetrieval] and `NBAttract` [@NaiveBayes]. There is also a reference implementation of our novel technique to create a textual representation of source code dependencies at an architectural level; Concrete Dependency Abstraction (CDA). It also contains a set of validated mappings between source code classes and architectural elements that are often used in software architecture erosion research. These systems have either been recovered from replication packages [@brunet2012evolutionary; @LenhardExploringSCMIndicatingArchInconsistency] or the [SAEroCon workshop repository](https://github.com/sebastianherold/SAEroConRepo). 

# The ``s4rdm3x`` Tool Suite

``S4rdm3x`` is an extensible suite of tools for source code analysis, architecture definition, mapping of source code modules to architecture elements, experiment definitions, and exploratory and visual analysis. The suite consists of an *extensible core* and two tools, a *graphical editor* to create and visualize mapping experiments and a *command-line tool to run experiment* at scale. 

![Overview of the s4rdm3x Core showing the implementation metamodel, and the mapper and experiment subsystem](classes.pdf)

The *core* provides Java bytecode analysis to extract a dependency graph (and naming information) as well as loading an architectural definition and source to module mapping. Figure 1 provides an overview of the important classes in the core and their main dependencies; the three leftmost classes represent source code modules, and their implemented dependencies contained as a graph and the architecture is represented by components and their allowed dependencies. There is a rich set of dependencies that are extracted from (Java) byte code, including the possibility to include implicit dependencies found via hard-coded constants. The `MapperBase` class is used to implement different mapping strategies as subclasses. `MapperExperiment` provides functionality to set up and run mapping experiments using combinations of random parameters at different intervals. An experiment is implemented as a subclass of `MapperExperiment` that instantiates the corresponding subclasses of `MapperBase`. This means that the mappers are not exclusive for experiments and can be reused in other situations, e.g., in a tool that performs semi-automatic mapping as part of a reflexion modeling approach. 

The *graphical editor* is used to define, visualize, analyze, and compare how well mapping algorithms perform with different parameters and initial sets of known mappings. It supports a range of visualizations and can be extended with new ones. The editor uses an Immediate Mode GUI approach, where the application renders the graphical primitives it needs (e.g., lines, rectangles, and points) every frame, an approach often used in computer games and tools used for computer game development since it offers fine-grained control over the visualization. This fine-grained control makes it possible to extend the editor with custom visualizations. The GUI uses [OpenGL](https://opengl.org) to provide hardware-accelerated rendering. 

The graphical editor can be used to run experiments, but these generally require a large number of combinations of, e.g., parameters, initial sets, and systems, so they can take a long time to run. The suite includes a command-line tool that runs these combinations in parallel on many-core machines. The command-line tool can read experiment definitions in XML exported from the graphical editor and save the results in a format that can be imported and visualized.
 


``S4rdm3x`` is implemented in Java and depends on [ASM](https://asm.ow2.io), [Weka](https://www.cs.waikato.ac.nz/ml/weka), and [Dear JVM ImGui](https://github.com/kotlin-graphics/imgui).

# Applications

The ``S4rdm3x`` tool suite has been used in research studies on orphan adoption [@ImprovedHuGMe; @NaiveBayes] and as a continuous integration tool-chain for static architecture conformance checking of student project submissions. 

# Acknowledgments

This work is supported by the [Linnaeus University Centre for Data Intensive Sciences and Applications (DISA)](https://lnu.se/forskning/sok-forskning/linnaeus-university-centre-for-data-intensive-sciences-and-applications) High-Performance Computing Center.  

# References
![Core Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/Core%20Build/badge.svg)
![Core Tests](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/Core%20Tests/badge.svg)
[![](https://jitpack.io/v/tobias-dv-lnu/s4rdm3x.svg)](https://jitpack.io/#tobias-dv-lnu/s4rdm3x)

![v3xt Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/v3xt%20Build/badge.svg)
![CmdExRunner Build](https://github.com/tobias-dv-lnu/s4rdm3x/workflows/CmdExRunner%20Build/badge.svg)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/tobias-dv-lnu/s4rdm3x/blob/master/LICENSE)

A tool suite to perform experiments in automatic mapping of source code to modular architecure definitions, also called the orphan adoption problem. It consists of a reusable base code (core) and two tools (v3xt & CMDExRunner).

## core
The base code provides Java bytecode analysis to extract a dependency graph (and naming information) as well as loading an architectural definition and source to module mapping. Furthermore it implements the HuGMe method and four attraction functions to map a source code file to an architectural module. The attraction functions are CountAttract, IRAttract, LSIAttract and NBAttract.

## v3xt
A tool that provides a GUI to define and run small scale experiments as well as visualize the results in real-time. This can be used to quickly try and asses new ideas and define larger experiments. Supports loading and saving of experiments definitions as experiments.

## CMDExRunner
A command line tool for executing experiments in parrallell. It reads an experiment definition xml-file and distributes the experiments over a number of threads. Typically useful for running experiments in multicore computing clouds.

# Installation
Download the latest [Release](https://github.com/tobias-dv-lnu/s4rdm3x/releases) from GitHub, it includes precompiled versions of the tools as jar applications for all platforms.

## Data Systems
S4rdm3x is distributed with architectural models and source code to implementation mappings for a number of systems (see the data/systems directory). Models and mappings are based on work during the [SAEroCon workshop](https://saerocon.wordpress.com/) and a replication package provided by Joao Brunet et. al. [On the Evolutionary Nature of Architectural Violations](https://code.google.com/archive/p/on-the-nature-dataset/wikis/ReplicabilityOfTheStudy.wiki)

The .sysmdl files are text files and should be viewable/editable in any text editor

To complement the models and mappings, the actual compiled systems are also needed (i.e. jar files). These are not included in this distribution as this would create a problematic licensing situation. However, the jar file dependencies are documented in the respective model file and should be available by either looking through the links from the sources above, or finding them in the actual official distribution of the systems.

You place the jar files in the same directory as the corresponding sysmdl file.

### Some Official Distributions
* [JabRef-3.7.jar](https://github.com/JabRef/jabref/releases/tag/v3.7)
* [ProM 6.9](http://www.promtools.org/doku.php?id=prom69)

# Prerequisites

- Java 11 or superior

# Running
Run via commandline:

`java -jar v3xt.jar`

`java -jar cmdexrunner.jar`

For OSX users the `-XstartOnFirstThread` JVM option needs to be supplied when running the `v3xt.jar`. Also note that the OSX version is highly unstable, you may need to try to start it several times.

`java -XstartOnFirstThread -jar v3xt.jar`

Check the included readme for further details.




## How To Use It
The release contains a ready to use test experiment as explained in the release README

# Contributing
Please report any bugs or anomalies in the github repository. 

# Development
* [Learn about the needed dependencies](dependencies "Dependencies")
* [How to set up a development environment](devenv "DevEnv")
* [How to add a new mapper](add_new_mapper)
* API [JavaDoc](https://tobias-dv-lnu.github.io/s4rdm3x/api "JavaDoc")


# Adding a new type of Mapper
This documentation describes the steps of adding a new type of mapper. The job of a mapper is to assign an orphan to an architectural module. To aid in this work there are a number of non-orphans that are already mapped to the modules.

## Experiment Subsystem
The most important aspect to understand is the experiment sub-system and how it works. The following UML diagram gives an overview of the current design. The starting class is the [ExperimentRunner]. This is where systems are loaded, experiments executed and results saved.
The basic steps to implement a new mapper would be to
1. Create the mapper, included examples are [HuGMe], [NBAttractMapper], [IRAttractMapper], and [LSIAttractMapper].
2. Create the mapper specific [ExperimentRun] subclass. If the mapper is based on information retrieval you can consider inheriting from [IRExperimentRunBase].
3. Create the specific data [ExperimentRunData] subclass and modify [RundDataCSVFileSaver] to save the data in csv format.
4. Modify [ExperimentXMLPersistence] to save mapper and mapper parameter settings to XML

![alt text](img/ex_ss.png "Experiment Subsystem")

## GUI
To be able to use a new mapper in the graphical user interface a few modifications is needed.
1. Add the new constant for the mapper in [ExperimentRunnerViewThread]
2. Add the parameters of the new mapper to [MapperView] and implement the visualization in [MapperView.doExperiment()] and mapper creation in [MapperView.createExperiment()]


[ExperimentRunner]:api/se/lnu/siq/s4rdm3x/experiments/ExperimentRunner.html
[ExperimentRun]:api/se/lnu/siq/s4rdm3x/experiments/ExperimentRun.html
[IRExperimentRunBase]:api/se/lnu/siq/s4rdm3x/experiments/IRExperimentRunBase.html
[ExperimentRunData]:api/se/lnu/siq/s4rdm3x/experiments/ExperimentRunData.html
[RundDataCSVFileSaver]:api/se/lnu/siq/s4rdm3x/experiments/RundDataCSVFileSaver.html
[ExperimentXMLPersistence]:api/se/lnu/siq/s4rdm3x/experiments/ExperimentXMLPersistence.html
[HuGMe]:api/se/lnu/siq/s4rdm3x/model/cmd/mapper/HuGMe.html
[NBAttractMapper]:api/se/lnu/siq/s4rdm3x/model/cmd/mapper/NBAttractMapper.html
[IRAttractMapper]:api/se/lnu/siq/s4rdm3x/model/cmd/mapper/IRAttractMapper.html
[LSIAttractMapper]:api/se/lnu/siq/s4rdm3x/model/cmd/mapper/LSIAttractMapper.html

[ExperimentRunnerViewThread]:api/experimenting/ExperimentRunnerViewThread.html
[MapperView]:api/experimenting/MapperView.html
[MapperView.doExperiment()]:api/experimenting/MapperView.html#doExperiment(gui.ImGuiWrapper,boolean)
[MapperView.createExperiment()]:/api/experimenting/MapperView.html#createExperiment()

# Development Environment
Core, v3xt and CmdExRunner are all distributed as independent gradle projects. Sources include gradlew scripts to build each projects. In addition the [core is distributed via jitpack](https://jitpack.io/#tobias-dv-lnu/s4rdm3x) for easy inclusion in new project, see cmdexrunnder/build.gradle, for an example of how to use it.

One caveat is that currenly asm, weka, snowball and bounce are distributed using copied files available in the lib folder. These must be specifically included into the gradle build script (see cmdexrunnder/build.gradle) as runtime libraries.
# Dependencies

## Compile-Time
## WEKA 3.8.3
available in lib

## ASM 6.2.1
available in lib

## Dear JVM IMGui IMGUI (v3xt only)
https://github.com/kotlin-graphics/imgui available via v3xt/gradle.build

## Lightweight Java Game Library (v3xt only)
https://lwjgl.org available via v3xt/gradle.build

## Test Dependencies
### JUnit 5
Automatic tests use JUnit5: https://junit.org/junit5

## Run-Time
### Bounce 0.18
available in lib

### Snowball
available in lib
