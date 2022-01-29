# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]

Workflows using v1 molviewer nodes will need to replace them with the v2 molviewer nodes.

## [2.0.2] - 2019-11-07

### Fixed

- Unable to show second pharmacophore ([#31](https://github.com/3D-e-Chem/knime-molviewer/issues/31))

## [2.0.1] - 2019-09-25

### Fixed

- Pharmacphore viewer gave syntax error

## [2.0.0] - 2019-09-24

### Changed

- Replaced Java based nodes via server with Dynamic js nodes ([#27](https://github.com/3D-e-Chem/knime-molviewer/issues/27))

## [1.1.3] - 2019-06-27

### Changed

- Compatible with KNIME 4 ([#25](https://github.com/3D-e-Chem/knime-molviewer/issues/25))

## [1.1.2] - 2017-12-11

### Added

- Centering on current scene (#21)

### Fixed

- Handle missing values (#19)

## [1.1.1] - 2017-04-10

### Fixed

- Executed node contains no data (#24)

## [1.1.0] - 2017-09-08

### Added

- Node to view pharmacophores with optional protein/ligand (#11)
- Node to view only ligands (#20)
- Node to view only proteins (#22)
- Color picker for ligands (#21)

## [1.0.2] - 2017-02-27

### Added

- Pocket as ball+stick toggle
- Added slider for adjusting pocket selection radius

## [1.0.1] - 2017-01-25

### Changed

* Use xdg-open on Linux to open url on web-browser (#16)

### Removed

* KNIME crash workaround in node description

## [1.0.0] - 2017-01-23

### Added

* KNIME crash workaround in node description (#16)

### Changed

* If KNIME is unable to open url in web browser then give a warning

## [0.1.5] - 2016-12-02

### Changed

* Allow proteins to be in mol2 format (#14)
* Switched from 3Dmol to NGL (https://github.com/3D-e-Chem/molviewer-tsx/issues/4)

## [0.1.4] - 2016-11-25

### Added

* Example workflow (#2)
* Test workflow (#4)
* Auto open web browser checkbox (#5)

### Changed

* Replace default node icon (#3)

### Fixed

* UI hosted as / (#1)
* Error if protein input is from Vernalis PDB Downloader node (#12) 
KNIME node which launches a web browser with a molecule viewer powered by [NGL](https://github.com/arose/ngl).

[![Build Status Linux](https://travis-ci.org/3D-e-Chem/knime-molviewer.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-molviewer)
[![Build status Windows](https://ci.appveyor.com/api/projects/status/595y9gh13d69y61q?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/knime-molviewer)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.knime.molviewer%3Anl.esciencecenter.e3dchem.knime.molviewer&metric=alert_status)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.knime.molviewer%3Anl.esciencecenter.e3dchem.knime.molviewer)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.knime.molviewer%3Anl.esciencecenter.e3dchem.knime.molviewer&metric=coverage)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.knime.molviewer%3Anl.esciencecenter.e3dchem.knime.molviewer)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597231.svg)](https://doi.org/10.5281/zenodo.597231)

If you are using KNIME workflows with large molecules and you want to view them in 3D then the molviewer nodes are able to handle this.
* Provides cheminformatics trying to model proteins within KNIME a way to view them
* Adds nodes to KNIME to visualize proteins, small molecules and pharmacophores
* Adds support to KNIME to  handle viewing big molecules
* Used to compare pharmacophore generation tools
* Used to check ligand repurposing

This project uses a web user interface based on https://github.com/3D-e-Chem/molviewer-tsx .

![From KNIME launch web based molviewer](https://raw.githubusercontent.com/3D-e-Chem/knime-molviewer/master/docs/molviewer-composite.png)

# Installation

Requirements:

* KNIME, https://www.knime.org, version 4.0 or higher

Steps to get the MolViewer KNIME node inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/updates`
4. Select --all sites-- in `work with` pulldown
5. Select the node called `MolViewer nodes for KNIME`
6. Install software
7. Restart KNIME

# Usage

1. Create a new KNIME workflow.
2. Find node in Node navigator panel under Community Nodes > 3D-e-Chem > Molviewer.
3. Drag node to workflow canvas.
4. Give molviewer nodes Proteins, Ligands and/or Pharmacophores as input. 
5. Open the view of the molviewer node, this will launch a web browser with the molecules visualized.
6. In web browser you can

   * Rotate/translate/zoom molecules with mouse
   * Toggle which molecules are visible
   * Transmit visible molecules as HiLite selection back to KNIME or do the reverse from KNIME to molviewer
   * and more

See example workflow in `examples/` directory.

# Build

```
mvn verify
```

An Eclipse update site will be made in `p2/target/repository` repository.
The update site can be used to perform a local installation.

# Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.molviewer.targetplatform/KNIME-AP-4.0.target` target definition.
6. A KNIME Analytics Platform instance can be started by right clicking on the `targetplatform/KNIME\ Analytics\ Platform.launch` file and selecting `Run As → KNIME Analytics Platform`. The KNIME instance will contain the target platform together with all extensions defined in the workspace.

During import the Tycho Eclipse providers must be installed.

## Web interface

The web interface in the `plugin/src/main/java/js-lib/molviewer` directory. Is a distribution from the https://github.com/3D-e-Chem/molviewer-tsx repository.

## Tests

Tests for the node are in `tests/src` directory.
Tests can be executed with `mvn verify`, they will be run in a separate KNIME environment.
Test results will be written to `test/target/surefire-reports` directory.

### Workflow tests

See https://github.com/3D-e-Chem/knime-testflow#3-add-test-workflow

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Test node by installing it from local update site
5. Append new release to an update site
  1. Make clone of an update site repo
  2. Append release to the update site with `mvn install -Dtarget.update.site=<path to update site>`
6. Add files, commit and push changes of update site repo.
7. Create a Github release
8. Update Zenodo entry, fix authors
9. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
10. Update DOI in CITATION.cff

# Technical architecture

The nodes are implemented as Dynamic JavaScript nodes also known as a org.knime.dynamic.node.generation.dynamicNodes extension.

The viewer has two way binding for the molecule visibility.

This project uses [Eclipse Tycho](https://www.eclipse.org/tycho/) to perform build steps.
