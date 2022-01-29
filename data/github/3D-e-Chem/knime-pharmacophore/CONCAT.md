# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [1.0.4] 2019-09-24

### Fixed

- Pharmacophore from molecule node parsed SDF incorrectly ([#12](https://github.com/3D-e-Chem/knime-pharmacophore/issues/12))

## [1.0.3] 2019-06-27

### Changed

- Compatible with KNIME 4 [#10](https://github.com/3D-e-Chem/knime-pharmacophore/issues/10)

## [1.0.2] 2018-07-05

### Fixed

- PharmacophorePoint#toArray method uses current locale [#9]

## [1.0.1] 2017-11-21

### Fixed

- Transformation matrix performs mirroring [#8]

## [1.0.0] 2017-09-26

### Added

- Node extract points from pharmacophore and node that does the reverse [#2]
- Nodes to read/write phar formatted files
- Nodes to convert between pharmacophore and molecule (sdf/mol) [#3]
- Node to align one pharmacophore to another [#4]

### Fixed

- RMSD very high for a good fit [#6]
[KNIME](https://www.knime.com) plugin with nodes to convert and align pharmacophores.

[![Travis-CI Build Status](https://travis-ci.org/3D-e-Chem/knime-pharmacophore.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-pharmacophore)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/0d27c4nhkjopy69r/branch/master?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/knime-pharmacohore)
[![SonarCloud Gate](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.knime.pharmacophore%3Anl.esciencecenter.e3dchem.knime.pharmacophore&metric=alert_status)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.knime.pharmacophore:nl.esciencecenter.e3dchem.knime.pharmacophore)
[![SonarCloud Coverage](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.knime.pharmacophore%3Anl.esciencecenter.e3dchem.knime.pharmacophore&metric=coverage)](https://sonarcloud.io/component_measures/domain/Coverage?id=nl.esciencecenter.e3dchem.knime.pharmacophore:nl.esciencecenter.e3dchem.knime.pharmacophore)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.997332.svg)](https://doi.org/10.5281/zenodo.997332)

A pharmacophore is an abstract description of molecular features that are necessary for molecular recognition of a ligand by a biological macromolecule.
Nodes in this plugin allow for converting pharmacophores, from and to molecules, by mapping elements to pharmacophore type and, reading from or writing to the phar file format used by the [Silicos IT align-it tool](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html#pharmacophores).
This plugin adds the Pharmacophore (Phar) data type to KNIME, allowing nodes to read/write/manipulate pharmacophores inside KNIME like the [Silicos-it align-it](https://github.com/3D-e-Chem/knime-silicos-it), [Kripo pharmacophore retrieval](https://github.com/3D-e-Chem/knime-kripodb) and [molviewer pharmacophore viewer](https://github.com/3D-e-Chem/knime-molviewer) nodes.

This project uses [Eclipse Tycho](https://www.eclipse.org/tycho/) to perform build steps.

# Installation

Requirements:

* KNIME, https://www.knime.org, version 4.0 or higher

Steps to get the Pharmacophore KNIME nodes inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/updates`
4. Select --all sites-- in `work with` pulldown
5. Select the `Pharmacophore KNIME nodes`
6. Install software
7. Restart KNIME

# Usage

1. Create a new KNIME workflow.
2. Find node in Node navigator panel.
3. Drag node to workflow canvas.

# Examples

The `examples/` folder contains example KNIME workflows.

# Build

```
mvn verify
```

An Eclipse update site will be made in `p2/target/repository` directory.
The update site can be used to perform a local installation.

## Continuous Integration

Configuration files to run Continuous Integration builds on Linux (Travis-CI), OS X (Travis-CI) and Windows (AppVeyor) are present.

See `./.travis.yml` file how to trigger a Travis-CI build for every push or pull request.
Also see `./.travis.yml` file how to upload coverage to https://www.codacy.com .

See `./appveyor.yml` file how to run on https://www.appveyor.com .

# Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.pharmacophore.targetplatform/KNIME-AP-4.0.target` target definition.
6. A KNIME Analytics Platform instance can be started by right clicking on the `targetplatform/KNIME\ Analytics\ Platform.launch` file and selecting `Run As → KNIME Analytics Platform`. The KNIME instance will contain the target platform together with all extensions defined in the workspace.

During import the Tycho Eclipse providers must be installed.

## Tests

Tests for the node are in `tests/src` directory.
Tests can be executed with `mvn verify`, they will be run in a separate KNIME environment.
Test results will be written to `test/target/surefire-reports` directory.
Code coverage reports (html+xml) can be found in the `tests/target/jacoco/report/` directory.

### Unit tests

Unit tests written in Junit4 format can be put in `tests/src/java`.

### Workflow tests

See https://github.com/3D-e-Chem/knime-testflow#3-add-test-workflow

## Speed up builds

Running mvn commands can take a long time as Tycho fetches indices of all p2 update sites.
This can be skipped by running maven offline using `mvn -o`.

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Create package with `mvn package`, will create update site in `p2/target/repository`
3. Run tests with `mvn verify`
4. Optionally, test node by installing it in KNIME from a local update site
5. Append new release to an update site
  1. Make clone of an update site repo
  2. Append release to the update site with `mvn install -Dtarget.update.site=<path to update site>`
6. Commit and push changes in this repo and update site repo.
7. Create a Github release
8. Update Zenodo entry
  1. Correct authors
9. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
