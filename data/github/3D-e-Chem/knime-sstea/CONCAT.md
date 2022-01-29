# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/

## [Unreleased]

## [1.0.6] - 2019-06-26

### Changes

- Requires KNIME 4.0 [#12](https://github.com/3D-e-Chem/knime-sstea/issues/12)

## [1.0.5] - 2016-07-18

### Added

- Integration tests, by running workflows

### Changed

- Nested Sygma node under 3D-e-chem category (#10)

### Removed

- Node views

## [1.0.4] - 2016-06-03

### Changed

- Improved node documentation (#3)

## [1.0.3] - 2015-05-19

### Fixed

- Throw exception when sequence is empty

## [1.0.2] - 2015-05-19

### Fixed

- exception when sequences are not same length (#2)

## [1.0.1] - 2015-05-04

Initial release

### Added

- Compute ss-TEA score on sequence alignment and subfamily members.
KNIME node to calculate subfamily specific two entropy analysis (ss-TEA) score.

The ss-TEA can identify specific ligand binding residue positions for any receptor, predicated on high quality sequence information.

See reference at https://doi.org/10.1186/1471-2105-12-332 for a description of the score.

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-sstea.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-sstea)
[![DOI](https://zenodo.org/badge/19641/3D-e-Chem/knime-sstea.svg)](https://zenodo.org/badge/latestdoi/19641/3D-e-Chem/knime-sstea)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.sstea%3Anl.esciencecenter.e3dchem.sstea&metric=alert_status)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.sstea%3Anl.esciencecenter.e3dchem.sstea)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.sstea%3Anl.esciencecenter.e3dchem.sstea&metric=coverage)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.sstea%3Anl.esciencecenter.e3dchem.sstea)

# Installation

Requirements:

* KNIME, https://www.knime.org

Steps install ss-TEA KNIME node:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/updates`
4. Select --all sites-- in work with pulldown
5. Open KNIME 3D-e-Chem Contributions folder
6. Select ss-TEA
7. Install software & restart

# Usage

See example workflow at [examples/ss-TEA-example.zip](examples/ss-TEA-example.zip).

It can be run by importing it into KNIME as an archive.

# Build

```
mvn verify
```

Jar has been made in `/target` folder.
An Eclipse update site will be made in `p2/target/repository` repository.

# Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.sstea.targetplatform/KNIME-AP-4.0.target` target definition.

During import the Tycho Eclipse providers must be installed.

## Tests

Tests for the node are in `tests/src` directory.
Tests can be executed with `mvn verify`, they will be run in a separate Knime environment.

### Unit tests

Unit tests written in Junit4 format can be put in `tests/src/java`.

### Workflow tests

See https://github.com/3D-e-Chem/knime-testflow#3-add-test-workflow

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Append new release to 3D-e-Chem update site
  1. Make clone of https://github.com/3D-e-Chem/3D-e-Chem.github.io repo
  2. Append release to 3D-e-Chem update site with `mvn install -Dtarget.update.site=<3D-e-Chem repo/updates>`
5. Commit and push changes in this repo and 3D-e-Chem.github.io repo
6. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
