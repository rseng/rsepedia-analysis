# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [2.0.0] - 2019-07-03

### Added

* CITATION.cff and .zenodo.json files, for better citation support

### Changed

* Requires KNIME version 4 ([#2](https://github.com/3D-e-Chem/knime-python-node-archetype/issues/2))
* Switched from KNIME SDK to Eclipse + target platform
* Source jars for plugin have been replaced with source reference in MANIFEST.MF

## [1.4.0] - 2018-02-07

Synced with [tycho-knime-node-archetype v1.5.0](https://github.com/3D-e-Chem/tycho-knime-node-archetype/releases/tag/v1.5.0).

### Changed

* Upgraded to v2.0.1 of nl.esciencecenter.e3dchem.python.plugin, which allows for running nodes in either Python version 2 or 3

## [1.3.0] - 2017-03-09

Synced with [tycho-knime-node-archetype v1.4.0](https://github.com/3D-e-Chem/tycho-knime-node-archetype/releases/tag/v1.4.0).

## [1.2.2] - 2017-01-24

### Added

* Instruction to add generated repo to 3D-e-Chem KNIME feature

## [1.2.1] - 2016-07-12

### Changed

* Templatify KNIME test workflow
* Process Travis-CI config file by template engine
* Improved development instructions

## [1.2.0] - 2016-07-11

Initial release

Forked from https://github.com/3D-e-Chem/tycho-knime-node-archetype

### Changed

* Refactored KNIME node to a Python Wrapper node (https://github.com/3D-e-Chem/knime-python-wrapper)
# KNIME Python node archetype

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-python-node-archetype.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-python-node-archetype)
[![Build status](https://ci.appveyor.com/api/projects/status/5dory9qjycepcmcn/branch/master?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/knime-python-node-archetype/branch/master)
[![Download](https://api.bintray.com/packages/nlesc/knime-python-node-archetype/knime-python-node-archetype/images/download.svg) ](https://bintray.com/nlesc/knime-python-node-archetype/knime-python-node-archetype/_latestVersion)
[![DOI](https://zenodo.org/badge/63080247.svg)](https://zenodo.org/badge/latestdoi/63080247)

Generates [KNIME](http://www.knime.org) workflow node skeleton repository with sample code.
The node executes a Python script which is included in the skeleton.
The script uses dictionary for dialog options and [Pandas](http://pandas.pydata.org/) DataFrames as input and output.

This archetype was made because the instructions to create KNIME nodes at https://tech.knime.org/developer-guide, requires interaction with Eclipse wizards. We wanted a way to start and perform node development from the command line and headless.
KNIME nodes are Eclipse plugins. The [Tycho](https://eclipse.org/tycho/) Maven plugin is used to build and handle dependencies of Eclipse plugins, so we use Tycho for KNIME node building.

The [Maven archetype](https://maven.apache.org/guides/introduction/introduction-to-archetypes.html) will generate a multi-module project with the following structure:

* / - parent Maven project
* /plugin/ - code for KNIME node
* /tests/ - tests of KNIME node
* /feature/ - eclipse feature
* /p2/ - eclipse update site

The KNIME node will execute a Python script called `/plugin/src/<package>/<python script file name>.py`.

See https://github.com/3D-e-Chem/knime-python-wrapper for more information how the Python Wrapper node works.

## Requirements

* Java ==8
* Maven >=3.0

The archetype is hosted on a [BinTray repository](https://dl.bintray.com/nlesc/knime-python-node-archetyp).
Maven does not resolve to this BinTray repository by default so it must be added.

The ~/.m2/settings.xml should contain the following profile:

```xml
<?xml version="1.0" encoding="UTF-8" ?>
<settings xmlns="http://maven.apache.org/SETTINGS/1.0.0"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://maven.apache.org/SETTINGS/1.0.0
                          https://maven.apache.org/xsd/settings-1.0.0.xsd">
  <profiles>
    <profile>
      <id>pythonknimearchetype</id>
      <repositories>
        <repository>
          <id>python-knime-archetype</id>
          <url>https://dl.bintray.com/nlesc/knime-python-node-archetype</url>
        </repository>
      </repositories>
    </profile>
  </profiles>
</settings>
```

## Generate

The following command will generate a skeleton project

```sh
mvn archetype:generate -DarchetypeGroupId=nl.esciencecenter \
  -DarchetypeArtifactId=knime-python-node-archetype \
  -DarchetypeVersion=2.0.0 \
  -P pythonknimearchetype
```

The command will ask the following questions:

1. Enter the groupId
2. Enter the artifactId
3. Enter the name of the package under which your code will be created
4. Enter the version of your project, use `x.y.z-SNAPSHOT` format (for example `1.2.3-SNAPSHOT`), where x.y.z is [semantic versioning](http://semver.org/).
5. Enter the GitHub organization name or GitHub username
6. Enter the GitHub repository name
7. Enter the KNIME node name
8. Enter the Python script file name (must be given without .py extension)
9. Enter the required Python package name (The presence of this Python package will be checked before executing the node)
10. Confirm

The skeleton has been generated in a sub-directory named after the artifactId in the current working directory.

The following steps are needed to get a ready to use project.

11. Change directory to generated code.
12. Make skeleton git aware, by running `git init`.
13. Fill in all placeholders (`[Enter ... here.]`) in

    * plugin/src/java/**/*.xml
    * feature/feature.xml

14. Commit all changes and push to GitHub
15. Optionally, setup Continuous Integration as described in the project README.md file.

Further instructions about generated project can be found in it's README.md file.

## Generate from inside KNIME SDK

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Register the archetype catalog which contains this archetype

      1. Goto Window > Preferences > Maven > Archetypes
      2. Add a remote catalog with `https://github.com/3D-e-Chem/knime-python-node-archetype/raw/master/archetype-catalog.xml`

5. Create a new project based on archetype

      1. Goto File > New > Project ... > Maven
      2. Select Maven Project & press Next button
      3. Use default location & press Next button
      4. Select Catalog you added in step 4.2
      5. Select the archetype with artifact id `knime-python-node-archetype` & press Next button
      6. Fill in the form & press Finish button

Further instructions about generated project can be found in it's README.md file.

## New release

1. Adjust version in pom.xml
2. Update CHANGELOG.md & README.md & archetype-catalog.xml & CITATION.cff
3. Commit & push
4. Test archetype by running `mvn verify`
5. Create GitHub release
6. Deploy to Bintray, see Deploy chapter below

### Deploy

To deploy current version to Bintray.

1. Add bintray API key to [~/.m2/settings.xml](https://maven.apache.org/settings.html)

```
<servers>
    <server>
        <id>bintray-nlesc-knime-python-node-archetype</id>
        <username>************</username>
        <password>********************************</password>
    </server>
<servers>
```

2. Run `mvn deploy`

## Attribution

The https://github.com/3D-e-Chem/tycho-knime-node-archetype was used as a starting point of this archetype.
# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]
#set( $symbol_pound = '#' )
#set( $symbol_dollar = '$' )
#set( $symbol_escape = '\' )
Skeleton for Python based KNIME nodes with sample code as described [here](https://tech.knime.org/developer-guide).

[![Build Status](https://travis-ci.org/${github_organization}/${github_repository}.svg?branch=master)](https://travis-ci.org/${github_organization}/${github_repository})
[![SonarCloud Gate](https://sonarcloud.io/api/badges/gate?key=${groupId}:${artifactId})](https://sonarcloud.io/dashboard?id=${groupId}:${artifactId})
[![SonarCloud Coverage](https://sonarcloud.io/api/badges/measure?key=${groupId}:${artifactId}&metric=coverage)](https://sonarcloud.io/component_measures/domain/Coverage?id=${groupId}:${artifactId})

This project uses [Eclipse Tycho](https://www.eclipse.org/tycho/) to perform build steps.

${symbol_pound} Installation

Requirements:

* KNIME, https://www.knime.org, version ${knime_version} or higher

Steps to get the ${node} KNIME node inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with url of update site which contains this node.
4. Select --all sites-- in `work with` pulldown
5. Select the node
6. Install software
7. Restart KNIME

${symbol_pound} Usage

1. Create a new KNIME workflow.
2. Find node in Node navigator panel.
3. Drag node to workflow canvas.

${symbol_pound} Build

To build the node extension and verify the tests run with the following command:
```
mvn verify
```

Make sure all code is commited as the snapshot version is determined by git commit timestamp.

An Eclipse update site will be made in `p2/target/repository` directory.
The update site can be used to perform a local installation.

${symbol_pound}${symbol_pound} Continuous Integration

Configuration files to run Continuous Integration builds on Linux (Travis-CI), OS X (Travis-CI) and Windows (AppVeyor) are present.

See `./.travis.yml` file how to trigger a Travis-CI build for every push or pull request.
Also see `./.travis.yml` file how to perform a [SonarCloud](https://sonarcloud.io/) analysis and code coverage.

See `./appveyor.yml` file how to run on https://www.appveyor.com .

To cite the KNIME node, a [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier) can be generated when a GitHub release is made. To enable, the GitHub repository must be connected on https://zenodo.org/account/settings/github/ . The connection must be made before creating a GitHub release.
To [cite the software](https://research-software.org/citation/developers/) a human and computer readable file called `CITATION.cff` is included.

${symbol_pound} Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (${knime_version}) - ${artifactId}.targetplatform/KNIME-AP-${knime_version}.target` target definition.

During import the Tycho Eclipse providers must be installed.

${symbol_pound}${symbol_pound} Tests

Tests for the node are in `tests/src` directory.
Tests can be executed with `mvn verify`, they will be run in a separate KNIME environment.
Test results will be written to `test/target/surefire-reports` directory.
Code coverage reports (html+xml) can be found in the `tests/target/jacoco/report/` directory.

The tests can be run against a different KNIME version using `mvn verify -Dtarget.file=KNIME-AP-4.1` where `4.1` is the major.minor version of KNIME and `KNIME-AP-4.1` is a target platform definition file called `targetplatform/KNIME-AP-4.1.target`.

${symbol_pound}${symbol_pound}${symbol_pound} Unit tests

Unit tests written in Junit4 format can be put in `tests/src/java`.

${symbol_pound}${symbol_pound}${symbol_pound} Workflow tests

See https://github.com/3D-e-Chem/knime-testflow${symbol_pound}3-add-test-workflow

The tests expect a `python3` executable in PATH with numpy, pandas and protobuf packages installed.

${symbol_pound}${symbol_pound} Speed up builds

Running mvn commands can take a long time as Tycho fetches indices of all p2 update sites.
This can be skipped by running maven offline using `mvn -o`.

${symbol_pound} New release

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
  2. Correct license
9. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
10. Update `CITATION.cff` file with new DOI, version, release date
