# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [2.0.4] - 2020-01-14

### Fixed

* Support KNIME 4.3 ([#9](https://github.com/3D-e-Chem/knime-python-wrapper/pull/9) + [#10](https://github.com/3D-e-Chem/knime-python-wrapper/pull/10))

### Removed

* Support for KNIME versions older than 4.3

## [2.0.3] - 2019-07-02

### Changed

- In tests call `PythonWrapperTestUtils.init()` instead of `PythonWrapperTestUtils.materializeKNIMEPythonUtils()`.

### Fixed

- Force usage of python3 in PATH during tests ((#8)[https://github.com/3D-e-Chem/knime-python-wrapper/issues/8])

## [2.0.2] - 2019-06-27

### Changed

- Compatible with KNIME 4 (#6)[https://github.com/3D-e-Chem/knime-python-wrapper/issues/6]

## [2.0.1] - 2018-02-07

### Changed

* Replaced PythonKernel from org.knime.python to org.knime.python2 (#4)

### Fixed

* Give proper exception when Python file is missing (#1)

## [1.1.0] - 2017-02-23

### Added

* Set node warning message from Python
* Allow changing of number and names of input and output tables,
  default is 1 input table called `input_table` in Python
  and 1 output table called `output_table` in Python.
* Test coverage
* More tests
* Codacy integration

### Fixed

* Give proper exception when Python file is missing (#1)

## [1.0.0] - 2016-07-11

### Added

* Abstract PythonWrapperNode classes
* Test utility to run tests which call PythonKernel execute

[Unreleased]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v2.0.4...HEAD
[2.0.4]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v2.0.3...v2.0.4
[2.0.3]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v2.0.2...v2.0.3
[2.0.2]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v1.1.0...v2.0.1
[1.1.0]: https://github.com/3D-e-Chem/knime-python-wrapper/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/3D-e-Chem/knime-python-wrapper/releases/tag/v1.0.0
Abstract Python wrapper KNIME node and helpers.
Used for development of KNIME nodes calling Python scripts.

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-python-wrapper.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-python-wrapper)
[![Build status](https://ci.appveyor.com/api/projects/status/y7u4n23sjo25pyg8/branch/master?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/knime-python-wrapper/branch/master)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.python%3Anl.esciencecenter.e3dchem.python&metric=alert_status)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.python%3Anl.esciencecenter.e3dchem.python)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.python%3Anl.esciencecenter.e3dchem.python&metric=coverage)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.python%3Anl.esciencecenter.e3dchem.python)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4537256.svg)](https://doi.org/10.5281/zenodo.4537256)

The nodes in Scripting>Python folder of the node repository (nodes part of the `KNIME Python integration` plugin) the end-user needs to paste Python code in a text area in the node dialog.
Nodes derived from this repo will have a Python script included in their jar file and the dialog of the node will contain no source code text area.
The included Python script is not editable by the end-user, but can read options from dialog like the input column name.

# Usage

Instructions for KNIME node developers that want to call a Python script.
Several steps must be performed:

1. [Add update site](#1-add-update-site)
2. [Add dependency](#2-add-dependency)
3. [Implement](#3-implement)
4. [Write tests](#4-write-tests)

## Add update site

The releases of this repository are available in the `https://3d-e-chem.github.io/updates` update site.

Configure target platform by adding the `https://3d-e-chem.github.io/updates` update site with `Abstract Python wrapper KNIME node and helpers` software.

To make use of it in a [Tycho based project](https://github.com/3D-e-Chem/tycho-knime-node-archetype/), add to `targetplatform/KNIME-AP-4.0.target` file the following:

```
<location includeAllPlatforms="false" includeConfigurePhase="false" includeMode="planner" includeSource="true" type="InstallableUnit">
		<unit id="nl.esciencecenter.e3dchem.python.plugin" version="0.0.0"/>
		<repository location="https://3d-e-chem.github.io/updates"/>
</location>
```

## Add dependency

To implement the node a dependency is needed for the plugin add tests.
To do this add `nl.esciencecenter.e3dchem.knime.python` as a required plugin to the `plugin/META-INF/MANIFEST.MF` and `tests/META-INF/MANIFEST.MF` file.

## Implement

### Config

Create your node config class extended from the `nl.esciencecenter.e3dchem.python.PythonWrapperNodeConfig` class.

Overwrite the constructor to add required Python modules with the `addRequiredModule("<module name>")` method.

PythonWrapperNodeConfig class are configured for a single input table called `input_table` and a single output table called `output_table`.
To change the the number and names of input and/or output tables, override the constructor.

### Dialog

In your nodes dialog the Python options panel should be added by adding the following to the dialog constructor

```java
pythonOptions = new PythonOptionsPanel<PredictMetabolitesConfig>();
addTab("Python options", pythonOptions);
```

To save the Python options to disk you must call the `pythonOptions.saveSettingsTo(config)` followed by `config.saveToInDialog(settings)` in the `save*To()` method of the dialog.
To load the Python options from disk you must call the `config.loadFromInDialog(settings)` followed by `pythonOptions.loadSettingsFrom(config)` in the `load*From()` methods of the dialog.

### Python script

Inside Python script the following variables are special:

- `options`, dictionary filled from Java with PythonWrapperNodeConfig.getOptionsValues() method, to read from
- `input_table`, Pandas Dataframe with input data table, to read from
- `output_table`, Pandas Dataframe with output data table, to assign with value
- `flow_variables`, dictionary of flow variables, to get or put key/value pairs
  - `flow_variables['warning_message']`, if key exists then value will be set as warning message of node

The Python script should be located in same directory as the model.

### Model

Create your node model class extended from the `nl.esciencecenter.e3dchem.python.PythonWrapperNodeModel` class.
Overwrite the `python_code_filename` fields in the constructor, this is the name of the Python script.

## Write tests

To run tests which execute the node a call to `PythonWrapperTestUtils.init()` is required
this will add the KNIME-Python utility scripts to the Python path and configure it to use the `python3` executable in current PATH.

# Build

To build the plugin and run the tests run the following command:

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
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.python.targetplatform/KNIME-AP-4.0.target` target definition.
6. A KNIME Analytics Platform instance can be started by right clicking on the `targetplatform/KNIME\ Analytics\ Platform.launch` file and selecting `Run As → KNIME Analytics Platform`. The KNIME instance will contain the target platform together with all extensions defined in the workspace.

During import the Tycho Eclipse providers must be installed.

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Append new release to an update site
5. Make clone of an update site repo
6. Append release to the update site with `mvn install -Dtarget.update.site=<path to update site>`
7. Commit and push changes in this repo and update site repo.
8. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
9. Create a GitHub release
10. Update CITATION.cff file to reflect DOI which was generated by GitHub release on https://zenodo.org
