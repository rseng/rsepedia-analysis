# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## Unreleased

## [3.0.0] - 2021-02-12

## Removed

- Python based nodes in favour of web service nodes ([#25](https://github.com/3D-e-Chem/knime-kripodb/issues/25))

## [2.5.0] - 2019-09-25

### Deprecated

- Python based nodes in favour of web service nodes ([#25](https://github.com/3D-e-Chem/knime-kripodb/issues/25))

## [2.4.2] - 2019-07-02

### Changed

- Requires KNIME 4.0 ([#29](https://github.com/3D-e-Chem/knime-kripodb/issues/29))

## [2.4.1] - 2018-04-04

### Fixed

- Last update of web service in node description [#28]

## [2.4.0] - 2018-02-07

### Changed

- Allow KripoDB local nodes to use Python 2 or 3
- Require KNIME 3.5

## [2.3.0] - 2017-07-21

### Added

- Node to fetch pharmacophore from local db file (#12)
- Node to fetch pharmacophore from webservice (#12)

### Changed

- Require KNIME 3.3

## [2.2.1] - 2017-03-07

### Fixes

- Handle nulls in server responses (#23)

## [2.2.0] - 2017-03-05

### Added

- Pure Java nodes to fetch fragment information and similarites from web service (#17)

## Changed

- Python nodes only work for local files, no longer for web service urls

## [2.1.5] - 2017-03-01

Requires KripoDB v2.2.0 or higher.

### Added

- Fragment information node fails completely when just one of the fragment ids is wrong (#11)
- Fragment similarity node gives warnings for query fragment ids that can not be found

## [2.1.4]

Retracted due to unresolved compilation problems.

## [2.1.3] - 2017-01-20

### Added

- kripodb Python installation instruction in node description (#16)
- code coverage (#19)

### Changed

- Allow fragments db to be a webservice url (#15)
- Tests run against mocked web service (#18)

## [2.1.2] - 2016-11-22

### Fixed

- Correct fragments.sqlite url (#14)

## [2.1.1] - 2016-07-18

### Changed

- Nest Kripo nodes under /community/3D-e-Chem (#7)

## [2.1.0] - 2016-07-14

### Added

- Explained fragments db file choices (#6)
- Explained similar fragments matrix file options (#5)

### Changed

- Renamed distance to similarity (#9)

### Removed

- Python templates

## [2.0.1] - 2016-07-11

### Changed

- Moved PythonWrapper classes to own repo (https://github.com/3D-e-Chem/knime-python-wrapper)

## [2.0.0] - 2016-07-06

Version <2.0.0 used Python templates which could be selected as source code and adjusted in a text area in the `Python script` node.
Version >= 2.0.0 uses KripoDB Knime node which have there own dialog with combo boxes and file pickers.

### Added

- Workflow tests
- Run workflow tests on Travis-CI
- Codacy badge
- Check that Python packages are available before executing

### Changed

- Use KripoDB Knime node instead of Python node with a KripoDB template (#3)

## [1.0.3] - 2016-06-21

### Added

- bioisosteric replacement workflow

### Fixed

- Example 'Add fragment for hits' node returns too many rows (#2)

## [1.0.2] - 2016-05-25

### Added

- Use webservice in template and example (#1)

## [1.0.1] - 2016-04-18

### Added

- Distance matrix can be local file or kripodb webservice base url.

### Removed

- DOI for this repo, see https://github.com/3D-e-Chem/kripodb for DOI.

## [1.0.0] - 2016-02-12

### Added

- Python templates to use KripoDB package
- Example workflow on Github repo.

[unreleased]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.5.0...HEAD
[2.5.0]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.4.2...v2.5.0
[2.4.2]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.3.1...v2.4.2
[2.4.1]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.3.0...v2.3.1
[2.3.0]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.2.1...v2.3.0
[2.2.1]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.2.0...v2.2.1
[2.2.0]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.5...v2.2.0
[2.1.5]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.4...v2.1.5
[2.1.4]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.3...v2.1.4
[2.1.3]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.2...v2.1.3
[2.1.2]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.1...v2.1.2
[2.1.1]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.1.0...v2.1.1
[2.1.0]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.0.1...v2.1.0
[2.0.1]: https://github.com/3D-e-Chem/knime-kripodb/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/3D-e-Chem/knime-kripodb/compare/v1.0.3...v2.0.0
[1.0.3]: https://github.com/3D-e-Chem/knime-kripodb/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/3D-e-Chem/knime-kripodb/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/3D-e-Chem/knime-kripodb/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/3D-e-Chem/knime-kripodb/releases/tag/v1.0.0
# KripoDB KNIME nodes

KRIPO stands for [Key Representation of Interaction in POckets](http://dx.doi.org/10.1186/1758-2946-6-S1-O26).

[KNIME](http://www.knime.org) nodes for KripoDB (https://github.com/3D-e-Chem/kripodb).

[![Java CI with Maven](https://github.com/3D-e-Chem/knime-kripodb/workflows/Java%20CI%20with%20Maven/badge.svg)](https://github.com/3D-e-Chem/knime-kripodb/actions?query=workflow%3A%22Java+CI+with+Maven%22)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.kripodb%3Anl.esciencecenter.e3dchem.kripodb&metric=alert_status)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.kripodb%3Anl.esciencecenter.e3dchem.kripodb)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=nl.esciencecenter.e3dchem.kripodb%3Anl.esciencecenter.e3dchem.kripodb&metric=coverage)](https://sonarcloud.io/dashboard?id=nl.esciencecenter.e3dchem.kripodb%3Anl.esciencecenter.e3dchem.kripodb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597262.svg)](https://doi.org/10.5281/zenodo.597262)

# Installation

Requirements:

- KNIME, https://www.knime.org, version 4.0 or higher

Optionally:

- KripoDB Python package & data files, https://github.com/3D-e-Chem/kripodb,
  required to start local running kripo web service instead of the public [one](http://3d-e-chem.vu-compmedchem.nl/kripodb/ui).

Steps to get KripoDB nodes inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/updates`
4. Select --all sites-- in work with pulldown
5. Open KNIME 3D-e-Chem Contributions folder
6. Select KripoDB
7. Install software & restart

# Usage

See a minimal example workflow at [examples/Knime-KripoDB-example.zip](examples/Knime-KripoDB-example.zip).
The workflow can be run by importing it into KNIME as an archive.

Other workflows using the KripoDB nodes can be found at https://github.com/3D-e-Chem/workflows

# Development

Development requirements:

- Maven, https://maven.apache.org

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.kripodb.targetplatform/KNIME-AP-4.0.target` target definition.
6. A KNIME Analytics Platform instance can be started by right clicking on the `targetplatform/KNIME\ Analytics\ Platform.launch` file and selecting `Run As → KNIME Analytics Platform`. The KNIME instance will contain the target platform together with all extensions defined in the workspace.

During import the Tycho Eclipse providers must be installed.

## Build

```shell
mvn package
```

Jar has been made in `plugin/target` directory.
An Eclipse update site will be made in `p2/target/repository` repository.

## Tests

```shell
mvn verify
```

Tests in `tests` module will have been run with results in `test/target/surefire-reports` directory.
There are unit tests and workflow tests both are executed in the KNIME eclipse application.
See https://github.com/3D-e-Chem/knime-testflow for more information about workflow tests.

# Create web service client

The web service client is generated using [Swagger Code Generator](https://github.com/swagger-api/swagger-codegen) and stored inside `plugin/src/java/nl/esciencecenter/e3dchem/kripodb/ws/client/` directory.

1. Start [KripoDB webservice](https://github.com/3D-e-Chem/kripodb#web-service)

```shell
kripodb serve data/similarities.frozen.h5 data/fragments.sqlite data/pharmacophores.h5
```

2. Download swagger code generator

```shell
wget http://repo1.maven.org/maven2/io/swagger/swagger-codegen-cli/2.2.3/swagger-codegen-cli-2.2.3.jar
```

3. Generate a client for KripoDB web service

```shell
java -jar swagger-codegen-cli-2.2.3.jar generate \
--input-spec http://localhost:8084/kripo/swagger.json \
--output client \
--lang java \
--config swagger-codegen.config.json
```

4. Compile client

```shell
cd client
mvn package
```

5. Populate plugin with client source code and dependencies

```shell
mkdir ../plugin/lib
cp target/lib/gson-* target/lib/logging-interceptor-* target/lib/ok* target/lib/swagger-annotations-* ../plugin/lib/
rm -r ../plugin/src/java/nl/esciencecenter/e3dchem/kripodb/ws/client
cp -r src/main/java/nl/esciencecenter/e3dchem/kripodb/ws/client ../plugin/src/java/nl/esciencecenter/e3dchem/kripodb/ws/
```

6. Update plugin/META-INF/MANIFEST.MF, plugin/build.properties files to reflect contents of lib/

## New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Test node by installing it from local update site
5. Append new release to 3D-e-Chem update site
6. Make clone of https://github.com/3D-e-Chem/3D-e-Chem.github.io repo
7. Append release to 3D-e-Chem update site with `mvn install -Dtarget.update.site=<3D-e-Chem repo/updates>`
8. Commit and push changes in this repo and 3D-e-Chem.github.io repo
9. Create a Github release
10. Update Zenodo entry
11. Correct authors
12. To Related/alternate identifiers section add http://dx.doi.org/10.1186/1758-2946-6-S1-O26 as `is cited by this upload` entry.
13. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
14. Update CITIATION.cff with new DOI

# Create stub recordings for integration tests

The test workflow are tested against a mocked web server and not the actual http://3d-e-chem.vu-compmedchem.nl/kripodb site.
The mock server is called [WireMock](http://WireMock.org/) and normally gives empty responses.
To have WireMock server return filled responses, stubs stored in `tests/src/test/resources/` directory must be provided.
The stubs can be recorded by starting a WireMock server in recording mode by:

```
java -jar tests/lib/wiremock-standalone-2.5.0.jar --proxy-all="http://3d-e-chem.vu-compmedchem.nl/" \
--port=8089 --record-mappings --verbose --root-dir=tests/src/test/resources/
java -jar tests/lib/wiremock-standalone-2.5.0.jar --proxy-all="http://localhost:8084/" \
--port=8089 --record-mappings --verbose --root-dir=tests/src/test/resources/
```

Then in a KNIME workflow in the KripoDB nodes set the base path to http://localhost:8089.
Executing the workflow will fetch data from http://3d-e-chem.vu-compmedchem.nl/kripodb via the WireMock server and cause new stubs to be recorded in the `tests/src/test/resources/` directory.

To run the test workflows from inside KNIME desktop enviroment start the WireMock server in mock mode by:

```shell
java -jar tests/lib/wiremock-standalone-2.5.0.jar --port=8089 --verbose --root-dir=tests/src/test/resources/
```

Then import the test workflows in `tests/src/knime/` directory, select the workflow in the KNIME explorer and in the context menu (right-click) select `Run as workflow test`.
