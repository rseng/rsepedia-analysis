# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.7.2] - 2021-02-01
### Changed
- Changed the url of Swagger UI from `/ui` to the base url `/`

## [0.7.1] - 2020-09-15
### Added
- PUT method used to update metadata

### Changed
- Replaced `flask-restplus` with `connexion` package for OpenAPI generation
- Refactored code based on OpenAPI specification
- Removed ending '/' in endpoints, e.g. '/catalog/' was changed to '/catalog'.
- Removed Bottle and Paste dependecies
- Removed docker compose file for production run


## [0.7.0] - 2020-06-12
### Added
- POST and DELETE methods
- SHACL validator to validate different layers of metadata
- Support for SPARQL database, such as Virtuoso RDF Triple Store
- Docker compose file
- Running in production mode

### Changed
- Improved GET and POST methods
- Updated support for common RDF serializations
- Removed loading metadata for the command tool `fdp-run`
- Removed support for INI-based configuration files
- Updated docker file
- Improved unit testing
- Updated README
- Use Bottle and Paste for production run

## [0.6.0] - 2020-02-10
### Added
- Loading meta-data from RDF/Turtle file (TTL).
- Automatic Swagger API (OpenAPI) generation (via `flask-restplus`).

### Changed
- Documentation in README
- Use Flask instead of Bottle

### Fixed
- Test syntax
# FAIR Data Point (FDP)

[![PyPI](https://img.shields.io/pypi/v/fairdatapoint)](https://pypi.org/project/fairdatapoint/)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/nlesc/fairdatapoint?label=Docker)](https://hub.docker.com/r/nlesc/fairdatapoint)
[![DOI](https://zenodo.org/badge/37470907.svg)](https://zenodo.org/badge/latestdoi/37470907)
[![Research Software Directory](https://img.shields.io/badge/RSD-FAIRDataPoint-red)](https://research-software.nl/software/fairdatapoint)
[![Build_Test](https://github.com/fair-data/fairdatapoint/actions/workflows/build_test.yml/badge.svg)](https://github.com/fair-data/fairdatapoint/actions/workflows/build_test.yml)
[![Coverage Status](https://coveralls.io/repos/github/fair-data/fairdatapoint/badge.svg?branch=master)](https://coveralls.io/github/fair-data/fairdatapoint?branch=master)

## Overview
Python implementation of FAIR Data Point.

FDP is a RESTful web service that enables data owners to describe and to expose their datasets (metadata) as well as data users to discover more information about available datasets according to the [FAIR Data Guiding Principles](http://www.force11.org/group/fairgroup/fairprinciples). In particular, FDP addresses the findability or discoverability of data by providing machine-readable descriptions (metadata) at four hierarchical levels:

*FDP -> catalogs -> datasets -> distributions*

FDP software specification can be found [here](https://github.com/FAIRDataTeam/FAIRDataPoint-Spec/blob/master/spec.md).
Other implementations are also available, e.g. [Java implementation](https://github.com/DTL-FAIRData/FAIRDataPoint)

### Demo server
A demo server of this Python implementation is http://fdp.fairdatapoint.nl/

## Installation

To install FDP, do

From pypi
```bash
pip install fairdatapoint
```

Or from this repo, but note that the in-development version might be unstable,
```bash
git clone https://github.com/fair-data/fairdatapoint.git
cd fairdatapoint
pip install .
```

## Running
```bash
fdp-run localhost 80
```

The [Swagger UI](https://swagger.io/tools/swagger-ui/) is enabled for FDP service, and you can have a try by visiting http://localhost.

## Unit testing
Run tests (including coverage) with:

```bash
pip install .[tests]
pytest
```

## Deploy with Docker

Check [fairdatapoint-service](https://github.com/CunliangGeng/fairdatapoint-service).

## Deploy without Docker

Before deploying FDP, it's necessary to first have a running SPARQL database which can be used to store metadata.

```
pip install fairdatapoint

# fdp-run <host> <port> --db=<sparql-endpoint>
# Let's assume your <host> is 'example.com' and <sparql-endpoint> is 'http://example.com/sparql', then
fdp-run example.com 80 --db='http://example.com/sparql'
```

## Web API documentation

FAIR Data Point (FDP) exposes the following endpoints (URL paths):

| Endpoint |  GET  | POST |  PUT | DELETE     |
|--------------|:--------------:|:-----------------:|:--------------:|:--------------:
| fdp | Output fdp metadata | Create new fdp metadata | Update fdp metadata | Not Allowed |
| catalog     | Output all catalog IDs   | Create new catalog metadata| Not Allowed | Not Allowed |
| dataset     | Output all dataset IDs   | Create new dataset metadata| Not Allowed | Not Allowed |
| distribution  | Output all distribution IDs  | Create new distribution metadata| Not Allowed | Not Allowed |
| catalog/\<catalogID\> | Output \<catalogID\> metadata | Not Allowed | Update \<catalogID\> metadata | Remove \<catalogID\> metadata |
| dataset/\<datasetID\> | Output \<datasetID\> metadata | Not Allowed | Update \<datasetID\> metadata | Remove \<datasetID\> metadata |
| distribution/\<distributionID\> | Output \<distributionID\> metadata | Not Allowed | Update \<distributionID\> metadata | Remove \<distributionID\> metadata |


### Access endpoints to request metadata programmatically

FDP: `curl -iH 'Accept: text/turtle' [BASE URL]/fdp`

Catalog: `curl -iH 'Accept: text/turtle' [BASE URL]/catalog/catalog01`

Dataset: `curl -iH 'Accept: text/turtle' [BASE URL]/dataset/dataset01`

Distribution: `curl -iH 'Accept: text/turtle' [BASE URL]/distribution/dist01`

### FDP supports the following RDF serializations (MIME-types):
* Turtle: `text/turtle`
* N-Triples: `application/n-triples`
* N3: `text/n3`
* RDF/XML: `application/rdf+xml`
* JSON-LD: `application/ld+json`


## Issues and Contributing
If you have questions or find a bug, please report the issue in the
[Github issue channel](https://github.com/fair-data/fairdatapoint/issues).

If you want to contribute to the development of FDP, have a look at the
[contribution guidelines](CONTRIBUTING.rst).