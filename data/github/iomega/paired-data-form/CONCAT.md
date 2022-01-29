# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.2] - 2021-03-03

## Added

- Conform to Bioschema profile ([#177](https://github.com/iomega/paired-data-form/issues/177))
- Allow notifications to be disabled with `SLACK_TOKEN=false` environment variable
- Mention CC BY 4.0 license on add and download page
- Health check API endpoint ([#171](https://github.com/iomega/paired-data-form/issues/171))

## Fixed

- URL in notification message is internal ([#168](https://github.com/iomega/paired-data-form/issues/168))
- Return 404 http error when path is not /api, a file or part of router paths ([#174](https://github.com/iomega/paired-data-form/issues/174))

## [0.9.1] - 2021-01-21

### Security

- Axios vulnerability ([CVE-2020-28168](https://github.com/advisories/GHSA-4w2v-q235-vp99))
- node-notifier vulnerability ([CVE-2020-7789](https://github.com/advisories/GHSA-5fw9-fq32-wv5p))

## [0.9.0] - 2020-12-07

### Added

- Notify admins on Slack when new project is ready for review ([#155](https://github.com/iomega/paired-data-form/issues/155))
- Notify admins on Slack when data archive on Zenodo is updated ([#163](https://github.com/iomega/paired-data-form/issues/163))

## [0.8.2] - 2020-11-05

### Fixed

- Use BioSample DB link to NCBI ([#162](https://github.com/iomega/paired-data-form/issues/162))

## [0.8.1] - 2020-09-28

### Added

- Links to methods page ([#159](https://github.com/iomega/paired-data-form/issues/159))

## [0.8.0] - 2020-09-22

### Added

- Column check for uploading of genome links ([#154](https://github.com/iomega/paired-data-form/pull/154))
- Methods page ([#158](https://github.com/iomega/paired-data-form/pull/158))

### Fixed

- Render publications delimited by spaces ([#152](https://github.com/iomega/paired-data-form/pull/152))

## [0.7.0] - 2020-07-10

### Added

- SEO optimizations like static sitemap, dynamic sitemap for projects and each project has structured data ([#148](https://github.com/iomega/paired-data-form/pull/148))

### Changed

- Titles synced with draft paper ([#146](https://github.com/iomega/paired-data-form/pull/146))

### Fixed

- Correct title for acetonitrile search example ([#141](https://github.com/iomega/paired-data-form/issues/141))
- Consistent casing of platform name ([#143](https://github.com/iomega/paired-data-form/issues/143))
- Smiles render with own optional scrollbar ([#147](https://github.com/iomega/paired-data-form/issues/147))

## [0.6.2] - 2020-04-23

### Added

- Paging projects ([#137](https://github.com/iomega/paired-data-form/issues/137))

### Changed

- Sort projects moved from web application to elastic search ([#138](https://github.com/iomega/paired-data-form/issues/138))
- Allow search and filter to be combined

## [0.6.1] - 2020-04-16

### Added

- Added ionization modes to stats page ([#132](https://github.com/iomega/paired-data-form/issues/132))
- Search query examples ([#132](https://github.com/iomega/paired-data-form/issues/132))

### Fixed

- Enrichments cause Limit of total fields exceeded error in elastic search ([#131](https://github.com/iomega/paired-data-form/issues/131))

## [0.6.0] - 2020-04-16

Search functionality using elastic search has been added.

### Added

- Full text search functionality ([#123](https://github.com/iomega/paired-data-form/issues/123))
- Filter projects on statistic functionality ([#124](https://github.com/iomega/paired-data-form/issues/124))

### Changed

- Replaced http with https ([#126](https://github.com/iomega/paired-data-form/issues/126))

### Fixed

- Schema version not constant ([#127](https://github.com/iomega/paired-data-form/issues/127))

## [0.5.0] - 2020-04-02

### Added

- Mention JSON schema on add form ([#115](https://github.com/iomega/paired-data-form/issues/115))
- Project list can be sorted by clicking on column header ([#117](https://github.com/iomega/paired-data-form/issues/117))
- Show software version on about page ([#109](https://github.com/iomega/paired-data-form/issues/109))
- Download page, shows DOI of dataset archive ([#109](https://github.com/iomega/paired-data-form/issues/109))

### Changed

- Project list sorted on metabolite id ([#117](https://github.com/iomega/paired-data-form/issues/117))

### Fixed

- Submitter email is used where PI email is expected ([#118](https://github.com/iomega/paired-data-form/issues/118))
- Download button on pending page yields incorrectly formatted json files ([#120](https://github.com/iomega/paired-data-form/issues/120))

## [0.4.0] - 2020-03-20

This version requires following migration steps.

- JSON schema changed to version 2. To migrate all projects in data/ dir from 1 to 2 run

  ```shell
  # Backup
  cp -a data/ backup-$(date -I)/
  # Perform migration
  docker-compose exec api npm run migrate
  # Validate projects
  docker-compose exec api npm run validateall
  # Fix any validation errors and rerun validation until all projects are valid
  # For example use VS Code extension https://marketplace.visualstudio.com/items?itemName=tiibun.vscode-docker-ws to edit file in Docker container
  # Edit CTRL-SHIFT-p, select dockerws command, select`paired-data-form_api_1` as Docker container and `/data` as path to open.
  code .
  # Restart api so new updated files are reindexed
  docker-compose restart api
  ```

- The enrichment of projects has been improved. To recreate enrichments of all projects run

  ```shell
  # Drop existing enrichment with
  docker-compose exec redis sh
  redis-cli --scan --pattern keyv:enrichment:* | xargs redis-cli del
  exit
  # Recreate all enrichments
  docker-compose exec api npm run enrich
  ```

### Added

- About page ([#86](https://github.com/iomega/paired-data-form/issues/86))
- Limit log size of Docker containers ([#89](https://github.com/iomega/paired-data-form/issues/89))
- Denied projects moved to thrash dir ([#95](https://github.com/iomega/paired-data-form/issues/95))
- Stats page ([#64](https://github.com/iomega/paired-data-form/issues/64))
- Fetch species from BioSample, ENA and JGI for each genome
- Download button for pending project in review section ([#98](https://github.com/iomega/paired-data-form/issues/98))
- Submitter name column to project lists ([#101](https://github.com/iomega/paired-data-form/issues/101))
- Commands to validate one or all projects ([#100](https://github.com/iomega/paired-data-form/issues/100))
- Command to perform data migrations ([#110](https://github.com/iomega/paired-data-form/pull/110))
- Second submitter ([#97](https://github.com/iomega/paired-data-form/issues/97))
- OpenAPI specification ([#112](https://github.com/iomega/paired-data-form/issues/112))
- Manual for developers who wish to consume the api

### Fixed

- Unable to submit large project ([#88](https://github.com/iomega/paired-data-form/issues/88))
- Spelling errors ([#87](https://github.com/iomega/paired-data-form/issues/87))
- Render error when growth medium is not set ([#92](https://github.com/iomega/paired-data-form/issues/92))
- Unable to download project from add form([#111](https://github.com/iomega/paired-data-form/issues/111))

### Changed

- Dropped Caddy web server from docker-compose, use nginx from app and external reverse proxy for https
- Download project directly using web service instead of data-url
- BGC number to BGC accession aka 1234 to BGC0001234 ([#94](https://github.com/iomega/paired-data-form/issues/94))
- Require more fields in Gene cluster - Mass spectra links ([#94](https://github.com/iomega/paired-data-form/issues/94))
- Increased JSON schema version to 2 due to issue #94

## [0.3.0] - 2019-12-11

### Added

- Confirmation of submission to review ([#74](https://github.com/iomega/paired-data-form/issues/74))
- Resins field to extraction method ([#76](https://github.com/iomega/paired-data-form/issues/76))
- erlemeyer flask option to aeration vessel ([#76](https://github.com/iomega/paired-data-form/issues/76))
- Links to description of fields and check list ([#76](https://github.com/iomega/paired-data-form/issues/76))
- Warning to not include spaces in urls ([#75](https://github.com/iomega/paired-data-form/issues/75))

### Fixed

- Do validation on labels when selected in links sections ([#73](https://github.com/iomega/paired-data-form/issues/73))
- Validate uploaded JSON documents ([#78](https://github.com/iomega/paired-data-form/issues/78))
- GNPS task id link broken ([#81](https://github.com/iomega/paired-data-form/issues/81))

## [0.2.0] - 2019-07-05

### Added

- Intro page ([#45](https://github.com/iomega/paired-data-form/issues/45))
- Password protected review section to review pending projects ([#45](https://github.com/iomega/paired-data-form/issues/45))
- Page with list of projects ([#45](https://github.com/iomega/paired-data-form/issues/45))
- Page to show single project ([#45](https://github.com/iomega/paired-data-form/issues/45))
- Web service to store project on disk as JSON documents ([#46](https://github.com/iomega/paired-data-form/issues/46))
- Task queue to enrich projects ([#46](https://github.com/iomega/paired-data-form/issues/46))
- Enrich project by fetching organism name based on genome identifier ([#46](https://github.com/iomega/paired-data-form/issues/46))

### Changed

- original form is now for adding a project for review ([#45](https://github.com/iomega/paired-data-form/issues/45))
- Added metabolights study id to genome ([#54](https://github.com/iomega/paired-data-form/issues/54))
- Made which fields are required more clear ([#42](https://github.com/iomega/paired-data-form/issues/42))
- Replaced run command from `yarn` to `docker-compose`

## [0.1.0] - 2019-05-01

Initial release.

[unreleased]: https://github.com/iomega/paired-data-form/compare/v0.9.2...HEAD
[0.9.2]: https://github.com/iomega/paired-data-form/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/iomega/paired-data-form/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/iomega/paired-data-form/compare/v0.8.2...v0.9.0
[0.8.2]: https://github.com/iomega/paired-data-form/compare/v0.8.1...v0.8.2
[0.8.1]: https://github.com/iomega/paired-data-form/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/iomega/paired-data-form/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/iomega/paired-data-form/compare/v0.6.2...v0.7.0
[0.6.2]: https://github.com/iomega/paired-data-form/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/iomega/paired-data-form/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/iomega/paired-data-form/compare/v0.5.9...v0.6.1
[0.5.0]: https://github.com/iomega/paired-data-form/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/iomega/paired-data-form/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/iomega/paired-data-form/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/iomega/paired-data-form/compare/v0.0.1...v0.2.0
[0.1.0]: https://github.com/iomega/paired-data-form/releases/tag/v0.1.0
# Pairing Omics Data Platform

Linking mas spectra and genomic information to discover new chemistry.

* Links MS/MS mass spectra with genome, sample preparation, extraction method and instrumentation method
* Links biosynthetic gene cluster with MS^2 mass spectra

A web application for storing paired omics data projects.

The [JSON schema (app/public/schema.json)](app/public/schema.json) describes the format of an project.

[![Node.js CI](https://github.com/iomega/paired-data-form/workflows/CI/badge.svg)](https://github.com/iomega/paired-data-form/actions?query=workflow%3A%22CI%22CI)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=iomega_paired-data-form&metric=alert_status)](https://sonarcloud.io/dashboard?id=iomega_paired-data-form)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=iomega_paired-data-form&metric=coverage)](https://sonarcloud.io/dashboard?id=iomega_paired-data-form)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3757/badge)](https://bestpractices.coreinfrastructure.org/projects/3757)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2656630.svg)](https://doi.org/10.5281/zenodo.2656630)

## Documentation

Developer and admin manuals can be found in [manuals/](manuals/) directory.

## Contributing

If you want to contribute to the platform, have a look at the [contribution guidelines](CONTRIBUTING.md).

## Architecture

The Pairing Omics Data Platform consists of:

1. A Web application, user interface, see [app/](app/) directory
2. An API web service, service responsible for storing projects, see [api.](api/) directory

[![Architecture diagram](app/src/pages/methods/architecture.svg)](app/src/pages/methods/architecture.svg)

The platform is implemented using Javascript based web service and a React (v16.13.1) based web application. The web application renders the submission form from the JSON schema. The web service stores each project as a file on disk. The application offers full text search functionality via web services using an elastic search (v7.6.2) index. The web service uses a redis queue (v5.0.5) to schedule jobs to fetch more information about the public identifiers and to upload the projects to Zenodo each month. For example, the scientific species name is fetched from GenBank using the public genome identifiers in the project. The web service has an OpenAPI (v3.0.3) specification ([https://www.openapis.org/](https://www.openapis.org/)) which can be used to submit and retrieve projects in a programmatic manner. The platform runs using Docker Compose (v1.25.4) with containers for the web application, web service and redis queue.

## Run using Docker compose

The application can be configured using environment variables:

* PORT, http port application is running on. Default is 8443.
* SHARED_TOKEN, token required to login to review area.
* ZENODO_ACCESS_TOKEN, Zenodo access token used for uploading database to Zenodo.
* ZENODO_DEPOSITION_ID, Zenodo deposition identifier used for uploading database to Zenodo. Set to -1 to disable scheduled uploading.
* SLACK_TOKEN, Token of Slack app with chat:write permission in workspace of channel
* SLACK_CHANNEL, Slack channel in which service should post messages

The environment variables can be set in the terminal or be put in a `.env` file.

```shell
docker-compose up -d --build
```

Starts application, api webservice and reverse proxy on [http://localhost:8443](http://localhost:8443).
Project JSON files are stored in a `./data/` directory.

To run on production put application behind a reverse proxy web server with a proper domain and secure transfer with https.

## New release

This chapter is for developers of the platform.

To make a new release of the platform do:

1. Determine new version of release, using semantic versioning (x.y.z)
2. Add version to [CHANGELOG.md](CHANGELOG.md)
    * Create a new `##` chapter for the new version
    * Update version links at bottom of CHANGELOG
3. Set new version of api web service by

    ```shell
    cd api
    npm version x.y.z
    ```

4. Set new version of web application by

    ```shell
    cd app
    npm version x.y.z
    ```

5. Commit & push changes
6. Create a GitHub release
7. On [https://doi.org/10.5281/zenodo.2656630](https://doi.org/10.5281/zenodo.2656630)
    * Update author list
    * Add `https://doi.org/10.5281/zenodo.3736430`, `is compiled/create by this upload` as `Dataset` in related identifiers section.
# Security Policy

## Supported Versions

Only the latest minor version will be supported with security updates.

## Reporting a Vulnerability

To report a vulnerability please create an [issue](https://github.com/iomega/paired-data-form/issues).
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at s.verhoeven@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq# Contributing guidelines

We welcome any kind of contribution to our platform, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you want to contribute a project;
1. you want to use the projects in the platform;
1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You want to contribute a project

You are free to add projects on [https://pairedomicsdata.bioinformatics.nl/add](https://pairedomicsdata.bioinformatics.nl/add).
After submission, the project will be reviewed and if approved will appear in the [list of projects](https://pairedomicsdata.bioinformatics.nl/projects).

## You want to use the data in the platform

Each project page on [https://pairedomicsdata.bioinformatics.nl](https://pairedomicsdata.bioinformatics.nl) has a download button, which gives you the JSON document of the project.

To access the platform in a programmatic manner see the [developer manual](manuals/developers.md).

Please inform us if you found a nice way to use the data.

## You have a question

1. use the search functionality [here](https://github.com/iomega/paired-data-form/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/iomega/paired-data-form/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. setup development environment by following instructions in [api/README.md](api/README.md) and [app/README.md](app/README.md);
1. make sure the existing tests still work by running ``npm run test`` in `api/` and/or `app/` directory;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the Python Template repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
# Developers

This manual is for developers who want to interact with the Paired Omics Data Platform in a programmatic manner.

How to contribute as a developer to the platform is described in [../CONTRIBUTING.md](../CONTRIBUTING.md).

## Web service

To consume the [web service](https://pairedomicsdata.bioinformatics.nl/api) some examples are provided and a specification.

### Examples

Using [httpie](https://httpie.org) as a http client.

```bash
# Install httpie
pip install httpie
# Submit a project
https -j pairedomicsdata.bioinformatics.nl/api/projects < ../app/public/examples/paired_datarecord_MSV000078839_example.json
# List pending projects
https pairedomicsdata.bioinformatics.nl/api/pending/projects 'Authorization: Bearer ashdfjhasdlkjfhalksdjhflak'
# Approve pending projects
https POST pairedomicsdata.bioinformatics.nl/api/pending/projects/a6f87bdf-6998-4336-8fb1-eca5b4fdb882.1 'Authorization: Bearer ashdfjhasdlkjfhalksdjhflak'
# List summary of projects
https pairedomicsdata.bioinformatics.nl/api/projects
# Get single project
https pairedomicsdata.bioinformatics.nl/api/projects/a6f87bdf-6998-4336-8fb1-eca5b4fdb882.1
# Get stats
https pairedomicsdata.bioinformatics.nl/api/stats
```

Change `a6f87bdf-6998-4336-8fb1-eca5b4fdb882.1` project id to what previous requests return.
Change `ashdfjhasdlkjfhalksdjhflak` token to what you got as an authentication token.

### OpenAPI specification

The web services ships with an [OpenAPI specification](https://www.openapis.org/).
The specification is available at [https://pairedomicsdata.bioinformatics.nl/openapi.yaml](https://pairedomicsdata.bioinformatics.nl/openapi.yaml).

To try out the specification you can use the [Swagger UI](https://swagger.io/tools/swagger-ui/) at [/api/ui/?url=/openapi.yaml](https://pairedomicsdata.bioinformatics.nl/api/ui/?url=/openapi.yaml). When trying out a local deployment don't forget to pick the right server in the Swagger UI.
# Manual for platform administrators

## Reviewing

When a project has been submitted it has to be reviewed and either be approved or denied.

To review a project

1. Goto [https://pairedomicsdata.bioinformatics.nl/pending](https://pairedomicsdata.bioinformatics.nl/pending) for the list of pending projects
1. Login with supplied password if not already logged in
1. Perform review by either
    1. Clicking on the `Metabolomics project identifier` to see the project rendered as if it was public.
    2. Download the JSON file and look at it offline
1. Click the Approve/Deny button
    1. Clicking `Approve` button will immediately make the project public for everyone
    2. Clicking `Deny` button will immediately move the project to the recycle-bin

The recycle-bin is the `data/thrash/` directory on the server.

## Deployment

Tasks for the person deploying platform.

### Schema change

When schema has changed the documents in `./data/` directory need to be migrated.

```shell
# Create a backup before starting migration
cp -a data/ backup-$(date -I)/
# Perform migration
docker-compose exec api npm run migrate
# Validate projects
docker-compose exec api npm run validateall
# Fix any validation errors by editing the project files and rerun validation until all projects are valid
```

### Rebuild

When the code has changed the Docker images have to be rebuild and restarted with

```shell
docker-compose stop
docker-compose up -d --build
```

### Enrich

Normally all submitted projects are enriched with for the species scientific name and other things.
Those enrichments are stored in the Redis database.
It can happen that, for example when enrichment code is improved, the Redis database does not have enrichments for a project.

To enrich existing projects run

```shell
docker-compose exec api npm run enrich
```
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Web application for paired data for mapping between genomic and metabolomic (mass spectra) projects

This project was bootstrapped with [Create React App](https://github.com/facebookincubator/create-react-app).

Some information on how to perform common tasks is available in the [guide](https://github.com/facebookincubator/create-react-app/blob/master/packages/react-scripts/template/README.md).

## Install

Atfer cloning the repo, it's dependencies must be installed.

Requires [nodejs](https://nodejs.org) v10 or higher.

```shell
npm install
```

## Development

```shell
npm start
```

Runs the app in development mode. Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

The page will automatically reload if you make changes to the code. You will see the build errors and lint warnings in the console.

## Docker

Build with

```shell
docker build -t iomega/paired-data-form app
```

Run with

```shell
docker run -d -p 8887:80 iomega/paired-data-form
```

Goto [http://localhost:8887](http://localhost:8887).

## Propagate schema changes

When `public/schema.json` changes the `src/schema.ts` must also be updated using

```shell
npm run schema2ts
```
# Api web service for Pairing Omics Data Platform

## Architecture

The Pairing Omics Data Platform api web service use a directories to store project json files.

* pending/, any projects added via POST to `/api/projects` or `/api/projects/:id` will end up here. Ready for reviewing
* approved/, all public visible projects. Once a project is approved it is moved from the `pending/` directory to the `approved/` directory.
* archive/, when an existing project has been edited and approve, the previously approved project file is moved here.
* thrash/, when an pending project has been denied it is moved from the `pending/` directory to the `thrash/` directory.

Extra information is gathered for some of the identifiers and urls in a project json files. This extra information is called an enrichment. [Redis](https://redis.io/) is used to store enrichments and as job queue.

[Elastic search](https://www.elastic.co/elasticsearch/) is used to perform full text search and filtering of projects.

Consumption of the web service is explained in the [developers manual](../manuals/developers.md).

## Install

```shell
npm install
```

## Configure

Create `./.env` file, use `./.env.example` as an example.

### Zenodo access token

To publish the data collection to Zenodo, a [Zenodo personal access token](https://zenodo.org/account/settings/applications/tokens/new/) is needed.
During token generation check the `deposit:actions` and `deposit:write` scopes.

### Slack integration

The service can post messages to a Slack channel when a project is submitted or pushed to Zenodo.

1. [Create an Slack app](https://api.slack.com/apps), pick any name / workspace you want
1. Optionally. Add collaborators so your are not the only one who can edit the Slack app
1. Add `chat:write` to bot token scopes
1. Install to workspace and allow bot to write chats
1. Copy token to `../.env:SLACK_TOKEN` value
1. Create channel to post messages to and set it's id (last bit of channel url) as `../.env:SLACK_CHANNEL` value
1. In channel select `Add an app` and add the created app

People that should be reviewing incoming projects should be invited to the Slack channel and be supplied with the review credentials.

## Build & Run

Redis can be started using docker:

```shell
docker run --name some-redis -d -p 6379:6379 redis
```

Or using [docker-compose](../README.md#run-using-docker-compose) to run the whole stack.

Elastic search can be started using docker:

```shell
docker run --name some-elasticsearch -d -p 9200:9200 -e "discovery.type=single-node" elasticsearch:7.6.2
```

Or using [docker-compose](../README.md#run-using-docker-compose) to run the whole stack.

Build with

```shell
npm run build
```

Run web service

```bash
npm run serve
```

## Docker

Build with

```shell
cd ..
docker build -t iomega/podp-api -f api/Dockerfile .
```

Run using `./data` dir as datadir with

```shell
docker run -d -p 8886:3001 --user $(id -u) -v $PWD/data:/data iomega/podp-api
```

Will start api web service on [http://localhost:8886](http://localhost:8886).

## Propagate schema changes

When `../app/public/schema.json` changes the `src/schema.ts` must also be updated using

```shell
npm run schema2ts
```

## Propagate api changes

When the api changes the OpenAPI specification at [../app/public/openapi.yaml](../app/public/openapi.yaml) should also be updated. Changes to schema (`../app/public/schema.json`) do not require updates to the OpenAPI spec as it is imported into the spec.

## Schema migration and validation

Whenever the schema is changed the projects in the data/ directory must be migrated and validated.

The projects can be validated with

```shell
npm run validateall
```

## Publish collection on Zenodo

The approved projects can be published to Zendo by running

```shell
npm run publish2zenodo -- --sandbox --deposition_id xxxx --access_token xxxxxxxxxxxxx
```

Or when running with docker-composose in production use

```shell
docker-compose exec api npm run publish2zenodo
```
