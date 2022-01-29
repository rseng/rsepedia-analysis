# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.1] - 2021-05-10

### Fixed

* Dependabot alerts of hosted-git-info, lodash and handlebars

## [0.5.0] - 2020-04-02

### Fixed

* Publishing date not updated ([#13](https://github.com/iomega/zenodo-upload/issues/13))
* Unable to use concept DOI ([#12](https://github.com/iomega/zenodo-upload/issues/12))

## [0.4.1] - 2020-04-01

### Added

* Code examples in doc strings

### Fixed

* DraftDiscardedError and others not exported in entrypoint
* Vague error messages

## [0.4.0] - 2020-03-31

### Added

* Command line tool ([#1](https://github.com/iomega/zenodo-upload/issues/1))

### Fixed

* Example does not work ([#10](https://github.com/iomega/zenodo-upload/issues/10))

## [0.3.1] - 2020-03-27

### Fixed

* Error: Digest method not supported ([#7](https://github.com/iomega/zenodo-upload/issues/7))

## [0.3.0] - 2020-03-27

### Added

* Perform checksum check and discard on match ([#5](https://github.com/iomega/zenodo-upload/issues/5))

## [0.2.3] - 2020-03-25

### Removed

* Publishing of package to [https://github.com/iomega/zenodo-upload/packages](https://github.com/iomega/zenodo-upload/packages)

## [0.2.2] - 2020-03-25

### Fixes

* html url is undefined ([#2](https://github.com/iomega/zenodo-upload/issues/2))

## [0.2.1] - 2020-03-25

Tying to get automated publishing after a GitHub release is created to work.

## [0.2.0] - 2020-03-25

Initial release.

[Unreleased]: https://github.com//iomega/zenodo-upload/compare/v0.5.1...HEAD
[0.5.1]: https://github.com//iomega/zenodo-upload/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com//iomega/zenodo-upload/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com//iomega/zenodo-upload/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com//iomega/zenodo-upload/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com//iomega/zenodo-upload/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com//iomega/zenodo-upload/compare/v0.2.3...v0.3.0
[0.2.3]: https://github.com//iomega/zenodo-upload/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com//iomega/zenodo-upload/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com//iomega/zenodo-upload/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com//iomega/zenodo-upload/releases/tag/v0.2.0# Zenodo upload

Uploads file to Zenodo.

[![npm version](https://badge.fury.io/js/%40iomeg%2Fzenodo-upload.svg)](https://badge.fury.io/js/%40iomeg%2Fzenodo-upload)
[![CI](https://github.com/iomega/zenodo-upload/workflows/CI/badge.svg)](https://github.com/iomega/zenodo-upload/actions?query=workflow%3ACI)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=iomega_zenodo-upload&metric=alert_status)](https://sonarcloud.io/dashboard?id=iomega_zenodo-upload)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=iomega_zenodo-upload&metric=coverage)](https://sonarcloud.io/dashboard?id=iomega_zenodo-upload)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3726851.svg)](https://doi.org/10.5281/zenodo.3726851)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3805/badge)](https://bestpractices.coreinfrastructure.org/projects/3805)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8F-yellow)](https://fair-software.eu)

A JavaScript library to create a new version of a [Zenodo](https://zenodo.org) upload with a file.

Makes a draft copy of an existing Zenodo upload.
After overwriting the file and version the upload is published.

Can be used to create a [DOI](https://doi.org) for a updated data file.
A Zenodo upload must already exist using this library.

[API documentation](https://iomega.github.io/zenodo-upload/modules/_index_.html#).

## Install

```shell
npm install @iomeg/zenodo-upload
```

## Usage

### As a command line tool

```shell
npx --package @iomeg/zenodo-upload zenodo_upload [--sandbox] [--no-checksum] <deposition_id> <file> <version> <access_token>
```

To create new version (`1.2.3`) of [https://zenodo.org/record/1234567](https://zenodo.org/record/1234567) by uploading a local file called `somefile`.

```shell
npx --package @iomeg/zenodo-upload zenodo_upload 1234567 somefile 1.2.3 sometoken
```

The `sometoken` string has to be replaced with a valid [Zenodo access token](https://sandbox.zenodo.org/account/settings/applications/tokens/new/).

### As a library

To create new version of [https://zenodo.org/record/1234567](https://zenodo.org/record/1234567).

Example usage using NodeJS:

```javascript
const fs = require('fs');
const zenodo_upload = require('@iomeg/zenodo-upload').default;

const deposition_id = 1234567;
const filename = 'somefile.txt';
fs.writeFileSync(filename, 'sometext', 'utf8');
const version = '1.2.3';
const access_token = 'sometoken';

zenodo_upload(deposition_id, filename, version, access_token)
    .then(r => console.log(`New zenodo entry ${r.doi} created`))
    .catch(e => console.error(e))
;
```

Or in modern javascript

```javascript
import fs from 'fs';
import zenodo_upload from '@iomeg/zenodo-upload';

const deposition_id = 1234567;
const filename = 'somefile.txt';
await fs.promises.writeFile(filename, 'sometext', 'utf8');
const version = '1.2.3';
const access_token = 'sometoken';

const result = await zenodo_upload(deposition_id, filename, version, access_token);
console.log(`New zenodo entry ${result.doi} created`);
```

To run the example code you will need a valid Zenodo access token and a deposition id that can be written to by that token.

## Development

To install dependencies:

```shell
yarn install
```

To run the project in development/watch mode. Your project will be rebuilt upon changes.

```shell
yarn start
```

To bundle the package to the dist folder.

```shell
yarn build
```

To runs the test watcher (Jest) in an interactive mode. By default, runs tests related to files changed since the last commit.

```shell
yarn test
```

To run linter and fix fixable errors.

```shell
yarn lint --fix
```

## Credits

This project was bootstrapped with [TSDX](https://github.com/jaredpalmer/tsdx).

This project follows the [fair-software-nl](https://fair-software.nl) recommendations.
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

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/iomega/zenodo-upload/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/iomega/zenodo-upload/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of concensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing unit tests still work by running ``yarn test``;
1. make sure that no lint errors are present work by running ``yarn lint``;
1. make sure the API documentation is in-sync by running ``yarn apidocs``;
1. add your own unit tests and integration tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the Xenon repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request (have a look at some of our old pull requests to see how this works, for example [#1](https://github.com/iomega/zenodo-upload/pull/1)).
