# Developer documentation

This [action](https://docs.github.com/en/actions/creating-actions/creating-a-javascript-action) is written in
[Typescript](https://www.typescriptlang.org)
and makes use of the
[Github @actions packages](https://github.com/actions/toolkit/blob/master/README.md#packages)

The underlining design can be found in [DESIGN.md](DESIGN.md).

## Requirements

This tool relies on the availability of [Node.js](https://nodejs.org/) and
[Docker](https://docs.docker.com/get-docker/).

Please verify that you have `Node.js` and the related package manager `npm`, and `docker` available on your
system. Make sure that the version of `Node.js` is at least `12`.

```bash
$ node --version
v14.17.0
$ npm --version
6.14.13
$ docker --version
Docker version 20.10.6, build 370c289
```

`Node.js` and `npm` can be downloaded in one package from [nodejs.org](https://nodejs.org/en/). And here are
[instructions for upgrading `Node.js`](https://phoenixnap.com/kb/update-node-js-version#ftoc-heading-3).

Install the dependencies

```bash
$ npm install
```

## Build

Build the typescript and package it for distribution

```bash
$ npm run build && npm run package
```

## Run unit test

The tests are stored in the directory `__tests__` and are written using
[jestjs](https://jestjs.io/).

Run the tests :heavy_check_mark:

```bash
$ npm test

 PASS  ./index.test.js
  ✓ throws invalid number (3ms)
  ✓ wait 500 ms (504ms)
  ✓ test runs (95ms)

...
```

To get information about the test coverage, run the tests with
`coverage npm test -- --coverage` and
examine the file `coverage/lcov-report/index.html`

## Linting

The code in the `src` directory can be linted with:

```bash
npm run lint
```

## Formatting

Some of the linting error can be fixed with formatting:

```bash
npm run format
```

## Run the analysis

### On the current repository

The tool will analyze the license dependencies in current Github
repository and store reports of the analyses in the `.tortellini/out/`
directory.

```shell
export INPUT_REPOSITORIES=
export INPUT_CURATIONS=''
export INPUT_CLASSIFICATIONS=https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/license-classifications.yml
export INPUT_RULES=https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/rules.kts
npm install
npm run build
npm run package
node dist/index.js
```

### On other repositories

You can also analyze other repositories on Github by storing their addresses in
a file and running node on the file, e.g.:

```shell
echo 'https://github.com/tortellini-tools/action' > urls.txt
echo 'https://github.com/fair-software/howfairis' >> urls.txt
```

The analysis expects a few environment variables. Here are their names and suggested values:

```
export INPUT_REPOSITORIES=urls.txt
export INPUT_CURATIONS=''
export INPUT_CLASSIFICATIONS=https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/license-classifications.yml
export INPUT_RULES=https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/rules.kts
sudo rm -r .tortellini
mkdir .tortellini
node dist/index.js
```

The analyses will be stored in the directories
`.tortellini/out/<owner>/<repository>/` .

## How to create a release

1. Update the citation metadata in `CITATION.cff`. Afterwards, follow the instructions from the `cffconvert` workfow to sync the information in `.zenodo.json` with that in `CITATION.cff`
1. Actions are run from GitHub repos so we need to generate the Javascript files in the `dist` folder and push the results:

    ```bash
    $ cd $(mktemp --directory --tmpdir tortellini-prep-release.XXXXXX)
    $ git clone https://github.com/tortellini-tools/action .
    $ npm install
    $ npm run all
    $ git add dist
    $ git commit --message "prod dependencies"
    $ git push origin main
    ```

1. Next, check if the workflows of the lastest commit on the main branch are green on the [action page](https://github.com/tortellini-tools/action/actions?query=branch%3Amain).
1. Create a release on the Github page via
   [Create a new release](https://github.com/tortellini-tools/action/releases/new).
1. On the new release page, for `Tag version` use `v` and the next version number, for example `v3`.
   See the [versioning documentation](https://github.com/actions/toolkit/blob/master/docs/action-versioning.md)
   for more information.
1. Make sure that usage workflows are using the new version tag and the examples in README.md are updated.

Your action is now published! :rocket:

Check if the new version has been published on the [Github Marketplace](https://github.com/marketplace/actions/tortellini-action).

You can now validate the action by going to
[this workflow](https://github.com/tortellini-tools/action/actions/workflows/usage-current-repository.yml) and [this workflow](https://github.com/tortellini-tools/action/actions/workflows/usage-multiple-repositories.yml)
and then clicking on the button `Run workflow`.
<p align="center">
  <a href="https://github.com/tortellini-tools/action/actions"><img alt="typescript-action status" src="https://github.com/tortellini-tools/action/workflows/build-test/badge.svg"></a>
  <a href="https://github.com/tortellini-tools/action/actions"><img alt="linting-action status" src="https://github.com/tortellini-tools/action/workflows/linting/badge.svg"></a>
  <a href="https://github.com/tortellini-tools/action/actions/workflows/usage-current-repository.yml"><img alt="tortellini-action status" src="https://github.com/tortellini-tools/action/actions/workflows/usage-current-repository.yml/badge.svg"></a>
  <a href="https://github.com/tortellini-tools/action/actions/workflows/usage-multiple-repositories.yml"><img alt="tortellini-action status" src="https://github.com/tortellini-tools/action/actions/workflows/usage-multiple-repositories.yml/badge.svg"></a>
  <a href="https://doi.org/10.5281/zenodo.4956072"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4956072.svg" alt="DOI"></a>
  <a href="https://www.research-software.nl/software/tortellini-github-action">
  <img src="https://img.shields.io/badge/rsd-tortellini-00a3e3.svg" alt="Research Software Directory"></a>
</p>

# Tortellini GitHub Action

This action checks dependency license issues using [ort](https://github.com/oss-review-toolkit/ort).

This GitHub action can
* Run license analysis on your repository
* Run license analysis on list of repositories given
* Detect licensing violations
* Summarize potential licensing issues
* Generate a report with summary

## Inputs

### `repositories`

A file containing list of GitHub repository urls. Format is a single url (`https://github.com/<owner>/<repo>`) on each line. If set then action will run check on each url and generate an overview HTML page.
If not set the action runs check on currently checked out repository (`.` directory).
By default the `repositories` input is not set.

### `rules`

**Required** A file or URL containing [ort rules](https://github.com/oss-review-toolkit/ort/blob/master/docs/file-rules-kts.md) to detect license violations. Default is [https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/rules.kts](https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/rules.kts).

### `classifications`

**Required** A file or URL containing classes for each license. For format see [ort documentation](https://github.com/oss-review-toolkit/ort/blob/master/docs/config-file-license-classifications-yml.md). Default is [https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/license-classifications.yml](https://github.com/NLeSC/tortellini-on-rsd/raw/main/config/license-classifications.yml).

### `curations`

A file or URL containing curations correct invalid or missing package metadata and set the concluded license for packages. See [ort documentation](https://github.com/oss-review-toolkit/ort/blob/master/docs/config-file-curations-yml.md) for format. If not set then action will not use any curations.

## Outputs

No action outputs are defined for this actions.
The action will write files to `.tortellini/out/` directory.

## Usage

### Own repository

```yaml
on:
    # Allows you to run this workflow manually from the Actions tab
    workflow_dispatch:

jobs:
    tortellini:
        runs-on: ubuntu-latest
        steps:
            - uses: actions/checkout@v2
            - uses: tortellini-tools/action@v3
            - uses: actions/upload-artifact@v2
              with:
                  name: tortellini-result
                  path: .tortellini/out
```

Tortellini action will generate `.tortellini/out/scan-report-web-app.html` file.
The HTML file can be downloaded from the workflow page (see the documentation of the [`Upload a Build Artifact` GitHub action](https://github.com/actions/upload-artifact#where-does-the-upload-go)). After unzipping the `scan-report-web-app.html` can be viewed in a web browser.

### Multiple repositories

```yaml
on:
    # Allows you to run this workflow manually from the Actions tab
    workflow_dispatch:

jobs:
    tortellini:
        runs-on: ubuntu-latest
        steps:
            - name: Create list of GitHub urls to perform check on
              run: |
                  echo 'https://github.com/tortellini-tools/action' > urls.txt
                  echo 'https://github.com/fair-software/howfairis' >> urls.txt
            - uses: tortellini-tools/action@v3
              with:
                  repositories: urls.txt
            - uses: actions/upload-artifact@v2
              with:
                  name: tortellini-results
                  path: .tortellini/out
```

Tortellini action will generate `.tortellini/out/index.html` file and `.tortellini/out/<GitHub user or organization>/<GitHub repository>/scan-report-web-app.html` files.

The HTML files can be downloaded from the workflow page (see the documentation of the [`Upload a Build Artifact` GitHub action](https://github.com/actions/upload-artifact#where-does-the-upload-go)). After unzipping the `index.html` can be viewed in a web browser.

## Developer documentation

See [README.dev.md](README.dev.md)

## Disclaimer

`tortellini` aims at providing insights into the license depencies of your package based on available information on open source licenses. We hope this information is helpful. However, we are not lawyers and we do make mistakes. Therefore `tortellini` provides information on an "as-is" basis and does not make warranties regarding this information. For any questions regarding the issues related to the licenses in your code, please consult a professional.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
s.verhoeven@esciencecenter.nl.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/tortellini-tools/action/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/tortellini-tools/action/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of concensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest main commit. While working on your feature branch, make sure to stay up to date with the main branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. see [developer documentation](README.dev.md) for hints how to make changes;
1. make sure the existing unit tests still work by running ``npm test``;
1. make sure that no lint errors are present work by running ``npm lint``;
1. add your own unit tests and integration tests (if necessary);
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the this repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request (have a look at some of our old pull requests to see how this works).
# Design

This document describes the goal, requirements and plan of this project.

## Goal

Gain insight in potential problems regarding the licensing of software packages we develop at eScienceCenter

## Purpose

This tool should produce a report about the license compatibility for 1 or more pieces of software and their dependencies.

## Background

To use open source software, a license needs to be applied to it.
Most software have dependencies which themselves have licenses.
Not all licenses can be combined with each other.
It is often difficult to find out the licenses of the dependencies, and to gain insight into whether there are any conflicting licenses in use in a software package.

## Personas

At the Netherlands eScience Center (NLeSC) there is a need for more attention to licensing issues. Personas in the center are

-   program manager, someone accountable for software adhering to guidelines and best practices.
-   research software engineer, someone writing software

## User stories

### Does this project have software with license problems

As a program manager at NLeSC I would like have a tool that can find potential license problems in the software created in a project.

### Does my piece of software have potential license problems

As a research software engineer at NLeSC I would like to find out if my piece of software has any license violations.

## Theoretical steps

Regardless of whether we choose to use an existing tool/service or make something ourselves, the solution needs to have the following elements:

1. input is a list of URLs to our repositories on GitHub, e.g. from Research Software Directory
1. visit each, look in the root of the repository for a list of runtime (?) dependencies, for example from dependency files such as

    1. requirements.txt
    1. setup.cfg
    1. setup.py
    1. Pipfile
    1. Pipfile.lock
    1. pyproject.toml
    1. environment.yml
    1. yarn.lock
    1. package.json
    1. package-lock.json
    1. DESCRIPTION

    when repositories contain subdirectories with any of these dependency files, we need to do more nested/recursing evaluation

1. for each of the project dependencies,
    1. figure out whose copy is being used (e.g. PyPI, GitHub, crates.io, npm, maven, Anaconda cloud, bower, etc. more [here](https://en.wikipedia.org/wiki/List_of_software_package_management_systems#Application-level_package_managers))
    1. identify the license as stated on the platform whose copy we're using
1. yields a one-to-many graph (a tree) with relations ("If I have a GPL-3.0 here, can I have Apache-2.0 one level higher?")
1. figure out if there is any license of a dependency that yields a conflict with the top level license (should be Apache-2.0) (based on which rules?)
1. choose which result best represents the situation:
    1. didn't find a dependency file
    1. found a dependency file,
        1. every name was resolvable
            1. no conflicts
            1. had conflicts (with list of conflicts found)
        1. not every name was resolvable
            1. needs human evaluation

## High level implementation plan

After [evaluating existing tools](https://github.com/tortellini-tools/action/issues/2) we decided to wrap [ort](https://github.com/oss-review-toolkit/ort) in a script. As ort can detect licenses of dependences of a freshly cloned repository and produce a nice HTML report and machine readable reports. Ort is also widely used and supported by the OSS community. The runner up was [https://github.com/pivotal/LicenseFinder](https://github.com/pivotal/LicenseFinder), but it requires that dependencies are already installed while ort installs them for you.

1. Write a script to regularly
    - Clone the repository
    - For each repo generate a report using [ort](https://github.com/oss-review-toolkit/ort)
2. Define curation files
    - Start with some predefined curations
    - If users want to change or overwrite these curations they can
3. Add index page to list repos
    - Each repo should have a detailed report (html webapp)
4. Instructions for program managers and engineers
    - How to run analysis if you are a program manager
    - How to run analysis if you are an engineer
    - How to update curations

## Technical implementation details

The tool will be written as a GitHub action so it can be called in a GitHub workflow
The action has 2 modes:

### List of repos

This mode targets program managers.

Steps in a Github workflow:

1. Weekly scheduled (eg [fairtally-test](https://github.com/jmaassen/fairtally-test/blob/main/.github/workflows/fairtally.yml))
2. Action fetches repos from RSD
3. Use ort Docker image from https://hub.docker.com/r/philipssoftware/ort/
4. Run https://github.com/tortellini-tools/action/blob/rsd-software-vs-ort/ort/batch.sh, replace shell script with Typescript
5. Create `index-<timestamp>.json` with stats of all repos using Typescript
6. Create symlink/copy `index-<timestamp>.json` to index-latest.json
7. To S3 Upload index-latest.json, `index-<timestamp>.json` and for each repo
    - scan-report-web-app.html,
    - out.txt
    - evaluation-result.yml
8. Vue app index.html which shows index-latest.json in table
9. Upload files to S3, the S3 bucket can be visted with web browser

Steps 2 .. 8 could be captured in Github Actions like

```yaml
on:
    schedule:
        - cron: '0 0 * * 4' # Every thursday

    # Allows you to run this workflow manually from the Actions tab
    workflow_dispatch:
jobs:
    tortellini:
        runs-on: ubuntu-latest
        steps:
            - name: Get the data dump from the RSD
              run: curl https://research-software.nl/api/software > software.json

            - name: Extract the list of URLs
              run: cat software.json | jq -r '[.[].repositoryURLs.github] | flatten | .[]' > urls.txt

            - name: Run ort on urls
              uses: tortellini-tools/action@v3
              with:
                  repositories: urls.txt
                  output-dir: results
                  curations: conf/curations.yml
                  rules: conf/rules.kts
                  classifications: conf/license-classifications.yml
            - uses: jakejarvis/s3-sync-action@master
              with:
                  args: --acl public-read --follow-symlinks
              env:
                  AWS_S3_BUCKET: ${{ secrets.AWS_S3_BUCKET }}
                  AWS_ACCESS_KEY_ID: ${{ secrets.AWS_ACCESS_KEY_ID }}
                  AWS_SECRET_ACCESS_KEY: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
                  AWS_REGION: 'us-west-1' # optional: defaults to us-east-1
                  SOURCE_DIR: 'results'
```

### Single repo

This mode targets engineers.

Using the GitHub Action for a single repo is very similar to using it for multiple repositories, but works on the currently checked out repository instead of on a list of URLs.

```yaml
on:
    release:
    # Allows you to run this workflow manually from the Actions tab
    workflow_dispatch:
jobs:
    tortellini:
        runs-on: ubuntu-latest
        steps:
            - uses: actions/checkout@v2
            - name: Run ort on .
              uses: tortellini-tools/action@v3
            - uses: actions/upload-artifact@v2
              with:
                  name: tortellini-results
                  path: results/**
```

## Action structure

Action repo will have:

-   action.yml, GitHub Action definition file
-   src/index.ts which will clone repos, start ort container, save output to repos/<OWNER>/<REPO>, render index-latest.json
-   package.json with docker library
-   README/LICENSE

Repo in NLeSC org with ort default config files

Start repo from boilerplate <https://github.com/actions/typescript-action>
Use same repo as a test
<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

    While the checklist is intended to be filled in by the
    reviewers, it is the responsibility of the author of the pull request to make
    sure all items on it are properly implemented.
-->

## Description

<!--
    Please describe your changes here, especially focusing on why this pull request makes tortellini better and what problem it solves.

    Before you start, please read our contribution guidelines: https://github.com/tortellini-tools/action/blob/main/CONTRIBUTING.md

    Please fill in the GitHub issue that is closed by this pull request, e.g. Refs #1903
-->

Refs: #issue_number

* * *

## Before you get started

<!--
    Please discuss your idea with the development team before getting started,
    to avoid disappointment or unnecessary work later. The way to do this is
    to open a new issue on GitHub.
-->

- [ ] I read and I followed [contributing guidelines](https://github.com/tortellini-tools/action/blob/main/CONTRIBUTING.md)
- [ ] [☝ Create an issue](https://github.com/NLeSC/licenseguard/issues/new/choose) to discuss what you are going to do

## Checklist

This section should be filled in by the creator of this pull request to make sure the pull request is ready to review.

- [ ] This pull request has a descriptive title
- [ ] Code is written according to the code quality guidelines
- [ ] Documentation if applicable is available
- [ ] Tests run successfully
- [ ] All checks below this pull request were successful

## Instructions to review the pull request

<!--
   Please describe how to test and what the expected behavior is.
-->
