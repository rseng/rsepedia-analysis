# Project Lead:

    - Vanessa Sochat <vsochat@stanford.edu> (@vsoch)


# Contributors:

    - David Trudgian <dave@trudgian.net> (@dctrud)
    - Victor <victorsv@gmail.com> (@victorsndvg)
# CHANGELOG

This is a manually generated log to track changes to the repository for each release. 
Each section should include general headers such as **Implemented enhancements** 
and **Merged pull requests**. All closed issued and bug fixes should be 
represented by the pull requests that fixed them. Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes


## [master](https://github.com/singularityhub/sregistry/tree/master) (master)
 - add: auto set "verify" attribute of s3 and s3_external obj in minio.py for SSL use (1.1.39) 
 - fix issues with psycopg2-binary and saml auth (1.1.38)
   - Pin psycopg2-binary 2.8.6 to deal with UTC errors
   - change formatting of the login URL to fix saml auth
 - bump django version update (1.1.37)
 - update base container to Python 3.6.13 (1.1.36)
 - fix google build deprecated djangorestframework function
 - add notes in docs/docker-compose.yaml to pin versions
 - adding GitHub enterprise backend for social auth (1.1.35)
 - remove un-needed lib PyYaml (1.1.34)
 - updating Django and Django Restframework (1.1.33)
 - API endpoint to create a collection (1.1.32)
 - allowing for Bearer token to have any casing (1.1.31)
 - adding minio environment file to https docker-compose (1.1.3)
 - enforcing usernames to be all lowercase (1.1.29)
 - Added ability to specify Minio direct download from interface (1.1.28)
 - Adding cleanup for Minio images no longer referenced in sregistry (1.1.27)
 - Django various updates, version bump 2.2.10 to 2.2.13 (1.1.26)
 - Minio environment needs to be added to scheduler and worker (1.1.25)
 - Adding minio backend for library endpoint (1.1.24)
 - GitHub API is deprecating use of GET parameters, must provide token in header (1.1.23)
 - multipart upload added to scs-library-client, needs to return 404 (1.1.22)
 - fixed bug push of new container (for same tag) does not update binary (1.1.21)
 - adding POSTGRES_HOST_AUTH_METHOD trust to account for changes postgres 9.6.17 (1.1.20)
 - bumping django version to fix two CVEs (1.1.19)
 - pinning verison of Django to not yet upgrade (1.1.18)
 - broken API and documentation links (1.1.17)
 - refactored collections treemap to only show collection container counts (1.1.16)
 - adding logs for cron jobs and fix their execution (1.1.15)
 - black linting and removing default upload_id (1.1.14)
 - push/pull fixes for library API (1.1.13)
 - bug fix with parsing valid token (1.1.12)
 - support for robot users (1.1.11)
 - library endpoint with push/auth (1.1.10)
 - adding key server (1.1.09)
 - fixing bug with select for team permissions (1.1.08)
 - adding backup as cron job (1.1.07)
 - collection settings are viewable by registry staff/superusers (1.1.06)
 - library pull needs to minimally check if container is private (1.1.05)
 - accidental removed of user profile prefix for custom profile (1.1.04)
 - adding django-ratelimit to all views, customized via settings (1.1.03)
   - button in profile to delete account
   - API throttle with defaults for users and anon
 - setting API schema backend to use coreapi.AutoSchema (1.1.02)
 - documentation fixes
 - application migrations added back in run_uwsgi.sh (1.1.01)
   - major update of documentation theme
 - addition of Google Cloud Build, versioning, tags to collections (1.1.0)
 - adding BitBucket authentication backend
 - updating sregistry-cli to 0.0.97, catching OSError earlier
 - updating sregistry-cli to 0.0.96, and Singularity download url to use sylabs organization
 - increasing length of name limit to 500, and catching error (with message and cleanup)
 - adding Globus integration
 - updating sregistry-cli to version 0.0.74
 - superusers and admins (global) can now create collections via a button in the interface
 - demos and customizated supplementary content removed - not used
 - user customization possible by superusers in the admin panel
 - adding teams and basic permissions to view and edit collections
 - changed internal functions to use sregistry client (instead of singularity)
 - authentication check added to sregistry pull, so private images can be pulled given correct credentials
 - updating the token to have format with registry at top level of dictionary (to support other sregistry clients).
 - from the *sregistry client* provided by Singularity Python, to support use of squashfs images and singularity 2.4, the default upload is not compressed, assuming squashfs, and the default download is not decompressed. To still compress an image add the `--compress` flag on push, and the `--decompress` flag on pull.
 - the generation of the date used for the credential has been fixed, done via updating singularity-python.
# Singularity Registry Server

[![status](https://joss.theoj.org/papers/050362b7e7691d2a5d0ebed8251bc01e/status.svg)](http://joss.theoj.org/papers/050362b7e7691d2a5d0ebed8251bc01e)
[![GitHub actions status](https://github.com/singularityhub/sregistry/workflows/sregistry-ci/badge.svg?branch=master)](https://github.com/singularityhub/sregistry/actions?query=branch%3Amaster+workflow%3Asregistry-ci)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1012531.svg)](https://doi.org/10.5281/zenodo.1012531)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-20-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

- [Documentation](https://singularityhub.github.io/sregistry)

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://vsoch.github.io"><img src="https://avatars0.githubusercontent.com/u/814322?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Vanessasaurus</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=vsoch" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=vsoch" title="Code">💻</a></td>
    <td align="center"><a href="tschoonj.github.io"><img src="https://avatars0.githubusercontent.com/u/65736?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tom Schoonjans</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=tschoonj" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=tschoonj" title="Code">💻</a></td>
    <td align="center"><a href="antoinecully.com"><img src="https://avatars3.githubusercontent.com/u/6448924?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Antoine Cully</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=Aneoshun" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=Aneoshun" title="Code">💻</a></td>
    <td align="center"><a href="https://dctrud.sdf.org"><img src="https://avatars1.githubusercontent.com/u/4522799?v=4?s=100" width="100px;" alt=""/><br /><sub><b>David Trudgian</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=dctrud" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/serlophug"><img src="https://avatars3.githubusercontent.com/u/20574493?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Sergio López Huguet</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=serlophug" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=serlophug" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/jbd"><img src="https://avatars2.githubusercontent.com/u/169483?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jbd</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=jbd" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=jbd" title="Code">💻</a></td>
    <td align="center"><a href="http://alex.hirzel.us/"><img src="https://avatars3.githubusercontent.com/u/324152?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Alex Hirzel</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=alhirzel" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=alhirzel" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://tangiblecomputationalbiology.blogspot.com"><img src="https://avatars0.githubusercontent.com/u/207407?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Steffen Möller</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=smoe" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=smoe" title="Code">💻</a></td>
    <td align="center"><a href="www.onerussian.com"><img src="https://avatars3.githubusercontent.com/u/39889?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Yaroslav Halchenko</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=yarikoptic" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=yarikoptic" title="Code">💻</a></td>
    <td align="center"><a href="http://sourceforge.net/u/victorsndvg/profile/"><img src="https://avatars3.githubusercontent.com/u/6474985?v=4?s=100" width="100px;" alt=""/><br /><sub><b>victorsndvg</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=victorsndvg" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=victorsndvg" title="Code">💻</a></td>
    <td align="center"><a href="arfon.org"><img src="https://avatars1.githubusercontent.com/u/4483?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Arfon Smith</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=arfon" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=arfon" title="Code">💻</a></td>
    <td align="center"><a href="https://ransomwareroundup.com"><img src="https://avatars3.githubusercontent.com/u/9367754?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brie Carranza</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=bbbbbrie" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=bbbbbrie" title="Code">💻</a></td>
    <td align="center"><a href="https://orcid.org/0000-0002-6178-3585"><img src="https://avatars1.githubusercontent.com/u/145659?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Fornika</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=dfornika" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=dfornika" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/RonaldEnsing"><img src="https://avatars2.githubusercontent.com/u/8299064?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ronald Ensing</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=RonaldEnsing" title="Documentation">📖</a> <a href="https://github.com/singularityhub/sregistry/commits?author=RonaldEnsing" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/vladdoster"><img src="https://avatars.githubusercontent.com/u/10052309?v=4?s=100" width="100px;" alt=""/><br /><sub><b>vladdoster</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=vladdoster" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/pini-gh"><img src="https://avatars.githubusercontent.com/u/1241814?v=4?s=100" width="100px;" alt=""/><br /><sub><b>pini-gh</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=pini-gh" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/0nebody"><img src="https://avatars.githubusercontent.com/u/26727168?v=4?s=100" width="100px;" alt=""/><br /><sub><b>0nebody</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=0nebody" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/dtrudg"><img src="https://avatars.githubusercontent.com/u/4522799?v=4?s=100" width="100px;" alt=""/><br /><sub><b>dtrudg</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=dtrudg" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/craigwindell"><img src="https://avatars.githubusercontent.com/u/44250868?v=4?s=100" width="100px;" alt=""/><br /><sub><b>craigwindell</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=craigwindell" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/hashkeks"><img src="https://avatars.githubusercontent.com/u/34633191?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Cedric</b></sub></a><br /><a href="https://github.com/singularityhub/sregistry/commits?author=hashkeks" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

## What is Singularity Registry

Singularity Registry Server is a server to provide management and storage of 
Singularity images for an institution or user to deploy locally. 
It does not manage building but serves endpoints to obtain and save containers. 

## Images Included

Singularity Registry consists of several Docker images, and they are integrated 
to work together using [docker-compose.yml](docker-compose.yml).

The images are the following:

 - **vanessa/sregistry**: is the main uWSGI application, which serves a Django (python-based) application.
 - **nginx**: pronounced (engine-X) is the webserver. The starter application is configured for HTTP. However, you should follow our [instructions](https://singularityhub.github.io/sregistry/docs/install/server#ssl) to set up HTTPS properly. Note that we build a custom NGINX image that takes advantage of the [nginx-upload-module](https://www.nginx.com/resources/wiki/modules/upload/).
 - **worker**: is the same uWSGI image, but with a running command for queueing jobs and processing them in the background. These jobs run via [django-rq](https://github.com/rq/django-rq) backed by a
 - **redis**: database to organize the jobs themselves.
 - **scheduler** jobs can be scheduled using the scheduler.

For more information about Singularity Registry Server, please reference the
[docs](https://singularityhub.github.io/sregistry). If you have any issues,
please [let me know](https://github.com/singularityhub/sregistry/issues)!

## License

This code is licensed under the MPL 2.0 [LICENSE](LICENSE).
---
title: 'Singularity Registry: Open Source Registry for Singularity Images'
tags:
  - containers
  - singularity
  - linux
  - registry
authors:
 - name: Vanessa Sochat
   orcid: 0000-0002-4387-3819
   affiliation: 1
affiliations:
 - name: Stanford University Research Computing
   index: 1
date: 26 September 2017
bibliography: paper.bib
---

# Summary

Singularity Registry is a non-centralized free and Open Source infrastructure to facilitate management and sharing of institutional or personal containers.

A container is the encapsulation of an entire computational environment that can be run consistently if the platform supports it. It is an aid in reproducibility [@Moreews2015-dy, @Belmann2015-eb, @Boettiger2014-cz, @Santana-Perez2015-wo, @Wandell2015-yt] because different researchers can run the exact same software stack on different underlying (Linux Intel) systems. Docker [@Merkel2014-da] has become popular as a general container system because it allows software to be bundled with administrator privileges, however,  it poses great security risks if installed on a multi-user, shared computational resource, and thus Docker has not had admission to these high performance computing (HPC) environments. Singularity [@Kurtzer2017-xj] offers similar features to Docker for software deployment but does not not allow the user to escalate to root, and so it has been easy to accept in HPC environments. Since its introduction, Singularity has been deployed in over 50 HPC resources across the globe.

## Sharing of Containers for Reproducible Science
Essential to the success of Singularity is not just creation of images, but sharing of them. To address this need, a free cloud service, Singularity Hub [@SingularityHub], was developed to build and share containers for scientists simply by way of building containers from a specification file in a version controlled Github repository. This setup is ideal given a small number of containers, each belonging in one Github repository, but is not optimal for institutions that wanted to build at scale using custom strategies.

Singularity Registry (sregistry) [@SingularityRegistryGithub] was developed to empower an institution or individual to build at scale, and push images that can be private or public to their own hosted registry.  It is the first user and institution installable, non-centralized Open Source infrastructure to faciliate the sharing of containers. While Singularity Hub works with cloud builders and object storage,  Singularity Registry is optimized for storage on a local filesystem and any choice of builder (e.g., continuous integration ([Travis](https://www.travis-ci.org), [Circle](https://circleci.com/)), cluster or private node, or separate server). A Registry is also customized with a center's name, and links to appropriate help contacts.

![Singularity Registry Home](registry-home.png)

The Registry, along with native integration into the Singularity software, includes several tools for organization, analysis, and logging of image metrics and usage. Administrators can control the ability for users or the larger community to create accounts, and give finely tuned access (e.g., an expiring token) to share containers. An application programming interface (API) exposes metadata such as container sizes, versioning, and build times. Every time a container is used, or starred by a user, the Registry keeps a record. This kind of metadata not only about containers but about their usage is highly useful to get feedback about highly used containers and general container use.

Public images are available for programmatic usage for anyone, and private images to authenticated users with the Singularity command line software. Sregistry allows for numerous authentication backends, tracks downloads and starring of images, tagging and versioning, search, and provides an interactive treemap visualization to assess size of image collections relative to one another.

![Singularity Collection Sizes](sizes.png)

Importantly, behind sregistry is a growing and thriving community of scientists, high performance computing (HPC) admins, and research software engineers that are incentivized to generate and share reproducible containers. Complete documentation including setup, deployment, and usage, is available [@SRegistry], and the developers welcome contribution and feedback in any form. Singularity Registry empowers the larger scientific community to build reproducible containers on their cluster or local resource, push them securely to the application, and share them toward transparency and reproducibility for discovery in science.


# References
# Backup

This README contains information about backing up containers and the database.

## Containers

If you use Singularity Registry Server with the default storing containers
on the filesystem, your containers will be preserved on the host in [images](../images).
This means that you can restore the registry with previous images, given that
this folder is on the host (bound at `/code/images` in the container).

## Database

This directory will be populated with a nightly backup, both for the previous
day and the one before that. It's run via a cron job that is scheduled in the
main uwsgi container:

```bash
RUN echo "0 1 * * * /bin/bash /code/scripts/backup_db.sh" >> /code/cronjob
```

It's a fairly simple strategy that will minimally allow you to restore your
registry database given that you accidentally remove and then recreate the db container.
For example, here we've started the containers, created some collections, uploaded
containers, and then we (manually) create a backup:

```bash
$ docker exec -it sregistry_uwsgi_1_210f9e7fc042 bash
$ /bin/bash /code/scripts/backup_db.sh
$ exit
```

Then stop and remove your containers.

```bash
$ docker-compose stop
$ docker-compose rm
```

Bring up the containers, and check the interface if you like to verify your
previous collection is gone. Again shell into the container,
but this time restore data.

```bash
$ docker exec -it sregistry_uwsgi_1_a5f868c10aa3 bash
$ /bin/bash /code/scripts/restore_db.sh
Loading table users
Installed 1 object(s) from 1 fixture(s)
Loading table api
Installed 1 object(s) from 1 fixture(s)
Loading table main
Installed 2 object(s) from 1 fixture(s)
Loading table logs
Installed 0 object(s) from 1 fixture(s)
```

If you browse to the interface, your collections should then be restored.
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
reported by contacting the project team leader @vsoch. All
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
# Contributor's Agreement
This code is licensed under the MPL 2.0 [LICENSE](LICENSE).

# Contributing

When contributing to Singularity Registry, it is important to properly communicate the
gist of the contribution. If it is a simple code or editorial fix, simply
explaining this within the GitHub Pull Request (PR) will suffice. But if this
is a larger fix or Enhancement, it should be first discussed with the project
leader or developers.

Please note we have a code of conduct, described below. Please follow it in
all your interactions with the project members and users.

## Pull Request Process

1. Bug fix PRs should be sent to both the master and development branches.
   Feature enhancements should only be submitted against the development
   branch.
2. Follow the existing code style precedent. We use black for linting, and you
   can install it `pip install black` and then run `make lint` to format the code.
   This will be tested.
3. Test your PR locally, and provide the steps necessary to test for the
   reviewers.
4. The project's default copyright and header have been included in any new
   source files.
5. All (major) changes to Singularity Registry must be documented in
   [docs](docs). If your PR changes a core functionality, please 
   include clear description of the changes in your PR so that the docs 
   can be updated, or better, submit another PR to update the docs directly.
6. If necessary, update the README.md.
7. The pull request will be reviewed by others, and the final merge must be
   done by the Singularity project lead, @vsoch (or approved by her).


# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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

### Our Responsibilities

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
reported by contacting the project leader @vsoch. All
complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. The project
team is obligated to maintain confidentiality with regard to the reporter of
an incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers, contributors and users who do not follow or enforce the
Code of Conduct in good faith may face temporary or permanent repercussions 
with their involvement in the project as determined by the project's leader(s).

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
---
name: Request for Documentation or Tutorial
about: request improvements or changes to docs or a tutorial

---


---
name: Bug report
about: Describe a bug you want fixed

---

**Describe the bug**

**To Reproduce**

**Expected behavior**

If applicable, add versions and screenshots to help explain your problem.
---
name: Feature request
about: Suggest an idea for this project

---

**Is your feature request related to a problem? Please describe.**

**Describe the solution you'd like**

**Describe alternatives you've considered**

Anything else?
---
name: Question
about: What's on your mind?

---


---
layout: registry
title:  "{{ REGISTRY_NAME }}"
base: {{ DOMAIN_NAME }}
date:   {{ DATETIME_NOW }}
author: {{ AUTHOR }}
categories:
- registry
img: robots/robot{{ NUMBER }}.png
thumb: robots/robot{{ NUMBER }}.png # wget https://vsoch.github.io/robots/assets/img/robots/robot15570.png
tagged: {{ REGISTRY_URI }}
institution: {{ REGISTRY_NAME }}
---

{{ REGISTRY_NAME }} is a Singularity Registry to provide institution level Singularity containers.
# Library Module

This module is dependent on the Containers and Collections models defined
in [shub/main](../main). It is kept separate from the original API in the
case that the two ever need to be separated.

## Development

There is no documentation for the library API, so I created it based on
what I could see from the interface.
---
title: Singularity Registry Server
permalink: /
excluded_in_search: true
---

<div style="float:right; margin-bottom:50px; color:#666">
</div>

<div>
    <img src="assets/img/logo.png" style="float:left">
</div><br><br>


# Singularity Registry Server

Hello there! It's so great that you want to use Singularity Registry Server. Let's get started. 

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/singularityhub/lobby)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00426/status.svg)](https://doi.org/10.21105/joss.00426)

<br>

## Introduction
This section covers rationale, background, and frequently asked questions.

 - [Introduction](docs/introduction): Covers some background and basic information.
 - [Use Cases](docs/use-cases): a few examples of when deploying a Singularity Registry Server is useful
 - [Frequenty Asked Questions](docs/faq): Quick answers to some questions you might have on your mind.

## Deploy a Registry Server
This section is going to cover installation, which means configuration of a Docker image, building of the image, and then starting your image with other services (docker-compose) to run your registry server.

 - [install](docs/install): configure, build, and deploy your registry server.
 - [setup](docs/setup): setting up and registering the running application.
 - [accounts](docs/accounts/credentials): User accounts, teams, and credentials for using the client.
 - [plugins](docs/plugins): Details about available plugins, and how to configure them.

## Use a Registry

 - [Interface](docs/interface): interacting with your collections in the interface
 - [sregistry Client](docs/client): The `sregistry-cli` to push, pull, list, and delete.
 - [Singularity Client](docs/singularity-client): A client provided by Singularity natively to pull.

Do you want a new feature? Is something not working as it should? @vsoch wants [your input!](https://www.github.com/singularityhub/sregistry/issues) This registry is driven by desire and need for clusters small and large, so if you don't tell me what you want, how will I know that you want it?

The robots thank you for giving them spirit and energy to provide this for you!
# Changelog

This is a manually generated log to track changes to the repository for each release. 
Each section should include general headers such as **Implemented enhancements** 
and **Merged pull requests**. All closed issued and bug fixes should be 
represented by the pull requests that fixed them.
Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes
 - migration guidance
 - changed behaviour

## [master](https://github.com/vsoch/mkdocs-jekyll/tree/master)
 - getting search working (0.0.1)
 - start of theme  (0.0.0)
# MkDocs Jekyll Theme

[![CircleCI](https://circleci.com/gh/vsoch/mkdocs-jekyll/tree/master.svg?style=svg)](https://circleci.com/gh/vsoch/mkdocs-jekyll/tree/master)

![assets/img/mkdocs-jekyll.png](assets/img/mkdocs-jekyll.png)

This is a [starter template](https://vsoch.github.com/mkdocs-jekyll/) for a mkdocs jekyll theme, based on these two
previous arts:

 - [alexcarpenter/material-jekyll-theme](http://alexcarpenter.github.io/material-jekyll-theme)
 - [squidfunk/mkdocs-material](https://github.com/squidfunk/mkdocs-material)

## Usage

### 1. Get the code

You can clone the repository right to where you want to host the docs:

```bash
git clone https://github.com/vsoch/mkdocs-jekyll.git docs
cd docs
```

### 2. Customize

To edit configuration values, customize the [_config.yml](_config.yml).
To add pages, write them into the [pages](pages) folder. 
You define urls based on the `permalink` attribute in your pages,
and then add them to the navigation by adding to the content of [_data/toc.myl](_data/toc.yml).

### 3. Options

Most of the configuration values in the [_config.yml](_config.yml) are self explanatory,
and for more details, see the [about page](https://vsoch.github.io/mkdocs-jekyll/about/)
rendered on the site.

### 4. Serve

Depending on how you installed jekyll:

```bash
jekyll serve
# or
bundle exec jekyll serve
```
---
title: About
permalink: /about/
---

# About

Singularity Registry Server is an open source registry for Singularity
containers. To read more about it's background, see {% include doc.html name="the introduction" path="introduction" %}.
For policy, including license, and suggested usage agreement, see [the policy page]({{ site.baseurl }}/policy).

## What infrastructure is needed?

Singularity Registry Server needs a web accessible server that can install and run Docker and Docker Compose. I (@vsoch) originally started developing it within a Singularity container, but decided to user Docker. Why?

 1. Feedback was strongly that Docker would be ok for such an application. More institutions are able to support separate servers for applications like this.
 2. Docker is great and optimized for orchestration of services. Singularity is not close to that, and with the default image type being an (unchangeable) squashfs, it does not seem to be moving in a direction to be optimized for service containers.
 3. Since images are pushed and pulled via `PUT` and `POST`, it isn't the case that the registry needs to be part of the cluster itself. We don't have the same security concerns as we do with running containers.

The components of the application include databases, a web server, worker, and application:

 - **quay.io/vanessa/sregistry**: is the main uwsgi application, which serves a Django (python-based) application.
 - **nginx**: pronounced (engine-X) is the webserver. The starter application is configured for http, however you should follow the instructions to set up https properly. Note that we build a custom nginx image that takes advantage of the [nginx upload module](https://www.nginx.com/resources/wiki/modules/upload/).
 - **worker**: is the same uwsgi image, but with a running command that is specialized to perform tasks. The tasks are run via [django-rq](https://github.com/rq/django-rq) that uses a
 - **redis**: database to organize the jobs themselves.
 - **scheduler** jobs can be scheduled using the scheduler.

This means that, given a pretty basic server to run the application, and enough space connected to it to store the images, you can bring the entire thing up relatively quickly. Awesome! Let's get started and talk about first steps of [install]({{ site.baseurl }}/docs/install). Or read about [use cases first]({{ site.baseurl }}/docs/use-cases)

## Support

If you need help, please don't hesitate to [open an issue](https://www.github.com/{{ site.github_repo }}/{{ site.github_user }}).

---
layout: page
title: Articles
permalink: /archive/
---
# News Archive

{% for post in site.posts  %}{% capture this_year %}{{ post.date | date: "%Y" }}{% endcapture %}{% capture next_year %}{{ post.previous.date | date: "%Y" }}{% endcapture %}

{% if forloop.first %}<h2 class="c-archives__year" id="{{ this_year }}-ref">{{this_year}}</h2>
<ul class="c-archives__list">{% endif %}
<li class="c-archives__item">
  {{ post.date | date: "%b %-d, %Y" }}: <a href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
  </li>{% if forloop.last %}</ul>{% else %}{% if this_year != next_year %}
</ul>
<h2 class="c-archives__year" id="{{ next_year }}-ref">{{next_year}}</h2>
<ul class="c-archives__list">{% endif %}{% endif %}{% endfor %}
---
title: Policy
permalink: /policy/
---

# Policy

## License

Usage of Singularity Registry Server is governed by the terms of the [license](https://github.com/singularityhub/sregistry/blob/master/LICENSE), which you should read and agree to before using the software.

## Usage Agreement

When deployed as a service, the following usage agreement is provided by any 
Singularity Registry server, which requires the user to agree with the terms,
and a record is kept on the server. If you do not wish to use these terms,
you are free to change the link under `shub/base/templates/terms/usage-agreement.html`.
It is also recommended that you ensure your registry is [GDPR Compliant](https://en.wikipedia.org/wiki/General_Data_Protection_Regulation).

<div class="content">
   <iframe frameBorder="0" width="100%" height="900px" src="https://docs.google.com/document/d/1np84Pd36nPWQGnrN2ZsH7mID__fyeBllSMjQ5vc9Z3s/pub?embedded=true"></iframe>
</div>

---
title: News
permalink: /news/
---

# News

<p>Subscribe with <a href="{{ site.baseurl }}/feed.xml">RSS</a> to keep up with the latest news.
For site changes, see the <a href="https://github.com/{{ site.github_user }}/{{ site.github_repo }}/blob/master/CHANGELOG.md">changelog</a> kept with the code base.</p>

{% for post in site.posts limit:10 %}
   <div class="post-preview">
   <h2><a href="{{ site.baseurl }}{{ post.url }}">{{ post.title }}</a></h2>
   <span class="post-date">{{ post.date | date: "%B %d, %Y" }}</span><br>
   {% if post.badges %}{% for badge in post.badges %}<span class="badge badge-{{ badge.type }}">{{ badge.tag }}</span>{% endfor %}{% endif %}
   {{ post.content | split:'<!--more-->' | first }}
   {% if post.content contains '<!--more-->' %}
      <a href="{{ site.baseurl }}{{ post.url }}">read more</a>
   {% endif %}
   </div>
   <hr>
{% endfor %}

Want to see more? See the <a href="{{ site.baseurl }}/archive/">News Archive</a>.
---
title:  "Rate Limits Added with 1.1.03"
date:   2019-08-20 01:51:21 -0500
categories: security update
badges:
 - type: info
   tag: security
 - type: success
   tag: update
---

This is an announcement for the 1.0.03 release of Singularity Registry Server!

## Updates

### Collection and Container Limits

We can never anticipate a malicious user making requests to an API or views for a server,
so we've added limits for pulling containers, customizable by the registry.

<!--more-->

### View and API Rate Limits

We've also added [view rate limits](/sregistry/docs/install/settings#view-rate-limits)
and API throttling as another component to customize your registry server. You can be
stringent or flexible in how you allow your users to interact with the server.

### Global Disables

In the case that you need to turn a component off to investigate or debug, we've 
added global variables to disabling building, receiving builds, or GitHub webhooks (in
the case of the [Google Build](/sregistry/docs/plugins/google-build) integration.

### Account Control

If users want to delete their account, there is now the the addition of a "Delete Account" button in the User profile
Deleting an account corresponds with deleting all associated containers and collections.

### Private Signed Urls

Previously, we served public containers from Google Storage for the Google Build plugin.
This can be risky if a URL is shared, and then maliciously used to incur charges for
the egress. To better control access, containers are now kept private in storage and accessed via signed urls.
---
title: collection interface
pdf: true
toc: false
---

# Interfaces

## Home

When you browse to your registry's root, it's a fairly expected interface with
links to collections and documentation. Here is what we see for the default "Tacosaurus Computing Center."

![teams.png](../assets/img/interface/home.png)

When you sign in with social authentication (e.g., GitHub), you first must sign the standard terms of service (TOS).
The TOS asserts that the software is provided for you, free to use, and you are responsible
for taking care of your registry and are accountable for the containers within it.

![teams.png](../assets/img/interface/terms.png)

When you agree, the robots welcome you!

![teams.png](../assets/img/interface/terms-welcome.png)

## Collections 

A collection is a set of containers under the same namespace. For example, `dinosaur/avocado` and `dinosaur/banana` could be two containers in the `dinosaur` collection. You can browse all collections in the main collection view:

![collections.png](../assets/img/collections.png)

and browse the containers within a collection by clicking on it. Here is a collection with one container:

![collection.png](../assets/img/interface/collection.png)

And a collection with multiple.

![collection.png](../assets/img/collection.png)

Each collection also has usage instructions.

![usage.png](../assets/img/interface/usage.png)

### Add a Container

A view has been added for users with permission to upload a container to a collection directly! You
can do this by clicking the "+" in the menu above the container table. here is what the upload page looks like:

![upload.png](../assets/img/upload.png)

This uploads directly to NGINX via the [nginx-upload-module](https://www.nginx.com/resources/wiki/modules/upload/), so it should be pretty speedy.

## User Profile

From the user settings menu, you can quickly grab your user token, which will allow you
to do authenticated pushes using Singularity.

![teams.png](../assets/img/interface/token.png)

You can also quickly access your user profile, where you can see your starred collections, a table of your collections,

![teams.png](../assets/img/interface/profile.png)

## Collection Settings

The most important control panel for your collections is the Settings page, which we
see as the first link in the menu at the top of the table. The most likely action you will want to do is add other users from your teams (as described above). You do this on each Collection settings page:

![team-settings.png](../assets/img/team-settings.png)

For example, if my lab has a set of users on sregistry and we intend to build images together, we would make the team for our lab and easily find one another to manage access to images.

### Badges

Recently added, you can get a badge to link to your collection.

![assets/img/badges.png](../assets/img/badges.png)

### Users

You might want to give other users control of your collection (to push and pull and generally manage), and these are called **owners**. You might also want to give some users pull access, most relevant if your collection is private. You can do that in the "Contributors" tab of the settings page:

![assets/img/team-settings.png](../assets/img/team-settings.png)

Remember that you can only choose to add individuals that are part of one of your teams. This means that you should generally make the team first.

### Danger Zone

And of course, if you need to delete, the settings page has a Danger Zone. Be careful!

![assets/img/danger.png](../assets/img/danger.png)


## Application Programming Interface

The application programming interface (API) used by the Singularity and sregistry clients
is available for browsing via two views. A standard one:

![teams.png](../assets/img/interface/api.png)

And a Swagger schema one, if that's how you roll:

![teams.png](../assets/img/interface/swagger.png)

## Tools

Along with the API, there are tools to visualize your registry and otherwise interact.

![teams.png](../assets/img/interface/tools.png) 

For example, here is the treemap.

![teams.png](../assets/img/container_treemap.png)

## Teams

Singularity Registry Server allows registry staff (and if the administrators allow it) authenticated users to create teams or groups of users that want to collaborate on container collections together.

![teams.png](../assets/img/teams.png)

If you are allowed to create and manage teams (see the setup page section about [teams](/sregistry/setup#teams) for information about this), the team permission level determines how others are added to the team.  If a team is **open**, then anyone can join. If it's **invite** only, then you need to generate an invitation. To do this, you can navigate to your Team page and click the button to "Invite User":

![team-invite.png](../assets/img/team-invite.png)

The interface will give you a link to send to your colleague to invite them to the team. Once used, it will expire.

![team-invite-link.png](../assets/img/team-invite-link.png)

Membership in teams is essential because when you add another user as a collaborator to one of your collections (either an owner or member), they must be part of one of your teams.

### Admin Control of Teams

While the Singularity Registry server doesn't directly allow administrators to add any users to be part of a
collection of contributors or owners (this we believe should be up to the collection owners), it is possible to do this
programmatically if it's necessary. Here is an example:

```python
# $ docker exec -it sregistry_uwsgi_1 bash
# python manage.py shell
from shub.apps.main.models import Collection
from shub.apps.users.models import User

# Target a collection and user
collection = Collection.objects.get(name="collection")
user = User.objects.get(username='vsoch')

# Who are current contributors, owners?
collection.contributors.all()                                                                                                          
collection.owners.all()      

# Add to contributors and/or owners
collection.contributors.add(user)
collection.owners.add(user)
```

It is the philosophy of this developer that collection owners should be responsible for this,
and thus the "Teams" feature is advocated for use as it explicitly states, "I am creating
this team of trusted users to add to my collection."
---
title: Singularity Client
pdf: true
toc: false
---

# Singularity Pull

Support for Singularity Registry Server pull via the Singularity software was added to the development branch in  [this pull request](https://github.com/singularityware/singularity/pull/889), and will be in the release of Singularity 2.4 [demo](https://asciinema.org/a/134694?speed=3).

## Pull

and then, given an image called "tacos" in collection "dinosaur" with tag "latest" (`shub://127.0.0.1/dinosaur/tacos:latest`) and a registry serving on your localhost (`127.0.0.1`) you would pull the image using Singularity as follows:

```
singularity pull shub://127.0.0.1/dinosaur/tacos:latest
```

### How is it different? 

This is different than the traditional shub pull command (which uses Github as a namespace) because of the inclusion of the registry uri. On Github, if you imagine that there was a repo named "tacos" owned by user "dinosaur" that pull would look different (and still work):

```
singularity pull shub://dinosaur/tacos:latest
```
---
title: Introduction
pdf: true
toc: false
---

# Introduction

## A Need for Reproducible Science

Using computational methods to answer scientific questions of interest is an important task to increase our knowledge about the world. Along with careful assembly of protocol and relevant datasets, the scientist must also write software to perform the analysis, and use the software in combination with data to answer the question of interest. When a result of interest to the larger community is found, the scientist writes it up for publication in a scientific journal. This is what we might call a single scientific result.

Replication of a result would increase our confidence in the finding. The extent to which a published finding affords a second scientist to repeat the steps to achieve the result is called reproducibility. Reproducibility, in that it allows for repeated testing of an interesting question to validate knowledge about the world, is a foundation of science. While the original research can be an arduous task, often the culmination of years of work and commitment, attempts to reproduce a series of methods to assess if the finding replicates is equally challenging. The researcher must minimally have enough documentation to describe the original data and detailed methods to put together an equivalent experiment. A comparable computational environment must then be used to look for evidence to assert or reject the original hypothesis.

Unfortunately, many scientists are not able to provide the minimum product to allow others to reproduce their work. It could be an issue of time - the modern scientist is burdened with writing grants, managing staff, and fighting for tenure. It could be an issue of education. Graduate school training is heavily focused on a particular domain of interest, and developing skills to learn to program, use version control, and test is outside the scope of the program. It might also be entirely infeasible. If the experiments were run on a particular supercomputer and/or with a custom software stack, it is a non trivial task to provide that environment to others. The inability to easily share environments and software serves as a direct threat to scientific reproducibility.

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00426/status.svg)](https://doi.org/10.21105/joss.00426)

{% include alert.html title="Citation" type="info" content="Sochat, (2017), Singularity Registry: Open Source Registry for Singularity Images<br>Journal of Open Source Software, 2(18), 426, doi:10.21105/joss.00426" %}

## Encapsulation of Environments with Containers

The idea that the entire software stack, including libraries and specific versions of dependencies, could be put into a container and shared offered promise to help this problem. Linux containers, which can be thought of as encapsulated environments that host their own operating systems, software, and flie contents, were a deal breaker when coming onto the scene in early 2015. Like Virtual Machines, containers can make it easy to run a newer software stack on an older host, or to package up all the necessary software to run a scientific experiment, and have confidence that when sharing the container, it will run without a hitch. In early 2015, an early player on the scene, an enterprise container solution called Docker, started to be embraced by the scientific community. Docker containers were ideal for enterprise deployments, but posed huge security hazards if installed on a shared resource.

It wasn't until the introduction of the Singularity software that these workflows could be securely deployed on local cluster resources. For the first time, scientists could package up all of the software and libraries needed for their research, and deliver a complete package for a second scientist to reproduce the work. Singularity took the high performance computing world by storm, securing several awards and press releases, and within a year being installed at over 45 super computing centers across the globe. 

## Singularity Registry Server

Singularity Registry Server is a Dockerized web application that an institution or individual can deploy to organize and manage Singularity images. After you [install](install) and [setup](setup) your registry, you are welcomed with the home screen. In this case, our institution is called "Tacosaurus Computing Center":

<br>
![../assets/img/registry-home.png](../assets/img/registry-home.png)

<br>
You can log in via the social backends that you've configured, in this case, the default is Twitter because it has the easiest setup:
<br>

![../assets/img/login.png](../assets/img/login.png)

<br>
And your registry "About" page is specific to your group, meaning a customized contact email and help link:
<br>

![../assets/img/about.png](../assets/img/about.png)

<br>
And you can quickly glimpse at the names, links, and relative sizes for all containers in the registry:
<br>

![../assets/img/sizes.png](../assets/img/sizes.png)
Enough screen shots! Let's get familiar first with some of the basics.


## How are images named?
When you deploy your registry, it lives at a web address. On your personal computer, this would be localhost (127.0.0.1), and on an institution server, it could be given its own domain or subdomain (eg, `containers.stanford.edu`). This means that, for example, if you had a container called `science/rocks` with tag `latest`, and if you wanted to pull it using the [Singularity software](singularity-client.md), the command would be:

```bash
$ singularity pull shub://127.0.0.1/science/rocks:latest
```

If you use the [sregistry software](/sregistry/client) (the main controller that is configured for a specific registry) then you don't need to use the domain, or the `shub://` uri.

```bash
$ sregistry pull science/rocks:latest
```

The name space of the uri (e.g., `/science/rocks:latest`) is completely up to you to manage. Here are a few suggestions for a larger cluster or institution:

 - `[ cluster ]/[ project ]`
 - `[ group ]/[ project ]`
 - `[ user ]/[ project ]`

For a personal user, you could use software categories or topics:

 - `[ category ]/[ software ]`
 - `[ neuroimaging ]/[ realign ]`

Singularity Hub, based on its connection with Github, uses `[ username ]/[ reponame ]`. If you manage repositories equivalently, you might also consider this as an idea. The one constaint on naming is that only the special character `-` is allowed, and all letters are automatically made lowercase. There are fewer bugs that way, trust us.

## How are images shared?

Akin to Singularity Hub, you share your images by making your registry publicly accessible (or some images in it) and then others can easily download your images with the pull command.

```bash
$ singularity pull shub://127.0.0.1/science/rocks:latest # localhost
```

You can also generate an expiring link for a user to download the image equivalently:
<br>

![../assets/img/share.png](../assets/img/share.png)


## Terms
Let's now talk about some commonly used terms.

### registry
The registry refers to this entire application. When you set up your registry, you will fill out some basic information in settings, and send it to Singularity Hub. When we have a few registries running, we will have a central location that uses endpoints served by each to make images easily findable.


### collections
Each container image (eg, `shub://fishman/snacks`) is actually a set of images called a `collection`. This is the view looking at all collections in a registry:

<br>
![../assets/img/collection.png](../assets/img/collection.png)

<br>
Within a collection you might have different tags or versions for images. For example:

 - `milkshake/banana:pudding`
 - `milkshake/chocolate:pudding`
 - `milkshake/vanilla:pudding`

All of these images are derivations of `milkshake`, and so we find them in the same collection. I chose above to change the name of the container and maintain the same tag, you could of course have more granular detail, different versions of the same container:

 - `milkshake/banana:v1.0`
 - `milkshake/banana:v2.0`

### Labels

A `label` is an important piece of metadata (such as version, creator, or build variables) that is carried with a container to help understand it's generation, and organize it. When you create a Singularity image, you might bootstrap from Docker, and or add a "labels" section to your build definition:

```bash
%labels
MAINTAINER vanessasaur
```

The registry automatically parses these labels, and makes them searchable for the user.
<br>
![../assets/img/labels.png](../assets/img/labels.png)
<br>
When you investigate an individual label, you can see all containers in the registry that have it! For example, here I reveal that my testing images are in fact the same image named differently, or at least they have two unique sizes:
<br>
![../assets/img/label_single.png](../assets/img/label_single.png)

### Topic Tags

Topic tags are a way for users to (after build, within the Registry) add "topic words" to describe containers. They are searchable in the same way that labels are:
<br>
![../assets/img/topic_tags.png](../assets/img/topic_tags.png)
<br>
While a label may be something relevant to the build environment, a topic tag is more something like "biology" or the operating system. For example, if you look at the single container view, the Singularity Registry automatically parses the "From" statement to create a topic for each operating system type:
<br>
![../assets/img/single_tag.png](../assets/img/single_tag.png)
<br>
Notably, if you are logged in, you can dynamically click and write in a new tag for the container, and it is automatically saved.


### Favorites
Do you have a favorite collection? You can star it! Each Singularity Registry keeps track of the number of downloads (pull) of the containers, along with stars! Here we see the number of stars for our small registry:
<br>
![../assets/img/stars.png](../assets/img/stars.png)
<br>

and equally, container downloads:
<br>
![../assets/img/downloads.png](../assets/img/downloads.png)
<br>
If there is other metadata you would like to see about usage, please let me know.
---
title: Singularity Clients
pdf: true
toc: false
---

# Singularity

## Singularity Pull

Singularity Registry Server implements a basic version of the Sylabs Library API, 
meaning that you can pull a container with Singularity directly. 

**Important** you must be using Singularity 3.3.0 or greater for this to work! 
If not, you should use [Singularity Registry Client](#singularity-registry-client)
 (which also has examples below of using the `shub://` endpoint with singularity).

For example, let's say that I have a collection with a container called `collection/container:tag`. and my registry is served at `containers.page`. I could pull it as follows:

```bash
$ singularity pull --library https://containers.page collection/container:tag
```

You can also pull a container using Singularity natively with the `shub://` uri:

```bash
$ singularity pull shub://containers.page/collection/container:tag
```

## Singularity Push

As of version `1.1.10`, Singularity Registry Server offers a library endpoint
to authenticate and then push containers. First, create an endpoint for your
registry:

```bash
$ singularity remote add --no-login DinosaurCloud cloud.dinosaur.io
```

Verify it's there:

```bash
$ singularity remote list
NAME           URI              	GLOBAL
DinosaurCloud  cloud.dinosaur.io        NO
[SylabsCloud]  cloud.sylabs.io   	YES
```

**Important** Sylabs requires these endpoints to have https, for obvious reasons. If you want to test with localhost, you'll need to edit [this file](https://github.com/sylabs/singularity/blob/5e483be4af2e120e646d33f0757e855c8d3be2da/internal/pkg/remote/remote.go#L237)
and re-compile Singularity to set the url to have http. The example above is a hypothetical "cloud.dinosaur.io" however in actual development, I used 127.0.0.1 (and you'll see it 
in various spots below). This is how I developed this set of endpoints.

Once you add the remote, then you'll first need to login and get your token at the /token endpoint, for example:

```bash
1eb5bc1daeca0f5a215ec242c9690209ca0b3d71
```

And then provide it (via copy paste) to the Singularity client to create a remote endpoint for your registry:

```bash
$ singularity remote login DinosaurCloud
INFO:    Authenticating with remote: DinosaurCloud
Generate an API Key at https://127.0.0.1/auth/tokens, and paste here:
API Key: 
INFO:    API Key Verified!
```

If you paste a token that isn't valid, you'll get a different message

```bash
$ singularity remote login DinosaurCloud
INFO:    Authenticating with remote: DinosaurCloud
Generate an API Key at https://127.0.0.1/auth/tokens, and paste here:
API Key: 
FATAL:   while verifying token: error response from server: Invalid Token
```

In case you are wondering, the token is kept in plaintext at `/home/vanessa/.singularity/remote.yaml`
so once you specify to use an endpoint, it knows the token. If you are having
issues copy pasting the token into your terminal (I had some when I wanted to
re-create it) you can also just open up this file and edit the text manually:

```
$ cat /home/vanessa/.singularity/remote.yaml Active: DinosaurCloud
Remotes:
  DinosaurCloud:
    URI: 127.0.0.1
    Token: 1eb5bcrdaeca0f5a215ef242c9690209ca0b3d71
    System: false
  SylabsCloud:
    URI: cloud.sylabs.io
    Token: hahhaayeahrightdude
    System: true
```

The easiest thing to do is now to set your remote to be DinosaurCloud (or whatever
you called it) so you don't need to specify the push command with `--library`:

```bash
singularity remote use DinosaurCloud
```

Now that we have a token, let's try a push! For security purposes, the collection
should already exist, and be owned by you. Collaborators are not allowed to push.

```bash 
                                         # library://user/collection/container[:tag]
$ singularity push -U busybox_latest.sif library://vsoch/dinosaur-collection/another:latest
```

If you don't do "remote use" then you can specify the library on the fly:

```bash
$ singularity push -U --library http://127.0.0.1 busybox_latest.sif library://vsoch/dinosaur-collection/another:latest
WARNING: Skipping container verifying
INFO:    http://127.0.0.1
INFO:    0a4cb168d3dabc3c21a15476be7f4a90396bc2c1
INFO:    library://vsoch/dinosaur-collection/another:latest
INFO:    [latest]
 656.93 KiB / 656.93 KiB [===================================================================] 100.00% 15.36 MiB/s 0s
```

We use `-U` for unsigned.

### Push Logic

 - If you push an existing tag, if the container is unfrozen, it will be replaced
 - If you push an existing tag and the container is frozen (akin to protected) you'll get permission denied.
 - If you push a new tag, it will be added.
 - If you push a new image, it will also be added.
 - If you push an image to a non existing collection, the collection will be created first, then the image will be added (version `1.1.32`).

Unlike the Sylabs API, when the GET endpoint is made to `v1/containers` and the image doesn't exist,
we return a response for the collection (and not 404). In other words, [this response](https://github.com/sylabs/scs-library-client/blob/acb520c8fe6456e4223af6fbece956449d790c79/client/push.go#L140) is always returned. We do this because
the Sylabs library client has a strange logic where it doesn't tag images until after the fact,
and doesn't send the user's requested tag to any of the get or creation endpoints. This means
that we are forced on the registry to create a dummy holder tag (that is guaranteed to be unique)
and then to find the container at the end to [set tags](https://github.com/sylabs/scs-library-client/blob/acb520c8fe6456e4223af6fbece956449d790c79/client/push.go#L187) based on the id of the image
that is created with the [upload request](https://github.com/sylabs/scs-library-client/blob/acb520c8fe6456e4223af6fbece956449d790c79/client/push.go#L174). I didn't see a logical way to create the container using the POST endpoint to
"v1/containers" given that we do not know the tag or version, and would need to know the exact container id
to return later when the container push is requested.

### Push Size

The push (as of this version) can now handle large images! Here is the largest that I've tested:

```bash
$ singularity push -U rustarok_latest.sif library://vsoch/dinosaur-collection/rustarok:latest
WARNING: Skipping container verifying
INFO:    http://127.0.0.1
INFO:    0a4cb168d3dabc3c21a15476be7f4a90396bc2c1
INFO:    library://vsoch/dinosaur-collection/rustarok:latest
INFO:    [latest]
 8.09 GiB / 8.09 GiB [====================================================================] 100.00% 91.31 MiB/s 1m30s
```

Of course, do this at your own risk! That is a *CHONKER*!

<hr>
<br>

# Singularity Registry Client

Singularity Registry Global Client, or [sregistry-cli](https://github.com/singularityhub/sregistry-cli),
is a general client to interact with Singularity images at remote endpoints, and it provides
such an endpoint for Singularity Registry Server. We will provide
basic instructions here, and for the full documentation, please see the [getting started guide here](https://singularityhub.github.io/sregistry-cli/client-registry). Note that you will need to [export your credentials](https://singularityhub.github.io/sregistry/credentials) in order to have authenticated interaction with sregistry.


## Install

### sregistry Installation

`sregistry` is the client for Singularity Registry server. To install, you can do the following:

```bash
git clone https://github.com/singularityhub/sregistry-cli
cd sregistry-cli
python setup.py install
```

To check your install, run this command to make sure the `sregistry` client is found.

which sregistry


### Container Install

We have provided a Singularity build definition for you, for which you can use to build a container that serves as the sregistry client (and this will likely be provided on Singularity Hub so you don't even need to do that.) To build, do the following:

```bash
cd sregistry/

# Singularity 2.4 and up
sudo singularity build sregistry Singularity

# For Singularity earlier than 2.4 (deprecated)
singularity create --size 2000 sregistry
sudo singularity bootstrap sregistry Singularity
```

If you install via this option, you will want to make sure the container itself is somewhere on your path, with appropriate permissions for who you want to be able to use it.


## Commands
This brief tutorial assumes that you have [Singularity installed](https://singularityware.github.io/install-linux).

### Pull
Not shown in the demo above is the pull command, but it does the same thing as the singularity pull.

```bash
sregistry pull banana/pudding:milkshake
Progress |===================================| 100.0% 
Success! banana-pudding-milkshake.img
```

This is useful so that you can (locally from your registry) pull an image without needing to specify the registry url. It's also important because registry support will only be added to Singularity when the entire suite of compoenents are ready to go!


### Push

If you don't have an image handy, you can pull one from Singularity Hub:

```bash
singularity pull shub://vsoch/hello-world
```

And then a push to your registry looks like this:

```bash
sregistry push vsoch-hello-world-master.img --name dinosaur/avocado --tag delicious
sregistry push vsoch-hello-world-master.img --name meowmeow/avocado --tag nomnomnom
sregistry push vsoch-hello-world-master.img --name dinosaur/avocado --tag whatinthe
```

If you don't specify a tag, `latest` is used. If you have authentication issues,
remember that you need to [export a token](https://singularityhub.github.io/sregistry/credentials) for your user, and ensure that the user is either an admin/manager, or
that you have set the `USER_COLLECTIONS` variable to true. You can read [more about roles here](https://singularityhub.github.io/sregistry/setup-roles), and [more about teams](https://singularityhub.github.io/sregistry/setup-teams) to manage groups of people.

### List

List is a general command that will show a specific container, a specific collection, optionally with a tag. Examples are provided below:

```bash
# All collections
sregistry list

# A particular collection
sregistry list dinosaur

# A particular container name across collections
sregistry list /avocado

# A named container, no tag
sregistry list dinosaur/avocado

# A named container, with tag
sregistry list dinosaur/avocado:delicious
```

In addition to listing containers, `sregistry` can show you metadata! It does this by issuing an inspect command at upload time, so that no processing is needed on the server side. Singularity Registry is a Dockerized application, so it would require --privileged mode, which is a bad idea. Anyway, we can look at environment (`--env/-e`), runscript (`--runscript/-r`), tests (`--test/-t`), or `Singularity` definition recipe (`--deffile/-d`):

```bash
# Show me environment
sregistry list dinosaur/tacos:delicious --env

# Add runscript
sregistry list dinosaur/tacos:delicious --e --r

# Definition recipe (Singularity) and test
sregistry list dinosaur/tacos:delicious --d --t

# All of them
sregistry list dinosaur/tacos:delicious --e --r --d --t
```

### Delete
Delete requires the same authentication as push, and you will need to confirm with `yes/no`

```bash
sregistry delete dinosaur/tacos:delicious
sregistry list
```

if you want to force it, add `--force`

```bash
sregistry delete dinosaur/tacos:delicious --force
```

### Labels
Labels are important, and so they are represented as objects in the database for index, query, etc. Akin to containers, we can list and search:

```bash
# All labels
sregistry labels

# A specific key
sregistry labels --key maintainer

# A specific value
sregistry labels --value vanessasaur

# A specific key and value
sregistry labels --key maintainer --value vanessasaur
```

# Curl

Like with the Sylabs Library API, it is possible to interact with Singularity Registry Server 
using Curl. You can browse the API schema via the `/api` path of your server.

## Authentication

Authentication is done presenting the token in an `Authorization` header:

```bash
$ curl -s -H 'Authorization: Bearer <token>' http://127.0.0.1/<api_endpoint>
```

The token can be found in the navigation in the top right of the registry interface
after you log in.

## Create a collection

As of version `1.1.32` it is possible to create a new collection via the API. It requires authentication.
First retrieve the numeric `id` associated with your username with a GET request to the endpoint `/v1/entities/<username>`.

```bash
$ curl -s -H 'Authorization: Bearer <token>' /v1/entities/<username>
```
Here is a response made pretty by piping into json_pp:
```
{
   "data" : {
      "collections" : [],
      "createdAt" : "2021-02-21T05:20:18.454003Z",
      "createdBy" : "",
      "customData" : "",
      "defaultPrivate" : false,
      "deleted" : false,
      "deletedAt" : "0001-01-01T00:00:00Z",
      "description" : "vsoch",
      "id" : "1",
      "name" : "vsoch",
      "quota" : 0,
      "size" : 0,
      "updatedAt" : "2021-02-21T05:20:18.479251Z",
      "updatedBy" : ""
   }
}
```

Notice that the id is 1? Great! We will use this to create a collection. We next issue
a POST request to the endpoint `/v1/collections` and this json payload:

```
{
  "entity": "<user_numeric_id>"
  "name": "<new_collection_name>"
  "private": true|false
}
```

Here is an example with our user id of 1:

```bash
$ curl -X POST -H 'Authorization: Bearer <token>' -H "Content-Type: application/json" --data '{"entity": 1, "name": "dinosaurs"}' http://127.0.0.1/v1/collections 
```

You can then see the response that the collection was created, and it will appear in the interface:

```bash
{
   "data" : {
      "containers" : [],
      "createdAt" : "2021-02-21T05:35:36.491446Z",
      "createdBy" : "1",
      "customData" : "",
      "deleted" : false,
      "deletedAt" : "0001-01-01T00:00:00Z",
      "description" : "Dinosaurs Collection",
      "entity" : "1",
      "entityName" : "vsoch",
      "id" : "2",
      "name" : "dinosaurs",
      "owner" : "1",
      "private" : false,
      "size" : 0,
      "updatedAt" : "2021-02-21T05:35:36.505902Z",
      "updatedBy" : "1"
   }
}
```

The `private` key is optional. If not provided, it defaults to the servers's configured default for collection creation.
In case of a `singularity push` to a non existing collection, the client triggers the collection creation first, using this endpoint, then pushes the image.
---
title: Use Cases
pdf: true
toc: false
---

# Use Cases 

Why would you need a Singularity Registry Server? The answer is when you would want to be able to manage, organize, and share images, and you don't want to be limited in the number of builds, or how you build the images.

## Personal Use
In this use case, I am an individual user, or share a computer resource with a small group. I decide to run a local Singularity Registry Server (on localhost or 127.0.0.1) that anyone that uses the computer can push images after being built. The images are then accessible on my computer based on the unique resource identifier, without me needing to keep track of where they were built. I can get a summary of image sizes, and search over images to find images of interest. I keep the entire Registry private so it cannot be discovered externally, but create an expiring sharing link to send securely to a collaborator to test one of my pipelines.

## Managed Cluster Registry
My university runs a shared computational resource manages a registry on a server next to it. Akin to supplying software modules, the administrators keep a version controlled repo of build recipes, and when software needs to be updated, create a new image with a tag for the version. The users can then use the images by way of specifying the unique resource identifier. 

## Collaborative Cluster Registry
It's often the case that pipelines are maintained internally within labs, or eventually discarded after papers are published and graduate students finish. In this use case, a large cluster wants to provide a central, organized resource for the scientific containers generated by its researchers. Perhaps alongside or instead of the core software and tools, this cluster decides to build final or published containers for its users. Building might lead to a private image for use at the institution, or a public image that can be referenced in a publication and easily disseminated. To build, the researcher simply might submit a pull request to a Github repo associated with the registry, it can be built and tested and discussed, and when ready, pushed to the resource from the continuous integration, or by the cluster's particular build server. Either way, the final upload is an authenticated, single line call to push the image with an appropriate name and tag. If you 
add [plugins](plugins) you can also have custom authentication and builds (e.g., GitHub webhooks + Google Cloud Build).

If you are a single user and looking for an image management tool, perhaps to work with images in multiple locations beyond a Singularity Registry, Server then you will be interested in the [Singularity Global Client](https://singularityhub.github.io/sregistry-cli).
---
title: Frequently Asked Questions
pdf: true
toc: false
---

# Frequently Asked Questions

 - [What is Singularity Registry Server](#what-is-singularity-registry-server)
 - [What is a Linux Container?](#what-is-a-linux-container)
 - [How is Singularity Registry Server different from Singularity Hub?](#how-is-singularity-registry-server-different-from-singularity-hub)
 - [Why do we need Singularity containers?](#why-do-we-need-singularity-containers)
 - [Why isn't the storage backed up?](#why-isnt-the-storage-backed-up)
 - [How is a Singularity Registry different from a Docker Registry?](#how-is-a-singularity-registry-different-from-a-docker-registry)
 - [Are there features of Singularity that are particularly supported by Singularity Hub?](#are-there-features-of-singularity-that-are-particularly-supported-by-singularity-hub)


### What is Singularity Registry Server

Singularity Registry Server is an open source registry for [Singularity](https://www.github.com/sylabs/singularity)
containers. It is optimized to be flexible for deployment with different plugins for
authentication, storage, and building, and for community contribution. 

### What is a Linux Container?

A container image is an encapsulated, portable environment that is created to distribute a scientific analysis or a general function. Containers help with reproducibility of such content as they nicely package software and data dependencies, along with libraries that are needed. Thus, the core of Singularity Hub are these Singularity container images, and by way of being on Singularity Hub they can be easily built, updated, referenced with a url for a publication, and shared. This small guide will help you to get started building your containers using Singularity Hub and your Github repositories.

### How is Singularity Hub related to Sylabs and Singularity?

[Sylabs](https://sylabs.io) is the company that maintains the open source [Singularity](https://www.github.com/sylabs/singularity) code base, along with providing services for paying customers that use Singularity containers. The creator and maintainer of Singularity Hub was one of the original Singularity (open source software) developers, however she is not affiliated with the company Sylabs. Singularity Hub continues to be maintained by this individual at Stanford University, with generous support from Google Cloud. The two now co-exist peacefully, both passionate about using and supporting users to build Singularity containers. The distinction between Sylabs Library and Singularity Hub (and [Singularity Registry Server](https://www.github.com/singularityhub/sregistry)) comes down to the intended communities that are served. Singularity Hub and Registry is a non-enterprise solution that is catered for research scientists.


### How is Singularity Registry Server different from Singularity Hub?

**Singularity Hub**

is the predecessor to Singularity Registry, and while it also serves as an image registry, in addition it provides a cloud build service for users. Singularity Hub also takes advantage of Github for version control of build recipes. The user pushes to Github, a builder is deployed, and the image available to the user. Singularity Hub would allow a user to build and run an image from a resource where he or she doesn't have sudo simply by using Github as a middleman.

**Singularity Registry Server** 

is similarly an image registry that plugs in natively to the singularity software, but it places no dependencies on Github, and puts the power of deciding how to build in the hands of the user. This could mean building after tests pass with a "push" command in a Github repository, building via a SLURM job, or on a private server. While Singularity Hub is entirely public and only allows for a minimum number of private images, a Singularity Registry Server can be entirely private, with expiring tokens that can be shared. The administrator can choose
to include Singularity Hub like features via plugins, but only if they are desired.

Importantly, both deliver image manifests that plug seamlessly into the Singularity command-line software, so a registry (or hub) image can be pulled easily:

```bash
singularity pull shub://vsoch/hello-world
```

We are hoping to add builder plugins to Singularity Registry Server so that it can (eventually) replace Singularity Hub
as a truly community maintained, open source resource. Please contribute!


**Singularity Global Client**

`sregistry` is a general client for the single user to interact with containers at different endpoints. An endpoint could be a Singularity Registry Server, Singularity Hub, Google Drive/Cloud, Dropbox, or Globus.  You can think of it as a glue or fabric between all these different endpoints and APIs, as it allows you to create a local database to manage images. For this reason, the executable is called [sregistry](https://singularityhub.github.io/sregistry-cli).

### Why do we need Singularity containers?

Singularity containers allow you to package your entire scientific analysis, including dependencies, libraries, and environment, and run it anywhere. Inside a Singularity container you are the same user as outside the container, so you could not escalate to root and act maliciously on a shared resource. For more information on containers, see [the Singularity site](https://singularityware.github.io).


### Why isn't the storage backed up?

At the initial release of the software, because there are many different options for storage, enforcing a particular backup strategy would possibly make the registry less flexible to fit into many different use cases. In the same way that the institution is able to decide how to build, it is also under their decision for how to backup. For the database, django has different options for backup (for example [django-backup](https://github.com/django-backup/django-backup)), along with a proper mirror (called a [hot standby](https://cloud.google.com/community/tutorials/setting-up-postgres-hot-standby)) of the database itself. An institution might simply want to mirror the filesystem, or to create freezes at consistent timepoints. The [InterPlanetary File System](https://en.wikipedia.org/wiki/InterPlanetary_File_System) has also been suggested, and we hope to have discussion and testing with the larger community to either provide a default or suggest top choices.

### How is a Singularity Registry different from a Docker Registry?

Both are similar, and in fact might be friends! We can easily talk about things they have in common:

 - both serve image manifests that link to relevant image files
 - both can be deployed for a single user or institution
 - both are Dockerized themselves

And things that are different:

 - A Docker registry serves layers (`.tar.gz`) and not entire images. A Singularity Registry serves entire images. 
 - Any Docker image from a Docker Registry can be immediately converted to Singularity. As far as I know, the other way around is not developed. 
 - Singularity Registry images can be pulled to a cluster resource. Docker images cannot.

There are really great use cases for both, and the decision of which to use is up to the goals of the user. 

### Are there features of Singularity that are particularly supported by Singularity Hub?

Singularity aims to support scientific containers, so Singularity Hub and Registry take an extra step to serve metadata about the containers via the API. It's important to know about usage (downloads and stars) but also software, environment variables, labels, and runscripts. This supports being able to do more research analytics across containers to better understand how containers (and more broadly software) help answer scientific questions. Given the issues we have with reproducibility, this is essential.

---
title: Credentials
pdf: true
toc: false
---

# Accounts

### User Accounts

If you remember from the [setup](../setup/teams#create-accounts) you created your first superuser account by giving permission to your account from the command line inside the container. For doing this again, you don't need to interact with the container in this way, but instead, you can manage admins (and other superusers) from the web interface. When logged in to your superuser account, you will see an "Admin" link in your profile in the top right:

![admin-option.png](../../assets/img/admin-option.png)

This will take you to the administrative panel. Once there, you can click on "Users" at the bottom of the list, and then see a list of users to edit, with filters in the right panel.

![admin-users.png](../../assets/img/admin-users.png)

Once you select a user, there will be checkboxes to give staff or superuser status, along with other options.


### Secrets

After you create a user, you will need a way to communicate to the registry, and validate your identity. This can be done by defining the `SREGISTRY_REGISTRY_BASE`, `SREGISTRY_REGISTRY_USERNAME` and `SREGISTRY_REGISTRY_TOKEN` environment variables. Each time we use the client, the secrets is used to encrypt a call and time-specific token that the registry can un-encrypt with the same key, and verify the payload. After creating your account in [setup](../setup/teams#create-accounts), making yourself a superuser and admin and logging in (remember this part?)

```bash
NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
docker exec $NAME python /code/manage.py add_superuser --username vsoch
docker exec $NAME python /code/manage.py add_admin --username vsoch
```

You will want to go to [http://127.0.0.1/token](http://127.0.0.1/token) and use the contents of the json object to define the necessary environment variabes. 

{% include alert.html title="Important!" content="You must be a superuser <strong>and</strong> admin to build images." %}

If you don't add yourself as an admin, the menu looks like this:

![/assets/img/without-superuser.png](../../assets/img/without-superuser.png)

As an admin, you see the button for "token":

![/assets/img/with-superuser.png](../../assets/img/with-superuser.png)

Here is the token page - note the button on the left will copy the text to your clipboard.

![/assets/img/token.png](../../assets/img/token.png)

Define the following variables, by using the corresponding keys from the json object:

```bash
export SREGISTRY_REGISTRY_BASE=http://127.0.0.1
export SREGISTRY_REGISTRY_USERNAME=vsoch
export SREGISTRY_REGISTRY_TOKEN=5027225717bf0030465db1e7099e496022c42181
```

Now when we try to communicate with [the client](/sregistry-cli/client-registry), it finds the token and can identify us.

```bash
export SREGISTRY_CLIENT=registry
sregistry list
No container collections found.
```

Next, see if you are interested in activating any additional [plugins](../plugins) for your Singularity Registry Server.
---
title: Accounts and Credentials
pdf: true
toc: false
permalink: docs/accounts/
---

# Accounts

Singularity Registry server creates user accounts by way of login with OAuth2, and
then interaction with the API via user specific tokens. These sections describe
how to generate accounts, teams, as well as robot users for interacting with 
continuous integration.

 - [credentials](credentials): overview of basic accounts and credentials
 - [robot users](robots) how to create a simple robot user account (admin required)
---
title: Robot Users
pdf: true
toc: false
---

# How to Generate a Robot User

In the case that you want to generate a robot user, or an account not associated with
a real person that can push from some CI service, you'll need to do the following.

{% include alert.html type="info" title="Important!" content="You must be an admin of the server to generate a robot user." %}

 1. Use the Django Administration site to add this user and a token for this user will be automatically created. If you need to refresh the token, you can do so here.
 2. Update the 'last login' field for this user in the Django Administration site to some value, e.g. to the current time.
 3. Since the collection must already exist and you cannot log in as the robot user, you should create the collection with your user account, and then you can add the robot user as an owner (shown below). Only owners are allowed to push to collections.

First, enter the uwsgi container (the name of your container may be different)

```bash
docker exec -it sregistry_uwsgi_1 bash
python manage.py shell
```

Find your robot user, and your collection:

```python
from shub.apps.users.models import User
from shub.apps.main.models import Collection
user = User.objects.get(username="myuser")
collection = Collection.objects.get(name="mycollection")
```

And then add the robot user as an owner to it.

```python
collection.owners.add(user)
collection.save()
```

As an alternative, if you intend to add this robot user to more than one collection,
you can create a [Team]({{ site.url }}{{ site.baseurl }}/docs/setup/teams) in the interface, 
add the robot user to it also via the console:

```python
from shub.apps.users.models import Team, User
team = Team.objects.get(name="myteam")
user = User.objects.get(username="myuser")
team.members.add(user)
```

And then in any collection interface where you can see the team, you can add
the robot user directly as an owner.

And that's it!
---
title: Plugins
pdf: true
toc: false
permalink: docs/plugins/
---

# Plugins

Singularity Registry Server supports added functionality through plugins. Plugins allow complex features,
such as container scanning, LDAP authentication, to be added without complicating the core code of
sregistry.

Plugins distributed with `sregistry` are found in the `shub/plugins` directory. 

## Included Plugins

The following plugins are included with sregistry, and can be enabled by adding them to the
`PLUGINS_ENABLED` entry in `shub/settings/config.py`. Plugins may require further configuration in
your registries' local `shub/settings/secrets.py` file.

 - [LDAP-Auth](ldap): authentication against LDAP directories
 - [PAM-Auth](pam): authentication using PAM (unix host users)
 - [Globus](globus): connect and transfer using Globus
 - [SAML](saml): Authentication with SAML
 - [Google Build](google-build) provides build and storage on Google Cloud.
 - [Keystore](pgp) provides a standard keystore for signing containers

The Dockerfile has some build arguments to build the Docker image according to the plugins software requirements. These variables are set to false by default:

```bash
ARG ENABLE_LDAP=false
ARG ENABLE_PAM=false
ARG ENABLE_GOOGLEBUILD=false
ARG ENABLE_GLOBUS=false
ARG ENABLE_SAML=false
```

Therefore, if you want to install the requirements of all current supported plugins, you can build the image as follows: 
```bash
docker build --build-arg ENABLE_LDAP=true --build-arg ENABLE_PAM=true  --build-arg ENABLE_GOOGLEBUILD=true --build-arg ENABLE_GLOBUS=true --build-arg ENABLE_SAML=true -t quay.io/vanessa/sregistry .
```


## Writing a Plugin

An sregistry plugin is a Django App, that lives inside `shub/plugins/<plugin-name>`.
Each plugin:

 - Must provide a `urls.py` listing any URLs that will be exposed under `/plugin-name`
 - Can provide additional, models, views, templates, static files.
 - Can register an additional `AUTHENTICATION_BACKEND` by specifying `AUTHENTICATION_BACKEND` in
   its `__init.py__`
 - Can register additional context processors by defining a tuple of complete paths to the relevant processors by specifying `CONTEXT_PROCESSORS` in its `__init.py__`
 - Must provide a documentation file and link in this README.

Plugins are loaded when the plugin name is added to `PLUGINS_ENABLED` in `shub/settings/config.py`.
A plugin mentioned here is added to `INSTALLED_APPS` at runtime, and any `AUTHENTICATION_BACKEND`
and `CONTEXT_PROCESSORS` listed in the plugin `__init.py__` is merged into the project settings.

More documentation will be added as the plugin interface is developed. For now, see plugins
distributed with sregisty under `shub/plugins` for example code.

Besides, if your plugin has any specific software requirements that are not currently available in the Docker image and **those requirements are compatible with the current software**, you can set a new build argument `ENABLE_{PLUGIN_NAME}` and add the corresponding installation commands in the `PLUGINS` section of the Dockerfile with the following format:
```bash
RUN if $ENABLE_{PLUGIN_NAME}; then {INSTALLATION_COMMAND}; fi;
```
## Writing Documentation
Documentation for your plugin is just as important as the plugin itself! You should create a subfolder under
`docs/pages/plugins/<your-plugin>` with an appropriate README.md that is linked to in this file.
Use the others as examples to guide you.

---
title: "pgp keystore"
pdf: true
toc: false
permalink: docs/plugins/pgp
---

# PGP KeyStore

The `pgp` plugin adds necessary endpoints for your Singularity Registry server to store
and interact with keys. It uses the [OpenPGP Server](https://tools.ietf.org/html/draft-shaw-openpgp-hkp-00)
protocol, meaning that activating the plugin will expose "lookup" and "add" endpoints.

To enable the pgp plugin you must:

  * Add `pgp` to the `PLUGINS_ENABLED` list in `shub/settings/config.py`
  * Build the docker image with the build argument ENABLE_PGP set to true:
    ```bash
    $ docker build --build-arg ENABLE_PGP=true -t quay.io/vanessa/sregistry .
    ```

The keystore, unlike other plugins, requires no further setup. Brief interactions
for usage are shown below.
  
## Quick Start

Singularity has it's own little key storage in your home, at `$HOME/.singularity/sypgp`. It also
has a set of client functions, `singularity key` to interact with local and remote keys. 
For the client commands we will be using this client.

**Important**: If you are deploying a server, you are required to have https for the
Sylabs Singularity client to interact with the keystore to work. if you don't,
you will see this message:

```
ERROR: push failed: TLS required when auth token provided
```

You can either secure your server with https, or you can test using localhost,
which will work without TLS.

### 1. Create a Key

The first thing that you want to do is likely to create a local key. Again, all of these
interactions from the client side will use the Singularity software.

```bash
$ singularity key newpair
```

Follow the prompts to generate a key, and don't choose yes to send to sylabs cloud (unless you want to).

### 2. List Keys

When you finish, you should be able to list the key:

```bash
$ singularity key list
Public key listing (/home/vanessa/.singularity/sypgp/pgp-public):

0) U: Vanessasaurus (dinosaurs) <vsochat@stanford.edu>
   C: 2019-10-01 12:23:48 -0400 EDT
   F: CFA6763B11637E52404A25F5DE565315F5198C71
   L: 4096
   --------
```

Note that I chose "dinosaurs" as a keyword so it would be searchable by that. The long string (F)
is the unique id.

### 3. Push Key

Now that we have a key, and we know it's identifier, let's push it to our registry server!
For the example below, we have the server running at http://localhost.
This command is exactly as you would do with Singularity, but we need to provide an additional `--url`
to designate the server base.

```bash
$ singularity key push --url http://localhost CFA6763B11637E52404A25F5DE565315F5198C71
public key `CFA6763B11637E52404A25F5DE565315F5198C71' pushed to server successfully
```

### 4. Search for Key

First you can try searching and providing the URL - and I'm not actually sure if this searches locally and remote, but you see the result:

```bash
$ singularity key search --url http://localhost  dinosaur
Showing 1 results

KEY ID    BITS  NAME/EMAIL
f5198c71  4096  Vanessasaurus (dinosaurs) <vsochat@stanford.edu>  

```

Let's delete all of our local keys to verify that we are really searching the Singularity Registry Server:

```bash
$ rm /home/vanessa/.singularity/sypgp/pgp-public 
```

Now let's do the search again. Since it still shows up, it must be from the registry!

```bash
$ singularity key search --url http://localhost  dinosaur
Showing 1 results

KEY ID    BITS  NAME/EMAIL
f5198c71  4096  Vanessasaurus (dinosaurs) <vsochat@stanford.edu>  
```

Now pull it! We can see that 1 key is added.

```bash
$ singularity key pull --url http://localhost CFA6763B11637E52404A25F5DE565315F5198C71
1 key(s) added to keyring of trust /home/vanessa/.singularity/sypgp/pgp-public
```

If you had tried to add the key without deleting the local one first, you already would have had it 
(note below that zero keys are added):

```bash
$ singularity key pull --url http://localhost CFA6763B11637E52404A25F5DE565315F5198C71
0 key(s) added to keyring of trust /home/vanessa/.singularity/sypgp/pgp-public
```
---
title: "Plugin: Custom Builder and Storage"
pdf: true
toc: true
permalink: docs/plugins/google-build
---

# Plugin: Google Cloud Build and Storage

The Singularity Registry client allows for [a large set](https://singularityhub.github.io/sregistry-cli/clients) of options for external storage endpoints. Specifically, this plugin uses storage and build provided by Google, meaning:

 - [Google Build](https://singularityhub.github.io/sregistry-cli/client-google-build)
 - [Google Storage](https://singularityhub.github.io/sregistry-cli/client-google-storage)

Other cloud vendors have been included with sregistry client (AWS, S3, Minio) and equivalent
build and storage pairs can be added here. If you would like to discuss adding a builder
and storage pair, please [open an issue](https://www.github.com/singularityhub/sregistry).

Don't forget to go back to the [install docs](https://singularityhub.github.io/sregistry/install-settings) where you left off. This quick start will walk through setting up custom storage using 
[Google Cloud Build](https://singularityhub.github.io/sregistry-cli/client-google-build)
and [Google Storage](https://singularityhub.github.io/sregistry-cli/client-google-storage) as
an endpoint.

## Configure sregistry

By default, google build is disabled. To configure sregistry to 
use Google Cloud build and Storage, in settings/config.py you can enable the plugin by 
uncommenting it from the list here:

```bash
PLUGINS_ENABLED = [
#    'ldap_auth',
#    'saml_auth',
#    'globus',
     'google_build'
]
```
You will need to build the image locally with, at least, the build argument ENABLE_GOOGLEBUILD set to true:

```bash
$ docker build --build-arg ENABLE_GOOGLEBUILD=true -t quay.io/quay.io/vanessa/sregistry .
```

## Secrets

Next, set the following variables in `shub/settings/secrets.py`, 
that you can create from `dummy_secrets.py` in the shub/settings folder.
The first two speak for themselves, your project name and path to your
Google Application Credentials.

### Project Identifiers

```python
# =============================================================================
# Google Cloud Build + Storage
# Configure a custom builder and storage endpoint
# =============================================================================

# google-storage, s3, google-drive, dropbox
GOOGLE_APPLICATION_CREDENTIALS=/path/to/credentials.json
SREGISTRY_GOOGLE_PROJECT=myproject-ftw

```

You can create custom [Google Application Credentials](https://cloud.google.com/docs/authentication/getting-started) for your server in the browser, and it will be enough to make the service account
a project owner. To allow for signed urls, you will need to also add [iam.serviceAccounts.signBlob](https://cloud.google.com/iam/credentials/reference/rest/v1/projects.serviceAccounts/signBlob?authuser=1)
to the permissions. This is associated with the role "Service Account Token Creator" and I've added
it in the past by going to IAM and Admin -> IAM and then selecting the service account and 
searching for that role. Yes, it's sort of annoying to get working the first time. ;/

If you are on a Google Cloud instance you can scp (with gcloud) using the command line as follows:

```bash
$ gcloud compute scp [credentials].json $USER@[INSTANCE]:/tmp --project [PROJECT]
```

Keep in mind that the path to the Google credentials file must be
within the container (/code is the root folder that is bound to the filesystem).

### Build Caching

```python
SREGISTRY_GOOGLE_BUILD_CACHE="true"
```

If you set this variable (to anything), it means that after build, you will not
delete intermediate dependencies in cloudbuild bucket (keep them as cache for rebuild if needed).
This defaults to being unset, meaning that files are cleaned up. If you define this as anything, 
the build files will be cached.

### Build Limit

```python
SREGISTRY_GOOGLE_BUILD_LIMIT=100
```

To prevent denial of service attacks on Google Cloud Storage, you should
set a reasonable limit for the number of active, concurrent builds. This
number should be based on your expected number of users, repositories, and
recipes per repository.


### Singularity Version

By default, we use the default version that is set by the [Google Build](https://singularityhub.github.io/sregistry-cli/client-google-build#environment) client that belongs to Singularity Registry Client.
However, as this value is subject to be updated, we recommend that you set it in your
secrets and can then decide when you want to update.

```python
SREGISTRY_GOOGLE_BUILD_SINGULARITY_VERSION="v3.2.1-slim"
```

The version must coincide with a container tag hosted under [singularityware/singularity](https://hub.docker.com/r/singularityware/singularity/).

### Storage Bucket Name

By default, the bucket name will be called `sregistry-gcloud-build-[hostname]`, and since
your host is a docker container, that will resolve to a random set of numbers and 
letters. For this reason, we *strongly recommend you set a consistent hostname*.
If you do not and need to remove and bring up the containers again, the bucket
metadata will not match the new bucket name. Here is an example of how to set a custom name:

```python
SREGISTRY_GOOGLE_STORAGE_BUCKET="taco-singularity-registry"
```

Additionally, a temporary bucket is created with the same name ending in `_cloudbuild`. This bucket is for build time dependencies, and is cleaned up after the fact. If you are having trouble getting a bucket it is likely because the name is taken, 
and we recommend creating both `[name]` and `[name]_cloudbuild` in the console and then setting the name here.

### Build Timeout

The number of seconds for the build to timeout. If set to None, will be 10 minutes. If
unset, it will default to 3 hours. This time should be less than the `SREGISTRY_GOOGLE_BUILD_EXPIRE_SECONDS`. If
you want to use the default, don't define this variable in your secrets.

```python
# SREGISTRY_GOOGLE_BUILD_TIMEOUT_SECONDS=None
```

### Signed URL Expiration

By default, the containers are made private, and then granted access via signed urls.
You can optionally adjust the time that the URL will expire in, although it's recommended
that this is kept small:

```python
CONTAINER_SIGNED_URL_EXPIRE_SECONDS=10
```

This can be much smaller than 10, as we only need it to endure for the POST. I've tested and found
that 3-5 seconds is about right.

### Build Expiration 

You must define the number of seconds that your build expires in, meaning that it would no
longer be accepted by the server.

```python
SREGISTRY_GOOGLE_BUILD_EXPIRE_SECONDS=28800
```

The default provided in the dummy secrets, shown above, would indicate 8 hours.


### Disable Github

If you need to globally disable GitHub, meaning that users cannot make new
collections and webhooks are disabled, you can do that:

```python
DISABLE_GITHUB=True
```


### Disable Building

Disable all building, including pushing of containers and recipes. By
default, for a working registry, this should be False.

```python
DISABLE_BUILDING=True
```

This setting is also in the main settings page, but mentioned here for 
Google Cloud Build.

### Disable Build Receive

Prevent responses from being received from Google Cloud Build (returns permission
denied).

```python
DISABLE_BUILD_RECEIVE=True
```

These variables are written in detail in the dummy_secrets.py file. 
If you need more information, you can read [the Google Cloud Build page](https://singularityhub.github.io/sregistry-cli/client-google-build).

If you are missing some variable, there will be an error message
on interaction with the Google Cloud Build API since you won't be able to 
authenticate. Once your settings are ready to go, you will want to continue
with the [install docs](https://singularityhub.github.io/sregistry/install-server#storage) where you left off,
and you can continue here after you've done:

```
$ docker-compose up -d
```

and confirmed the registry running at localhost, and also have logged in
(so you have an account with permission to push containers and recipes.)

## Singularity Registry Client

If you haven't yet, you will need the [sregistry client](https://singularityhub.github.io/sregistry-cli/) in order to push recipes to build with Google Cloud Build. The minimum version that supports this
is `0.2.19`. An easy way to install is any of the following:

```bash
$ pip install sregistry[google-build]
$ pip install sregistry[google-build-basic] # without local sqlite database
```

Next, export the client to be your registry.

```bash
$ export SREGISTRY_CLIENT=registry
```

If you are reading here from the installation docs, you likely haven't
brought up your registry and should [return there](https://singularityhub.github.io/sregistry/install-settings) where you left off.

## Building Containers

There are two ways to trigger builds:

 1. Automated trigger from GitHub webhooks
 2. Manual push of a recipe

The recommended approach is to enable GitHub authentication and then
have pushes to your repository trigger builds. For the second approach,
while you can upload a recipe directly, it is not recommended
as it doesn't have the recipe kept under any version control.

### Trigger from Github

You will first need to log in with GitHub, and then navigate to the
container collections page (the "Containers" link in the navigation):

![/sregistry/assets/img/google-build-new-collection.png](../../assets/img/google-build-new-collection.png)

If the Google Build plugin is correctly enabled, you'll see a second option on the 
right:

![/sregistry/assets/img/google-build-connect-github.png](../../assets/img/google-build-connect-github.png)

Select this, and your repositories (and organizations) that you granted
permission to connect to will show up. You can select one:

![/sregistry/assets/img/google-build-repos.png](../../assets/img/google-build-repos.png)

Once you've connected the repository, an initial build will build
the latest version of the recipes that are discovered. Any recipe that 
is in the format `Singularity.<tag>` or just `Singularity` (tag defaults
to latest) will be built.

![assets/img/google-build-collection.png](../../assets/img/google-build-collection.png)

If you have two recipes named equivalently in different folders, the
recipe that was more recently updated will be used. 

### Push a Recipe

When the server is started and the client is ready, it's time to push a recipe
to build! By default, you will need to specify the name of the collection and
container, and to include the fact that you want to use Google Cloud Build.
You'll need to install Singularity Registry Client version 0.2.21 or later:

````bash
$ pip install sregistry[google-build]>=0.2.21
$ pip install sregistry[google-build-basic]>=0.2.21 # without local database
```

Then to submit a build, you'll need to grab your credentials from https://<registry>/token.
You can write them to your Singularity Registry secrets at `$HOME/.sregistry`. Once your
token and registry base are defined, you will need to create the collection
in the web interface first to establish yourself as an owner. **You cannot
push to a collection that does not exist**. Once the collection is
created (for example, below I created the collection "collection"), you can push like this:

```bash
$ sregistry build --name registry://collection/container:tag Singularity --builder google_build
```

Notice that we specify the builder to be "google_build." Also notice 
that the command simply requires a name for your collection (it doesn't
need to exist, but you need push access and to have [exported your token](https://singularityhub.github.io/sregistry/credentials) to your local machine.

If you get this error:

```bash
[================================] 0/0 MB - 00:00:00
Recipe upload failed: 403 Client Error: Forbidden for url: https://containers.page/google_build/build/.
```

you forgot to create a collection called "collection" and need to make it in the interface before
proceeding.

## Pull Containers

Once you have a container, you of course want to pull it! You can use
the Singularity Client to do this. Let's say that our server is at `https://www.containers.page`:

```bash
$ singularity pull shub://containers.page/singularityhub/hello-registry:latest
 760.00 KiB / 760.00 KiB [=========================================================================================] 100.00% 5.22 MiB/s 0s
```

And there you have it!

```bash
$ ls
hello-registry_latest.sif

$ singularity run hello-registry_latest.sif 
Tacotacotaco!
``` 

Note that having a custom registry name (containers.page, in the above example)
was a bug in early versions of Singularity 3.x. if you have trouble with 
this command, you will need to upgrade Singularity.

You can technically also just pull it with simple bash commands, if you
don't want to rely on Singularity.

```bash
$ wget $(curl https://containers.page/api/container/singularityhub/hello-registry:latest | jq --raw-output .image)
```

If you want to pull with Singularity (but get the error) you can also do this:

```bash
$ singularity pull $(curl https://containers.page/api/container/singularityhub/hello-registry:latest | jq --raw-output .image)
```

Finally, it should be pointed out that you can use the Google Builder integration
from your command line without having a registry at all. [Singularity Registry Client](https://singularityhub.github.io/sregistry-cli/client-google-build) can serve to build and then pull the image on its own.
---
title: "globus - Authentication with Globus"
pdf: true
toc: false
permalink: docs/plugins/globus
---

# Globus Connect

The `globus` plugin allows a logged in user to connect their Globus account to allow for transfer of images from the registry to a Globus endpoint. To use the plugin, you want to take the following steps:


## Setup

In your `shub/settings/secrets.py` file you need to add a client id and secret generated at [https://developers.globus.org/](https://developers.globus.org/). Navigate to the site and do the following:

 - Click on the first option, "Register your app with Globus"
 - In the top right click "Add --> New App"
 - Don't check any of the boxes at the bottom.
 - Choose the following scopes:

```
profile (Know some details about you.)
email (Know your email address.)
openid (Know who you are in Globus.)
urn:globus:auth:scope:transfer.api.globus.org:all (Transfer files using Globus Transfer)
```

And finally, the redirect URIs should include the following, where localhost is appropriate if your container is running on your local machine, and some other uri can be used if hosting on a server with a domain name.

```
http://localhost/globus_auth/login
http://localhost/globus_auth/login/
http://localhost/complete/globus/
http://localhost/complete/globus
http://localhost/globus/login/
http://localhost/globus/login
```
For reference, we are following [these steps](http://globus-sdk-python.readthedocs.io/en/stable/tutorial/#step-1-get-a-client).
Then click "Create app." Once you have the application created, you should copy the client secret and id, and add to your `shub/settings/secrets.py` file like so:

```
SOCIAL_AUTH_GLOBUS_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxx"
SOCIAL_AUTH_GLOBUS_USERNAME="xxxxxxxxxxxxxxxxxxxyyyy@clients.auth.globus.org"
SOCIAL_AUTH_GLOBUS_SECRET="xxxxxxxxxxxxxxx"
```

If you don't yet have a `secrets.py` you can copy the `dummy_secrets.py` in the same folder, which has commented out examples of the above.
Then, you must build the container setting the build argument ENABLE_GLOBUS to true.
```
docker build --build-arg ENABLE_GLOBUS=true -t quay.io/vanessa/sregistry .
```

Once the build is done, start the container.

```
docker-compose up -d
```

and after you've started it, run the script to generate the endpoint (in the example below, the container is named `sregistry_uwsgi_1` and we figured this out with `docker ps`).

```
docker exec -it sregistry_uwsgi_1 /bin/bash /code/scripts/globus/globus-setup.sh
docker-compose restart
```

The script above will ask you to open your browser to authenticate. This first step is with regard to the endpoint. You, as the admin, are
the owner of the endpoint. The user will have to further authenticate from the application interface to interact with it.


## Usage

### Authentication
Once your server is setup with an endpoint, any logged in user must authenticate with Globus in order to issue a refresh token to make transfers. To add this integration, do the following:

 1. When you are logged in, you can access your Integrations under Settings in your user menu:

![pages/plugins/globus/img/settings.png](pages/plugins/globus/img/settings.png)

 2. Go to the "Integrations" tab of the user profile and click "Connect Globus"

![pages/plugins/globus/img/connect-globus.png](pages/plugins/globus/img/connect-globus.png)

 3. You will be redirected to the Globus Login page. Do that.

![pages/plugins/globus/img/globus-login.png](pages/plugins/globus/img/globus-login.png)

 4. And when you login, you will be redirected to the Globus Transfer page:

![pages/plugins/globus/img/transfer.png](pages/plugins/globus/img/transfer.png)


### Search
The Globus Transfer portal shows endpoints with scope `shared-with-me` and 
`my-endpoints`. To see endpoints that you might actually want to transfer to,
you should do a search using the box at the top:

![pages/plugins/globus/img/search.png](pages/plugins/globus/img/search.png)

Here is the search result for the query "Stanford" that shows Stanford's various endpoints.

![pages/plugins/globus/img/search-result.png](pages/plugins/globus/img/search-result.png). 

Any endpoint that says "ACTIVATE" you need to follow the link to activate it from Globus. Any endpoint that has
a green label "activated" should be ready for use.


### Transfer
To do an actual transfer, you can do this directly from your container collections. First, don't forget that you
need to have your [credentials set up](https://singularityhub.github.io/sregistry/credentials#secrets) to push a container
using the [sregistry client](https://singularityhub.github.io/sregistry/client). Let's pull an image from
singularity hub, and push it to our registry:

```bash
singularity pull --name avocados.simg shub://vsoch/hello-world
SREGISTRY_CLIENT=registry sregistry push --name tacos/avocados.simg avocados.simg
```

We can now navigate to our "tacos" collection and see the container! We also see that there is a little transfer
button (the font awesome exchange icon) next to it:

![pages/plugins/globus/img/exchange.png](pages/plugins/globus/img/exchange.png).

This will take you back to the endpoints page, and you are again free to search for an endpoint that you want to
transfer to. This time when you see a green label that indicates an endpoint is activated, you will also see the option
to transfer:

![pages/plugins/globus/img/transfer-option.png](pages/plugins/globus/img/transfer-option.png).

A notification will pop up that the task is underway!

![pages/plugins/globus/img/verify.png](pages/plugins/globus/img/verify.png).

If you click "view task" in the notification you will be taken again to Globus to view the task status. It will hopefully
succeed.

![pages/plugins/globus/img/complete-transfer.png](pages/plugins/globus/img/complete-transfer.png).

Finally, when the time is needed (and if you want to, you don't have to) you can return to
the user profile "Integrations" tab and click "Disconnect Globus" to log out.

![pages/plugins/globus/img/logout.png](pages/plugins/globus/img/logout.png)


## Background

### How does Singularity Registry work with Globus?
Singularity Registry works with Globus by serving its filesystem with images as a read only Globus endpoint. The user authenticates through Globus and is then able to initiate a transfer of any container in the registry to an endpoint they have write permission for, such as one on their local resource.


<div>
    <a href="/sregistry/plugin-pam"><button class="previous-button btn btn-primary"><i class="fa fa-chevron-left"></i> </button></a>
    <a href="/sregistry/plugin-saml"><button class="next-button btn btn-primary"><i class="fa fa-chevron-right"></i> </button></a>
</div><br>
---
title: "ldap-auth - Authentication against LDAP directories"
pdf: true
toc: false
permalink: docs/plugins/ldap
---

# LDAP Authentication

The `ldap_auth` plugin allows users to login to sregistry using account information stored in an
LDAP directory. This supports logins against [Microsoft Active Directory](https://msdn.microsoft.com/en-us/library/bb742424.aspx), as well open-source
[OpenLDAP](https://www.openldap.org/) etc.

To enable LDAP authentication you must:

  * Build the docker image with the build argument ENABLE_LDAP set to true
  * Add `ldap_auth` to the `PLUGINS_ENABLED` list in `shub/settings/config.py`
  * Configure the details of your LDAP directory in `shub/settings/secrets.py`. See
    `shub/settings/dummy_secrets.py` for an example OpenLDAP configuration. A good start is to do the following:

```
cp shub/settings/dummy_secrets.py shub/settings/secrets.py
```
  
Because no two LDAP directories are the same, configuration can be complex and there are no
standard settings. The plugin uses `django-auth-ldap`, which provides more [detailed documentation
at Read the Docs](https://django-auth-ldap.readthedocs.io/en/1.2.x/authentication.html).

To test LDAP authentication you may wish to use a docker container that provides an OpenLDAP
directory. `mwaeckerlin/openldap` [(GitHub)](https://github.com/mwaeckerlin/openldap) [(Docker
Hub)](https://hub.docker.com/r/mwaeckerlin/openldap/) is a useful container configured 
with unencrypted, StartTLS, and SSL access to an OpenLDAP directory.

## Quick Start
This quick start is intended to demonstrate basic functionality of the LDAP server, and you should
review the links referenced above for more detail.

### What is LDAP?

LDAP (Lightweight Directory Access Protocol) is a common protocol used by
organizations to hold user information and perform authentication. An LDAP
directory hold records for each user, and groups which users may belong to.

LDAP directories are implemented by many different directory servers. The most
commonly encountered are OpenLDAP on Linux, and Microsoft Active Directory on Windows platforms.

To test sregistry LDAP authentication we can use a Dockerized OpenLDAP server.


#### Create the server
As instructed in [https://github.com/mwaeckerlin/openldap](https://github.com/mwaeckerlin/openldap) and [(here!)](https://marc.xn--wckerlin-0za.ch/computer/setup-openldap-server-in-docker)  let's bring 
up a dummy LDAP server:

```bash
docker run -d --restart unless-stopped \
              --name openldap \
              -p 389:389 \
              -e DOMAIN=my-company.com \
              -e ORGANIZATION="Tacosaurus" \
              -e PASSWORD=avocados \
              mwaeckerlin/openldap
```

With this command we are:

  - Allowing access on port 389 (unecrypted LDAP / StartTLS encrypted)
  - Creating a directory that will have a `basedn: dc=my-company,dc=com`. The
    basedn is the root of the LDAP directory tree. It is usually created by
    breaking your domain name into domain components (dc).
  - Creating an admin account, which will have the dn (distinguished name)
    `cn=admin,dc=my-company,dc=com` and password `avocados`.

The `-d` means "detached" so you won't see the output in the terminal. If you need to see output, remove the `-d`. Here is the 
running container:

```
docker ps

CONTAINER ID        IMAGE                  COMMAND                  CREATED             STATUS              PORTS                           NAMES
398b6297d6ff        mwaeckerlin/openldap   "/bin/sh -c /start.sh"   3 minutes ago       Up 3 minutes        0.0.0.0:389->389/tcp, 636/tcp   openldap
```

#### Interact with it
Here is a way to get familiar with the executables inside the image for ldap:

```bash
docker exec -it openldap bash

root@docker[72b21bd3c290]:/# which ldapadd
/usr/bin/ldapadd

root@docker[72b21bd3c290]:/# which ldapwhoami
/usr/bin/ldapwhoami
```

Note that the long string with cn= through dc= is your username! The password is the one you set for the image.

```bash
root@docker[4ec2c4f2737a]:/# ldapwhoami -x -D 'cn=admin,dc=my-company,dc=com' -W
Enter LDAP Password:
dn:cn=admin,dc=my-company,dc=com
```

For the password above you would enter the one we set in the environment for the image. In our case this was `avocados`.

#### Add a user and group to the directory

If all has gone well we can check the content of the directory with:

```bash
$ ldapsearch -x -b 'dc=my-company,dc=com'

# extended LDIF
#
# LDAPv3
# base <dc=my-company,dc=com> with scope subtree
# filter: (objectclass=*)
# requesting: ALL
#

# my-company.com
dn: dc=my-company,dc=com
objectClass: top
objectClass: dcObject
objectClass: organization
o: Tacosaurus
dc: my-company

# admin, my-company.com
dn: cn=admin,dc=my-company,dc=com
objectClass: simpleSecurityObject
objectClass: organizationalRole
cn: admin
description: LDAP administrator

# search result
search: 2
result: 0 Success

# numResponses: 3
# numEntries: 2
```

We now need to add some test users and groups to our directory. Create a file
called 'example.ldif' with the following content:

```ldif
##
## Basic LDIF for LDAP simulation setup
## D. C. Trudgian Nov 2014
##

# ----------------------------------------------------------------------------
# STRUCTURE
# ----------------------------------------------------------------------------

dn: ou=users,dc=my-company,dc=com
ou: users
objectClass: top
objectClass: organizationalUnit

dn: ou=groups,dc=my-company,dc=com
ou: Group
objectClass: top
objectClass: organizationalUnit

# ----------------------------------------------------------------------------
# Test Groups
# ----------------------------------------------------------------------------

dn: cn=test,ou=groups,dc=my-company,dc=com
cn: test
description: Test User Group
gidNumber: 1000
objectClass: posixGroup
memberUid: testuser
memberUid: testadmin

dn: cn=admin,ou=groups,dc=my-company,dc=com
cn: admin
description: Test Admin Users
gidNumber: 1002
objectClass: posixGroup
memberUid: testadmin

# ----------------------------------------------------------------------------
# Test Users
# ----------------------------------------------------------------------------

dn: uid=testuser,ou=users,dc=my-company,dc=com
displayName: Test User
cn: Test User
title: Tester
objectClass: inetOrgPerson
objectClass: posixAccount
objectClass: shadowAccount
loginShell: /bin/bash
uidNumber: 1000
gecos: Test User, Test Lab, Test Dept
sn: User
homeDirectory: /home2/testuser
mail: testuser@localhost
givenName: Test
employeeNumber: 1000
shadowExpire: 99999
shadowLastChange: 10000
gidNumber: 1000
uid: testuser
userPassword: {SSHA}IZQSiHPR9A/xKUPTKAM82EJoejtb70vD

dn: uid=testadmin,ou=users,dc=my-company,dc=com
displayName: Test Admin
cn: Test Admin
title: Tester
objectClass: inetOrgPerson
objectClass: posixAccount
objectClass: shadowAccount
loginShell: /bin/bash
uidNumber: 1001
gecos: Test Admin, Test Lab, Test Dept
sn: Admin
homeDirectory: /home2/testadmin
mail: testadmin@localhost
postalAddress: NL5.136
givenName: Test
employeeNumber: 1000
shadowExpire: 99999
shadowLastChange: 10000
gidNumber: 1000
uid: testadmin
userPassword: {SSHA}84O5yFQxQwvQc1Dluc5fJrehrucmCFdH
```

This will create a directory with:

  - Two organizational units (ou=users, ou=groups) to hold our users and groups
  - Two groups, *test* and *admin.*
  - A user *testuser* with password *testuser* who belongs to the *test* group.
  - A user *testadmin* with password *testadmin* who belongs to the *test* and
    *admin* groups.

We use the `ldapadd` command to import this ldif file into our directory (note this will prompt for the password):

```bash
cat example.ldif  | ldapadd -x -H ldap://localhost -D 'cn=admin,dc=my-company,dc=com' -W -v
```

The variables mean the following:

  - `-x` uses simple authentication
  - `-H` specifies the LDAP server to connect to
  - `-D` specifies we want to bind as our admin account
  - `-W` prompts for the password for that account

**Important** You need to get the ip-address of your ldap server. Since we aren't using docker-compose,
the containers won't magically see one another. You can inspect the container's networking as follows:

```bash
docker inspect openldap | grep IPAddress
            "SecondaryIPAddresses": null,
            "IPAddress": "172.17.0.2",
            "IPAddress": "172.17.0.2",
```

The IPAddress thus is `172.17.0.2`. Note that you will need this address in the next step for `AUTH_LDAP_SERVER_URI`.

#### Configure sregistry

To configure sregistry to authenticate against our LDAP directory we need to set
the following options in `shub/settings/secrets.py`:

```python
import ldap
from django_auth_ldap.config import LDAPSearch, PosixGroupType

# The URI to our LDAP server (may be ldap_auth:// or ldaps://)
AUTH_LDAP_SERVER_URI = "ldap://172.17.0.2"

# Any user account that has valid auth credentials can login
AUTH_LDAP_USER_SEARCH = LDAPSearch("ou=users,dc=my-company,dc=com",
                                   ldap.SCOPE_SUBTREE, "(uid=%(user)s)")

AUTH_LDAP_GROUP_SEARCH = LDAPSearch("ou=groups,dc=my-company,dc=com",
                                    ldap.SCOPE_SUBTREE, "(objectClass=posixGroup)"
                                    )
AUTH_LDAP_GROUP_TYPE = PosixGroupType()


# Populate the Django user model from the LDAP directory.
AUTH_LDAP_USER_ATTR_MAP = {
    "first_name": "givenName",
    "last_name": "sn",
    "email": "mail"
}

# Map LDAP group membership into Django admin flags
AUTH_LDAP_USER_FLAGS_BY_GROUP = {
    # Anyone in this group is a superuser for the app
    "is_superuser": "cn=admin,ou=groups,dc=my-company,dc=com"
}
```

Also ensure 'ldap_auth' is listed in `PLUGINS_ENABLED` inside `shub/settings/config.py`.

Finally, you must build the Docker image with the build argument ENABLE_LDAP set to true:
```bash
docker build --build-arg ENABLE_LDAP=true -t quay.io/vanessa/sregistry . 
```

It's recommended to have the uwsgi logs open so any issue with ldap is shown clearly. You can do that with:

```bash
docker-compose logs -f uwsgi
```

For example, if you put in an incorrect credential, you would see the following in the logs:

```bash
uwsgi_1   | [pid: 56|app: 0|req: 4/4] 172.17.0.1 () {42 vars in 1025 bytes} [Thu Oct 26 07:18:10 2017] GET /ldap_auth/login/?next=http://127.0.0.1/login/ => generated 13475 bytes in 26 msecs (HTTP/1.1 200) 7 headers in 382 bytes (1 switches on core 0)
uwsgi_1   | search_s('ou=users,dc=my-company,dc=com', 2, '(uid=%(user)s)') returned 0 objects: 
uwsgi_1   | Authentication failed for adminuser: failed to map the username to a DN.
```

Once you have set these options, startup sregistry and you should be able to see the ldap option on the login page:

![ldap.png](ldap.png)

and login with the username/password pairs *testuser/testuser* and *testadmin/testadmin*. As a final note, if you choose this method to deploy an actual ldap server, you might consider adding the container to the docker-compose. If you've done this and need help, or want to contribute what you've learned, please submit a Pull Request to update these docs.
---
title: "pam-auth - Authentication with PAM"
pdf: true
toc: false
permalink: docs/plugins/pam
---

# PAM Authentication

The `pam_auth` plugin allows users to login to sregistry using the unix accounts on 
the host system.

To enable PAM authentication you must:
  * Add `pam_auth` to the `PLUGINS_ENABLED` list in `shub/settings/config.py`
  * Uncomment binds to /etc/shadow and /etc/passwd in `docker-compose.yml`
  * Build the docker image with the build argument ENABLE_PAM set to true
More detailed instructions are below.

## Permissions

The rules with respect to user collections still hold true - each user is given
push access given that they are added to a team, or `USER_COLLECTIONS` is true,
and each user will still each need to export their token to push.  You can read [more about roles here](https://singularityhub.github.io/sregistry/setup-roles), and [more about teams](https://singularityhub.github.io/sregistry/setup-teams) to manage groups of people.


## Getting Started

This is the detailed walkthough to set up the PAM AUthentication plugin. 

First, uncomment "pam_auth" at the bottom of `shub/settings/config.py` to 
enable the login option.

```bash
PLUGINS_ENABLED = [
#    'ldap_auth',
    'pam_auth',
#    'globus',
#    'saml_auth'
]
```

Since we need to get access to users from the host,
you need to edit the `docker-compose.yml` and uncomment binds to your host:

```bash
uwsgi:
  restart: always
  image: quay.io/vanessa/sregistry
  volumes:
    - .:/code
    - ./static:/var/www/static
    - ./images:/var/www/images
    # uncomment for PAM auth
    #- /etc/passwd:/etc/passwd 
    #- /etc/shadow:/etc/shadow
  links:
    - redis
    - db
```

If you do this, we lose the user added in the container for nginx! 
You also need to add the nginx user to your host:

```bash
$ sudo addgroup --system nginx
$ sudo adduser --disabled-login --system --home /var/cache/nginx --ingroup nginx nginx
```

Note that this solution [would require restarting the container](https://github.com/jupyterhub/jupyterhub/issues/535) for changes on the host to take effect (for example,
adding new users). If you find a better way to do this, please test and open an issue to add to this documentation.

Finally, you must build the docker image with the build argument ENABLE_PAM set to true:
```bash
$ docker build --build-arg ENABLE_PAM=true -t quay.io/vanessa/sregistry .
```
---
title: "saml-auth - Shibboleth Authentication"
pdf: true
toc: false
permalink: docs/plugins/saml
---

# SAML Authentication

The `saml_auth` plugin allows users to authentication with your [SAML provider](https://en.wikipedia.org/wiki/Security_Assertion_Markup_Language) of choice.

To enable SAML authentication you must:

  * Add `saml_auth` to the `PLUGINS_ENABLED` list in `shub/settings/config.py`
  * Add some configuration details to `shub/settings/config.py`
  * Configure the details of your SAML provider in in `shub/settings/secrets.py` per instructions provided [here](http://python-social-auth.readthedocs.io/en/latest/backends/saml.html).
  * Build the docker image with the build argument ENABLE_SAML set to true:
    ```bash
    $ docker build --build-arg ENABLE_SAML=true -t quay.io/vanessa/sregistry .
    ```


If you haven't yet created a secrets.py, a good start is to do the following:

```
cp shub/settings/dummy_secrets.py shub/settings/secrets.py
```


## Quick Start
This quick start is intended to demonstrate basic functionality of the SAML authentication. 


#### Edit Config.py

In the file `shub/settings/config.py` you should add the name of your institution (used to render the button)
along with the idp (the unique identifier for your SAML server request). That means uncommenting these lines.

```bash
# AUTH_SAML_IDP = "stanford"
# AUTH_SAML_INSTITUTION = "Stanford University"
```

so they appear like:


```bash
AUTH_SAML_IDP = "stanford"
AUTH_SAML_INSTITUTION = "Stanford University"
```

#### Setting up SAML Auth

In `secrets.py` you will need to define the variables specified [here](http://python-social-auth.readthedocs.io/en/latest/backends/saml.html), and that includes generating your certificate, which looks something like:

```bash
openssl req -new -x509 -days 3652 -nodes -out saml.crt -keyout saml.key
cat saml.key
mv saml.crt /etc/ssl/certs
mv saml.key /etc/ssl/private
```

and then generate the `metadata.xml` by going to `http://localhost/saml_auth/saml.xml`. Usually institutions have different portals for submitting metadata / getting information about SAML, for Stanford the information is via the [SAML SP Service Provider Database](https://spdb.stanford.edu/).

---
title: "Setup: Image Interaction"
pdf: true
toc: false
---

# Image Interaction

Before we work with accounts and other application setup, you need to know how to interact with the application, meaning Docker images. Here are the basic commands you should get comfortable with as an administrator of your registry. Note that these are really great for debugging too:

```
docker ps  # gets you the image ids
docker-compose logs worker  # show me the logs of the worker instance
docker-compose logs nginx  # logs for the main application
docker-compose logs -f nginx  # keep the logs open until I Control+C
```

and to shell in to an image, you can do the following:

```
NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
docker exec -it ${NAME} bash
```

Next, learn how to [create and manage](roles) different user roles.

---
title: Setup
pdf: true
toc: false
permalink: docs/setup/
---

# Setup

By the time you get here, you have added all required secrets and settings to your [settings](https://github.com/singularityhub/sregistry/tree/master/shub/settings) folder, and you've built and started the image. Next, you should navigate to [http://127.0.0.1](http://127.0.0.1) (localhost) to make sure the registry server is up and running. These sections will
detail setup of your registry, including user roles, teams, and interaction.

 1. [Interact](interact) with your images, an essential step to controlling and using it.
 2. [User Roles](roles) to give you control over who can create and access images.
 3. [Teams](teams) of users can be created to give group level access to image collections.
 4. [Registration](registration) of your registry so others can easily find it.

And don't forget to read about [accounts and credentials](../accounts/credentials) to govern the above.
---
title: Teams
pdf: true
toc: false
---

# Teams

To add a level of organization of users, sregistry has created loose groups of users called Teams. A registry admin can create a team, or if `USER_COLLECTIONS` is True, the an authenticated user can also create them. Creating a team means that the creator (admin or authenticated user) becomes the Owner of the team that can add and remove users. If an admin creates a team for a group of users, he or she must manage it or add a user to the list of owners to do the same. To create a team:

 1. Click on the "teams" tab in the navigation bar
 2. Choose a name, team name, and image.
 3. Decide if your team is "open" or "invite" only

There are two kinds of teams:

 - **invite** only means that an owner must send an invitation link.
 - **open** means that anyone can join the team that is authenticated in the registry.

![team-edit.png](../../assets/img/team-edit.png)

The default setting is "invite." Teams are important because when you add individuals as collaborators to your collections, they must come from one of your teams, and you do this on each Collection settings page:

![team-settings.png](../../assets/img/team-settings.png)

For example, if my lab has a set of users on sregistry and we intend to build images together, we would make a team for our lab, and then easily find one another to manage access to images.


## Create Accounts
To create your admin account, the first thing you need to do (in the web interface) is click Login on the top right. You should see the social account options that you configured in the [install](../install) step. You can now log in with your social account of choice. This will create a user account for yourself that you can use to manage the application, and when you are logged in you should see your username in the top right. It usually corresponds with the one from your social account.


### The Application Manager
At this point, you've started the application, and created a user with your social auth. Your username is in the top right, and it's usually the same as your social account. Keep this in mind because you will need it when you shell into the image to make yourself a superuser. Let's first shell inside:

```bash
NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
docker exec -it ${NAME} bash
```

you will find yourself in a directory called `/code`, which is where the main application lives. For administration you will be using the file [manage.py](https://github.com/singularityhub/sregistry/blob/master/manage.py) to interact with the registry. If you want to see all the different options, type `python manage.py` and it will show you.

Let's first make yourself a superuser and an admin, meaning that you are an administrator **and** have godlevel control of the registry. Just by way of being inside the Docker image you already have that. You will be able to set other people as admins. In summary:

 - `superuser`: you are an admin that can do anything, you have all permissions.
 - `admin`: you can push images, but not have significant impact on the registry application.

Of course anyone that shells into your Docker image could just explode everything - it's up to you to secure and manage the server itself! Let's say my username is `vsoch`. Here is the command I would run to make myself a superuser and admin. Note that we are inside the Docker image `sregistry_uwsgi_1`:

```bash
$ python manage.py add_superuser --username vsoch
DEBUG Username: vsoch
DEBUG vsoch is now a superuser.

$ python manage.py add_admin --username vsoch
DEBUG Username: vsoch
DEBUG vsoch can now manage and build.

# And from outside the Docker image
NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
docker exec $NAME python /code/manage.py add_superuser --username vsoch
docker exec $NAME python /code/manage.py add_admin --username vsoch
```

You can also choose to remove a superuser or admin at any time. This means he or she will no longer be able to build and access the token and secret to do so.


```bash
# Inside the image
$ python manage.py remove_superuser --username vsoch
$ python manage.py remove_admin --username vsoch

# Outside
NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
docker exec $NAME python /code/manage.py remove_superuser --username vsoch
docker exec $NAME python /code/manage.py remove_admin --username vsoch
```

Guess what! You don't need to continue to manage admins (and other superusers) from the command line after this step. When logged in to your superuser account, you will see an "Admin" link in your profile in the top right:

![admin-option.png](../../assets/img/admin-option.png)

This will take you to the administrative panel. Once there, you can click on "Users" at the bottom of the list, and select one or more checkboxes to assign other users to roles:

![admin-users.png](../../assets/img/admin-users.png)

[Register your registry](registration) next!
---
title: "Setup: User Roles"
pdf: true
toc: false
---

# Roles

Before we make accounts, let's talk about the different roles that can be associated with a registry. The core of Django supports some pre-defined roles, and we use those to the greatest extent that we can to keep things simple.

## superuser
You can think of as an application master. They have the highest level of permissions (akin to root) to shell into the application, add and remove users and roles, and do pretty much whatever they want. In that you are reading this and setting up the registry, you are going to be a superuser.

## admin
An admin corresponds with Django's "staff" role. An admin is designated by the superuser to have global ability to manage collections. This also gives permission to create teams, or groups of one or more users that can be added as contributors to a collection. An admin has a credential file to push images. An admin is a manager, but only of container collections, not the application.

## authenticated user
is a user that creates an account via the interface, but does not have a global ability to push images. Instead, the authenticated user can edit and manage collections that he or she contributes to.
  - If the variable `USER_COLLECTIONS` is set to True, the authenticated user can create and manage collections, and create teams. Each collection can have one or more owners and contributors, and both can push and pull images. Only owners can delete the collection or containers within it.
  - If the variable `USER_COLLECTIONS` is False, the authenticated user cannot create his or her own collections, but can still be added as a contributor to collections managed by admins. In this case, the admin is also in charge of creating teams.

Since you get to choose your authentication backend (e.g., LDAP, Twitter) you get to decide who can become an authenticated user. Here are a couple of scenarios:

 - You can keep container management tightly controlled by setting `USER_COLLECTIONS` to False, and then making a small set of individuals `admin`, meaning they manage public and/or private collections and teams for all users. In the case that a collection is private, the authenticated users must be added as contributors in the settings to view and pull images.
 - You can allow your users to manage their own collections teams by setting `USER_COLLECTIONS` to True. You can still have `admin` roles to be global managers, but put users in charge of managing their own images. The same rules apply with public and private collections - if a collection is private, the user would need to add collaborators to give pull ability.

## visitors
is an anonymous user of the registry. In the case of a private registry, this individual cannot pull images. In the case of a public registry, there is no huge benefit to being authenticated.

Based on the above and granted that you are setting up the server and reading this, you will be a **superuser** because you have permissions to control the Docker images and grant other users (and yourself) the ability to push with the role **admin**.

### Google Build + GitHub

If you have enabled the [Google Build+Github]({{ site.baseurl }}/plugin-google-build) plugin,
then your users will be able to log in with GitHub, and build collections that are
linked to GitHub repositories. In this case, permissions for the registry interaction
do not extend to GitHub. For example, if you build from a repository that you own,
adding a collaborator or another owner will not change anything on GitHub.

Speaking of collaborators, next, learn how users can be a part of [teams](teams)
---
title: Deployment backup
description: Deployment backup
tags: 
 - docker
---

# Backup

With changes to models and heavy development, there can be mistakes and errors
related to losing data. We can only do our best to back up the data locally,
back up the containers, and take snapshots of the server. This guide will provide 
detail to that.

## Snapshots

If you are using Google Cloud, Google Cloud makes it easy to generate a snapshot schedule,
and then edit a Disk to associate it. Full instructions are [here](https://cloud.google.com/compute/docs/disks/scheduled-snapshots), and basically it comes down to:

 - Creating a snapshot schedule under Compute -> Snapshots. I typically chose daily, with snapshots expiring after 14 days
 - Editing the Disk under Compute -> Disks to use it.


## Containers

Since this is run directly on the server using Docker, it must be run with
the local cron. You can use crontab -e to edit the cron file, and crontab -l
to list and verify the edits. Specifically, you need to add this line:

```cron
0 1 * * * docker commit sregistry_uwsgi_1 quay.io/vanessa/sregistry:snapshot
0 1 * * * docker commit sregistry_db_1 quay.io/vanessa/sregistry-postgres:snapshot
0 1 * * * docker commit sregistry_redis_1 quay.io/vanessa/sregistry-redis:snapshot
0 1 * * * docker commit sregistry_scheduler_1 quay.io/vanessa/sregistry-scheduler:snapshot
0 1 * * * docker commit sregistry_worker_1 quay.io/vanessa/sregistry-worker:snapshot
```

This will run a docker commit at 1:00am, daily, using the container name
"sregistry_uwsgi_1" and saving to "quay.io/vanessa/sregistry:snapshot". When the snapshot
is saved for the disk (at 3-4am per the previous instructions) then
this container should be included. A few notes:

 - make sure that your containers are named as specified above - if you start from a different folder or use a different version of compose, you might see differences from the commands above.
 - You don't need to backup any database (db) container if you are using a non-container database.
 - The worker, scheduler, and base container are the same, so you can remove the last two lines if desired, and then if you need to restore, just tag the base image with the other names.

The above will commit the containers to your host, but of course doesn't help if you lose the host (hence why an additional
backup strategy for your host is recommended, as suggested above with the Google Cloud backup strategy).

## Database

For some hosts, I've found that the scheduled cron jobs to backup the container will work interactively, but not
successfully with cron. Thus, it's a better idea to run the database backup scripts from the host
via cron:

```
0 3 * * * docker exec sregistry_uwsgi_1 /code/scripts/backup_db.sh
```

It's also good practice to test running this script manually from inside the container,
along with running it via docker exec (outside of the cron job) after that. Finally,
once you've checked those two, you should also check the timestamp on the backup files
one day later to ensure that it ran on its own.
---
title: "Setup: Registration"
pdf: true
toc: false
---

# Registration

We maintain a "registry of registries" ([https://singularityhub.github.io/containers](https://singularityhub.github.io/containers) (one registry to rule them all...) where you can easily have your registry's public images available for others to find. Adding your registry is easy - it comes down to automatically generating a file, adding it to the repo, and then doing a pull request (PR). Specifically:


## 1. Fork the repo
Fork the repo, and clone to your machine. That might look like this, given a username `vsoch`:

```bash
git clone https://www.github.com/vsoch/containers
cd containers
```

## 2. Generate your Metadata
Then use the manager to generate a markdown file for your registry. In this example, we have 
interactively shelled inside the `uwsgi` image for our sregistry like `docker run -it <container> bash` 
and are sitting in the `/code` folder.

```bash
# Inside the image
$ python manage.py register

Registry template written to taco-registry.md!

Your robot is at https://vsoch.github.io/robots/assets/img/robots/robot5413.png
1. Fork and clone https://www.github.com/singularityhub/containers
2. Add taco-registry.md to the registries folder
3. Download your robot (or add custom institution image) under assets/img/[custom/robots]
4. Submit a PR to validate your registry.
```

Specifically, this produces a markdown file in the present working directory (which is mapped to your host) that can be plopped into a folder. It is named based on your registry `uri`, and looks like this:

```bash
$ cat taco-registry.md 
---
layout: registry
title:  "Tacosaurus Computing Center"
base: http://127.0.0.1
date:   2017-08-30 17:45:44
author: vsochat
categories:
- registry
img: robots/robot5413.png
thumb: robots/robot5413.png # wget https://vsoch.github.io/robots/assets/img/robots/robot15570.png
tagged: taco
institution: Tacosaurus Computing Center
---

Tacosaurus Computing Center is a Singularity Registry to provide institution level Singularity containers.

```

At this point, you can send this file to `@vsoch` and she will be happy to add your 
registry to the... registry! If you want to customize your robot image, or submit
the file yourself via a pull request, continue reading!

## 3. Choose your image
For the image and thumbnail, we have a [database of robots](https://vsoch.github.io/robots) that we have randomly selected a robot for you. If you don't like your robot, feel free to browse and choose a different one. Importantly, you will need to add the robot to the github repo:

```bash
cd containers/assets/img/robots
wget https://vsoch.github.io/robots/assets/img/robots/robot15570.png
```

If you have some other custom image, add it to the "custom" folder. If it's not created yet, make it.

```bash
cd containers/assets/img
mkdir -p custom
cd custom
mv /path/to/institution/logo/taco-logo.png
```

Then for each of the `thumb` and `img` fields you would want to look like this:

```bash
img: custom/taco-logo.png
thumb: custom/taco-logo.png
```

## 4. Submit a PR
You can then add your files, and submit a PR to the main repo. We will have tests that ping your registry to ensure correct naming of files and registry address, along with a preview of the content that is added. If you want to prevew locally, you can run `jekyll serve`.


Great! Now that you have your accounts, you probably want to learn about how to build and push images! 
To push directly, you will first need to generate a [credential](../accounts/credentials). If you 
have enabled the [Google Build+Github](../plugins/google-build) plugin,
then you will be able to log in with GitHub, and connect GitHub repositories to build 
on commit. Either way, you should next read about the [client](../client).
---
title: "Installation: Web Server and Storage"
pdf: true
---

# Installation: Web Server and Storage

Before doing `docker-compose up -d` to start the containers, there are some specific things that need to be set up.

## Release Version

If you've downloaded a [release](https://github.com/singularityhub/sregistry/releases/),
you'll want to update the `docker-compose.yaml` instances of the quay.io/vanessa/sregistry
image to be tagged with the release version that matches your clone. E.g.:

```yaml
# Change instances of
  image: quay.io/vanessa/sregistry
# to
  image: quay.io/vanessa/sregistry:1.1.34
```

If you have cloned master, then the current master should coincide with 
latest and you don't need to do these updates.

## Under Maintenance Page

If it's ever the case that the Docker images need to be brought down for maintenance, a static fallback page should be available to notify the user. If you noticed in the [prepare_instance.sh](https://github.com/singularityhub/sregistry/blob/master/scripts/prepare_instance.sh) script, one of the things we installed is nginx (on the instance). This is because we need to use it to get proper certificates for our domain (for https). Before you do this, you might want to copy the index that we've provided to replace the default (some lame page that says welcome to Nginx!) to one that you can show when the server is undergoing maintenance.

```bash
cp $INSTALL_ROOT/sregistry/scripts/nginx-index.html /var/www/html/index.html
rm /var/www/html/index.nginx-debian.html
```

If you want your page to use the same SSL certificates, a nginx-default.conf is also
provided that will point to the same certificates on the server (generation discussed later).
Please execute this command only after the certificates have been generated, or you won't be able to generate them:

```bash
cp $INSTALL_ROOT/sregistry/scripts/nginx-default.conf /etc/nginx/conf.d/default.conf
```

If you don't care about user experience during updates and server downtime, you can just ignore this.

## Custom Domain

In the [config settings file](https://github.com/singularityhub/sregistry/blob/master/shub/settings/config.py#L30)
you'll find a section for domain names, and other metadata about your registry. You will need to update
this to be a custom hostname that you use, and custom names and unique resource identifiers for your
registry. For example, if you have a Google Domain and are using Google Cloud, you should be able to set it up using [Cloud DNS](https://console.cloud.google.com/net-services/dns/api/enable?nextPath=%2Fzones&project=singularity-static-registry&authuser=1). Usually this means
creating a zone for your instance, adding a Google Domain, and copying the DNS records for
the domain into Google Domains. Sometimes it can take a few days for changes to propagate.
You are strongly encouraged to register both `your.domain.com`, as well as `www.your.domain.com` and have them point to the same IP address.
We will discuss setting up https in a later section.

## Storage

For Singularity Registry versions prior to 1.1.24, containers that were uploaded to your registry
were stored on the filesystem, specifically at `/var/www/images` that was bound to the host
at `images`. We did this by way of using a custom nginx image with the nginx upload module
enabled (see [this post](https://vsoch.github.io/2018/django-nginx-upload/) for an example).

There is also the other option to use [custom builders]({{ site.url }}/install-builders) 
that can be used to push a recipe to Singularity Registry Server, and then trigger a 
cloud build that will be saved in some matching cloud storage.

### Default

However for versions 1.1.24 and later, to better handle the Singularity `library://`
client that uses Amazon S3, we added a [Minio Storage](https://docs.min.io/docs/python-client-api-reference.html#presigned_get_object) 
backend, or another container (minio) that is deployed alongside Singularity Registry server.
If you look in the [docker-compose.yml](https://github.com/singularityhub/sregistry/blob/master/docker-compose.yml) that looks something like this:

```
minio:
  image: minio/minio
  volumes:
    - ./minio-images:/images
  env_file:
   - ./.minio-env
  ports:
   - "9000:9000"  
  command: ["server", "images"]
```


At the time of development we are using this version of minio:

```
/ # minio --version
minio version RELEASE.2020-04-02T21-34-49Z
```

which you can set in the docker-compose.yml file to pin it. Notice that we bind the
folder "minio-images" in the present working directory to /images in the container,
which is where we are telling minio to write images to the filesystem. This means
that if your container goes away, the image files will still be present on the host.
For example, after pushing two images, I can see them organized by bucket, collection,
then container name with hash.

```
$ tree minio-images/
minio-images/
└── sregistry
    └── test
        ├── big:sha256.92278b7c046c0acf0952b3e1663b8abb819c260e8a96705bad90833d87ca0874
        └── container:sha256.c8dea5eb23dec3c53130fabea451703609478d7d75a8aec0fec238770fd5be6e
```

### Configuration

For secrets (the access and secret key that are used to create the container) 
we are reading in environment variables for the server in `.minio-env`
that looks like this:

```bash
# Credentials (change the first two)
MINIO_ACCESS_KEY=newminio
MINIO_SECRET_KEY=newminio123
MINIO_ACCESS_KEY_OLD=minio
MINIO_SECRET_KEY_OLD=minio123

# Turn on/off the Minio browser at 127.0.0.1:9000?
MINIO_BROWSER=on

# Don't clean up images in Minio that are no longer referenced by sregistry
DISABLE_MINIO_CLEANUP = False
```

If you [read about how minio is configured](https://github.com/minio/minio/tree/master/docs/config)
you will notice that since the container was built previously with the default credentials,
to change them we are required to define the same variables with `*_OLD` suffix. This
means that you should be able to keep the last two lines (minio and minio123, respectively)
and change the first two to your new access key and secret. And if it's not clear:

> **You should obviously change the access key and secret in the minio-env file!**

The `.minio-env` file is also bound to the uwsgi container, so that the generation of the minio
storage can be authenticated by the uwsgi container, which is the interface between
the Singularity client and minio. For variables that aren't secrets, you can look
in `shub/settings/config.py` and look for the "Storage" section with various
minio variables:

```python
MINIO_SERVER = "minio:9000"  # Internal to sregistry
MINIO_EXTERNAL_SERVER = (
    "127.0.0.1:9000"  # minio server for Singularity to interact with
)
MINIO_BUCKET = "sregistry"
MINIO_SSL = False  # use SSL for minio
MINIO_SIGNED_URL_EXPIRE_MINUTES = 5
MINIO_REGION = "us-east-1"
MINIO_MULTIPART_UPLOAD = True
```

Since the container networking space is different from what the external
Singularity client interacts with, we define them both here. If you deploy
a minio server external to the docker-compose.yml, you can update both of
these URLs to be the url to access it. The number of minutes for the signed
url to expire applies to single PUT (upload), GET (download), and upload Part (PUT) requests.
Finally, the logs that you see with `docker-compose logs minio` are fairly limited,
it's recommended to install the client [mc](https://docs.minio.io/docs/minio-client-quickstart-guide)
to better inspect:

```bash
wget https://dl.min.io/client/mc/release/linux-amd64/mc
chmod +x mc
./mc --help
```

You'll still need to add the host manually:

```bash
./mc config host add myminio http://127.0.0.1:9000 $MINIO_ACCESS_KEY $MINIO_SECRET_KEY
Added `myminio` successfully.
```

You can then list the hosts that are known as follows (and make sure yours appears)

```bash
/ # ./mc config host ls
gcs  
  URL       : https://storage.googleapis.com
  AccessKey : YOUR-ACCESS-KEY-HERE
  SecretKey : YOUR-SECRET-KEY-HERE
  API       : S3v2
  Lookup    : dns

local
  URL       : http://localhost:9000
  AccessKey : 
  SecretKey : 
  API       : 
  Lookup    : auto

myminio
  URL       : https://127.0.0.1:9000
  AccessKey : YOUR-ACCESS-KEY-HERE
  SecretKey : YOUR-SECRET-KEY-HERE
  API       : S3v4
  Lookup    : auto


play 
  URL       : https://play.min.io
  AccessKey : YOUR-ACCESS-KEY-HERE
  SecretKey : YOUR-SECRET-KEY-HERE
  API       : S3v4
  Lookup    : auto

s3   
  URL       : https://s3.amazonaws.com
  AccessKey : YOUR-ACCESS-KEY-HERE
  SecretKey : YOUR-SECRET-KEY-HERE
  API       : S3v4
  Lookup    : dns
```

And then getting logs like this:

```bash
./mc admin trace -v myminio
```

When the trace is running, you should be able to do some operation and see output.
(It will just hang there open while it's waiting for you to try a push).
This is really helpful for debugging!

```bash
127.0.0.1 [REQUEST s3.PutObjectPart] 22:59:38.025
127.0.0.1 PUT /sregistry/test/big:sha256.92278b7c046c0acf0952b3e1663b8abb819c260e8a96705bad90833d87ca0874?uploadId=a1071852-3407-4c2b-9444-6790cfafae51&partNumber=1&Expires=1586041182&Signature=2dN0tY%2F0esPKVDDD%2B%2F1584I0qqQ%3D&AWSAccessKeyId=minio
127.0.0.1 Host: 127.0.0.1:9000
127.0.0.1 Content-Length: 928
127.0.0.1 User-Agent: Go-http-client/1.1
127.0.0.1 X-Amz-Content-Sha256: 2fc597b42f249400d24a12904033454931eb3624e8c048fe825c360d9c1e61bf
127.0.0.1 Accept-Encoding: gzip
127.0.0.1 <BODY>
127.0.0.1 [RESPONSE] [22:59:38.025] [ Duration 670µs  ↑ 68 B  ↓ 806 B ]
127.0.0.1 403 Forbidden
127.0.0.1 X-Xss-Protection: 1; mode=block
127.0.0.1 Accept-Ranges: bytes
127.0.0.1 Content-Length: 549
127.0.0.1 Content-Security-Policy: block-all-mixed-content
127.0.0.1 Content-Type: application/xml
127.0.0.1 Server: MinIO/RELEASE.2020-04-02T21-34-49Z
127.0.0.1 Vary: Origin
127.0.0.1 X-Amz-Request-Id: 1602C00C5749AF1F
127.0.0.1 <?xml version="1.0" encoding="UTF-8"?>
<Error><Code>SignatureDoesNotMatch</Code><Message>The request signature we calculated does not match the signature you provided. Check your key and signing method.</Message><Key>test/big:sha256.92278b7c046c0acf0952b3e1663b8abb819c260e8a96705bad90833d87ca0874</Key><BucketName>sregistry</BucketName><Resource>/sregistry/test/big:sha256.92278b7c046c0acf0952b3e1663b8abb819c260e8a96705bad90833d87ca0874</Resource><RequestId>1602C00C5749AF1F</RequestId><HostId>e9ba6dec-55a9-4add-a56b-dd42a817edf2</HostId></Error>
127.0.0.1 
```
The above shows an error - the signature is somehow wrong. It came down to not specifying the signature type in a config
when I instantiated the client, and also needing to customize the `presign_v4` function to allowing
sending along the sha256sum from the scs-library-client (it was using an unsigned hash).
For other configuration settings (cache, region, notifications) you should consult
the [configuration documentation](https://github.com/minio/minio/tree/master/docs/config).
SSL instructions for minio are included in the next section.

## SSL

### Server Certificates

Getting https certificates is really annoying, and getting `dhparams.pem` takes forever. But after the domain is obtained, it's important to do. Again, remember that we are working on the host, and we have an nginx server running. You should follow the instructions (and I do this manually) in [generate_cert.sh](https://github.com/singularityhub/sregistry/blob/master/scripts/generate_cert.sh). 

 - starting nginx
 - installing certbot
 - generating certificates
 - linking them to where the docker-compose expects them
 - add a reminder or some other method to renew within 89 days

With certbot, you should be able to run `certbot renew` when the time to renew comes up. There is also an [older
version](https://github.com/singularityhub/sregistry/blob/master/scripts/generate_cert_tiny-acme.sh) that uses tiny-acme instead of certbot. For this second option, it basically comes down to:

 - starting nginx
 - installing tiny acme
 - generating certificates
 - using tinyacme to get them certified
 - moving them to where they need to be.
 - add a reminder or some other method to renew within 89 days

Once you have done this (and you are ready for https), you should use the `docker-compose.yml` and the `nginx.conf` provided in the folder [https](https://github.com/singularityhub/sregistry/blob/master/https/). So do something like this:

```bash
mkdir http
mv nginx.conf http
mv docker-compose.yml http

cp https/docker-compose.yml .
cp https/nginx.conf.https nginx.conf
```

If you run into strange errors regarding any kind of authentication / server / nginx when you start the images, likely it has to do with not having moved these files, or a setting about https in the [settings](https://github.com/singularityhub/sregistry/tree/master/shub/settings). If you have trouble, please post an issue on the [issues board](https://www.github.com/singularityhub/sregistry/issues) and I'd be glad to help.

### Minio

Minio has detailed instructions for setting up https [here](https://github.com/minio/minio/tree/master/docs/tls),
and you can use existing or self signed certificates. The docker-compose.yml in the https folder takes the
strategy of binding the same SSL certs to `${HOME}/.minio/certs` in the container.
Note that inside the certs directory, the private key must by named private.key and the public key must be named public.crt.

## Build the Image (Optional)
If you want to try it, you can build the image. Note that this step isn't necessary as the image is provided on [Quay.io]({{ site.registry }}). This step is optional. However, if you are developing you likely want to build the image locally. You can do:


```bash
cd sregistry
docker build -t quay.io/vanessa/sregistry .
```

## Nginx

This section is mostly for your FYI. The nginx container that we used to rely on for uploads is a custom compiled
nginx that included the [nginx uploads module](https://www.nginx.com/resources/wiki/modules/upload/).
This allowed us to define a server block that would accept multipart form data directly, and 
allowed uploads directly to the server without needing to stress the uwsgi application. 
The previous image we used is still provided on Docker Hub at
[quay.io/vanessa/sregistry_nginx](https://hub.docker.com/r/quay.io/vanessa/sregistry_nginx)
While we don't use this for standard uploads, we still use it for the web interface upload,
and thus have not removed it. To build this image on your own change this section:


```bash
nginx:
  restart: always
  image: quay.io/vanessa/sregistry_nginx
  ports:
    - "80:80"
  volumes:
    - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
    - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
  volumes_from:
    - uwsgi
  links:
    - uwsgi
    - db
``` 
to this, meaning that we will build from the nginx folder:

```bash
nginx:
  restart: always
  build: nginx
  ports:
    - "80:80"
  volumes:
    - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
    - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
  volumes_from:
    - uwsgi
  links:
    - uwsgi
    - db
```

Next, why don't you [start Docker containers](containers) to get your own registry going.
---
title: Custom Builders and Storage
pdf: true
toc: false
---

Currently, we support custom installation of the following builder and storage pairs. Each of these is provided
as a plugin, so you can enable them in the same way. Instructions are included with the links below.

 - [Google Cloud Build + Storage]({{ site.baseurl }}/plugin-google-build)

Don't forget to go back to the [install docs](https://singularityhub.github.io/sregistry/install-server#storage) where you left off.

<div>
    <a href="/sregistry/install"><button class="previous-button btn btn-primary"><i class="fa fa-chevron-left"></i> </button></a>
</div><br>
---
title: "Installation: Host"
pdf: true
toc: false
---

# Host Installation

## Dependencies

If starting from scratch, you need the following dependencies on the host:

 - [Docker](https://docs.docker.com/install/): a container engine
 - [docker-compose](https://docs.docker.com/compose/install/): an orchestration tool for Docker images.
 - python: docker compose requires some additional python libraries, `ipaddress` and `oauth2client`

{% include alert.html type="info" title="Important!" content="If you are just installing Docker <strong>you will need to log in and out after adding your user to the Docker group</strong>. " %}

## Steps

For a record of the installation procedure that I used for a Google Cloud host, I've provided the [basic commands](https://github.com/singularityhub/sregistry/blob/master/scripts/prepare_instance.sh). This script was run manually for an instance. This was done on a fairly large fresh ubuntu:16.04 instance on Google Cloud. This setup is done only once, and requires logging in and out of the instance after installing Docker, but before bringing the instance up with `docker-compose`. A few important points:

- The `$INSTALL_BASE` is set by default to `/opt`. It is recommended to choose somewhere like `/opt` or `/share` that is accessible by all those who will maintain the installation. If you choose your home directory, you can expect that only you would see it. If it's for personal use, `$HOME` is fine.
- Anaconda3 is installed for python libraries needed for `docker-compose`. You can use whatever python you with (system installed or virtual environment)
- Make sure to add all users to the docker group that need to maintain the application, and log in and out before use.

For the rest of the install procedure, you should (if you haven't already) clone the repository:

```bash
git clone https://github.com/singularityhub/sregistry
cd sregistry
```

For the files linked below, you should find the correspoinding file in the Github repository that you cloned. If you are setting this up for the first time, it's recommended to try locally and then move onto your production resource.

Next, why don't you [configure settings](settings) to customize your installation.
---
title: "Installation: Settings"
pdf: true
toc: false
---

# Settings

See that folder called [settings](https://github.com/singularityhub/sregistry/blob/master/shub/settings)? inside are a bunch of different starting settings for the application. We will change them in these files before we start the application. There are actually only two files you need to poke into, generating a `settings/secrets.py` from our template [settings/dummy_secrets.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/dummy_secrets.py) for application secrets, and [settings/config.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/config.py) to configure your database and registry information.

## Secrets

There should be a file called `secrets.py` in the shub settings folder (it won't exist in the repo, you have to make it), in which you will store the application secret and other social login credentials.

An template to work from is provided in the settings folder called `dummy_secrets.py`. You can copy this template:

```bash
cp shub/settings/dummy_secrets.py shub/settings/secrets.py
```

Or, if you prefer a clean secrets file, create a blank one as below:

```bash
touch shub/settings/secrets.py
```

and inside you want to add a `SECRET_KEY`. You can use the [secret key generator](http://www.miniwebtool.com/django-secret-key-generator/) to make a new secret key, and call it `SECRET_KEY` in your `secrets.py` file, like this:

```      
SECRET_KEY = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
```

### Authentication Secrets

One thing I (@vsoch) can't do for you in advance is produce application keys and secrets to give your Registry for each social provider that you want to allow users (and yourself) to login with. We are going to use a framework called [python social auth](https://python-social-auth-docs.readthedocs.io/en/latest/configuration/django.html) to achieve this, and in fact you can add a [number of providers](http://python-social-auth-docs.readthedocs.io/en/latest/backends/index.html) (I have set up a lot of them, including SAML, so please <a href="https://www.github.com/singularityhub/sregistry/isses" target="_blank">submit an issue</a> if you want one added to the base proper.). Singularity Registry uses OAuth2 with a token--> refresh flow because it gives the user power to revoke permission at any point, and is a more modern strategy than storing a database of usernames and passwords. You can enable or disable as many of these that you want, and this is done in the [settings/config.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/config.py):


```python
# Which social auths do you want to use?
ENABLE_GOOGLE_AUTH=False
ENABLE_TWITTER_AUTH=False
ENABLE_GITHUB_AUTH=True
ENABLE_GITLAB_AUTH=False
ENABLE_BITBUCKET_AUTH=False
ENABLE_GITHUB_ENTERPRISE_AUTH=False
```

and you will need at least one to log in. I've found that GitHub works the fastest and easiest, and then Google.
Twitter now requires an actual server name and won't work with localhost, but if you are deploying on a server with a proper domain go ahead and use it. All avenues are extremely specific with regard to callback urls, so you should be very careful in setting them up. If you want automated builds from a repository
integration with Google Cloud Build, then you must use GitHub.

## Plugins

Other authentication methods, such as LDAP, are implemented as [plugins](https://singularityhub.github.io/sregistry/docs/plugins/) to sregistry. See the [plugins documentation](https://singularityhub.github.io/sregistry/docs/plugins/) for details on how to configure these. You should also now look here to see which plugins you will
want to set up (and then build into your container).

For authentication plugins, we will walk through the setup of each in detail here. 
For other plugins, you should look at the [plugins](https://singularityhub.github.io/sregistry/docs/plugins/) documentation now before proceeding. For all of the below, you should put the content in your `secrets.py` under settings. Note that if you are deploying locally, you will need to put localhost (127.0.0.1) as your domain, and Github is now the only one that worked reliably without an actual domain for me.

### Google OAuth2

You first need to [follow the instructions](https://developers.google.com/identity/protocols/OpenIDConnect) and setup an OAuth2 API credential. The redirect URL should be every variation of having http/https, and www. and not. (Eg, change around http-->https and with and without www.) of `https://www.sregistry.org/complete/google-oauth2/`. Google has good enough debugging that if you get this wrong, it will give you an error message with what is going wrong. You should store the credential in `secrets.py`, along with the complete path to the file for your application:

```python
GOOGLE_CLIENT_FILE='/code/.grilledcheese.json'

# http://psa.matiasaguirre.net/docs/backends/google.html?highlight=google
SOCIAL_AUTH_GOOGLE_OAUTH2_KEY = 'xxxxxxxxxxxxxxxxxxxxxx.apps.googleusercontent.com'
SOCIAL_AUTH_GOOGLE_OAUTH2_SECRET = 'xxxxxxxxxxxxxxxxx'
# The scope is not needed, unless you want to develop something new.
#SOCIAL_AUTH_GOOGLE_OAUTH2_SCOPE = ['https://www.googleapis.com/auth/drive']
SOCIAL_AUTH_GOOGLE_OAUTH2_AUTH_EXTRA_ARGUMENTS = {
    'access_type': 'offline',
    'approval_prompt': 'auto'
}
```

Google is great in letting you specify multiple acceptable callback urls, so you should set every version of `http://127.0.0.1/complete/google-oauth2` (I did with and without http/https, along with the ending and without the ending slash, just in case). Note that `1.` extra arguments have been added to ensure that users can refresh tokens, and `2.` in testing I was using `http` and not `https`, and I eventually added `https` (and so the url was adjusted accordingly). Next, we need to follow instructions for [web applications](https://developers.google.com/identity/protocols/OAuth2WebServer). 


### Setting up Github OAuth

For users to connect to Github, you need to [register a new application](https://github.com/settings/applications/new), and add the key and secret to your `secrets.py` file like this: 

```python
# http://psa.matiasaguirre.net/docs/backends/github.html?highlight=github
SOCIAL_AUTH_GITHUB_KEY = ''
SOCIAL_AUTH_GITHUB_SECRET = ''

# If you want to use the google_build plugin, you will need to include the following:
SOCIAL_AUTH_GITHUB_SCOPE = ["admin:repo_hook",
                            "repo:status",
                            "user:email",
                            "read:org",
                            "admin:org_hook",
                            "deployment_status"]
```

The callback url should be in the format `http://127.0.0.1/complete/github`, and replace the localhost address with your domain. See the [Github Developers](https://github.com/settings/developers) pages to browse more information on the Github APIs.


### Setting up Github Enterprise OAuth

The GitHub Exterprise [docs are here](https://python-social-auth.readthedocs.io/en/latest/backends/github_enterprise.html).  You will want to register a new application on your instance of GitHub Enterprise in Developer Settings, set the callback URL to "http://example.com/complete/github-enterprise/" replacing example.com with your domain, and then the following environment variables should be defined in your secrets.

```python
# The URL for your GitHub Enterprise appliance:
SOCIAL_AUTH_GITHUB_ENTERPRISE_URL = "https://git.example.com/"

# Set the API URL for your GitHub Enterprise appliance:
SOCIAL_AUTH_GITHUB_ENTERPRISE_API_URL = "https://git.example.com/api/v3/"

# Fill the Client ID and Client Secret values from GitHub in the settings:
SOCIAL_AUTH_GITHUB_ENTERPRISE_KEY = ""
SOCIAL_AUTH_GITHUB_ENTERPRISE_SECRET = ""
```

### Gitlab OAuth2

Instructions are provided [here](https://github.com/python-social-auth/social-docs/blob/master/docs/backends/gitlab.rst). Basically:

1. You need to [register an application](https://gitlab.com/profile/applications), be sure to add the `read_user` scope. If you need `api`, add it to (you shouldn't).
2. Set the callback URL to `http://registry.domain/complete/gitlab/`. The URL **must** match the value sent. If you are having issues, try adjusting the trailing slash or http/https/. 
3. In your `secrets.py` file under settings, add:

```
SOCIAL_AUTH_GITLAB_SCOPE = ['api', 'read_user']
SOCIAL_AUTH_GITLAB_KEY = ''
SOCIAL_AUTH_GITLAB_SECRET = ''
```
Where the key and secret are replaced by the ones given to you. If you have a private Gitlab, you need to add it's url too:

```
SOCIAL_AUTH_GITLAB_API_URL = 'https://example.com'
```


### Bitbucket OAuth2

We will be using the [bitbucket](https://python-social-auth.readthedocs.io/en/latest/backends/bitbucket.html) backend for Python Social Auth.

First, register a new OAuth Consumer by following the instructions in the [Bitbucket documentation](https://confluence.atlassian.com/bitbucket/oauth-on-bitbucket-cloud-238027431.html). Overall, this means registering a new consumer, and making sure to add the "account" scope to it. You can find the button to add a consumer in your BitBucket profile (click your profile image from the bottom left of [the dashboard](https://bitbucket.org/dashboard/overview).

![../../assets/img/bitbucket-consumer.png](../../assets/img/bitbucket-consumer.png)

After clicking the button, fill in the following values:

 - Name: choose a name that will be easy to link and remember like "Singularity Registry Server"
 - Callback URL: should be `http://[your-domain]/complete/bitbucket` For localhost, this is usually `http://127.0.0.1/complete/bitbucket`
 - Keep the button "This is a private consumer" checked.
 - Under Permissions (the scope) click on Account (email, read, write).


Then, when you click to add the consumer, it will take you back to the original pacge. To get the key and secret, you should click on the name of the consumer. Then add the following variables to your `secrets.py` file under settings:

```python
SOCIAL_AUTH_BITBUCKET_OAUTH2_KEY = '<your-consumer-key>'
SOCIAL_AUTH_BITBUCKET_OAUTH2_SECRET = '<your-consumer-secret>'
```

 3. Optionally, if you want to limit access to only users with verified e-mail addresses, add the following:

```python
SOCIAL_AUTH_BITBUCKET_OAUTH2_VERIFIED_EMAILS_ONLY = True
```

Finally, don't forget to enable the bitbucket login in settings/config.py:

```python
ENABLE_BITBUCKET_AUTH=True
```

### Setting up Twitter OAuth2
You can go to the [Twitter Apps](https://apps.twitter.com/) dashboard, register an app, and add secrets, etc. to your `secrets.py`:

```bash
SOCIAL_AUTH_TWITTER_KEY = ''
SOCIAL_AUTH_TWITTER_SECRET = ''
```

Note that Twitter now does not accept localhost urls. Thus, 
the callback url here should be `http://[your-domain]/complete/twitter`.



## Config

In the [config.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/config.py) you need to define the following:

### Google Analytics

If you want to add a Google analytics code, you can do this in the settings/config.py:

```python
GOOGLE_ANALYTICS = "UA-XXXXXXXXX"
```

The default is set to None, and doesn't add analytics to the registry.


### Domain Name
A Singularity Registry Server should have a domain. It's not required, but it makes it much easier for yourself and users to remember the address. The first thing you should do is change the `DOMAIN_NAME_*` variables in your settings [settings/config.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/config.py#L30).

For local testing, you will want to change `DOMAIN_NAME` and `DOMAIN_NAME_HTTP` to be localhost. Also note that I've set the regular domain name (which should be https) to just http because I don't have https locally:

```python
DOMAIN_NAME = "http://127.0.0.1"
DOMAIN_NAME_HTTP = "http://127.0.0.1"
#DOMAIN_NAME = "https://singularity-hub.org"
#DOMAIN_NAME_HTTP = "http://singularity-hub.org"
```

It's up to the deployer to set one up a domain or subdomain for the server. Typically this means going into the hosting account to add the A and CNAME records, and then update the DNS servers. Since every host is a little different, I'll leave this up to you, but [here is how I did it on Google Cloud](https://cloud.google.com/dns/quickstart).


### Registry Contact
You need to define a registry uri, and different contact information:

```python
HELP_CONTACT_EMAIL = 'vsochat@stanford.edu'
HELP_INSTITUTION_SITE = 'https://srcc.stanford.edu'
REGISTRY_NAME = "Tacosaurus Computing Center"
REGISTRY_URI = "taco"
```

The `HELP_CONTACT_EMAIL` should be an email address that you want your users (and/or visitors to your registry site, if public) to find if they need help. The `HELP_INSTITUTION_SITE` is any online documentation that you want to be found in that same spot. Finally, `REGISTRY_NAME` is the long (human readable with spaces) name for your registry, and `REGISTRY_URI` is a string, all lowercase, 12 or fewer characters to describe your registry.


### User Collections
By default, any authenticated user in your Registry can create collections, and decide to make them public or private. If you would prefer to revoke this permission (meaning that only administrators can create and manage collections) then you would want to set this variable to `False`.

```python
# Allow users to create public collections
USER_COLLECTIONS=True
```

Setting `USER_COLLECTIONS` to False also means that users cannot create [Teams](../setup/teams), which are organized groups of users that then can be added as contributors to individual collections. With this setting as True, any authenticated user, staff, or administrator can create and manage new collections and teams, and this is done by issuing a token.

Finally, you can also allow users to create collections, but limit the number created.

```python
# Limit users to N collections (None is unlimited)
USER_COLLECTION_LIMIT = None
```

The default is None, meaning that users can create unlimited collections, given that `USER_COLLECTIONS`
is True. If you set this to a non-zero positive integer, user collections will be limited to
this number. If a user is staff or an administrator, they are not subject to this limit.


### Registry Private
By default Singularity Registry will provide public images, with an option to set them to private. If you are working with sensitive data and/or images, you might want all images to be private, with no option to make public. You can control that with the variable `PRIVATE_ONLY`.

```python
PRIVATE_ONLY=True
```

The above would eliminate public status and make private the default. Alternatively, if you want to allow for public images but make the default private (and collection owners can make collections of their choice public) set `DEFAULT_PRIVATE` to True.

```python
DEFAULT_PRIVATE=True
```

`PRIVATE_ONLY` takes preference to `DEFAULT_PRIVATE`. In other words, if you set `PRIVATE_ONLY` to True, the default has to be private, the change to `DEFAULT_PRIVATE` is meaningless, and a user cannot change a collection to be public.


### Collections Page Display

On the main server's `<domain>/collections` page, users will be shown
some limited set of collections, plus those that are private that they own.
Since this could slow down the page if the number is huge, you are given
control of this number:

```python
# The number of collections to show on the /<domain>/collections page
COLLECTIONS_VIEW_PAGE_COUNT=250
```

For larger registries, it's recommended to disable this view all together, and
encourage users to find containers via "search." If you think this should
be a default, please open an issue to discuss.


### View Rate Limits

While it's unlikely someone would be malicious to request a view, we can't
disregard it completely. For all views, we use django-ratelimit to limit
views to a certain number per day based on the ip address. For most views,
you can define the variables:

```python
VIEW_RATE_LIMIT="50/1d"  # The rate limit for each view, django-ratelimit, "50 per day per ipaddress)
VIEW_RATE_LIMIT_BLOCK=True # Given that someone goes over, are they blocked for the period?
```

In the example above, we limit each ip address to 50/day. We block any addresses
that go over, until the next period begins.

### Container GET Limits

Too many get requests for any particular container, whether stored locally or in
Google Storage (Google Cloud Build + GitHub plugin) could lead to a DoS for the server.
Thus, we have a limit on the total number of weekly GET requests per container:

```python
# The maximum number of downloads allowed per container, per week
CONTAINER_WEEKLY_GET_LIMIT=100
```

The `Container` object has a get_limit and get_count that are adjusted when a
user downloads a container. A weekly cron job will reset the get_count.

### Collection GET Limits

It could be the case that a user continually rebuilds containers to get around the
single container get limit, in which case we've also added a collection
weekly get limit.

```python
# The maximum number of downloads allowed per collection, per week
COLLECTION_WEEKLY_GET_LIMIT=100
```

The `Collection` object also has a get_limit and get_count that are adjusted when a
user downloads a container, reset by the same cron task.


### Disable Building

Disable all building, including pushing of containers and recipes. By
default, for a working registry, this should be False.

```python
DISABLE_BUILDING=True
```

For other disable and limit arguments (for GitHub, creating, or receiving builds) see the
[Google Build Plugin](/sregistry/docs/plugins/google-build).

### Database

By default, the database itself will be deployed as a postgres image called `db`. You probably don't want this for production (for example, I use a second instance with postgres and a third with a hot backup, but it's an ok solution for a small cluster or single user. Either way, we recommend backing it up every so often.

When your database is set up, you can define it in your `secrets.py` and it will override the Docker image one in the `settings/main.py file`. It should look something like this

```python
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'dbname',
        'USER': 'dbusername',
        'PASSWORD':'dbpassword',
        'HOST': 'localhost',
        'PORT': '5432',
    }
}
```

### Logging
By default, Singularity Registry keeps track of all requests to pull containers, and you have control over the level of detail that is kept. If you want to save complete metadata (meaning the full json response for each call) then you should set `LOGGING_SAVE_RESPONSES` to True. If you expect heavy use and want to save the minimal (keep track of which collections are pulled how many times) the reccomendation is to set this to False. 

```python
LOGGING_SAVE_RESPONSES=True
```

## API

Take a look in [settings/api.py](https://github.com/singularityhub/sregistry/blob/master/shub/settings/api.py)
to configure your restful API. You can uncomment the first block to require authentication:

```python
  'DEFAULT_PERMISSION_CLASSES': (
      'rest_framework.permissions.IsAuthenticated',
   ),
```

And also choose throttle rates for users and anonymous API requests.

```python
    # You can also customize the throttle rates, for anon and users
    'DEFAULT_THROTTLE_CLASSES': (
        'rest_framework.throttling.AnonRateThrottle',
    ),
    # https://www.django-rest-framework.org/api-guide/throttling/
    'DEFAULT_THROTTLE_RATES': {
        'anon': '100/day',
        'user': '1000/day',
    },
```

These are important metrics to ensure that your server isn't subject to a DoS attack.


Great job! Let's now [configure your web server and storage](server).
---
title: Installation
pdf: true
toc: false
permalink: docs/install/
---

# Installation

Installation will consist of the following steps:

 1. [Install host Dependencies](host) including Docker and compose to run the application.
 2. [Configure Settings](settings) to customize your installation.
 3. [Web Server and Storage](server) configuration based on your host needs.
 4. [Containers](containers) are ready to be started!

When you are done with this guide, you will have a running server! You can then move on to the [setup](../setup) guide to better understand how to interact with your fresh Singularity Registry.
---
title: "Installation: Containers"
pdf: true
toc: true
---

# Installation: Start Containers

Whether you build or not, the compose command will bring up the application (and
download containers provided on Quay.io, previously on Docker Hub, if they
aren't in your cache).

## What containers are provided?

Singularity Registry Server uses the following images, all provided on Quay.io
(or you can build the registry-specific ones locally):

 - [quay.io/vanessa/sregistry]({{ site.registry }}): is the core application
   image, generated from the Dockerfile in the base of the repository.
 - [quay.io/vanessa/sregistry_nginx]({{ site.registry }}_nginx/): Is the NGINX
   container installed with the NGINX upload module, intended for use with
   speedy uploads. It is generated from the subfolder "nginx" in the repository.

To use these images provided, you can bring up the containers like so:

## Start Containers

```bash
$ docker-compose up -d
```

The `-d` means detached, and that you won't see any output (or errors) to the
console. You can easily restart and stop containers, either specifying the
container name(s) or leaving blank to apply to all containers. Note that these
commands must be run in the folder with the `docker-compose.yml` :

```bash
$ docker-compose restart uwsgi worker nginx
$ docker-compose stop
```

When you do `docker-compose up -d` the application should be available at
`http://127.0.0.1/` , and if you've configured https, `https://127.0.0.1/` . If
you need to shell into the application, for example to debug with `python
manage.py shell` you can get the container id with `docker ps` and then do:

```bash
$ NAME=$(docker ps -aqf "name=sregistry_uwsgi_1")
$ docker exec -it ${NAME} bash
```

## Debugging Containers

Sometimes you might want to start containers and debug. The first thing to do is
to stop and remove old containers, and if necessary, remove old images.

```bash
$ docker-compose stop
$ docker-compose rm
```

If you want to re-pull (or for other reason, remove) the core images, do that
too:

```bash
$ docker rmi quay.io/vanessa/sregistry
$ docker rmi quay.io/vanessa/sregistry_nginx
```

You can inspect any container by looking at its logs:

```bash
$ docker-compose logs uwsgi

# Only the last 30 lines
$ docker-compose logs --tail=30 uwsgi

# Hanging
$ docker-compose logs --tail=30 -f uwsgi
```

It's also helpful to (after stopping and removing) bring up the containers but
leave out the `-d` This can commonly show issues related to starting up (and
ordering of it).

```bash
$ docker-compose up
```

And then press Control+C to kill the command and continue.

## Interactive Debugging

This is my preference for debugging - because you can shell into any of the
containers and inspect things in real-time (interactively). For example, let's
say uWSGI is running, and the container has `CONTAINER_ID` of `37b0f7d1332a` (do
`docker ps` to get the identifiers). We could shell into it via:

```bash
$ docker exec -it 37b0f7d1332a bash
```

and then find the codebase at `/code` . You can get an interactive Django shell:

```bash
$ cd /code
$ python manage.py shell
```

or an interactive database shell.

```bash
$ python manage.py dbshell
```

## Restart Containers

If you modify the container code (the Python, or a configuration value, etc.), 
you should restart the container for changes to take effect. It's good to be
conservative and only restart those containers that are needed (e.g., usually
NGINX and uWSGI).

```bash
$ docker-compose restart uwsgi nginx    # uwsgi and nginx
$ docker-compose restart                # all containers
```

## Build Containers

If you make changes to either of the images locally (or have other reasons to
build them on your own), you can do this!  In the base of the repository:

```bash
$ docker build -t quay.io/vanessa/sregistry .
```

And then to build NGINX:

```bash
$ cd nginx
$ docker build -t quay.io/vanessa/sregistry_nginx .
```

That's it! Likely the easiest thing to do is `docker-compose up -d` and let the
containers be pulled and started, and debug only if necessary. Once you have
issued the commands to generate and start your containers, it's time to read the
[setup](../setup) guide to better understand how to configure and interact with
your Singularity Registry.
---
title: Configure https for your server
pdf: true
toc: true
---

# Configure HTTPs

There are two strategies we discuss here to get https. The first is for development, meaning
you can make a faux https certificate to test it locally, and the second is for
production.

## Faux https

We are going to be using [FiloSottile/mkcert](https://github.com/FiloSottile/mkcert)
for this case. You can following the instructions in the README to generate your
"certificates." You will need [Go installed](https://golang.org/doc/install) on your system first.
If you want to install from a package manager, there are instructions in the
[repository README](https://github.com/FiloSottile/mkcert).

```bash
git clone https://github.com/FiloSottile/mkcert
cd mkcert
go build
```

I also needed to install:

```bash
$ sudo apt install libnss3-tools
```

So the certificates could be installed to my browsers.
This places `mkcert` in your current working directory!

```bash
$ ./mkcert -h
Usage of mkcert:

	$ mkcert -install
	Install the local CA in the system trust store.

	$ mkcert example.org
	Generate "example.org.pem" and "example.org-key.pem".

	$ mkcert example.com myapp.dev localhost 127.0.0.1 ::1
	Generate "example.com+4.pem" and "example.com+4-key.pem".

	$ mkcert "*.example.it"
	Generate "_wildcard.example.it.pem" and "_wildcard.example.it-key.pem".

	$ mkcert -uninstall
	Uninstall the local CA (but do not delete it).

For more options, run "mkcert -help".
```

I then did:

```bash
$ ./mkcert -install
The local CA is already installed in the system trust store! 👍
The local CA is now installed in the Firefox and/or Chrome/Chromium trust store (requires browser restart)! 🦊
```

Let's pretend we are generating a certificate for singularity-registry.org.

```bash
./mkcert singularity-registry.org "*.singularity-registry.org" singularity-registry.test localhost 127.0.0.1 ::1

Created a new certificate valid for the following names 📜
 - "singularity-registry.org"
 - "*.singularity-registry.org"
 - "singularity-registry.test"
 - "localhost"
 - "127.0.0.1"
 - "::1"

Reminder: X.509 wildcards only go one level deep, so this won't match a.b.singularity-registry.org ℹ️

The certificate is at "./singularity-registry.org+5.pem" and the key at "./singularity-registry.org+5-key.pem" ✅

It will expire on 29 August 2023 🗓
```

Then I moved them into the registry root, and updated my shub/settings/config.py to use
https on localhost.

```python
DOMAIN_NAME = "https://127.0.0.1"
DOMAIN_NAME_HTTP = "https://127.0.0.1"
DOMAIN_NAKED = DOMAIN_NAME_HTTP.replace("https://", "")
```

Finally, we need to make sure that we are using the docker-compose file for https,
the nginx.conf for https, and that the certificates are correctly bound.

```bash
mv docker-compose.yml docker-compose.yml.http
mv https/docker-compose.yml .
mv nginx.conf nginx.conf http
mv https/nginx.conf.https nginx.conf
```

In the docker-compose.yml that is newly copied, change the binds of the
paths to use the files in your present working directory.

```yaml
nginx:
  restart: always
  image: quay.io/vanessa/sregistry_nginx
  ports:
    - "80:80"
    - "443:443"
  volumes:
    - ./nginx.conf:/etc/nginx/conf.d/default.conf:ro
    - ./uwsgi_params.par:/etc/nginx/uwsgi_params.par:ro
    - ./singularity-registry.org+5-key.pem:/code/domain.key
    - ./singularity-registry.org+5.pem:/code/cert.pem
#    - /etc/ssl/certs:/etc/ssl/certs:ro
#    - /etc/ssl/private:/etc/ssl/private:ro
```

We're almost done! In the newly moved nginx.conf, make sure to comment out the dhparam.pem,
the old certificate paths, and add the new ones we just created:

```
    ssl on;
#    ssl_certificate /etc/ssl/certs/chained.pem;
#    ssl_certificate_key /etc/ssl/private/domain.key;
    ssl_certificate /code/cert.pem;
    ssl_certificate_key /code/domain.key;
    ssl_session_timeout 5m;
    ssl_protocols TLSv1 TLSv1.1 TLSv1.2;
    ssl_ciphers ECDHE-RSA-AES256-GCM-SHA384:ECDHE-RSA-AES128-GCM-SHA256:DHE-RSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-SHA384:ECDHE-RSA-AES128-SHA256:ECDHE-RSA-AES256-SHA:ECDHE-RSA-AES128-SHA:DHE-RSA-AES256-SHA:DHE-RSA-AES128-SHA;
    ssl_session_cache shared:SSL:50m;
#    ssl_dhparam /etc/ssl/certs/dhparam.pem;
    ssl_prefer_server_ciphers on;
```

You can then do `docker-compose up -d`. As instructed by the tool, we probably
need to restart our browsers.

```bash
docker-compose up -d
```

And amazingly (I opened Chrome freshly) my site has https!

![teams.png](../../assets/img/https-works.png)


Once you are done, you can (along with deleting the certs you created) you can
uninstall mkcert and delete the certs you created:

```bash
$ ./mkcert -uninstall
...
The local CA is now uninstalled from the system trust store(s)! 👋
```

You can see where the root certs were located:
 
```bash
$ ./mkcert -CAROOT
/home/vanessa/.local/share/mkcert
```

And remove them:

```bash
rm -rf  /home/vanessa/.local/share/mkcert/
```

I also deleted the mkcert directory and the original certificates,
and undid the changes above to settings, the docker-compose and nginx.conf files

## Production https

### Get a hostname

Recall that the first step to get https is to have a hostname. You can use Google
Domains, or you can also create an account on [https://www.dynu.com](https://www.dynu.com) (It’s free).
Log into your account and under the Control Panel go to DDNS Services. On the next page, click the **+ Add** button, and fill out the Host and Top Level fields under Option 1 using whatever you like. 
This will be how users access your server (e.g., `sregistry.dynu.net)`. Click + Add.

 - On the next page, change the IPv4 Address to the IP address for your droplet. Change the TTL to 60. Click Save.
 - With a few minutes, you should be able to access your server using that hostname.

### Test Nginx

In the install script, we installed nginx. Now, you merely need to start it (it might
already be started).

```bash
$ sudo service nginx start
```

For this next step, we are still working on the host where you will run your container. What we first need to do is generate certificates, start a local web server, and ping "Let's Encrypt" to verify that we own the server, and then sign the certificates.

### SSL Certificates

We'll use "certbot" to install and renew certificates.

#### Step 1. Set some variables

First we'll set some variables that are used in later steps.

```bash
EMAIL="youremail@yourdomain.com"
DOMAIN="containers.page"
```

The email you set here will be used to send you renewal reminders at 20 days, 10 days, and 1 day before expiry (super helpful!)

#### Step 2. Install certbot

Certbot automates certificate generation and renewal. In other words, it makes it really easy to setup SSL.

```bash
sudo add-apt-repository ppa:certbot/certbot
sudo apt-get update
sudo apt-get install python-certbot-nginx
```

#### Step 3. Get certificates with certbot

Now obtain a certificate by running this command.  Note that if you aren't using a container, or you aren't the root user, you might need to add `sudo`.

```bash
certbot certonly --nginx -d "${DOMAIN}" -d "www.${DOMAIN}" --email "${EMAIL}" --agree-tos --redirect
```

Equivalently, if your domain doesn't have `www.` you can remove the second `-d` argument.

#### Step 4. Stop nginx

Now we need to stop nginx because we have what we need from it!

```bash
sudo service nginx stop
```

#### Step 5. Copy certs to a new location

Now we'll move the certs to where they're expected later.

```bash
sudo cp /etc/letsencrypt/live/$DOMAIN/fullchain.pem /etc/ssl/certs/chained.pem
sudo cp /etc/letsencrypt/live/$DOMAIN/privkey.pem /etc/ssl/private/domain.key
sudo cp /etc/letsencrypt/ssl-dhparams.pem /etc/ssl/certs/dhparam.pem
```

#### Step 6. Renewal (and remembering to renew!)

Certificates expire after 90 days. You'll get reminders at 20 days, 10 days, and 1 day before expiry to the email you set before. Before the cert expires, you can run this command to renew:

```bash
sudo certbot renew
```

{% include alert.html title="Important!" content="Before renewing you need to stop the docker container running expfactory and start nginx outside of docker." %}

The commands to stop the nginx container and
renew the certificates might look like this (this is for typical Ubuntu or similar).

```bash
docker-compose stop nginx
sudo service nginx start
sudo certbot renew
sudo service nginx stop
docker-compose start nginx
```

And then issue the command to start your container.

Importantly, when you start the container (that will be generated in the next steps)
you will need to bind to these files on the host, and
expose ports 80 and 443 too.
---
title: "Installation: Using Ansible"
pdf: true
toc: true
---

# Installation using Ansible

{% include alert.html type="info" title="Important!" content="Important: this deployment works for Singularity Registry Server v1.1.23 and earlier, before Minio was added" %}

It is possible to automate the deployment of your Singularity Registry by using Ansible. There is an  Ansible role that can be used for deploying your own Singularity Registry Server and installing [sregistry-cli](https://github.com/singularityhub/sregistry-cli). This role is available [here](https://galaxy.ansible.com/grycap/singularity_registry). The complete documentation can be found [here](https://github.com/grycap/ansible-role-singularity-registry).

For example, if you want to install [sregistry](https://github.com/singularityhub/sregistry) with GitHub and PAM authentication, and also the client [sregistry-cli](https://github.com/singularityhub/sregistry-cli), you can create the following ansible playbook:

```yml
- become: true
  hosts: localhost
  vars:
    # Variables to configure GITHUB authorization
    sregistry_secrets_vars:
    - { option: 'SOCIAL_AUTH_GITHUB_KEY', value: "XXXXXXXXXX" }
    - { option: 'SOCIAL_AUTH_GITHUB_SECRET', value: "XXXXXXXXXX" }

    sregistry_config_vars:
      - { option: 'ENABLE_GITHUB_AUTH', value: True }
      - { option: 'HELP_CONTACT_EMAIL', value: 'email@email.com' }
      - { option: 'HELP_INSTITUTION_SITE', value: 'https://www.yourinstitution.com'}
      - { option: 'REGISTRY_NAME', value: 'My Singularity Registry' }
      - { option: 'REGISTRY_URI', value: 'mysreg' }

    # Use PAM authorization
    sregistry_plugins_enabled:
      - pam_auth

    # sregistry-cli in Docker
    sregistry_cli_use_docker: true

  roles:
    - { role: grycap.singularity_registry }
```

Once you have completed all the configuration of the role in your recipe, it is time to run the playbook (in this example, it is stored in `sregistry.yml`). It should be pointed out that root privileges are required to install it:

```bash
sudo ansible-playbook sregistry.yml
```

Want to learn more? Check out the [official documentation](https://github.com/grycap/ansible-role-singularity-registry) 
for the entire story.
