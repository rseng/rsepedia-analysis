# 3.0.1

- Fixed bug in migration script (there was no key `tags` in the project collection, so the `$rename` operation failed)
- Updated the sample data with the latest production data (grabbed Dec 8, 2020)
- Updated the data entry notes for project collection

# 3.0.0

commit comparison [2.0.2...3.0.0](https://github.com/research-software-directory/research-software-directory/compare/2.0.2...3.0.0)


- Added project page accessible through new route `/projects/<id>` and `/projects/<slug>`
- Added project index page accessible through new route `/projects`
- Added 404 error page
- Clean up and improve scss files
- Cleaned up the file structure for templates
- Improve documentation
- `docker-compose` now uses `.env` for environment variables instead of `export` command
- This version works on Windows and updated the developer documentation
- Added error messages if environment variables are undefined
- Migrated from Travis CI to GitHub Actions
    - Added GitHub Superlinter
    - Added OSSAR tests
    - Added Markdown link checker
    - Added `backend` tests
    - Added `frontend` tests
    - Added `harvesting` tests
    - Added integration tests
- Certificate microservice is now called `https` in docker-compose.yml instead of `nginx_ssl`
- Use static docker volumes instead of 
    - `docker-volumes/letsencrypt` 
    - `docker-volumes/cert`
- Python docker containers now based on Python 3.8
- Node docker containers now use Node 14.x
- Updated security and other dependencies
- Use `Caddy` in favor of `letsencrypt`
- API changes
    - removed `corporateUrl` and `principalInvestigator` from `project` collection
    - `image` now uses blobs like `{data: string, mimeType: string}`
    - added required properties to `project` collection
        - `callUrl`
        - `codeUrl`
        - `dateEnd`
        - `dateStart`
        - `description`
        - `grantId`
        - `isPublished`
        - `related.organizations`
        - `related.projects`
        - `related.software`
        - `slug`
        - `topics`
        - `team`
    - added optional properties to `project` collection
        - `dataManagementPlanUrl`    
        - `homeUrl`
        - `imageCaption`
    - renamed required properties to `project` collection
        - `tags` to `technologies`
- the sample data from `database/db-init` was updated to schema version 3.0.0. For notes on
  how to migrate your own data, please refer to the
  [data migration notes](/data-migration/2.0-to-3.0/README.md).
- the sample data also include data harvested from esciencecenter.nl/projects such as logos,
  related people, related organizations, descriptive text, and a hero image; furthermore the
  project start dates and end dates were added based on data from Exact. 

# 2.0.2

- Had to rewrite history due to a copyright violation on the included Akkurat font files.

# 2.0.1

- fixed a bug with retrieving images from the corporate site when they have http addresses instead of https

# 2.0.0

<!-- - Bugfix | Change | Feature | Documentation | Security -->
## Summary of changes

- no need for starting ``docker-compose`` with ``--project-name`` anymore
- updated documentation
- updated ``admin`` style for better legibility
- projects can now be harvested from the corporate site, and users can add mentions as 'impact' or 'output' via the admin interface
- security updates
- fixed bug where server logging would fill the entire disk space
- now providing feedback to page maintainers via frontend
- ``/graphs`` is now available via a button on the frontend
- added codemeta and CITATION.cff files as downloadable files
- added the content type of downloadable citation manager files
- downloadable citation manager files now have filenames consistent with their respective standard
- full commit diff [1.2.0...2.0.0](https://github.com/research-software-directory/research-software-directory/compare/e1e10fc781089d19aedc32824ffe4641f746baa2...2be41cb88be237700f60feb03fd4702e7bee9cff)

For the data migration instructions, go [here](/data-migration/1.x-to-2.x/README.md).

# 1.2.0

- added github issue templates
- replaced latest codemeta fields in releases with latest schema.org; the schema.org data is also used in the frontend to help discovery by Google et al.
- replaced the prominent header in the frontend with a more subtle one that helps make explicit that there are multiple instances of the Research Software Directory
- added new /graphs page showing metrics as well as their distribution over the software packages
- removed deprecated code from harvesting scrapers
- fixed error with blog scraping after Medium site changed its layout
- added more detailed control of the harvesting of dois and of zotero items
- added more informative logging messages for harvesters
- added throttling of queries to Zenodo
- cleaned up docker-compose file, simplified building (removed required ``-p`` option)
- added documentation for maintainers, e.g. how to make a release, how to update a production instance; added notes on security aspects

# 1.1.0

- added a simple OAI-PMH interface to allow harvesting of metadata about the 
items in the Research Software Directory in datacite4 format. The interface
implementation is not complete; at the moment, the OAI-PMH verbs ``ListRecord``
and ``GetRecords`` are implemented (but without any time based slicing such as
using ``from`` or ``until``, and without subsetting based on ``set``)
- added a service that visualizes the state of the Research Software Directory 
instance as simple graphs, e.g. histograms of how many contributors there are 
per software package.
- no longer using git submodules, this should make installing a lot easier. 
Having a monolithic repo also means that it is easier to see diffs between an 
'upstream' instance of the Research Software Direcotry, and one of its
descendants. Finally, archiving of the software via Zenodo works better for a
monorepo than for a repo with multiple submodules.
- bugfixes
- added a lot of documentation

# 1.0.0

First stable release.

Names below refer to repositories within the https://github.com/research-software-directory/ GitHub organization.

submodules SHAs:
```
 6856e70204eac7e0bb5b63cd0c391b13917966d9 admin (v1-140-g6856e70)
 fa3e0e62033dba864437af8c3b929e639428a4ac auth-github (heads/master)
 97c26c5f2825396f03fd029c32c3e21404c3b3a5 backend (v1-122-g97c26c5)
 f553e6d076471eab390b44d39c915e48f3b27ba4 db-dump (heads/master)
 1d82f242c3b2d15973b50f9a4096a16e72bcb0be frontend (v1-79-g1d82f24)
 362da73a4d594972318d5d11e3d6fb81402a46da readthedocs (heads/master)
 b4a99564651a547488c3a00efbb0b0ccbd75a7fe tasks-nlesc (heads/master)
```

# 0.0.0

First release, mostly testing the (Zenodo| GitHub) infrastructure at this point.
[![Research Software Directory](https://img.shields.io/badge/rsd-Research%20Software%20Directory-00a3e3.svg)](https://www.research-software.nl/software/research-software-directory)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1154130.svg)](https://doi.org/10.5281/zenodo.1154130)
[![Admin tests](https://github.com/research-software-directory/research-software-directory/workflows/Admin%20tests/badge.svg)](https://github.com/research-software-directory/research-software-directory/actions?query=workflow%3A%22Admin+tests%22)
[![Frontend tests](https://github.com/research-software-directory/research-software-directory/workflows/Frontend%20tests/badge.svg)](https://github.com/research-software-directory/research-software-directory/actions?query=workflow%3A%22Frontend+tests%22)
[![Backend tests](https://github.com/research-software-directory/research-software-directory/workflows/Backend%20tests/badge.svg)](https://github.com/research-software-directory/research-software-directory/actions?query=workflow%3A%22Backend+tests%22)
[![Integration Tests](https://github.com/research-software-directory/research-software-directory/workflows/Integration%20Tests/badge.svg)](https://github.com/research-software-directory/research-software-directory/actions?query=workflow%3A%22Integration+tests%22)
[![Check Markdown links](https://github.com/research-software-directory/research-software-directory/workflows/Check%20Markdown%20links/badge.svg)](https://github.com/research-software-directory/research-software-directory/actions?query=workflow%3A%22Check+Markdown+links%22)

# Research Software Directory

## What is it?

The Research Software Directory is a content management system that is tailored to research software.

The idea is that institutes for whom research software is an important output, can run their own instance of the Research Software Directory. The system is designed to be flexible enough to allow for different data sources, database schemas, and so on. By default, the Research Software Directory is set up to collect data from GitHub, Zenodo, Zotero, as well as Medium blogs.

For each software package, a _product page_ can be created on the Research Software Directory if the software is deemed useful to others.

## What the Research Software Directory can do for you

The Research Software Directory:

1. presents software packages alongside the context necessary for visitors to understand how the software can help them
1. makes scientific impact of research software visible in a qualitative way
1. provides automatically generated citation metadata in a variety of reference manager file formats, for easy citation
1. improves findability of software packages by applying Search Engine Optimization techniques such as schema.org metadata. This helps search engines understand what a given software package is about, thus improving ranking of search results
1. provides aggregated insights through a metrics dashboard, helping to make more accurate and more timely business decisions
1. provides metadata about its software packages via [OAI-PMH](https://www.openarchives.org/pmh/), the standard protocol for metadata harvesting. Digital libraries and other services can use this feature to automatically update their records with data about the software packages published in the Research Software Directory.
1. provides all of its data via a JSON API
1. integrates with third-party services such as [Zotero](http://zotero.org/) (reference manager), [Zenodo](https://zenodo.org/) (archiving), GitHub (code development platform)

## Examples

1. [https://research-software.nl](https://research-software.nl)
1. [https://software.process-project.eu](https://software.process-project.eu)


## Try it out

### Requirements

1. You'll need a minimum of about 3 GB free disk space to
store the images, containers and volumes that we will be making.
1. Linux OS (we use Ubuntu 18.04)
1. [docker](https://docs.docker.com/install/) (v19.03 or later)
1. [docker-compose](https://docs.docker.com/compose/install/) (v1.26 or later)
1. [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) (v2.21 or later)

### Instructions

To quickly get a running Research Software Directory up and running on your local machine do the following

1. Fork this repo to your own GitHub organization or GitHub profile and clone it
1. Start the complete stack using

    ```shell
    cp rsd-secrets.env.example .env
    docker-compose build
    docker-compose up
    ```
Go to [http://localhost](http://localhost) (disregard certificate warning) to see the Research Software Directory
website. You should be able to see all non-authenticated pages, but editing data or harvesting data from external
sources won't work. To bring up the website with all bells and whistles, refer to selected resources from the list
below.

## Further resources

<!-- see also docs/README.md -->

1. [Entering data about your software in an existing instance](./docs/entering-data.md)
1. [Configuring your instance to use your own data sources](./docs/configuring-your-instance-to-use-your-own-data-sources.md)
1. [Changing the look and feel](./docs/changing-the-look-and-feel.md)
1. [Hosting your instance online](./docs/hosting.md)
1. [Running an instance of the Research Software Directory in production](./docs/production.md)
1. [Finding your way: Research Software Directory services](./docs/services-overview.md)
1. [Documentation for developers](./docs/documentation-for-developers.md)
1. [Documentation for maintainers](./docs/documentation-for-maintainers.md)
1. [Security concerns](./docs/security.md)
1. [Contributing](./.github/CONTRIBUTING.md)
# `harvesting` service for the Research Software Directory

Periodic harvesting of external data

## Install

```shell
pip install -r requirements.txt
```

## Configuration

Harvesting requires environment variables.
See `.env.example` for the required variables.

### Harvesting mentions data from Zotero

To fetch mentions of software from Zotero you need the Zotero group identifier
to search and an API key to access it.

The group identifier can be found in the url of a group on the
[https://www.zotero.org/groups/](https://www.zotero.org/groups/) page. For example for the Netherlands eScience
center the group url is
[https://www.zotero.org/groups/1689348/netherlands_escience_center](https://www.zotero.org/groups/1689348/netherlands_escience_center) and the group
identifier is 1689348. The Zotero group identifier must be set as value for the
`ZOTERO_LIBRARY` environment variable.

The API key can be created at [https://www.zotero.org/settings/keys](https://www.zotero.org/settings/keys) The API key
must be granted read only permission on the group. The Zotero API key must be
set as value for the `ZOTERO_API_KEY` environment variable.

To save the mentions into the database the location of the backend server and a
jwt secret must be configured. A jwt secret can be found in the `JWT_SECRET`
env var of the backend server. The jwt secret must be set as value for the
`JWT_SECRET` environment variable. The location of backend server must be set
as value for the `BACKEND_URL` environment variable.

### Harvesting commits data from Github

To fetch GitHub commits for lots of repositories an access token is required.
The token can be generated at [https://github.com/settings/tokens](https://github.com/settings/tokens), no scopes need
to be selected. The token must be set as value for the `GITHUB_ACCESS_TOKEN`
environment variable.

The harvester directly injects into the database so the `DATABASE_PORT` and
`DATABASE_NAME` environment variables should be set.

### Harvesting citations data from Zenodo

The harvester retrieves all the data to be able to cite a software item.

The harvester directly injects into the database so the `DATABASE_HOST`, the
`DATABASE_PORT` and `DATABASE_NAME` environment variables should be set.

## Usage

Refer to the help like this:

```shell
python app.py --help
```

The mentions can be fetched from Zotero using

```shell
python app.py harvest mentions
python app.py harvest mentions --help
python app.py harvest mentions --since-version VERSION
python app.py harvest mentions --keys STRING
```

The Github commits can be fetched using

```shell
python app.py harvest commits
```

The releases of each software can be fetched using

```shell
python app.py harvest citations
python app.py harvest citations --help
python app.py harvest citations --dois STRING
```

The metadata of each software can be fetched using

```shell
python app.py harvest metadata
python app.py harvest metadata --help
python app.py harvest metadata --dois STRING
```

If you want to fetch data from all sources, use:

```shell
python app.py harvest all
```

Once the data have been fetched, data from different collections need to be
combined into one document such that it can be fed to the `frontend`'s template.
Combining is done with the `resolve` task, as follows:

```shell
# resolve references in all documents of type project
python app.py resolve projects
```

```shell
# resolve references in all documents of type software
python app.py resolve software
```

```shell
# resolve references in all documents of type project and in all
# documents of type software
python app.py resolve all
```


## Dockerized

Build with

```shell
docker build --tag rsd/harvesting .
```
# Admin interface

This is the admin user interface for the Research Software Directory. It can be used to add and update software items to
the directory. The admin interface requires credentials for users to see the admin interface.

TODO:

- add notes on how to keep the dependencies updated, e.g. using `yarn audit` and `yarn outdated` ([#367](https://github.com/research-software-directory/research-software-directory/issues/367))

## Development setup

This project was bootstrapped with [Create React App](https://github.com/facebookincubator/create-react-app) using
scripts package [@nlesc/react-scripts](https://github.com/NLeSC/create-react-app). Original documentation
[here](https://github.com/NLeSC/create-react-app/blob/master/packages/react-scripts/template/README.md).

### Prerequisites

1. Install node version `^14.0`, e.g. using `nvm` (Node Version Manager, see [https://github.com/nvm-sh/nvm#install-script](https://github.com/nvm-sh/nvm#install-script)):
    1. `curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.36.0/install.sh | bash`
    1. restart terminal
    1. run `nvm install 14`
1. Install yarn, `npm install -g yarn` (Note: there exists a similarly named package, Apache Hadoop YARN --please ignore
    that, it's something else)

### Getting a local development build running

1. To login in to admin interface you need to set the authentication callback to `http://localhost:8000/auth/get_jwt` on [https://github.com/settings/developers](https://github.com/settings/developers)
1. Start the services required by admin interface:

    ```shell
    cd research-software-directory
    docker-compose -f docker-compose.yml -f admin/docker-compose.dev.yml up
    ```

1. New terminal
1. `cd admin`
1. Install the `admin` service's dependencies: `yarn install`
1. Start the admin service in a development server: `yarn start`. It will tell you where to go to check the `admin`
    interface ([http://localhost:8000](http://localhost:8000)).

The Docker containers started with `docker-compose` is a complete Research Software Directory instance on [https://localhost](https://localhost) including an admin interface. Any changes made to software/projects/etc. in the admin interface running on [http://localhost:8000](http://localhost:8000) will be saved in the instance and can be viewed there.

## Production setup: non-dockerized

To build the app for production run:

```bash
yarn build
```

The deployable app will build to the `./build/` directory.

## Production setup: dockerized

The Docker image should not be used on its own, as the code expects the [backend server](/backend) to be running at
`/api` and the [auth server](/auth-github) to be running at `/auth`.

The `rsd/admin` image should be used as part of a `docker-compose`, see
[https://github.com/research-software-directory/research-software-directory](https://github.com/research-software-directory/research-software-directory):

```bash
cd research-software-directory
docker-compose build admin
```
# `backup` service for the Research Software Directory
# `graphs` service for the Research Software Directory

Present integrated statistics and distributions about this instance of the Research Software Directory.
# Placeholder directory for the Research Software Directory's bind mounts

This directory is a placeholder for any data volumes that Docker creates upon
running the software.

After running the complete stack, this directory should look like this:

```shell
docker-volumes/
├── db
├── oaipmh-cache
└── README.md
```

Note that in addition to these bind mounts, other data are stored in static
volumes. Refer to [docker-compose.yml](/docker-compose.yml) to see which volumes
exist.
# data migration

When upgrading, refer to the instructions in the subdirectories below; you may need to follow multiple sets of instructions in sequence.

- [2.0-to-3.0](2.0-to-3.0)
- [1.x-to-2.x](1.x-to-2.x)
# Updating data from 2.0.2 to 3.0.x

The following code snippet emulates the situation where the code is in version 3.0.x while the data is in version 2.0.2.

```shell
# update the references
git fetch --tags

# give me the code state for 3.0.0 ...
git checkout 3.0.0
```

If you don't have a database filled with 2.x version data, it can be filled with sample data from version 2.0.2 using

```shell
git checkout 2.0.2 -- database/db-init
```

but for the 2.0.2 sample data to get used, you'll need to remove the current database files from `docker-volumes/db` as
well, see instructions [here](https://github.com/research-software-directory/research-software-directory/blob/3.0.0/docs/dev.md#removing-local-state). Note this has the potential for **THE LOSS OF
DATA**.

It should now be possible to upgrade the sample data as follows:

```shell
docker-compose build
docker-compose up -d
docker-compose log --follow
```

If you opted to fill the database with 2.0.2 sample data by doing the `git checkout 2.0.2 -- database/db-init` above,
you can roll back those changes once the Research Software Directory is up (and has had time to restore the database
files from `db-init`) using:

```shell
git reset HEAD database/db-init
git checkout -- database/db-init
```

Verify that the API is serving v2.0.2 data by visiting [http://localhost/api/project/764](http://localhost/api/project/764).

Then, in a new terminal,

```shell
# copy the migrate script to inside the running database service
docker cp data-migration/2.0-to-3.0/migrate.js $(docker-compose ps -q database):/tmp

# run the migrate script
docker-compose exec database mongo rsd /tmp/migrate.js

# update the cache
docker-compose exec harvesting python app.py resolve all
```

The data you get from the API should now be according to the 3.0.x schema, e.g.
[http://localhost/api/project/764](http://localhost/api/project/764), and all aspects of the site should now work. You
should verify if everything works by doing the checks mentioned in section [Verifying the local
installation](https://github.com/research-software-directory/research-software-directory/blob/3.0.0/docs/dev.md#verifying-the-local-installation).

## Optional: get project data

Optionally update the project data from [esciencecenter.nl corporate website](https://esciencecenter.nl) by doing following (continue in the
terminal from the previous section):

```shell
python3 -m virtualenv -p python3 venv3
source venv3/bin/activate
pip install -r ./data-migration/2.0-to-3.0/requirements.txt
source .env
export PYTHONWARNINGS="ignore:Unverified HTTPS request"
export PYTHONPATH=$PYTHONPATH:`pwd`/harvesting
python3 data-migration/2.0-to-3.0/harvest_project_info_nlesc.py

# again, update the cache
docker-compose exec harvesting python app.py resolve all
```

The data you get from the API should now include richer data for the projects, e.g.
[http://localhost/api/project/764](http://localhost/api/project/764), which means that the corresponding project pages
are also richer, e.g. [http://localhost/projects/764](http://localhost/projects/764). All aspects of the site should
now work. You should verify if everything works by doing the checks mentioned in section [Verifying the local
installation](https://github.com/research-software-directory/research-software-directory/blob/3.0.0/docs/dev.md#verifying-the-local-installation).

## Optional: add project start and end dates

```shell
# copy the script to inside the running database service
docker cp data-migration/2.0-to-3.0/add-project-dates.js $(docker-compose ps -q database):/tmp

# run the migrate script
docker-compose exec database mongo rsd /tmp/add-project-dates.js

# update the cache
docker-compose exec harvesting python app.py resolve all
```
# Updating data from 1.x to 2.x

In version 2.0.0, the ``project`` collection is partly filled by harvesting from
an external data source, and partly filled by means of users making edits in the
admin interface. This means that version 2.0.0 of the Research Software
Directory requires changes to the database. Below are the steps to migrate data
from 1.2.0 to 2.0.0. Furthermore, the frontend now shows information for page
maintainers, for which a new MongoDB collection ``logging`` is needed.

When migrating data there is always the possibility of **LOSS OF DATA**. Review the
notes on how to make a backup of the Mongo data [here](https://github.com/research-software-directory/research-software-directory/tree/2.0.2/README.md#updating-a-production-instance).

```shell
source rsd-secrets.env
docker-compose exec database mongo rsd
```

Create collection "logging":

```shell
db.createCollection("logging")
```

**Remove** all ``release`` documents, ``project`` documents, and ``project_cache`` documents entirely:

```shell
db.release.deleteMany({})
db.project.deleteMany({})
db.project_cache.deleteMany({})
```
Then, update the project identifiers as used in the ``software`` collection by
copy-pasting the contents of the [data migration script](https://github.com/research-software-directory/research-software-directory/blob/2.0.2/data-migration-1.x-to-2.js) into the Mongo shell.

Exit the Mongo shell with Ctrl-d or ``exit``, then run the harvester:

```shell
docker-compose exec harvesting python app.py harvest all
```

See if it all worked by running:

```shell
docker-compose exec harvesting python app.py resolve all
```

The ``resolve`` command should list only INFO messages, not ERROR messages.
# Data API for Research Software Directory

Data API for Research Software Directory

## Requirements

- Python 3.8+ (or run through Docker)

## Setup

This service depends on the following services:
- MongoDB (3.6) (`docker pull mongo:3.6 && docker run -p 27017:27017 mongo:3.6`)


## Configuration

Configuration consists of two parts:

### Environmental variables

```shell
JWT_SECRET=[hidden]                - JSON web token secret to generate/verify tokens]
DATABASE_HOST=localhost            - MongoDB host
DATABASE_PORT=27017                - MongoDB tcp port
DATABASE_NAME=rsd                  - MongoDB database to use
SCHEMAS_PATH=./schemas             - Path where the schema files can be found
```

### Schema files

Schema files should be in `SCHEMAS_PATH`. See the directory `schemas_example`
for the schemas we use at the eScience Center.

## Run unit tests

```shell
mkvirtualenv data-api -p `which python3`
source data-api/bin/activate
pip install -r requirements.txt
PYTHONPATH=`pwd` pytest
```

## Run API server

Make sure that your Python is up to date and requirements are installed (same as under unit tests).
Set environmental variables (eg. `export $(cat .env.example | xargs)`).

```shell
python entry.py                                                               # starts server
FLASK_APP=`pwd`/entry.py flask generate_jwt --sub test_user -p write -p read  # generates a JWT for read+write
```

Or run through Docker:

```shell
docker build -t rsd/backend .
docker run --env-file ./.env -p 5001:8000 -it --name rsd-backend rsd/backend
```

## Usage

### `GET /[resource_type]`

Get list of resources of type `resource_type`. Arguments:

- ```shell
  ?sort=[field]
  ?sort=[field.subfield]
  ```
  Sorts by field `field`

- ```shell
  ?direction=desc
  ```
  Sort descending (default is ascending)

- ```shell
  ?skip=[skip]&limit=[limit]
  ```
  Skip first `skip` records, show max `limit` results.

### `GET /[resource_type]/[id_or_slug]`

Returns record with `slug` or `primaryKey/id` of `id_or_slug`

### `POST /[resource_type](?test=1)`

Create a new record. Body should contain the contents in JSON format.
Argument `test=1` means it won't be really saved, just tested if it can be saved.
**Requires `write` permission, use `Authorization: Bearer [JWT]` header with a valid JWT with write permissions**.

### `PATCH /[resource_type]/[id_or_slug]`

Updates record. Body should contain the (partial) contents in JSON format.
A field is not modified if it is not set in the JSON body.
**Requires `write` permission.**
# Contributor Covenant Code of Conduct

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
reported by contacting the project team at [rsd@esciencecenter.nl](mailto:rsd@esciencecenter.nl). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# How to contribute

Thank you for considering to contribute to the Research Software Directory!

There are many ways to contribute, e.g.:

- Write blog posts, tutorials, etc. about the Research Software Directory. Heck, even tweeting about it helps!
- Review any and all parts of the code base including documentation, schemas or anything else you find in there
- Submit bug reports, feature requests; make suggestions
- Correct typos (a one-character diff in a Pull Request is perfectly fine)
- Fork the repository and set up your own custom instance.

However, when participating, we ask that you adhere to some **ground rules**.

Your contribution to the Research Software Directory is valued, and it should be an enjoyable experience.
To ensure this, there is the Research Software Directory
[Code of Conduct](/.github/CODE_OF_CONDUCT.md), which you are required to follow.

If you have never contributed to an open source project, you may find this tutorial helpful:
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github).
# Pull request details

## List of related issues or pull requests

Refs: #ISSUE_NUMBER


## Describe the changes made in this pull request

<!-- include screenshots if that helps the review -->


## Instructions to review the pull request
---
name: Bug report
about: Something doesn't work like I expected.
title: ''
labels: bug
assignees: ''

---

Use your best judgment to provide a useful level of information. Depending on the nature of the issue, consider including, e.g.

- Which operating system you're running
- Which version of the software you're running
- Any logging or error output you see
---
name: Feature request
about: I'd like to see something added or changed to make the Research Software Directory
  even better.
title: ''
labels: enhancement
assignees: ''

---


---
name: I found a security issue!
about: There is an issue with the software's security
title: ''
labels: ''
assignees: ''

---

Please do not proceed to file this issue publicly; instead send an email to ``rsd at esciencecenter dot nl``.
# auth-github

Github authentication service for the Research Software Directory. It should be configured as the callback endpoint of a
GitHub OAuth app.

This verifies the user is part of Github organization `AUTH_GITHUB_ORGANIZATION`,
then forwards the user to `[AUTH_CALLBACK_URL]?jwt=[GENERATED_JWT]` (or an
error message is shown).

The JWT is generated using secret `JWT_SECRET`.

Claims in the generated token:

```json
{
  "sub": [github_username],
  "subType": "GITHUB",
  "permissions": ["read", "write"],
  "iat": [current timestamp (in seconds)],
  "user": {
    "name": [github_username],
    "image": [github_users_avatar_image_url],
  }
}
```

## Install

```shell
mkvirtualenv data-api -p `which python3` # create virtual python environment
source data-api/bin/activate             # activate environment
pip install -r requirements.txt          # install python dependencies
```

## Run

Copy `rsd-secrets.env.example` to `.env` and check the values.

```shell
export $(cat .env | xargs)       # load settings as env vars
python app.py                    # run debug server on default port (5002)
```

or use docker (you can skip the `install` instructions):

```shell
docker build . -t rsd/auth
docker run --env-file ./.env -p 5002:8000 -it --name auth-github auth-github
```

## Configuration

Configuration is done through environmental variables

```shell
# The secret used to generate jwt token that is used by the api backend to authenticate/authorize users from the admin site
JWT_SECRET=
# The id and secret of a GitHub OAuth app, for authorization callback URL in the GitHub application use `<url of server>/get_jwt`
AUTH_GITHUB_CLIENT_ID=
AUTH_GITHUB_CLIENT_SECRET=
# The GitHub organization that users must be a member of to be authorized
AUTH_GITHUB_ORGANIZATION=nlesc
# The url to redirect the user to when the authentication has been completed
AUTH_CALLBACK_URL=http://localhost:3000
```

## Docker

The Docker image can be built using

```shell
docker build -t rsd/auth .
``
# Make your instance available to others by hosting it online (deployment)

Amazon Web Services (AWS) is a online service provider that offers all kinds of
services relating to compute, storage, and hosting. The Netherlands eScience
Center uses AWS to run their instance of the Research Software Directory. This
section describes how to deploy your own customized instance of the Research
Software Directory to AWS.

Go to [https://aws.amazon.com/console/](https://aws.amazon.com/console/). Once there, you'll see something like:

[![AWS Management Console login](/docs/images/aws-management-console-login.png)](/docs/images/aws-management-console-login.png)

Create a free account if you don't already have one, and subsequently click
``Sign In to the Console``.

Once in the console, you'll be presented with an overview of all the services
that Amazon Web Services has to offer:

[![AWS Management Console Services Overview](/docs/images/aws-management-console-services-overview.png)](/docs/images/aws-management-console-services-overview.png)

It's easy to get lost in this plethora of services, but for running an instance
of the Research Software Directory, you'll only need 3 of them:

1. **EC2**: this is where we will run your customized instance of the Research
Software Directory and host it online; [jump to the EC2 section](#configuring-ec2)
1. **IAM**: we use this to create a user with limited privileges, so we don't
have to use root credentials when we don't have to; [jump to the IAM section](#configuring-iam)
1. **S3**: this is where we will store our daily backups; [jump to the S3 section](#configuring-s3)

## Configuring EC2

In the ``All Services`` overview, click ``EC2`` or use this link
[https://console.aws.amazon.com/ec2](https://console.aws.amazon.com/ec2).

<!-- TODO how to configure default zone -->

1. Click the blue ``Launch instance`` button
1. Scroll down to where it says ``Ubuntu Server 18.04 LTS``, click ``Select``
1. Choose instance type ``t2.small``
1. Proceed in the wizard until you get to 'Add storage'. Set the storage to 10GB.
1. Proceed in the wizard by clicking ``Next`` until you get to ``Configure
Security Group``. It should already have one rule listed. However, its security
settings should be a bit more secure, because currently it allows SSH
connections from any IP. Click the ``Source`` dropdown button, select ``My IP``.
1. Now click the blue ``Review and Launch`` button in the lower right corner
1. In the ``Review`` screen,  click the blue ``Launch`` button in the lower
right corner to bring the instance up
1. In the ``Keypair`` popup, select ``Create a new key pair``, try to give it a
meaningful name, e.g. ``rsd-instance-on-aws`` or something
1. Click ``Download Key Pair``, save the ``*.pem`` file in ``~/.ssh`` on your
local machine, then click ``Launch Instances`` (it takes a moment to
initialize).
1. On your local machine, open a terminal and go to ``~/.ssh``. Change the
permissions of the key file to octal 400 (readable only by user):

    ```shell
    chmod 400 <the keyfile>
    ```
1. Verify that the ``.ssh`` directory itself has octal permission 700 (readable,
writable, and executable by user only).
1. Go back to Amazon, click ``View instances``
1. Make a note of your instance's public IPv4, e.g. ``3.92.182.176``
1. On your own machine use a terminal to log in to your instance
1. ``ssh -i path-to-the-keyfile ubuntu@<your-instance-public-ip>``
1. Once logged in to the remote machine, update the package manager's list of
   software packages and their versions:

    ```shell
    sudo apt update
    ```

1. Upgrade any software packages to a higher version if available:

    ```shell
    sudo apt upgrade
    ```

1. Install ``docker`` and
``docker-compose``, then add user ``ubuntu`` to the group ``docker``, same as
before (see section _Documentation for developers_
[above](/README.md#documentation-for-developers)).
1. Make a new directory and change into it:

    ```shell
    cd ~
    mkdir rsd
    cd rsd
    ```
1. The machine should have ``git`` installed, use it to ``git clone`` your
customized Research Software Directory instance into the current directory as
follows:

    ```shell
    git clone https://github.com/<your-github-organization>/research-software-directory.git .
    ```
    (Note the dot at the end)

1. Open a new terminal and secure-copy your local ``rsd-secrets.env`` file to
the Amazon machine as follows:

    ```shell
    cd <where rsd-secrets.env is>
    scp -i path-to-the-keyfile ./rsd-secrets.env \
    ubuntu@<your-instance-public-ip>:/home/ubuntu/rsd/rsd-secrets.env
    ```
1. On the remote machine, create the symlink named `.env` and have it point to the secrets file:

    ```shell
    ln -s rsd-secrets.env .env
    ```

1. Follow the instructions
[above](/README.md#auth_github_client_id-and-auth_github_client_secret) to make
a second key pair ``AUTH_GITHUB_CLIENT_ID`` and ``AUTH_GITHUB_CLIENT_SECRET``.
However, let this one's ``Authorization callback url`` be ``https://`` plus your
instance's IPv4 plus ``/auth/get_jwt``. Update the Amazon copy of
``rsd-secrets.env`` according to the new client ID and secret.
1. Start the Research Software Directory instance with:

    ```shell
    cd ~/rsd
    docker-compose build
    docker-compose up -d
    ```
1. On your local machine, open a new terminal. Connect to the Amazon instance,
run the harvesters, and resolve the foreign keys:

    ```shell
    ssh -i path-to-the-keyfile ubuntu@<your-instance-public-ip>
    cd ~/rsd
    docker-compose exec harvesting python app.py harvest all
    docker-compose exec harvesting python app.py resolve all
    ```

At this point we should have a world-reachable, custom instance of the Research
Software Directory running at ``https://<your-instance-public-ip>/``. However,
if we go there using a browser like Firefox or Google Chrome, we get a warning
that the connection is not secure.

To fix this, we need to configure the security credentials, but this in turn
requires us to claim a domain and configure a DNS record. There are free
services available that you can use for this, e.g. [https://noip.com](https://noip.com). Here's how:

1. Go to [https://noip.com](https://noip.com), sign up and log in.
1. Under My services, find ``DNS Records``
1. Click the ``Add a hostname`` button
1. Choose your free (sub)domain name, e.g. I chose ``myrsd.ddns.net``
1. Fill in the IP address of your Amazon machine. In my case,
[https://myrsd.ddns.net](https://myrsd.ddns.net) will serve as an alias for [https://3.92.182.176](https://3.92.182.176)
1. Once you have the (sub)domain name, update ``DOMAIN`` and ``SSL_DOMAINS`` in the file
``rsd-secrets.env`` on your Amazon instance (leave out the ``https://`` part, as
well as anything after the ``.com``, ``.nl``, ``.org`` or whatever you may
have).
1. Fill in your e-mail for ``SSL_ADMIN_EMAIL``.
1. Finally, revisit your OAuth app here [https://github.com/settings/developers](https://github.com/settings/developers),
replace the Amazon IP address in the ``Authorization callback url`` with
your freshly minted domain name.
1. Now, stop the Research Software Directory if it is still running with Ctrl-c
or ``docker-compose stop``.
1. Start the Research Software Directory back up

    ```shell
    cd ~/rsd
    docker-compose up -d
    ```
1. Pointing your browser to your (sub)domain name should now show your instance
of the Research Software Directory (although be aware that sometimes it takes a
while before the domain name resolves to the IP address.

## Configuring IAM

1. In the ``All Services`` overview, click ``IAM`` or use this link
[https://console.aws.amazon.com/iam](https://console.aws.amazon.com/iam).
1. In the menu on the left, click ``Groups``.
1. Click the ``Create New Group`` button.
1. Name the group ``s3-users``.
1. When asked to attach a (security) policy, use the search bar to find
``AmazonS3FullAccess`` and check its checkbox.
1. Click the ``Next step`` button in the lower right corner.
1. Review your group, go back if need be. When you're ready, click the ``Create
Group`` button in the lower right corner.
1. Now you should be presented with a group, but the group is still empty; there
are no users.
1. In the menu on the left, click ``Users``.
1. Click the ``Add user`` button in the top left corner.
1. Choose your user name. I chose to call mine ``rsd-backup-maker``. For this
user, check the checkbox labeled ``Programmatic access``. This user won't need
``AWS Management Console access``, so leave that box unchecked.
1. In the lower right corner, click the ``Next: Permissions`` button.
1. Select ``Add user to group``, and make user ``rsd-backup-maker`` a member of
group ``s3-users``.
1. In the lower right corner, click the ``Next: Tags`` button. We don't need to
assign any tags, so proceed to the next page by clicking ``Next: Review``. Go
back if you need to, but if everything looks OK, click ``Create User``. You will
be presented with the new user's credentials. Download the CSV file now; we'll
use the ``Access key ID`` and the ``Secret access key`` later to set up the
backup mechanism.

## Configuring S3

In the ``All Services`` overview, click ``S3`` or use this link
[https://console.aws.amazon.com/s3](https://console.aws.amazon.com/s3).

1. create a bucket with a random name (bucket names must be globally unique;
websites like [https://www.random.org/strings/](https://www.random.org/strings/) are useful to get a random string)
1. in that bucket, make a directory, e.g. ``rsd-backups``
1. The backup service contains a program
([xenon-cli](https://github.com/xenon-middleware/xenon-cli)) that can copy to a range of
storage providers. You can use it to make daily backups of the MongoDB database,
and store the backups on Amazon's S3. For this, configure the environmental
variable ``BACKUP_CMD`` as follows (naturally, you'll need to use a different
location, username, and password; see explanation below):

    ```shell
    BACKUP_CMD='xenon filesystem s3 \
    --location http://s3-us-west-2.amazonaws.com/nyor-yiwy-fepm-dind/ \
    --username AKIAJ52LWSUUKATRQZ2A \
    --password xQ3ezZLKN7XcxIwRko2xkKhV9gdJ5etA4OyLbXN/ \
    upload rsd-backup.tar.gz /rsd-backups/rsd-backup-$BACKUP_DATE.tar.gz'
    ```

    - The bucket name is ``nyor-yiwy-fepm-dind``. It is physically located in
    zone ``us-west-2``.
    - We access the bucket as a limited-privileges IAM user, for whom we
    created an access key (it has been deactivated since). The Access key ID is
    ``AKIAJ52LWSUUKATRQZ2A``, and its corresponding Secret access key is
    ``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``.
    - The variable ``BACKUP_DATE`` is set by the backup script (see
    [``/backup/backup.sh``](/backup/backup.sh)); no need to change this for your application.
    - ``rsd-backup.tar.gz`` is the name of the backup archive as it is called
    inside the container; no need to change this for your application.
    - ``/rsd-backups/rsd-backup-$BACKUP_DATE.tar.gz`` is the path inside
    the bucket. It includes the date to avoid overwriting previously existing
    archives; no need to change this for your application.
1. Test the setup by stopping the Research Software Directory on Amazon, by

    ```shell
    # ssh into the remote machine
    cd rsd
    docker-compose stop
    # update BACKUP_CMD by editing the rsd-secrets.env file
    docker-compose up -d
    ```

    Wait until the Research Software Directory is up and running again, then

    ```shell
    docker-compose exec backup /bin/sh
    /app # /bin/sh backup.sh
    ```
# How to populate the Research Software Directory with data

This document is an instruction to Research Software Engineers working at the
Netherlands eScience Center. It describes what they need to do in order to help
make our Research Software Directory an attractive website.

First read [this
blog](https://blog.esciencecenter.nl/the-research-software-directory-and-how-it-promotes-software-citation-4bd2137a6b8) for the
overall picture of what a Research Software Directory is and what we are trying
to accomplish with it.

## Accounts

To populate the Research Software Directory with your content, you'll need
accounts for GitHub and for Zotero.

1. On GitHub, you need to be a member of the NLeSC organization, and that
   membership should be public. [This website](https://github.com/orgs/NLeSC/people)
   will tell you if you meet those criteria.
1. On Zotero, you'll need to be a member of the Netherlands eScience Center group,
   which you can request here [https://www.zotero.org/groups/1689348](https://www.zotero.org/groups/1689348).
   A group administrator needs to grant you access, so be aware it may take a while.

## Basic workflow

The blog mentions that the Research Software Directory harvests various types of
data, which are then combined using its Admin interface, as follows:

![/docs/images/data-sources.svg](/docs/images/data-sources.svg)

So, now the question is: what does this all mean for how you do your day-to-day
work? The answer is, as long as you use public source code repositories on
GitHub and make a release every so often following the guidelines outlined
[here](https://guide.esciencecenter.nl/#/citable_software/making_software_citable),
you'll be OK. By simply adhering to these best practices, the source code
related boxes (``GitHub``, ``Zenodo`` and ``CITATION.cff``) are already covered.

However, product pages are not only about source code; they also aim to provide
context. An important part of this are the _mentions_ that you'll see on most
product pages. You can control which mentions should appear on a product page by
selecting them via the dropdown list in the Admin interface of the Research Software Directory.

### Output and impact

The items on the dropdown list are harvested periodically from the
[Netherlands eScience Center group on
Zotero](https://www.zotero.org/groups/1689348). For reference, Zotero is the
place where we keep track of NLeSC's _output_ and _impact_:

- Output includes all the things we produced, be it code, data, videos,
documentation, papers, posters, demos, etc.
- Impact includes things we did not produce but that makes use of things we
produced, e.g. a paper, data or software created by an external party which uses software
developed by us

So, if you have a new blog post, new paper, or a new screencast that you want to
add as a mention to a product page, you need to add the new item to Zotero (and
then wait a bit, see section [_How/when do I get to see my
changes?_](#howwhen-do-i-get-to-see-my-changes) below). For detailed guidance on
how to add various types of records to Zotero, go [here](zotero.md).

Hopefully this clarifies where the various data are coming from. The next step
is to establish links between the parts, e.g. linking a paper to a software.
This is done by hand via the Admin interface (see sections below).

## What software should be added?

As a rule of thumb, software should be added to the Research Software Directory if

1. we can realistically claim (co-)ownership of the software AND
2. the software could be useful to somebody else

Under this definition, the following is not recorded in the Research Software Directory:

- A script that was just meant to wrangle some data unique to a problem we were having (violates 2).
- A contribution to an existing package (violates 1).
- A package that is very useful but that NLeSC wasn't involved with (violates 1).

## How do I add a new ``Software``?

1. Go to the Research Software Directory's admin interface
([https://www.research-software.nl/admin/](https://www.research-software.nl/admin/)).
1. In the left pane, select ``Software``
1. Check that the software you want to add doesn't exist yet by searching via the
search box. If the search comes up empty, add the new software by clicking the
blue ``+`` symbol.
1. You should now see a form. Each field has a brief explanation of what you should
fill in. Additionally, there is formal form validation, e.g. to make sure that
URLs are valid, and to make sure you don't accidentally skip required fields.
Complete the form. Don't worry about mistakes, you can always come back later
and fix it.
1. Set the ``isPublished`` slider near the top of the page to the right to have it
included in [research-software.nl](https://research-software.nl).
1. **Don't forget** to click ``Save`` when you're done.
1. Refer to section [_Delays_](#delays) and section [_The query
trick_](#the-query-trick) to know when you get to see your changes.

## How do I add a new ``Person``?

1. Go to the Research Software Directory's admin interface
([https://www.research-software.nl/admin/](https://www.research-software.nl/admin/)).
1. In the left pane, select ``Person``
1. Check that the person you want to add doesn't exist yet by searching via the
search box. If the search comes up empty, add the new person by clicking the
blue ``+`` symbol.
1. Fill the various parts of the person's name
1. Optionally, fill in the person's email address
1. Optionally, provide an image of the person
1. **Don't forget** to click ``Save`` when you're done.
1. Refer to section [_Delays_](#delays) and section [_The query
trick_](#the-query-trick) to know when you get to see your changes.

## How do I add a new ``Mention``?

You can't. The list of Mentions is harvested via Zotero's API. Check
[https://research-software.nl/schedule](https://research-software.nl/schedule) and look for ``python app.py harvest
mentions`` to determine when the harvester is scheduled to run.

## How do I add a new ``Project``?

1. Go to the Research Software Directory's admin interface
([https://www.research-software.nl/admin/](https://www.research-software.nl/admin/)).
1. In the left pane, select ``Project``
1. Check that the project you want to add doesn't exist yet by searching via the
search box. If the search comes up empty, add the new project by clicking the
blue ``+`` symbol.
1. Fill the form
   1. Please respect the required image size of 1024x600
   1. Add mentions to a project, either as "output" (the mention was
      produced with NLeSC funding), or as "impact" (the mention was produced without
      NLeSC funding, for example when a paper cites code that NLeSC created).
1. **Don't forget** to click ``Save`` when you're done.
1. Refer to section [_Delays_](#delays) and section [_The query
trick_](#the-query-trick) to know when you get to see your changes.

## How do I add a new ``Organization``?

1. Go to the Research Software Directory's admin interface
([https://www.research-software.nl/admin/](https://www.research-software.nl/admin/)).
1. In the left pane, select ``Organization``
1. Check that the organization you want to add doesn't exist yet by searching via
the search box. If the search comes up empty, add the new organization by
clicking the blue ``+`` symbol.
1. Fill the name of the organization and provide a URL
1. Optionally, provide the organization's logo as an image.
   1. Make sure to use high quality logos, preferably with a transparent background.
   1. Organizations often offer the source files for their logos via a link on
      their website, usually named something like "Press" or "House style".
1. **Don't forget** to click ``Save`` when you're done.
1. Refer to section [_Delays_](#delays) and section [_The query
trick_](#the-query-trick) to know when you get to see your changes.

## How/when do I get to see my changes?

That depends. There are delays involved, and you'll need to use [_the query
trick_](#the-query-trick).

### Delays

Your changes/additions to the Admin interface, Zotero, GitHub or other places do
not show up immediately. The length of delay depends on the frequency at which
data is harvested from the external source. The schedule for the harvesting is
published at [https://research-software.nl/schedule](https://research-software.nl/schedule).

If your change only involves the Research Software Directory's Admin interface,
for example when you add a pre-existing ``Mention`` to a ``Software`` and press
``Save``, the corresponding database collection is updated immediately but the
data needed for populating the product page template is collected at intervals.
Refer to the schedule at
[https://research-software.nl/schedule](https://research-software.nl/schedule) and look for the
``python app.py resolve`` task.

If your change involves an external data source, e.g. you have new commits on
GitHub, a new release on Zenodo, or you've added an entry to Zotero, you will only see the resulting data show up
in the Research Software Directory after the corresponding data is harvested.
Check [https://research-software.nl/schedule](https://research-software.nl/schedule) and look for
the ``python app.py harvest commits`` and ``python app.py harvest citations``
tasks, respectively.

### The query trick

Second, note that both your browser and the Research Software Directory server
employ data caching for better performance. This means that when you ask to see
a product page, you will likely get to see an old version of that page, because
the browser knows it recently downloaded the data for that page, so it assumes
it doesn't need to download it again. You can easily circumvent this behavior
using the so-called _query_ part of URLs. Let's say you want to see
[https://www.research-software.nl/software/xenon](https://www.research-software.nl/software/xenon) but you want to make sure that
it really is the latest version, not some cached version. You can do so by
appending a question mark followed by some random characters, for example
[https://www.research-software.nl/software/xenon?somerandomstring](https://www.research-software.nl/software/xenon?somerandomstring)
(``?somerandomstring`` is the query part of the URL). This way, both the browser
and server consider it a previously unseen web page because the URL is
different, and will therefore show you the latest data. You'll need to come up
with a different value for ``somerandomstring`` every time you reload the data.
# Customize your instance of the Research Software Directory

Let's say you followed the steps above, and have a running instance of the
Research Software Directory. Now it is time to start customizing your Research
Software Directory. We have prepared some FAQs for customizations that are
common. For example, you can read up on the following topics:

1. [How do I change the colors?](faq/how-do-i-change-the-colors.md)
1. [How do I change the font?](faq/how-do-i-change-the-font.md)
1. [How do I change the logo?](faq/how-do-i-change-the-logo.md)
1. [How do I change when data collection scripts run?](faq/how-do-i-change-when-data-collection-scripts-run.md)
1. [How do I empty the database?](faq/how-do-i-empty-the-database.md)
1. [How do I make changes to the admin interface?](faq/how-do-i-make-changes-to-the-admin-interface.md)
1. [How do I add properties to the data schema?](faq/how-do-i-add-properties-to-the-data-schema.md)

It is suggested that you first do one or more of:

1. [How do I change the colors?](faq/how-do-i-change-the-colors.md)
1. [How do I change the font?](faq/how-do-i-change-the-font.md)
1. [How do I change the logo?](faq/how-do-i-change-the-logo.md)

Then, learn how to add properties to the schema:

1. [How do I add properties to the data schema?](faq/how-do-i-add-properties-to-the-data-schema.md)

Finally, learn how to empty the database, such that you can replace the sample
data with your own:

1. [How do I empty the database?](faq/how-do-i-empty-the-database.md)


## General workflow when making changes

After making your changes, here's how you get to see them:

1. Go to the terminal where you started `docker-compose`
1. Use Ctrl+C to stop the running instance of Research Software Directory
1. Check which docker containers you have with:

    ```shell
    docker-compose ps
    ```

    For example, mine says:

    ```shell
    docker-compose ps
           Name                     Command                State     Ports
    ----------------------------------------------------------------------
    rsd-admin            sh -c rm -rf /build/* && c ...   Exit 0
    rsd-authentication   /bin/sh -c gunicorn --prel ...   Exit 0
    rsd-backend          /bin/sh -c gunicorn --prel ...   Exit 0
    rsd-database         /mongo.sh --bind_ip 0.0.0.0      Exit 137
    rsd-frontend         /bin/sh -c sh -c "mkdir -p ...   Exit 0
    rsd-nginx-ssl        /bin/sh -c /start.sh             Exit 137
    rsd-reverse-proxy    /bin/sh -c nginx -g 'daemo ...   Exit 137
    rsd-harvesting       /bin/sh -c crond -d7 -f          Exit 137
    ```

    Use `docker-compose rm` to delete container by their **service name**, e.g. the `rsd-frontend` container:

    ```shell
    docker-compose rm frontend
    ```

    List all docker images on your system:

    ```shell
    docker images
    ```

    Note that image names consist of the environment variable `COMPOSE_PROJECT_NAME`, followed by `/`,
    followed by the service name. Remove as follows:

    ```shell
    docker rmi rsd/frontend
    ```
1. Make changes to the source code of the service whose container and image you just removed
1. Rebuild containers as necessary, using:

    ```shell
    docker-compose build frontend
    docker-compose up -d frontend
    ```
# Notes on security

The Research Software Directory is set up as a collection of services such as `backend`, `frontend`, `harvesting`,
etc. To avoid one service interfering with another, each service is dockerized. In a sense, docker is a bit like object
oriented programming: you have your data and methods together, and other methods don't have access to data unless you
specifically said that is OK.

Let's say that an attacker succeeds in somehow escaping the containment of the docker environment. If you set up your
instance on Amazon EC2/S3 as described in above, that may mean that they then have access to:

1. the Research Software Directory software
1. the collections in the Mongo database
1. the plaintext keys that are stored in `rsd-secrets.env`

Note that it does not mean they will have access to any of the rest of your institute's web site, since that content is
hosted on physically different machines, in a physically different location, with different networks, different
credentials, and probably a login procedure that is more challenging than just username/password.

With regard to (1), that information is public anyway. Just `git clone` would be a much easier way to get that
information. I guess the worst they could do here is make a change to the code and break the site. Or possibly, keep the
website like it is but use the Amazon machine to start mining bitcoins in the background. If that would happen though,
the usage graphs that Amazon provides would clearly show a change in behavior (from a spikey pattern directly related to
the crontab file to a uniform pattern).

With regard to (2), that information is (by default) harvested from public sources, so not much to be gained there. A
possible risk would be if the attacker aims to change the information displayed on the website, for example, pointing
links to the wrong place, or changing the description of a software package. Other risks might be that they empty the
database, or change data in such a way that the site no longer works.

With regard to (3), having access to some keys matters more than having access to others. Keys that don't matter so much
are `DOMAIN`, `SSL_ADMIN_EMAIL`, `SSL_DOMAINS`, `ZOTERO_LIBRARY`. Those are not really secret anyway, they are
more of a configuration value.

`AUTH_GITHUB_CLIENT_ID` and `AUTH_GITHUB_CLIENT_SECRET` are only useful if the attacker is a member of
`AUTH_GITHUB_ORGANIZATION`. If they in fact are a member, having access to the id and secret does not give them much
extra, because they already have direct access to the database service at that point (item 2 from the list above).

`GITHUB_ACCESS_TOKEN` provides readonly access to GitHub's publicly available data, it's just used to increase the
rate limit of GitHub's API.

`ZOTERO_API_KEY` provides readonly access to your Zotero library, which is probably public information anyway. Again
its main purpose is to increase the rate limit of allowed API usage.

If an attacker had access to `BACKUP_CMD`, that could potentially lead to the loss of your backups. They could use the
username and password to throw away any backups you have in that particular S3 bucket. (Note that you could make copies
to another bucket if you wanted to, or set up a different backup mechanism altogether; it _might_ be possible to
configure your S3 bucket such that you can write a backup file with the credentials from `BACKUP_CMD`, but not delete
them).

`JWT_SECRET` is only used to have services talk to each other, but doesn't give an attacker any abilities that they
would not already have, given that we assumed they have access to every service already.

A couple more remarks:

- Mongo has been in the news for mongo instances running on the internet without authentication (the default
  installation) leaking information. The Research Software Directory runs the Mongo instance in a private network
  wrapped by the token-protected backend service.

- Whoever is in charge of the Amazon machine needs to do the security updates of the host machine, in particular those
  updates that relate to the docker/docker-compose installation. Furthermore it's a good idea to also rebuild the docker
  images for each service, because then they get their respective updates.

- Also be aware that a service can have dependencies which may not be completely up to date, for example if the
  `requirements.txt` is outdated. This can have security implications.

- Regarding DDOS attacks, this is possible of course but not very likely in our opinion. However in such a case you
  would be charged more because there is more outbound traffic. You can mitigate it by setting a "budget exceeded" alarm
  on your usage.
# Using your own data sources

The research software directory is configured using a file with environment
variables called `.env`. An example config file
(`rsd-secrets.env.example`) is available, use it as a starting point.

```bash
cd research-software-directory
cp rsd-secrets.env.example .env
```

The config file has some placeholder values (`changeme`); they must be set by
editing the `rsd-secrets.env` file. Below are instructions on how to get the
different tokens and keys.

## `COMPOSE_PROJECT_NAME`

This is a prefix that docker-compose uses in naming its images, containers, and
volumes in order to avoid name clashes. Its default value is `rsd`.

## `AUTH_GITHUB_CLIENT_ID` and `AUTH_GITHUB_CLIENT_SECRET`

These environment variables are used for authenticating a user, such that they
can be granted access to the admin interface to create, read, update, and delete
items in the Research Software Directory.

These are the steps to assign values:

1. Go to [https://github.com/settings/developers](https://github.com/settings/developers)
1. Click the `New OAuth App` button
1. Under `Application name`, write something like _The Research Software
Directory's admin interface on localhost_
1. Under `Homepage URL` fill in some URL, for example, let it point to this
readme on GitHub. Despite the fact that it is a required field, its value is
not used as far as I can tell.
1. Optionally add a description. This is especially useful if you have multiple OAuth apps
1. The most important setting is the value for `Authorization callback url`.
Set it to [http://localhost/auth/get_jwt](http://localhost/auth/get_jwt) for now. We will revisit
`AUTH_GITHUB_CLIENT_ID` and `AUTH_GITHUB_CLIENT_SECRET` in the section about
deployment
[below](#make-your-instance-available-to-others-by-hosting-it-online-deployment)
1. Click `Register application`
1. Assign the `Client ID` as value for `AUTH_GITHUB_CLIENT_ID` and assign
the `Client Secret` as value for `AUTH_GITHUB_CLIENT_SECRET`

## `AUTH_GITHUB_ORGANIZATION`

Data is entered into the Research Software Directory via the admin interface.
Set `AUTH_GITHUB_ORGANIZATION` to the name of the GitHub organization whose
members should be allowed access to the admin interface. Most likely, it is the
name of the organization where you forked this repository to.

Note: members should make their membership of the GitHub organization public. Go
to
[https://github.com/orgs/&lt;your-github-organization&gt;/people](https://github.com/orgs/your-github-organization/people)
to see which users are a member of &lt;your-github-organization&gt;, and whether
their membership is public or not.

## `GITHUB_ACCESS_TOKEN`

To query GitHub's API programmatically, we need an access token. Here's how you can get one:

1. Go to [https://github.com/settings/tokens](https://github.com/settings/tokens)
1. Click `Generate new token`
1. Under `Token description`, fill in something like _Key to programmatically retrieve information from GitHub's API_
1. Verify that all scopes are unchecked
1. Use token as value for `GITHUB_ACCESS_TOKEN`

## `ZENODO_ACCESS_TOKEN`

To query Zenodo's API programmatically, we need an access token. Here's how you can get one:

1. Go to [https://zenodo.org/account/settings/applications/tokens/new/](https://zenodo.org/account/settings/applications/tokens/new/)
1. For name, fill in something like _Key to retrieve data from Zenodo_
1. Make sure all scopes are unselected
1. Click Create
1. Fill in the long string you get as value for `ZENODO_ACCESS_TOKEN`

## `ZOTERO_LIBRARY`

When getting the references data from Zotero, this environment variable
determines which library on Zotero is going to be harvested. Go to
[https://www.zotero.org/groups/](https://www.zotero.org/groups/) to see which Zotero groups you are a member of.
If you click on the `Group library` link there, the URL will change to
something like
[https://www.zotero.org/groups/1689348/netherlands_escience_center/items](https://www.zotero.org/groups/1689348/netherlands_escience_center/items), where
`1689348` is the value you need to assign to `ZOTERO_LIBRARY`.


## `ZOTERO_API_KEY`

To query Zotero's API programmatically, we need an API key. Here's how
you can get one:

1. [https://www.zotero.org/settings/keys](https://www.zotero.org/settings/keys)
1. Click `Create new private key`
1. Type a description of the key, e.g. _API key to access library X on Zotero_
1. Under `Personal library`, make sure only `Allow library access` is checked.
1. Under `Default group permissions`, choose `None`
1. Under `Specific groups`, check `Per group permissions`
1. Set `Read only` for the group that you want to harvest your references data from; verify that any other groups are set to `None`
1. Click the `Save Key` button at the bottom of the page.
1. On the `Key Created` page, you will see a string of random character,
something like `bhCJSBCcjzptBvd3fvliYOoE`. This is the key; assign it to
`ZOTERO_API_KEY`

## `BACKUP_CMD`

This environment variable is used for making a daily backup of the database with
software, people, projects, etc. As it is typically only used during deployment,
leave its value like it is for now; we will revisit it in the page about
[deployment](hosting.md).


## `JWT_SECRET`

<!-- This environment variable is used for ... TODO -->

The `JWT_SECRET` is simply a string of random characters. You can generate one
yourself using the `openssl` command line tool, as follows:

```bash
openssl rand -base64 32
```

Assign the result to `JWT_SECRET`.

## `DOMAIN` and `SSL_DOMAINS`

These environment variables are not relevant when you're running your instance
locally. Leave their values like they are in `rsd-secrets.env.example` for the
time being. We will revisit them in the page about [deployment](hosting.md).
# Creating a Zotero account and installing the Zotero client

1. Create an account with Zotero at
[zotero.org/user/register](https://www.zotero.org/user/register/). Note:
Lastpass will **not** popup and ask to save the password, so **save it
beforehand**!
1. Login.
1. Go to the Netherlands eScience Center group on
[zotero.org/groups](https://www.zotero.org/groups/1689348/netherlands_escience_center).
1. Click `Join`. Be aware that a group administrator needs to grant you access
to the group so it may take a while, depending on people's availability.
1. From here on out, use the offline client. Download it from
[zotero.org](https://www.zotero.org/download/).
1. Install Zotero.

## Configure Zotero client

Once the Zotero application is running:

1. Go to `Edit`/`Preferences`/`Sync`
1. Enter username & password -> OK
1. Click the `sync` button:

    ![sync-button](/docs/images/zotero-sync-button.png)

1. Netherlands eScience Center should appear under `Group Libraries`.

## General workflow for adding items

Make sure the item you want to add has some sort of identifier such as a DOI or
a URL; without it, whatever you're adding is just hearsay.

1. Synchronize with the Zotero server.
1. On the left-hand side, select the ``Miscellaneous`` folder (or, you could try to find a more appropriate folder, but
   know that that won't matter for anything related to the Research Software Directory -- it does not use the folder
   structure from Zotero, just the items that are in it).
1. Click the `Add by identifier` button.

    ![add-by-identifier](/docs/images/zotero-add-by-identifier.png)

    Fill in your item's DOI, e.g. ``10.5281/zenodo.1299523``

1. Wait, you're not done yet! Verify that the metadata in the `Info` tab is
correct, as follows:
   1. In the right-hand panel, select the `Info` tab
   1. Verify that the `Item Type` (top of the list) is correct.
   1. For some `Item Type`s, you need to fill in `Type` as well (see the list below).
   1. Verify that any names have been entered correctly, for instance in the `Author` or
   `Programmer` fields. A person's first name and last name each have their
   own input field. You can switch between single and two field entry by
   pressing the small button next to the names. Name particles (`de`, `van der`, etc.)
   should be included in the lastname so for `Jan de Groot` use lastname `de
   Groot`.
   1. Verify that the dates have been entered correctly. Use a single date (don't
   use date ranges). Because dates are tricky, Zotero shows a string such as
   `y m d` or `d m y` next to each date, to show how it has interpreted each number.
1. If everything looks good, synchronize with the Zotero server again.


## `Item type`-specific information

Below is a list of output we would like to keep track of, with a short
description. Pick the one that best describes your output, and fill out the
metadata required. By default Zotero shows a much larger list of metadata,
use your judgment to fill out the other fields as necessary.

If you have an item that doesn't fit, please [open an issue
here](https://github.com/research-software-directory/research-software-directory/issues) and we'll figure it
out and update this document.

### Software

You don't need to add `Software` items to Zotero. We keep track of our software
output via the Research Software Directory. If you have a software package that
you want to add, use the Research Software Directory's Admin interface as
explained [here](entering-data.md).

### Papers

If you have a DOI, use the 'Add by identifier' button, and check that the ``Item type`` in the ``Info`` tab is correct.

### Datasets

If you don't have a DOI for the dataset yet, make sure to [upload a copy of the item
to Zenodo](https://zenodo.org/deposit/new), FigShare or an other place that provides a DOI.
Set `Item Type` to `Web page`, and set the item's URL to [https://doi.org/&lt;yourdoi&gt;](https://doi.org/<yourdoi>).

### Conference poster or presentation slides

Use `Item Type` `Presentation`, set the field `Type` to either `Conference
Poster` or `Conference Presentation`. Note that you need to [upload the poster or
slides to Zenodo](https://zenodo.org/deposit/new), FigShare, or some other place that provides a DOI. Since
Zotero does not have a DOI field for `Presentation`s, use the URL field for this
purpose (e.g. [https://doi.org/&lt;yourdoi&gt;](https://doi.org/<yourdoi>)).

### Workshop, lecture, or demonstration

Use `Item Type` `Presentation`, set the field `Type` to `Workshop`, `Lecture`,
or `Demonstration`.

### Report

For items that have not been peer reviewed, such as internal reports, white
papers, etc., use `Item Type` `Report`. Make sure to [upload a copy of the item
to Zenodo](https://zenodo.org/deposit/new), FigShare or an other place that provides a DOI. The item's URL should
point to [https://doi.org/&lt;yourdoi&gt;](https://doi.org/<yourdoi>).

### Thesis

A PhD, Master, or Bachelor thesis. Set `Item Type` to `Thesis` and fill in
`Bachelor`, `Master`, `PhD` in the `Type` field.

### Other types

Please choose the most appropriate type from `Blogpost`, `Book`, `Book Section`,
`Interview`, `Magazine Article`, `Newspaper Article`, `Podcast`, `Radio
Broadcast`, `TV Broadcast`, `Video Recording`, `Webpage`. And make a best effort
attempt at filling out the other metadata on the `Info` tab.
# Contents

<!-- see also bottom of main readme -->

1. [Entering data about your software in an existing instance](entering-data.md)
1. [Configuring your instance to use your own data sources](configuring-your-instance-to-use-your-own-data-sources.md)
1. [Changing the look and feel](changing-the-look-and-feel.md)
1. [Hosting your instance online](hosting.md)
1. [Running an instance of the Research Software Directory in production](production.md)
1. [Finding your way: Research Software Directory services](services-overview.md)
1. [Documentation for developers](documentation-for-developers.md)
1. [Documentation for maintainers](documentation-for-maintainers.md)
1. [Security concerns](security.md)
1. [Contributing](../.github/CONTRIBUTING.md)
# Overview of services of the Research Software Directory

The Research Software Directory is made of the following services

- [frontend](/frontend): Python web application which renders HTML pages for normal visitors
- [backend](/backend): Python web service for programmatic access to the directory data (software, projects, persons, organizations). Used by other services to fetch and set data.
- [reverse-proxy](/reverse-proxy): Web server responsible for combining all web based services behind a single domain and port. Also hosts the static files of other services for best performance and caching.
- [admin](/admin/): React application for editing the data in the directory. Hosted by `reverse-proxy` service.
- [auth-github](/auth-github): Protects the `admin` service by forcing authentication with a GitHub account and authorization using GitHub organization membership.
- [https](/https): Responsible for encrypting (HTTPS) traffic from `reverse-proxy` service.- [backup](/backup): For backup, copies database to an S3 bucket every day. Only runs when configured.
- [database](/database): A Mongo database. Used by `backend` service to store data. Initializes with sample data when database is empty.
- [graphs](/graphs): Web page which shows metrics of directory. Hosted by `reverse-proxy` service.
- [harvesting](/harvesting): Scheduled jobs which periodicaly fetch external data. For example commits from GitHub and mentions from Zotero.

All these services are started in the [Docker compose file](/docker-compose.yml). All services have their own directory in the repository.
# Documentation for developers

## Try it out, step 1/3: Fork and clone

Click the ``Fork`` button on
[https://github.com/research-software-directory/research-software-directory/](https://github.com/research-software-directory/research-software-directory/)
to fork to your own GitHub organization or GitHub profile, then:

```shell
git clone https://github.com/<your-github-organization>/research-software-directory.git
```

## Try it out, step 2/3: Configure

See section [Configuring your instance to use your own data sources](./configuring-your-instance-to-use-your-own-data-sources.md).

## Try it out, step 3/3: Start the complete stack using [docker-compose](https://docs.docker.com/compose/)

```shell
# build all containers:
docker-compose build

# start the full stack using docker-compose:
docker-compose up -d

# see logging from all services with
docker-compose logs --follow

# or from a specific service only, e.g. backend
docker-compose logs --follow backend
```

After the Research Software Directory instance is up and running, we want to
start harvesting data from external sources such as GitHub, Zotero, Zenodo, etc.
To do so, open a new terminal and run

```shell
docker-compose exec harvesting python app.py harvest all
```

You should see some feedback in the newly opened terminal.

After the ``harvest all`` task finishes, several database collections should
have been updated, but we still need to use the data from those separate
collections and combine them into one document that we can feed to the frontend.
This is done with the ``resolve all`` task, as follows:

```shell
docker-compose exec harvesting python app.py resolve all
```

By default, the ``resolve all`` task runs every 10 minutes anyway, so you could just wait for a bit, until you see some output scroll by that is generated by the ``rsd-harvesting`` container, something like:

```shell
rsd-harvesting     | 2018-07-11 10:30:02,990 cache [INFO] processing software Xenon command line interface
rsd-harvesting     | 2018-07-11 10:30:03,013 cache [INFO] processing software Xenon gRPC server
rsd-harvesting     | 2018-07-11 10:30:03,036 cache [INFO] processing software xtas
rsd-harvesting     | 2018-07-11 10:30:03,059 cache [INFO] processing software boatswain
rsd-harvesting     | 2018-07-11 10:30:03,080 cache [INFO] processing software Research Software Directory
rsd-harvesting     | 2018-07-11 10:30:03,122 cache [INFO] processing software cffconvert
rsd-harvesting     | 2018-07-11 10:30:03,149 cache [INFO] processing software sv-callers
```

## Verifying the local installation

Open a web browser to verify that everything works as it should. Below are some things to check:

### Frontend

- [``http://localhost``](http://localhost) should show the software index page to the local instance of the Research Software Directory
- [``http://localhost/projects``](http://localhost/projects) should show the project index page to the local instance of the Research Software Directory
- [``http://localhost/projects/``](http://localhost/projects/)  should show the project index page to the local instance of the Research Software Directory
- [``http://localhost/projects/764``](http://localhost/projects/764) should show a project page (here: ABC-MUSE) in the local instance of the Research Software Directory
- [``http://localhost/software``](http://localhost/software) should show the software index page to the local instance of the Research Software Directory
- [``http://localhost/software/``](http://localhost/software/) should show the software index page to the local instance of the Research Software Directory
- [``http://localhost/software/xenon``](http://localhost/software/xenon) should show a product page (here: Xenon) in the local instance of the Research Software Directory
- [``http://localhost/graphs``](http://localhost/graphs) should show you some integrated statistics of all the packages in the local instance of the Research Software Directory
- [``http://localhost/about``](http://localhost/about) should show the about page in the local instance of the Research Software Directory

### Admin interface

- [``http://localhost/admin``](http://localhost/admin) should show the Admin interface to the local instance of the Research Software Directory

### API

- [``http://localhost/api/mention``](http://localhost/api/mention) should show a JSON representation of all mentions in the local instance of the Research Software Directory
- [``http://localhost/api/organization``](http://localhost/api/organization) should show a JSON representation of all organizations in the local instance of the Research Software Directory
- [``http://localhost/api/project_cache``](http://localhost/api/project_cache) should show a JSON representation of all projects in the local instance of the Research Software Directory, with all references resolved
- [``http://localhost/api/project``](http://localhost/api/project) should show a JSON representation of all projects in the local instance of the Research Software Directory
- [``http://localhost/api/release``](http://localhost/api/release) should show a JSON representation of all releases in the local instance of the Research Software Directory
- [``http://localhost/api/schema``](http://localhost/api/schema) should show the schema for the local instance of the Research Software Directory
- [``http://localhost/api/software_cache``](http://localhost/api/software_cache) should show a JSON representation of all software in the local instance of the Research Software Directory, with all references resolved
- [``http://localhost/api/software/xenon``](http://localhost/api/software/xenon) should show a JSON representation of a product (here: Xenon) in the local instance of the Research Software Directory
- [``http://localhost/api/software``](http://localhost/api/software) should show a JSON representation of all software in the local instance of the Research Software Directory

The api endpoints also support the following query parameters:

- ``sort`` (e.g. ``sort=updatedAt``)
- ``direction`` (e.g. ``direction=desc``)
- ``limit`` (e.g. ``limit=3``)

Which can be combined in the usual way, e.g.

- [``http://localhost/api/mention?limit=3&direction=desc&sort=updatedAt``](http://localhost/api/mention?limit=3&direction=desc&sort=updatedAt) should return the 3 mentions that were updated most recently.

### Citation

- [``http://localhost/cite/xenon?version=3.0.4&format=bibtex``](http://localhost/cite/xenon?version=3.0.4&format=bibtex) should return a reference manager file for software package Xenon version 3.0.4 in BibTeX format.

### OAI-PMH

- [``http://localhost/oai-pmh?verb=ListRecords&metadataPrefix=datacite4``](http://localhost/oai-pmh?verb=ListRecords&metadataPrefix=datacite4) should return an XML document with metadata about all the packages that are in the local instance of the Research Software Directory, in DataCite 4 format.

### Harvesting schedule

- [``http://localhost/schedule``](http://localhost/schedule) should return the cron job describing when each harvester is scheduled to run.

## Removing local state

The Research Software Directory stores its state in a couple of places. While
doing development, sometimes you need to clear the local state, therefore this
section lists some ways to clear such state. Be aware that running these
commands results in the **LOSS OF DATA**.

- Remove a docker container:

    ```shell
    # remove a container associated with a specific service from docker-compose.yml
    docker-compose rm <service name>

    # remove any container corresponding to any of the services defined in docker-compose.yml
    docker-compose rm
    ```

- Remove a docker image:

    ```shell
    # remove a specific image
    docker rmi <image name>
    ```

- Docker bind mounts store data in ``<project directory>/docker-volumes``, remove with:

    ```shell
    sudo rm -rf docker-volumes/db/ docker-volumes/oaipmh-cache/
    ```

- Docker static volumes store data. Refer to [/docker-compose.yml](/docker-compose.yml) to see which services use which volumes. Remove a volume with:

    ```shell
    # remove volumes that are not in use by any containers
    docker volume prune

    # or remove a specific volume
    docker volume rm <volume>
    ```

- Docker networks. By default, the services in [/docker-compose.yml](/docker-compose.yml) share a network named ``rsd_default``. Remove a network with

    ```shell
    # remove networks that are not in use by any containers
    docker network prune

    # or remove a specific network
    docker network rm <network>
    ```

- To remove Docker container, images, static volumes and networks in single step (bind mounts still need to be removed separately) with

    ```shell
    docker-compose down --rmi all -v
    ```

## Checking if there's any documentation with invalid links

The repository comes with documentation spread out over multiple MarkDown files.
[This workflow file](./../.github/workflows/markdown-link-checker.yml) is set up to check whether there are broken links in
any of them. If you want to check this locally, you can do so with:

```shell
# from the repository root directory
npm install
npm run mlc
```

## Running the superlinter locally

We use GitHub's [Super-Linter](https://github.com/github/super-linter) to lint all directories using a variety of
linters. You can run the superlinter using `docker`, as follows:

```shell
# get the linter
docker pull github/super-linter:latest

# run the linter
docker run \
   -e RUN_LOCAL=true \
   -v ${PWD}:/tmp/lint \
   github/super-linter
```

The superlinter can be a bit slow if you run all checks on all directories, but you can run just one check on one file
with a [specific linter](https://github.com/github/super-linter#environment-variables) if needed, like so:

```shell
# run the linter with only pylint check enabled
docker run \
   -e RUN_LOCAL=true \
   -e VALIDATE_PYTHON_PYLINT=true \
   -v ${PWD}/harvesting/app.py:/tmp/lint/app.py \
   github/super-linter
```

or evaluate a whole subdirectory, all checks:

```shell
cd harvesting
docker run \
   -e RUN_LOCAL=true \
   -v ${PWD}:/tmp/lint \
   github/super-linter
```

or evaluate the same whole subdirectory, but do just one check:

```shell
cd harvesting
docker run \
   -e RUN_LOCAL=true \
   -e VALIDATE_PYTHON_PYLINT=true \
   -v ${PWD}:/tmp/lint \
   github/super-linter
```

By default, the GitHub Super-Linter generates a log file inside the container, which is subsequently mapped to your file
system where it appears as a file named `super-linter.log` with root permissions. You can disable this behavior by
setting the bind mount as read-only (`:ro`), as follows:

```shell
docker run \
   -e RUN_LOCAL=true \
   -v ${PWD}:/tmp/lint:ro \
   github/super-linter
```

### The output on GitHub looks different

The workflow file that we use for our continuous integration on GitHub Actions has a per-directory configuration with a
custom list of linters for each directory (see [/.github/workflows/linting.yml](/.github/workflows/linting.yml)). If you
want your local setup to reflect exactly what runs on GitHub, it may therefore be convenient to use
[`act`](https://github.com/nektos/act).

See also [this
comment](https://github.com/research-software-directory/research-software-directory/pull/624#pullrequestreview-528215446).

## Visualizing ``docker-compose.yml``

It is sometimes helpful to visualize the structure in the ``docker-compose.yml`` file.
Use [https://github.com/pmsipilot/docker-compose-viz](https://github.com/pmsipilot/docker-compose-viz) to generate a png image.

```shell
docker run --rm -it --name dcv -v $(pwd):/input pmsipilot/docker-compose-viz \
   render -m image --output-file=docs/images/docker-compose.png docker-compose.yml
```

For example,

![images/docker-compose.png](images/docker-compose.png)
# Running an instance of the Research Software Directory in production

## Temporarily disabling the admin interface

You can stop users from making changes to the database by disabling the authentication to the admin interface, as follows:

```shell
docker-compose stop auth
```

This can be useful when upgrading the data or software to a newer version.

Note that even with a stopped `auth` service, determined users can still access `backend` directly, and users who were
logged in before you disable `auth` will still be able to use the admin interface to make changes.

It it recommended that you post a message that users see when they try to use the admin interface. The easiest way to do
this is as follows:

```shell
# log in to admin service
docker-compose run --workdir="/build" admin /bin/sh

# rename index.html
mv index.html index.html.disabled

# put a message in index.html
echo "<html>" \
"<head></head>" \
"<body>Sorry, we're doing maintenance right now. " \
"Hopefully be back soon.</body></html>" > index.html
```

You can check if it works by using your browser to navigate to the admin interface. Instead of the normal interface, you
should now see your message.

To enable admin interface again, do the previous instructions in reverse:

```shell
# inside admin container
mv /build/index.html.disabled /build/index.html

# Enable logins again
docker-compose start auth
```


## Updating a production instance

Every now and then, the production instance needs to be updated, so the server
can get the latest security patches, and the Research Software Directory
software itself can be updated to include the latest features.

The steps below differentiate between the old and the new instance of the Research
Software Directory; the old instance has IP ``35.156.38.208``, the new one has
IP ``3.122.233.225``. Your IP addresses will likely be different.

1. Make a new Amazon instance by following the notes above. Some things to think about:
    - Reuse the existing security group.
    - Reuse the existing key pair.
    - Verify that you're allowed to ssh into the new instance.
1. Transfer the ``rsd-secrets.env`` file from the old instance to the new instance.

    ```shell
    cd $(mktemp -d)
    scp -i ~/.ssh/rsd-instance-for-nlesc-on-aws.pem \
       ubuntu@35.156.38.208:/home/ubuntu/rsd/rsd-secrets.env .
    scp -i ~/.ssh/rsd-instance-for-nlesc-on-aws.pem \
       ./rsd-secrets.env \
       ubuntu@3.122.233.225:/home/ubuntu/rsd/rsd-secrets.env
    ```
1. On the remote, create the symlink `.env` and let it point to `rsd-secrets.env`:

    ```shell
    cd ~/rsd
    ln -s rsd-secrets.env .env
    ```

1. Stop new additions to the database in the old Research Software Directory instance by following the notes from
   [_Temporarily disabling the admin interface_](#temporarily-disabling-the-admin-interface). This will effectively
   disable the ``rsd-admin`` service.

1. Create the backup files in the old Research Software Directory instance:

    ```shell
    # start an interactive shell in the backup container
    docker-compose exec backup /bin/sh

    # create the backup files in the container's /dump directory
    /app # mongodump \
      --host ${DATABASE_HOST} \
      --port ${DATABASE_PORT} \
      --db ${DATABASE_NAME} \
      --out /dump

    # leave the backup container
    exit

    # Copy the dump directory out of the docker container
    docker cp $(docker-compose ps -q backup):/dump/rsd /home/ubuntu/rsd/dump
    ```

1. Transfer the dumped json and bson files from the old to the new instance

    ```shell
    scp -r -i ~/.ssh/rsd-instance-for-nlesc-on-aws.pem \
       ubuntu@35.156.38.208:/home/ubuntu/rsd/dump .

    scp -r -i ~/.ssh/rsd-instance-for-nlesc-on-aws.pem \
       ./dump/* ubuntu@3.122.233.225:/home/ubuntu/rsd/database/db-init/

    ```

1. Start the new Research Software Directory instance.

    ```shell
    ssh -i ~/.ssh/rsd-instance-for-nlesc-on-aws.pem ubuntu@3.122.233.225
    cd /home/ubuntu/rsd

    docker-compose build
    docker-compose up -d
    ```

1. Check [/CHANGELOG.md](/CHANGELOG.md) to see if you need to run any command to
   migrate data, e.g. when a collection has changed its schema.

1. Next, harvest all the data from external sources using:

    ```shell
    docker-compose exec harvesting python app.py harvest all
    docker-compose exec harvesting python app.py resolve all
    ```

1. In case the old instance had problems with harvesting of the mentions, you
   may need to retrieve all mentions, as follows:

    ```shell
    docker-compose exec harvesting python app.py harvest mentions --since-version 0
    ```

1. Check if the instance works correctly using a browser to navigate to
   the new instance's IP address.
1. If everything looks good, stop the Research Software Directory in the old instance

    ```shell
    docker-compose stop
    ```

1. Disassociate the ElasticIP address from the old instance.
1. Associate the ElasticIP address with the new instance.

As a final step, use the Amazon EC2 management console to ``Stop`` (not
``Terminate``) the old instance. This way, the old instance can still be
reactivated in case you need to get back to the old version.
# Documentation for repository maintainers

## Making a release

1. Write the release notes
1. Update CITATION.cff
1. Generate the metadata file for Zenodo using [cffconvert](https://pypi.org/project/cffconvert/).

    ```shell
    pip install --user cffconvert
    cffconvert --outputformat zenodo --ignore-suspect-keys --outfile .zenodo.json
    ```

    ```shell
    # git add, commit, and push everything
    ```

1. Make sure that everything is pushed

    ```shell
    cd $(mktemp -d)
    git clone https://github.com/research-software-directory/research-software-directory.git
    cd research-software-directory
    ```

1. Follow the notes from the ['For developers'](#documentation-for-developers) section above, and verify that it all works as it should.
1. Use GitHub's ``Draft a new release`` button [here](https://github.com/research-software-directory/research-software-directory/releases) to make a release.

## Pulling in changes from upstream using a three-way merge

Set ``UPSTREAM`` and ``DOWNSTREAM`` to the different sources you want to
three-way merge between, e.g.

```shell
UPSTREAM=https://github.com/research-software-directory/research-software-directory.git
DOWNSTREAM=https://github.com/process-project/research-software-directory.git
```

Then:

```shell
cd $(mktemp -d)
mkdir left middle right
cd left && git clone $UPSTREAM . && cd -
cd middle && git clone $DOWNSTREAM . && git branch develop && git checkout develop && cd -
cd right && git clone $DOWNSTREAM . && cd -
meld left middle right &
```

You should only make changes to the ``middle`` one. When you're done making your changes,

```shell
git add <the files>
git commit
git push origin develop
```
# How do I change the logo?

## Relevant files

- `frontend/templates/layout_template.html`
- `frontend/templates/logo.html`
- Section `#header_logo` in `frontend/style/components/_header.scss`


By default the site uses some SVG trickery combined with Jinja templating to
show an SVG logo on the frontend. The logo's appearance is controlled from the
CSS, which in turn is generated from SCSS source files. While this offers full
control over the logo's appearance, it is a little bit tricky to set up. Details
are described in [The vector image way](#the-vector-image-way) below.

However, if your institute already has a regular website that includes a logo in
a raster format like PNG or JPG, there is a much easier way to do things. See
[The raster image way](#the-raster-image-way) below.

## The vector image way

This part is a stub.

## The raster image way


In `layout_template.html`, find the `div` with `id=header_logo`:

```html+jinja
<div id="header_logo">
    <a href="/">
        {% with position="header" %}
        {% include 'logo.html' %}
        {% endwith %}

        {% block rsd_home %}
        <div id="header_text">
            research software directory
        </div>
        {% endblock %}
    </a>
</div>
```

and replace this part:

```html+jinja
{% with position="header" %}
{% include 'logo.html' %}
{% endwith %}
```

with the link to the raster image as used on your institute website, e.g.:

```html+jinja
<img src="https://www.esciencecenter.nl/img/cdn/logo_escience_center.png" />
```

You may need to either use a `style` tag (not preferred)

```html+jinja
<img src="https://www.esciencecenter.nl/img/cdn/logo_escience_center.png" style="some style here" />
```

or include some styling in `_header.scss` in order to set some of the `img`'s style properties such as its height.

Note that for the latter option, you'll need a program that can compile SCSS
files into CSS. For example, on Ubuntu you can use `pysassc`:

```shell
sudo apt install pysassc
```

From the `frontend` root directory, run:

```shell
pysassc --style=compressed --sourcemap style/rsd.scss static/style/rsd.scss.css
```

Afterwards, your new `img` style should be included in `static/style/rsd.scss.css`.

Refer to the [general workflow when making changes](/README.md#general-workflow-when-making-changes) to update the Docker container with the new content.
# How do I make changes to the admin interface?

The admin interface is derived automatically from the schema, so read
[How do I add properties to the data schema?](/docs/faq/how-do-i-add-properties-to-the-data-schema.md)
to learn how to make changes to the admin interface.
# How do I add properties to the data schema?

## Relevant files

- `backend/schemas/*.json`
- `admin/public/settings.json`

The data schema is the document that describes what the data in the database
looks like. You can find those documents here: `backend/schemas/`. As you can
see, there are multiple documents, one for each collection in the MongoDB
database.

We can make different kinds of changes to the schema, for example, we may change
a property, remove a property or add a property. Some of these are more complex
than others. Because the simplest change is making an addition to the schema,
that's where we'll start.

Let's say we want to add a key `grants` to the `software` schema, such that
we can record the grant numbers of the project in which certain software was
developed.

Insert `grants` as a property of the existing top-level property `properties`:

```json
"grants": {
    "type": "array",
    "items": {
        "type": "string"
    }
},
```

This extends the software schema with an array of strings. Additionally, we need
to find the corresponding part in `admin/public/settings.json`, such that the
admin interface reflects our adding of `grants`. To do so, first find the
`software` property under `resources`, then look for its property `properties`.

Add `grants`, as follows:

```json
"grants": {
},
```

This minimal addition will make the `grants` property show up in the admin
interface. However, it's probably a good idea to provide a small description to
whoever is going to fill in the grants value in the admin interface. For this,
you can add a `label`, like so:

```json
"grants": {
    "label": "Specify which research grants were used to develop this software."
},
```

Finally, you can also control the order in which properties appear using
`sortIndex`; if you don't define a `sortIndex`, it will be added at the
bottom of the form. Given that the `sortIndex` of most other properties are
around 100 or even higher, a value of 50 will place it near the top of the form:

```json
"grants": {
    "label": "Specify which research grants were used to develop this software.",
    "sortIndex": 50
},
```

Refer to the [general workflow when making
changes](/README.md#general-workflow-when-making-changes) to update the Docker
container with the new content.

## Migrating pre-existing data to the new schema

If you go to [http://localhost/admin](http://localhost/admin), then select Software from the menu on the left and click
the blue plus sign, there should be a section in the form labeled `grants`, with the description that you provided.
However, pre-existing software documents do not have any grant numbers yet, so if you try to load e.g.
[http://localhost/admin/software/xenon](http://localhost/admin/software/xenon) the admin interface will complain about
missing data.

To avoid those errors, we need to update the pre-existing data from the old
schema to the new schema. We will use the MongoDB console for this.

Make sure the Research Software Directory is still up and running with

```shell
docker ps -a
```

Start a new terminal and go to the project's root directory.

Start a shell in the `database` container, as follows:

```shell
docker-compose exec database /bin/sh
```

Once in the shell, run `mongo rsd` to gain access to the `rsd` database via
the MongoDB shell.

In the MongoDB shell, run (see explanation below):

```shell
db.software.update({}, {$set: {"grants": []}}, {"multi": true})
```

`db.software.update` takes three arguments here. The first `{}` selects all
documents from the collection; the second `{$set: {"grants": []}}` specifies
that the `grants` property should be set to an empty array; the third argument
`{"multi": true}` specifies that if there are multiple documents in the
selection, that they should all be updated.

For more details on how to use the MongoDB shell, please refer to the documentation:
[https://docs.mongodb.com/manual/reference/method/db.collection.update/](https://docs.mongodb.com/manual/reference/method/db.collection.update/).

Go to any software package in the admin interface and verify that the errors
which were there previously have now gone.

# How do I change the font?

Which font you see on the frontend is governed by the style file `frontend/static/style/rsd.scss.css`. This file is generated from SCSS source
files. In a nutshell, these are the steps to change to a different font (more
details below):

1. Add static assets.
1. Update the SCSS variables that point to the fonts.
1. Recompile the CSS from its SCSS source files.


## Add static assets

Relevant files:

- `frontend/static/style/fonts/<new directory>/YourFontLight.ttf`
- `frontend/static/style/fonts/<new directory>/YourFontRegular.ttf`
- `frontend/static/style/fonts/<new directory>/YourFontBold.ttf`

Go to [https://fonts.google.com/](https://fonts.google.com/) to select your
favorite font, and download the TrueType font files (`*.ttf`) via the popup
menu. Font files are part of a website's so-called _static assets_.

Make a new directory in `frontend/static/style/fonts`, e.g. I downloaded
the RobotoCondensed font, so I made `frontend/static/style/fonts/roboto-condensed`.
I then copied `RobotoCondensed-Regular.ttf`, `RobotoCondensed-Light.ttf`, and `RobotoCondensed-Bold.ttf` into the new directory.

## Update the SCSS variables

Relevant files:

- `frontend/style/base/_fonts.scss`
- `frontend/style/settings/_variables.scss`

In order to use the new static assets, we need to update some SCSS variables. First,
we need to make an `@font` entry in `frontend/style/base/_fonts.scss` for every font we want to use. Currently, the `_fonts.scss` file uses the Akkurat font, but we can't use that due to its restricted usage rights. So, replace

```scss
@font-face {
    font-family: "Akkurat-bold";
    src: url('fonts/akkurat/lineto-akkurat-bold.eot');
    src: url('fonts/akkurat/lineto-akkurat-bold.eot?#iefix') format('embedded-opentype'),
         url('fonts/akkurat/lineto-akkurat-bold.woff') format('woff'),
         url('fonts/akkurat/lineto-akkurat-bold.ttf') format('truetype');
    font-weight: normal;
    font-style: normal;
}
```

with


```scss
@font-face {
    font-family: "RobotoCondensedBold";
    src: url('fonts/roboto-condensed/RobotoCondensed-Bold.tff') format('truetype');
    font-weight: normal;
    font-style: normal;
}
```

Note: there are different font file formats. Since we're not using any `eot` or `woff` files, we only need one url per font.

Also update the Regular and Light sections in `_fonts.scss`.

Now, open the file `frontend/style/settings/_variables.scss`, find the
defintion of the `$primaryFont`, `$primaryFontBold`, and
`$primaryFontLight` variables, and for all three, replace the 'Akkurat' string
with the corresponding name you used for `font-family` in `_fonts.scss`.

My font section in `_variables.scss` now looks like this:

```scss
$primaryFont:  'RobotoCondensedRegular', Helvetica, arial, sans-serif;
$primaryFontBold:  'RobotoCondensedBold', Helvetica, arial, sans-serif;
$primaryFontLight:  'RobotoCondensedLight', Helvetica, arial, sans-serif;
```

## Recompile the CSS

Relevant files:

- `frontend/static/style/rsd.scss.css`
- `frontend/static/style/rsd.scss.css.map`

For this step, you'll need a program that can compile SCSS files into CSS. For
example, on Ubuntu you can use `pysassc`:

```shell
sudo apt install pysassc
```

From the `frontend` root directory, run:

```shell
pysassc --style=compressed --sourcemap style/rsd.scss static/style/rsd.scss.css
```

Afterwards, your new font should be included in `static/style/rsd.scss.css`.

Refer to the [general workflow when making changes](../changing-the-look-and-feel.md#general-workflow-when-making-changes) to update
the Docker container with the new content.
# How do I empty the database?

Obviously, this part of the documentation can lead to the **LOSS OF DATA**. Make
sure you have copies of all data that you care about.

Assuming that `docker images` shows the `rsd/` images, and `docker ps -a`
shows the `rsd-` docker containers, add the environment variables to the
terminal:

```shell
docker-compose up -d && docker-compose logs --follow
```

In a new terminal,

```shell
docker-compose exec database /bin/sh
```

Run the `mongo` command inside the `database` container to start the Mongo
shell there.

```shell
mongo
```

In Mongo shell, tell Mongo you want to use the `rsd` database:

```shell
use rsd
```

Ask for the list of collections that Mongo knows about:

```shell
show collections
```

For every collection that you want to delete, e.g. `commit` and `project`:

```shell
db.commit.deleteMany({})
db.project.deleteMany({})
```

For reference, here is the Link to the Mongo shell documentation:
[https://docs.mongodb.com/manual/reference/method/db.collection.deleteMany/#db.collection.deleteMany](https://docs.mongodb.com/manual/reference/method/db.collection.deleteMany/#db.collection.deleteMany)

Type Ctrl-D to exit the Mongo shell.

After you are done making changes to the collections, you will want to update
the data in `/database/db-init`. The data in this directory is part of the
GitHub repo and serves as sample data for when people do a `git clone`. Now
that you have emptied some collections, that initial data needs to be updated,
as follows:

Dump the contents of the `rsd` database to a directory by running (still from
within the `database` container):

```shell
mongodump --db rsd
```

The dump files should be located at `/dump/rsd/` (inside the `database`
container). Verify that the files are there and then leave the `database`
container with:

```shell
exit
```

You should now be back in the original terminal. From there, copy the database dump files from inside the
container to outside the container:

```shell
docker cp rsd-database:/dump/rsd/ database/db-init/
```

Move the data to the appropriate place (`./database/db-init`) and delete the `rsd` directory:

```shell
cd database/db-init
mv rsd/* .
rm -r rsd
cd ../..
```

Update your git repository with:

```shell
git branch updated-data
git checkout updated-data
git add database/db-init/*
git commit
git push origin updated-data
```

After that, people who do a new `git clone` of your fork of the Research Software
Directory, should get the updated sample data.
# How do I change when data collection scripts run?

## Relevant files

- ``harvesting/crontab``

The various scripts run according to the schedule as defined by the crontab file.

You can check the crontab file for a running instance at ``https://<hostname>/schedule``.
# How do I change the colors?

Which colors you see on the frontend is governed by the style file `frontend/static/style/rsd.scss.css`. This file is
generated from SCSS source files. In a nutshell, these are the steps to change the colors (more details below):

1. Update the SCSS variables.
1. Recompile the CSS from its SCSS source files.
1. Update the activity graph line color

## 1. Update the CSS variables

Relevant files:

- `frontend/style/settings/_variables.scss`

## 2. Recompile the CSS

Relevant files:

- `frontend/static/style/rsd.scss.css`
- `frontend/static/style/rsd.scss.css.map`

For this step, you'll need a program that can compile SCSS files into CSS. For example, on Ubuntu you can use `pysassc`:

```shell
sudo apt install pysassc
```

From the ``frontend`` root directory, run:

```shell
pysassc --style=compressed --sourcemap style/rsd.scss static/style/rsd.scss.css
```

Afterwards, your new color should be included in ``static/style/rsd.scss.css``.

## 3. Update the activity graph line color

Relevant files:

- ``frontend/static/scripts/software_scripts.js``

Update the line color and marker color to the color you want.

Refer to the [general workflow when making changes](/README.md#general-workflow-when-making-changes) to update the Docker
container with the new content.

# `frontend` service of the Research Software Directory

Visitors of the Research Software Directory website will request pages from this web application.

## Installation

Setup a virtual Python 3 environment with

```shell
pip3 install -r requirements.txt
```

## Configuration

The frontend requires the [api backend server](/backend) to be running.
The url of the backend server must be set as value for the BACKEND_URL environment variable.

Setup `.env` file, see `.env.example` as an example.

Setup environmental variables: `export $(cat .env | xargs)`

## Test

Run unit tests with fixtures:

```shell
PYTHONPATH=. pytest -m "not live"
```

You can also test against live backend server, it will check if all pages render:

```shell
BACKEND_URL=https://www.research-software.nl/api PYTHONPATH=. pytest -m live -v
```

## Run

Before running the installation and configuration must be completed.

Run in development mode with:

```shell
python entry.py
```

Run in production mode with Docker using:

```shell
docker build -t rsd/frontend .
docker run --env-file .env --rm -it --name test -p5004:5004 rsd/frontend
```

## Develop

Changes to the sass style files should be followed up by css generation with

```shell
sassc --style=compressed --sourcemap style/rsd.scss static/style/rsd.scss.css
```
