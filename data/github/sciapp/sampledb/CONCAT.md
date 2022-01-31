# SampleDB

[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![DOI](https://zenodo.org/badge/221237572.svg)](https://zenodo.org/badge/latestdoi/221237572)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02107/status.svg)](https://doi.org/10.21105/joss.02107)

SampleDB is a web-based sample and measurement metadata database.

## Documentation

You can find the documentation for the current release at https://scientific-it-systems.iffgit.fz-juelich.de/SampleDB/.

## Getting Started

We recommend using our pre-built Docker images for setting up `SampleDB`. You will need two containers, one for a PostgreSQL database and another for SampleDB itself, and a directory to store all files in.

If you would like to set up a development version of SampleDB instead, please see the [contribution guide](https://github.com/sciapp/sampledb/blob/develop/CONTRIBUTING.md).

If you do not have Docker installed yet, please [install Docker](https://docs.docker.com/engine/install/).

### Using docker-compose

First, get the [docker-compose.yml](https://raw.githubusercontent.com/sciapp/sampledb/develop/docker-compose.yml) configuration file. You can git clone this repo or just get the file:

```bash
curl https://raw.githubusercontent.com/sciapp/sampledb/develop/docker-compose.yml --output docker-compose.yml
```

Then simply bring everything up with:

```bash
docker-compose up -d
```

### Using docker commands

First, start your database container:

```bash
docker run \
    -d \
    -e POSTGRES_PASSWORD=password \
    -e PGDATA=/var/lib/postgresql/data/pgdata \
    -v `pwd`/pgdata:/var/lib/postgresql/data/pgdata:rw \
    --restart=always \
    --name sampledb-postgres \
    postgres:12
```

Next, start the SampleDB container:

```bash
docker run \
    -d \
    --link sampledb-postgres \
    -e SAMPLEDB_CONTACT_EMAIL=sampledb@example.com \
    -e SAMPLEDB_MAIL_SERVER=mail.example.com \
    -e SAMPLEDB_MAIL_SENDER=sampledb@example.com \
    -e SAMPLEDB_ADMIN_PASSWORD=password \
    -e SAMPLEDB_SQLALCHEMY_DATABASE_URI=postgresql+psycopg2://postgres:password@sampledb-postgres:5432/postgres \
    -e SAMPLEDB_FILE_STORAGE_PATH=/home/sampledb/files/ \
    -v `pwd`/files:/home/sampledb/files:rw \
    --restart=always \
    --name sampledb \
    -p 8000:8000 \
    sciapp/sampledb:0.19.3
```

### Once it's started

This will start a minimal SampleDB installation at `http://localhost:8000` and allow you to sign in with the username `admin` and the password `password` (which you should change immediately after signing in).

To learn how to further set up SampleDB, please follow the rest of the [Getting Started guide](https://scientific-it-systems.iffgit.fz-juelich.de/SampleDB/administrator_guide/getting_started.html).

## Contributing

If you want to improve SampleDB, please read the [contribution guide](https://github.com/sciapp/sampledb/blob/develop/CONTRIBUTING.md) for a few notes on how to report issues or submit changes.

## Support

If you have any questions about SampleDB or run into any issues setting up or running SampleDB, please [open an issue on GitHub](https://github.com/sciapp/sampledb/issues/new).
# Contributing to SampleDB

SampleDB is open source and we hope that ideas and improvements from different fields and institutions will help to make it as useful as possible. Thank you for considering contributing to SampleDB!

## Reporting Issues

- Describe the situation to help us reproduce it. This might include providing relevant schemas and/or objects as JSON files, scripts for WebAPI issues and/or a modified `sampledb/scripts/set_up_demo.py` that can recreate the situation.
- Describe what you expected to happen.
- Describe what actually happened. Please include logs and screenshots or screen recordings for UI issues.
- List your SampleDB version. If possible, try to reproduce your issue using the current `develop` branch.

## Setting up a Development Installation

- To work on SampleDB you will need a PostgreSQL database. On Linux, your distribution will likely have a package for this and on macOS, you can use [Postgres.app](https://postgresapp.com/).
- Set up a virtual environment using `python3 -m venv env` and activate it using `source env/bin/activate`
- Install the requirements, using `pip install -r requirements.txt`
- Set [configuration environment variables](https://scientific-it-systems.iffgit.fz-juelich.de/SampleDB/developer_guide/configuration.html). At the very least you will need to set a mail server and sender, e.g. by using `export SAMPLEDB_MAIL_SERVER=mail.example.com`, `export SAMPLEDB_MAIL_SENDER=sampledb@example.com` and `export SAMPLEDB_CONTACT_EMAIL=sampledb@example.com`. Depending on how you set up your database, you may have to set the `SAMPLEDB_SQLALCHEMY_DATABASE_URI`.
- Start an instance using demo data from the `set_up_demo` script, using `python demo.py`. This way, you will have some example instruments, actions, objects and users. If you try to access a route that requires a user account, you will automatically be signed in.

## Translating SampleDB

- SampleDB was primarily developed in English, but can be translated using [Babel](http://babel.pocoo.org/en/latest/cmdline.html) and the *gettext* internationalization and localization system.
- A new translation can be created using `pybabel extract -F babel.cfg -k lazy_gettext -o sampledb/messages.pot sampledb` followed by `pybabel init -i sampledb/messages.pot -d sampledb/translations -D extracted_messages -l <locale>`.
- The locale for this translation also needs to be added to `sampledb.logic.locales.SUPPORTED_LOCALES`, so users will be able to select it.
- Additionally, a language should be created for it, so user-generated content can be translated as well.
- If you would like to translate SampleDB to a different language, please open up an [issue](https://github.com/sciapp/sampledb/issues/new/choose) to coordinate with us.

## Submitting Changes

- We aim to support a wide variety of use cases from different fields. Please keep this in mind and prefer generic solutions to specialized ones.
- Adhere to the code style. You can check this using `python3 -m pycodestyle --ignore=E402,E501,W504 sampledb` and `python3 -m pyflakes sampledb`.
- Include tests if your change introduces a new feature or fixes a bug. Make sure that the test fails
  without your change. You can run all tests using `python3 -m pytest -s --cov=sampledb/ tests`. As the tests require an LDAP server, you will need to set the following [configuration environment variables](https://scientific-it-systems.iffgit.fz-juelich.de/SampleDB/developer_guide/configuration.html):
  - `SAMPLEDB_LDAP_SERVER`
  - `SAMPLEDB_LDAP_USER_BASE_DN`
  - `SAMPLEDB_LDAP_UID_FILTER`
  - `SAMPLEDB_LDAP_NAME_ATTRIBUTE`
  - `SAMPLEDB_LDAP_MAIL_ATTRIBUTE`
  - `SAMPLEDB_LDAP_OBJECT_DEF`
  - `SAMPLEDB_TESTING_LDAP_LOGIN`
  - `SAMPLEDB_TESTING_LDAP_PW`
- If you introduce or significantly improve a feature, please document it and add it to the changelog in `docs/changelog.rst`. You can build the documentation using `python -m sphinx docs/ build/`.
- Once you are done, push your changes to your fork of SampleDB and open up a [pull request](https://github.com/sciapp/sampledb/compare).
---
title: 'SampleDB: A sample and measurement metadata database'
tags:
  - Python
  - sample management
  - research data management
authors:
  - name: Florian Rhiem
    orcid: 0000-0001-6461-9433
    affiliation: 1
affiliations:
  - name: PGI/JCNS-TA, Forschungszentrum Jülich
    index: 1
date: 10 December 2019
bibliography: paper.bib
---

# Summary

One of the key aspects of good scientific practice is the handling of research data [@dfg]. Archiving research data and making it findable, accessible, interoperable and re-usable by other researchers is crucial for reproducibility and can also yield new findings [@fair]. This is not only true for research data resulting from experiments or simulations, but particularly for information on how a sample was created, how a measurement was performed and which parameters were used for a simulation.

``SampleDB`` is a web-based sample and measurement metadata database developed at Jülich Centre for Neutron Science (JCNS) and Peter Grünberg Institute (PGI). Researchers can use ``SampleDB`` to store and retrieve information on samples, measurements and simulations, analyze them using Jupyter notebooks, track sample storage locations and responsibilities and view sample life cycles.

The application was designed to support the wide variety of instruments and processes found at PGI, JCNS and other institutes by allowing users to define the metadata of new processes using a graphical editor or a JSON-based schema language. These schemas can contain common datatypes such as booleans, texts and arrays, but also more specialized datatypes such as physical quantities and references to existing samples. Using these schemas, the application can generate forms for entering the information and then validate the information, including support for using regular expressions as constraints for text input. By storing the information using the ``JSONB`` datatype of ``PostgreSQL`` [@postgresql], the built-in search function can support complex expressions and search for quantities with differing units but equal dimensionality.

In addition to the process-specific metadata, users can also store other information, including the location of a sample, auxiliary files, publications and additional comments. To provide a complete history of the metadata on an object, all previous versions are stored and can be accessed by users.

A Web API allows automated data entry using already existing information, e.g. by integrating it into an instrument control system or by monitoring log files, thereby archiving the information without additional overhead for the scientists. The Web API can also be used for automated data retrieval, e.g. for accessing measurement parameters during data analysis.

The latter is further simplified if ``SampleDB`` is used in conjunction with a ``JupyterHub`` server. Instrument scientists can provide Jupyter notebook templates for analyzing data from their instruments, along with a list of schema entries containing the necessary metadata. Users can then use these templates to create Jupyter notebooks from within ``SampleDB``, which will copy the required metadata from its ``SampleDB`` entry into the notebook, preparing a process-specific data analysis.

Using a schema system for process-specific metadata allows ``SampleDB`` to support a wide variety of processes and presents an advantage over alternative applications, such as the JuliaBase project [@juliabase], which define instruments and processes as part of the application code. That approach inherently limits the group of users that can define or alter processes to those users with administrative access. Keeping application logic and process-specific metadata separate has the disadvantage that ``SampleDB`` cannot provide built-in data analysis tools the way other applications such as the JuliaBase project or BikaLIMS [@bikalims] / Senaite [@senaite] can, however the support for Jupyter notebook templates can be used as an alternative to these data analysis tools. Using pre-defined calculations as they are offered by BikaLIMS / Senaite could also be explored for ``SampleDB`` in the future.

``SampleDB`` was developed using Python 3, the Flask web framework [@flask] and the SQLAlchemy package [@sqlalchemy]. To run ``SampleDB``, we recommend using the provided [Docker images](https://hub.docker.com/r/sciapp/sampledb/) along with a PostgreSQL container, though it can also be run directly from source. Information on setting up a ``SampleDB`` instance can be found on the [GitHub project site](https://github.com/sciapp/sampledb/) and a more in-depth guide can be found in the [project documentation](https://scientific-it-systems.iffgit.fz-juelich.de/SampleDB/).

# Acknowledgements

We thank Dorothea Henkel for her contributions, Daniel Kaiser for code review, and Paul Zakalek and Jörg Perßon for their feedback during the development of ``SampleDB``.

# References
.. raw:: html

    <div class="table_of_contents_row">
    <div class="table_of_contents_block">

User Guide
----------

.. toctree::
    :maxdepth: 2

    user_guide/users.rst
    user_guide/groups.rst
    user_guide/projects.rst
    user_guide/instruments.rst
    user_guide/actions.rst
    user_guide/schemas.rst
    user_guide/objects.rst
    user_guide/search.rst
    user_guide/export.rst
    user_guide/citations.rst

.. raw:: html

    </div>
    <div class="table_of_contents_block">

Administrator Guide
-------------------

.. toctree::
    :maxdepth: 1

    administrator_guide/getting_started.rst
    administrator_guide/configuration.rst
    administrator_guide/tls_termination.rst
    administrator_guide/backup_and_restore.rst
    administrator_guide/upgrading.rst
    administrator_guide/administration_scripts.rst
    administrator_guide/jupyterhub_support.rst
    administrator_guide/dataverse_export.rst
    administrator_guide/languages.rst

Developer Guide
---------------

.. toctree::
    :maxdepth: 2

    developer_guide/api.rst
    developer_guide/contributing.rst

.. raw:: html

    </div>
    </div>

.. toctree::
    :hidden:
    :maxdepth: 2

    changelog.rst

.. todolist::

.. note::
    Both |service_name| and this documentation are still a work in progress. If you come across any issues or want to request a feature, please `let us know`_.
Changelog
=========

Version 0.20
------------

Currently in development.

- Added support for any, all and not conditions
- Add schema templates

Version 0.19.3
--------------

Released on January 19th, 2022.

- Fix schema upgrade for multi language choices

Version 0.19.2
--------------

Released on January 7th, 2022.

- Fix editing notes in schema editor

Version 0.19.1
--------------

Released on December 20th, 2021.

- Fix missing object type and ID on object page when using inline edit mode

Version 0.19
------------

Released on December 9th, 2021.

- Allow filtering instrument log entries by author
- Allow sorting instrument log entries by author
- Added event datetime for instrument log entries
- Added internationalization features
- Added german localization
- Store file contents in database by default
- Allow setting a publicly visible user role
- Added support for configurable user fields
- Added label for administrators in user list
- Allow individual exemptions for Use as Template
- Allow setting a default number of items for arrays
- Improved user interface
- Added support for a custom CSS file
- Added support for conditional properties
- Allow filtering object references by action
- Implemented TOTP-based two factor authentication
- Added tree view for instrument log entries
- Allow editing individual fields
- Allow hiding object type and id on object page

Version 0.18
------------

Released on May 7th, 2021.

- Moved example_data functionality to set_up_data script
- Allow administrators to enforce user names to be given as surname, given names
- Added plotly_chart data type
- Improved search page
- Improved object version HTTP API
- Improved action HTTP API
- Improved user interface

Version 0.17
------------

Released on February 10th, 2021.

- Added Dataverse export using the EngMeta "Process Metadata" block
- Added short descriptions to actions and instruments
- Added array style "horizontal_table"
- Improved handling of optional text input
- Allow linking to headers in Markdown content
- Allow disabling of "Use in Measurement" button for samples
- Added markdown support to object metadata
- Added markdown support to instrument log
- Reimplemented PDF export
- Added configuration variables to allow only administrators to create groups or projects
- Added asterisks to mark required fields when editing objects
- Project permissions can be set when inviting a user
- Allow default value "self" for user fields
- Allow searching for tags in dropdown object selection fields
- Renamed projects to project groups and groups to basic groups to avoid ambiguity
- Allow disabling of subprojects / child project groups
- Allow giving basic or project groups initial permissions
- Allow configuring the Help link
- Allow linking project groups to objects
- Fixed action ID filtering when loading objects in the background
- Added action permissions to user interface
- Improved handling of quantities for the HTTP API

Version 0.16.1
--------------

Released on January 27th, 2021.

- Fixed object name escaping when loading objects in the background

Version 0.16
------------

Released on December 9th, 2020.

- Allow restricting object references to specific action id
- Improved performance of object lists
- Allow setting display properties as part of the object list URL
- Improved performance of instrument pages
- Added image upload via drag and drop to Markdown editors
- Added support for placeholder texts for text and quantity schemas
- Added additional options to the HTTP API objects endpoint
- Display projects based on parent-child relationship
- Improved "View Objects" for users, groups and projects
- Added object comments to the HTTP API

Version 0.15
------------

Released on November 6th, 2020.

- Added versioning to instrument log entries
- Added user to metadata types
- Allow setting instrument log entry order
- Allow custom action types
- Allow administrators to deactivate users
- Allow disabling group deletion by non-administrators
- Fixed pagination for viewing objects of a project
- Added Docker Compose configuration file
- Ensure that file storage path is owned by sampledb user in docker container
- Added ``SAMPLEDB_LOAD_OBJECTS_IN_BACKGROUND`` option to load object select options using ajax
- Added "list" array style
- Added Markdown editor for editing instrument and action Markdown content

Version 0.14.1
--------------

Released on October 13th, 2020.

- Upgraded dependencies

Version 0.14
------------

Released on September 23rd, 2020.

- Allow restricting location management to administrators
- Do not show hidden users as instrument scientists
- Added setting for admin permissions
- Allow hiding instruments and actions
- Added object name to properties of publications linked to an object
- Improved invitation token handling
- Made invitation time limit configurable
- Show pending group and project invitations to members
- Show all group and project invitations to administrators
- Allow copying permissions from another object
- Improved user interface

Version 0.13.1
--------------

Released on September 9th, 2020.

- Fixed a user interface issue

Version 0.13
------------

Released on September 2nd, 2020.

- Added Dublin Core metadata in RDF/XML format
- Added fullscreen image preview of object and instrument log images
- Added instrument log to HTTP API
- Allow filtering instrument log by month
- Allow setting a publicly visible user affiliation

Version 0.12
------------

Released on July 29th, 2020.

- Added data export as PDF document, .zip or .tar.gz archive
- Allow adding a logo to object export PDF documents
- Allow setting a publicly visible ORCID iD
- Added instrument log
- Added instrument scientist notes

Version 0.11
------------

Released on June 18th, 2020.

- Allow usage of Markdown in instrument and action descriptions
- Added configuration values for creating an admin user during initial setup
- Added administrator guide to documentation

Version 0.10
------------

Released on May 11th, 2020.

- Allow configuring label formats
- Added search filters to objects API

Version 0.9
-----------

Released on March 10th, 2020.

- Allow creating and editing instruments using the web frontend
- Allow referencing measurements as object properties
- Added readonly users
- Allow hiding users
- Added API tokens
- Added administration functions to the web frontend
- Fixed various minor bugs

Version 0.8.1
-------------

Released on December 10th, 2019.

- Simplified deployment

Version 0.8
-----------

Released on November 12th, 2019.

- Added search to group and project dialogs
- Fixed various minor bugs


Version 0.7
-----------

Released on September 13th, 2019.

- Allow deleting groups and projects
- Allow group and project member removal
- Allow users to accept responsibility assignments
- Fixed various minor bugs


Version 0.6
-----------

Released on August 30th, 2019.

- Added JupyterHub notebook templates
- Added list of tags
- Fixed various minor bugs


Version 0.5
-----------

Released on April 15th, 2019.

- Added publications
- Removed activity log
- Added files to HTTP API
- Improved user interface


Version 0.4
-----------

Released on February 13th, 2019.

- Added object pagination
- Added posting of external links for objects
- Added schema editor
- Added 'Use in Measurement' button to samples
- Fixed various minor bugs


Version 0.3.1
-------------

Released on January 21st, 2019.

- Improved performance of object permissions


Version 0.3
-----------

Released on January 16th, 2019.

- Added custom actions
- Added locations
- Added notifications
- Added search by user name
- Added users and object permissions to HTTP API
- Improved documentation
- Improved email design
- Improved user interface
- Fixed various minor bugs


Version 0.2
-----------

Released on November 30th, 2018.

- Added documentation
- Added HTTP API
- Added *Related Objects* to objects' pages
- Added PDF export for objects
- Added label generation for objects
- Added GHS hazards as optional metadata
- Added error messages during object creation and editing
- Changed advanced search to be automatic for some queries
- Added sorting to object tables
- Added favorites for actions and instruments
- Improved user interface
- Fixed various minor bugs

Version 0.1
-----------

First stable release.
.. _contributing:

Contributing
============

You can find information on how to report issues, set up a development installation or submit changes in the `Contribution Guide <https://github.com/sciapp/sampledb/blob/develop/CONTRIBUTING.md>`_ in the `SampleDB GitHub repository <https://github.com/sciapp/sampledb>`_.
.. _http_api:

HTTP API
========

Authentication
--------------

The |service_name| HTTP API either uses `Basic Authentication <https://tools.ietf.org/html/rfc7617>`_ using normal user credentials (e.g. using the header :code:`Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=`) or `Bearer Authentication <https://tools.ietf.org/html/rfc6750>`_ using the API token (e.g. using the header :code:`Authorization: Bearer bf4e16afa966f19b92f5e63062bd599e5f931faeeb604bdc3e6189539258b155`). API tokens are meant as an alternative method for authentication for individual scripts and allow you to monitor the requests made with the token. You can create an API token when editing your :ref:`preferences`. If you have a two factor authentication method enabled, you cannot use your user credentials to use the API and will have to use an API token instead.

Please make sure to use HTTPS when accessing the API.

Objects
-------

Reading a list of all objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/

    Get a list of all objects visible to the current user.

    The list only contains the current version of each object. By passing the parameter :code:`q` to the query, the :ref:`advanced_search` can be used. By passing the parameters :code:`action_id` or :code:`action_type` objects can be filtered by the action they were created with or by their type (e.g. :code:`sample` or :code:`measurement`).

    Instead of returning all objects, the parameters :code:`limit` and :code:`offset` can be used to reduce to maximum number of objects returned and to provide an offset in the returned set, so allow simple pagination.

    If the parameter :code:`name_only` is provided, the object data and schema will be reduced to the name property, omitting all other properties and schema information.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "object_id": 1,
                "version_id": 0,
                "action_id": 0,
                "schema": {
                    "title": "Object Information",
                    "type": "object",
                    "properties": {
                        "name": {
                            "title": "Object Name",
                            "type": "text"
                        }
                    }
                },
                "data": {
                    "name": {
                        "_type": "text",
                        "text": "Example Object"
                    }
                }
            },
            {
                "object_id": 2,
                "version_id": 3,
                "action_id": 0,
                "schema": {
                    "title": "Object Information",
                    "type": "object",
                    "properties": {
                        "name": {
                            "title": "Object Name",
                            "type": "text"
                        }
                    }
                },
                "data": {
                    "name": {
                        "_type": "text",
                        "text": "Other Object"
                    }
                }
            }
        ]

    :statuscode 200: no error


Getting the current object version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)

    Redirect to the current version of an object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 302 Found
        Location: /api/v1/objects/1/versions/0

    :statuscode 302: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading an object version
^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/versions/(int:version_id)

    Get the specific version (`version_id`) of an object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/versions/0 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "object_id": 1,
            "version_id": 0,
            "action_id": 0,
            "user_id": 1,
            "utc_datetime": "2021-04-29 12:34:56",
            "schema": {
                "title": "Object Information",
                "type": "object",
                "properties": {
                    "name": {
                        "title": "Object Name",
                        "type": "text"
                    }
                }
            },
            "data": {
                "name": {
                    "_type": "text",
                    "text": "Example Object"
                }
            }
        }

    :>json number object_id: the object's ID
    :>json number version_id: the object version's ID
    :>json number action_id: the action's ID
    :>json object action: the action (if the parameter embed_action is set to a non-empty value)
    :>json number user_id: the ID of the user who created this version
    :>json object user: the user (if the parameter embed_user is set to a non-empty value)
    :>json string utc_datetime: the time and date when this version was created in UTC
    :>json object schema: the object's schema
    :>json object data: the object's data
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object/version combination does not exist


Creating a new object
^^^^^^^^^^^^^^^^^^^^^

.. http:post:: /api/v1/objects/

    Create a new object.

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/objects/1/versions/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Content-Type: application/json
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "action_id": 0,
            "schema": {
                "title": "Object Information",
                "type": "object",
                "properties": {
                    "name": {
                        "title": "Object Name",
                        "type": "text"
                    }
                }
            },
            "data": {
                "name": {
                    "_type": "text",
                    "text": "Example Object"
                }
            }
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: /api/v1/objects/1/versions/0

    :<json number version_id: the object version's ID (optional, must be 0)
    :<json number action_id: the action's ID
    :<json object schema: the object's schema (optional, must equal current action's schema)
    :<json object data: the object's data
    :statuscode 201: no error
    :statuscode 400: invalid data


Updating an object / Creating a new object version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:post:: /api/v1/objects/(int:object_id)/versions/

    Create a new version of an object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/objects/1/versions/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Content-Type: application/json
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "data": {
                "name": {
                    "_type": "text",
                    "text": "Example Object"
                }
            }
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: /api/v1/objects/1/versions/1

    :<json number object_id: the object's ID (optional, must equal `object_id` in URL)
    :<json number version_id: the object version's ID (optional, must equal new version's ID)
    :<json number action_id: the action's ID (optional, must equal previous `action_id`)
    :<json object schema: the object's schema (optional, must equal previous `schema` or current action's schema)
    :<json object data: the object's data
    :statuscode 201: no error
    :statuscode 400: invalid data
    :statuscode 403: the user does not have WRITE permissions for this object
    :statuscode 404: the object does not exist


Object Permissions
------------------


Reading whether an object is public
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/public

    Get whether or not an object is public.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/public HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        true

    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Setting whether an object is public
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:put:: /api/v1/objects/(int:object_id)/permissions/public

    Get whether or not an object is public.

    **Example request**:

    .. sourcecode:: http

        PUT /api/v1/objects/1/permissions/public HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        false

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        false

    :statuscode 200: no error
    :statuscode 403: the user does not have GRANT permissions for this object
    :statuscode 404: the object does not exist


Reading all users' permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/users/

    Get a mapping of user IDs to their permissions.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/users/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "1": "read",
            "2": "grant"
        }

    :queryparam include_instrument_responsible_users: If given, permissions from being an instrument responsible user will be included (optional)
    :queryparam include_groups: If given, permissions from basic group memberships will be included (optional)
    :queryparam include_projects: If given, permissions from project group memberships will be included (optional)
    :queryparam include_admins: If given, permissions from being an administrator will be included (optional)
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading a user's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/users/(int:user_id)

    Get the permissions of a user for an object.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/users/2 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "grant"

    :queryparam include_instrument_responsible_users: If given, permissions from being an instrument responsible user will be included (optional)
    :queryparam include_groups: If given, permissions from basic group memberships will be included (optional)
    :queryparam include_projects: If given, permissions from project group memberships will be included (optional)
    :queryparam include_admins: If given, permissions from being an administrator will be included (optional)
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object or user does not exist


Setting a user's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:put:: /api/v1/objects/(int:object_id)/permissions/users/(int:user_id)

    Set the permissions of a user for an object.

    **Example request**:

    .. sourcecode:: http

        PUT /api/v1/objects/1/permissions/users/2 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        "write"

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "write"

    :statuscode 200: no error
    :statuscode 400: invalid data (should be "read", "write", "grant" or "none")
    :statuscode 403: the user does not have GRANT permissions for this object
    :statuscode 404: the object or user does not exist


Reading all basic groups' permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/groups/

    Get a mapping of basic group IDs to their permissions.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/groups/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "4": "write"
        }

    :queryparam include_projects: If given, permissions from project group memberships will be included (optional)
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading a basic group's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/groups/(int:group_id)

    Get the permissions of a basic group for an object.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/groups/4 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "write"

    :queryparam include_projects: If given, permissions from project group memberships will be included (optional)
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object or basic group does not exist


Setting a basic group's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:put:: /api/v1/objects/(int:object_id)/permissions/groups/(int:group_id)

    Set the permissions of a basic group for an object.

    **Example request**:

    .. sourcecode:: http

        PUT /api/v1/objects/1/permissions/groups/2 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        "read"

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "read"

    :statuscode 200: no error
    :statuscode 400: invalid data (should be "read", "write", "grant" or "none")
    :statuscode 403: the user does not have GRANT permissions for this object
    :statuscode 404: the object or basic group does not exist


Reading all project groups' permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/projects/

    Get a mapping of project group IDs to their permissions.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/projects/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "7": "read"
        }

    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading a project group's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/permissions/projects/(int:project_id)

    Get the permissions of a project group for an object.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/permissions/projects/7 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "read"

    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object or project group does not exist


Setting a project group's permissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:put:: /api/v1/objects/(int:object_id)/permissions/projects/(int:project_id)

    Set the permissions of a project group for an object.

    **Example request**:

    .. sourcecode:: http

        PUT /api/v1/objects/1/permissions/projects/2 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        "read"

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        "read"

    :statuscode 200: no error
    :statuscode 400: invalid data (should be "read", "write", "grant" or "none")
    :statuscode 403: the user does not have GRANT permissions for this object
    :statuscode 404: the object or project group does not exist


Instruments
-----------


Reading a list of all instruments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/

    Get a list of all instruments.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "instrument_id": 1,
                "name": "Example Instrument",
                "description": "This is an example instrument",
                "is_hidden": false,
                "instrument_scientists": [1, 42]
            }
        ]

    :statuscode 200: no error


Reading an instrument
^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)

    Get the specific instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "instrument_id": 1,
            "name": "Example Instrument",
            "description": "This is an example instrument",
            "is_hidden": false,
            "instrument_scientists": [1, 42]
        }

    :>json number instrument_id: the instrument's ID
    :>json string name: the instruments's name
    :>json string description: the instruments's description
    :>json bool is_hidden: whether or not the instrument is hidden
    :>json list instrument_scientists: the instrument scientists' IDs
    :statuscode 200: no error
    :statuscode 404: the instrument does not exist


Instrument Log Entries
----------------------

Reading a list of all log entries for an instrument
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/

    Get a list of all log entries for a specific instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "log_entry_id": 1,
                "utc_datetime": "2020-08-19T12:13:14.123456",
                "author": 1,
                "content": "Example Log Entry 1",
                "categories": []
            },
            {
                "log_entry_id": 2,
                "utc_datetime": "2020-08-19T13:14:15.123456",
                "author": 1,
                "content": "Example Log Entry 2",
                "categories": [
                    {
                        "category_id": 1
                        "title": "Error Report"
                    },
                    {
                        "category_id": 7
                        "title": "Maintenance Log"
                    }
                ]
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument does not exist


Reading an instrument log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/(int:log_entry_id)

    Get the specific log entry (`log_entry_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries/2 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "log_entry_id": 2,
            "utc_datetime": "2020-08-19T13:14:15.123456",
            "author": 1,
            "content": "Example Log Entry 2",
            "categories": [
                {
                    "category_id": 1
                    "title": "Error Report"
                },
                {
                    "category_id": 7
                    "title": "Maintenance Log"
                }
            ]
        }

    :>json number log_entry_id: the log entry's ID
    :>json string utc_datetime: the date and time of the log entry in UTC in ISO format
    :>json string content: the log entry's content
    :>json number author: the user ID of the log entry's author
    :>json list categories: the log entry's categories
    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument or the log entry do not exist


Reading a list of all log categories for an instrument
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_categories/

    Get a list of all log categories for a specific instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_categories HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "category_id": 1
                "title": "Error Report"
            },
            {
                "category_id": 7
                "title": "Maintenance Log"
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument does not exist


Reading an instrument log category
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_categories/(int:category_id)

    Get the specific log category (`category_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_categories/7 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "category_id": 7
            "title": "Maintenance Log"
        }

    :>json number category_id: the log category's ID
    :>json string title: the log category's title
    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument or the log category do not exist


Reading a list of all file attachments for a log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/(int:log_entry_id)/file_attachments/

    Get a list of file attachments for a specific log entry (`log_entry_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries/2/file_attachments HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "file_attachment_id": 1,
                "file_name": "example.txt",
                "content": "RXhhbXBsZSBDb250ZW50"
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument or the log entry do not exist


Reading a file attachment for a log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/(int:log_entry_id)/file_attachments/(int:file_attachment_id)

    Get a specific file attachment (`file_attachment_id`) for a log entry (`log_entry_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries/2/file_attachments/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "file_attachment_id": 1,
            "file_name": "example.txt",
            "content": "RXhhbXBsZSBDb250ZW50"
        }

    :>json string file_attachment_id: the file attachment's ID
    :>json string file_name: the original file name
    :>json string content: the base64 encoded file content
    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument, the log entry or the file attachment do not exist


Reading a list of all object attachments for a log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/(int:log_entry_id)/object_attachments/

    Get a list of object attachments for a specific log entry (`log_entry_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries/2/object_attachments HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "object_attachment_id": 1,
                "object_id": 1
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument or the log entry do not exist


Reading an object attachment for a log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/instruments/(int:instrument_id)/log_entries/(int:log_entry_id)/object_attachments/(int:object_attachment_id)

    Get a specific object attachment (`object_attachment_id`) for a log entry (`log_entry_id`) for an instrument (`instrument_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/instruments/1/log_entries/2/object_attachments/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "object_attachment_id": 1,
            "object_id": 1
        }

    :>json string object_attachment_id: the object attachment's ID
    :>json string object_id: the object ID
    :statuscode 200: no error
    :statuscode 403: the instrument log can only be accessed by instrument scientists
    :statuscode 404: the instrument, the log entry or the object attachment do not exist


Creating an instrument log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:post:: /api/v1/instruments/(int:instrument_id)/log_entries/

    Create a log entry for an instrument (`instrument_id`) and optionally attach files and objects to it.

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/instruments/1/log_entries/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "content": "Example Log Entry Text",
            "category_ids": [1, 7],
            "file_attachments": [
                {
                    "file_name": "example.txt",
                    "base64_content": "RXhhbXBsZSBDb250ZW50"
                }
            ],
            "object_attachments": [
                {
                    "object_id": 1
                },
                {
                    "object_id": 2
                }
            ]
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: https://iffsamples.fz-juelich.de/api/v1/instruments/1/log_entries/1

    :<json string content: the log entry's content
    :<json list category_ids: an optional list of category IDs for the log entry
    :<json list file_attachments: an optional list of file attachments as json objects with file_name and base64_content attributes
    :<json list object_attachments: an optional list of object attachments as json objects with an object_id attribute
    :statuscode 201: the log entry and optional attachments have been created successfully
    :statuscode 400: there was an error in the given json data
    :statuscode 403: only instrument scientists can write to the instrument log
    :statuscode 404: the instrument does not exist


Actions
-------


Reading a list of all actions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/actions/

    Get a list of all actions.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/actions/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "action_id": 1,
                "instrument_id": null,
                "user_id": null,
                "type": "sample",
                "type_id": -99,
                "name": "Example Sample Creation",
                "description": "This is an example action",
                "is_hidden": false,
                "schema": {
                    "title": "Example Sample",
                    "type": "object",
                    "properties": {
                        "name": {
                            "title": "Sample Name",
                            "type": "text"
                        }
                    },
                    "required": ["name"]
                }
            },
            {
                "action_id": 2,
                "instrument_id": 1,
                "user_id": null,
                "type": "measurement",
                "type_id": -98,
                "name": "Example Measurement",
                "description": "This is an example action",
                "is_hidden": false,
                "schema": {
                    "title": "Example Measurement",
                    "type": "object",
                    "properties": {
                        "name": {
                            "title": "Measurement Name",
                            "type": "text"
                        }
                    },
                    "required": ["name"]
                }
            }
        ]

    :statuscode 200: no error


Reading an action
^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/actions/(int:action_id)

    Get the specific action (`action_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/actions/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "action_id": 1,
            "instrument_id": null,
            "user_id": null,
            "type": "sample",
            "type_id": -99,
            "name": "Example Sample Creation",
            "description": "This is an example action",
            "is_hidden": false,
            "schema": {
                "title": "Example Sample",
                "type": "object",
                "properties": {
                    "name": {
                        "title": "Sample Name",
                        "type": "text"
                    }
                },
                "required": ["name"]
            }
        }

    :>json number action_id: the action's ID
    :>json number instrument_id: the action's instrument's ID or null
    :>json number user_id: the action's user ID, if it is a user-specific action, or null
    :>json string type: the action's type ("sample", "measurement", "simulation" or "custom")
    :>json number type_id: the ID of the action's type
    :>json string name: the action's name
    :>json string description: the action's description
    :>json bool is_hidden: whether or not the action is hidden
    :>json object schema: the action's schema
    :statuscode 200: no error
    :statuscode 404: the action does not exist


Action Types
------------


Reading a list of all action types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/action_types/

    Get a list of all action types.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/action_types/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "type_id": -99,
                "name": "Sample Creation",
                "object_name": "sample",
                "admin_only": false
            },
            {
                "type_id": -98,
                "name": "Measurement",
                "object_name": "measurement",
                "admin_only": false
            },
            {
                "type_id": -97,
                "name": "Simulation",
                "object_name": "simulation",
                "admin_only": false
            }
        ]

    :statuscode 200: no error


Reading an action type
^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/action_types/(int:type_id)

    Get the specific action type (`type_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/action_types/-99 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "type_id": -99,
            "name": "Sample Creation",
            "object_name": "sample",
            "admin_only": false
        }

    :>json number type_id: the action type's ID
    :>json string name: the action type's name
    :>json string object_name: the name of objects created with this action type
    :>json bool admin_only: whether or not actions with this type can only be created by administrators
    :statuscode 200: no error
    :statuscode 404: the action does not exist


Users
-----


Reading a list of all users
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/users/

    Get a list of all users.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/users/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "user_id": 1,
                "name": "Example User",
                "orcid": null,
                "affiliation": null,
                "role": null
            }
        ]

    :statuscode 200: no error


Reading a user
^^^^^^^^^^^^^^

.. http:get:: /api/v1/users/(int:user_id)

    Get the specific user (`user_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/users/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "user_id": 1,
            "name": "Example User",
            "orcid": null,
            "affiliation": null,
            "role": null
        }

    :>json number user_id: the user's ID
    :>json string name: the user's name
    :>json string orcid: the user's ORCid ID (optional)
    :>json string affiliation: the user's affiliation (optional)
    :>json string role: the user's role (optional)
    :>json string email: the user's email (only for API requests by administrators)
    :statuscode 200: no error
    :statuscode 404: the user does not exist


Reading the current user
^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/users/me

    Get the current user.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/users/me HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "user_id": 1,
            "name": "Example User",
            "orcid": null,
            "affiliation": null,
            "role": null
        }

    :>json number user_id: the user's ID
    :>json string name: the user's name
    :>json string orcid: the user's ORCid ID (optional)
    :>json string affiliation: the user's affiliation (optional)
    :>json string role: the user's role (optional)
    :>json string email: the user's email (only for API requests by administrators)
    :statuscode 200: no error


Locations
---------


Reading a list of all locations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/locations/

    Get a list of all locations.

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/locations/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "location_id": 1,
                "name": "Example Location",
                "description": "This is an example location",
                "parent_location_id": null
            }
        ]

    :statuscode 200: no error


Reading a location
^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/locations/(int:location_id)

    Get the specific location (`location_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/locations/1 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "location_id": 1,
            "name": "Example Location",
            "description": "This is an example location",
            "parent_location_id": null
        }

    :>json number location_id: the location's ID
    :>json string name: the locations's name
    :>json string description: the locations's description
    :>json number parent_location_id: the parent location's ID
    :statuscode 200: no error
    :statuscode 404: the location does not exist


Reading a list of an object's locations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/locations/

    Get a list of all object locations assignments for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/locations/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "object_id": 1,
                "location_id": 3,
                "responsible_user_id": 6,
                "user_id": 17,
                "description": "Shelf C",
                "utc_datetime": "2018-12-11 17:50:00"
            }
        ]

    :statuscode 200: no error


Reading an object's location
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/locations/(int:index)

    Get a specific object location assignment (`index`) for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/locations/0 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "object_id": 1,
            "location_id": 3,
            "responsible_user_id": 6,
            "user_id": 17,
            "description": "Shelf C",
            "utc_datetime": "2018-12-11 17:50:00"
        }

    :>json number object_id: the object's ID
    :>json number location_id: the location's ID
    :>json number responsible_user_id: the ID of the user who is responsible for the object
    :>json number user_id: the ID of the user who assigned this location to the object
    :>json string description: the description of the object's position
    :>json number utc_datetime: the datetime when the object was stored
    :statuscode 200: no error
    :statuscode 404: the object or the object location assignment does not exist


Files
-----


Reading a list of an object's files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/files/

    Get a list of all files for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/files/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "object_id": 1,
                "file_id": 0,
                "storage": "url",
                "url": "https://iffsamples.fz-juelich.de"
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading information for a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/files/(int:file_id)

    Get a specific file (`file_id`) for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/files/0 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "object_id": 1,
            "file_id": 0,
            "storage": "url",
            "url": "https://iffsamples.fz-juelich.de"
        }

    :>json number object_id: the object's ID
    :>json number file_id: the file's ID
    :>json string storage: how the file is stored (local, database or url)
    :>json string url: the URL of the file (for url storage)
    :>json string original_file_name: the original name of the file (for local or database storage)
    :>json string base64_content: the base64 encoded content of the file (for local or database storage)
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object or the file does not exist


Uploading a file
^^^^^^^^^^^^^^^^

.. http:post:: /api/v1/objects/(int:object_id)/files/

    Create a new file with local storage for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/objects/1/files/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "storage": "local",
            "original_file_name": "test.txt",
            "base64_content": "dGVzdA=="
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: https://iffsamples.fz-juelich.de/api/v1/objects/1/files/0

    :<json string storage: how the file is stored (local)
    :<json string original_file_name: the original name of the file
    :<json string base64_content: the base64 encoded content of the file
    :statuscode 201: the file has been created successfully
    :statuscode 403: the user does not have WRITE permissions for this object
    :statuscode 404: the object does not exist


Posting a link
^^^^^^^^^^^^^^

.. http:post:: /api/v1/objects/(int:object_id)/files/

    Create a new file with url storage for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/objects/1/files/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "storage": "url",
            "url": "https://iffsamples.fz-juelich.de"
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: https://iffsamples.fz-juelich.de/api/v1/objects/1/files/0

    :<json string storage: how the file is stored (url)
    :<json string url: the URL of the file
    :statuscode 201: the file has been created successfully
    :statuscode 403: the user does not have WRITE permissions for this object
    :statuscode 404: the object does not exist

Comments
--------


Reading a list of an object's comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/comments/

    Get a list of all comments for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/comments/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        [
            {
                "object_id": 1,
                "user_id": 1,
                "comment_id": 0,
                "content": "This is an example comment"
                "utc_datetime": "2020-12-03T01:02:03.456789"
            }
        ]

    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object does not exist


Reading information for a comment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. http:get:: /api/v1/objects/(int:object_id)/comments/(int:comment_id)

    Get a specific comment (`comment_id`) for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        GET /api/v1/objects/1/comments/0 HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 200 OK
        Content-Type: application/json

        {
            "object_id": 1,
            "user_id": 1,
            "comment_id": 0,
            "content": "This is an example comment"
            "utc_datetime": "2020-12-03T01:02:03.456789"
        }

    :>json number object_id: the object's ID
    :>json number user_id: the ID of the user who posted the comment
    :>json number comment_id: the comment's ID
    :>json string content: the comment's content
    :>json string utc_datetime: the time the comment was posted in UTC formatted in ISO 8601 format
    :statuscode 200: no error
    :statuscode 403: the user does not have READ permissions for this object
    :statuscode 404: the object or the comment does not exist


Posting a comment
^^^^^^^^^^^^^^^^^

.. http:post:: /api/v1/objects/(int:object_id)/comments/

    Create a new comment for a specific object (`object_id`).

    **Example request**:

    .. sourcecode:: http

        POST /api/v1/objects/1/comments/ HTTP/1.1
        Host: iffsamples.fz-juelich.de
        Accept: application/json
        Authorization: Basic dXNlcm5hbWU6cGFzc3dvcmQ=

        {
            "content": "This is an example comment"
        }

    **Example response**:

    .. sourcecode:: http

        HTTP/1.1 201 Created
        Content-Type: application/json
        Location: https://iffsamples.fz-juelich.de/api/v1/objects/1/comments/0

    :<json string content: the (non-empty) content for the new comment
    :statuscode 201: the comment has been created successfully
    :statuscode 403: the user does not have WRITE permissions for this object
    :statuscode 404: the object does not exist
.. _export:

Export
======

|service_name| can export the data for all objects a :ref:`User <users>` has **READ** permissions for, either as a PDF file or as an archive containing a machine-readable JSON file and the files uploaded for these objects.

The PDF file consists of the :ref:`PDF files <pdf_export>` for the individual objects. As a result, only the current metadata versions and only those uploaded files that can be represented in the PDF are included.

The archive contains a file ``data.json`` with information on objects in JSON format, including all metadata versions, comments, publications, location history and basic information on uploaded files, alongside information on relevant actions, instruments, users and locations. The uploaded files themselves are in subdirectories, grouped by the individual objects, so the first file uploaded for an object with the ID #42 would be located in the directory ``objects/42/files/0`` under its original upload name. It also contains a ``.rdf`` file containing Dublin Core metadata in RDF/XML format for each object.
.. _citations:

Citation Guide
==============

Code
----

To cite the code of SampleDB, please use:

- Florian Rhiem et al., SampleDB. [Online; accessed <today>]. 2017-2021. DOI: 10.5281/zenodo.4012175. URL: https://github.com/sciapp/sampledb

Here is an example BibTeX entry:

.. code-block:: bibtex

    @Misc{SampleDB,
      author = {Florian Rhiem and others},
      title = {{SampleDB}},
      year = {2017--2021},
      url = "https://github.com/sciapp/sampledb",
      note = {[Online; accessed YYYY-MM-DD]},
      doi = {10.5281/zenodo.4012175}
    }

JOSS article
------------

To cite the `article about SampleDB published in the Journal of Open Source Software <https://doi.org/10.21105/joss.02107>`_, please use:

- Rhiem, F., (2021). SampleDB: A sample and measurement metadata database. Journal of Open Source Software, 6(58), 2107, https://doi.org/10.21105/joss.02107

.. code-block:: bibtex

    @article{Rhiem2021,
      doi = {10.21105/joss.02107},
      url = {https://doi.org/10.21105/joss.02107},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {58},
      pages = {2107},
      author = {Florian Rhiem},
      title = {SampleDB: A sample and measurement metadata database},
      journal = {Journal of Open Source Software}
    }

Objects
-------

If you wish to reference a specific sample, measurement or simulation in an instance of SampleDB, you can use its URL, e.g. ``https://iffsamples.fz-juelich.de/objects/0``. Bear in mind, however, that users need to be signed in even to view public objects.

An alternative to this would be to publish the information using :ref:`Dataverse <dataverse_export>` or uploading the information to a service such as `zenodo <https://zenodo.org/>`_.
Search
======

Users can find objects using the |service_name| search. There are two modes for the object search:

- A *simple* text-based search, and
- an *advanced* search using property comparisons

Simple Search
-------------

To use the simple search, users can enter words or phrases into the search field and will find all objects containing these.


.. _advanced_search:

Advanced Search
---------------

The advanced search allows a more fine grained search by performing comparisons on objects' properties and supporting Boolean algebra. Users can enter a query into the search field and select 'Advanced Search' in the adjacent dropdown, though using operators like ``=`` will automatically enable the advanced search mode. Another way to perform an advanced search is by clicking on the search icon next to an object's property that will search for other objects having the same value.

.. figure:: ../static/img/generated/advanced_search_by_property.png
    :alt: Search button for finding objects with an equal property

    Search button for finding objects with an equal property

When an advanced search is used, |service_name| will show the search tree that the query has been parsed into, to clearly show which filters were used.

.. figure:: ../static/img/generated/advanced_search_visualization.png
    :alt: Advanced search tree for the query '"Sb" in substance and (temperature < 110degC or temperature > 120degC)'

    Advanced search tree for the query ``"Sb" in substance and (temperature < 110degC or temperature > 120degC)``

Property Names
``````````````

To search for objects which have property fulfilling some condition, the internal name of that property must be known. Property names are set in an action's schema and the easiest way to find a property is to use the search button shown above. When searching for properties inside an object or array, dots are used as separators and a question mark can be used as wildcard for an array index, e.g. ``layers.?.name == "Base Layer"``.

Boolean operators
`````````````````

Boolean properties or comparisons of other properties can be combined with the boolean operators ``and``, ``or`` and ``not``. ``not`` is first in the order of operations, followed by ``and`` and ``or``. For a different order, parentheses can be used.

Text comparisons
````````````````

Text properties can currently be used either for direct comparisons, e.g. ``name = "Sample"``, or by checking whether a property contains a text, e.g. ``"MBE" in name``.

Quantity comparisons
````````````````````

Quantity properties can be compared using the basic mathematical comparison operators ``<``, ``<=``, ``>``, ``>=``, ``=`` and ``!=``. Comparisons will be performed in the quantities' base units.

Date comparisons
````````````````

Datetime properties can be compared using the basic mathematical comparison operators or the operators ``before``, ``after`` and ``on``. Dates to compare with must be entered using the format ``YYYY-MM-DD``.

Tag search
``````````

Objects with certain tags can be found using ``#`` and the name of the tag, e.g. ``#hbs``.
.. _instruments:

Instruments
===========

Instruments in |service_name| map real instruments to :ref:`actions` performed with them.

You can view a list of instruments at |service_instruments_url|. To make navigating the growing list of instruments easier, users can select **favorites** by clicking the star next to an instrument's name.

.. note::
    At this time, instruments can only be created by the |service_name| administrators. If you would like your instrument or action to be included, please `let us know`_.

.. _instrument_scientists:

Instrument Scientists
---------------------

Most instruments are associated with one or more **instrument scientists**, who will automatically be granted permissions for all objects created with their instruments. For more information on how permissions are handled, see :ref:`permissions`.

.. _instrument_scientist_notes:

Instrument Scientist Notes
^^^^^^^^^^^^^^^^^^^^^^^^^^

Instrument Scientists can save notes, e.g. for internal information on  instrument maintenance or configuration.
These notes will not be shown to regular users, though they should not be used to store sensitive information such as passwords.

.. _instrument_log:

Instrument Log
--------------

Instrument scientists can use the instrument log to keep track of problems, maintenance or other events. They can decide whether the log can be seen by other users, and whether other users can also create log entries, e.g. to report issues. When a new log entry is created, a notification will be sent to the instrument scientists.
.. _metadata:

Metadata and Schemas
====================

A schema specifies what metadata should be recorded when performing an action. The metadata should contain all the information required to reproduce what a user did.

Simple schemas can be created and edited using the graphical schema editor:

.. figure:: ../static/img/generated/schema_editor.png
    :alt: Editing an action using the schema editor

    Editing an action using the schema editor

To access more advanced features like arrays and conditions, users can instead edit the schema directly in a text-based editor by setting the *Use Schema Editor* setting to *No* in their preferences.

.. figure:: ../static/img/generated/disable_schema_editor.png
    :alt: Disabling the graphical schema editor

    Disabling the graphical schema editor

Schemas are defined using `JSON <https://www.json.org/>`_. They consist of nested JSON objects, with each object having a ``title`` and a ``type`` attribute, and more attributes depending on the ``type`` (note: when discussing JSON, what is called an "attribute" here is usually called a "property", this would be ambiguous however, as the metadata fields of an object are also called properties in SampleDB). The root object of each schema must have the ``object`` type and contain at least a text attribute called ``name``, as shown here:

.. code-block:: json
    :caption: A minimal schema containing only the object name

    {
      "type": "object",
      "title": "Object Information",
      "properties": {
        "name": {
          "title": "Name",
          "type": "text"
        }
      },
      "propertyOrder": ["name"],
      "required": ["name"]
    }

Data types
----------

Currently, the following basic data types are supported for metadata:

- Texts
- Booleans
- Quantities
- Datetimes

These can be used to form the following composite data types:

- Arrays
- Objects

Additionally, there are special data types:

- :ref:`Tags <tags>`
- :ref:`Hazards <hazards>`
- *plotly* Charts
- User References
- Object References
    - Sample References
    - Measurement References
    - Generic Object References
- Schema Templates

In the following, each data type and the attributes in a schema of each type are listed.

Objects
```````

Objects represent complex composite data types containing named properties. All SampleDB schemas start with an object (the root object), which consists of various other properties, mapping the name of each property to its subschema. In the minimal example above, the root object contains only a name, but you can add many more properties, as long as they all have a unique name.

.. code-block:: json
    :caption: A minimal schema containing only the object name

    {
      "type": "object",
      "title": "Example Object Information",
      "properties": {
        "name": {
          "title": "Name",
          "type": "text"
        }
        "created": {
          "title": "Creation Datetime",
          "type": "datetime"
        }
      },
      "propertyOrder": ["name", "created"],
      "required": ["name"]
    }

Object instances are JSON objects mapping the property names to the property data.

.. code-block:: json
    :caption: An object instance for the schema above.

    {
      "name": {
        "_type": "text",
        "text": "Demo Object"
      }
      "created": {
        "_type": "datetime",
        "utc_datetime": "2021-07-22 01:23:45"
      }
    }


type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``object``.

title
^^^^^

The title for the object as a JSON string or object, e.g. ``"Sample Information"`` or ``{"en": "Simulation Parameters"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given object property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

properties
^^^^^^^^^^

This JSON object maps names to the subschemas for other properties. The names may consist of lowercase characters, digits and underscores, but must not begin or end with an underscore. These names are, for example, used for the advanced search.

.. code-block:: json
    :caption: The properties attribute from the example above.

    "properties": {
      "name": {
        "title": "Name",
        "type": "text"
      }
      "created": {
        "title": "Creation Datetime",
        "type": "datetime"
      }
    }

.. note:: As mentioned above, the root object must have a required property called ``name`` with the type ``"text"``. This is the object name used on SampleDB to represent the object. Even though it is not process-specific, it might have process-specific restrictions, which is why it needs to be included in the schema.

.. note:: Try to use consistent property names between schemas, as this can greatly simplify searches, automated data entry or data analysis.

propertyOrder
^^^^^^^^^^^^^

As the ``properties`` JSON objects does not necessarily preserve the order of properties when processed by SampleDB, this attribute can set the desired order of properties when creating or displaying objects created with this schema. It is optional, though recommended to ensure consistent behavior. The property names are given as JSON strings in a JSON array, e.g. ``["name", "created"]``.

required
^^^^^^^^

This JSON array lists the names of all properties which must be set for an object to be valid, e.g. ``["name"]`` if only the ``name`` property must be set.

.. note:: For the root object, the ``name`` property must be required. If a ``hazards`` property exists in the root object, it must also be required.

.. note:: Sometimes, the behavior of required properties of type ``text`` may appear confusing, as even an empty text (``""``) is technically a text. If you want a text property to be non-empty, you can set a ``minLength`` for it in addition to setting it as required. See the text data type below for more information.

default
^^^^^^^

An object instance may be provided as the ``default`` attribute, e.g. for creating a new object. This should be a JSON object mapping each property name to its default data. The default must be a valid instance of the object schema, so the properties in it must fulfill all restrictions from their individual subschemas.

.. code-block:: json
    :caption: A default attribute for the example above.

    "default": {
      "name": {
        "_type": "text",
        "text": "Demo Object"
      }
      "created": {
        "_type": "datetime",
        "utc_datetime": "2021-07-22 01:23:45"
      }
    }

displayProperties
^^^^^^^^^^^^^^^^^

This attribute can be set to a JSON array containing the names of properties that should be displayed in a list of objects for the action this schema belongs to.

.. note:: This attribute may only be set for the root object.

.. note:: For some data types, it may be impossible to display them in the table, e.g. due to size restrictions. If you encounter issues with a property that should be possible to display but isn't shown correctly, you can `report it on GitHub <https://github.com/sciapp/sampledb/issues/new>`_.


batch
^^^^^

This attribute is a boolean that sets whether or not the action for this root object should create a batch of objects. If set to ``true``, the user will be able to enter how many objects should be created during object creation, and that number of objects will be created with identical data except for the name. By default, it is set to ``false``.

.. note:: This attribute may only be set for the root object.

batch_name_format
^^^^^^^^^^^^^^^^^

If the ``batch`` attribute is set to ``true``, this string attribute sets the format for the suffix that will be attached to the name of the objects created as a batch. It must follow the Python string format syntax and will be provided with the index of the individual object in that batch (starting with 1). If no ``batch_name_format`` is provided, the index will be used by itself. If the user set the name as ``Demo`` and were to create three items in a batch, then the default would result in the names ``Demo1``, ``Demo2`` and ``Demo3``, while a ``batch_name_format`` set to ``"-{:03d}"`` would result in the names ``Demo-001``, ``Demo-002`` and ``Demo-003``.

.. note:: This attribute may only be set for the root object.

notebookTemplates
^^^^^^^^^^^^^^^^^

A JSON array containing information about Jupyter notebook templates. For more information, see :ref:`jupyterhub_support`.

.. note:: This attribute may only be set for the root object.

Arrays
``````

Array properties represent a list of properties of the same type called ``items``. While each property in an object must have an individual subschema, all items in an array share their subschema.

.. code-block:: json
    :caption: An array property containing texts, with a default and length restrictions

    {
      "title": "Notes",
      "type": "array",
      "items": {
        "title": "Note",
        "type": "text"
      },
      "minItems": 1,
      "maxItems": 10,
      "default": [
        {
          "_type": "text",
          "text": "First default note"
        },
        {
          "_type": "text",
          "text": "Second default note"
        }
      ]
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``array``.

title
^^^^^

The title for the array as a JSON string or object, e.g. ``"Preparation Steps"`` or ``{"en": "Notes"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given array property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

items
^^^^^

This JSON object contains the subschema for the items in this array. Arrays may contain all other data types (aside from the special types ``tags`` and ``hazards``, which may only occur in the root object).

.. code-block:: json
    :caption: The items attribute from the example above.

    "items": {
      "title": "Note",
      "type": "text"
    }

minItems
^^^^^^^^

A number that sets how many items must at least be present in the array for it to be valid, e.g. ``1``. By default, there is no minimum number of items.

maxItems
^^^^^^^^

A number that sets how many items must at most be present in the array for it to be valid, e.g. ``10``. By default, there is no limit to the number of items.

default
^^^^^^^

A JSON array containing the default data for the individual items. See also the ``defaultItems`` attribute below.

.. code-block:: json
    :caption: The default attribute from the example above.

    "default": [
      {
        "_type": "text",
        "text": "First default note"
      },
      {
        "_type": "text",
        "text": "Second default note"
      }
    ]

defaultItems
^^^^^^^^^^^^

If the ``default`` attribute is not set, this number can be used to set how many items should be present by default, e.g. if it is common to have at least one item, but this is not a strict requirement, ``defaultItems`` could be set to ``1``.

style
^^^^^

This attribute is a string indicating how the array should be displayed. By default, the items will be shown one after another, but sometimes a different behavior may be desired. If the items are objects, using the ``table`` style may be useful to create a table with the items as rows and their properties in the columns. Alternatively, if the items should be in the columns and their properties should be in the rows, the ``horizontal_table`` style can be used. If the items are neither objects nor arrays, the ``list`` style may be useful to create a simple list.

.. note:: Using a style other than the default may lead to issues when entering or viewing object data. Please test the action and how its objects are displayed. If you encounter issues with a style, you can `report it on GitHub <https://github.com/sciapp/sampledb/issues/new>`_.

Texts
`````

Text properties represent various types of textual data:

- Single line texts
- Multi line texts
- Rich text using Markdown
- A selection from a list of predefined texts (displayed as a dropdown field)


.. code-block:: json
    :caption: A sample name as a text property with a default, a pattern and length restrictions

    {
      "title": "Sample Name",
      "type": "text",
      "minLength": 1,
      "maxLength": 100,
      "default": "Sample-",
      "pattern": "^.+$"
    }

.. code-block:: json
    :caption: A sample description allowing multiple lines of text

    {
      "title": "Description",
      "type": "text",
      "multiline": true
    }

.. code-block:: json
    :caption: A sample description allowing Markdown content

    {
      "title": "Description",
      "type": "text",
      "markdown": true
    }

.. code-block:: json
    :caption: A measurement option using predefined choices

    {
      "title": "Measurement Option",
      "type": "text",
      "choices": [
        "Option A",
        "Option B"
      ]
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``text``.

title
^^^^^

The title for the text as a JSON string or object, e.g. ``"Description"`` or ``{"en": "Substrate"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Please describe the process in detail."`` or ``{"en": "Can be filled in later."}``.

placeholder
^^^^^^^^^^^

The placeholder for the text when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Description"`` or ``{"en": "Substrate"}``.

default
^^^^^^^

The default value for this property, as a JSON string or object, e.g. ``"Example"`` or ``{"en": "Example"}``. If there are ``choices`` defined for this property, then the default must be one of the choices.


.. code-block:: json
    :caption: The default attribute from one of the examples above

    "default": "Sample-"

minLength
^^^^^^^^^

This attribute sets the minimum number of characters for the value of this property, e.g. ``1``. By default, there is no minimum length.

maxLength
^^^^^^^^^

This attribute sets the maximum number of characters for the value of this property, e.g. ``1``. By default, there is no maximum length.

pattern
^^^^^^^

A JSON string containing a `regular expression <https://docs.python.org/3/library/re.html#regular-expression-syntax>`_ limiting what values are valid for this property, e.g. ``^Sample-[0-9]{4}$`` to ensure only values starting with ``Sample-`` followed by a four digit number will be valid.

languages
^^^^^^^^^

Either a JSON array containing the allowed language codes for this property, e.g. ``["en", "de"]`` or the JSON string ``"all"`` to allow all languages enabled for user input. By default, this attribute is set to ``["en"]`` only allowing english language input.

choices
^^^^^^^

A JSON array of acceptable values, either as JSON objects or JSON strings. If choices are provided, the value for this property must be one of the choices and a dropdown menu will be used to let the user select the choice. If this property is not required, not selecting a choice at all and therefore not providing a value for this property will also be valid.

.. code-block:: json
    :caption: The choices attribute from one of the examples above

    "choices": [
      "Option A",
      "Option B"
    ]

.. note:: For properties with ``choices`` set, you cannot provide a ``placeholder`` value and should not set a ``minLength``, ``maxLength`` or ``pattern``. Setting ``choices``, ``multiline`` and ``markdown`` are all mutually exclusive.

multiline
^^^^^^^^^

This attribute is a boolean that sets whether or not the value of this property may contain multiple lines. By default, this is ``false`` and the field when creating or editing an object using this schema will be for a single line only.

.. note:: Setting ``choices``, ``multiline`` and ``markdown`` are all mutually exclusive.

markdown
^^^^^^^^

This attribute is a boolean that sets whether or not the value of this property should be rich text based on the Markdown syntax. If this attribute is true, users will be able to input multiple lines and use a Markdown editor to include formatting, images and other rich text elements in the value of this property. By default, this is ``false``.

.. note:: Setting ``choices``, ``multiline`` and ``markdown`` are all mutually exclusive.

Booleans
````````

Booleans represent a value that is either true or false.

.. code-block:: json
    :caption: A boolean property with a default

    {
      "title": "Lid Open?",
      "type": "bool",
      "default": true
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``bool``.

title
^^^^^

The title for the boolean as a JSON string or object, e.g. ``"Pressurization"`` or ``{"en": "Target Set"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Set if chamber was pressurized."`` or ``{"en": "Check box if a target was set"}``.

default
^^^^^^^

The default value for this property as a boolean, so ``true`` or ``false``.

Quantities
``````````

Properties of the ``quantity`` type represent physical quantities and unitless numbers. The ``units`` attribute is mandatory, so for unitless numbers it must be set to ``1``.

.. code-block:: json
    :caption: A temperature property with a default of 25°C (298.15K)

    {
      "title": "Temperature",
      "type": "quantity",
      "units": "degC",
      "default": 298.15
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``quantity``.

title
^^^^^

The title for the quantity as a JSON string or object, e.g. ``"Temperature"`` or ``{"en": "Detector Distance"}``.

placeholder
^^^^^^^^^^^

The placeholder for the text when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Temperature in K"`` or ``{"en": "Detector Distance (horizontal)"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Temperature in measurement chamber."`` or ``{"en": "Horizontal distance between sample and detector"}``.

default
^^^^^^^

The default value for this property as a number. This should be the value in base units, so if ``units`` is set to ``nm`` and you want to set a default of 10nm, you need to set ``default`` to ``0.00000001`` as it will be interpreted in meters.

units
^^^^^

A JSON string containing the units for this property, e.g. ``nm`` or ``degC``.

.. note:: These units will be parsed using the `pint Python Package <https://pint.readthedocs.io/en/latest/index.html>`_ with additional `units defined by SampleDB <https://github.com/sciapp/sampledb/blob/develop/sampledb/logic/unit_definitions.txt>`_.

Datetimes
`````````

Datetime properties represent points in time. They are stored using ``YYYY-MM-DD hh:mm:ss`` notation and UTC, though users may enter and display them in differing timezones.

.. code-block:: json
    :caption: A datetime property with a default value

    {
      "title": "Creation Datetime",
      "type": "datetime",
      "default": "2018-12-05 15:38:12"
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``datetime``.

title
^^^^^

The title for the datetime as a JSON string or object, e.g. ``"Creation Date"`` or ``{"en": "Use Before"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Use experiment starting time"`` or ``{"en": "Include cool down time in estimate"}``.

default
^^^^^^^

A default value for the property, as a JSON string using ``YYYY-MM-DD hh:mm:ss`` notation and UTC, e.g. ``"2021-07-23 08:00:00"``. If no default is given, the current date and time when creating or editing an object using this schema will be used as the default.

Tags
````

Tags are keywords that can be used to categorize and quickly find objects relating to a specific topic. They may only be used as a property of the root object with the name ``tags``. The tag values themselves may only consist of lowercase characters, digits and underscores.

.. code-block:: json
    :caption: A tags property with default tags

    {
      "title": "Tags",
      "type": "tags",
      "default": ["tag1", "tag2"]
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``tags``.

title
^^^^^

The title for the tags as a JSON string or object, e.g. ``"Tags"`` or ``{"en": "Keywords"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

default
^^^^^^^

A JSON array containing default tags as strings, e.g. ``[]`` or ``["demo", "documentation"]``. There must be no duplicates in the array and as noted above, tags are limited to lowercase characters, digits and underscores.

Hazards
```````

Hazards allow users to declare whether or not the substance represented by the object poses any hazards by selecting the relevant GHS pictograms. Hazards may only be used as a property of the root object with the name ``hazards``. If such a property exists, it must be required to avoid any ambiguity, so that users have to explicitly declare that a substance poses no hazards instead of just not entering any.

.. code-block:: json
    :caption: A hazards property

    {
      "title": "GHS hazards",
      "type": "hazards"
    }

type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``hazards``.

title
^^^^^

The title for the hazards as a JSON string or object, e.g. ``"GHS hazards"`` or ``{"en": "Hazards (GHS)"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

note
^^^^

A note to display below the hazards selection when creating or editing an object using this schema, as a JSON string or object, e.g. ``"See lab guidelines"`` or ``{"en": "Please provide additional information in the description."}``.

plotly Charts
`````````````

Properties of this type allow users to store JSON data that can be rendered by `plotly <https://plotly.com/>`_. This is most useful in combination with automated data entry as opposed to manually creating and entering the JSON data.

.. code-block:: json
    :caption: A plotly chart from the plotly documentation

    {
      "data": [
        {
          "x": [
            "giraffes",
            "orangutans",
            "monkeys"
          ],
          "y": [
            20,
            14,
            23
          ],
          "type": "bar"
        }
      ]
    }

For more information on the plotly JSON format, see the `JSON chart schema <https://plotly.com/chart-studio-help/json-chart-schema/>`_.

.. code-block:: json
    :caption: A plotly chart property

    {
      "title": "Temperature",
      "type": "plotly_chart"
    }


type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``plotly_chart``.

title
^^^^^

The title for the plotly chart as a JSON string or object, e.g. ``"Temperature"`` or ``{"en": "Z Distance Movement"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the JSON field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Will be filled by bot"`` or ``{"en": "Upload raw log file as well"}``.

User References
```````````````

Properties of this type allow you to reference SampleDB users.

.. code-block:: json
    :caption: A user reference property

    {
      "title": "Client",
      "type": "user"
    }


type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``user``.

title
^^^^^

The title for the property as a JSON string or object, e.g. ``"Client"`` or ``{"en": "Principal Investigator"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"For external users, leave blank and fill in information below"`` or ``{"en": "Remember to set as responsible user as well"}``.

default
^^^^^^^

A JSON number containing the user ID to be used as default selection, or a JSON string ``"self"`` to denote that the user who is currently creating or editing the object should be the default.

Object References
`````````````````

Properties of this type allow you to reference other objects, e.g. to denote a precursor material or a dataset used for a simulation. Using ``action_type_id`` or ``action_id`` you can limit which objects may be referenced using this property.

.. code-block:: json
    :caption: An object reference property

    {
      "title": "Measured Object",
      "type": "object_reference"
    }


type
^^^^

This sets the type for this subschema as a JSON string and must be set to ``object_reference``.

title
^^^^^

The title for the property as a JSON string or object, e.g. ``"Precursor"`` or ``{"en": "Calibration Measurement"}``.

may_copy
^^^^^^^^

This attribute is a boolean that sets whether or not the data for the given property may be copied when using the **Use as template** functionality in SampleDB. By default, it is set to ``true``.

dataverse_export
^^^^^^^^^^^^^^^^

This attribute is a boolean that controls whether this property should be exported as part of a :ref`dataverse_export` or not, although the exporting user will still have the choice to enable or disable this property during the export. By default, it is set to ``false``.

conditions
^^^^^^^^^^

This attribute is a JSON array containing a list of conditions which need to be fulfilled for this property to be available to the user. By default, no conditions need to be met. For examples and more information, see :ref:`conditions`.

note
^^^^

A note to display below the field when creating or editing an object using this schema, as a JSON string or object, e.g. ``"Leave blank if no precursor was used."`` or ``{"en": "Select the associated calibration measurement"}``.

action_type_id
^^^^^^^^^^^^^^

This attribute is a number that sets the ID of an action type to limit which actions an object referenced by this property may have been created with, e.g. ``-99`` to limit the property to samples.

action_id
^^^^^^^^^

This attribute is a number that sets the ID of an action to limit that only objects created with that action may be referenced by this property, e.g. ``1``.

Sample References
^^^^^^^^^^^^^^^^^

Properties of this type are a special case of object reference, limited to referencing samples. The same can be achieved using an object reference with ``action_type_id`` set to -99. These properties support the same attributes as those of type ``object_reference``, aside from ``action_id`` and ``action_type_id``. Their type must be ``sample``.

.. code-block:: json
    :caption: A sample reference property

    {
      "title": "Previous Sample",
      "type": "sample"
    }

Measurement References
^^^^^^^^^^^^^^^^^^^^^^

Properties of this type are a special case of object reference, limited to referencing measurements. The same can be achieved using an object reference with ``action_type_id`` set to -98. These properties support the same attributes as those of type ``object_reference``, aside from ``action_id`` and ``action_type_id``. Their type must be ``measurement``.

.. code-block:: json
    :caption: A measurement reference property

    {
      "title": "Preparatory Measurement",
      "type": "measurement"
    }

Schema Templates
```````````````

Schema Templates offer a way to easily reuse action schemas.

If an *action_type* is marked as includable into other actions it's possible to reuse the schema.

The schema for a template action could look like the following:

.. code-block:: json
   :caption: Minimal schema template

    {
      "title": "test",
      "type": "object",
      "properties": {
        "name": {
          "title": "Name",
          "type": "text"
        },
        "value": {
          "title": "Value",
          "type": "text"
        }
      },
      "required": [
        "name"
      ],
      "propertyOrder": [
        "name",
        "value"
      ]
    }

There is generally no difference to the schemas of other actions.

Schema templates can be included into other actions by providing a ``template`` for a property of type ``object``

.. code-block:: json
   :caption: Action with included schema template

    {
      "title": "Action with included Schema Template",
      "type": "object",
      "properties": {
        "name": {
          "title": "Name",
          "type": "text"
        },
        "included": {
          "title": "Included Schema Template",
          "type": "object",
          "template": 15
        }
      },
      "required": [
        "name"
      ],
      "propertyOrder": [
        "name",
        "included"
      ]
    }

Internally, this will then be treated as if the schema template were used for the property ``included`` there, except that the ``name`` property will be removed to avoid redundancies. The resulting action will be equivalent to:

.. code-block:: json
   :caption: Action with schema template

    {
      "title": "Action with included Schema Template",
      "type": "object",
      "properties": {
        "name": {
          "title": "Name",
          "type": "text"
        },
        "included": {
          "title": "Included Schema Template",
          "type": "object",
          "properties": {
            "value": {
              "title": "Value",
              "type": "text"
            }
          },
          "required": [],
          "propertyOrder": ["value"]
        }
      },
      "required": [
        "name"
      ],
      "propertyOrder": [
        "name",
        "included"
      ]
    }

When the schema template action is updated, all actions using it will be updated as well, as long as the resulting schema is still valid.

.. _conditions:

Conditional Properties
----------------------

Some properties might only sometimes be needed, based on some conditions, such as a particular setting of an instrument. Properties can contain conditions like this, consisting of a JSON object with a ``type`` and additional information depending on the type of condition.

.. code-block:: javascript
    :caption: A schema with a conditional property

    {
      "title": "Example Object",
      "type": "object",
      "properties": {
        "name": {
          "title": "Object Name",
          "type": "text",
          "languages": ["en", "de"]
        },
        "dropdown": {
          "title": "Dropdown",
          "type": "text",
          "choices": [
            {"en": "A"},
            {"en": "B"},
          ],
          "default": {"en": "A"}
        },
        "conditional_text": {
          "title": "Conditional Text",
          "type": "text",
          "markdown": true,
          "conditions": [
            {
              "type": "choice_equals",
              "property_name": "dropdown",
              "choice": {"en": "B"}
            }
          ]
        }
      },
      "required": ["name"]
    }

In the example schema above the property ``conditional_text`` will only be enabled if its ``choice_equals`` condition is fulfilled, that is if the ``dropdown`` property has the value ``{"en": "B"}`` selected.

The following types of conditions are supported by SampleDB:

choice_equals
`````````````

For this type of condition, the ``property_name`` attribute must be the name of another property, in the same object as the property the condition is for. The property of that name must be a property of type ``text`` with the ``choices`` attribute set. The condition must have a ``choice`` attribute that must be one of those choices, and for the condition to be fulfilled that choice must be selected.

.. code-block:: javascript
    :caption: A choice_equals condition

    {
      "type": "choice_equals",
      "property_name": "dropdown",
      "choice": {"en": "B"}
    }

user_equals
```````````

For this type of condition, the ``property_name`` attribute must be the name of another property, in the same object as the property the condition is for. The property of that name must be a property of type ``user``. The condition must have a ``user_id`` attribute that must be the ID of a user, and for the condition to be fulfilled that user must be selected.

.. code-block:: javascript
    :caption: A user_equals condition

    {
      "type": "user_equals",
      "property_name": "client",
      "user_id": 1
    }

If the ``user_id`` is set to ``null`` instead, the condition will be fulfilled if no user has been selected.

.. code-block:: javascript
    :caption: A user_equals condition for not having selected a user

    {
      "type": "user_equals",
      "property_name": "client",
      "user_id": null
    }

bool_equals
```````````

For this type of condition, the ``property_name`` attribute must be the name of another property, in the same object as the property the condition is for. The property of that name must be a property of type ``bool``. The condition must have a ``value`` attribute that must be either ``true`` or ``false``, and for the condition to be fulfilled the property must also be true or false, correspondingly.

.. code-block:: javascript
    :caption: A bool_equals condition

    {
      "type": "bool_equals",
      "property_name": "heating_on",
      "value": true
    }

object_equals
`````````````

For this type of condition, the ``property_name`` attribute must be the name of another property, in the same object as the property the condition is for. The property of that name must be a property of type ``object_reference``, ``sample`` or ``measurement``. The condition must have a ``object_id`` attribute that must be the ID of an object, and for the condition to be fulfilled that object must be selected.

.. code-block:: javascript
    :caption: An object_equals condition

    {
      "type": "object_equals",
      "property_name": "precursor",
      "object_id": 1
    }

If the ``object_id`` is set to ``null`` instead, the condition will be fulfilled if no user has been selected.

.. code-block:: javascript
    :caption: An object_equals condition for not having selected an object

    {
      "type": "object_equals",
      "property_name": "precursor",
      "object_id": null
    }

any / all
`````````

To denote that either only one or all of a list of conditions need to be fulfilled, the ``any`` or ``all`` condition type can be used, containing other conditions. An ``any`` condition is fulfilled, if any one of the conditions in it is fulfilled. If it does not contain any conditions, it will be considered as not being fulfilled. An ``all`` condition is fulfilled, if all of the conditions in it are fulfilled. If it does not contain any conditions, it will be considered as being fulfilled.

.. code-block:: javascript
    :caption: An any condition

    {
      "type": "any",
      "conditions": [
        {
          "type": "bool_equals",
          "property_name": "example_bool_1",
          "value": true
        },
        {
          "type": "bool_equals",
          "property_name": "example_bool_2",
          "value": true
        }
      ]
    }

.. code-block:: javascript
    :caption: An all condition

    {
      "type": "all",
      "conditions": [
        {
          "type": "bool_equals",
          "property_name": "example_bool_1",
          "value": true
        },
        {
          "type": "bool_equals",
          "property_name": "example_bool_2",
          "value": true
        }
      ]
    }

not
```

To denote that a certain condition must not be met, the ``not`` condition type can be used together with that other condition.

.. code-block:: javascript
    :caption: A not condition

    {
      "type": "not",
      "condition": {
        "type": "object_equals",
        "property_name": "example_object",
        "object_id": null
      }
    }

.. note:: If you need a new type of conditions, please `open an issue on GitHub <https://github.com/sciapp/sampledb/issues/new>`_ to let us know.
.. _objects:

Objects
=======

Objects in |service_name| represent **samples**, **measurements** or **simulations**. They store their metadata, allow users to track their life cycle and provide easy access to related objects, like measurements performed on a sample.

To create an object, users can select an :ref:`Action <actions>` and click on the button labeled **Create Sample**, **Perform Measurement** or **Perform Simulation**.

.. figure:: ../static/img/generated/action.png
    :alt: Sample Creation using an Example Action and Instrument

    Sample Creation using an Example Action and Instrument

Alternatively, an existing object can be used as a template for new objects, copying its metadata.

Upon its creation, each object in |service_name| is assigned a unique integral ID. While users can use text based names to quickly identify an object, only its ID will be used internally.

Metadata
--------

Each object has a set of metadata, specific to the way the object was created. A sample created using one instrument may have completely different relevant information than a measurement performed on it with another. As such, each :ref:`Action <actions>` has a so-called schema that defines what metadata should be stored for objects created with it.

In general, users will need to enter the information required to reproduce what they did. For more information on metadata, see :ref:`metadata`.

.. _tags:

Tags
````

Tags or Keywords serve as a means to organize objects using short, lower case words. A tag may consist of characters *a* to *z* (including *ä*, *ö*, *ü* and *ß*), digits *0* to *9*, *-* and *_*. When entering tags, previously used tags will be suggested to the users to encourage a common tag vocabulary without restricting users.


.. figure:: ../static/img/generated/tags.png
    :alt: Tags in the Object Metadata

    Tags in the Object Metadata

.. _hazards:

Hazards
```````

Samples which pose hazards to human health or the environment should be labelled using the `Globally Harmonized System of Classification and Labelling (GHS) <https://www.unece.org/trans/danger/publi/ghs/ghs_welcome_e.html>`_. For actions used to creating such samples, these hazards can be part of the object metadata and can be set during object creation or modification.


.. figure:: ../static/img/generated/hazards_input.png
    :alt: Input of GHS Hazards during Sample Creation

    Input of GHS Hazards during Sample Creation

Versioning
``````````

Once an object has been created, changes to its metadata will be versioned in |service_name|. Users can edit the object using the **Edit Object** button, and once they save their changes, a new version of the object will be stored internally. Although only the latest version is shown by default, users can view older versions from the :ref:`Activity Log <activity_log>`. If a change has been made erroneously, users can restore a previous version without any changes being lost.

Files
-----

Users can upload files related to the object, like sketches or notes generated by measuring software. As |service_name| only serves as a database for metadata, storage is limited and this function is **not** for uploading large datasets or measurement results. For those, a user could leave a :ref:`Comment <comments>` containing the files' path on an archival storage system, for example.

At this time, files can either be uploaded directly from the browser or from a smartphone or other mobile device by scanning a QR code. To upload a file click **Browse** and select either **Local Files** or **Smartphone**.

.. figure:: ../static/img/generated/files.png
    :alt: Files

    Files

File Information
````````````````

Users can view additional information on a file by clicking on the **i** icon on the right side of the file table. There they can edit a file's title and description, view its history and hide it.

Due to its nature as an archive, files uploaded to |service_name| cannot be deleted. If, however, a wrong file was uploaded accidentally or for some other reason a file's content should be hidden, clicking the **Hide** button will allow users to hide a file. It will not be deleted, but regular users will be unable to access it afterwards.

.. figure:: ../static/img/generated/file_information.png
    :alt: File Information

    File Information

.. _comments:

Comments
--------

Users can leave comments on objects as free form text, e.g. to share additional information that does not fit the predefined metadata fields. These comments are displayed chronologically on the object's page and will be included in data exports.

.. figure:: ../static/img/generated/comments.png
    :alt: Comments

    Comments

.. _activity_log:

Activity Log
------------

The activity log shows a timeline of the object's life cycle, containing events like its creation, file uploads and when it was used for another object.

.. figure:: ../static/img/generated/activity_log_dontblock.png
    :alt: Activity Log

    Activity Log

.. _locations:

Location
--------

To indicate where a sample is stored, a location and/or a responsible user can be assigned to it. When a user is assigned responsibility for an object, they can confirm this either on the object's page or using the :ref:`Notification <notifications>` they received for the assignment.

The location log shows where an object has been stored and when it was moved.

.. figure:: ../static/img/generated/locations.png
    :alt: Location

    Location

.. _permissions:

Permissions
-----------

By default, samples, measurements and simulations are visible only to the user who created them and to the instrument scientists of the instrument the objects were created with. Additionally, administrators of |service_name| have access to the database the information is stored in. Object permissions can be used to share access to these objects with other :ref:`users`, :ref:`groups` or :ref:`projects`.

The object permissions built into |service_name| fall into three categories:

- **Read**: The permission to **view objects** and their properties, files and comments.
- **Write**: The permission to **edit objects** and their properties and add files and comments.
- **Grant**: The permission to **grant permissions** to other users.

Each of these categories is built on top of the other, with **Write** permissions including **Read** permissions and **Grant** permissions including **Write** permissions.

.. figure:: ../static/img/generated/object_permissions.png
    :alt: Object Permissions

    Object Permissions

To modify the permissions of an object, any user with **Grant** permissions can click the **Edit permissions** button on the object's page. They can then view the existing permissions, modify them or add new permissions for users, basic groups or project groups.

Although administrators are shown to have **Grant** permissions for all objects, this only reflects their access to the database mentioned above. At this time, administrators do not automatically have **Grant** permissions for all objects.

.. _default_permissions:

Default Permissions
```````````````````

When an object is created, its creator, any associated instrument scientists and the administrators will have **Grant** permissions. They can then allow other users to access the data by granting them permissions. To make this more convenient, each user has a set of **default permissions** in the :ref:`preferences`, which will be applied to all objects they create in the future.

.. figure:: ../static/img/generated/default_permissions.png
    :alt: Default Permissions

    Default Permissions in the User Preferences

.. _pdf_export:

Data Export
-----------

Users can export object information to a PDF file, e.g. for printing or offline usage. Note that the exported object information will not be fully complete, e.g. only files of some formats will be included in the PDF and only the current metadata version will be shown.

Alternatively, users can export object information as an archive, which contains the full object information as a JSON file and all files uploaded for the object.

Along with the current object, related objects can be exported along with it, e.g. a sample can be exported together with all measurements performed with it.

Users can also export information for all objects which the user has **READ** permissions for (see :ref:`export`).

Labels
------

|service_name| can be used to create labels for newly created samples. These labels will contain the object's ID, name, creator and creation date, along with :ref:`hazards` if those were specified as part of the object's metadata.

.. note::
    If you require a label format that isn't covered by the ones generated at this time, please `let us know`_.

.. figure:: ../static/img/generated/labels.png
    :alt: Generated Labels

    Generated Labels
.. _projects:

Project Groups
==============

Project Groups are one of the two ways for organizing :ref:`users` in |service_name|. For :ref:`Object Permissions <permissions>`, a project group acts as a single entity with members of the group sharing the permissions granted to it up to the permissions they have for the project group itself.

Users can be members of any number of project groups, leave the projects they are in or invite other users if they have **GRANT** permissions for the project group. Members with **GRANT** permissions can also remove other members or delete the project group as a whole.

A project group can contain :ref:`groups` as well as users, though it cannot contain other project groups. Instead, a parent-child relationship can be established between project groups, optionally allowing users with **GRANT** permissions in a child project group to invite other users to both the child and parent project groups.

A project group may be linked to an object, which shows that the project group is meant for managing the permissions for that object or objects related to it.

Project groups may be best suited for those who want to model their object permissions in |service_name| according to an organizational hierarchy. For a simpler approach, see :ref:`groups`.

.. note::
    Administrators can disable creation of project groups by regular users. Administrators can also disable subproject support.

.. note::
    Project groups were formerly known only as projects. The name was changed to distinguish them from objects containing information about research projects.

.. list-table:: Comparing Basic and Project Groups
   :header-rows: 1

   * - Feature
     - Basic Groups
     - Project Groups
   * - Members
     - Users
     - Users and Basic Groups
   * - Hierarchy
     - Flat
     - Nested (Project Groups can have Child Project Groups)
   * - Object Permissions
     - All members get the Object Permissions granted to the Basic Group
     - Each Project Group member has a permissions level (**Read**/**Write**/**Grant**) for the Project Group and will get Object Permissions granted to the Project Group up to that level
   * - Member Invitation
     - Any member of a Basic Group can invite new members
     - Only members of a Project Group with **Grant** permissions can invite new members
.. _groups:

Basic Groups
============

Basic Groups are one of the two ways for organizing :ref:`users` in |service_name|. For :ref:`Object Permissions <permissions>`, a basic group acts as a single entity with all members of the group sharing any permissions granted to it.

Users can be members of any number of basic groups, leave the basic groups they are in or invite other users. Any member of a basic group can remove other members or delete the basic group as a whole.

As all users in a basic group are equal, basic groups may be best suited for those who want to freely share objects with their colleagues without fine-grained permissions. For a more detailed or hierarchical approach, see :ref:`projects`.

.. note::
    Administrators can disable creation of basic groups by regular users.

.. note::
    Basic Groups were formerly known only as groups. The name was changed along with the renaming of :ref:`projects`.

.. list-table:: Comparing Basic and Project Groups
   :header-rows: 1

   * - Feature
     - Basic Groups
     - Project Groups
   * - Members
     - Users
     - Users and Basic Groups
   * - Hierarchy
     - Flat
     - Nested (Project Groups can have Child Project Groups)
   * - Object Permissions
     - All members get the Object Permissions granted to the Basic Group
     - Each Project Group member has a permissions level (**Read**/**Write**/**Grant**) for the Project Group and will get Object Permissions granted to the Project Group up to that level
   * - Member Invitation
     - Any member of a Basic Group can invite new members
     - Only members of a Project Group with **Grant** permissions can invite new members
.. _actions:

Actions
=======

Processes like **creating a sample**, **performing a measurement** or **running a simulation** are represented as actions in |service_name| and whenever such an action is performed a new :ref:`Object <objects>` should be created in |service_name|.

Either generic or associated with an :ref:`Instrument <instruments>`, each action contains a name, a description and a :ref:`Schema <metadata>`.

You can view a list of actions at |service_actions_url|. Similar to instruments, users can select **favorites** by clicking the star next to an action's name.

Action Types
------------

The type of an Action describes the general kind of process it represents, such as **Sample Creation**. There are four built-in action types for creating samples, performing a measurement, running a simulation, and for use as schema templates. Administrators can create additional custom action types.
Schema templates have a special role, as they are usually not used like other actions. Instead, actions marked for use as a schema template can be included in the schemas of other actions.

Custom Actions
--------------

Users can create custom actions to represent their own processes or instruments that are not yet part of |service_name|. These actions can be private, only usable by the users who created them, or public, so all everyone can see and use them, or the user can grant permissions to specific users, basic groups or project groups.

To create a custom action, users can either use an existing action as a template or write a :ref:`Schema <metadata>` from scratch.

.. note::
    Custom Actions are an advanced feature that most users of |service_name| will not need. If you would like your instrument or action to be included without writing your own schema, please `let us know`_.

    .. only:: iffSamples

        If you would like to try working with custom actions, please `use the development and testing deployment of iffSamples <https://docker.iff.kfa-juelich.de/dev-sampledb/>`_.
.. _users:

Users
=====

.. _authentication:

Authentication
--------------

Users at facilities which use LDAP for user management can use their **LDAP username** with the corresponding password to sign in. An account will be created automatically.

Guests, or users at facilities without LDAP, should ask another user of |service_name| for an **invitation**, e.g. the scientist responsible for the instrument they will be using. Once they have confirmed their email address by clicking the confirmation link in the invitation email, they can then set a password for their new |service_name| account.

.. figure:: ../static/img/generated/guest_invitation.png
    :alt: User Invitation Form

    User Invitation Form

Users can find the invitation form at |service_invitation_url|.

.. _preferences:

Preferences
-----------

Users can edit their preferences by clicking on their name in the top right corner and selecting preferences.

The preferences are split into the following sections:

User Information
````````````````

Users can update their user name displayed on |service_name|, e.g. in the event of a marriage. They can also change their email address, which will be updated once the new address has been confirmed, and set a publicly visible ORCID iD, affiliation and role.

Authentication Methods
``````````````````````

Users can have multiple ways of signing in to |service_name|, for example using their LDAP account or using an email address. This section of the user preferences can be used to add, modify or remove such authentication methods, e.g. for users leaving their institute but still requiring access to their sample data.

API Tokens
``````````

API tokens are meant as an alternative method for authentication when using the :ref:`HTTP API <http_api>` for individual scripts and allow you to monitor the requests made with the token. When you create a new API token, it will be shown to you once and you will be asked to save it. If you lose access to a token, simply delete it and create a new one as replacement.

While API tokens pose less of a risk than using your own user credentials in a script, please keep your API tokens private. Do not commit them with a version control system like git and do not share them with others!

Default Permissions
```````````````````

To automatically set permissions for future objects, users can set **default permissions** in their preferences. These will be applied whenever an object like a sample or measurement is created afterwards.

For more information, see :ref:`default_permissions`.

.. _notifications:

Notifications
-------------

Users will receive notifications whenever they need to be informed about an activity on |service_name|. Whenever a user has unread notifications, a bell with the number of unread notifications is shown in the navigation bar.

.. figure:: ../static/img/generated/unread_notification_icon.png
    :alt: Unread Notification Icon

    Unread Notification Icon

Bot Users
---------

Tasks like object creation can be automated by using the :ref:`HTTP API <http_api>`. When this is done in connection to an instrument or a facility instead of an individual user, it may be better to create a dedicated user account solely for this purpose.

Readonly Users
--------------

Users can be limited to READ permissions, e.g. for former employees who should still have access to their data but should not be able to create new SampleDB entries.

Hidden users
------------

Users can also be hidden from users lists, which may be useful in similar use cases as when marking a user as readonly. These users can still be seen as part of an object's history or as members of basic and project groups, but they will not be shown in the central users list, when granting permissions, inviting a user to a group, etc.

Deactivated users
------------------

Users can also be deactivated. These users will be unable to sign in to their account or use the API until they have been reactivated by an administrator. As they will be unable to access their own data, this should only be used if marking a user as readonly will not suffice.
.. _tls_termination:

TLS Termination
===============

The SampleDB container comes with a built-in, production ready CherryPy web server which serves SampleDB via HTTP on port 8000. This is great for local testing, but if you run SampleDB in production, you should use HTTPS instead.

The recommended way of doing this is to set up an nginx container as a so-called *reverse proxy* that will handle TLS termination. Clients will send requests to nginx using HTTPS and nginx will internally communicate to the SampleDB container using HTTP.

You may have existing infrastructure which you can integrate SampleDB into, this guide is merely an example in case you are starting on a completely fresh system.

Please follow general system administration best practices.

Requirements
------------

To follow this guide and set up TLS termination, you will need:

- a domain name that is mapped to your host
- a TLS certificate, saved as ``certificate.crt``
- the corresponding private key, saved as ``certificate.key``

Your facility is likely to already have a domain name, and you may be able to use one of its subdomains, e.g. ``sampledb.yourdomain.tld``. If so, there likely is an established certficate authority your facility uses. If not, consider using a certificate authority like https://letsencrypt.org/.

Docker network
--------------

The first step is to make sure the three containers run in their own Docker network, so that they can communicate with each other while port 8000 of the SampleDB container will only be reachable internally:

.. code-block:: bash

    docker network create --driver=bridge --subnet=172.24.0.0/16 sampledb-network

Next, assign IP addresses to the individual containers. The rest of this guide assumes that you use:

=================  ==========
container          IP address
=================  ==========
sampledb           172.24.0.1
sampledb-postgres  172.24.0.2
sampledb-nginx     172.24.0.3
=================  ==========

When creating the containers with ``docker run``, pass the network and IP address to docker, by using the option:

.. code-block:: bash

    --network=sampledb-network --ip=<ip_address>

Also make sure to remove the ``-p 8000:8000`` option from the SampleDB container.

nginx Configuration
-------------------

You can then create a ``sampledb.conf`` file for nginx:

.. code-block::

    upstream sampledb {
       server 172.24.0.1:8000;
    }

    ssl_session_cache shared:ssl_session_cache:10m;

    server {
        listen              80;
        server_name         <domain_name>;
        return              301 https://$server_name$request_uri;
    }
    server {
        listen              443 ssl;
        server_name         <domain_name>;
        ssl_certificate     certificate.crt;
        ssl_certificate_key certificate.key;
        ssl_protocols       TLSv1.2;
        ssl_ciphers         "EECDH+AESGCM:EDH+AESGCM:!aNULL";
        ssl_prefer_server_ciphers on;
        ssl_dhparam         dhparam.pem;

        add_header Strict-Transport-Security max-age=63072000;
        add_header X-Frame-Options DENY;
        add_header X-Content-Type-Options nosniff;

        client_max_body_size 64M;

        location / {
            proxy_redirect     off;
            proxy_set_header   Host $host;
            proxy_set_header   X-Real-IP $remote_addr;
            proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header   X-Forwarded-Proto $scheme;
            proxy_pass         http://sampledb;
        }
    }

Note the usage of ``<domain_name>`` and fill in your domain name there.

dhparam.pem
-----------

This configuration file requires the existance of a file ``dhparam.pem`` containing the parameters for the Diffie-Helmann key exchange. To generate this file, run:

.. code-block:: bash

    openssl dhparam -out dhparam.pem 4096

This may take a very long time to run.

nginx Container
---------------

With these files in place, you can start the nginx container:

.. code-block::

    docker run \
        -d \
        -v `pwd`/certificate.crt:/etc/nginx/certificate.crt:ro \
        -v `pwd`/certificate.key:/etc/nginx/certificate.key:ro \
        -v `pwd`/dhparam.pem:/etc/nginx/dhparam.pem:ro \
        -v `pwd`/sampledb.conf:/etc/nginx/conf/default.conf:ro \
        --network=sampledb-network \
        --ip=172.24.0.3 \
        --restart=always \
        --name sampledb-nginx \
        -p 80:80 \
        -p 443:443 \
        nginx

You should now be able to access SampleDB using your domain name and HTTPS.
.. _jupyterhub_support:

JupyterHub Support
==================

JupyterHub is a web application that can provide Jupyter notebook servers to multiple users. This can allow your users to run Jupyter notebooks on one of your systems, giving them access to their data and a software environment configured for your facility. For general information on JupyterHub, see `the JupyterHub documentation <https://jupyterhub.readthedocs.io/>`_

SampleDB can interface with JupyterHub with the help of a notebook templating server. This way metadata stored for a sample or measurement in SampleDB can be inserted into a Jupyter notebook, which can then be run on the JupyterHub, allowing your users to start a data analysis with a few clicks.

If you do not wish to run a JupyterHub instance, you can still run the notebook templating server and configure it so that the generated notebooks can be downloaded by the users. In that case, you can skip the following step.

JupyterHub
----------

You can find information on how to set up JupyterHub in the `Get Started guide <https://jupyterhub.readthedocs.io/en/stable/getting-started/index.html>`_ for JupyterHub. During this setup, you will pick an ``Authenticator`` and a ``Spawner``. Depending on the ``Spawner``, you will also have to configure how persistent storage is handled for the JupyterHub. The rest of this guide assumes that there is some form of persistent storage, which the notebook templating server can write notebooks to.

Notebook Templating Server
--------------------------

Depending on your specific configuration of JupyterHub, you will need to write and set up a notebook templating server. SampleDB will redirect the user to a URL containing the name of the desired template, e.g.: ``https://jupyterhub.example.com/templates/t/some_instrument/data_analysis.ipynb`` In this example:

- ``https://jupyterhub.example.com/templates/`` is the base URL of a notebook templating server,
- ``/t/`` is a separator, and
- ``/some_instrument/data_analysis.ipynb`` is the name of the template.

The selected metadata will be sent as an urlencoded form parameter in the request body, e.g. ``{"measurement_name": "Experiment-12345"}`` would be sent as ``params=%7B%measurement_name%22%3A+%22Experiment-12345%22%7D``.

The templating server will check if the requested template exists and if so, ask the user for confirmation and will insert the properties in ``params`` into the notebook as a new cell. The resulting notebook is then stored on the persistent storage for your JupyterHub and the user is redirected to it.

If the name of a template contains a parameter name in braces, e.g. ``Data_Analysis_{measurement_name}.ipynb``, then the value of that parameter is inserted into the name when generating the notebook.

For more information and a Flask blueprint for writing a server for your JupyterHub instance, see `the notebook_templates project page <https://github.com/sciapp/notebook_templates>`_.

Notebook Templates
------------------

The notebook templates themselves are simply regular Jupyter notebooks, which expect some variables to be defined in the cell inserted by the templating server. A template might, for example, rely on the path to a data file being passed to it as a parameter and then open that file for analysis.

SampleDB configuration
----------------------

To configure SampleDB to use your JupyterHub and the notebook templating server, you should set the SAMPLEDB_JUPYTERHUB_NAME, SAMPLEDB_JUPYTERHUB_URL and SAMPLEDB_JUPYTERHUB_TEMPLATES_URL :ref:`configuration values <jupyterhub_configuration>`. For actions with a notebook template, SampleDB will then add a button for creating a Jupyter notebook.

Actions
-------

Notebook templates are defined for each action, as they depend on the type of metadata that's available. Notebook templates are added to the schema as the top-level property ``notebookTemplates``, which contains a list of objects.

Each of the objects represents one notebook template, with a ``title``, its relative ``url`` and a its ``params``. Each template parameter in ``params`` is mapped to the path to the actual data in the schema as a list of keys, as shown in the following example:

.. code-block:: javascript

    "notebookTemplates": [
      {
        "title": "Data Analysis",
        "url": "some_instrument/data_analysis.ipynb",
        "params": {
          "measurement_name": [
            "name",
            "text"
          ],
          "created": [
            "created",
            "utc_datetime"
          ]
        }
      }
    ]

Here the template parameter ``measurement_name`` will be filled with the ``text`` of the ``name`` property of an existing object, and the parameter ``created`` will contain the value of the ``created`` property as a string containing the UTC datetime.
.. _configuration:

Configuration
=============

When running a SampleDB installation, you can set environment variables to configure it. The following sections will describe groups of such variables and their effects.

E-Mail
------

.. list-table:: E-Mail Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_CONTACT_EMAIL
     - An email address for users to contact
   * - SAMPLEDB_MAIL_SENDER
     - The email address used for outbound emails
   * - SAMPLEDB_MAIL_SERVER
     - The mail server used for outbound emails
   * - SAMPLEDB_MAIL_PORT
     - The port to use for connections to the mail server (default: ``25``)
   * - SAMPLEDB_MAIL_USE_TLS / SAMPLEDB_MAIL_USE_SSL
     - Whether to use TLS or SSL for connections to the mail server (default: ``False``)
   * - SAMPLEDB_MAIL_USERNAME
     - The username sent to the mail server
   * - SAMPLEDB_MAIL_PASSWORD
     - The password sent to the mail server

While the ``SAMPLEDB_CONTACT_EMAIL``, ``SAMPLEDB_MAIL_SENDER`` and ``SAMPLEDB_MAIL_SERVER`` variables are required, you may need to set one or more of the other variables to connect to your mail server, depending on its configuration.

.. _ldap_configuration:

LDAP
----

.. list-table:: LDAP Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_LDAP_NAME
     - The name of the LDAP server shown to the users
   * - SAMPLEDB_LDAP_SERVER
     - The ldaps-URL of the LDAP server
   * - SAMPLEDB_LDAP_USER_BASE_DN
     - The LDAP base DN to search users with
   * - SAMPLEDB_LDAP_UID_FILTER
     - The filter to use for identifying a user in LDAP as python template, e.g. ``(uid={})``
   * - SAMPLEDB_LDAP_NAME_ATTRIBUTE
     - The name of the attribute containing a user's name in LDAP, e.g. ``cn``
   * - SAMPLEDB_LDAP_MAIL_ATTRIBUTE
     - The name of the attribute containing a user's email address in LDAP, e.g. ``mail``
   * - SAMPLEDB_LDAP_OBJECT_DEF
     - The object def to use for looking up user attributes, e.g. inetOrgPerson
   * - SAMPLEDB_LDAP_USER_DN
     - The DN of an LDAP user to use when searching for other users
   * - SAMPLEDB_LDAP_PASSWORD
     - The password for the user identified by SAMPLEDB_LDAP_USER_DN
   * - SAMPLEDB_TESTING_LDAP_LOGIN
     - The uid of an LDAP user (only used during tests)
   * - SAMPLEDB_TESTING_LDAP_PW
     - The password for the ldap user identified by SAMPLEDB_TESTING_LDAP_LOGIN (only used during tests)

If you use LDAP for user management, you can use these variables to configure how SampleDB should connect to your LDAP server.


Customization
-------------

.. list-table:: Customization Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_SERVICE_NAME
     - The name of the service
   * - SAMPLEDB_SERVICE_DESCRIPTION
     - A short description of the service
   * - SAMPLEDB_SERVICE_IMPRINT
     - The URL to use for the imprint link
   * - SAMPLEDB_SERVICE_PRIVACY_POLICY
     - The URL to use for the privacy policy link
   * - SAMPLEDB_PDFEXPORT_LOGO_URL
     - A file, http or https URL for a PNG or JPEG logo to be included in object export PDF documents
   * - SAMPLEDB_PDFEXPORT_LOGO_ALIGNMENT
     - The alignment (left, center or right) of the logo, if SAMPLEDB_PDFEXPORT_LOGO_URL is set (default: right)
   * - SAMPLEDB_HELP_URL
     - The URL to use for the help link

You can use these variables to customize how your SampleDB instance is called, described and which links are included in the footer. The logo at the given PDFEXPORT_LOGO_URL will be fetched when SampleDB is started and cached afterwards. To refresh the logo, you will need to restart SampleDB.

.. _jupyterhub_configuration:

JupyterHub Support
------------------

.. list-table:: JupyterHub Support Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_JUPYTERHUB_NAME
     - The name of your JupyterHub server (default: ``JupyterHub``)
   * - SAMPLEDB_JUPYTERHUB_URL
     - The base URL of your JupyterHub server
   * - SAMPLEDB_JUPYTERHUB_TEMPLATES_URL
     - The URL of a Jupyter notebook templating server (default: SAMPLEDB_JUPYTERHUB_URL + ``/templates``, if SAMPLEDB_JUPYTERHUB_URL is set)

For more information on JupyterHub support and Jupyter notebook templates, see :ref:`jupyterhub_support`.

.. _dataverse_configuration:

Dataverse Export
----------------

.. list-table:: Dataverse Export Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_DATAVERSE_NAME
     - The name of the Dataverse server (default: ``Dataverse``)
   * - SAMPLEDB_DATAVERSE_URL
     - The base URL of the Dataverse server
   * - SAMPLEDB_DATAVERSE_ROOT_IDS
     - A comma seperated list of IDs of Dataverses, which objects may be exported to  (default: ``:root``)

For more information on the Dataverse export, see :ref:`dataverse_export`.

Administrator Account
---------------------

.. list-table:: Administrator Account Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_ADMIN_PASSWORD
     - The password for the admin account.
   * - SAMPLEDB_ADMIN_USERNAME
     - The username for the admin account (default: ``admin``)
   * - SAMPLEDB_ADMIN_EMAIL
     - See email address for the admin account (default: SAMPLEDB_CONTACT_EMAIL)


If no users exist yet and the ``SAMPLEDB_ADMIN_PASSWORD`` variable is set, a new user will be created with this password. This user will be a SampleDB admin. The username for this user will be set to value of ``SAMPLEDB_ADMIN_USERNAME`` and the email address for this user will be set to the value of ``SAMPLEDB_ADMIN_EMAIL``.

If another user already exists, these variables will have no effect. It is meant for creating an administrator account as part of the initial setup.

Miscellaneous
-------------

.. list-table:: Miscellaneous Configuration Environment Variables
   :header-rows: 1

   * - Variable Name
     - Description
   * - SAMPLEDB_FILE_STORAGE_PATH
     - A path to the directory that uploaded files should be stored in
   * - SAMPLEDB_SERVER_NAME
     - The server name for Flask. See: https://flask.palletsprojects.com/en/1.1.x/config/#SERVER_NAME
   * - SAMPLEDB_SQLALCHEMY_DATABASE_URI
     - The database URI for SQLAlchemy. See: https://flask-sqlalchemy.palletsprojects.com/en/2.x/config/
   * - SAMPLEDB_SECRET_KEY
     - The secret key for Flask and Flask extensions. See: https://flask.palletsprojects.com/en/1.1.x/config/#SECRET_KEY
   * - SAMPLEDB_WTF_CSRF_TIME_LIMIT
     - The time limit for WTForms CSRF tokens in seconds. See: https://flask-wtf.readthedocs.io/en/stable/config.html
   * - SAMPLEDB_INVITATION_TIME_LIMIT
     - The time limit for invitation links in seconds.
   * - SAMPLEDB_ONLY_ADMINS_CAN_MANAGE_LOCATIONS
     - If set, only administrators will be able to create and edit locations.
   * - SAMPLEDB_ONLY_ADMINS_CAN_CREATE_GROUPS
     - If set, only administrators will be able to create basic groups.
   * - SAMPLEDB_ONLY_ADMINS_CAN_DELETE_GROUPS
     - If set, only administrators will be able to delete non-empty basic groups.
   * - SAMPLEDB_ONLY_ADMINS_CAN_CREATE_PROJECTS
     - If set, only administrators will be able to create project groups.
   * - SAMPLEDB_LOAD_OBJECTS_IN_BACKGROUND
     - If set, object selections will be loaded in the background using AJAX.
   * - SAMPLEDB_DISABLE_USE_IN_MEASUREMENT
     - If set, the "Use in Measurement" button will not be shown.
   * - SAMPLEDB_DISABLE_SUBPROJECTS
     - If set, project groups cannot have child project groups assigned to them.
   * - SAMPLEDB_ENFORCE_SPLIT_NAMES
     - If set, force names to be entered as "surname, given names". **Note:** this will prevent users with a mononym from setting their name correctly!
   * - SAMPLEDB_PYBABEL_PATH
     - The path to the pybabel executable (default: ``pybabel``)
   * - SAMPLEDB_EXTRA_USER_FIELDS
     - A JSON-encoded dict containing extra user fields, e.g. ``{"phone": {"name": {"en": "Phone No."}, "placeholder": {"en": "Phone No."}}}`` (default: ``{}``)
   * - SAMPLEDB_SHOW_PREVIEW_WARNING
     - If set, a warning will be shown indicating that the instance is a preview installation and that data will be deleted.
   * - SAMPLEDB_DISABLE_INLINE_EDIT
     - If set, the inline edit mode will be disabled and users will not be able to edit individual fields.
   * - SAMPLEDB_SHOW_OBJECT_TITLE
     - If set, object schema titles will be shown when viewing metadata by default. Users may override this setting in their preferences.
   * - SAMPLEDB_HIDE_OBJECT_TYPE_AND_ID_ON_OBJECT_PAGE
     - If set, the object type and id, e.g. "Sample #4" will not be shown on the object page.
   * - SAMPLEDB_MAX_BATCH_SIZE
     - Maximum number of objects that can be created in one batch (default: 100)

There are other configuration values related to packages used by SampleDB. For more information on those, see the documentation of the corresponding packages.
.. _backup_and_restore:

Backup and Restore
==================

SampleDB stores all its information in:

- the PostgreSQL database, and
- the file directory

So to create a backup of SampleDB, you will need to create backups of these.

It is recommended that you stop the SampleDB container before creating backups and start it again afterwards.

While you yourself will need to decide when and how exactly you want to create backups, the following sections show examples of how backups of these two sources of information can be created.

Please follow general system administration best practices.

The PostgreSQL database
-----------------------

One way of creating a backup of a PostgreSQL database is to create an SQL dump using the `pg_dump` tool:

.. code-block:: bash

    docker exec sampledb-postgres pg_dump -U postgres postgres > backup.sql

The resulting ``backup.sql`` file can then be copied to another system.

To restore the PostgreSQL database from such an SQL dump, you should first remove the existing database:

.. code-block:: bash

    docker stop sampledb-postgres
    docker rm sampledb-postgres
    rm -rf pgdata

You can then recreate the database container and restore the backup using the ``psql`` tool:

.. code-block:: bash

    docker run \
        -d \
        -e POSTGRES_PASSWORD=password \
        -e PGDATA=/var/lib/postgresql/data/pgdata \
        -v `pwd`/pgdata:/var/lib/postgresql/data/pgdata:rw \
        --restart=always \
        --name sampledb-postgres \
        postgres:12
    docker exec -i sampledb-postgres psql -U postgres postgres < backup.sql

If you have set different options for the database container before, e.g. setting it in a specific network and giving it a fixed IP, you should also set these options here.

For more information on backing up a PostgreSQL database and restoring a backup, see the `PostgreSQL documentation on Backup and Restore <https://www.postgresql.org/docs/current/backup.html>`_

The file directory
------------------

If you have followed the Getting Started guide, you will have set the ``SAMPLEDB_FILE_STORAGE_PATH`` variable and mounted a local directory ``files`` into the SampleDB container so that SampleDB can store uploaded files in this directory.

Files are uploaded there over time and should not change, so this directory is especially suited for incremental backups.

One way to copy it to another system is to use the ``rsync`` tool:

.. code-block:: bash

    rsync -a files <hostname>:<backup_directory>

Or, if you wish to create a backup locally, you can simply use ``cp``:

.. code-block:: bash

    cp -an files <backup_directory>

Here the ``-n`` option prevents copying files which already exist in the backup directory.

To restore the backup, simply copy the backup to your local file directory.

.. note::

    By default, new files will be stored in the database instead of in the file directory, so the file directory may be empty... _setup:

Getting Started
===============

Step 1: Minimal Installation
----------------------------

We recommend using our pre-built Docker images for setting up ``SampleDB``. You will need two containers, one for a PostgreSQL database and another for SampleDB itself, and a directory to store all files in.

If you would like to set up a development version of SampleDB instead, please see the `contribution guide <https://github.com/sciapp/sampledb/blob/develop/CONTRIBUTING.md>`_.

If you do not have Docker installed yet, please `install Docker <https://docs.docker.com/engine/install/>`_.

First, start your database container:

.. code-block:: bash

    docker run \
        -d \
        -e POSTGRES_PASSWORD=password \
        -e PGDATA=/var/lib/postgresql/data/pgdata \
        -v `pwd`/pgdata:/var/lib/postgresql/data/pgdata:rw \
        --restart=always \
        --name sampledb-postgres \
        postgres:12


Next, start the SampleDB container:

.. code-block:: bash

    docker run \
        -d \
        --link sampledb-postgres \
        -e SAMPLEDB_CONTACT_EMAIL=sampledb@example.com \
        -e SAMPLEDB_MAIL_SERVER=mail.example.com \
        -e SAMPLEDB_MAIL_SENDER=sampledb@example.com \
        -e SAMPLEDB_ADMIN_PASSWORD=password \
        -e SAMPLEDB_SQLALCHEMY_DATABASE_URI=postgresql+psycopg2://postgres:password@sampledb-postgres:5432/postgres \
        -e SAMPLEDB_FILE_STORAGE_PATH=/home/sampledb/files/ \
        -v `pwd`/files:/home/sampledb/files:rw \
        --restart=always \
        --name sampledb \
        -p 8000:8000 \
        sciapp/sampledb:0.19.3

This will start a minimal SampleDB installation at ``http://localhost:8000`` and allow you to sign in with the username ``admin`` and the password ``password``.

Step 2: Basic Configuration
---------------------------

After visiting ``http://localhost:8000`` to verify that SampleDB runs, you should sign in as administrator and update your password.

Although you have a minimal SampleDB installation now, you will need to perform some further configuration. The code block starting the SampleDB Docker container contains several environment variables with placeholder values. After stopping the container, please start it with valid values for:

- SAMPLEDB_CONTACT_EMAIL
- SAMPLEDB_MAIL_SERVER
- SAMPLEDB_MAIL_SENDER

Your mail server might require additional settings and you can find descriptions of these variables in the :ref:`list of configuration environment variables<configuration>`.

You can use this moment to set other configuration variables, e.g. for the service name, service description, or imprint and privacy policy links.

You can restart the SampleDB container to change configuration variables at any time, so you can always return to this step later on.

Step 3: Creating Actions
------------------------

Your first step in actually using SampleDB should be to set up your first :ref:`actions`. Actions represent processes that create samples, perform measurements or running simulations.

1. We recommend that you start by setting up an Action for creating a generic sample in SampleDB with:
    - a name,
    - a creation date,
    - a multiline text description, and
    - an optional input sample.

   This will allow you to import samples as you move on, until all your processes for sample creation, preparation or import are modelled as individual Actions.
2. Next, we recommend that you set up two equally generic Actions for performing a measurement and for running a simulation. This is to allow your users to enter any measurements or simulations into your SampleDB installation which have not yet been modelled as individual Actions.
3. With these three generic Actions in place, you should then pick one process and model it as an Action. Try to capture all the information that would be required to reproduce a sample creation, measurement or simulation, and add it as properties to this Action's schema.

You can then improve your Actions' schemas and add new Actions as you become more experienced using SampleDB and gather feedback from your users.

Instruments
```````````

As you add more Actions, you may want to group some Actions by the instrument they are performed with and give the instrument scientists control over these Actions. To do so:

- create a new :ref:`Instrument <instruments>`,
- assign :ref:`instrument_scientists`, and
- create :ref:`actions` for this instrument.

Step 4: Preparing SampleDB for Production
-----------------------------------------

After the previous steps, you can fully evaluate SampleDB locally using the admin user. At this stage, however, you might want to make your SampleDB installation available to others and run SampleDB in production. We **strongly** recommend that you set up :ref:`TLS Termination<tls_termination>` and that you regularly create :ref:`backups <backup_and_restore>`.

Step 5: User Management
-----------------------

At this time, SampleDB users can either sign in using a username and password specific to SampleDB, or by using LDAP if it has been enabled using the :ref:`LDAP configuration variables<ldap_configuration>`.

If your facility already has an LDAP system for user management, we recommend that you configure LDAP in SampleDB so that users can use their existing credentials.

Otherwise, you can invite your users using the :ref:`User Invitation Form<authentication>`.

Next Steps
----------

- You might want to create :ref:`groups` or :ref:`projects` to model your existing team structures. While this can be useful, it is completely optional as users can set these up themselves.
- You might want to create a basic hierarchy of :ref:`locations`. Like groups, users can create these themselves so this is optional.
- If you already have a JupyterHub installation or want to set up one, you might want to enable SampleDB :ref:`JupyterHub support <jupyterhub_support>`.
- SampleDB is still under active development. When a new version is released, you should consider :ref`upgrading your SampleDB installation <upgrading>`.
- If you have any questions about SampleDB or run into any issues setting up or running SampleDB, please `create an issue on GitHub <https://github.com/sciapp/sampledb/issues/new>`_.
.. _upgrading:

Upgrading
=========

To update an existing SampleDB installation, please create a :ref:`backup <backup_and_restore>` and then re-start SampleDB using the Docker image of the new version.

There is no supported path for downgrading SampleDB at this time.
.. _languages:

Languages
=========

Although SampleDB was primarily developed in English, it now supports two ways of using it with other languages:

- SampleDB itself has been translated to German and can be translated in other languages as well in the future.
- Object metadata fields, actions, instruments, groups and other user-generated content can contain translations for various languages, defined by the administrator.

If you wish to contribute to a translation of SampleDB, see the `SampleDB contribution guide <https://github.com/sciapp/sampledb/blob/develop/CONTRIBUTING.md>`_.

By default, a SampleDB instance will have English and German as built-in languages for user input, with only English enabled. This way all text fields will appear as they were before this change. If you enable German input or create a new language that is enabled for input, many fields for which a translation would be reasonable will offer users to add a translation in the newly enabled language. When a user uses a language for which some content does not provide a translation, the content will be shown in English instead.

.. figure:: ../static/img/generated/translations.png
    :alt: Fields for providing translations for an Action name

    Fields for providing translations for an Action name

To add a new language for user input, administrators can go to *More* → *Languages* and click *Add new language*. The field *Datetime format for datetime* expects a format in the syntax used by the Python `datetime module <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`_ and the field *Datetime format for moment* expects a format in the syntax used by the JavaScript `moment.js package <https://momentjs.com/docs/#/parsing/string-format/>`_. These formats will be used for user input, and formats based on the language code will be use used for displaying datetimes.
.. _administration_scripts:

Administration Scripts
======================

Most administration features are available to administrators through the user interface. Some features which are rarely used or which would be useful for scripting purposes are available as scripts, so that they can be run using the SampleDB command line interface.

If you are running SampleDB using Docker, you can execute these scripts using ``docker exec``. To get a list of all available scripts, you can run:

.. code-block:: bash

    docker exec sampledb env/bin/python -m sampledb help

To then get a list of all actions, for example, you can run:

.. code-block:: bash

    docker exec sampledb env/bin/python -m sampledb list_actions

To get information on how to run one of these scripts, you can run pass the `help` parameter to it:

.. code-block:: bash

    docker exec sampledb env/bin/python -m sampledb list_actions help
.. _dataverse_export:

Dataverse Export
================

Dataverse is web application for managing research data, including its publication. In contrast to |service_name|, which focuses on managing flexible, process-specific metadata within an organization, Dataverse has an administrator-managed schema system and aims at the sharing and visibility of research data.

If the :ref:`Dataverse Export configuration variables <dataverse_configuration>` have been set, users with **GRANT** permissions can export object information and files to a Dataverse server. A draft dataset is then created there, which the user can edit, review and publish according to the Dataverse permissions.

To be able to represent the process-specific metadata in a format compatible with Dataverse, |service_name| relies both on the Dataverse *Citation Metadata* block for generic metadata and the *Process Metadata* block developed by the University of Stuttgart as part of `EngMeta <https://www.izus.uni-stuttgart.de/fokus/engmeta/>`_.

By default, none of the files or process-specific metadata will be exported. Users can enable individual files and properties, and the schema of a property may include a ``dataverse_export`` boolean to set whether a property should be exported by default. This way, e.g. instrument scientists can provide a suggestion for what should be exported for a particular action.
