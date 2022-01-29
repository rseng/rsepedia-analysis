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
