# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

* BASE_PATH env var to overwrite base path (#15)

## [0.2.0] - 2018-12-06

### Fixed

* workspace is already in use in another JupyterLab window (#11)

### Changed

* Upgraded to Connexion v2
* Upgraded to OpenAPI v3

## [0.1.0] - 2018-10-09

Initial release
# experiment-launcher

[![Python application](https://github.com/eWaterCycle/experiment-launcher/actions/workflows/python-app.yml/badge.svg)](https://github.com/eWaterCycle/experiment-launcher/actions/workflows/python-app.yml)
[![SonarCloud quality gate](https://sonarcloud.io/api/project_badges/measure?project=experiment-launcher&metric=alert_status)](https://sonarcloud.io/dashboard?id=experiment-launcher)
[![SonarCloud coverage](https://sonarcloud.io/api/project_badges/measure?project=experiment-launcher&metric=coverage)](https://sonarcloud.io/component_measures?id=experiment-launcher&metric=coverage)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1453264.svg)](https://doi.org/10.5281/zenodo.1453264)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

eWaterCycle Experiment Launcher a webservice to generate and launch a Jupyter notebook.

The API of the webservice is described in OpenAPI specification at [openapi.yaml](https://github.com/eWaterCycle/experiment-launcher/blob/main/ewatercycle_experiment_launcher/swagger.yaml) and can be seen in [Swagger UI](http://petstore.swagger.io/?url=https://raw.githubusercontent.com/eWaterCycle/experiment-launcher/main/ewatercycle_experiment_launcher/openapi.yaml)

# Install

Instructions below have been tested on Linux, but should also work on OSX and Windows Subsystem for Linux.

## JupyterHub server

The experiment launcher requires a JupyterHub server.

JupyterHub can be installed with the following commands
```bash
pip3 install jupyterhub jupyterlab
sudo npm install -g configurable-http-proxy
```

JupyterHub must accept calls from the experiment launcher service to start a notebook server for any hub user and upload a notebook.
In the JupyterHub configuration file (jupyterhub_config.py) you must register the experiment launcher as a service with admin permissions and a token which the launcher must use to communicate with JupyterHub.

The token (shared secret) for the launcher can be generated with

```bash
openssl rand -hex 32
```

A JupyterHub config file can be made using the [./jupyterhub_config.py.example](jupyterhub_config.py.example) file with

```bash
cp jupyterhub_config.py.example jupyterhub_config.py
```

Put the generated token in the config file by editing it with your favourite editor

```bash
nano jupyterhub_config.py
```

JupyterHub can be started with

```bash
jupyterhub
```

Test JupyterHub by going to http://localhost:8000 and login with your OS credentials.

## Installation for production

```bash
pip install ewatercycle_experiment_launcher
```

## Installation for development

To install the launcher in development mode clone the repo and run

```bash
pip install -r requirements_dev.txt
```

# Run

The launcher must be given the same token as configured in the JupyterHub config file as `JUPYTERHUB_TOKEN` environment variable. To get it you can use following oneliner:

```bash
# Use token from jupyterhub_config.py
export JUPYTERHUB_TOKEN=$(python3 -c "from traitlets.config import Application;\
    print([s['api_token'] for s in \
    next(Application._load_config_files('jupyterhub_config'))[0]['JupyterHub']['services'] \
    if s['name'] == 'experiment-launcher'][0])")
```

To start launcher use

```bash
# JUPYTERHUB_URL is URL where JupyterHub is running. If path like `/jupyter` then origin header is appended.
export JUPYTERHUB_URL=http://172.17.0.1:8000
gunicorn -w 4 -b 0.0.0.0:8888 ewatercycle_experiment_launcher.serve:app
```

Goto http://localhost:8888/ui/ for Swagger UI.

The JupyterHub and Experiment Launcher both use local OS accounts for authentication and authorization.

In the Swagger UI you must authorize before trying an operation.

When running on Internet make sure https is enforced so the authentication is secure.

The webservice by default runs on `/` base path. This can be changed by setting the `BASE_PATH` environment variable.
For example `export BASE_PATH=/launcher` will host the Swagger UI on http://localhost:8888/launcher/ui/ .
* eWaterCycle Experiment Launcher version:
* Python version:
* Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/ewatercycle/ewatercycle_experiment_launcher/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

eWaterCycle Experiment Launcher could always use more documentation, whether as part of the
official eWaterCycle Experiment Launcher docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/ewatercycle/ewatercycle_experiment_launcher/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)


Adding a new notebook type
~~~~~~~~~~~~~~~~~~~~~~~~~~

The web service has an path for each type of notebook.

To add a new type of notebook the following steps must be performed:

1. In `ewatercycle_experiment_launcher/openapi.yaml` file create a new path
    * The http method should be `post`
    * The requestBody should be a json object which includes a `notebook` property of schema type `NotebookRequest`
    * The 200 response should of response type NotebookResponse
    * The default response should of response type ErrorResponse
2. In `ewatercycle_experiment_launcher/api` directory create a file with same name as the chosen path +'.py'
    * Create a `post()` function, using the following template

.. code-block:: python

        from ewatercycle_experiment_launcher.process import process_notebook

        def post(body):
            """Generate notebook and launch it

            Args:
                body: The json POST body as a Python dictionary
            """
            nb = ... # <Add code that generates a nbformat.NotebookNode object>
            return process_notebook(body['notebook'], nb)

3. Write unit tests in `tests/api/` directory

Make a Pull Request after the new notebook type has been implemented and tested.

Get Started!
------------

Ready to contribute? Here's how to set up `ewatercycle_experiment_launcher` for local development.

1. Fork the `ewatercycle_experiment_launcher` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/ewatercycle_experiment_launcher.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv ewatercycle_experiment_launcher
    $ cd ewatercycle_experiment_launcher/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 ewatercycle_experiment_launcher tests
    $ python setup.py test or py.test
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.7, 3.4, 3.5 and 3.6, and for PyPy. Check
   https://travis-ci.org/ewatercycle/ewatercycle_experiment_launcher/pull_requests
   and make sure that the tests pass for all supported Python versions.

Release
-------

A reminder for the maintainers on how to release a new version.

1. Make sure tests pass by running::

    $ pytest

2. Bump the version by running::

    $ bumpversion patch # possible: major / minor / patch

3. Update or create an entry for the new version in the `CHANGELOG.md` file
4. Make sure all your changes are committed and pushed
5. Publish to pypi with::

    $ rm -rf dist
    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*

6. Create GitHub release
7. Update DOI in `CITATION.cff` file
.. include:: ../CONTRIBUTING.rst
.. include:: ../HISTORY.rst
=====
Usage
=====

To use eWaterCycle Experiment Launcher in a project::

    import ewatercycle_experiment_launcher
.. include:: ../README.rst
.. highlight:: shell

============
Installation
============


Stable release
--------------

To install eWaterCycle Experiment Launcher, run this command in your terminal:

.. code-block:: console

    $ pip install ewatercycle_experiment_launcher

This is the preferred method to install eWaterCycle Experiment Launcher, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for eWaterCycle Experiment Launcher can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/ewatercycle/ewatercycle_experiment_launcher

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/ewatercycle/ewatercycle_experiment_launcher/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/ewatercycle/ewatercycle_experiment_launcher
.. _tarball: https://github.com/ewatercycle/ewatercycle_experiment_launcher/tarball/master
Welcome to eWaterCycle Experiment Launcher's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   readme
   installation
   usage
   modules
   contributing
   history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
