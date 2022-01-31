============
Introduction
============

Texcavator is a text mining web application.
Its back-end is written in the Django_ framework, while the front-end relies heavily upon the Dojo_ toolkit.
The web application allows users to fire queries to an Elasticsearch index and will convert the responses into visualizations (for which we use D3.js_).
For larger requests (i.e. word cloud generation over multiple documents) we use Celery_ as a background task queue.

======
Python
======

Texcavator uses the following Python packages (see ``requirements.in``):

- Django_: web framework.
- elasticsearch_: elasticsearch connector for Python.
- Celery_: distributed task queue.
- redis_: key-value store, connects Django with Celery.
- MySQL-python_: MySQL connector for Python. Can be replaced with a connector to any other RDMS (e.g. pyscopg2)
- sqlparse_: allows raw SQL operations on databases.
- requests_: makes API requests more easy.
- DAWG_: allows creation of direct acyclic word graphs, which allows us to efficiently save frequencies of all words in the corpus.
- dicttoxml_: converts dictionaries to XML; used for generating XML output in the download functionality.
- Sphinx_: documentation generator.
- nose_: test enhancement suite.
- django-nose_: nose connector for Django.

.. _Django: https://www.djangoproject.com/
.. _Celery: http://www.celeryproject.org/
.. _redis: http://redis.io/
.. _elasticsearch: http://elasticsearch-py.readthedocs.io/en/master/
.. _MySQL-python: https://pypi.python.org/pypi/MySQL-python
.. _sqlparse: https://sqlparse.readthedocs.io/en/latest/
.. _requests: http://docs.python-requests.org/en/latest/
.. _DAWG: https://github.com/kmike/DAWG
.. _dicttoxml: https://pypi.python.org/pypi/dicttoxml
.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _nose: http://nose.readthedocs.io/en/latest/
.. _django-nose: https://django-nose.readthedocs.io/en/latest/

==========
JavaScript
==========

Texcavator heavily uses the Dojo_ toolkit for its frontend.
jQuery_ is used in places where Dojo doesn't have ready solutions available.

For our visualizations we leverage D3.js_ (word clouds and time lines), nvd3_ (metadata visualizations) and Cal-HeatMap_ (heatmaps).

To aid users in using correct syntax when querying the elasticsearch index, we have created a parsing expression grammar in PEG.js_.

.. _Dojo: https://dojotoolkit.org/
.. _jQuery: http://jquery.com/
.. _D3.js: https://d3js.org/
.. _nvd3: http://nvd3.org/
.. _Cal-HeatMap: http://cal-heatmap.com/
.. _PEG.js: http://pegjs.org/
==========
Deployment
==========

Texcavator is currently battle-tested on a Debian-based distribution with Apache_ (`mod_wsgi`), MySQL_, Redis_ and Postfix_.
Starting from a clean install, this guide should get you going.
Feel free to replace one or all of these components, but then the usual caveats apply.

.. _Apache: https://httpd.apache.org/
.. _MySQL: https://www.mysql.com/
.. _Redis: http://redis.io/
.. _Postfix: http://www.postfix.org/

************
Dependencies
************

First, let's install some dependencies::

    sudo apt-get update
    sudo apt-get upgrade
    sudo apt-get install apache2
    sudo apt-get install mysql-server
    sudo apt-get install redis-server
    sudo apt-get install postfix
    sudo apt-get install libapache2-mod-wsgi libmysqlclient-dev libxml2-dev libxslt-dev

During installation of ``mysql-server`` and ``postfix``, make sure to jot down the details to set these in Texcavator's settings later on.

*************
elasticsearch
*************

You can run elasticsearch_ on the same server as Texcavator, but we would propose setting up a separate server.
See https://www.elastic.co/guide/en/elasticsearch/guide/current/deploy.html for some instructions.

.. _elasticsearch: https://www.elastic.co/

**********************
Retrieving the sources
**********************

This should be easy, as you can retrieve the sources from GitHub_::

    sudo apt-get install git
    git clone https://github.com/UUDigitalHumanitieslab/texcavator.git

.. _GitHub: https://github.com/

*******************
Virtual environment
*******************

It's customary to install the Python dependencies in a virtual environment, so::

    sudo apt-get install python-pip
    pip install virtualenv

Then, activate the virtual environment and install the requirements using the ``requirements.txt`` provided in the source directory::

    pip install -r requirements.txt

********************
Application settings
********************

Copy ``texcavator/settings_local_default.py`` to ``texcavator/settings_local.py``.
The latter file is not kept under version control.
Then, change the settings in ``texcavator/settings_local.py`` and (if necessary) ``texcavator/settings.py`` to reflect your setup.

Run the commands in the section **Prerequisite commands** from ``README.rst``.

********
Database
********

Create a database (and user) in your RDMS. After that, run::

    python manage.py migrate

from your virtual environment to migrate the database to the latest and greatest version.

******
Celery
******

For Celery_, it's (probably) best to use the ``celery multi`` command, from your virtual environment, with a specific celery user.
The following command would start 3 workers with 8 concurrent processes, more than enough for Texcavator::

   ${virtualenv}/bin/celery multi restart worker1 worker2 worker3 --app=texcavator.celery:app --workdir=${app_source} --time-limit=300 --concurrency=8 --logfile=${logging_location}/celery/%N.log --pidfile=${celery_root}/%N.pid --uid=${celery_user} --gid=${celery_group}

.. _Celery: http://www.celeryproject.org/

******
Apache
******

Setting up Apache and mod_wsgi to serve Django is quite a challenge.
Luckily, there is some excellent documentation available online.
Be sure to update settings.py and settings_local.py according to your settings.
You can collect static files using the ``collectstatic`` management command::

    python manage.py collectstatic

======
Puppet
======

We have a Puppet module available that will install the application:
you'll only have to install the dependencies, and afterwards the complete application environment will be created.
Drop us an e-mail if you're interested in that.
Mapping
=======

The elasticsearch mapping we currently use is as follows::

    {
      "settings": {
        "analysis" : {
          "analyzer" : {
            "dutch_analyzer" : {
              "type" : "custom",
              "tokenizer": "standard",
              "filter" : ["standard", "lowercase", "dutch_stemmer"]
            }
          },
          "filter" : {
            "dutch_stemmer" : {
              "type" : "stemmer",
              "name" : "dutch_kp"
            }
          }
        }
      },
      "mappings": {
        "doc": {
          "properties" : {
            "article_dc_subject": {
              "type": "string",
              "include_in_all": "false",
              "index": "not_analyzed"
            },
            "article_dc_title": {
              "type": "string",
              "term_vector": "with_positions_offsets_payloads",
              "fields": {
                "stemmed": {
                  "type": "string",
                  "analyzer": "dutch_analyzer",
                  "term_vector": "with_positions_offsets_payloads"
                }
              }
            },
            "identifier": {
              "type": "string",
              "include_in_all": "false",
              "index": "not_analyzed"
            },
            "paper_dc_date": {
              "format": "dateOptionalTime",
              "type": "date"
            },
            "paper_dc_title": {
              "type": "string",
              "term_vector": "with_positions_offsets_payloads",
              "fields": {
                "raw": {
                  "type": "string",
                  "index": "not_analyzed"
                }
              }
            },
            "paper_dcterms_spatial": {
              "type": "string",
              "include_in_all": "false",
              "index": "not_analyzed"
            },
            "paper_dcterms_temporal": {
              "type": "string",
              "include_in_all": "false",
              "index": "not_analyzed"
            },
            "text_content": {
              "type": "string",
              "term_vector": "with_positions_offsets_payloads",
              "fields": {
                "stemmed": {
                  "type": "string",
                  "analyzer": "dutch_analyzer",
                  "term_vector": "with_positions_offsets_payloads"
                }
              }
            }
          }
        }
      }
    }


An example document would then be::

    {
        "article_dc_subject": "newspaper",
        "article_dc_title": "Test for Texcavator",
        "identifier": "test1",
        "paper_dc_date": "1912-04-15",
        "paper_dc_title": "The Texcavator Test",
        "paper_dcterms_spatial": "unknown",
        "paper_dcterms_temporal": "daily",
        "text_content": "This is a test to see whether Texcavator works!"
    }
::

     _____                             _             
    |_   _|____  _____ __ ___   ____ _| |_ ___  _ __ 
      | |/ _ \ \/ / __/ _` \ \ / / _` | __/ _ \| '__|
      | |  __/>  < (_| (_| |\ V / (_| | || (_) | |   
      |_|\___/_/\_\___\__,_| \_/ \__,_|\__\___/|_|   


Copyright '`Digital Humanities lab` @ Utrecht University',`Netherlands eScience Center`_, `University of Amsterdam`_. 

From 2015 onwards developed by the `Digital Humanities lab`_, Utrecht University.

.. _`Netherlands eScience Center`: https://www.esciencecenter.nl/
.. _`University of Amsterdam`: http://www.uva.nl/en/
.. _`Digital Humanities lab`: http://dig.hum.uu.nl/

Developer quick-start
=====================

************
Dependencies
************

Before installing Texcavator, make sure your packages are up-to-date and
a relational database (we use MySQL_) and Redis_ server are present on the system.
In apt-based Linux distros like Ubuntu/Debian, issue::

    sudo apt-get update
    sudo apt-get upgrade
    sudo apt-get install mysql-server redis-server

Make sure both servers are running. Furthermore, you will need a few development packages::

    sudo apt-get install libmysqlclient-dev libxml2-dev libxslt-dev python-dev

.. _MySQL: https://www.mysql.com/
.. _Redis: http://redis.io/

************
Installation
************

To install Texcavator, clone the repository (using git_) in your home directory
and make a virtualenv_, activate it, and install the requirements::

    sudo apt-get install git python-pip
    pip install virtualenv
    cd ~
    git clone https://github.com/UUDigitalHumanitieslab/texcavator.git
    mkdir .virtualenvs
    virtualenv .virtualenvs/texc
    source .virtualenvs/texc/bin/activate
    pip install -r texcavator/requirements.txt

In ``texcavator/settings.py``, you can change the path to the log file, if you like.

Copy ``texcavator/settings_local_default.py`` to ``texcavator/settings_local.py``.
The latter file is not kept under version control.

In ``texcavator/settings_local.py``, set up the database; for a quick test, you can use SQLite::

    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': os.path.join(PROJECT_PARENT, 'db.sqlite3')
        }
    }

.. _git: https://git-scm.com/
.. _virtualenv: http://virtualenv.readthedocs.io/

*************
Elasticsearch
*************

The next step is to load your data into an elasticsearch_ index.
To get started using elasticsearch see the quickstart_.

In ``texcavator/settings_local.py``, you can specify the elasticsearch host and port
(typically elasticsearch runs on ``localhost:9200``).
Texcavator assumes by default the data is in an index called ``kb`` and
that the documents are stored in a type ``doc`` that has at least the following fields:

* article_dc_subject
* article_dc_title
* identifier
* paper_dc_date
* paper_dc_title
* paper_dcterms_spatial
* paper_dcterms_temporal
* text_content

We use the mapping specified in ``mapping.rst``.

.. _elasticsearch: https://www.elastic.co/
.. _quickstart: https://www.elastic.co/guide/en/elasticsearch/reference/current/getting-started.html

*********************
Prerequisite commands
*********************

Texcavator requires you to install some external packages and management commands in order to function correctly.
Before issuing the commands below, make sure Elasticsearch, MySQL and Redis are still running at the specified ports.

* Install ShiCo, which allows for visualizing shifting concepts over time::

    python install-shico.py

* Populate the database::

    python manage.py migrate

* Create a Django superuser. The username and password you pick will be the administrator account for Texcavator::

    python manage.py createsuperuser

* Run the management command ``gatherstatistics`` to be able to display timelines::

    python manage.py gatherstatistics

* Run the management command ``add_stopwords`` to add a default list of (Dutch) stop words::

    python manage.py add_stopwords stopwords/nl.txt

* Run the management command ``gathertermcounts`` to be able to create word clouds normalized for inverse document frequency::

    python manage.py gathertermcounts

* Run the management command ``add_guest_user`` to add a guest environment (with limited options)::

    python manage.py add_guest_user

.. _Dojo: http://dojotoolkit.org/

******************
Development server
******************

First, make sure Elasticsearch, MySQL and Redis are still running at the specified ports.
Then, start Celery and Django's built-in webserver::

    celery --app=texcavator.celery:app worker --loglevel=info
    # In a separate terminal
    python manage.py runserver

Texcavator is now ready for use at ``http://127.0.0.1:8000``.

Downloading of query data requires a running SMTP server; you can use Python's built-in server for that::

    # In a separate terminal
    python -m smtpd -n -c DebuggingServer localhost:1025

Deployment
==========

You can find instructions for deploying Texcavator in ``deployment.rst``

Documentation
=============

The documentation for Texcavator is in Sphinx_. You can generate the documentation by running::

    make html

in the /doc/ directory.

.. _Sphinx: http://sphinx-doc.org/index.html


License
=======

Texcavator is distributed under the terms of the Apache2 license. See ``LICENSE`` for details.
.. Texcavator documentation master file, created by
   sphinx-quickstart on Mon Feb 23 18:41:39 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Texcavator's documentation!
======================================

Contents:

.. toctree::
   :maxdepth: 2

   modules/models
   modules/views
   modules/other
   modules/commands


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Additional Modules
==================

Texcavator
----------

Utils
+++++

.. automodule:: texcavator.utils 
    :members:

Query
-----

Burstsdetector
++++++++++++++

.. automodule:: query.burstsdetector.bursts 
    :members:

Download
++++++++

.. automodule:: query.download 
    :members:

Utils
+++++

.. automodule:: query.utils 
    :members:

Services
--------

Elasticsearch_biland
++++++++++++++++++++

.. warning::
   
   Must be replaced as soon as possible (however, this requires a new user
   interface).

.. automodule:: services.elasticsearch_biland
    :members:

Es
++

.. note::

   This is the ElasticSearch functionality that should be used.

.. automodule:: services.es
    :members:

Tasks
+++++

.. automodule:: services.tasks
    :members:
Management Commands
===================

Query
-----

export_queries
++++++++++++++

.. automodule:: query.management.commands.export_queries
    :members:

export_stopwords
+++++++++++++++++

.. automodule:: query.management.commands.export_stopwords
    :members:

gatherstatistics
+++++++++++++++++

.. automodule:: query.management.commands.gatherstatistics
    :members:

add_stopwords
+++++++++++++

.. automodule:: query.management.commands.add_stopwords
    :members:

Services
--------

gatherdocids
+++++++++++++

.. automodule:: services.management.commands.gatherdocids
    :members:

gatherqueryterms
++++++++++++++++

.. automodule:: services.management.commands.gatherqueryterms
    :members:

esperformance
+++++++++++++

.. automodule:: services.management.commands.esperformance
    :members:

tag_wordclouds
++++++++++++++

.. automodule:: services.management.commands.tag_wordclouds
    :members:

termvector_wordclouds
+++++++++++++++++++++

.. automodule:: services.management.commands.termvector_wordclouds
    :members:

tv_wc_no_client
+++++++++++++++++++++

.. automodule:: services.management.commands.tv_wc_no_client
    :members:

weighted_queries
+++++++++++++++++++++

.. automodule:: services.management.commands.weighted_queries
    :members:
Models
======

Query 
-----
.. automodule:: query.models
    :members:

Services
--------
.. automodule:: services.models
    :members:
Views
=====

Texcavator
----------
.. automodule:: texcavator.views 
    :members:

Query
-----
.. automodule:: query.views 
    :members:

Services
--------
.. automodule:: services.views 
    :members:


