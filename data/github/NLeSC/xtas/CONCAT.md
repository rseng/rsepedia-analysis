.. image:: https://api.travis-ci.org/NLeSC/xtas.png?branch=master
   :target: https://travis-ci.org/NLeSC/xtas

xtas
====

Distributed text analysis suite based on Celery.

Copyright University of Amsterdam, Netherlands eScience Center and
contributors, distributed under the Apache License (see ``AUTHORS.txt``,
``LICENSE.txt``). Parts of xtas use GPL-licensed software, such as the
Stanford NLP tools, and datasets that may incur additional restrictions.
Check the documentation for individual functions.


Quickstart
----------

Install::

    pip install xtas

Start worker::

    python -m xtas.worker --loglevel=info

Start web frontend::

    python -m xtas.webserver

For full documentation, please visit http://nlesc.github.io/xtas/.
REST API
========

When run as a webserver, xtas exposes a simple REST API to its tasks.
To tokenize a string using this API, you need to make two HTTP requests,
one to start a task and one to wait until it has completed.
Assuming default settings for the webserver, you can try::

    $ curl -H "Content-type: text/plain" -X POST -d 'Hello, world!' \
        http://127.0.0.1:5000/run/tokenize
    11d0f158-abdb-4a0a-860d-b365456122f4

The last line shows an id for the job that you submitted.
Since tokenization is very fast, you can immediately query the REST endpoint
again to get the results from the tokenization::

    $ curl http://127.0.0.1:5000/result/11d0f158-abdb-4a0a-860d-b365456122f4
    ["Hello", ",", "world", "!"]

Passing arguments to the task is also possible. To do so, you have to send the
API request as JSON, as follows::

    $ curl -H "Content-type: application/json" -X POST -d \
        '{"data": "Hello, world!", "arguments": {"output": "rank"}}'

All single-document tasks can be called in this way.
The arguments should be JSON-wrapped versions
of the arguments ordinarily passed to the Python functions
(see :ref:`api`).
.. _setup:

Getting started
===============

xtas can be run locally or on a cluster; the :ref:`overview` explains this
in more detail. This guide shows first how to get xtas running locally,
then how to run it as a distributed service.


Installation
------------

xtas runs on Linux.
It depends on Python 2.7, SciPy, the Java Runtime (JRE) and RabbitMQ.
On a Debian/Ubuntu/Linux Mint system, these dependencies can be installed with::

    sudo apt-get install libatlas-dev liblapack-dev rabbitmq-server \
        python-scipy openjdk-7-jre python-virtualenv build-essential \
        python-pip libxslt-dev

For CentOS/Fedora/Red Hat-style systems (this list is incomplete)::

    sudo yum install atlas-devel java-1.7.0-openjdk lapack-devel \
        libxslt-devel numpy python-devel rabbitmq-server scipy

(RabbitMQ is an `EPEL <https://fedoraproject.org/wiki/EPEL>`_ package.)

Next, set up a virtualenv for xtas::

    virtualenv --system-site-packages /some/where
    . /some/where/bin/activate

Use `pip <https://pypi.python.org/pypi/pip/1.1>`_ to install xtas.
To get the latest release::

    pip install xtas

To get the bleeding edge version from GitHub::

    pip install git+https://github.com/NLeSC/xtas.git

For installing from local source (if you want to modify xtas),
see :doc:`extending`.

Try it out in a Python shell::

    >>> from xtas.tasks import tokenize
    >>> tokenize("Hello, world!")
    [u'Hello', u',', u'world', u'!']

If you only want to use xtas as a library on your local machine, you're done.
Check the :ref:`api` for what xtas can do.

Note that xtas will download, install and run several free, open source
programs on your computer as needed. See LICENSE.rst for (legal) details.


Distributed xtas
----------------

To run xtas in distributed mode, you need to have RabbitMQ
and, optionally, Elasticsearch running on some machine.
xtas by default assumes that these are running locally on their standard ports.
If they are not, run::

    python -m xtas.make_config

and edit the generated configuration file, ``xtas_config.py``,
to point xtas to RabbitMQ ``BROKER_URL`` (see `Celery configuration
<http://docs.celeryproject.org/en/latest/configuration.html>`_ for details)
and Elasticsearch.
Make sure this file is somewhere in your ``PYTHONPATH``
(test this with ``python -c 'import xtas_config'``).

Then start an xtas worker::

    python -m xtas.worker --loglevel=info &

If you want to use the xtas REST API, also start the webserver::

    python -m xtas.webserver &

Verify that it works::

    curl http://localhost:5000/tasks | python -m json.tool

You should see a list of supported tasks.

Now to perform some actual work, make sure Elasticsearch is populated with
documents, and visit a URL such as

    http://localhost:5000/run_es/morphy/20news/post/1/text

This runs the Morphy morphological analyzer on the "text" field of "post" 1
in ES index "20news". After some time, the results from Morphy are written to
a child document of this post, that can be obtained using::

    curl http://localhost:9200/20news/post__morphy/1?parent=1

You can now run the unittest suite using::

    python setup.py test

in the source directory (``pip install nose`` if needed). This requires a
running worker process and Elasticsearch. Running the tests first is a good
idea, because it will fetch some dependencies (e.g. NLTK models) that will
otherwise be fetched on demand.

To learn more about using xtas as a distributed text analysis engine,
see the :ref:`tutorial`.


Running as a service
--------------------

xtas can be run as a service on Linux. See the directory ``init.d`` in the
xtas source distribution for example init scripts.
.. _extending:

Extending xtas
==============

This is a short guide to extending xtas to suit specific needs
and (optionally) contributing code back.
It describes how to write new tasks and tie them in with the package.


Writing new tasks
-----------------

If you have a custom task that you want xtas to perform,
then you can add it as follows.
Suppose you want to perform sentiment analysis on French text
using the `Pattern <http://www.clips.ua.ac.be/pages/pattern>`_ toolkit.
First, install Pattern::

    pip install pattern

Now define an xtas task that uses Pattern to process a French text.
Put the following in a file, say, ``pattern_tasks.py``::

    import pattern.fr
    from xtas.core import app

    @app.task
    def fr_sentiment(text):
        """Perform sentiment analysis on French text.

        Returns the text (for reference), a subjectivity score,
        and a positive/negative score.
        """
        return (text,) + pattern.fr.sentiment(text)

Make sure this file can be imported::

    >>> from pattern_tasks import fr_sentiment

If the above gave an error, adjust your ``PYTHONPATH`` (in the shell)::

    export PYTHONPATH=.:$PYTHONPATH

Now adjust the xtas configuration. Run ``python -m xtas.make_config`` to get
a configuration file ``xtas_config.py`` in the current directory. At the bottom
of the file is an empty list called ``EXTRA_MODULES``. Put your module in it::

    EXTRA_MODULES = [
        'pattern_tasks',
    ]

Now restart the worker. It should report ``pattern_tasks.fr_sentiment``
as its first task, followed by all the built-in tasks.
You can now run your task function asynchronously, e.g. in a Python shell::

    >>> from pattern_tasks import fr_sentiment
    >>> result = fr_sentiment.apply_async(["Bon!"])
    >>> result.get()
    ['Bon!', 0.875, 0.7]

If you also restart the webserver, you should see the new task in the list of
single-document tasks::

    $ curl -s http://localhost:5000/tasks | python -m json.tool | grep pattern
        "pattern_tasks.fr_sentiment",

To use the custom task from the REST API, e.g. with ``/run_es``, give its
fully qualified name (``pattern_tasks.fr_sentiment``).
Only built-in tasks have their name abbreviated to not include the module name.

.. note::
   The webserver will currently assume your task is a single-document one,
   rather than a batch task. This is a `known defect
   <https://github.com/NLeSC/xtas/issues/40>`_.


Contributing code
-----------------

If you have code that is reusable for others, and you want to and are legally
able to distribute it, then we're happy to consider it for inclusion in
xtas. Make sure you have copyright to the code you write or your employer
gives you permission to contribute under the terms of the Apache License
(``LICENSE.txt`` in the main source directory).

Fork the main repository on GitHub, then install this instead of the released
version. First make a new virtualenv::

    virtualenv --system-site-packages /some/where/xtas-work
    . /some/where/xtas-work/bin/activate

Then, in the xtas source directory, issue::

    pip install .

When you make changes, issue this command to update the virtualenv::

    pip install --upgrade --no-deps .

To contribute code back, commit your changes to a separate branch.
Push this branch to GitHub and do a pull request. Your code will be reviewed
before pulling.

In the case of new tasks, put them in either ``xtas/tasks/single.py`` or
``xtas/tasks/cluster.py``, depending on the type of task. Follow the
conventions laid out in the docstrings of the modules. All xtas code should
conform to `PEP 8 <http://legacy.python.org/dev/peps/pep-0008/>`_, the style
guide for the Python standard library. Use the `pep8
<http://pep8.readthedocs.org/en/latest/>`_ and `pyflakes
<https://pypi.python.org/pypi/pyflakes>`_ tools to check code for compliance.
Also, be sure to use names starting with an underscore for private helper
functions.


Writing documentation
---------------------

Make sure to document your tasks.  Documentation is primarily written in the
form of docstrings, and we tend to follow the `NumPy docstring conventions
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

To tie docstrings into the HTML documentation, edit the ``api.rst`` file
in the directory ``doc``. To generate HTML, make sure you have Sphinx,
numpydoc and Celery 3.1.10 or later::

    pip install -U sphinx numpydoc celery sphinx_bootstrap_theme

then type ``make html`` inside the ``doc`` directory. HTML will be generated
in ``doc/_build/html``.
Frequently anticipated questions
================================

* If xtas downloads optional dependencies at runtime, where will it put those?

By default, in ``~/xtas_data``. You can override this by setting the
``XTAS_DATA`` environment variable.

In addition, xtas uses NLTK extensively, and that will download resource files
to ``~/nltk_data``.


* I get ``SystemError: error return without exception set`` when starting a
  Celery worker

Check if RabbitMQ is running and xtas is properly configured to talk to it.
.. _overview:

Architecture overview
=====================

xtas can be used in one of three modes,
characterized by the following combinations of API and execution model:

* Python, standalone
* Python, distributed
* REST, distributed

In the first mode, xtas is simply a Python library of NLP tasks
that are called from your Python script or interactive interpreter,
and they are executed synchronously on your local machine.
Tasks are the functions listed in the :ref:`api`.

.. graphviz::

    digraph {
        rankdir=LR
        node [shape=rect, color=lightblue]

        "Python script" -> "xtas API"
    }

In the second mode, xtas workers run on some cluster (or a big machine),
and the local xtas library sends tasks to those workers
rather than executing them locally.
Tasks are submitted by calling tasks asynchronously::

    >>> result = xtas.tasks.guess_language.apply_async(['Welke taal zou dit zijn?'])
    >>> result
    <AsyncResult: 8b774f33-512e-41fe-90fa-33c09a7fc9c2>
    >>> result.get()
    ['nl', 0.999999995851766]

Tasks are distributed using `Celery <http://www.celeryproject.org>`_,
a Python wrapper around the RabbitMQ task queuing middleware.

.. graphviz::

    digraph {
        rankdir=LR
        node [shape=rect, color=lightblue]

        q [label="Task queue (Celery)"]

        "Python script" -> "xtas API"
        "xtas API" -> q
        q -> "xtas worker 1"
        q -> "xtas worker N"
    }

In the final mode, your code communicates with xtas workers through a
REST API that in turn communicates with the workers.

.. graphviz::

    digraph {
        rankdir=LR
        node [shape=rect, color=lightblue]

        rest [label="xtas REST API"]
        q [label="Task queue (Celery)"]

        "Program in any language" -> rest -> q
        q -> "xtas worker 1"
        q -> "xtas worker N"
    }

It follows that xtas consists of three parts:
a Python library, a worker program, and a web server program.
To :ref:`get started <setup>`, you only need the Python library and
no workers need to be running.
The web server is not need to use xtas from Python.
.. _tutorial:

xtas + Elasticsearch tutorial
=============================

Assuming you've properly configured and started xtas as described in the
section :ref:`setup`, here's how to do interesting work with it.

First, you need a document collection. If you don't have one already, download
the 20newsgroups dataset::

    curl http://qwone.com/~jason/20Newsgroups/20news-bydate.tar.gz | tar xzf -

Store the documents in Elasticsearch::

    >>> from elasticsearch import Elasticsearch
    >>> import os
    >>> from os.path import join
    >>> es = Elasticsearch()
    >>> files = (join(d, f) for d, _, fnames in os.walk('20news-bydate-train')
    ...          for f in fnames)
    ...
    >>> for i, f in enumerate(files):
    ...     body = {'text': open(f).read().decode('utf-8', errors='ignore')}
    ...     es.create(index='20news', doc_type='post', body=body, id=i)
    ...

Now, we can run named-entity recognition on the documents. Let's try it on one
document::

    >>> from xtas.tasks import es_document, stanford_ner_tag
    >>> doc = es_document('20news', 'post', 1, 'text')
    >>> tagged = stanford_ner_tag(doc)
    >>> [token for token, tag in tagged if tag == 'PERSON']
    ['Dane', 'C.', 'Butzer', 'Dane']

We just fetched the document from ES to run Stanford NER locally. That's not
the best we can do, so let's run it remotely. We can do so by running the
``stanford_ner_tag`` tasks asynchronously. First, observe that ``doc`` isn't
really the document: it's only a handle on the ES index::

    >>> doc
    {'index': '20news', 'type': 'post', 'id': 1, 'field': 'text'}

This handle can be sent over the wire to make Stanford NER run in the worker::

    >>> result = stanford_ner_tag.apply_async([doc])
    >>> result
    <AsyncResult: 44128ea3-970c-427c-80cc-2af499426b33>
    >>> [token for token, tag in result.get() if tag == 'PERSON']
    [u'Dane', u'C.', u'Butzer', u'Dane']

We have the same result, but now from a worker process. The ``result`` object
is an ``AsyncResult`` returned by Celery; see
`its documentation <http://docs.celeryproject.org/en/latest/>`_ for full
details.


Batch tasks
-----------

Some tasks require a batch of documents to work; an example is topic modeling.
Such tasks are available in the ``xtas.tasks.cluster`` package,
so named because most of the tasks can be considered a form of clustering.
Batches of documents are addressed using Elasticsearch queries,
which can be performed using xtas.
For example, to search for the word "hello" in the 20news collection::

    >>> from xtas.tasks.es import fetch_query_batch
    >>> hello = fetch_query_batch('20news', 'post',
    ...                           {'term': {'text': 'hello'}}, 'text')
    >>> len(hello)
    10

This fetches the ``'text'`` field of the documents that match the query.
(``'text'`` appears twice since you might want to match on the title,
but retrieve the body text, etc.)

Now we can fit a topic model to these document. (You need the gensim package
for this, ``pip install gensim``.) Try::

    >>> from xtas.tasks.cluster import lda
    >>> from pprint import pprint
    >>> pprint(lda(hello, 2))
    >>> pprint(lda(hello, 2))
    [[(u'and', 0.12790937338676811),
      (u'am', 0.1152182699301362),
      (u'any', 0.11363680997981712),
      (u'application', 0.10684945575768415),
      (u'algorithm', 0.10684176173416),
      (u'advance', 0.096258764014596071),
      (u'be', 0.089261082314667867),
      (u'anyone', 0.088502698282068443),
      (u'an', 0.087148958488685924),
      (u'ac', 0.068372826111416152)],
     [(u'and', 0.1150624380922104),
      (u'application', 0.11312568874979338),
      (u'anyone', 0.10728409178444523),
      (u'advance', 0.10695761539379559),
      (u'algorithm', 0.10218031622874352),
      (u'be', 0.10193138214362385),
      (u'any', 0.098412728034106681),
      (u'am', 0.097021607665297285),
      (u'an', 0.089350803758697001),
      (u'ac', 0.06867332814928688)]]

Here, the ``lda`` task returns (term, weight) pairs for two topics.
Admittedly, the topics aren't very pretty on this small set.

Of course, fetching the documents and running the topic model locally isn't
optimal use of xtas. Instead, let's set up a *chain* of tasks that runs the
query and fetches the results on a worker node, then runs the topic model
remotely as well. We'll use Celery syntax to accomplish this::

    >>> from celery import chain
    >>> fetch = fetch_query_batch.s('20news', 'post',
    ...                             {'term': {'text': 'hello'}}, 'text')
    >>> fetch_lda = chain(fetch, lda.s(k=2))    # make a chain
    >>> result = fetch_lda()                    # run the chain
    >>> pprint(result.get())                    # get results and display them
    [[[u'application', 0.11542413644453535],
      [u'am', 0.11459672375838384],
      [u'and', 0.11376035386021534],
      [u'algorithm', 0.11359529150248926],
      [u'advance', 0.10468087522675153],
      [u'be', 0.10361386971376114],
      [u'any', 0.10321250189311466],
      [u'anyone', 0.08927608350583244],
      [u'an', 0.08631073215334073],
      [u'ac', 0.055529431941575814]],
     [[u'and', 0.12924727078744289],
      [u'any', 0.1088955461011441],
      [u'anyone', 0.10640790208811389],
      [u'application', 0.10453772359410955],
      [u'advance', 0.09849716348992137],
      [u'am', 0.09774308133723099],
      [u'algorithm', 0.09546989905871668],
      [u'an', 0.09017462431916073],
      [u'be', 0.08754428312668586],
      [u'ac', 0.0814825060974741]]]

More details on creating chains can be found in the `Celery userguide
<http://celery.readthedocs.org/en/latest/userguide/canvas.html#chains>`_.


Storing results
---------------

We just saw how to run jobs remotely, fetching documents from an Elasticsearch
index. What is even more interesting is that we can also store results back to
ES, so we can use xtas as preprocessing for a semantic search engine.

We can use the ``store_single`` task to run NER on a document from the index
and store the result back, if we append it to our chain::

    >>> from celery import chain
    >>> from xtas.tasks.es import store_single
    >>> doc = es_document('20news', 'post', 3430, 'text')
    >>> ch = chain(stanford_ner_tag.s(doc, output="names"),
    ...            store_single.s('ner', doc['index'], doc['type'], doc['id']))
    >>> result = ch()
    >>> pprint(result.get())
    [[u'School of Computer Science', u'ORGANIZATION'],
     [u'McGill University Lines', u'ORGANIZATION'],
     [u'Sony', u'ORGANIZATION'],
     [u'SONY', u'ORGANIZATION'],
     [u'Invar Shadow Mask', u'ORGANIZATION'],
     [u'NEC', u'ORGANIZATION'],
     [u'NEC', u'ORGANIZATION'],
     [u'Tony', u'PERSON'],
     [u'McGill University', u'ORGANIZATION'],
     [u'Floyd', u'PERSON']]


``result.get()`` will now report the output from the NER tagger, but getting it
locally is not what we're after. The ``store_single`` task has also stored the
result back into the document, as you can verify with::

    >>> from xtas.tasks import get_all_results, get_single_result
    >>> pprint(get_all_results(doc['index'], doc['type'], doc['id']))
    {u'ner': [[u'Christopher Taylor', u'PERSON'],
              [u'Bradley University Distribution', u'ORGANIZATION'],
              [u'NHL', u'ORGANIZATION']],
    }
    >>> pprint(get_single_result('ner', doc['index'], doc['type'], doc['id']))
    [[u'Christopher Taylor', u'PERSON'],
     [u'Bradley University Distribution', u'ORGANIZATION'],
     [u'NHL', u'ORGANIZATION']]

``get_all_results`` returns all results for a document,
while ``get_single_result`` only returns the results for a specific taskname
(in our case specified as ``"ner"``).
But what if we had the task run forever ago and can't remember the tasks
that were run on a specific index?

::
    >>> pprint(get_tasks_per_index(doc['index'], doc['type']))
    set([u'ner'])

We can now actually query the xtas results.
Let's say we are interested in all documents that contain a PERSON name,
as identified by the named entity recognition::

    >>> from xtas.tasks import fetch_documents_by_task
    >>> query = {"match" : { "data" : {"query":"PERSON"}}}
    >>> pprint(fetch_documents_by_task('20news', 'post', query, 'ner', full=True))
    [[u'3430',
      {u'_id': u'3430',
       u'_index': u'20news',
       u'_score': 1.0,
       u'_source': {u'text': u"From: nittmo@camelot.bradley.edu (Christopher Taylor)\nSubject: Anyone Have Official Shorthanded Goal Totals?\nNntp-Posting-Host: camelot.bradley.edu\nOrganization: Bradley University\nDistribution: na\nLines: 4\n\nDoes anyone out there have the shorthanded goal totals of the NHL players\nfor this season?  We're trying to finish our rotisserie stats and need SHG\nto make it complete.\n\n"},
       u'_type': u'post',
       u'ner': [[u'Christopher Taylor', u'PERSON'],
                [u'Bradley University Distribution', u'ORGANIZATION'],
                [u'NHL', u'ORGANIZATION']]}]]

We now query only the results of the NER task!
The parameter ``full=True`` returns the full documents,
i.e., includes the results of the task as well.
The results of the xtas task are always stored in the ``data`` field:
make sure to take that into account when building your queries.

We can see that NHL is classified as an organisation.
What if we would like to query all occurences of NHL
to see if it is consistently classified as organisation?

::
    >>> from xtas.tasks import fetch_results_by_document
    >>> query = {"match" : { "text" : {"query":"NHL"}}}
    >>> pprint(fetch_results_by_document('20news', 'post', query, 'ner'))
    [[u'3430',
      {u'_id': u'3430',
       u'_index': u'20news',
       u'_score': 1.0,
       u'_source': {u'data': [[u'Christopher Taylor', u'PERSON'],
                              [u'Bradley University Distribution',
                               u'ORGANIZATION'],
                              [u'NHL', u'ORGANIZATION']],
                    u'timestamp': u'2014-12-18T10:57:37.259356'},
       u'_type': u'post__ner'}]]

Those two functions are convenience wrappers for the actual query function
``fetch_query_details_batch``:
here we can fire any elastic search query and get all results in a batch.
If we would want to have a simple query for all documents containing NHL,
we would say::

    >>> from xtas.tasks import fetch_query_details_batch
    >>> query = {"match" : { "text" : {"query":"NHL"}}}
    >>> pprint(fetch_query_details_batch('20news', 'post', query, True))
    [[u'3362',
      {u'_id': u'3362',
       u'_index': u'20news',
       u'_score': 0.8946465,
       u'_source': {u'text': u'From: mmilitzo@scott.skidmore.edu (matthew militzok)\nSubject: 1992 - 1993 FINAL NHL PLAYER STATS\nOrganization: Skidmo
    .... cut ....
    ]]

When multiple tasks were performed on the index,
``tasknames`` restricts the returned annotation results
to the ones in that list.

More complex queries can be built,
including relations between the between the different tasks and the documents.
Keep in mind that the type of a task is (currently)
a concatenation of the type of the original document and the taskname.
The query build to fetch a document by its task is then::

    >>> query = {'has_child': {'query': {'match': {'data': {'query': 'PERSON'}}}, 'type': 'post__ner'}}
    >>> pprint(fetch_query_details_batch('20news', 'post', query, True))
    [[u'3430',
      {u'_id': u'3430',
       u'_index': u'20news',
       u'_score': 1.0,
       u'_source': {u'text': u"From: nittmo@camelot.bradley.edu (Christopher Taylor)\nSubject: Anyone Have Official Shorthanded Goal Totals?\nNntp-Posting-Host: camelot.bradley.edu\nOrganization: Bradley University\nDistribution: na\nLines: 4\n\nDoes anyone out there have the shorthanded goal totals of the NHL players\nfor this season?  We're trying to finish our rotisserie stats and need SHG\nto make it complete.\n\n"},
       u'_type': u'post',
       u'ner': [[u'Christopher Taylor', u'PERSON'],
                [u'Bradley University Distribution', u'ORGANIZATION'],
                [u'NHL', u'ORGANIZATION']]}]]
.. xtas documentation master file, created by
   sphinx-quickstart on Wed Apr  2 14:39:31 2014.

.. raw:: html

	<div class="jumbotron">

xtas, the eXtensible Text Analysis Suite
========================================

xtas is a collection of natural language processing and text mining tools,
brought together in a single software package
with built-in distributed computing
and support for the Elasticsearch document store.

xtas functionality consists partly of wrappers for existing packages,
with automatic installation of software and data;
and partly of custom-built modules coming out of research.
Currently offered are various parsers for Dutch and English
(Alpino, CoreNLP, Frog, Semafor),
named entity recognizers (Frog, Stanford and custom-built ones),
a temporal expression tagger (Heideltime)
and a sentiment tagger based on SentiWords.

A basic installation of xtas works like a Python module.
Built-in package management and a simple, uniform interface
take away the hassle of installing, configuring and using
many existing NLP tools.

xtas's open architecture makes it possible to include custom code,
run this in a distributed fashion and have it communicate with Elasticsearch
to provide document storage and retrieval.
See :ref:`extending` for details.

.. raw:: html

		<a href="setup.html" class="btn btn-primary btn-lg">Getting started</a>
		<a href="overview.html" class="btn btn-primary btn-lg">Overview</a>
		<a href="api.html" class="btn btn-primary btn-lg">API</a>
	</div>

Contents
========

.. toctree::
   :maxdepth: 2

   setup
   tutorial
   overview
   api
   rest
   extending
   changelog
   faq


Developed by
============

.. TODO find a better way of displaying these (in the theme?)

.. raw:: html

	<div class="col-md-6">
		<a href="http://esciencecenter.nl">
			<img src="_static/logo_nlesc.png" title="Netherlands eScience Center" />
		</a>
	</div><div class="col-md-6" style="padding-top:6px;">
		<br /><br />
		<a href="http://ilps.science.uva.nl">
			<img src="_static/logo_uva.png" title="ILPS, University of Amsterdam" />
		</a>
	</div>

Changelog/release notes
=======================

3.5
---

* Stanford CoreNLP is now downloaded automatically.

3.4
---

* Added: emotion classifier for film reviews (and maybe other English text).
* Documentation updates.
* Latent dirichlet allocation (LDA, topic modeling) is now performed using
  scikit-learn instead of Gensim. This removes one dependency.
* xtas now requires NLTK 3.1. The previously required version, 3.0a2,
  conflicted with some versions of six, causing cryptic ImportErrors when
  loading xtas in some setups.

3.3
---

* New fast NER tagger for Dutch. Contributed by Daan Odijk of UvA, and based on the CoNLL dataset.
* Fixed a bug that prevented parsimonious language models from returning top-k words for large k. Patch contributed by Sicco van Sas.
* Various minor fixes.

3.2
---

The xtas REST server can now run within an application container such as
uwsgi. By default, it uses Tornado for improved throughput.

A wrapper for the Heideltime temporal tagger has been added.

3.1
---

xtas 3.1 has improved integration with Elasticsearch compared to its
predecessor, storing processing results as child documents. (Previously, they
were stored as fields.) Documentation has also been greatly improved to
demonstrate the use of xtas as a preprocessing engine for semantic search.

Other new features are:

* Stemming support for English, German, Norwegian, Italian, Dutch,
  Portuguese, French, Swedish (via PyStemmer).
* TCP port for communication with Frog (POS tagger/NER/parser for Dutch) is
  now configurable.
* Various bugfixes and small optimizations.

3.0
---

xtas 3.0 is a complete rewrite of the University of Amsterdam's xTAS system
(retroactively, xtas 2) and its predecessor, Fietstas (xtas 1). No attempt
has been made to retain backward compatibility, so that code could be
redesigned in a clean, simple, scalable way.

Features introduced in xtas 3 include:

* Communication with Elasticsearch
* REST API to single-document tasks
* Python API for both synchronous and asynchronous execution (with Celery)

And the following NLP tasks:

* Tokenization, based on NLTK
* Named-entity recognition with Stanford NER (``stanford_ner_tag``)
* Wrapper for the Frog Dutch POS tagger/NER tagger/dependency parser
* Wrapper for the Alpino parser for Dutch
* Language guessing
* English lemmatization with morphy
* Movie review polarity classification
* English POS tagging with NLTK
* Wrapper for the Semafor semantic parser
* Word-level sentiment polarity tagging with SentiWords
* Document clustering with :math:`k`-means (``kmeans``, ``big_kmeans``)
* Topic modeling with latent semantic analysis (LSA), latent dirichlet
  allocation (LDA), based on scikit-learn and Gensim
* Parsimonious language modeling, useful for discriminative word clouds

Preliminary versions of 3.0 were named 2.99.
.. _api:

API reference
=============

xtas.core
---------

.. automodule:: xtas.core

.. autofunction:: configure


xtas.tasks.single
-----------------

.. automodule:: xtas.tasks.single

.. It seems we have to list all the tasks here. automodule doesn't pick them
   up, probably because of the decorator.

.. autotask:: alpino
.. autotask:: corenlp
.. autotask:: corenlp_lemmatize
.. autotask:: dbpedia_spotlight
.. autotask:: frog
.. autotask:: guess_language
.. autotask:: morphy
.. autotask:: movie_review_polarity
.. autotask:: pos_tag
.. autotask:: semafor
.. autotask:: semanticize
.. autotask:: sentiwords_tag
.. autotask:: stanford_ner_tag
.. autotask:: stem_snowball
.. autotask:: tokenize
.. autotask:: untokenize


xtas.tasks.cluster
------------------

.. automodule:: xtas.tasks.cluster

.. autotask:: big_kmeans
.. autotask:: kmeans
.. autotask:: lda
.. autotask:: lsa
.. autotask:: parsimonious_wordcloud
