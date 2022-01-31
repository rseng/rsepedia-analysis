# Changelog

## v0.7.5 2022-01-17

- Drop support for Python 3.6 [[#108]](https://github.com/aiidateam/kiwipy/pull/108)
- `RmqCommunicator`: add the `server_properties` property [[#107]](https://github.com/aiidateam/kiwipy/pull/107)
- Expose `aio_pika.Connection.add_close_callback` [[#104]](https://github.com/aiidateam/kiwipy/pull/104)

## v0.7.4 2021-03-02

- ♻️ REFACTOR: BroadcastFilter to extract filter conditions into a separate `is_filtered` method.

## v0.7.3 2021-02-24

- 👌 IMPROVE: Add debug logging for sending task/rpc/broadcast to RMQ.
- 👌 IMPROVE: Close created asyncio loop on RmqThreadCommunicator.close

## v0.7.2 2021-02-11

- 🐛 FIX: an aio-pika deprecation, to use async context managers when processing messages.

## v0.7.1

The default task message TTL setting was changed in `v0.5.4` but this breaks existing queues since RabbitMQ does not allow changing these parameters on existing queues.
Therefore the change was reverted which was released in `v0.5.5`.
However, since that was a patch release, it had not been merged back to `v0.6.0` as well, which therefore from the problem described.
The same revert is applied in this release to restore original functionality.

### Changes
- Revert "Increase the default TTL for task messages" [[#93]](https://github.com/aiidateam/kiwipy/pull/93)


## v0.7.0

### Changes
- Add support for Python 3.9 [[#87]](https://github.com/aiidateam/kiwipy/pull/87)
- Drop support for Python 3.5 [[#89]](https://github.com/aiidateam/kiwipy/pull/89)
- Replace old format string interpolation with f-strings [[#90]](https://github.com/aiidateam/kiwipy/pull/90)

### Bug fixes
- Fix warning caused by excepted task and no reply [[#83]](https://github.com/aiidateam/kiwipy/pull/83)

### Dependencies
- Dependencies: update upper limit requirement for `pytray>=0.2.2,<0.4.0` [[#80]](https://github.com/aiidateam/kiwipy/pull/80)
- Dependencies: update requirement `pytest-asyncio~=0.12` [[#82]](https://github.com/aiidateam/kiwipy/pull/82)
---
title: 'kiwiPy: Robust, high-volume, messaging for big-data and computational science workflows'
tags:
  - python
  - high-throughput computing
  - message broker
  - scientific workflows
  - fault-tolerant computing
authors:
  - name: Martin Uhrin
    orcid: 0000-0001-6902-1289
    affiliation: "1, 2, 3"
  - name: Sebastiaan P. Huber
    orchid: 0000-0001-5845-8880
    affiliation: "2, 3"
affiliations:
 - name: Department of Energy Conversion and Storage, Technical University of Denmark, 2800 Kgs. Lyngby, Denmark
   index: 1
 - name: National Centre for Computational Design and Discovery of Novel Materials (MARVEL), École Polytechnique Fédérale de Lausanne, CH-1015 Lausanne, Switzerland
   index: 2
 - name: Theory and Simulation of Materials (THEOS), Faculté des Sciences et Techniques de l’Ingénieur, École Polytechnique Fédérale de Lausanne, CH-1015 Lausanne, Switzerland
   index: 3
date: 19 April 2020
bibliography: paper.bib
---

# Summary

The computational sciences have seen a huge increase in the use of high-throughput, automated, workflows over the course of the last two decades or so.  Focusing on just our domain of computational materials science, there have been several large scale initiatives to provide high-quality results from standardised calculations [@Landis2012; @Curtarolo2012; @Jain2013; @Saal2013; @Draxl2019a; @Talirz2020].  Almost all of these repositories are populated using results from high-throughput quantum mechanical calculations that rely on workflow frameworks [@Jain2015a; @Mortensen2020], including our own [AiiDA](http://www.aiida.net/) [@Pizzi2016; @Huber2020] which powers the [Materials Cloud](https://www.materialscloud.org/).  One of the many challenges for such frameworks is maximising fault-tolerance whilst simultaneously maintaining high-throughput, often across several systems (typically the client launching the tasks, the supercomputer carrying out the computations, and the server hosting the database).

On the software level, these problems are perhaps best addressed by using messaging brokers that take responsibility for guaranteeing the durability (or persistence) and atomicity of messages and often enable event-based communication.
Indeed, solutions such as [RabbitMQ](https://www.rabbitmq.com/), see widespread adoption in industry. However, adoption in academia has been more limited, with home-made queue data structures, race condition susceptible locks and polling based solutions being commonplace.
This is likely due to message brokers typically having complex APIs (which reflect the non-trivial nature of the underlying protocol) as well as the lack of familiarity with event-based programming in general within the community.
[KiwiPy](https://kiwipy.readthedocs.io/en/latest/) was designed specifically to address both these issues, by providing a tool that enables building robust, event-based systems with an interface that is as simple as possible.

In kiwiPy, all messages are saved to disk by RabbitMQ, meaning that any or all systems involved in a workflow, including the broker, can be shut down (abruptly or gracefully), and the previous state can be recreated allowing the workflow to continue when the necessary resources are brought back online.  This is especially important for long-running HPC jobs that may take days or weeks.  This feature differentiates kiwiPy from protocols such as MPI, [ZeroMQ](https://zeromq.org/) or libraries such as [Dask](https://docs.dask.org/en/latest/), which do not persist their state.

A number of libraries for interacting directly with message brokers exist, including [Pika](https://pika.readthedocs.io/en/stable/), [aio-pika](https://aio-pika.readthedocs.io/en/latest/index.html), [py-amqp](https://barryp.org/software/py-amqplib), [kombu](https://github.com/celery/kombu) and others.  These tend to be rather low-level and focus on channels, exchanges, routing, sockets and so on.   A comparison of the difference in focus between Pika and kiwiPy can be found in the [documentation](https://kiwipy.readthedocs.io/en/latest/).  At the other end of the spectrum are libraries such as [Celery](https://docs.celeryproject.org/en/stable/getting-started/introduction.html), [RQ](https://python-rq.org/) and [others](https://www.fullstackpython.com/task-queues.html) that provide task queues and libraries such as [RPyC](https://rpyc.readthedocs.io/en/latest/), [Spyne](http://spyne.io), [Python-JRPC](https://github.com/alex-sherman/python-jrpc) and [others](https://stackoverflow.com/questions/1879971/what-is-the-current-choice-for-doing-rpc-in-python) that enable remote procedure calls.  In contrast, kiwiPy brings together three commonly used message types (task queues, Remote Procedure Calls (RPCs), and, broadcasts) in a single interface.

All messaging in kiwiPy is done in the `Communicator` class, which can be trivially constructed by providing a URI string pointing to the RabbitMQ server.  By default, kiwiPy creates a separate communication thread that the user never sees, allowing them to interact with the communicator using familiar Python syntax, without the need to be familiar with either coroutines or multithreading.  This has the additional advantage that kiwiPy will maintain heartbeats (a periodic check to make sure the connection is still alive) with the server whilst the user code can be doing other things. Heartbeats are an essential part of RabbitMQ's fault tolerance; two missed checks will automatically trigger the message to be requeued to be picked up by another client.

To demonstrate some of the possible usage scenarios, we briefly outline the way kiwiPy is used in AiiDA.  AiiDA, amongst other things, manages the execution of complex workflows made up of processes that may have checkpoints.

Task queues
-----------

As is common for high-throughput workflow engines, AiiDA maintains a task queue to which processes are submitted (typically from the user's workstation).  These tasks are then consumed by multiple daemon processes (which may also be on the user's workstation or remote) and will only be removed from the task queue once they have been acknowledged to be completed by the consumer.  The daemon can be gracefully or abruptly shut down and no task will be lost, since the task will simply be requeued by the broker once it notices that the consumer has died.  Furthermore, there are no worries about race conditions between multiple daemon processes, since the task queue is guaranteed to only distribute each task to, at most, one consumer at a time.

Remote Procedure Calls
----------------------

These are used to control live processes.  Each process has a unique identifier and can be sent a `pause`, `play`, or `kill` message, the response to which is optionally sent back to the initiator to indicate success or something else.

Broadcasts
----------

These currently serve two purposes: sending `pause`, `play`, or `kill` messages to all processes at once by broadcasting the relevant message, and controlling the flow between processes.  If a parent process is waiting for a child to complete, it will be informed of this via a broadcast message from the child saying that its execution has terminated. This enables decoupling, as the child need not know about the existence of the parent.


Together these three message types allow AiiDA to implement a highly-decoupled, distributed, yet, reactive system that has proven to be scalable from individual laptops to workstations, driving simulations on high-performance supercomputers with workflows consisting of varying durations, ranging from milliseconds up to multiple days or weeks.

It is our hope that by lowering the barriers to adoption, kiwiPy will bring the benefits of industry grade message brokers to academia and beyond, ultimately making robust scientific software easier to write and maintain.


# Acknowledgements

We would like to thank Giovanni Pizzi, Nicola Marzari, and the AiiDA team for their continuous coordination and development of the project.  We also thank Jason Yu for contributing the first version of the documentation.
This work is supported by the MARVEL National Centre for Competency in Research funded by the Swiss National Science Foundation (grant agreement ID 51NF40-182892) and the European Materials Modelling Council-CSA funded by the European Union‘s Horizon 2020 research and innovation programme under Grant Agreement No 723867.

# References
.. _AiiDA: https://www.aiida.net
.. _rmq tutorial: https://www.rabbitmq.com/getstarted.html
.. _documentation: https://kiwipy.readthedocs.io/en/latest/index.html


kiwiPy
======

.. image:: docs/source/_static/logo.svg
   :height: 64px
   :width: 64px
   :alt: kiwiPy

.. image:: https://codecov.io/gh/aiidateam/kiwipy/branch/develop/graph/badge.svg
    :target: https://codecov.io/gh/aiidateam/kiwipy
    :alt: Coveralls

.. image:: https://github.com/aiidateam/kiwipy/workflows/continuous-integration/badge.svg
    :target: https://github.com/aiidateam/kiwipy/actions?query=workflow%3Acontinuous-integration
    :alt: Github Actions

.. image:: https://img.shields.io/pypi/v/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/
    :alt: Latest Version

.. image:: https://img.shields.io/pypi/pyversions/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/

.. image:: https://img.shields.io/pypi/l/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/

.. image:: https://joss.theoj.org/papers/10.21105/joss.02351/status.svg
   :target: https://doi.org/10.21105/joss.02351



`kiwiPy`_ is a library that makes remote messaging using RabbitMQ (and possibly other message brokers) EASY.  It was
designed to support high-throughput workflows in big-data and computational science settings and is currently used
by `AiiDA`_ for computational materials research around the world.  That said, kiwiPy is entirely general and can
be used anywhere where high-throughput and robust messaging are needed.

Here's what you get:

* RPC
* Broadcast (with filters)
* Task queue messages

Let's dive in, with some examples taken from the `rmq tutorial`_.  To see more detail head over to the `documentation`_.

RPC
---

The client:

.. code-block:: python

    import kiwipy

    with kiwipy.connect('amqp://localhost') as comm:
        # Send an RPC message
        print(" [x] Requesting fib(30)")
        response = comm.rpc_send('fib', 30).result()
        print((" [.] Got %r" % response))

`(rmq_rpc_client.py source) <https://raw.githubusercontent.com/aiidateam/kiwipy/develop/examples/rmq_rpc_client.py>`_


The server:

.. code-block:: python

    import threading
    import kiwipy

    def fib(comm, num):
        if num == 0:
            return 0
        if num == 1:
            return 1

        return fib(comm, num - 1) + fib(comm, num - 2)

    with kiwipy.connect('amqp://127.0.0.1') as comm:
        # Register an RPC subscriber with the name 'fib'
        comm.add_rpc_subscriber(fib, 'fib')
        # Now wait indefinitely for fibonacci calls
        threading.Event().wait()

`(rmq_rpc_server.py source) <https://raw.githubusercontent.com/aiidateam/kiwipy/develop/examples/rmq_rpc_server.py>`_


Worker
------

Create a new task:

.. code-block:: python

    import sys
    import kiwipy

    message = ' '.join(sys.argv[1:]) or "Hello World!"

    with rmq.connect('amqp://localhost') as comm:
        comm.task_send(message)

`(rmq_new_task.py source) <https://raw.githubusercontent.com/aiidateam/kiwipy/develop/examples/rmq_new_task.py>`_


And the worker:

.. code-block:: python

    import time
    import threading
    import kiwipy

    print(' [*] Waiting for messages. To exit press CTRL+C')


    def callback(_comm, task):
        print((" [x] Received %r" % task))
        time.sleep(task.count(b'.'))
        print(" [x] Done")


    try:
        with kiwipy.connect('amqp://localhost') as comm:
            comm.add_task_subscriber(callback)
            threading.Event().wait()
    except KeyboardInterrupt:
        pass

`(rmq_worker.py source) <https://raw.githubusercontent.com/aiidateam/kiwipy/develop/examples/rmq_worker.py>`_

Citing
======

If you use kiwiPy directly or indirectly (e.g. by using `AiiDA`_) then please cite:

Uhrin, M., & Huber, S. P. (2020). kiwiPy : Robust , high-volume , messaging for big-data and computational science workflows, 5, 4–6. http://doi.org/10.21105/joss.02351

This helps us to keep making community software.

Versioning
==========

This software follows `Semantic Versioning`_

Contributing
============

Want a new feature? Found a bug? Want to contribute more documentation or a translation perhaps?

Help is always welcome, get started with the `contributing guide <https://github.com/aiidateam/kiwipy/wiki/Contributing>`__.

.. _Semantic Versioning: http://semver.org/

Development
===========

This package utilises `tox <https://tox.readthedocs.io>`__ for unit test automation, and `pre-commit <https://pre-commit.com>`__ for code style formatting and test automation.

To install these development dependencies:

.. code-block:: bash

    pip install tox pre-commit

To run the unit tests:

.. code-block:: bash

    tox

For the ``rmq`` tests you will require a running instance of RabbitMQ.
One way to achieve this is using Docker and launching ``test/rmq/docker-compose.yml``.

To run the pre-commit tests:

.. code-block:: bash

    pre-commit run --all

To build the documentation:

.. code-block:: bash

    tox -e docs-clean

Changes should be submitted as Pull Requests (PRs) to the ``develop`` branch.

Publishing Releases
===================

1. Create a release PR/commit to the ``develop`` branch, updating ``kiwipy/version.py`` and ``CHANGELOG.md``.
2. Fast-forward merge `develop` into the `master` branch
3. Create a release on GitHub (https://github.com/aiidateam/kiwipy/releases/new), pointing to the release commit on `master`, named ``v.X.Y.Z`` (identical to version in ``kiwipy/version.py``)
4. This will trigger the ``continuous-deployment`` GitHub workflow which, if all tests pass, will publish the package to PyPi. Check this has successfully completed in the GitHub Actions tab (https://github.com/aiidateam/kiwipy/actions).

(if the release fails, delete the release and tag)
.. _tasks example: examples/tasks.ipynb
.. _rpc example: examples/rpc.ipynb
.. _broadcast example: examples/broadcast.ipynb



Concepts
========

Here we introduce the meaning of task, rpc, and broadcast and their usage scenarios, since the meaning of a the task in
``kiwiPy`` differs slightly from the one in RabbitMQ's tutorial.

Task
----

Tasks are one to many messages.  This means that you sent out a task to a queue and there can be zero or more workers
attached one of which will process the task when it is ready.  The result of the task can optionally be delivered to the
sender.  See the `tasks example`_ to see this in action.

RPC
---

A remote procedure is one-to-one.  This is used when you want to call a particular remote function/method and (usually)
expect an immediate response. For example imagine asking a remote process to pause itself.  Here you would make a RPC
and wait to get confirmation that it has, indeed, paused.  See the `rpc example`_ to see this in action.


Broadcast
---------

One sender to zero or more consumers.  These are fire and forget calls that broadcast a message to anyone who is
listening.  Consumers may optionally apply a filter to only receive messages that match some criteria.
See the `broadcast example`_ to see this in action.


Coroutines vs Threads
=====================

KiwiPy's RabbitMQ backend has both coroutine and threaded version of the Communicator.  If you're in a coroutine
environment you may prefer to use the :class:`~kiwipy.rmq.communicator.RmqCommunicator` which uses coroutines and
:class:`asyncio.futures.Future` s while if you don't have an event loop just stick with the :class:`~kiwipy.rmq.threadcomms.RmqThreadCommunicator` which
runs an event loop in a separate thread to allow it to do asynchronous communication while your application logic gets on
with other stuff.

If all of this means nothing to you then stick with the :class:`~kiwipy.rmq.threadcomms.RmqThreadCommunicator`.

Examples
========

Here are some examples of the basic functionality.  Keep in mind that what we should is running on one system but
usually the communicator would be used on multiple, distributed, clients.

.. toctree::
    :maxdepth: 1
    :glob:

    examples/broadcast.ipynb
    examples/tasks.ipynb
    examples/rpc.ipynb
.. _pika: https://pika.readthedocs.io/en/stable/
.. _aio-pika: https://pika.readthedocs.io/en/stable/
.. _msgpack: https://github.com/msgpack/msgpack-python


Performance
===========

KiwiPy is designed primarily with ease-of-use and robustness in mind with less of a focus on performance.
Never the less, it is expected that depending on the configuration and network conditions kiwPy should be able to handle hundreds to thousands of messages per second.

The chart below shows benchmarks comparing kiwiPy to `pika`_ using the following configuration:

* Sending 100 task messages consisting of 32 to 2048 bytes in size (default 1024)
* Messages are sent and _then_ received (i.e. serialised) to be able to show split timings
* Everything running locally (server and client) on Dell XPS 15 (Core i9-9980H @ 2.3 x 16, 1TB SSD HDD), Ubuntu 20.04
* Using `msgpack`_ for encoding and decoding of messages
* Mandatory routing enabled meaning that the code will wait for messages to be routed by the broker (or raise an exception)
* Full code can be found in ``test/rmq/bench/test_benchmarks.py``

Tests are split into ``send``, ``get`` and ``send_get`` to show a breakdown of the full process of sending and receiving a message.

The ``send_get`` test is carried out for different message size (in bytes) as indicated in square brackets.


.. image:: /_static/bench.svg
  :width: 640
  :alt: Performance comparison of kiwiPy vs pika


As we can see, kiwiPy is faster in getting tasks (``test_kiwi_get``: 21.7 ms vs ``test_pika__get``: 50.6 ms) but significantly slower in sending them out (``test_kiwi_send``: 237.7 ms vs ``test_pika_send``: 31.5 ms) and this results in overall (send/get) times that are roughly 3.5x (or 200ms) slower than that of pika.
This may be due to the fact that kiwiPy is running a separate thread to perform all communications meaning that when a client makes a request to send a message (from the main thread) kiwiPy synchronises with the communications thread and waits until the message has been send before returning control the user.  This makes it easier to write correct code when using kiwiPy as any errors encountered in sending will be raised immediately as exceptions, as opposed to the user having to explicitly check a future.  Incoming messages, on the other hand, are delivered directly on the communications thread and it is the users responsibility to synchronise with their main thread if they need, hence we see no performance loss.

To run these benchmarks yourself, just clone the source, then install and run like so

.. code-block:: shell

    pip install -e .[rmq,tests]
    pytest test


At the end of the tests you should see a summary of all the timings.

.. _examples: examples.rst


Installation
============

Python
------

KiwiPy supports Python versions 3.7 and above.

RabbitMQ
--------

KiwiPy depends on RabbitMQ as the message broker.
On Ubuntu this is as simple as:

.. code-block:: shell

    apt install rabbitmq

For other platforms refer to `Downloading and Installing RabbitMQ <https://www.rabbitmq.com/download.html>`_

Basic Installation
------------------

    $ pip install kiwipy[rmq]

Now you're ready to run `examples`_!

Building from Source
--------------------

In order to develop kiwipy it's best to install kiwipy in ``editable`` mode. This allows changes you
make to kiwipy to be reflected immediately in your runtime environment.

First, clone the source:

.. code-block:: shell

   $ git clone https://github.com/aiidateam/kiwipy.git

Then, create and activate a virtualenv:

.. code-block:: shell

    virtualenv venv
    . venv/bin/activate
    pip install -e "kiwipy[rmq,pre-commit,tests]"

To run the tests, make sure the RabbitMQ server is up and running (see the RabbitMQ documentation on how to accomplish and/or verify this) and type:

.. code-block:: shell

    pytest test

.. _Pika: https://pika.readthedocs.io/en/stable/
.. _RabbitMQ tutorials: https://www.rabbitmq.com/getstarted.html
.. _work queues: https://www.rabbitmq.com/tutorials/tutorial-two-python.html
.. _topics: https://www.rabbitmq.com/tutorials/tutorial-five-python.html
.. _RPC: https://www.rabbitmq.com/tutorials/tutorial-six-python.html
.. _aio-pika: https://aio-pika.readthedocs.io/en/latest/index.html

Coming from Pika
================

KiwiPy comes natively with a RabbitMQ communicator (others can be added by extending the :py:class:`~kiwipy.Communicator` interface) and thus it may be useful to see how to achieve the same things in `Pika`_ (the standard library used in the `RabbitMQ tutorials`_) and kiwiPy.  This also shows some of the differences, particularly where kiwiPy is less verbose by setting sensible defaults for the user.


Work Queues
-----------

Let's start with RabbitMQ's `work queues`_ example.


Pika
++++

The code for sending a task:

.. literalinclude:: examples/pika/new_task.py
    :language: python


The code for running a worker:

.. literalinclude:: examples/pika/worker.py
    :language: python


KiwiPy
++++++

And now, in kiwiPy.  The code for sending the same task:

.. literalinclude:: examples/kiwipy/new_task.py
    :language: python

Here, compared to the Pika snippet, we see that we don't have to decide or know about an ``exchange``, ``routing_key``, ``properties``, ``channel`` or durability of the message.
Task queues in kiwiPy are always durable and the corresponding messages persistent.
The ``routing_key`` is simply the name of the queue that the user declared and an ``exchange`` is selected automatically for the user (set to ``kiwipy.tasks`` by default).
The ``exchange`` and some other defaults can be changed when constructing a :py:class:`~kiwipy.rmq.threadcomms.RmqThreadCommunicator`.

Here we have explicitly created a task queue called ``task_queue``, however, even this has a default version if the user doesn't need multiple queues.  In which case the code would simply be::

    with kiwipy.connect('amqp://localhost') as comm:
        comm.task_send(message)
        print(" [x] Sent %r" % message)

Now to run a worker it's just:

.. literalinclude:: examples/kiwipy/worker.py
    :language: python

This is fairly straightforward and the only decision the user has to make is about the prefetch count (the number of tasks and worker can be working on simultaneously).


Topics
------

Moving on to the `topics`_ example (known as broadcasts with filters in kiwiPy).


Pika
++++

The Pika code to emit a log topic is:

.. literalinclude:: examples/pika/emit_log_topic.py
    :language: python

And to receive, the code is:

.. literalinclude:: examples/pika/receive_logs_topic.py
    :language: python


KiwiPy
++++++

Emitting in kiwiPy:

.. literalinclude:: examples/kiwipy/emit_log_topic.py
    :language: python

As we've come to see in kiwiPy there's no need to worry about the channel, exchange, exchange type while the routing key is known as the subject of the message.

And to receive the log:

.. literalinclude:: examples/kiwipy/receive_logs_topic.py
    :language: python

Here we can see that the filtering is performed by using a :py:class:`~kiwipy.BroadcastFilter` which can match a string, that optionally includes wildcards ``*``.


RPC
---

Finally, let's end with the `RPC`_ example.


Pika
++++

The code to run an RPC server is:

.. literalinclude:: examples/pika/rpc_server.py
    :language: python


The client side is:

.. literalinclude:: examples/pika/rpc_client.py
    :language: python


KiwiPy
++++++

Now, in kiwiPy the server becomes:

.. literalinclude:: examples/kiwipy/rpc_server.py
    :language: python

As usual, the routing of messages is automatically set up by kiwiPy and the ``identifier`` of the RPC subscriber is the only thing needed to link the client and server.

The client code is:

.. literalinclude:: examples/kiwipy/rpc_client.py
    :language: python

Here, we send the RPC request to the identifier we used before.  The ``rpc_send`` call gives us back a future from which we can get the result once it has been received from the server.


Roundup
-------

As these examples demonstrate, kiwiPy tries to use sensible defaults and keep many of the routing details hidden from the user.
This has the advantage of bringing the power of a dedicated message broker such as RabbitMQ to less expert users.
Additionally, it allows messaging to be done using a simple interface which can be implemented by other messaging protocols, allowing clients to be protocol agnostic.  The trade-off, of course, is reduced control over the finer details of how the messaging is performed.
For users that need such fine-grained control, `pika`_ and the excellent `aio-pika`_ may be better alternatives.
.. kiwipy documentation master file, created by
   sphinx-quickstart on Tue May 14 13:41:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _kiwiPy: https://github.com/aiidateam/kiwipy
.. _AiiDA: https://www.aiida.net
.. _Pika: https://pika.readthedocs.io/en/stable/
.. _examples: examples.rst
.. _concepts: concepts.rst
.. _coming from pika: pika-comparison.rst
.. _installation: installation.rst
.. _performance: performance.rst
.. _API documentation: apidoc.rst


Welcome to kiwiPy's documentation!
==================================

.. image:: https://codecov.io/gh/aiidateam/kiwipy/branch/develop/graph/badge.svg
    :target: https://codecov.io/gh/aiidateam/kiwipy
    :alt: Coveralls

.. image:: https://github.com/aiidateam/kiwipy/workflows/continuous-integration/badge.svg
    :target: https://github.com/aiidateam/kiwipy/actions?query=workflow%3Acontinuous-integration
    :alt: Github Actions

.. image:: https://img.shields.io/pypi/v/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/
    :alt: Latest Version

.. image:: https://img.shields.io/pypi/pyversions/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/

.. image:: https://img.shields.io/pypi/l/kiwipy.svg
    :target: https://pypi.python.org/pypi/kiwipy/

.. image:: https://joss.theoj.org/papers/10.21105/joss.02351/status.svg
   :target: https://doi.org/10.21105/joss.02351


`kiwiPy`_ is a library that makes remote messaging using RabbitMQ (and possibly other message brokers) EASY.  It was
designed to support high-throughput workflows in big-data and computational science settings and is currently used
by `AiiDA`_ for computational materials research around the world.  That said, kiwiPy is entirely general and can
be used anywhere where high-throughput and robust messaging are needed.


Features
++++++++

* Support for `1000s of messages per second <performance_>`__
* Highly robust - no loss of messages on connection interruptions, etc., as messages are automatically persisted to disk
* Generic communicator interface with native support for RabbitMQ
* Supports task queues, broadcasts and RPC
* Support for both thread and coroutine based communication
* Python 3.7+ compatible.


Getting Started
+++++++++++++++

* To install kiwiPy follow the instructions in the `installation`_ section
* After you have successfully installed kiwipy, give in to some of the `examples`_ to see what kiwiPy can do.
* The design concepts behind kiwiPy can be found in `concepts`_ section
* If you're already familiar with `Pika`_ you might find the `coming from pika`_ section useful
* Finally check out the complete `API documentation`_

.. admonition:: Development Contributions
   :class: note

   Want a new feature? Found a bug? Want to contribute more documentation or a translation perhaps?
   Help is always welcome, get started with the `contributing guide <https://github.com/aiidateam/kiwipy/wiki/Contributing>`__.



Table Of Contents
+++++++++++++++++

.. toctree::
   :maxdepth: 2

   installation
   concepts
   examples
   pika-comparison
   performance
   API Reference <apidoc/kiwipy>


Citing
++++++

If you use kiwiPy directly or indirectly (e.g. by using `AiiDA`_) then please cite:

Uhrin, M., & Huber, S. P. (2020). kiwiPy : Robust , high-volume , messaging for big-data and computational science workflows, 5, 4–6. http://doi.org/10.21105/joss.02351

This helps us to keep making community software.

Versioning
++++++++++

This software follows `Semantic Versioning`_

.. _Semantic Versioning: http://semver.org/
