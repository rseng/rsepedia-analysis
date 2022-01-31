Generated using

```bash
ssh-keygen -t rsa
```

id1_rsa has no passphrase
id2_rsa has passphrase 'kingfisher'
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

1. use the search functionality `here <https://github.com/MD-Studio/cerulean/issues>`_ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

1. use the search functionality `here <https://github.com/MD-Studio/cerulean/issues>`_ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_ and `here <https://help.github.com/articles/syncing-a-fork/>`_);
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the Cerulean repository on GitHub;
1. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`_.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.3.7
*****

Fixed
-----

* Symlink resolving issue
* Tooling improvements


0.3.6
*****

Added
-----

* Support for redirecting scheduler/system output
* Path.remove() method

0.3.5
*****

Fixed
-----

* Detect a closed SSH socket and auto-reconnect.

0.3.4
*****

Fixed
-----

* Directory permissions when using mkdir(). This is a security issue, and you
  should upgrade as soon as possible.

Added
-----

* Equality comparison of FileSystems and Terminals, which makes equality
  comparison for Paths work better as well.
* Add `interval` parameter to scheduler.wait(), and improve default behaviour.

0.3.3
*****

Fixed
-----

* Copy silently ignored missing file, now raises FileNotFoundError

Added
-----

* FileSystem.root() to get a Path for the root
* Path.__repr__() for better debugging output

0.3.2
*****

Fixed
-----

* Various small things

Added
-----

* Support for Slurm 18.08 (worked already, now also part of the tests)
* Add command prefix for schedulers
* Add support for WebDAV

0.3.1
*****

Fixed
-----

* Extraneous slashes in paths
* Properly handle errors on Slurm submission
* Leftover print statement in Path.walk


0.3.0
*****

Added
-----

* New copy_permissions option to copy()
* New callback option to copy()
* New Path.walk() method

Fixed
-----

* Add missing EntryType and Permission classes to API
* SFTP-to-SFTP copy deadlock


0.2.0
*****

Added
-----

* Path.write_text()
* Scheduler.submit_job() is now Scheduler.submit()
* Scheduler.wait()

Fixed
-----

* Bugs in copy()
* Documentation for JobDescription


0.1.0
*****

Initial release
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

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

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at l.veen@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
.. image:: https://readthedocs.org/projects/cerulean/badge/?version=develop
    :target: https://cerulean.readthedocs.io/en/latest/?badge=develop
    :alt: Documentation Build Status

.. image:: https://github.com/MD-Studio/cerulean/actions/workflows/ci.yaml/badge.svg?branch=develop
    :target: https://github.com/MD-Studio/cerulean/actions
    :alt: Build Status

.. image:: https://app.codacy.com/project/badge/Grade/4909f6a7c0d94cc3a4b3c83230a20248
    :target: https://www.codacy.com/gh/MD-Studio/cerulean/dashboard
    :alt: Codacy Grade

.. image:: https://app.codacy.com/project/badge/Coverage/4909f6a7c0d94cc3a4b3c83230a20248
    :target: https://www.codacy.com/gh/MD-Studio/cerulean/dashboard
    :alt: Code Coverage

.. image:: https://requires.io/github/MD-Studio/cerulean/requirements.svg?branch=develop
    :target: https://requires.io/github/MD-Studio/cerulean/requirements/?branch=develop
    :alt: Requirements Status

################################################################################
Cerulean
################################################################################

Cerulean is a Python 3 library for talking to HPC clusters and supercomputers.
It lets you copy files between local and SFTP filesystems using a
``pathlib``-like API, it lets you start processes locally and remotely via SSH,
and it lets you submit jobs to schedulers such as Slurm and Torque/PBS.

Documentation and Help
**********************

Cerulean can be installed as usual using pip:

`pip install cerulean`

Instructions on how to use Cerulean can be found in `the Cerulean documentation
<https://cerulean.readthedocs.io>`_.

Code of Conduct
---------------

Before we get to asking questions and reporting bugs, we'd like to point out
that this project is governed by a code of conduct, as described in
CODE_OF_CONDUCT.rst, and we expect you to adhere to it. Please be nice to your
fellow humans.

Questions
---------

If you have a question that the documentation does not answer for you, then you
have found a bug in the documentation. We'd love to fix it, but we need a bit of
help from you to do so. Please do the following:

#. use the search functionality `here
   <https://github.com/MD-Studio/cerulean/issues>`__
   to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new
   issue;
#. apply the "Question" label; apply other labels when relevant.

We'll answer your question, and improve the documentation where necessary.
Thanks!

Bugs
----

Like most software, Cerulean is made by humans, and we make mistakes. If you
think you've found a bug in Cerulean, please let us know! Reporting bugs goes as
follows.

#. Use the search functionality `here
   <https://github.com/yatiml/yatiml/issues>`_
   to see if someone already filed the same issue.
#. If your issue search did not yield any relevant results, make a new issue.
   Please explain:
   - what you were trying to achieve,
   - what you did to make that happen,
   - what you expected the result to be,
   - what happened instead.
   It really helps to have the actual code for a simple example that
   demonstrates the issue, but excerpts and error messages and a
   description are welcome too.
#. Finally, apply any relevant labels to the newly created issue.

With that, we should be able to fix the problem, but we may ask for some more
information if we can't figure it out right away.

Development
-----------

More information for Cerulean developers may be found in `the Cerulean
documentation <https://cerulean.readthedocs.io>`_.

License
*******

Copyright (c) 2018, The Netherlands eScience Center and VU University Amsterdam

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
Testing Cerulean
================

Testing a library that is intended to work with external hardware resources is a
bit tricky. Of course you can run unit tests, but you want to test the whole
thing as well. Cerulean solves this by using a set of Docker containers that run
simulated clusters running a variety of schedulers. The containers also contain
an SSH server, so everything can be tested.

Running the Cerulean tests is done by calling ``python3 setup.py test`` from the
root directory. When you do this, the following happens:

- ``setup.py`` runs ``pytest``
- ``pytest`` collects tests, but only in the ``tests/`` directory, finding
  ``test_cerulean.py``
- ``pytest`` runs ``test_cerulean.py``
- ``test_cerulean.py`` runs ``docker-compose pull`` to pull the images for the
  target containers from DockerHub. ``tests/docker-compose.yml`` guides this
  process.
- ``test_cerulean.py`` runs ``docker-compose build`` to build the test image.
  This copies the whole Cerulean repository into a ``cerulean-test-container``
  image.
- ``test_cerulean.py`` runs ``docker-compose up``, which launches all the target
  containers and then the ``cerulean-test-container``.
- ``pytest`` runs inside ``cerulean-test-container``, collecting tests from
  ``cerulean/test``, and running everything except scheduler tests. This uses
  ``pytest.ini`` from ``tests/container-test``.
- ``pytest`` runs a second time, still inside the container, and runs the
  scheduler tests, with some settings to make it merely take a long time, rather
  than a crazy long time. Coverage is aggregated together with the coverage
  reports from the first ``pytest`` run, and output as ``coverage.xml``.
- ``test_cerulean.py`` copies ``coverage.xml`` out of the container and into the
  main directory.

If we're running on Travis CI, the ``.travis.yml`` file in the root of the
repository contains some additional instructions for setting up the environment,
and a final command to upload the extracted ``coverage.xml`` to Codacy.

The images for the target containers are built from a `separate repository
<https://github.com/MD-Studio/cerulean-test-docker-images>`_.
========
Tutorial
========

Welcome to the Cerulean tutorial. This tutorial demonstrates the basics of using
Cerulean: using local and remote file systems, running processes locally and
remotely, and using schedulers.

To install Cerulean, use

.. code-block:: bash

  pip install cerulean

If you're using Cerulean in a program, you will probably want to use a
virtualenv and install Cerulean into that, together with your other
dependencies.


Accessing files
===============

The file access functions of Cerulean use a ``pathlib``-like API, but unlike in
``pathlib``, Cerulean supports remote file systems. That means that there is no
longer just the local file system, but multiple file systems, and that Path
objects have a particular file system that they are on.

Of course, Cerulean also supports the local file system. To make an object
representing the local file system, you use this:

.. code-block:: python

  import cerulean

  fs = cerulean.LocalFileSystem()

And then you can make a path on the file system using:

.. code-block:: python

  import cerulean

  fs = cerulean.LocalFileSystem()
  my_home_dir = fs / 'home' / 'username'

In this example, ``my_home_dir`` will be a :class:`cerulean.Path` object,
which is very similar to a normal Python ``pathlib.PosixPath``. For example, you
can read the contents of a file through it:

.. code-block:: python

  import cerulean

  fs = cerulean.LocalFileSystem()
  passwd_file = fs / 'etc' / 'passwd'

  users = passwd_file.read_text()
  print(users)

Note that :class:`cerulean.Path` does not support ``open()``. Cerulean can copy
files and stream data from and to them, but it does not offer random access, as
not all remote file access protocols support this.

You can use the ``/`` operator to build paths from components as with
``pathlib``, and there's a wide variety of supported operations. See the API
documentation for :class:`cerulean.Path` for details.

Remote filesystems
------------------

Cerulean supports remote file systems through the SFTP protocol. (It uses the
Paramiko library internally for this.) Accessing a remote file system through
SFTP goes like this:

.. code-block:: python

  import cerulean

  credential = cerulean.PasswordCredential('username', 'password')
  with cerulean.SshTerminal('remotehost.example.com', 22, credential) as term
      with SftpFileSystem(term) as fs:
          my_home_dir = fs / 'home' / 'username'
          test_txt = (my_home_dir / 'test.txt').read_text()
          print(test_txt)

Since we are going to connect to a remote system, we need a credential.
Cerulean has two types of credentials, :class:`PasswordCredential` and
:class:`PubKeyCredential`. They are what you expect, one holds a username and
a password, the other a username, a local path to a public key file, and
optionally a passphrase for the key.

Once we have a credential, we can open a terminal. Like a terminal window on
your desktop, a :class:`Terminal` object lets you run commands. Cerulean
supports local terminals and remote terminals through SSH. Since the SFTP
protocol is an extension to the SSH protocol, we need an SSH terminal connection
first, so we make one, connecting to a host, on a port, with our credential.
This terminal holds an SSH connection, which needs to be closed when we are done
with it. :class:`SshTerminal` is therefore a context manager and needs to be
used in a ``with`` statement. Note that :class:`LocalTerminal` is not a context
manager, as it does not hold any resources.

Once we have the terminal, we can make an :class:`SftpFileSystem` object, and
from there it works just like a local file system. Just like
:class:`SshTerminal`, :class:`SftpFileSystem` is a context manager, so we need
another ``with``-statement.

Copying files
-------------

When running jobs on HPC machines, you often start with copying the input files
from the local system to the HPC machine, and finish with copying the results
back. Cerulean's :meth:`copy` function takes care of this for you, and works as
you would expect:

.. code-block:: python

  import cerulean


  local_fs = cerulean.LocalFileSystem()

  credential = cerulean.PasswordCredential('username', 'password')
  with cerulean.SshTerminal('remotehost.example.com', 22, credential) as term
      with SftpFileSystem(term) as remote_fs:
          input_file = local_fs / 'home' / 'username' / 'input.txt'
          job_dir = remote_fs / 'home' / 'username' / 'my_job'
          cerulean.copy(input_file, job_dir)

          # run job and wait for it to finish

          output_file = local_fs / 'home' / 'username' / 'output.txt'
          cerulean.copy(job_dir / 'output.txt', output_file)

Running commands
================

If you have read the above, then the secret is already out: running commands
using Cerulean is done using a :class:`Terminal`. For example, you can run a
command locally using:

.. code-block:: python

  import cerulean

  term = cerulean.LocalTerminal()

  exit_code, stdout_text, stderr_text = term.run(
          10.0, 'ls', ['-l'], None, '/home/username')

The first argument to :meth:`Terminal.run` is a timeout value in seconds,
which determines how long Cerulean will wait for the command to finish. The
second argument is the command to run, followed by a list of arguments. Next is
an optional string that, if you specify it, will be fed into the standard input
of the program you are starting. The final argument is a string specifying the
working directory in which to execute the command.

The function returns a tuple containing three values: the exit code of the
process (or `None` if it didn't finish in time), a string containing text
printed to standard output, and a string containing text printed to standard
error by the command you ran.

Running commands remotely through SSH of course works in exactly the same way,
except you use an :class:`SshTerminal`, as above:

.. code-block:: python

  import cerulean

  credential = cerulean.PasswordCredential('username', 'password')
  with cerulean.SshTerminal('remotehost.example.com', 22, credential) as term
      exit_code, stdout_text, stderr_text = term.run(
              10.0, 'ls', ['-l'], None, '/home/username')


Submitting jobs
===============

On High Performance Computing machines, you don't run commands directly.
Instead, you submit batch jobs to a scheduler, which will place them in a queue,
and run them when everyone else in line before you is done. The most popular
scheduler at the moment seems to be Slurm, but Cerulean also supports
Torque/PBS.

The usual way of working with a scheduler is to use ``ssh`` to connect to the
cluster, where you run commands that submit jobs and check on their status.
Cerulean works in the same way:

.. code-block:: python

  import cerulean
  import time

  credential = cerulean.PasswordCredential('username', 'password')
  with cerulean.SshTerminal('remotehost.example.com', 22, credential) as term
      sched = cerulean.SlurmScheduler(term)

      job = cerulean.JobDescription()
      job.name = 'cerulean_test'
      job.command = 'ls'
      job.arguments = ['-l']

      job_id = sched.submit_job(job)

      time.sleep(5)
      status = sched.get_status(job_id)

      if status == cerulean.JobStatus.DONE:
          exit_code = sched.get_exit_code()
          print('Job exited with code {}'.format(exit_code))

Of course, if you intend to run your submission script on the head node, then
the scheduler is local, and you want to use a :class:`LocalTerminal` with your
:class:`SlurmScheduler`. If your HPC machine runs Torque/PBS, use a
:class:`TorqueScheduler` instead.


More information
================

To find all the details of what Cerulean can do and how to do it, please refer
to the :doc:`API documentation<apidocs/cerulean>`.
.. Cerulean documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Cerulean's documentation!
==========================================================

Cerulean is a Python 3 library for connecting to HPC compute resources, such as
compute clusters and supercomputers. It lets you copy files between local and
SFTP filesystems using a ``pathlib``-like API, it lets you start processes
locally and remotely via SSH, and it lets you submit jobs to schedulers such as
Slurm and Torque/PBS.

Cerulean supports Python 3.4 and later.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial


API Reference
=============

.. toctree::
  :maxdepth: 1

   API reference <apidocs/cerulean.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
