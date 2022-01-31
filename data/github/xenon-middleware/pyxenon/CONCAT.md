Python interface to Xenon 3.0
=============================
|ZenodoBadge| |ReadTheDocsBadge| |Apache2License| |BuildStatus| |CodacyBadge|

Python interface to the `Xenon middleware library, v. 3.0
<http://xenon-middleware.github.io/xenon/>`__. Xenon provides a simple programming
interface to various pieces of software that can be used to access distributed
compute and storage resources.

Underneath it uses `GRPC <https://grpc.io>`__, to connect to the `Xenon-GRPC
<https://github.com/xenon-middleware/xenon-grpc>`__ service.
We've taken care to mirror the original Java API in this Python module as much
as possible.

Installing
----------
Clone this repository, and do::

    pip install .

The code will appear on PyPI when it is ready for release.

Documentation
-------------
The compiled documentation is hosted on `Read the Docs
<http://pyxenon.readthedocs.io/en/latest>`__. This includes a quick-start
guide.

Example
-------

.. code-block:: python

    import xenon
    from pathlib import Path
    import os

    xenon.init()

    # create a new job scheduler, using SSH to localhost to submit new jobs.
    with xenon.Scheduler.create(
            adaptor='ssh', location='localhost') as scheduler:

        # make a new job description. The executable must already be present on
        # the target host.
        target = Path('.') / 'stdout.txt'
        desc = xenon.JobDescription(
            executable='hostname',
            stdout=str(target.resolve()))

        # submit a job
        job = scheduler.submit_batch_job(desc)
        status = scheduler.wait_until_done(job, 1000)

        # read the standard output of the job. We can do this directly because
        # we ran on localhost, otherwise, we need to transfer the file first.
        with open(target) as f:
            print(f.read())

Development
-----------
PyXenon ships with the `Xenon-GRPC` jar-file and command-line executable. If
these need upgrading, build them manually, following instructions at
`Xenon-GRPC <https://github.com/nlesc/xenon-grpc>`__, and place the contents of the
``build/install/xenon-grpc-shadow`` folder (``lib`` and ``bin``) here.

To generate the `GRPC` code, run ``scripts/protoc.sh`` from the project root.

Steps for upgrading Xenon-GRPC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``xenon-grpc`` repository, run::

   ./gradlew shadowJar

The updated JAR file will be located in ``./build/libs/xenon-grpc-${version}.jar``.
Also make sure to copy updated binary files from ``./build/install/xenon-grpc/bin``.
The target files should go to the ``./bin`` and ``./lib`` folders in the ``pyxenon`` repository.

Update ``./xenon/versions.py``.

To update the GRPC Python bindings, run the ``./scripts/protoc.sh`` script.

To update the Xenon adaptor documentation, first ``pip install --upgrade .``, then run ``python ./scripts/print_adaptor_docs.py > docs/adaptors.rst``.

Run ``tox`` before pushing anything back to github.

Testing
-------
Unit tests all run against the `local` scheduler and the `file` adaptor for
filesystems. To run them, just do::

    $ pytest ./tests

For faster testing it may be useful to start the ``xenon-grpc`` daemon
manually; start it in a separate terminal as it may give useful output for
debugging.

For integration testing, run the following docker container to test against
remote slurm

.. code-block:: bash

    docker run --detach --publish 10022:22 nlesc/xenon-slurm:17

An example of some code running against this container is in
``examples/tutorial.py``.

Contributing
------------

Contributions can be made using GitHub pull requests. To add a feature,
first install the test requirements

::

    pip install -U tox

and then run

::

    tox

until all tests succeed. The command checks against flake8 code
standards and syntax errors on Python 3.5 and 3.6. Then commit, to make sure
the change didn't break any code. The pull request will be evaluated in
`Travis <https://travis-ci.org/NLeSC/pyxenon>`__.

.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.60929.svg
   :target: http://dx.doi.org/10.5281/zenodo.60929
.. |PyPi version| image:: https://img.shields.io/pypi/v/pyxenon.svg
   :target: https://pypi.python.org/pypi/pyxenon
.. |Apache2License| image:: https://img.shields.io/github/license/NLeSC/pyxenon.svg?branch=master
   :target: https://raw.githubusercontent.com/NLeSC/pyxenon/master/LICENSE
.. |PythonVersions| image:: https://img.shields.io/pypi/pyversions/pyxenon.svg
.. |BuildStatus| image:: https://travis-ci.org/xenon-middleware/pyxenon.svg?branch=master
   :target: https://travis-ci.org/NLeSC/pyxenon
.. |CodacyBadge| image:: https://api.codacy.com/project/badge/grade/35e155e3bb08459aa2c24622d5fdb0d3
   :target: https://www.codacy.com/app/NLeSC/pyxenon
.. |ReadTheDocsBadge| image:: https://readthedocs.org/projects/pyxenon/badge/?version=latest
   :target: http://pyxenon.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |ZenodoBadge| image:: https://zenodo.org/badge/47132292.svg
   :target: https://zenodo.org/badge/latestdoi/47132292
Adaptors
========
This section contains the adaptor documentation which is generated from the
information provided by the adaptors themselves.

.. contents::


File System
-----------

.. note:: Supported property names should be prefixed with
``"xenon.adaptors.filesystems"``.  We've left this prefix out to improve
readability of the tables.


S3
~~
The S3 adaptor uses Apache JClouds to talk to s3 and others. To
authenticate use PasswordCredential with access key id as username and
secret access key as password

+------------------------------------+----------------------+
| field                              | value                |
+====================================+======================+
| supports_third_party_copy          | False                |
+------------------------------------+----------------------+
| can_create_symboliclinks           | False                |
+------------------------------------+----------------------+
| can_read_symboliclinks             | False                |
+------------------------------------+----------------------+
| is_connectionless                  | True                 |
+------------------------------------+----------------------+
| supported_credentials              | `PasswordCredential` |
+------------------------------------+----------------------+
| can_append                         | False                |
+------------------------------------+----------------------+
| supports_reading_posix_permissions | False                |
+------------------------------------+----------------------+
| supports_setting_posix_permissions | False                |
+------------------------------------+----------------------+
| supports_rename                    | False                |
+------------------------------------+----------------------+
| needs_size_beforehand              | True                 |
+------------------------------------+----------------------+

location string:
    * `http[s]://host[:port]/bucketname[/workdir]`
    * `https://s3.region.amazonaws.com/bucketname[/workdir]`

supported properties:

+---------------+-------------------------------------------------------+-----------+---------+
| name          | description                                           | data_type | default |
+===============+=======================================================+===========+=========+
| s3.bufferSize | The buffer size to use when copying files (in bytes). | size      | `64K`   |
+---------------+-------------------------------------------------------+-----------+---------+

File
~~~~
This is the local file adaptor that implements file functionality for
local access.

+------------------------------------+---------------------+
| field                              | value               |
+====================================+=====================+
| supports_third_party_copy          | False               |
+------------------------------------+---------------------+
| can_create_symboliclinks           | True                |
+------------------------------------+---------------------+
| can_read_symboliclinks             | True                |
+------------------------------------+---------------------+
| is_connectionless                  | True                |
+------------------------------------+---------------------+
| supported_credentials              | `DefaultCredential` |
+------------------------------------+---------------------+
| can_append                         | True                |
+------------------------------------+---------------------+
| supports_reading_posix_permissions | True                |
+------------------------------------+---------------------+
| supports_setting_posix_permissions | True                |
+------------------------------------+---------------------+
| supports_rename                    | True                |
+------------------------------------+---------------------+
| needs_size_beforehand              | False               |
+------------------------------------+---------------------+

location string:
    * `(null)`
    * `(empty string)`
    * `[/workdir]`
    * `driveletter:[/workdir]`

supported properties:

+-----------------+-------------------------------------------------------+-----------+---------+
| name            | description                                           | data_type | default |
+=================+=======================================================+===========+=========+
| file.bufferSize | The buffer size to use when copying files (in bytes). | size      | `64K`   |
+-----------------+-------------------------------------------------------+-----------+---------+

Sftp
~~~~
The SFTP adaptor implements all file access functionality to remote
SFTP servers

+------------------------------------+----------------------------------------------------+
| field                              | value                                              |
+====================================+====================================================+
| supports_third_party_copy          | False                                              |
+------------------------------------+----------------------------------------------------+
| can_create_symboliclinks           | True                                               |
+------------------------------------+----------------------------------------------------+
| can_read_symboliclinks             | True                                               |
+------------------------------------+----------------------------------------------------+
| is_connectionless                  | False                                              |
+------------------------------------+----------------------------------------------------+
| supported_credentials              | `DefaultCredential`, `CertificateCredential`,      |
|                                    | `PasswordCredential`, `CredentialMap`              |
+------------------------------------+----------------------------------------------------+
| can_append                         | True                                               |
+------------------------------------+----------------------------------------------------+
| supports_reading_posix_permissions | True                                               |
+------------------------------------+----------------------------------------------------+
| supports_setting_posix_permissions | True                                               |
+------------------------------------+----------------------------------------------------+
| supports_rename                    | True                                               |
+------------------------------------+----------------------------------------------------+
| needs_size_beforehand              | False                                              |
+------------------------------------+----------------------------------------------------+

location string:
    * `host[:port][/workdir]`

supported properties:

+----------------------------+------------------------------------------------------------+-----------+---------+
| name                       | description                                                | data_type | default |
+============================+============================================================+===========+=========+
| sftp.strictHostKeyChecking | Enable strict host key checking.                           | boolean   | `true`  |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.loadKnownHosts        | Load the standard known_hosts file.                        | boolean   | `true`  |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.loadSshConfig         | Load the OpenSSH config file.                              | boolean   | `true`  |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.agent                 | Use a (local) ssh-agent.                                   | boolean   | `false` |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.agentForwarding       | Use ssh-agent forwarding when setting up a connection.     | boolean   | `false` |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.connection.timeout    | The timeout for creating and authenticating connections    | natural   | `10000` |
|                            | (in milliseconds).                                         |           |         |
+----------------------------+------------------------------------------------------------+-----------+---------+
| sftp.bufferSize            | The buffer size to use when copying files (in bytes).      | size      | `64K`   |
+----------------------------+------------------------------------------------------------+-----------+---------+

Ftp
~~~
The FTP adaptor implements file access on remote ftp servers.

+------------------------------------+-------------------------------------------+
| field                              | value                                     |
+====================================+===========================================+
| supports_third_party_copy          | False                                     |
+------------------------------------+-------------------------------------------+
| can_create_symboliclinks           | False                                     |
+------------------------------------+-------------------------------------------+
| can_read_symboliclinks             | True                                      |
+------------------------------------+-------------------------------------------+
| is_connectionless                  | False                                     |
+------------------------------------+-------------------------------------------+
| supported_credentials              | `DefaultCredential`, `PasswordCredential` |
+------------------------------------+-------------------------------------------+
| can_append                         | True                                      |
+------------------------------------+-------------------------------------------+
| supports_reading_posix_permissions | True                                      |
+------------------------------------+-------------------------------------------+
| supports_setting_posix_permissions | False                                     |
+------------------------------------+-------------------------------------------+
| supports_rename                    | True                                      |
+------------------------------------+-------------------------------------------+
| needs_size_beforehand              | False                                     |
+------------------------------------+-------------------------------------------+

location string:
    * `host[:port][/workdir]`

supported properties:

+----------------+-------------------------------------------------------+-----------+---------+
| name           | description                                           | data_type | default |
+================+=======================================================+===========+=========+
| ftp.bufferSize | The buffer size to use when copying files (in bytes). | size      | `64K`   |
+----------------+-------------------------------------------------------+-----------+---------+

Webdav
~~~~~~
The webdav file adaptor implements file access to remote webdav
servers.

+------------------------------------+-------------------------------------------+
| field                              | value                                     |
+====================================+===========================================+
| supports_third_party_copy          | False                                     |
+------------------------------------+-------------------------------------------+
| can_create_symboliclinks           | False                                     |
+------------------------------------+-------------------------------------------+
| can_read_symboliclinks             | False                                     |
+------------------------------------+-------------------------------------------+
| is_connectionless                  | True                                      |
+------------------------------------+-------------------------------------------+
| supported_credentials              | `DefaultCredential`, `PasswordCredential` |
+------------------------------------+-------------------------------------------+
| can_append                         | False                                     |
+------------------------------------+-------------------------------------------+
| supports_reading_posix_permissions | False                                     |
+------------------------------------+-------------------------------------------+
| supports_setting_posix_permissions | False                                     |
+------------------------------------+-------------------------------------------+
| supports_rename                    | True                                      |
+------------------------------------+-------------------------------------------+
| needs_size_beforehand              | False                                     |
+------------------------------------+-------------------------------------------+

location string:
    * `http://host[:port][/workdir]`
    * `https://host[:port][/workdir]`

supported properties:

+-------------------+-------------------------------------------------------+-----------+---------+
| name              | description                                           | data_type | default |
+===================+=======================================================+===========+=========+
| webdav.bufferSize | The buffer size to use when copying files (in bytes). | size      | `64K`   |
+-------------------+-------------------------------------------------------+-----------+---------+


Scheduler
---------

.. note:: Supported property names should be prefixed with
``"xenon.adaptors.schedulers"``.  We've left this prefix out to improve
readability of the tables.


Local
~~~~~
The local jobs adaptor implements all functionality by emulating a
local queue.

+-----------------------+---------------------+
| field                 | value               |
+=======================+=====================+
| is_embedded           | True                |
+-----------------------+---------------------+
| supports_interactive  | True                |
+-----------------------+---------------------+
| supports_batch        | True                |
+-----------------------+---------------------+
| uses_file_system      | True                |
+-----------------------+---------------------+
| supported_credentials | `DefaultCredential` |
+-----------------------+---------------------+

location string:
    * `[/workdir]`

supported properties:

+-------------------------------------+--------------------------------------------+-----------+---------+
| name                                | description                                | data_type | default |
+=====================================+============================================+===========+=========+
| local.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq.                                |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+

Ssh
~~~
The SSH job adaptor implements all functionality to start jobs on ssh
servers.

+-----------------------+---------------------------------------------------------------------------------+
| field                 | value                                                                           |
+=======================+=================================================================================+
| is_embedded           | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supports_interactive  | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supports_batch        | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| uses_file_system      | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supported_credentials | `DefaultCredential`, `CertificateCredential`, `PasswordCredential`,             |
|                       | `CredentialMap`                                                                 |
+-----------------------+---------------------------------------------------------------------------------+

location string:
    * `host[:port][/workdir][ via:otherhost[:port]]*`

supported properties:

+-----------------------------------+--------------------------------------------+-----------+---------+
| name                              | description                                | data_type | default |
+===================================+============================================+===========+=========+
| ssh.strictHostKeyChecking         | Enable strict host key checking.           | boolean   | `true`  |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadKnownHosts                | Load the standard known_hosts file.        | boolean   | `true`  |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadSshConfig                 | Load the OpenSSH config file.              | boolean   | `true`  |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.agent                         | Use a (local) ssh-agent.                   | boolean   | `false` |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.agentForwarding               | Use ssh-agent forwarding                   | boolean   | `false` |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.timeout                       | The timeout for the connection setup and   | long      | `10000` |
|                                   | authetication (in milliseconds).           |           |         |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                   | jobs (in milliseconds).                    |           |         |
+-----------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                   | the multiq..                               |           |         |
+-----------------------------------+--------------------------------------------+-----------+---------+

At
~~
The At Adaptor submits jobs to an at scheduler.  This adaptor uses
either the local or the ssh scheduler adaptor to run commands on the
machine running at,  and the file or the stfp filesystem adaptor to
gain access to the filesystem of that machine.

+-----------------------+---------------------------------------------------------------------------------+
| field                 | value                                                                           |
+=======================+=================================================================================+
| is_embedded           | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_interactive  | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_batch        | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| uses_file_system      | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supported_credentials | `DefaultCredential`, `CertificateCredential`, `PasswordCredential`,             |
|                       | `CredentialMap`                                                                 |
+-----------------------+---------------------------------------------------------------------------------+

location string:
    * `local://[/workdir]`
    * `ssh://host[:port][/workdir][ via:otherhost[:port]]*`

supported properties:

+-------------------------------------+--------------------------------------------+-----------+---------+
| name                                | description                                | data_type | default |
+=====================================+============================================+===========+=========+
| at.poll.delay                       | Number of milliseconds between polling the | long      | `1000`  |
|                                     | status of a job.                           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.strictHostKeyChecking           | Enable strict host key checking.           | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadKnownHosts                  | Load the standard known_hosts file.        | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadSshConfig                   | Load the OpenSSH config file.              | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agent                           | Use a (local) ssh-agent.                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agentForwarding                 | Use ssh-agent forwarding                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.timeout                         | The timeout for the connection setup and   | long      | `10000` |
|                                     | authetication (in milliseconds).           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.pollingDelay              | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.multi.maxConcurrentJobs   | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq..                               |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq.                                |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+

Slurm
~~~~~
The Slurm Adaptor submits jobs to a Slurm scheduler.  This adaptor
uses either the local or the ssh scheduler adaptor to run commands on
the machine running Slurm,  and the file or the stfp filesystem
adaptor to gain access to the filesystem of that machine.

+-----------------------+---------------------------------------------------------------------------------+
| field                 | value                                                                           |
+=======================+=================================================================================+
| is_embedded           | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_interactive  | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supports_batch        | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| uses_file_system      | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supported_credentials | `DefaultCredential`, `CertificateCredential`, `PasswordCredential`,             |
|                       | `CredentialMap`                                                                 |
+-----------------------+---------------------------------------------------------------------------------+

location string:
    * `local://[/workdir]`
    * `ssh://host[:port][/workdir][ via:otherhost[:port]]*`

supported properties:

+-------------------------------------+--------------------------------------------+-----------+---------+
| name                                | description                                | data_type | default |
+=====================================+============================================+===========+=========+
| slurm.disable.accounting.usage      | Do not use accounting info of slurm, even  | boolean   | `false` |
|                                     | when available. Mostly for testing         |           |         |
|                                     | purposes                                   |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| slurm.poll.delay                    | Number of milliseconds between polling the | long      | `1000`  |
|                                     | status of a job.                           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.strictHostKeyChecking           | Enable strict host key checking.           | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadKnownHosts                  | Load the standard known_hosts file.        | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadSshConfig                   | Load the OpenSSH config file.              | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agent                           | Use a (local) ssh-agent.                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agentForwarding                 | Use ssh-agent forwarding                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.timeout                         | The timeout for the connection setup and   | long      | `10000` |
|                                     | authetication (in milliseconds).           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.pollingDelay              | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.multi.maxConcurrentJobs   | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq..                               |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq.                                |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+

Gridengine
~~~~~~~~~~
The SGE Adaptor submits jobs to a (Sun/Oracle/Univa) Grid Engine
scheduler. This adaptor uses either the local or the ssh scheduler
adaptor to run commands on the machine running Grid Engine,  and the
file or the stfp filesystem adaptor to gain access to the filesystem
of that machine.

+-----------------------+---------------------------------------------------------------------------------+
| field                 | value                                                                           |
+=======================+=================================================================================+
| is_embedded           | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_interactive  | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_batch        | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| uses_file_system      | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supported_credentials | `DefaultCredential`, `CertificateCredential`, `PasswordCredential`,             |
|                       | `CredentialMap`                                                                 |
+-----------------------+---------------------------------------------------------------------------------+

location string:
    * `local://[/workdir]`
    * `ssh://host[:port][/workdir][ via:otherhost[:port]]*`

supported properties:

+-------------------------------------+--------------------------------------------+-----------+---------+
| name                                | description                                | data_type | default |
+=====================================+============================================+===========+=========+
| gridengine.ignore.version           | Skip version check is skipped when         | boolean   | `false` |
|                                     | connecting to remote machines. WARNING: it |           |         |
|                                     | is not recommended to use this setting in  |           |         |
|                                     | production environments!                   |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| gridengine.accounting.grace.time    | Number of milliseconds a job is allowed to | long      | `60000` |
|                                     | take going from the queue to the qacct     |           |         |
|                                     | output.                                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| gridengine.poll.delay               | Number of milliseconds between polling the | long      | `1000`  |
|                                     | status of a job.                           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.strictHostKeyChecking           | Enable strict host key checking.           | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadKnownHosts                  | Load the standard known_hosts file.        | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadSshConfig                   | Load the OpenSSH config file.              | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agent                           | Use a (local) ssh-agent.                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agentForwarding                 | Use ssh-agent forwarding                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.timeout                         | The timeout for the connection setup and   | long      | `10000` |
|                                     | authetication (in milliseconds).           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.pollingDelay              | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.multi.maxConcurrentJobs   | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq..                               |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq.                                |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+

Torque
~~~~~~
The Torque Adaptor submits jobs to a TORQUE batch system. This adaptor
uses either the local or the ssh scheduler adaptor to run commands on
the machine running TORQUE,  and the file or the stfp filesystem
adaptor to gain access to the filesystem of that machine.

+-----------------------+---------------------------------------------------------------------------------+
| field                 | value                                                                           |
+=======================+=================================================================================+
| is_embedded           | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_interactive  | False                                                                           |
+-----------------------+---------------------------------------------------------------------------------+
| supports_batch        | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| uses_file_system      | True                                                                            |
+-----------------------+---------------------------------------------------------------------------------+
| supported_credentials | `DefaultCredential`, `CertificateCredential`, `PasswordCredential`,             |
|                       | `CredentialMap`                                                                 |
+-----------------------+---------------------------------------------------------------------------------+

location string:
    * `local://[/workdir]`
    * `ssh://host[:port][/workdir][ via:otherhost[:port]]*`

supported properties:

+-------------------------------------+--------------------------------------------+-----------+---------+
| name                                | description                                | data_type | default |
+=====================================+============================================+===========+=========+
| torque.ignore.version               | Skip version check is skipped when         | boolean   | `false` |
|                                     | connecting to remote machines. WARNING: it |           |         |
|                                     | is not recommended to use this setting in  |           |         |
|                                     | production environments!                   |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| torque.accounting.grace.time        | Number of milliseconds a job is allowed to | long      | `60000` |
|                                     | take going from the queue to the accinfo   |           |         |
|                                     | output.                                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| torque.poll.delay                   | Number of milliseconds between polling the | long      | `1000`  |
|                                     | status of a job.                           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.strictHostKeyChecking           | Enable strict host key checking.           | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadKnownHosts                  | Load the standard known_hosts file.        | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.loadSshConfig                   | Load the OpenSSH config file.              | boolean   | `true`  |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agent                           | Use a (local) ssh-agent.                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.agentForwarding                 | Use ssh-agent forwarding                   | boolean   | `false` |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.timeout                         | The timeout for the connection setup and   | long      | `10000` |
|                                     | authetication (in milliseconds).           |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.pollingDelay              | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| ssh.queue.multi.maxConcurrentJobs   | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq..                               |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.pollingDelay            | The polling delay for monitoring running   | long      | `1000`  |
|                                     | jobs (in milliseconds).                    |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+
| local.queue.multi.maxConcurrentJobs | The maximum number of concurrent jobs in   | integer   | `4`     |
|                                     | the multiq.                                |           |         |
+-------------------------------------+--------------------------------------------+-----------+---------+

Advanced: Streaming & Interactive jobs
======================================

In several cases it is desireable to stream data from/to interactive jobs as
well as data to a remote filesystem. The GRPC API has build-in support for
asynchronous streaming through many simultaneous requests. In Python this API is
exposed in terms of generators.

Example: an online job
----------------------

In this example we'll show how to obtain
bi-directional communication with an online job. An online job is started with
:py:meth:`Scheduler.submit_online_job()`.

Streaming input, a.k.a. The Halting Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to stream input to the online job. In the :doc:`quick-start`, we saw
that we could send data to a stream by simply giving a list of bytes objects.
Here we aim a bit more advanced to play a kind of real-time ping-pong with a
remote process. We need to provide `PyXenon` with an generator that pulls its
messages from a queue. The GRPC module ensures that this generator is being run
asynchonously from the main thread.

The tricky part is that we need to be able to tell the generator when the work
is done and no more input is to be expected. We could have it recieve strings and
make it check for end-of-file messages in some way, but in essence we'll always
have to define a little protocol to deal with the finiteness of the generator's
life. To make this explicit we define a little 2-tuple micro-language:

+-------------------------------+--------------------------+
| message                       | action                   |
+===============================+==========================+
| ``('msg', <value: string>)``  | ``yield value.encode()`` |
+-------------------------------+--------------------------+
| ``('end', None)``             | ``return``               |
+-------------------------------+--------------------------+

Implementing this:

.. code-block:: python

    from queue import Queue

    def make_input_stream():
        input_queue = Queue()

        def input_stream():
            while True:
                cmd, value = input_queue.get()
                if cmd == 'end':
                    input_queue.task_done()
                    return
                elif cmd == 'msg':
                    yield value.encode()
                    input_queue.task_done()

        return input_queue, input_stream

Reading output
~~~~~~~~~~~~~~
The return-value of :py:meth:`submit_online_job()` is an iterator yielding
objects of type `SubmitOnlineJobResponse`. These objects have a ``stdout``
field containing (binary) data that the job wrote to standard output, as well
as a ``stderr`` field containing data written to standard error. For any message
either field may be empty or not. In this example we're only interested in data
from ``stdout``:

.. code-block:: python

    def get_stdout(stream):
        return stream.next().stdout.decode()

The "remote" script
~~~~~~~~~~~~~~~~~~~
For the purpose of this example, we have defined a small Python ``rot13``
program:

.. code-block:: python
    :caption: rot13.py

    import codecs

    try:
        while True:
            line = input()
            print(codecs.encode(line, 'rot_13'))

    except EOFError:
        pass

Defining the job
~~~~~~~~~~~~~~~~
Online job descriptions are the same as normal job descriptions.

.. code-block:: python

    # our input lines
    input_lines = [
        "Zlfgvp aboyr tnf,",
        "Urnil lrg syrrgvat sebz tenfc,",
        "Oyhr yvxr oheavat vpr."
    ]

    # the job description, make sure you run the script from the examples
    # directory!
    job_description = xenon.JobDescription(
        executable='python',
        arguments=['rot13.py'],
        queue_name='multi')

Putting it together
~~~~~~~~~~~~~~~~~~~

The rest is history.

.. code-block:: python

    import xenon

    # start the xenon-grpc server
    xenon.init()

    # on the local adaptor
    with xenon.Scheduler.create(adaptor='local') as scheduler:
        input_queue, input_stream = make_input_stream()

        # submit an interactive job, this gets us the job-id and a stream
        # yielding job output from stdout and stderr.
        job, output_stream = scheduler.submit_interactive_job(
            description=job_description, stdin_stream=input_stream())

        # next we feed the input_queue with messages
        try:
            for line in input_lines:
                print(" [sending]   " + line)
                input_queue.put(('msg', line + '\n'))
                msg = get_stdout(output_stream)
                print("[received]   " + msg)

        # make sure to close our end whatever may happen
        finally:
            input_queue.put(('end', None))
            input_queue.join()

        scheduler.wait_until_done(job)


Protocol definitions
--------------------
It can be instructive to see what the GRPC protocol with respect to interactive
jobs looks like.

.. code-block:: proto

    message SubmitInteractiveJobRequest {
        Scheduler scheduler = 1;
        JobDescription description = 2;
        bytes stdin = 3;
    }

    message SubmitInteractiveJobResponse {
        Job job = 1;
        bytes stdout = 2;
        bytes stderr = 3;
    }

    service SchedulerService {
        rpc submitInteractiveJob(
                stream SubmitInteractiveJobRequest)
            returns (stream SubmitInteractiveJobResponse) {}
    }

In `PyXenon` the remote procedure call ``submitInteractiveJob`` is wrapped to
the method :py:meth:`submit_interactive_job()` of the :py:class:`Scheduler`
class. Note that the ``SubmitInteractiveJobRequest`` specifies (next to the
scheduler, which is obtained from ``self`` in the method call) the job
description and ``bytes`` for standard input. Requests of this type are
streamed.  This means that GRPC expects to get an iterator of
``SubmitInteractiveJobRequest`` objets.

The `PyXenon` :py:meth:`submit_interactive_job()` method separates the
job-description and input-stream arguments. Sending the ``scheduler`` and
``description`` fields in the first request, followed up by a sequence of
requests where only the ``stdin`` field is specified. This latter sequence
is yielded from the ``stdin_stream`` argument.

Similarly, the first item in the output stream is guaranteed to only contain
the job-id, this first item is available immediately. Subsequent calls to
``next(output_stream)`` will block until output is available. The
:py:meth:`submit_interactive_job()` method takes the first item of the
iterator, and extracts the job-id. The user recieves a tuple with the
extracted job-id and the iterator.
.. PyXenon documentation master file, created by
   sphinx-quickstart on Wed Aug 30 10:12:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyXenon's documentation!
===================================

.. toctree::
   :maxdepth: 2

   self
   quick-start
   streaming
   adaptors
   api

The PyXenon module interfaces with Xenon-GRPC to get an interface to the Xenon
2.0 Java library. We kept this interface close to the original Java API.
PyXenon 2.0 only works on Python 3.

Installing
----------

.. code-block:: bash

    pip install pyxenon


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Quick Start
-----------

We like to test Xenon against a `Docker`_ image: `nlesc/xenon-slurm`_.
If you have docker all setup, you can run this image as follows::

    $ docker pull nlesc/xenon-slurm
    ...
    $ docker run --detach --publish 10022:22 nlesc/xenon-slurm

Try logging onto this image by `ssh`, to make sure everything works. The
username is `xenon`, the password `javagat`::

    $ ssh localhost -p 10022 -l xenon
    xenon@localhost's password: <javagat>
    $ exit
    Connection to localhost closed.

Starting the server
~~~~~~~~~~~~~~~~~~~
To get anything done in PyXenon, we need to start the GRPC server:

.. code-block:: python

    import xenon

    xenon.init()

Writing to a remote filesystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Next, let's try to copy a file to the container. We need credentials to
access anything on the remote side.

.. code-block:: python

    from xenon import PasswordCredential, FileSystem

    credential = PasswordCredential(
        username='xenon',
        password='javagat')

    remotefs = FileSystem.create(
        'sftp', location='localhost:10022',
        password_credential=credential)

We can write to a file by streaming. The second argument to
:py:meth:`write_to_file` should be an iterable. It will be read in a separate
thread, so it is allowed to be blocking. Here we'll do nothing so fancy:

.. code-block:: python

    from xenon import Path

    target = Path('hello.sh')

    if remotefs.exists(target):
        remotefs.delete(target)

    remotefs.write_to_file(
        target,
        [b'#!/bin/sh\n',
         b'echo "Hello, World!"\n'])

Running a script
~~~~~~~~~~~~~~~~
The remote machine runs a SLURM job scheduler. We describe a job in a
:py:class:`JobDescription` object. This seems a bit long-winded, but in
practice you'll be reusing the descriptions a lot.

.. code-block:: python

    from xenon import Scheduler, JobDescription

    scheduler = Scheduler.create(
        adaptor='slurm',
        location='ssh://localhost:10022',
        password_credential=credential)

    job_description = JobDescription(
        executable='/bin/sh',
        arguments=['hello.sh'],
        stdout='result.txt')

    job = scheduler.submit_batch_job(job_description)

    state = scheduler.wait_until_done(job)
    print(state)


Retrieving the result
~~~~~~~~~~~~~~~~~~~~~
Just as we can write data by sending an iterable, we can read data from a file
and recieve a generator yielding bytes objects. Here we realize the transfer by
joining the data chunks into a string:

.. code-block:: python

    text = ''.join(chunk.decode() for chunk in
        remotefs.read_from_file(Path('result.txt')))
    print(text)

.. _Docker: https://www.docker.com/
.. _nlesc/xenon-slurm: https://hub.docker.com/r/nlesc/xenon-slurm/
API
===

.. contents::

.. automodule:: xenon

The Server
----------

.. autofunction:: init

File Systems
------------

.. autoclass:: FileSystem
    :members:

.. autoclass:: Path
    :members:

Message classes
~~~~~~~~~~~~~~~
.. autoclass:: PosixFilePermission
    :members:
    :undoc-members:

.. autoclass:: CopyMode
    :members:
    :undoc-members:

.. autoclass:: CopyStatus
    :members:
    :undoc-members:

Schedulers
----------
.. autoclass:: Scheduler
    :members:

Message classes
~~~~~~~~~~~~~~~
.. autoclass:: Job
    :members:

.. autoclass:: JobDescription
    :members:

.. autoclass:: JobStatus
    :members:
    :undoc-members:

.. autoclass:: QueueStatus
    :members:
    :undoc-members:

Credentials
-----------
.. autoclass:: CertificateCredential
    :members:

.. autoclass:: PasswordCredential
    :members:

Exceptions
----------
.. automodule:: xenon.exceptions
    :members:
