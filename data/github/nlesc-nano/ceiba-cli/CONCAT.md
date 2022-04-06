##########
Change Log
##########

1.0.0 [22/03/2021]
******************
New
---
* Make the client generic to uses other that quantum chemistry


0.3.0 [22/02/2021]
******************

Changed
-------
* Make library more generic by removing references to the smiles (#26)

New
---
* Accept jobs as a list of JSON objects (#26)

0.2.0 [Unreleased]
******************

Added
-----
* Command to manage the jobs status
* Allow to pass the input via the command line (#10)
* Allow to retrieve the available collections (#13)
* Authentication functionality (#18)

Changed
-------
* Make input file optional (#10)
* Choose automatically the new collection name when adding jobs (#8)


  Changed
-------
* Running jobs must include username (#9)

0.1.0 [03/11/2020]
******************

Added
-----

* Interface to web service using `add`, `query`, `compute` or `report` actions (#1)
* Interface to [Openstack Swift](https://docs.openstack.org/swift/latest/) to store large objects
* Allow to query jobs based on their size (#4)
* Allow to report either in CSV or JSON format.

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
reported by contacting the project team at f.zapata@esciencecenter.nl. All
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
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/nlesc-nano/ceiba-cli/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/nlesc-nano/ceiba-cli/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest main commit. While working on your feature branch, make sure to stay up to date with the main branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the ceiba-cli repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
.. image:: https://github.com/nlesc-nano/ceiba-cli/workflows/build/badge.svg
   :target: https://github.com/nlesc-nano/ceiba-cli/actions
.. image:: https://readthedocs.org/projects/ceiba-cli/badge/?version=latest
   :target: https://ceiba-cli.readthedocs.io/en/latest/?badge=latest
.. image:: https://codecov.io/gh/nlesc-nano/ceiba-cli/branch/main/graph/badge.svg
  :target: https://codecov.io/gh/nlesc-nano/ceiba-cli
.. image:: https://zenodo.org/badge/296641388.svg
   :target: https://zenodo.org/badge/latestdoi/296641388

#########
ceiba-cli
#########

command line interface to interact with the `insilico web server <https://github.com/nlesc-nano/ceiba>`_.
See the `documentation <https://ceiba-cli.readthedocs.io/en/latest/>`_ and `blog post <https://blog.esciencecenter.nl/building-a-web-service-to-manage-scientific-simulation-data-using-graphql-a0bbf1c3f6e9>`_.


Installation
------------

To install ceiba-cli, do:

.. code-block:: console

python -m pip install git+https://github.com/nlesc-nano/ceiba-cli.git@main

Contributing
############

If you want to contribute to the development of ceiba-cli,
have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

License
#######

Copyright (c) 2020-2021, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.



Credits
#######

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.
.. include:: ../README.rst
Usage
#####
The **ceibacli** command line interface offers four actions to interact
with the `Ceiba web service <https://github.com/nlesc-nano/ceiba>`_.
You can check them by trying the following command in your terminal:
::

   user>  ceibacli --help

You should see something similar to:
::

  usage: ceibacli [-h] [--version] {login,add,compute,report,query,manage} ...

  positional arguments:
    {login,add,compute,report,query,manage}
                          Interact with the properties web service
      login               Log in to the Insilico web service
      add                 Add new jobs to the database
      compute             Compute available jobs
      report              Report the results back to the server
      query               Query some properties from the database
      manage              Change jobs status

  optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit


After running one of the previous commands a log file named ``ceibacli_output.log``
is generated.
.. ceibacli documentation main file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ceibacli's documentation!
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   includereadme
   authentication
   usage
   login
   add
   compute
   report
   query
   manage

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Login
#####
We certainly want to restrict who can access and modify the data. Therefore users are required
to login with the web service. For doing so, you should have a `GitHub account <https://github.com/>`_,
then you need to request a **read-only** token from `GitHub personal access token service <https://github.com/settings/tokens>`_.

Once you have a read-only GitHub token, you can login into the web service like:
::

  ceibacli login -w http://YourCeibaInstance:8080/graphql -t Your_token

How does it work?
#################
The *Ceiba server* will contact `GitHub <https://github.com/>`_ and will check if you are a known user there.
Compute
=======
The ``compute`` command ask the *web service* for available jobs that need to be run.
To run some jobs you need to type in the terminal:
::

   ceibacli compute -i input_compute.yml

Where the *input_compute.yml* is an file in `YAML format <https://en.wikipedia.org/wiki/YAML>`_ containing the :ref:`compute input` metadata.

The compute command takes the user's input, request some available job and :ref:`schedule` those jobs using the information
provided by the user.


.. _compute input:

Compute Input File
******************

The input file contains the following mandatory keywords:
::

   # Web service URL 
   web: "http://YourCeibaInstance:8080/graphql"
   
   # Name of the collection to compute
   collection_name: "simulation_name"

   # Command use to run the workflow
   command: compute_properties
   

Other optional keywords are:
::
   
   # Configuration of the job scheduler
   scheduler:
      "none"


   # Path to the directory where the calculations are going to run (default: workdir_ceibacli)
   workdir:
      /path/to/workdir

   # Number of jobs to request and run (default: 10)
   max_jobs:
      5
      
.. _schedule:

Job Scheduling
**************
Most of the scientific simulation are usually perform in supercomputers that use a
`job scheduler <https://en.wikipedia.org/wiki/Job_scheduler>`_. *ceiba-cli* supports two of the most popular ones: `SLURM <https://www.openpbs.org/>`_.
If you choose a *scheduler* different from ``none``, *ceiba-cli* will automatically contact
the job scheduler with the options that you have provided. Below you can find a description
of the available options:
::

   # Job scheduler. Of of "none" or "slurm" (default: none)
   scheduler:
      slurm
   
   # Number of computing nodes to request (default: 1)
   nodes:
      1

   # Number of CPUs per task (default: None)
   cpus_per_task:
      48

   # Total time to request ind "days:hours:minutes" format (default: 1day)
   walltime:
     "01:00:00"

   # Partion name (queue's name) where the job is going to run (default: None)
   partion_name:
     "short"

You can alternatively provide a string with all the options for the queue system like,
::

   scheduler:
     slurm
   
   # String with user's Configuration
   free_format: "#!/bin/bash
   #SBATCH -N 1
   #SBATCH -t 00:15:00
   ....
   "


.. _Job state:

Job State
*********
The user's requested jobs are initially marked as ``RESERVERED``, in the web service to
avoid conflicts with other users. Then, if the jobs are sucessfully scheduled they
are marked as `RUNNING`. If there is a problem during the scheduling or subsequent
running step the job would be marked as `FAILED`.

Report
======
The ``report`` command send the results of the jobs computed by the user to
the web service. You can also send data that is not associated to any job to the server.
In the last case, the results don't have all the metadata associated with a job in the server,
for example because it has been previously computed or computed in another facility.

To report the results you need to type in the terminal:
::

   ceibacli report -w http://yourCeibaInstance:8080/grapqhl

Or if you want to have more control over what is reported you can provide an input file like:
::

   ceibacli report -i input_report.yml

Where the *input_compute.yml* is an file in `YAML format <https://en.wikipedia.org/wiki/YAML>`_ containing the :ref:`report input` metadata.

You can also report results without associated jobs, follow the :ref:`report stand alone results`. 

   
.. _report input:

Report results from a job
*************************
If the results that you want to report where computed with the `ceibacli compute` command, you can
optionally provide the following input:
::

   # Path to the Folder where the jobs run (default "workdir_ceibacli")
   path_results: "workdir_ceibacli"

   # Pattern to search for the result files (default "results*csv")
   output: "results*csv"

   # Pattern to search for the input files (default "inputs*json")
   input: "inputs*json"

   # If the data is already in server you can either:
   # KEEP the old data
   # OVERWRITE and discard the old data
   # MERGE the new and the old data (Default)
   # APPEND new data at the end of the old data array
   # Default = KEEP 
   duplication_policy: "KEEP"

Check the :ref:`large objects data storage` for further information on
saving large output files.

.. _report stand alone results:

Report results without associated jobs
**************************************
Sometimes you have some results that you have previously computed and you want to share them with your colleagues.
You can upload those results into the database very similarly to the previous section, but you need to
provide an additional keyword:
::

   has_metadata: False


You also need to provide the ``path_results`` and the ``output`` to look for. The ``has_metadata``
indicates to *Ceiba-cli* that the results that you want to report don't have metadata about how the
results where computed.

.. _job metadata:

How does it work?
*****************
The library enters the ``path_results`` and search recursively all the files and
directories name like ``job_*``. In each subfolder, apart from the
computed data (specificied with the ``pattern`` keyword), the ``report`` command
would try to collect the metadata associated with the job in a files named
*metadata.yml* containing the following information:
::

   job_id: 1271269411
   property:
       collection_name: awesome_data
       id: 76950

*Without the metadata no data is reported back to the server*.


Reporting data without associated jobs
**************************************


.. _large objects data storage:

Large objects data storage
**************************
For many simulation it is desirable to store the output plain data and/or the binary checkpoints.
Those files can be used to retrieve data that is not available in the database or to restart
a calculation to perform further computations.

Those large objects are not suitable for storage in a database but fortunately there are
technologies like `swift openstack <https://docs.openstack.org/swift/latest/>`_ that allows
to store these kind of data in an efficient and safely way.



In order to storage large output you need to provide in the yaml file the following keywords:
::

     large_objects:
       # URL to the datastorage service
       web: "http://large_scientific_data_storage.pi"
       # The large file(s) to search for
       patterns:  ["output*hdf5"]
       

.. Note::
   * Installing, deploying an mantaining a `swift openstack data storage service <https://docs.openstack.org/swift/latest/getting_started.html>`_ 
     is a nontrivial task. Therefore it is recommended to request access to this service to a provider.
     Be aware that **IT COSTS MONEY** to maintain the service running in a server!
   * The large files and their corresponding metadata are going to be stored in the `swift collection <https://docs.openstack.org/swift/latest/api/object_api_v1_overview.html>`_.
     using the same ``collection_name`` that has been specified in the :ref:`job metadata`.
Add (for Administrators)
########################
The ``add`` command is an adminstrative action to add new jobs into the database.

To add jobs you need to run the following command in the terminal:
::

   ceibacli add -w http://yourCeibaInstance:8080/grapqhl -c collection_name  -j Path/to/jobs.json

Where the `-w` option is the web service URL. the `collection_name` is the collection where the data is going to be stored.
Finally, the `-j` options is the path to the *JSON* file containing the jobs as an array of JSON objects.
See the next :ref:`jobs file` section for further information.

.. _jobs file:

Jobs File
*********
The job file is a list of json objects, like:
::

  [
      {
          "type": "awesome_simulation_1",
          "parameters": {
              "value": 3.14
          }
      },
      {
          "type": "awesome_simulation_2",
          "parameters": {
              "value": 2.187
          }
      }
  ]

Each job is a JSON object with the parameters to perform the simulation.
  
How does it work?
*****************
The `add` command will read each job in the JSON jobs file. For each job
it will generate a unique identifiers. Then, the jobs and their identifier will
be stored a collection named `job_your_collection_name`.
Manage (For administrators)
###########################
The ``manage`` command is an adminstrative action to change the jobs status. For example,
jobs that have been marked as ``RESERVED`` or ``RUNNING`` for a long period of time
can be marked again as ``AVAILABLE`` if the user doesn't report the results.

To change the jobs status you need to type in the terminal:
::

   ceibacli manage -i input_manage.yml

Where the *input_manage.yml* is an file in `YAML format <https://en.wikipedia.org/wiki/YAML>`_ containing the :ref:`jobs metadata` specification.

.. _jobs metadata:

Manage Input File
*****************
The following snippet represent an input example for the *manage* action:
::

   # Web service URL 
   web: "http://YourCeibaInstance:8080/graphql"  

   # Target collection to change job status
   collection_name: "example_collection"

   # Metadata to change jobs status
   change_status:
    old_status: RUNNING
    new_status: AVAILABLE
    expiration_time: 24 # one day

How does it work?
*****************
ceiba-cli will research in the ``collection_name`` for all the jobs with ``old_status`` then
it will check if those jobs have been scheduled before the ``expiration_time``. If
the jobs have expired, ceiba-cli will marked the expired jobs with the ``new_status``.
 



Authentication
##############

Generate an OAuth token
-----------------------
You need to generate an **OAuth token from GitHub** in order to login into the application!
For doing so, you should:

1. Go to `github tokens <https://github.com/settings/tokens>`_ and click the ``Generate new token`` button.
2. Provide your GitHub password when prompted
3. Fill in a description for the token, for example, *ceiba access token*.
4. **Do not enable any scope** therefore the token will grant read-only access to the app
5. Click ``Generate`` at the bottom. Make sure to copy its value because you will need it to login!

Query
=====
The ``query`` actions requests some data from the web service
and writes the requested data in a csv file.

There are currently two possible query actions:
 * request what collections are available
 * request a single collection

To request what collections are available you just need to run the following command:
::

   ceibacli query -w http://yourCeibaInstance:8080/grapqhl

Previous command will ouput something similar to:
::

   Available collections:
  name size
  simulation1 3
  simulation2 42
  ....

In the previous ``name`` indicates the actual collections' names and ``size`` how many datasets are stored
in that particular collection.

To request all the datasets available in a given collection, you just need to run the following command:
::
   
   ceibacli query -w http://yourCeibaInstance:8080/grapqhl -c simulation2

That command will write into your current work directory a file called ``simulation2.csv``
containing the properties in the requested collection.
