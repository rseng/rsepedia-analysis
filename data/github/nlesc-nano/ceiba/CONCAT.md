---
title: 'Ceiba: A web service to handle scientific simulation data'
tags:
  - Python
  - web service
  - scientific simulations
  - high performance computing
authors:
  - name: Felipe Zapata
    orcid: 0000-0001-8286-677X
    affiliation: 1
  - name: Nicolas Renaud
    orcid: 0000-0001-9589-2694
    affiliation: 1
affiliations:
 - name: Netherlands eScience Center
   index: 1

date: March 2021
bibliography: paper.bib
---

# Summary

Safe and efficient handling of large and complex data has become a critical part of many research projects.
In particular, many datasets containing the results of computationally expensive simulations are being assembled
in different scientific fields ranging from biology to material science. However, many research teams lack the proper
tools to collaboratively create, store and access datasets that are crucial for their research [@Wilson2017].
This lack of suitable digital infrastructure can lead to data loss, unreproducible results, and inefficient resources
usage that hinders scientific progress. The *Ceiba* web service provides a technical solution for teams of researchers
to jointly run the simulations needed to create their dataset, organize the data and the associated metadata, and
immediately share the datasets with each other. With *Ceiba*, academic researchers can not only improve their
data-handling practices but also promote collaboration among independent teams in need of the same data set. 

# Statement of need

Many research projects require running a large number of computationally heavy but independent simulations. 
Those can be molecular dynamics simulations of proteins [@Proteins_benchmark],
material properties [@OQMD], fluid dynamics simulations with different initial conditions [@ERCOFTAC], etc.
Recent advances in machine-learning have stimulated the creation of databases containing these types of calculations
and have highlighted the
importance of data quality and provenance as described by the FAIR data principle [@Wilkinson2016].
As the independent simulations grow in number, their orchestration and execution require a collaborative effort among a
team of researchers. Several platforms have already been developed to orchestrate large-scale collaborative efforts leveraging 
local computing resources [@SETI_at_home], [@Folding_at_home]. These platforms are technically very impressive but
offer, unfortunately, limited opportunity for reuse in smaller scale initiatives.

Here we present *Ceiba*, a light-weight library that aims to enable collaborative database creation by small and medium-sized teams.
*Ceiba* is implemented in Python using the Tartiflete GraphQL server [@Graphql;@Tartiflette].
*Ceiba* orchestrates the interaction between 3 distinct components: the client, the server and the database.
The scheme in Figure \ref{fig:architecture} represents the architecture of the web service. 

The *Ceiba web service* has been designed as a two-level service: one for managing the jobs generating
the data and another for handling the actual data. This partition allows to keep a clear boundary between
the metadata and provenance of a given job, from the concrete datasets. Since these two layers are independent,
*Ceiba* users can also manage data without associated jobs; for example,
the data has been previously computed and users just want to share it among themselves.

We use docker [@Docker] to set up and run both the server and the database in their own isolated and independent
Linux containers. The server and database containers are deployed using docker-compose [@Dockercompose] and they
communicate with each other using their own internal network [@Dockernetwork]. The docker-compose tool makes sure that
the server is listening to client requests in a given port (e.g. 8080 by default) and the database is stored
on the host computer where the docker containers are running, so it can be periodically backed up.
*Ceiba* uses MongoDB [@Mongodb] as backing database. Using a non-SQL database like MongoDB helps to
manipulate semi-structured data, like JSON files, without having to impose a schema over the simulation data.

Since both the server and the database need some computational resources to run, we anticipate that both
the server and database can be deployed at a local/national or cloud computing infrastructure. 
Once the server is up and running, users can install the client (*ceiba-cli*) on their local computer,
national computing infrastructure, cloud, etc. 

Using the client (*ceiba-cli*) the user can interact with the server and perform actions like:

* store new jobs in the database
* request some jobs to compute
* report job results
* query some available data
* perform administrative tasks on the database

Notice that in order to keep the data safe, it is required for users to log in to the Ceiba web service.
Since managing our own authentication system takes considerable time and resources, we use the
GitHub authentication system [@Authentication] to authenticate users on behalf of the Ceiba Web service. Users just
need to have a GitHub account and request a personal access token [@Token].

Once the user has authenticated with the web service, she can add new jobs by calling the client (`ceibacli add`)
with a JSON input file specifying the parameters to run the simulation. Similarly, the user can request through
the command line interface the parameters to compute new data points (`ceibacli compute`). This last command
will fetch from the server the parameters to run a specific calculation and it will feed them to the executable
provided by the user, as part of the input for the `compute` command. Finally, the client will run the job locally
or on the resources specified by the user (cluster/cloud etc.). Notice that when a user requests to compute a job,
that job is no longer available for other users and will remain in an *in-progress* state until its corresponding
results are reported or a given amount of time has passed without receiving the results. This reservation mechanism
ensures that two users do not compute the same datapoint, saving computational resources and human time.

Having run the requested jobs, the user can easily upload (using `ceibacli report`) the results and their metadata in
the server. In addition, the user can retrieve available datapoints from the database 
(using `ceibacli query`) at all times. The example section will provide a hands-on ilustration of the aforementioned actions.

Optionally, *Ceiba* allows you to store large binary/text objects using the Swift OpenStack data
storage service [@Openstack]. Large objects are not suitable for storage in a database but the 
Swift service allows to handle these kinds of objects efficiently. The only drawback of this approach
is that users need to request (and pay) for the cloud infrastructure necessary to provide this extra
service.


![Diagram representing the Ceiba architecture.\label{fig:architecture}](architecture.jpg){ width=90% }


# Examples

We present in this section a simple example illustrating the use of *Ceiba*. For a more comprehensive discussion
about how to interact with the web service, see the Ceiba-CLI documentation [@Ceiba_CLI].

## Deploying the server and the database
Before using *Ceiba* the administrator of the server, Adam, must deploy the server and database.
While in a real application both the server and database will most likely be hosted on a cloud service, we will create for the sake of illustration
a couple of containers hosted locally.

In order to start the *Ceiba server*, Adam needs to install Docker [@Docker]
and Docker-compose [@Dockercompose]. Then he needs to clone the Ceiba repository [@Ceiba]
and go to the *provisioning* folder. Inside that folder he needs to define an environmental variable defining
the MongoDB password like:
```bash
export MONGO_PASSWORD="secure_password"
```
And now he can launch the server like:
```bash
docker compose up -d
```
The previous command will launch two containers to run in the background, one with the database and the other with the server.


## Adding jobs to the database
Once the server and the database are created, Adam must specify the jobs that his collaborators will
run to compute the different data points.

For this example, we will consider a simple case where we want to compute *Pi*; using the Monte-Carlo method.
To perform the simulations we will use [a Python code script called computepi.py](https://github.com/nlesc-nano/ceiba/blob/main/paper/computepi.py).
Each job parameter is the number of *samples* to estimate *Pi*;. 

Adam must define the jobs using a JSON file that looks like:
```json
  [
    { "samples": 100 },
    { "samples": 1000 },
    { "samples": 5000 },
    ........
  ]
```

Adam can then add the jobs to the database like:
```bash
ceibacli add -w http://localhost:8080/graphql -c monte_carlo -j jobs.json
```

## Requesting job and uploading the results

Now that the database is created, Julie, a collaborator of Adam wants to request 5 jobs to compute. But before requesting the jobs,
she must first log in in to the server:

```bash
ceibacli login -t ${LOGIN_TOKEN} -w "http://localhost:8080/graphql"
```
where `LOGIN_TOKEN` is a [read-only GitHub token to authenticate the user](https://ceiba-cli.readthedocs.io/en/latest/authentication.html#authentication).
Once the authentication step is complete, Julie can request jobs and run them with the following command:

```bash
ceibacli compute -i compute_input.yml
```

where ``compute_input.yml`` is a YAML file, containing the input to perform the computation. This input file looks like:
```yml
web: "http://localhost:8080/graphql"

collection_name: "monte_carlo"

command: computepi.py

max_jobs: 5
```

This command fetches 5 available jobs from the server that still need to be computed.
These jobs will now be marked as *in progress* in the server so that other collaborators cannot compute them.
By default, Julie's job are run locally but she can also provide a ``schedule``
[section in the input file](https://ceiba-cli.readthedocs.io/en/latest/compute.html),
if she wants to run the jobs using a job scheduler like slurm [@Jette2002].

After Julie invokes the ``compute`` command the jobs will be immediately run locally or remotely.
In the background, *Ceiba-CLI* takes each job's parameters and writes them down into YAML (or JSON) format.
*Ceiba* then calls the command that Julie has provided (here ``computepi.py``). Note that all these operations are 
orchestrated by *Ceiba* and they are invisible to users.

Once the jobs have finished, Julie can upload the freshly computed datapoints to the server by executing the following command:

```bash
 ceibacli report  -w http://localhost:8080/graphql -c monte_carlo
```
The jobs executed by Julie will now be marked as `Completed` on the server. Julie and other collaborators can keep on
requesting new jobs through the `ceibacli compute` command and report those results via `ceibacli report`.
Once there are no more jobs available, *ceiba-cli* will send a message declaring that there are not more jobs available
for running.

## Querying the database

At any point all the collaborators can obtain an overview of the current status in the server via:

```bash
 ceibacli query -w http://localhost:8080/graphql
```

This will return:

```
Available collections:
  name size
monte_carlo 7
```

indicating that there is currently one datasets called *monte_carlo* and it contains 7 data points.

If users want to retreive all the available data in *monte_carlo* they can use:
```bash
 ceibacli query -w http://localhost:8080/graphql --collection_name monte_carlo
```
that will create a `monte_carlo.csv` file containing the dataset.

The example presented above is, of course, trivial and does not necessitate the collaborative efforts of multiple people.
In real-life applications, for example, each job could be a type of computationally expensive calculation
like the quantum-mechanical simulation of the molecular properties of a given structure or
the molecular-dynamics-based simulation of the docking process between two large proteins.
We hope that for such cases, where each job can require up to several days of calculation on a super-computer,
*Ceiba* can provide an easy solution to orchectrate the creation of the database and ensure its consistency.


# Acknowledgements
Felipe would like to express his deepest gratitude to Stefan Verhoeven
for guiding him on the web developing world. We are also grateful to Jen Wehner
and Pablo Lopez-Tarifa for their support and feedback designing the Ceiba web service.

# References
# What is this repository for?
This recipe installs and configures the **ceiba** app using [ansible](https://www.ansible.com/).

## Step to install the app
1. Install [Ansible](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html) in your computer.
2. Clone [this repo](https://github.com/nlesc-nano/ceiba)
3. Edit the [inventory](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html) file with the address of the server(s) where you want to install the runner.
4. Edit the [playbook](https://docs.ansible.com/ansible/latest/user_guide/playbooks.html) file with the `remote_user` name for the hosts.
5. Make sure that you can ssh to your server(s).
6. Install the runner with the following command:
   ``ansible-playbook -i inventory playbook.yml``

## Note:
The script will ask you for the `Mongodb_password` for the app.

## Backup
To backup the data in the VM a [CEPH datablock](https://doc.hpccloud.surfsara.nl/create-datablocks)
must be created and then a [crontab job](https://crontab.guru/) must be run periodically to backup
the data into the disk. First,
add a new `cron` task:
```bash
sudo crontab -e  # sudo is required to copy the mongo_data folder
```

Then add the time and the task to execute the job

```bash
00 9 * * * rsync -a /root/mongo_data /data
```

The previous command will copy the `mongo_data` folder every day at 9AM to the `/data` folder.

Here we assume that `/home/ubuntu/mongo_data` is the docker volumen where the Mongodb containers stores
the database and `/data` is the directory where the external backup disk has been mounted.


### Supported OS
Currently the recipe only works for **Ubuntu** and **Debian**
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

#. use the search functionality `here <https://github.com/nlesc-nano/ceiba/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/nlesc-nano/ceiba/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the ceiba repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
##########
Change Log
##########

1.0.0 [22/03/2021]
******************
Added
-----
* First stable versions

Changed
-------
* Add Test for Python 3.9 and drop 3.7


0.3.0 [22/02/2021]
******************
CHANGED
-------
* Generalize server by removing references to smiles and adding a metadata field (#21)

0.2.0 [10/02/2021]
******************

Added
-----
* Add query resolver to retrieve the available collections
* Add mutation to authenticate user (#2)
* Add token to identify user (#6)
* Docker container recipe (#8)
* use `Caddy <https://caddyserver.com/>`_ to generate the certificate and the start the reverse-proxy (#12)
* Use `docker-compose to start the app <https://github.com/nlesc-nano/ceiba/issues/13>`_

CHANGED
-------
* Allow to mutate the jobs timestamps and user in the jobstatus resolver
* Request a token to mutate data

0.1.0 [03/11/2020]
******************

Added
-----

* Web service prototype to handle client requests to store/retrieve simulation data (#1)
* Use `GraphQL <https://graphql.org/>`_ query languages for the API (#1)
* Use `Tartiflette <https://github.com/tartiflette/tartiflette#tartiflette-over-http>`_ Python GraphQL server implementation (#1)
* Use `MongoDB <https://www.mongodb.com/>`_ as database (#1)
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
.. image:: https://github.com/nlesc-nano/ceiba/workflows/build/badge.svg
   :target: https://github.com/nlesc-nano/ceiba/actions
.. image:: https://readthedocs.org/projects/ceiba/badge/?version=latest
   :target: https://ceiba.readthedocs.io/en/latest/?badge=latest
.. image:: https://codecov.io/gh/nlesc-nano/ceiba/branch/main/graph/badge.svg?token=MTD70XNYEA
   :target: https://codecov.io/gh/nlesc-nano/ceiba
.. image:: https://zenodo.org/badge/297567281.svg
   :target: https://zenodo.org/badge/latestdoi/297567281

#####
ceiba
#####
🧬 🔭 🔬 Scientific simulations generate large volume of data that needs to be stored and processed
by multidisciplinary teams across different geographical locations. Distributing computational expensive
simulations among the available resources, avoiding duplication and keeping the data safe are challenges
that scientists face every day.

Ceiba and its command line interface `Ceiba-cli <https://github.com/nlesc-nano/ceiba-cli>`_
solve the problem of computing, storing and securely sharing
computationally expensive simulation results. Researchers can save significant time and resources by easily
computing new data and reusing existing simulation data to answer their questions.

See `documentation <https://ceiba.readthedocs.io/en/latest/>`_ and `blog post <https://blog.esciencecenter.nl/building-a-web-service-to-manage-scientific-simulation-data-using-graphql-a0bbf1c3f6e9>`_.


Installation
************
The `provisioning folder <https://github.com/nlesc-nano/ceiba/tree/main/provisioning>`_ contains the instructions
to deploy the web service using docker. Alternatively, you can deploy the web service locally using
the following instructions:

#. 🐳 Install `Docker <https://www.docker.com/>`_

#. 🚀 Define a environment variable `MONGO_PASSWORD` with the database password. Now you can run the following
   command to start both the server and the mongodb services:
   ::

      provisioning/start_app.sh


Contributing
************

If you want to contribute to the development of ceiba,
have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

License
*******

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
*******

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.
Schema Queries
##############

See `the available queries schema definitions <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Query.graphql>`_.
Schema Mutations
################

See `the available queries schema definitions <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Mutation.graphql>`_.
Data Layout
###########
The data layouts specifies how is the data and its metadata store in the database.
Since we use `Mongodb <https://www.mongodb.com/>`_  for the database, we will
also explain the data layout using Mongobd's terminology.


Properties collections
**********************
The following schema defines how the properties are stored:
::

  # Unique identifier
  _id: Integer

  # Name to which the property belongs. e.g. Theory level
  collection_name: String

  # Metadata associated with the given property
  metadata: String
  
  # Properties values as JSON
  data: Optional[String]

  # Input with which the property was computed encoded as JSON
  input: Optional[String]

Notice that the previous schema mirros the
`GraphQL definition of Property in the server <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Query.graphql>`_.


Jobs collections
****************
The following schema defines how jobs are defined:
::
   
  # Unique identifier
  id: Integer

  # compute Properties
  property: Ref[Property]
  
  # Input to perform the computation
  settings: String

  # Job status
  status: Status

  # User who es executing the job
  user: Optional[String]

  # Timestamp = datatime.timestamp()
  schedule_time: Optional[Float]

  # Timestamp = datatime.timestamp()
  report_time: Optional[Float]

  # platform where the job was run: platform.platform()
  platform: Optional[String]
   

Notice that the previous schema mirros the
`GraphQL definition of Job in the server <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Query.graphql>`_.


.. Note::
   * `Optional[T] <https://docs.python.org/3/library/typing.html#typing.Optional>`_  is a type that could be either ``None`` or some ``T``.
   * References ``ref`` are implemented as `Mongodb DBRefs <https://docs.python.org/3/library/typing.html#typing.Optional>`_.
   * Jobs are stored in a collections named like ``jobs_<property_collection_name>``.

Queries
#######
Queries are defined using the `schema definition language <https://graphql.org/learn/schema/>`_.
You can find the mutations definition `ceiba/sdl/Query.graphql <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Query.graphql>`_

Queries involved *read-only* interactions between the client and the database.


.. automodule:: ceiba.query_resolvers
Mutations
#########
Mutations are defined using the `schema definition language <https://graphql.org/learn/schema/>`_.
You can find the mutations definition at `ceiba/sdl/Mutation.graphql <https://github.com/nlesc-nano/ceiba/blob/main/ceiba/sdl/Query.graphql>`_

Mutations involved a change in the database state, requested by the client. **All the mutations required authentication**.


.. automodule:: ceiba.mutation_resolvers
.. ceiba documentation main file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ceiba's documentation!
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   includereadme
   intro
   sdl_queries
   sdl_mutations
   queries
   mutations
   datalayout


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. include:: ../README.rst

The Ceiba Web Service
########################
🧬🧪  Most of the scientific simulations are usually performed in supercomputer
or high tech facilities. Usually the data is kept on those facilities
stored in a raw format, 🗝  in contradiction with the
`scientific FAIR principles for data <https://www.go-fair.org/fair-principles/>`_.

This repo contains a library to create a web service to interact with a database
containing a set of numerical properties. All the interactions with the database are
defined by a `GraphQL API <https://graphql.org/>`_ and the service is developed using `tartiflette <https://tartiflette.io/>`_

Adding user to the web service
##############################
In the root folder of the *ceiba* repo there is a plain text file called `users.txt`. You can add users to the
web service by adding the Github's usernames in that file.

Interactions with the database
##############################
Using the `GraphQL query language <https://graphql.org/>`_  the service
define a set of rules to interact with the services: **queries** and **mutations**.

The queries are just read only actions against the database, while the mutations,
involved some change in the database state.
