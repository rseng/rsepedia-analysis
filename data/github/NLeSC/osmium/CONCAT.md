Osmium
======

|TravisCILink|_

.. |TravisCILink| image:: https://travis-ci.org/NLeSC/osmium.png
.. _TravisCILink: https://travis-ci.org/NLeSC/osmium
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1065436.svg
   :target: https://doi.org/10.5281/zenodo.1065436

Web service to submit jobs via a Xenon (https://nlesc.github.io/Xenon) supported scheduler.

Submit job to web service using a HTTP POST request.
Reports status of job to submitter using a callback url.

Requirements
------------

- JDK 7 (http://www.java.com)
- Maven 3 (http://maven.apache.org)

Install
-------

#. Make copy of 'joblauncher.yml-dist' to 'joblauncher.yml'

  #. Configure Xenon scheduler
  #. Configure Xenon sandbox root directory
  #. Configure optional MAC id/key

2. Build uber-jar or execute from maven.

2.1. Uber-jar, to start on other machine the `osmium-2.0.jar` and `joblauncher.yml` files must be copied.

.. code-block:: bash

   mvn package
   java -jar target/osmium-2.0.jar server joblauncher.yml


2.2. Execute from maven

.. code-block:: bash

   mvn compile exec:java

Usage
-----

Web service listens on http://localhost:9998 .

Create a directory with an input file and script

.. code-block:: bash

   mkdir myjob
   cd myjob
   echo 'Lorem ipsum' > input_file
   echo 'hostname;date;wc -l input_file > output_file' > runme.sh

Create a json file (query.json)

.. code-block:: json

   {
      "jobdir": "<absolute path to myjob directory>",
      "executable": "/bin/sh",
      "prestaged": ["runme.sh", "input_file"],
      "poststaged": ["output_file"],
      "stderr": "stderr.txt",
      "stdout": "stdout.txt",
      "arguments": ["runme.sh"],
      "status_callback_url": "<optional url where job status is PUT to>"
   }

The `status_callback_url` value is the url where the web service will send the job status with a http PUT request.
Can be used as an alternative for polling the job status.

Then submit it

.. code-block:: bash

   curl -H "Content-Type: application/json" -H 'Accept: application/json' -i -X POST -d @query.json http://localhost:9998/job

   HTTP/1.1 201 Created
   Date: Thu, 23 May 2013 11:50:28 GMT
   Location: http://localhost:9998/job/local-1234
   Content-Type: application/json
   Content-Length: 0

The submit response contains no content only headers.
The `Location` header value is the url where the job can be queried for it's status or where it can be canceled.

Callback authentication
^^^^^^^^^^^^^^^^^^^^^^^

The status callbacks uses MAC Access Authentication.
The MAC key indentifier and MAC key must be obtained from the provider.

Status
^^^^^^

In the submit response the url is a relative url to the job.

.. code-block:: bash

   curl -H "Content-Type: application/json" -H 'Accept: application/json' http://localhost:9998/job/local-1234

Example response when job is running:

.. code-block:: json

   {
       "request": {
           "jobdir": "/tmp/jobdir",
           "executable": "/bin/sh",
           "stderr": "stderr.txt",
           "stdout": "stdout.txt",
           "arguments": [
               "runme.sh"
           ],
           "prestaged": [
               "runme.sh", "input.dat"
           ],
           "poststaged": ["output.dat"],
           "status_callback_url": "http://localhost/status"
       },
       "status": {
         "state": "RUNNING",
         "exitCode": null,
         "exception": null,
         "running": true,
         "done": false,
         "schedulerSpecficInformation": null
      }
   }

Example response when job is done:

.. code-block:: json

   {
       "request": {
           "jobdir": "/tmp/jobdir",
           "executable": "/bin/sh",
           "stderr": "stderr.txt",
           "stdout": "stdout.txt",
           "arguments": [
               "runme.sh"
           ],
           "prestaged": [
               "runme.sh", "input.dat"
           ],
           "poststaged": ["output.dat"],
           "status_callback_url": "http://localhost/status"
       },
       "status": {
         "state": "DONE",
         "exitCode": 0,
         "exception": null,
         "running": false,
         "done": true,
         "schedulerSpecficInformation": null
      }
   }

Example response when job has been canceled (see below for cancel command):

.. code-block:: json

   {
      "request": {
         "jobdir": "/tmp/myjob",
         "status_callback_url": null,
         "poststaged": [
            "output_file"
         ],
         "stderr": "stderr.txt",
         "executable": "/bin/sh",
         "arguments": [
            "runme.sh"
         ],
         "prestaged": [
            "runme.sh",
            "input_file"
         ],
         "stdout": "stdout.txt"
      },
      "status": {
         "running": false,
         "done": true,
         "exception": "Process cancelled by user.",
         "schedulerSpecficInformation": null,
         "exitCode": null,
         "state": "KILLED"
      }
   }

Cancel
^^^^^^

Cancel a pending or running job.
Deletes any generated output in the sandbox where the job was running.

.. code-block:: bash

   curl -H "Content-Type: application/json" -H 'Accept: application/json' -X DELETE http://localhost:9998/job/local-1234

Documentation
-------------

A maven site can be generated with

.. code-block:: bash

   mvn site
   firefox target/site/index.html

Integration tests
-----------------

Run integration tests with

.. code-block:: bash

   mvn verify
