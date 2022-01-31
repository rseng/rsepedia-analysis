To contribute contact me via Github issue or pull request at https://github.com/NLeSC/MAGMa.


Installation
============

Each sub-project has it own installation instructions.

Import into Eclipse
-------------------

To get the subprojects as seperate projects into Eclipse do:

1. Import 'Projects from Git'
2. Use URI with 'git@github.com:NLeSC/MAGMa.git'
3. Import all branches
4. Use default destination '~/git/MAGMa'
5. Import 'Projects from Git'
6. Use Local with '~/git/MAGMa'
7. Select import as general project
8. Open working directory tree and select 'emetabolomics_site'
9. Repeat steps 5..8 for 'job' and 'web', rename projects by prefixing 'magma'
10. On 'job' and 'web' project folder right click and select 'PyDev>Set as PyDev Project'
11. On 'job' and 'web' project folder right click and select 'Properties>PyDev PYTHONPATH'
12. Press 'Add source folder' and select project folder
13. Install virtualenv python with packages required by job and web and run 'python setup.py develop' and setup python interpreter.
14. Import Maven>'Existing Maven Projects'
15. Set Root Directory to ~/git/MAGMa/joblauncher
16. On 'joblauncher' project folder right click and select 'Team>Share project'
17. Select git select magma git repo
The eMetabolomics research project
==================================

.. image:: https://travis-ci.org/NLeSC/MAGMa.svg?branch=master
    :target: https://travis-ci.org/NLeSC/MAGMa

.. image:: https://landscape.io/github/NLeSC/MAGMa/master/landscape.svg?style=flat
    :target: https://landscape.io/github/NLeSC/MAGMa/master
    :alt: Code Health

.. image:: https://coveralls.io/repos/NLeSC/MAGMa/badge.svg?branch=master
    :target: https://coveralls.io/r/NLeSC/MAGMa?branch=master

.. image:: https://img.shields.io/badge/docker-ready-blue.svg
    :target: https://hub.docker.com/r/nlesc/magma

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1043226.svg
   :target: https://doi.org/10.5281/zenodo.1043226

The eMetabolomics project is funded by the Netherlands eScience Center and is carried out at Wageningen University and the Netherlands eScience Center in collaboration with the Netherlands Metabolomics Centre. The project develops chemo-informatics based methods for metabolite identification and biochemical network reconstruction in an integrative metabolomics data analysis workflow.

Homepage: http://www.emetabolomics.org

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.

  .. image:: web/magmaweb/static/img/metabolites.png
     :alt: Screenshot MAGMa results page

Subprojects:

- emetabolomics_site - The http://www.emetabolomics.org website
- job - Runs MAGMa calculation
- joblauncher - Webservice to execute jobs
- pubchem - Processing of PubChem database, used to find mass candidates
- web - Web application to start jobs and view results

Subproject interdependencies
----------------------------

- The `emetabolomics_site` website can be used as starting pages for the `web` application.
- The `job` calculation requires a pubchem lookup database which can be made using the `pubchem` application.
- The `web` application starts `job` calculations via the `joblauncher` webservice.

Joblauncher submodule
---------------------

Use following command to initialize and fetch the joblauncher submodule:

.. code-block:: bash

    git submodule update --init

License
-------

MAGMa is released under the Apache License Version 2.0.
The MAGMa web application uses ExtJS GPLv3 with application exception.
MAGMaWeb

By Stefan Verhoeven, Lars Ridder and Marijn Sanders.

MAGMaWeb is the web application to start new MAGMa calculations and view the results.

It's a web application with the server-side written in `Python <http://www.python.org>`_ using the `Pyramid web framework <http://www.pylonsproject.org>`_
and the client-side is written in `ExtJS <http://www.sencha.com/products/extjs>`_.

Development installation
========================

0. Requires Python 2.6 with development libraries for compilations
1. Requires ExtJS 4 for user inteface.

  1. Download ExtJS from http://www.sencha.com/products/extjs/download/ext-js-4.2.1/2281 or directy http://cdn.sencha.com/ext/gpl/ext-4.2.1-gpl.zip
  2. Unzip it in in `web/magmaweb/static`.
  3. The ext.js is missing some bits, fix by using sencha cmd

.. code-block:: bash

   cd magmaweb/static/ext-4.2.1.883
   mv ext.js ext-cmd.js
   sencha fs minify -yui -f ext-debug.js -t ext.js

Or by using yui compressor

.. code-block:: bash

   cd magmaweb/static/ext-4.2.1.883
   mv ext.js ext-cmd.js
   wget https://github.com/yui/yuicompressor/releases/download/v2.4.8/yuicompressor-2.4.8.jar
   java -jar yuicompressor-2.4.8.jar ext-debug.js -o ext.js

2. Install MAGMa web and it's dependencies

.. code-block:: bash

   virtualenv env
   . env/bin/activate
   cd web
   python setup.py develop

3. Create users and register jobs see :ref:`User management <user>` or `docs/user.rst <docs/user.rst>`_.
4. Configure and start web server

.. code-block:: bash

   cp production.ini-dist development.ini # make copy of example configuration
   mkdir data # default location; to store data elsewhere modify development.ini
   pserve development.ini --reload # starts web server

Goto http://localhost:6543/magma

Production installation
=======================

Additional to the `Development installation`_ to make application more complete/robust/faster.

* `Minimize javascript`_
* Configure reverse http proxy webserver like `nginx`_ to host static content
* Use a faster wsgi python server like `gunicorn`_ or `uWSGI`_
* Install :ref:`job launcher <launcher>` or `docs/launcher.rst <docs/launcher.rst>`_ to perform calculations
* (Optionally) Make ExtJS installation smaller by removing it's `docs`, `builds` directories

nginx
-----

`nginx web server <http://www.nginx.org>`_ is an open source Web server and a reverse proxy server for HTTP,
 SMTP, POP3 and IMAP protocols, with a strong focus on high concurrency, performance and low memory usage.

.. code-block:: bash

   sudo apt-get install nginx-full

Edit /etc/nginx/sites-enabled/default to:

.. code-block:: nginx

   server {
       #listen   80; ## listen for ipv4; this line is default and implied
       #listen   [::]:80 default ipv6only=on; ## listen for ipv6

       server_name $hostname;

       location /magma {
           proxy_set_header        Host $host;
               proxy_set_header        X-Real-IP $remote_addr;
               proxy_set_header        X-Forwarded-For $proxy_add_x_forwarded_for;
               proxy_set_header        X-Forwarded-Proto $scheme;

               client_max_body_size    1000m;
               client_body_buffer_size 128k;
               proxy_connect_timeout   60s;
               proxy_send_timeout      90s;
               proxy_read_timeout      90s;
               proxy_buffering         off;
               proxy_temp_file_write_size 64k;
               proxy_pass http://127.0.0.1:6543;
               proxy_redirect          off;
       }

       location /magma/static/ {
           alias       /home/stefanv/workspace/MAGMa/web/magmaweb/static/;
           expires     30d;
           add_header  Cache-Control public;
           access_log  off;
       }
   }

gunicorn
--------

Gunicorn wsgi server (http://gunicorn.org/) is a Python WSGI HTTP Server for UNIX.

Edit `development.ini` file by commenting out the `server:main` section with `waitress`.
And remove comment in-front of the `server:main` section with `gunicorn`.

Then start gunicorn with:

.. code-block:: bash

   pip install gunicorn
   pserve development.ini

uWSGI
-----

uWSGI wsgi server (http://projects.unbit.it/uwsgi/)  is a fast,
self-healing and developer/sysadmin-friendly application container server coded in pure C.

The HttpUwsgiModule (http://wiki.nginx.org/HttpUwsgiModule) is required.

In `production.ini-dist` there is a section for uwsgi configuration.

Change /magma section in /etc/nginx/sites-enabled/default to:

.. code-block:: nginx

    location /magma {
        proxy_set_header        Host $host;
        proxy_set_header        X-Real-IP $remote_addr;
        proxy_set_header        X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header        X-Forwarded-Proto $scheme;

        client_max_body_size    1000m;
        client_body_buffer_size 128k;
        include uwsgi_params;
        uwsgi_pass unix:/tmp/magma.uwsgi.sock;
        uwsgi_param SCRIPT_NAME /magma;
        uwsgi_modifier1 30;
        uwsgi_param  UWSGI_SCHEME   $scheme;
    }

Then start uWSGI with:

.. code-block:: bash

   pip install uwsgi
   uwsgi -H env --ini-paste-logged development.ini

Minimize javascript
-------------------

Install Sencha SDK tools by following instructions at http://www.sencha.com/products/sencha-cmd use version 4.x.x.
Direct download https://cdn.sencha.com/cmd/4.0.5.87/SenchaCmd-4.0.5.87-linux-x64.run.zip.

Then concatenate and compress javascript with:

.. code-block:: bash

   cd magmaweb
   sencha build -d static/app -p magmaweb.results-4.2.1.jsb3
   ln -s magmaweb/static/app/resultsApp-all-4.2.1.js magmaweb/static/app/resultsApp-all.js

Now not hundreds of seperate javascript files are loaded, but a single javascript file.

Create magmaweb.results.jsb3 file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This only needs to be done if magmaweb.results*.jsb3 does not yet create.

The `sencha create` command does not work for our pages. So we role our own jsb3 writer.

1. Load result page.
2. Goto developers/firebug console
3. Enter `copy(Ext.Loader.history)`
4. Open file `myhistory` and paste clipboard (CTRL-p)
5. Load workspace page
6. Goto developers/firebug console
7. Enter `copy(Ext.Loader.history)`
8. Open file `myhistory` for appending and paste clipboard (CTRL-p)
9. Run `perl loader2jsb3.pl myhistory > magmaweb.results-4.2.1.jsb3`

loader2jsb3.pl looks like:

.. code-block:: perl

   #!/usr/bin/env perl

   use strict;
   use warnings;
   use JSON;

   my %paths = (
      'Ext' => 'static/ext-4.2.1.883/src',
      'Ux'  => 'static/ext-4.2.1.883/examples/ux',
      'Esc' => 'static/esc',
      'App' => 'static/app'
   );
   my @files;
   my %cache;

   while (<>) {
     my $line = $_;
     chomp($line);
     for my $dep (split(/,/,$line)) {
       if ($cache{$dep}) {
         next;
       } else {
         $cache{$dep}++;
       }
       my ($path, $name) = $dep =~ /(.*)\.(.*)/;
       $name .= '.js';
       $path =~ s/\./\//g;
       $path .= '/';
       if ($path=~/^Esc\/magmaweb/) {
           $path =~ s/^Esc\/magmaweb/$paths{App}/;
       } elsif ($path=~/^Esc/) {
           $path =~ s/^Esc/$paths{Esc}/;
       } elsif ($path=~/^Ext\/ux/) {
           $path =~ s/^Ext\/ux/$paths{Ux}/;
       } else {
   	       $path =~ s/^Ext/$paths{Ext}/;
       }
       push(@files, {'path'=> $path, 'name'=> $name});
     }
   }

   print to_json({
     'projectName'=> 'MAGMA web results',
     licenseText=> "Copyright(c) 2011 Netherlands eScience Center",
       "builds"=> [
           {
               "name"=> "All Classes",
               "target"=> "resultsApp-all-4.2.1.js",
               "compress"=> JSON::true,
               "files"=> \@files
   }
       ],
       "resources"=> []
   }, {pretty=>1});

Running tests
=============

Python
------

Python tests can be run with:

.. code-block:: bash

   pip install nose coverage
   nosetests

To run only unit tests:

.. code-block:: bash

   nosetests -a '!functional'

To run only functional tests:

.. code-block:: bash

   nosetests -a functional

Javascript
----------

The ExtJS tests can be run using karma runner (http://karma-runner.github.io/).
The tests require NodeJS v6.

.. code-block:: bash

    npm install -g karma-cli
    npm install
    karma start

It will generate JUnit XML files as `TEST-*.xml` and a coverage report in coverage/ directory.

Generate documentation
======================

Python
------

Generate Python documentation with

.. code-block:: bash

   pip install sphinx
   cd docs
   make html
   firefox _build/html/index.html

Javascript
----------

Javascript documentation generation with JSDuck.
See https://github.com/senchalabs/jsduck

.. code-block:: bash

   jsduck magmaweb/static/ext-4.2.1.883/src magmaweb/static/ext-4.2.1.883/examples/ux \
   magmaweb/static/d3/d3.min.js magmaweb/static/esc magmaweb/static/app --builtin-classes \
   --output jsdoc --images magmaweb/static/ext-4.2.1.883/docs/images
   firefox jsdoc/index.html

Database migration
==================

When `magmaweb/models.py` is changed then all the databases have to migrated to this new state.
Alembic (http://readthedocs.org/docs/alembic/) is used to perform database migrations.

When `models.py` has changed use ``alembic -x jobid=ff52323b-c49a-4387-b964-c6dafab5f0c4 revision --autogenerate -m "Added metabolize scenario"`` to make a migration script.
You might need to force the database to the head revision using ``alembic -x jobid=e1e4951e-30e1-4ce7-b0e1-e8af6b998580 stamp 185259a481ee``.

Upgrade all the job result databases with:

.. code-block:: bash

    for x in `ls data/jobs`
    do
    echo $x
    alembic -x jobid=$x upgrade head
    done

The migration version of a job db can be queried with ``alembic -x jobid=ff52323b-c49a-4387-b964-c6dafab5f0c4 current``.
====================
Web service consumer
====================

MAGMa has a nice web user interface.
The user interface recieves/sends data using a JSON web service.
The JSON web service can also be used by other applications.

The Web service api documentation is available at http://nlesc.github.io/MAGMA-apiary/ or at http://docs.magma.apiary.io/

Below is an example how to use the web service using `curl <http://curl.haxx.se/>`_, a command line HTTP client.

.. contents::

Start session
=============

.. code-block:: bash

   curl -L -c cookie.jar -b cookie.jar http://www.emetabolomics.org/magma

Returns the html for the MAGMa home page.
The output can be ignored. This was just to get a session cookie into the cookie.jar.

Submit job
==========

.. code-block:: bash

   curl -L -c cookie.jar -b cookie.jar \
   -F ms_data_file=@input.mtree -F ms_data_format=mass_tree \
   -F structure_database=pubchem \
   -F max_mz=1200 -F min_refscore=1 \
   -F ionisation_mode=-1 \
   -F max_broken_bonds=3 -F max_water_losses=1 \
   -F mz_precision=5 \
   http://www.emetabolomics.org/magma/start

Parameters:

- ``ms_data_file``, file upload of mass spectra data
- ``ms_data_format``, which format ``ms_data_file`` is in. See http://www.emetabolomics.org/magma/help for examples.
- ``structure_database``, Retrieve molecules from a database. Can be pubchem, kegg, hmdb.
- ``max_mz``, Mass filter for structure database
- ``min_refscore``, PubChem reference score filter for pubchem structure database.
- ``ionisation_mode``, -1 or 1
- ``max_broken_bonds``, Bond breaks
- ``max_water_losses``, Additional water losses
- ``mz_precision``, Relative precision (ppm)

Example response:

.. code-block:: javascript

   {
      "success": true,
      "id": "844bcea5-058b-4b7f-8d29-ba2cc131a568"
   }

Making job public
=================

By default, jobs can only be seen by the user that submitted it.
An additional command is needed to make it visible for anyone (who knows the url).

Query file (query.json):

.. code-block:: javascript

   {
      "description": "New description",
      "ms_filename": "New file name for MS data",
      "is_public": true
   }

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -d @query.json -X PUT http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568

http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568 can now be shared and shown in a web-browser.
When job is not yet completed it will show a status page, after completion the results will be automatically shown.

Poll status
===========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar http://www.emetabolomics.org/magma/status/844bcea5-058b-4b7f-8d29-ba2cc131a568.json

Where ``844bcea5-058b-4b7f-8d29-ba2cc131a568`` is the job identifier returned by the job submission.

Retry until job has status STOPPED.

Example response:

.. code-block:: javascript

   {
      "status" : "STOPPED",
      "jobid" : "844bcea5-058b-4b7f-8d29-ba2cc131a568"
   }

Fetching results
================

Molecules
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/molecules.json?start=0;limit=10'

Parameters:

- ``start``, Offset in list of molecules
- ``limit``, Maximum nr of molecules to return
- ``scanid``, only return molecules that have hits in scan with this identifier (optional)

Example response:

.. code-block:: javascript

   {
      "totalUnfiltered": 1,
      "total": 1,
      "rows": [{
         "name": "CHLOROGENIC ACID (1794427)",
         "smiles": "CWVRJTMFETXNAD",
         "refscore": 10234.0,
         "formula": "C16H18O9",
         "assigned": false,
         "reference": "hyperlink ....",
         "mol": "molblock ....",
         "reactionsequence": [],
         "predicted": true,
         "mim": 354.095082,
         "logp": -0.4,
         "level": 1,
         "molid": 1,
         "nhits": 1
      }],
      "scans": [{"rt": null, "id": 1}]}

Long values have been replaced with `....`.

Fields:

- ``molid`` is the molecule identifier.
- ``name`` is the name of the molecule.

Molecules for one scan
----------------------

To fetch a ranked list of molecules which are annotated for a certain scan.

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/molecules.json?start=0;limit=10;scanid=1;sort=%5B%7B%22property%22%3A%22score%22%2C%22direction%22%3A%22ASC%22%7D%5D'

``%5B%7B%22property%22%3A%22score%22%2C%22direction%22%3A%22ASC%22%7D%5D``
is the URL encoded (see http://www.faqs.org/rfcs/rfc3986) version of
``[{"property":"score","direction":"ASC"}]`` and orders the molecules with the highest Candidate score first.

Same response as above, but with additional ``score`` and ``deltappm`` fields.

Fragments
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/fragments/123/456.json?node=root'

Where ``123`` is the scan identifier and ``456`` is the molecule identifier.

Parameters:

- ``node``, The fragment identifier to fetch children fragments for.

Example response:

.. code-block:: javascript

   {
      "expanded" : true,
      "children" : [
         {
            "deltah" : -1,
            "deltappm" : -0.8824098991817264,
            "mol" : "molblock ....",
            "formula": "C16H17O9",
            "molid" : 23,
            "fragid" : 5,
            "score" : 3,
            "mass" : 370.1263823051,
            "scanid" : 1789,
            "expanded" : true,
            "mz" : 369.119262695312,
            "mslevel" : 1,
            "atoms" : "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15",
            "isAssigned" : false,
            "leaf" : false,
            "children" : [
               {
                  "deltah" : -2,
                  "deltappm" : -1.861685339415437,
                  "mol" : "molblock ....",
                  "formula" : "C7H11O6",
                  "molid" : 23,
                  "fragid" : 6,
                  "score" : 2,
                  "mass" : 115.039519091,
                  "scanid" : 1790,
                  "expanded" : true,
                  "mz" : 113.024360656738,
                  "mslevel" : 2,
                  "atoms" : "14,15,16,20,22,23,24,25",
                  "leaf" : true
               }
            ]
         }
      ]
   }

Fields:

- ``fragid`` is the fragment identifier.
- ``molid`` is the molecule identifier.
- ``scanid`` is the scan identifier.

Chromatogram
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/chromatogram.json'

Example response:

.. code-block:: javascript

   {
      "cutoff": 0.0,
      "scans": [{
         "rt": null,
         "ap": 0,
         "intensity": 69989984.0,
         "id": 1
      }]
   }

Fields:

- ``rt`` is the retention time.
- ``ap`` whether scan has molecules assigned to peaks
- ``id`` is the scan identifier.

Mass spectra
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/mspectra/1234.json'

Where ``1234`` is the scan identifier.

Example response:

.. code-block:: javascript

   {
      "precursor": {
         "id": 0,
         "mz": 0.0
      },
      "cutoff": 0.0,
      "peaks": [{
         "intensity": 69989984.0,
         "assigned_molid": null,
         "mz": 353.087494
      }],
      "mslevel": 1,
      "fragments": [{"mz": 1353.087494}]
   }

Fields:

- ``precursor``, The precursor scan identifier and mz of current scan.
- ``peaks``, list of peaks for current scan.
- ``fragments``, Peaks which have one or more fragments.

Extracted ion chromatogram
--------------------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/extractedionchromatogram/1234.json'

Where ``1234`` is the molecule identifier.

Example response:

.. code-block:: javascript

   {
      "chromatogram": [{
         "rt": null,
         "intensity": 69989984.0
      }],
      "scans": [{
         "rt": null,
         "id": 1
      }]
   }

Fields:

- ``chromatogram`` is list of intensities of selected molecule for each retention time.
- ``scans`` is list of scan identifiers where selected molecule was fragmented.

Delete job
==========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -X DELETE 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568'
.. _launcher:

Job launcher
============

To perform MAGMa calculations the web application offloads the job to a launcher deamon.

The job launcher is in the `joblauncher/` directory of the MAGMa repository.

The example configuration expects a joblauncher deamon running on `http://localhost:9998/job`.

The status of a job managed by the job launcher is pushed back to the web application using a callback url.

The joblauncher is authenticated against the web application using MAC Access Authentication.
See http://tools.ietf.org/html/draft-ietf-oauth-v2-http-mac-01 .

Installation
------------

1. Create a user (e.g. "joblauncher") for the joblauncher, see :ref:`User management <user>` or `user.rst <user.rst>`_.
2. Set 'monitor_user' key in web config file, so that user is authorized with to update statuses of jobs;

    monitor_user = joblauncher

3. Login as the joblauncher user on web site
4. On workspace page generate an access token. Example token object:

.. code-block:: javascript

    {
        "acesss_token": "bGFpYmFldGhlZXRoZWV0aGFwaGVpZm9vZmFlc2hvcmVpbW9oamluZ2lheG9jaGV6b3VwZW92YWVzaGVhYmFobmdvb3F1ZWlib2thaG5nZWV0ZWVwaG9odGhldXR=",
        "mac_key": "dW9oYml1cGllZ2hhaWNhdWZvaAo=",
        "expires_in": 30758400, # a year
        "token_type": "mac",
        "mac_algorithm": "hmac-sha-1"
    }

5. Configure joblauncher with supplied token.

* Use "access_token" from token object in joblauncher config as "id"
* Use "mac_key" from token object in joblauncher config as "key"
* Use "http(s):\\host:port" where MAGMa web is running at as "scope" in joblauncher config. Examples:
  * if MAGMaWeb is running at http://localhost:6543/magma use http://localhost:6543 as scope
  * if MAGMaWeb is running at https://www.emetabolomics.org/magma use https://www.emetabolomics.org/magma as scope

The joblauncher config will look like:

.. code-block:: yaml

   macs:
   - id: bGFpYmFldGhlZXRoZWV0aGFwaGVpZm9vZmFlc2hvcmVpbW9oamluZ2lheG9jaGV6b3VwZW92YWVzaGVhYmFobmdvb3F1ZWlib2thaG5nZWV0ZWVwaG9odGhldXR=
     key: dW9oYml1cGllZ2hhaWNhdWZvaAo=
     scope: https://www.emetabolomics.org
.. _user:

User & job management
=====================

Users and jobs can be added/edited/deleted using the `magma-web` script. For options, run:

.. code-block:: bash

   magma-web -h
   magma-web add -h

Anonymous and restricted mode
=============================

To make the MAGMa application usable by everyone,
MAGMa can be configured to run in anonymous and restricted mode.

Also for certain journals it is required that you can use the application without creating an account.

Anonymous mode
--------------

To enable put in the `*.ini` config file:

    auto_register = true

If you go to the application you will not be asked to login.
You still have a workspace and make results public if you want.

Restricted mode
---------------

To enable put in the `*.ini` config file:

    restricted = true

This applies restrictions to the calculations to make sure calculations run fast.
.. MAGMaWeb documentation master file, created by
   sphinx-quickstart on Tue Nov 15 15:35:20 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MAGMaWeb's documentation!
====================================

By Stefan Verhoeven, Lars Ridder and Marijn Sanders.

MAGMaWeb is the web application to start new MAGMa calculations and view the results.

It’s a web application with the server-side written in Python using the Pyramid web framework and the client-side is written in ExtJS.

.. toctree::
   :maxdepth: 1

   Installation <install>
   Job Launcher <launcher>
   User management <user>
   Web service consumer <ws_consumer>
   API Documentation <api>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: ../README.rstAPI Documentation
=================

.. automodule:: magmaweb
    :members:

:mod:`magmaweb.views`
---------------------

.. automodule:: magmaweb.views
    :members:

:mod:`magmaweb.rpc`
-------------------

.. automodule:: magmaweb.rpc
    :members:

:mod:`magmaweb.models`
----------------------

.. automodule:: magmaweb.models
    :members:

:mod:`magmaweb.user`
--------------------

.. automodule:: magmaweb.user
    :members:

:mod:`magmaweb.job`
-------------------

.. automodule:: magmaweb.job
    :members:

:mod:`magmaweb.jobquery`
------------------------

.. automodule:: magmaweb.jobquery
    :members:
Create local PubChem and HMDB structure databases
=================================================

PubChem is a large database of chemical molecules and their activities against biological assays.
HMDB is a database of human metabolites
For more information see http://pubchem.ncbi.nlm.nih.gov/ and http://www.hmdb.ca/

Create local PubChem and Kegg (=subset of PubChem) databases:
-------------------------------------------------------------

1) Download/update local copy of PubChem data files:

.. code-block::

   usage: process_pubchem.py update [-h] [-d DATA_DIR]
 
   Update downloads from PubChem server
   
   optional arguments:
     -h, --help            show this help message and exit
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where PubChem data is stored (default: ./)

2) Manual download of file with information of Kegg compounds in PubChem.
   Must be downloaded from http://pubchem.ncbi.nlm.nih.gov/ as follows:
   Substance => Advanced => SourceName: kegg => Search;
   Display Settings: Format = ID Map => Apply;
   Send to: File

3) Process the PubChem data twice, once for non-halogenated compounds and
   once for halogenated compounds (-f option)

.. code-block::

   usage: process_pubchem.py process [-h] [-k KEGG] [-s] [-f] [-d DATA_DIR]
                                     [-b DATABASE_DIR]
   
   Update local PubChem databases. This has to be run twice, with and without -f option
   to generate db's for halogenated and non-halogenated compounds respectively
   
   optional arguments:
     -h, --help            show this help message and exit
     -k KEGG, --kegg KEGG  File obtained in step 2
     -s, --skip_names      Skip update of PubChem names db (default: False)
     -f, --halogens        Generate database with halogenated compounds (default:
                           False)
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where PubChem data is stored (default: ./)
     -b DATABASE_DIR, --database_dir DATABASE_DIR
                           Directory where PubChem databases are stored (default:
                           ./)

Create local HMDB database:
---------------------------

.. code-block::

   usage: process_hmdb.py [-h] [-v] [-d DATA_DIR] [-b DATABASE_DIR]

   Update local HMDB database
   
   optional arguments:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where HMDB structures.zip file is stored
                           (default: None, attempt to read directly from HMDB server)
     -b DATABASE_DIR, --database_dir DATABASE_DIR
                           Directory where HMDB database is stored (default: ./)
   MAGMa
=====

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.
MAGMa subproject which performs actual calculations.

Docker
------

The MAGMa command-line script can be executed in a Docker container. For more information on installing and running Docker see https://www.docker.com/.
With the Docker image it is particularly straightforward to run the "light" subcommand which takes a single MS/MS spectrum as input, retrieves candidate structures from the online databases or a local file and writes the result as smiles or sdf to stdout.

Examples:

.. code-block:: bash

   $ docker run --rm nlesc/magma light -h
   $ cat >glutathione.mgf 
   BEGIN IONS
   TITLE=CASMI 2014, Challenge 9
   PEPMASS=308.0912 100.0
   116.0165 3.2
   144.0114 6.3
   162.0219 40.2
   179.0485 100.0
   233.0590 21.6
   290.0802 5.1
   END IONS
   ^d
   $ docker run --rm -v $PWD:/data nlesc/magma light -f mgf -s hmdb glutathione.mgf

Installation
------------------------

.. code-block:: bash

   # Install conda
   # For python 2:
   wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
   sh Miniconda2-latest-Linux-x86_64.sh
   # Or, for python 3:
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   
   # Optionally, create a dedicated conda environment and activate
   conda create -n magma
   source activate magma
   
   # Install dependencies
   conda install -c rdkit rdkit
   conda install cython lxml nose coverage
   # For python 2:
   pip install http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip
   # For python 3:
   pip install https://www.parallelpython.com/downloads/pp/pp-1.6.4.4.zip
   
   # If needed install C compiler
   sudo apt-get update && sudo apt-get install gcc
   
   # Install MAGMa
   git clone https://github.com/NLeSC/MAGMa.git
   cd MAGMa/job
   python setup.py develop

Configuration
-------------

A 'magma_job.ini' config file is read from users home directory (~/).

Exampe config file to read candidate molecules from the emetabolomics server:

.. code-block:: INI

   [magma job]
   # Retrieve candidate molecules from the emetabolomics server
   structure_database.online = True
   structure_database.service = http://www.emetabolomics.org/magma/molecules

Example config file to read candidate molecules from local databases (can be created by the scripts in MAGMa/pubchem):

.. code-block:: INI

   [magma job]
   # Location of structure database from which to retrieve candidate molecules locally
   structure_database.online = False
   structure_database.pubchem = /home/user/magma_databases/Pubchem_MAGMa.db
   structure_database.pubchem_halo = /home/user/magma_databases/Pubchem_MAGMa_halo.db
   structure_database.kegg = /home/user/magma_databases/Pubchem_MAGMa_kegg.db
   structure_database.kegg_halo = /home/user/magma_databases/Pubchem_MAGMa_kegg_halo.db
   structure_database.hmdb = /home/user/magma_databases/HMDB_MAGMa.db

   # MACS authentication, used for sending progress reports to MAGMa web application
   macs.id = <MAC key identifier>
   macs.key = <MAC key>

Usage
-----

Annotate a tree file using PubChem database:

.. code-block:: bash

   echo '353.087494: 69989984 (191.055756: 54674544 (85.029587: 2596121, 93.034615: 1720164, 109.029442: 917026, 111.045067: 1104891 (81.034691: 28070, 83.014069: 7618, 83.050339: 25471, 93.034599: 36300, 96.021790: 8453), 127.039917: 2890439 (57.034718: 16911, 81.034706: 41459, 83.050301: 35131, 85.029533: 236887, 99.045074: 73742, 109.029404: 78094), 171.029587: 905226, 173.045212: 2285841 (71.013992: 27805, 93.034569: 393710, 111.008629: 26219, 111.045029: 339595, 137.024292: 27668, 155.034653: 145773), 191.055725: 17000514), 353.087097: 4146696)' > example.tree
   magma read_ms_data --ms_data_format tree -l 5 -a 0  example.tree results.db
   magma annotate -p5 -q0 -c0 -d0 -b3 -i -1 -s pubchem -o ../pubchem/Pubchem_MAGMa_new.db,0,9999 -f results.db

Running on cluster
------------------

On the compute node not all dependencies of Magma will be installed.
By freezing the magma application on the head node we include all dependencies like rdkit.

On head node:

.. code-block:: bash

   pip install bbfreeze
   python setup.py bdist_bbfreeze
   cd dist
   chmod +x dist/Magma-<version>/Magma-<version>-py2.7.egg/magma/script/reactor
   tar -zcf Magma-<version>.tar.gz Magma-<version>

On compute node:

.. code-block:: bash

   tar -zxf Magma-<version>.tar.gz
   ./Magma-<version>/magma ...
