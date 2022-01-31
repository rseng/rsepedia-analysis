Install instructions
====================

Red hat 6/CentOS 6 dependencies
-------------------------------

  yum install postgresql-devel.x86_64 mod_perl ntp libXp java-1.6.0-openjdk

Octave, R, redis are available in epel repo::

  rpm -ivh http://mirror.1000mbps.com/fedora-epel/6/i386/epel-release-6-8.noarch.rpm
  yum install R redis octave
  chkconfig redis on
  /etc/init.d/redis start

Python-2.7 is needed for R integration, it is not-standard on rh6 so use compile it by hand.

As root::

  yum groupinstall "Development tools"
  yum install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel xz
  wget http://python.org/ftp/python/2.7.6/Python-2.7.6.tar.xz
  xz -d Python-2.7.6.tar.xz
  tar -xf Python-2.7.6.tar
  cd Python-2.7.6
  ./configure
  make
  make altinstall
  wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
  python2.7 ez_setup.py
  wget https://raw.github.com/pypa/pip/master/contrib/get-pip.py
  python2.7 get-pip.py
  pip install virtualenv

As application runner::

  cd $APPHOME
  virtualenv env
  . env/bin/activate
  python setup.py develop

Allow apache to reverse proxy and allow reading of static files in home dir with selinux::

  setsebool -P httpd_can_network_connect 1
  setsebool -P httpd_enable_homedirs true
  chcon -R -t httpd_user_content_t $APPHOME/trackertask/static

Appliction configuration
------------------------

  cp development.ini-dist development.ini

Or

  cp production.ini-dist production.ini

And edit it to fit needs.

Web application start/stop
--------------------------

Copy `extras/init/script-wrapper-web.conf` into `/etc/init/`.
Edit the application root and config filename variable.
Start with `start script-wrapper-web`.

Worker start/stop
-----------------

Copy `extras/init/script-wrapper-worker` into `/etc/init`.
Edit the application root and config filename variable.
Start with `start script-wrapper-worker`.

Monitoring web application and work
-----------------------------------

Nagios plugins are available in `extras/nagios`.

Apache configuration
--------------------

Create `/etc/httpd/conf.d/script-wrapper.conf` with following content or copy `extras/apache/script-wrapper-web.conf`::

  <Location /sw>
    AuthType basic
    AuthName "e-ecology credentials required"
    AuthBasicProvider file
    AuthUserFile /etc/httpd/script-wrapper.passwords
    Require valid-user
  </Location>

  ProxyPass /sw/static !
  ProxyPass /sw/ http://localhost:6543/sw/
  ProxyPassReverse /sw/ http://localhost:6543/sw/

  Alias /sw/static $APPHOME/trackertask/static

  <Directory $APPHOME/trackertask/static >
    Order allow,deny
    Allow from all
  </Directory>

Create `/etc/httpd/script-wrapper.passwords` with::

  htpasswd -c /etc/httpd/script-wrapper.passwords <username>

Enable it by

  chkconfig script-wrapper on
  /etc/init.d/script-wrapper start

Enable it by restarting apache with `/etc/init.d/httpd restart`.

See https://services.e-ecology.sara.nl/redmine/projects/uvagps/wiki/Apache_authentication_against_DB to add db authentication.

Script wrapper distributed workers
----------------------------------

To distribute work in cloud have one master machine with the web front-end and redis server.
Have multiple slaves with celery workers and local database. Script results will be shared over NFS.

1. Create cloud ubuntu machine with 25Gb disk.
2. Install dependencies
2.0 Update system + install git, build essentials + remove apache, add users myself+sw
2.1 NFS, /shared/Scratch/script-wrapper-sessions
2.2 postgresql + postgis
2.3 python, virtualenv
2.4 R, DBI, RPostgresql, stringr
2.5 Matlab MCR's
2.6 Octave
3. Install script-wrapper
3.1 Add ssh key to github + git clone
3.2 virtualenv, pip install numpy oct2py, python setup.py install
3.3 setup ini
3.4 add start/stop script

   description     "Script wrapper worker"

   start on (mounted MOUNTPOINT=/shared)
   stop on runlevel [!2345]

   setuid verhoes
   umask 022

   script
     cd /home/verhoes/eEcology-script-wrapper
     . env/bin/activate
     pceleryd development.ini
   end script

3.5 Redis server on master bind to all, open firewall `-A INPUT -i eth1 -j ACCEPT` for private network
4. Stop, rename template, start several instances.

Docker deployment
-----------------

Script wrapper consist of 3 container:
- web, instance of script-wrapper image
- redis
- worker, instance of script-wrapper image

Orchistration is done with fig (http://fig.sh).

1. Create docker image for script-wrapper
(Docker puts images in /var/lib/docker, this can be changed by starting docker deamon with `-g <graphdir>` option.)

    sudo docker build -t sverhoeven/eecology-script-wrapper:2.2.1 .
    sudo docker tag sverhoeven/eecology-script-wrapper:2.2.1 sverhoeven/eecology-script-wrapper:latest
    sudo docker login
    sudo docker push sverhoeven/eecology-script-wrapper:2.2.1
    sudo docker push sverhoeven/eecology-script-wrapper:latest

Or setup automated builds in docker registry hub, pushing commit will trigger a build. Version management needs to be done inside docker hub.

2. Setup jobs dir

    mkdir jobs
    chown www-data jobs

The script wrapper job results are shared between using the web and worker container using a host directory.
The `fig.yml` defaults to `jobs/` directory in current working directory.
Update `fig.yml` if jobs need to stored elsewhere.
The jobs dir should be writable by www-data (uid=33) user, as both web and worker service run as www-data user.

3. Start it, somewhere with docker(http://docker.com) and fig (http://fig.sh) installed

    fig -p script-wrapper up
Script wrapper
==============

Web application to make Matlab, R scripts available for data mining the UvA Bird tracking system (http://www.uva-bits.nl/ ).

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1034003.svg
   :target: https://doi.org/10.5281/zenodo.1034003
   
Available scripts:

* Tracker Calendar, Calendar overview with daily statistics of GPS-tracker
* KMZ Generator, open tracks in Google Earth
* Flight Generator, Generate a GPX file for a track. Can be used in Doarama visualization (http://www.doarama.com/).

Documentation
-------------

Documentation is available in docs/ folder.

Issue tracker
-------------

Issue tracker at https://services.e-ecology.sara.nl/redmine/projects/uvagps with category `Scriptwrapper`.

Copyrights & Disclaimers
------------------------

eEcology script wrapper is copyrighted by the Netherlands eScience Center and releases under
the Apache License, Version 2.0.

See <http://www.esciencecenter.nl> for more information on the Netherlands
eScience Center.

See the "LICENSE" and "NOTICE" files for more information.
=====
Usage
=====

Script wrapper is a web site with several pages.

.. contents:: Pages
    :local:

Home
====

List of available scripts with hyperlinks to the form of the script.

Form
====

On the form page the description of the script is given.

Each script can have a different fields, most common are a date range and a tracker selection.

To submit the form press the **Submit** button in the left bottom corner.

Save/Restore selection
----------------------

A selection is the collection of values for every field on the form.
To save a selection:

1. Press the **Save selection** button
2. Fill in a name and press **OK** button

To restore a selection:

1. Press the **Restored saved selection** button
2. In the list of selection, double-click on the name or press **Load** button next to the name

For convenience the last submitted selection is saved with the name **Last used**.

Result
======

While the script is running a status page is shown.
When the script is done the result page will be shown.
The result page consists out of a list of output files.
===============
Adding a script
===============

Adding a script is performed by the script wrapper administrator. As several files need to be added to the script wrapper application.

Adding a script consists out of several steps.

.. contents:: Steps
    :local:

1. Compiling Matlab script
==========================

To make Matlab script run without requiring a Matlab license token, it has to be compiled.

1. Copy script and all it dependencies like libraries (eg. openearth, googleearth) to a machine with a Matlab Compiler (``mcc``).
2. Load Matlab environment with ``module load matlab``. (Optional)
3. Compile it with ``mcc -mv -R -nodisplay -I <library> <script>.m [<dependency>.m]``.
4. Copy generated ``<script>`` and ``run_<script>.sh`` files back to script wrapper server.
5. The generated run script can not handle arguments with spaces. Fix ``run_<script>.sh`` by removing the lines from ``args=`` to ``done`` and use ``"$@"`` instead of ``$args``.

Example compilations:

.. code-block:: bash

    # Matlab script which uses openearth postgresql library
    mcc -mv -R -nodisplay -I openearth/io/postgresql -I openearth \
    -I openearth/general -I openearth/general/io_fun dbq.m

    # Matlab script like above, but also uses googleearth library and a dependency dist.m
    mcc -mv -R -nodisplay -I openearth/io/postgresql -I openearth \
    -I openearth/general -I openearth/general/io_fun -I googleearth gpsvis.m dist.m

2. Make skeleton
================

Depending on the language make a copy of ``script_wrapper/tasks/example_<language>`` to ``script_wrapper/tasks/<script_id>``.

In ``script_wrapper/tasks/<script_id>/__init__.py``:

* Replace ``class example_<language>`` with ``class <script_id>``.
* Fill in the fields at the beginning of the class (name, description, etc.).
* Enable script by setting ``autoregister`` to ``True``.
* By default the Matlab scripts will be compiled with Matlab 2012a, if the script has been compiled with a different version set the matlab_version property and make sure the Matlab Compile Runtime is installed and configured in the *.ini file.

3. Define form
==============

Edit ``script_wrapper/tasks/<script_id>/form.js``.

4. Validate form and map form result to script arguments
========================================================

In ``script_wrapper/tasks/<script_id>/__init__.py`` the ``formfields2taskargs`` function has to be customized.
This function recieves the form submission result.

After validation and mapping the script can be executed.

Examples of validation:

* Prevent script from running when there is no data or too much data to work with.
* Give hints which date range to choose.

Make use of colander (http://docs.pylonsproject.org/projects/colander/en/latest/index.html) for validation.

The script arguments must be JSON serializable, so don't pass objects like DateTime.

5. Run script
=============

The ``run`` function in ``script_wrapper/tasks/<script_id>/__init__.py`` has to be customized.
See documentation for ``run`` function and other scripts for examples.

The script has to be added to the ``script_wrapper/tasks/<script_id>/`` folder and ``script`` field in ``__init__.py`` has to set to it's filename.

Any binaries that the script calls should also be in ``script_wrapper/tasks/<script_id>/`` folder.

6. Run script response
======================

The script will generate one or more output files. The will be listed on the result page.

The response is a dict with the following keys:

1. query, the query object used as input for the script
2. return_code, the posix return code of a executable started. Implemented when task is subclassed from SubProcessTask or MatlabTask class.

7. Custom result page
=====================

By default the result page will consist out of a list of output files and a message whether script completed successfully.

The list of output files can be replaced with an string with html.
1. Set `result_template` in task class to a Mako template file.
2. Construct template, the variables available inside the template are:

  * `task`, that is task object or self
  * `celery` result object
  * `query`, same a run script response['query']
  * `files`, dictionary of with filename as key and url as value.

8. Add unit tests
=================

When script contains a lot of Python code write unit tests for it.
Deployment
==========

Requirements
------------

Requires a http reverse proxy which does basic authentication and passes HTTP_AUTHENTICATION and HTTP_REMOTE_USER environment variables.
By default uses a redis store to communicate between web application (pserve) and workers (pceleryd) so requires a redis store.

ExtJS (http://www.sencha.com/products/extjs) is used as widget library. 
ExtJS should be extracted in `script_wrapper/static/ext` folder. 

  cd script_wrapper/static
  wget http://cdn.sencha.com/ext/gpl/ext-4.2.1-gpl.zip
  unzip ext-4.2.1-gpl.zip
  ln -s ext-4.2.1.883 ext

Installation
------------

To install the script wrapper::

  python setup.py develop
  cp development.ini-dist development.ini

Run
---

  pserve developement.ini
  # in other shell
  pceleryd development.ini

To run as daemons::

  pserve --daemon --log-file=error.log development.ini
  nohup pceleryd development.ini --pidfile=$PWD/worker.pid -f worker.log &

Matlab compile runtime
----------------------

Install same version as used to compile.

Add postgresql jdbc (/usr/share/java/postgresql-jdbc3.jar) to
~/MATLAB/MATLAB_Compiler_Runtime/v717/toolbox/local/classpath.txt

Java is missing, add by::

   cd /home/stefanv/MATLAB/MATLAB_Compiler_Runtime/v711/sys/
   mkdir -p java/jre
   ln -s /usr/lib/jvm/java-6-openjdk-amd64 java/jre/glnxa64


R
-

To query database install from R prompt:

  install.packages('RPostgreSQL')
  install.packages('stringr')

Add cert for postgresql jdbc ssl connection
-------------------------------------------

See http://jdbc.postgresql.org/documentation/91/ssl-client.html
===================
Script requirements
===================

Requirements where the script that is added needs to adhere to.

.. contents:: Requirements
    :local:

General
=======

Several pieces of information about the script are required:

* **Id**, identifier of the script, will be used in package and url.
* **Label**, Human readable name.
* **Title**, Single line description of script.
* **Description**, Description of script.
* **Arguments** that are required to run script. Includes order, format, possible choices, db credentials and db connection string. An input form will be made based on the arguments.

Database
========

The script will be executed with the database privileges of the end-user.
This means only use table/views everyone can use eg. only use ``\*_limited`` views.

The database stores datetime with UTC timezone. Datetimes from the submit form are also in UTC timezone.
Inside script make sure you use **UTC timezone** when quering database and generating output.

The script will be called with database credentials, hostname and database name.

Matlab
======

* The script should start with a function which takes only **strings as arguments**. The reason is that the script will be compiled and started from command-line. The arguments can be converted to Matlab native variables using ``val = eval(valAsStr)``.
* Use the **OpenEarth PostgreSQL library** to perform database queries. For more information see https://services.e-ecology.sara.nl/redmine/projects/uvagps/wiki/Matlab_with_OpenEarth
* Standard out and error are saved as ``stdout.txt`` and ``stderr.txt`` resp.
* Script will be run inside the result directory so any generated output file should have no directory prefixed to it.
* Do not use hardcoded absolute paths in script, as the machine where it is compiled or being run may not have those paths.
* The script must be compiled for Linux so it can be run without a Matlab installation, to compile some information is required:

  * The version of Matlab
  * Required toolboxes
  * Additional Matlab files needed to run script.

Example Matlab script to find number of timepoints of a tracker in a certain date range:

.. code-block:: matlab

    function db_query(dbUsername, dbPassword, dbName, databaseHost,...
                      TrackerIdentifier, startTime, stopTime)
    conn = pg_connectdb(dbName, 'host', databaseHost,...
                        'user', dbUsername, 'pass', dbPassword);
    sql_tpl = ['SELECT device_info_serial, count(*) ',...
           'FROM gps.ee_tracking_speed_limited ',...
           'WHERE device_info_serial=%d ',...
           'AND date_time BETWEEN ''%s'' AND %s'' ',...
           'GROUP BY device_info_serial'];
    sql = sprintf(sql_tpl, TrackerIdentifier, starTime, stopTime)
    results = pg_fetch_struct(conn, sql);
    display(results);

R
=

* Python list variables have to converted to R vectors using eg. ``rpy2.robjects.IntVector([1,2,3])`` (For more information see http://rpy.sourceforge.net/rpy2/doc-2.2/html/introduction.html#r-vectors)
* To write output files to the ``output_dir``, it has to be passed a R function argument

Example R script to find number of timepoints of a tracker in a certain date range:

.. code-block:: r

    db_query <- function(dbUsername, dbPassword, dbName, databaseHost,
                         TrackerIdentifier, startTime, stopTime, outputDir) {
        library(RPostgreSQL)
        library(stringr)
        drv <- dbDriver("PostgreSQL")
        conn = dbConnect(drv, user=dbUsername, password=dbPassword,
                         dbname=dbName, host=databaseHost)

        tpl <- paste("SELECT device_info_serial, count(*) FROM gps.ee_tracking_speed_limited ",
           "WHERE device_info_serial="%s ",
           "AND date_time BETWEEN '%s' AND '%s' ",
           "GROUP BY device_info_serial")
        sql =  sprintf(tpl, TrackerIdentifier, startTime, stopTime)

        results <- dbGetQuery(conn, sql)

        dbDisconnect(conn)

        # Save as text
        dump('results', file=file.path(outputDir, "results.txt"))
    }

Python
======

Use SQLAlchemy models of e-ecology database.

Example Python run function to find number of timepoints of a tracker in a certain date range:

.. code-block:: python

    from script_wrapper.models import DBSession, Speed

    def run(self, db_url, tracker_id, start, end):
        # Perform a database query
        db_url = self.local_db_url(db_url)
        s = DBSession(db_url)()
        q = s.query(Speed)
        q = q.filter(Speed.device_info_serial==tracker_id)
        q = q.filter(Speed.date_time.between(start, end))
        count = q.count()

        s.close()

        # Write results to text files
        fn = os.path.join(self.output_dir(), 'result.txt')
        with open(fn, 'w') as f:
            f.write(count)
        return {'query': {'start': start,
                          'end': end,
                          'tracker_id': tracker_id,
                          }}


.. Script Wrapper documentation master file, created by
   sphinx-quickstart on Wed Sep  4 12:58:57 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Script Wrapper's documentation!
==========================================

Contents:

.. toctree::
   :maxdepth: 2

   usage
   script_requirements
   adding_script
   deployment
   api

A way make Matlab, Python and R scripts usable by wider community.
By wrapping the scripts in a web application.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

===
API
===

For developers and interested parties only.

Tasks
=====

.. automodule:: script_wrapper.tasks
    :members:

Models
======

.. automodule:: script_wrapper.models
    :members:

Validation
==========

.. automodule:: script_wrapper.validation
    :members:

Web application
===============

.. automodule:: script_wrapper
    :members:

.. automodule:: script_wrapper.views
    :members:
