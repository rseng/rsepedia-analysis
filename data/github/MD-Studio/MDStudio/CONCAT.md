# MDStudio
[![Build Status](https://travis-ci.org/MD-Studio/MDStudio.svg?branch=master)](https://travis-ci.org/MD-Studio/MDStudio)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aff6245274f44a7991a3a25976ad6472)](https://www.codacy.com/app/tifonzafel/MDStudio?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=MD-Studio/MDStudio&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/MD-Studio/MDStudio/branch/master/graph/badge.svg)](https://codecov.io/gh/MD-Studio/MDStudio)

![Configuration settings](docs/img/mdstudio-logo.png)

MDStudio is a stand-alone framework to build distributed applications using a broker-based ([Crossbar](https://crossbar.io)) microservice concept.
The only requirement for any piece of software or script to become a microservice is a formal definition of the functional endpoints that expose
the functionality the software has to offer. In MDStudio, this definition is a very thin layer available in nearly 12 different [languages](https://crossbar.io/about/Supported-Languages/)
in which communication is network based using JSON and the input/output of the endpoints is fully described using JSON schemas.



## Installation

To make development and deployment easier we have setup a docker environment hosting the microservice broker and support microservices in a self
consistent environment.

> **NOTE**: some of the Python dependencies of MDStudio are version fixed. To prevent accidental updates that may prevent MDStudio
  from working we advise to install MDStudio and it's dependencies in a virtual environment such as `virtual env` or Anaconda/miniconda
  virtual environment.

### Quick-start docker based installation

The only requirement MDStudio has is [Docker](https://www.docker.com/).
Next, run the builder and start the microservice docker environment under bash as:

	Clone the repository using the devel branch
	>> git clone -b devel git@github.com:MD-Studio/MDStudio.git
	Install the dependencies from the MDStudio folder
	>> pip install -r requirements-dev.txt
	Install mdstudio
	>> pip install -e mdstudio
    >> ./build.sh
    >> ./standalone.sh

The docker containers only need to be build once and for subsequent usage running
standalone.sh is sufficient.

### Manual (non-docker) based installation

The docker based environment is the most convenient way for developing microservices and deploying distributed applications.
If you are unable to use Docker or want to keep track of the whole environment for yourself, you should follow these installation
instructions.

The MDStudio application is written in Python and mostly self contained thanks to the use of an in-application Python virtual environment.
The application has currently been successfully tested with Python versions: 2.7

The only external dependencies are:

 * [MongoDB](https://www.mongodb.com) - A NoSQL database.
 * [Pipenv](https://github.com/kennethreitz/pipenv) - A python virtual environment manager.
 * [Redis](https://redis.io/) - A fast caching layer.

Install the virtual environment with:

    >> pipenv install --skip-lock --dev

The application is started on the command line as:

    >> pipenv shell
    >> python .

#### PyCharm IDE Integration
Go to `File > Project Settings > Project Interpreter`, and add a remote interpreter,
and make sure it matches this screen.

![Configuration settings](docs/img/pycharm-config.png)

Note specifically:

|                      |                                                            |
|----------------------|------------------------------------------------------------|
| Interpreter path     | `/root/mdstudio/app-4PlAip0Q/bin/python`   |


## Creating microservices

Now that you have successfully installed the MDStudio framework you can start having fun by writing some
microservices and using their endpoints in any number of different ways as briefly described in the intro.

For guidance and inspiration on how to write and use microservices, please have a look at the [MDStudio examples](https://github.com/MD-Studio/MDStudio_examples)
repository on GitHub.
# docker-redis-cluster

[![Docker Stars](https://img.shields.io/docker/stars/grokzen/redis-cluster.svg)](hub])
[![Docker Pulls](https://img.shields.io/docker/pulls/grokzen/redis-cluster.svg)](hub])

Docker image with redis built and installed from source.

The main usage for this container is to test redis cluster code. For example in https://github.com/Grokzen/redis-py-cluster repo.

The cluster is 6 redis instances running with 3 master & 3 slaves, one slave for each master. They run on ports 7000 to 7005.

It also contains 2 standalone instances that is not part of the cluster. They are running on port 7006 & 7007

This image requires at least `Docker` version 1.10 but the latest version is recommended.

Update 2018-03-06: All images/tags on dockerhub has been rebuilt with the latest code and re-uploaded to dockerhub.


# Available tags

The following tags with pre-built images is available on `docker-hub`.

- latest == 4.0.8

Redis 4.0.x versions:

- 4.0.8
- 4.0.7
- 4.0.6
- 4.0.5
- 4.0.4
- 4.0.3
- 4.0.2
- 4.0.1
- 4.0.0

Redis 3.2.x versions:

- 3.2.11
- 3.2.10
- 3.2.9
- 3.2.8
- 3.2.7
- 3.2.6
- 3.2.5
- 3.2.4
- 3.2.3
- 3.2.2
- 3.2.1
- 3.2.0

Redis 3.0.x versions:

- 3.0.7
- 3.0.6
- 3.0.5
- 3.0.4
- 3.0.3
- 3.0.2
- 3.0.1
- 3.0.0



# Usage

There is 2 primary ways of building and running this container


## docker build

To build your own image run:

    docker build -t <username>/redis-cluster .

    # or

    make build

And to run the container use:

    docker run -i -t -p 7000:7000 -p 7001:7001 -p 7002:7002 -p 7003:7003 -p 7004:7004 -p 7005:7005 -p 7006:7006 -p 7007:7007

    # or

    make run

    # and top stop the container run

    make stop

To connect to your cluster you can use the redis-cli tool:

    redis-cli -c -p 7000


## docker compose

It is also possible to build the container using docker-compose. The advantage with a compose build is that it simplifies container management.

To build the container run:

    docker-compose -f docker-compose.yml build

    # or

    make compose-build

To start the container run:

    docker-compose -f docker-compose.yml up

    # or

    make compose-up

    # and to stpo the container run

    make compose-stop

To connection to your cluster you can run redis-cli tool:

    redis-cli -c -p 7000


## Omit standalone redis instances

Set env variable CLUSTER_ONLY=true.

* Running with docker compose, modify docker-compose file

      version: '2'
      services:
        redis-cluster:
          build:
            context: .
            args:
              redis_version: '3.2.7'
          hostname: server
        environment:
          CLUSTER_ONLY: 'true'

* Running with docker directly:

      docker run -i -t -p 7000:7000 -p 7001:7001 -p 7002:7002 -p 7003:7003 -p 7004:7004 -p 7005:7005 -e CLUSTER_ONLY=true <image name/id>


## Build alternative redis versions

For a release to be buildable it needs to be present at this url: http://download.redis.io/releases/


### docker build

To set a different redis version use the argument --build-arg

    # Example
    docker build --build-arg redis_version=3.2.0 -t grokzen/redis-cluster .


### docker-compose

To set a different redis version you must change the variable 'redis_version' inside the docker-compose file.

Then you simply rebuild the compose file and it should be updated with the desired redis version.


# License

This repo is using the MIT LICENSE.

You can find it in the file [LICENSE](LICENSE)
###LIEStudio

This is the read me for the LIEStudio Angular based web application

####Installation
1) Install NodeJs (use NodeJS installer: https://nodejs.org/en/download/)
2) Install the Node Package Manager (npm)
3) Run 'npm install' to install all required packages listed in package.json
4) Run the application in interactive mode using 'gulp serve'

####Fonts
The LIEStudio app uses the Roboto font loaded from the server with fallback to Trebuchet MS, 
Helvetica or Arial respectivly.
Scalable vector icons are from the fontawesome project (http://fontawesome.io).

####UI-elements
The LIEStudio app mostly uses the reusable application UI widgets of the PrimeUI project packed as
Angular 2 directives in the PrimeNG package (https://github.com/primefaces/primeng).
Application specific styles are implemented in the '/shared/styles/_theme_ui_elements.scss'
stylesheet on top of the primeui-ng-all.css.

####Modifications made to third-party packages
- PrimeNG: added support for default active (unfolded) menu panel (PanelMenu component)
  Changed MenuItem interface in 'primeng/components/common.d.ts: added defaultActive?: boolean;
  
  export interface MenuItem {
      label?: string;
      icon?: string;
      command?: (event?: any) => void;
      url?: string;
      routerLink?: any;
      eventEmitter?: EventEmitter<any>;
      items?: MenuItem[];
      defaultActive?: boolean;
  }
  
  Changed PanelMenu.prototype.headerClick function in 'primeng/components/panelmenu/panelmenu.js':
  
  PanelMenu.prototype.headerClick = function (event, item) {
      var index = this.activeItems.indexOf(item);
      if (index == -1)
          this.activeItems.push(item);
      else
          this.activeItems.splice(index, 1);
      if (item.defaultActive) {
        for (var _i = 0, _a = this.model; _i < _a.length; _i++) {
          if (this.model[_i] === item) {
            this.model[_i].defaultActive = false;
          }
        }
      }
      event.preventDefault();
  };
  
  Changed PanelMenu.prototype.isActive function in 'primeng/components/panelmenu/panelmenu.js':
  
  PanelMenu.prototype.isActive = function (item) {
      var index = this.activeItems.indexOf(item);
      if (item.defaultActive && index == -1) {
        this.activeItems.push(item);
        return true;
      }
      return index != -1;
  };

- PrimeNG: Changed menu bar positioning in 'primeng/components/menubar/menubar.js':
 
  line 28 change: sublist.style.left = '0px'; 
          to:     sublist.style.left = '-' + this.domHandler.getOuterWidth(item.children[0]) / 2 + 'px';
  
- PrimeNG: Added class to p-menubarsub class of PrimeNG menubar component in 'primeng/components/menubar/menubar.js':
  
  template: "\n        <ul [ngClass]=\"{'ui-helper-reset':root, 'ui-widget-content ui-corner-all 
  ui-helper-clearfix ui-menu-child ui-shadow':!root}\" class=\"ui-menu-list\"\n            
  (click)=\"listClick($event)\">\n            <template ngFor let-child [ngForOf]=\"(root ? item : item.items)\">\n
  <li #item [ngClass]=\"{'ui-menuitem ui-widget ui-corner-all':true,'ui-menu-parent':child.items,
  'ui-menuitem-active':item==activeItem, 'ui-menuitem-notext':!child.label}\"\n
  (mouseenter)=\"onItemMouseEnter($event, item)\" (mouseleave)=\"onItemMouseLeave($event, item)\">\n
  <a #link [href]=\"child.url||'#'\" class=\"ui-menuitem-link ui-corner-all\" 
  [ngClass]=\"{'ui-state-hover':link==activeLink}\" (click)=\"itemClick($event, child)\">\n
  <span class=\"ui-submenu-icon fa fa-fw\" *ngIf=\"child.items\" [ngClass]=\"{'fa-caret-down':root,
  'fa-caret-right':!root}\"></span>\n <span class=\"ui-menuitem-icon fa fa-fw\" *ngIf=\"child.icon\" 
  [ngClass]=\"child.icon\"></span>\n <span class=\"ui-menuitem-text\">{{child.label}}</span>\n 
   </a>\n <p-menubarSub class=\"ui-submenu\" [item]=\"child\" *ngIf=\"child.items\"></p-menubarSub>\n
  </li>\n            </template>\n        </ul>\n    ",

- PrimeNG: Changed InputSwitch ui px width calculation to match rounded corner switcher in 'primeng/components/inputswitch/inputswitch.js':
  
  Line 81 from: this.onContainer.style.width = this.offset + 'px';
  Line 81 to: this.onContainer.style.width = this.offset + 4 + 'px';
  Line 84 from: this.handle.style.left = this.offset + 'px';
  Line 84 to: this.handle.style.left = this.offset - 4 + 'px';
  
####TODO
- Upon app init, check browser compatibility with respect to enabled javascript, flexbox support.# Manual typings folder

This folder contains all typings needed that have not been published by the community.

The file manual-typings.d.ts include all others typings.###Folder shared

This folder contains all shared informations: services, directives, styles

You can generate through your command prompt a service, a directive or a style with the following commands:

- yo angular2gen:service MyService
- yo angular2gen:directive MyDirective
- yo angular2gen:style MyStyle
###Folder directives

This folder contains all directives for your angular 2 project.

You can generate through your command prompt a directive with the following commands:

```
yo angular2gen:directive MyDirective
```

##### The directive name will be MyDirectiveDirective.
##### For instance, you run yo angular2gen:directive Draggable, the name of the class will be DraggableDirective

As you have seen in the folder architecture of the generator, the folder directives has two folders: one for the sources *src* and another for the tests *test*.
When you run the previous command, it will create two files as follow:
```
- src
         │_ my-directive.directive.ts : The main file of your directive
- test
         │_ my-directive.directive.spec.ts: The test file of your directive
```
###Folder services

This folder contains all services for your angular 2 project.

You can generate through your command prompt a service with the following commands:

```
yo angular2gen:service MyService
```

##### The service name will be MyServiceService.
##### For instance, you run yo angular2gen:service CallDataBase, the name of the class will be CallDataBaseService

As you have seen in the folder architecture of the generator, the folder services has two folders: one for the sources *src* and another for the tests *test*.
When you run the previous command, it will create two files as follow:
```
- src
         │_ my-service.service.ts : The main file of your service
- test
         │_ my-service.service.spec.ts: The test file of your service
```
###Folder styles

This folder contains all global CSS styles for the LIEStudio app as Sass .scss
files.

The main.scss file will serve as the final Sass compiled .css endpoint file 
containing all global app style definitions in one file.

Individual Angular 2 components main implement their own styles or global style
overwrites that operate in the scope of the component. These styles are defined
in .scss files in the component directory.

Sass style variables that determine the look and feel of the app theme are 
conveniatly gathered in the theme_variables.scss file. All substyle defenitions
inherit from this variables file.###Folder assets

This folder contains all resources you need for your project. For instance images, videos, json etc...# Gulp tasks folder# LIEStudio components

Most of the functionality of the LIEStudio application is arranged as a
number of independent Python modules that are installed in the LIEStudio
virtual environment during application setup by the installer.sh script.

The ~/LIEStudio/components directory contains the Python packages for 
these modules that are installed using pip by the LIEStudio installer
script. Any package dependencies will be resolved automatically.

## UnitTesting the components
Each component package contains a `tests` directory with a number of
Python unittests. After activating the virtual environment, they can be 
run as:

    python tests/
  
Some of the component functions require access to a (running) MongoDB
process via PyMongo. If such a process is not active or cannot be 
activated by the lie_db drivers, the unittests for these functions are
scipped.
Components
==========

.. toctree::
   :maxdepth: 4

   components/mdstudio_amber/modules.rst
   components/mdstudio_atb/modules.rst
   components/mdstudio_auth/modules.rst
   components/mdstudio_authorizer/modules.rst
   components/mdstudio_cli/modules.rst
   components/mdstudio_componentbase/modules.rst
   components/mdstudio_config/modules.rst
   components/mdstudio_db/modules.rst
   components/mdstudio_graph/modules.rst
   components/mdstudio_haddock/modules.rst
   components/mdstudio_logger/modules.rst
   components/mdstudio_gromacs/modules.rst
   components/mdstudio_smartcyp/modules.rst
   components/mdstudio_schema/modules.rst
   components/mdstudio_scriptrunner/modules.rst
   components/mdstudio_structures/modules.rst
   components/mdstudio_user/modules.rst
   components/mdstudio_workflow/modules.rst
Welcome to the LIEStudio documentation!
=======================================

Contents:

.. toctree::
   :maxdepth: 3

   about
   development/index
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`.. _about:

=====
About
=====

LIEStudio is an interactive application to explore molecular interactions between small ligands and proteins 
using powerful molecular simulation, data analysis and visualization techniques. The LIEStudio application 
has been constructed with a number of "design goals" in mind that make it:

Flexible
--------

LIEStudio uses an event driven programming design where individual tasks such as docking or simulations are
handeled by dedicated modules that communicate with one another using the WAMP_ messaging protocol (Web Application
Messaging Protocol) implemented using the Crossbar_ WAMP router. 

The functionality of these modules is available at the script level and via a graphical user interface (GUI)
implemented as a web browser application. Performing complicated tasks is now as easy as chaining together the
methods that the individual modules expose in a controled manner.

The LIEStudio API enables easy extention of the applications functionality by addition of new (custom) modules
that extend the functionality. The individual modules can be written in any of the 12 languages_ for wich a WAMP
API is available allowing you to build modules in your language of choice. Most of the default modules shipped 
with the application are written in Python_

.. _Python: http://www.python.org
.. _WAMP: http://wamp-proto.org
.. _Crossbar: http://crossbar.io
.. _languages: http://crossbar.io/about/Supported-Languages/.. _architecture:

Application architecture
========================

LIEStudio uses the Web Application Messaging Protocol (WAMP) to enable messaging between the components of the
application. WAMP communication is handeled by the Autobahn_ familly of WAMP protocol libraries and Crossbar_
operates as WAMP router. The benefits of this approach are plentyfull:

* Autobahn currently provides WAMP libraries for 13 different languages allowing developers to easily contribute
  components to LIEStudio written in their favourite language.
* The WAMP infrastructure enabled by Autobahn and Crossbar allow messaging over TCP and Unix domain sockets as
  transport layer connecting threads, seperated processes on the same machine or processes on physically
  different machines. This flexibility allows for processes on heterogeneous infrastructures such as different
  OS versions, Virtual Machines, Docker images or different servers to be connected together. This means efficient
  scale up- and out.
* The Crossbar router in particular is a mature product that offers additional functionality on top of basic 
  routing such as: communication realms, secure transport (TLS), authorization and authentication among others.
* WAMP bundels Remote Procedure Call (RPC) and Publish & Subscribe messaging patterns in one unified protocol

Directory Structure
-------------------
LIEStudio consists of the following directories:

 * app - The frontend Angularjs_ code.
 * components - The python lie components.
 * data - Configuration files needed to run the application.
 * docker - Docker configuration files.
 * docs - Application documentation folder.
 * scripts - Collection of scripts.



.. _Crossbar: http://crossbar.io
.. _Autobahn: http://autobahn.ws

.. _Angularjs: https://angularjs.org/.. _get-started:

Get Started
===========

.. toctree::
    :maxdepth: 2
    
    docker
    manual

To get started with LIEStudio, we should first install the application. 
First install git_, and after that clone LIEStudio:

.. code-block:: bash

    $ git clone https://github.com/NLeSC/LIEStudio.git


We can now choose two options:

    * A :ref:`docker` installation - The most simple and complete installation (recommended)
    * A complete :ref:`manual` installation - A bit harder but more control

Using LIEStudio
---------------
LIEStudio consists of three parts:

 * Backend code based on Crossbar_ written in python.
 * Frontend code based on Angularjs_ written in typescript.
 * The documentation.

Backend
~~~~~~~
LIEStudio runs both the backend, and will spin up a simple webserver to run the frontend.

On docker we run the application with:

.. code-block:: bash

    $ python .

Manually we have to use pipenv:

.. code-block:: bash

    $ pipenv run python .

.. tip::
   The LIEStudio application runs on  `http://localhost:8080/ <http://localhost:8080>`_ or   `http://localhost/ <http://localhost>`_ on docker.

Frontend
~~~~~~~~
The frontend is compiled using ``gulp compile`` in the ``installer.sh``. When there is no active development
on the GUI, this should suffice. However to manually recompile the frontend, you should change the working
directory to ``app``. After this we have two commands:

 * ``gulp serve`` - A live recompilation of Angularjs_, that allows for simple development.
 * ``gulp compile`` - A single compile to bring LIEStudio to the latest version of the GUI.

 .. tip::
   In the docker container you can also use the command ``serve``, which will run ``gulp serve`` from the correct directory for you.
   Also the command ``compile`` is available.

 .. tip::
   To see the live compilation you should go to `http://localhost:5000 <http://localhost:5000>`_.

Documentation
~~~~~~~~~~~~~
The documentation by default is build when LIEStudio is installed. However we can also update and live serve the documentation while
developing. First go to the ``docs`` directory. After this we can either run ``make html`` or ``make livehtml``.

We can also compile the documentation by using the installer:

.. code-block:: bash

    $ ./installer.sh -d

.. tip::
   In the docker container you can also use the command ``livedocs``, which will run the livedocs reloading.
   Also the command ``docs`` is available.

.. tip::
   To see the live documentation you should go to `http://localhost:8000 <http://localhost:8000>`_.


.. _git: https://git-scm.com/
.. _Crossbar: http://crossbar.io
.. _Angularjs: https://angularjs.org/===========
Development
===========

.. toctree::
    get-started
    architecture
    wamp

In the following chapters we will guide you through the proces of writting modules for the LIEStudio ecosystem.

.. _manual:

Manual Installation
===================
When you want to keep track of the whole environment for yourself, you should follow these
instructions.

Prerequisites
-------------
The only external dependencies are the MongoDB_ NoSQL database
making available the `mongod` MongoDB exacutable in the users path, and the ``nodejs`` package.
For a manual installation, make sure you install:

* MongoDB_ running locally.
* Python_
* nodejs_
* pipenv_

First we need to install npm global dendencies:

* gulp_
* gulp-cli_
* typescript_
* typings_

We install them by running:

.. code-block:: bash

    $ npm i -g gulp gulp-cli typescript typings


Installation
------------
Run the ``installer.sh`` script as:

.. code-block:: bash

    $ ./installer.sh -s

for a quick install using the default Python version. Use -h for more information on
customizing the installation.

A quick install will in sequence:

* Setup a python virtual environment (including installation of a local pip and pipenv_)
* Install required packages from the Python package repository.
* Install LIEStudio component Python packages and their dependencies
* Create a self-signed certificate for WAMP communication over TLS secured websockets.
  Certificate creation requires OpenSSL. If not available the default certificate
  shipped with the package will be used (liestudio/data/crossbar).
  It is recommended to replace the certificate with a personal one signed by a offical
  certificate authority when using the application in a production environment.
* Compile API documentation available from the browser when the program is running at
  http://localhost/help.
  
Usage
-----
The application is started on the command line as:

.. code-block:: bash

    $ pipenv run python .

.. _gulp: http://gulpjs.com/
.. _gulp-cli: https://github.com/gulpjs/gulp-cli
.. _typescript: https://www.typescriptlang.org/
.. _typings: https://github.com/typings/typings


.. _Docker: https://www.docker.com/
.. _MongoDB: https://www.mongodb.com
.. _pipenv: https://github.com/kennethreitz/pipenv_
.. _Python: https://www.python.org/download/releases/2.7/
.. _nodejs: https://nodejs.org/en/.. _docker:

Docker Installation
===================
To make development easier we have setup a docker environment.

Prerequisites
-------------

 * Docker_, that's it.

Usage
-----
To use the docker environment you have to start the container:

.. code-block:: bash

    $ ./start.sh

After this we login to our docker container and we can start LIEStudio with:

.. code-block:: bash

    $ python .

This command will spin up the complete environment including MongoDB, and ssh's into the 
container. When you want to exit this mode just use `>> exit` to exit. Containers can be
stopped using:

.. code-block:: bash

    $ ./stop.sh

SSH Access
----------
When you want to use your private shell you can login using ssh using the following setttings:

+--------+----------------------------------+
| Host   | ``127.0.0.1``                    |
+--------+----------------------------------+
| Key    | ``docker/insecure_id_rsa.ppk``   |
+--------+----------------------------------+
| User   | ``root``                         |
+--------+----------------------------------+
| Port   | ``65432``                        |
+--------+----------------------------------+

IDE Integration
---------------

Pycharms
--------

Go to `File > Project Settings > Project Interpreter`, and add a remote interpreter,
and make sure it matches this screen.

.. image:: ../img/pycharm-config.png

Note specifically:

+--------------------+------------------------------+
| Interpreter path   | ``/app/.venv/bin/python``    |
+--------------------+------------------------------+
| Pycharm helpers    | ``/app/.pycharm_helpers``    |
+--------------------+------------------------------+

Debug Hook
----------
While we now support breakpoints and the likes natively, Pycharm still fails to do post morten
debugging in components. Fixing this is easy; We go to `Run > View Breakpoints`. We add a 
python exception breakpoint. 

.. image:: ../img/pycharm-breakpoint.png

After that we select the runpy._error exception:

.. image:: ../img/pycharm-error.png

Make sure `On Raise` is selected:

.. image:: ../img/pycharm-raise.png.. _wamp:

WAMP messaging
==============

The WAMP based API is the only formal interface of a LIEStudio component and the application as a whole.
The code that enables the functionality of a component is to great extend shielded from the user or other 
components. To ensure that all components together form a functional whole requires the message exchange between
them to be reliable, predictable and complete.

Messages are therefor wrapped in an envelope that defines metadata on the WAMP session, the component dealing with
the request and user information on behalf of which components are contacted. 

* **_realm**:          WAMP realm the component is connected to
* **_package_name**:   Component name
* **_authid**:         user name
* **_class_name**:     name of the WAMP API class of the component
* **_authrole**:       user role
* **_authmethod**:     method used to authenticate the user
* **_session**:        WAMP session ID 
* **_status**:         status message of the task being executed as: submitted, waiting, ready, scheduled, running,
  done, aborted, cancelled or cleared.
* **_init_time**:      time message was first created.
* **_update_time**:    time message was last updated.
* **_status_message**: human readable status message
* **_id**:             unique ID identifing the message. 

**messaging and database storage**
The above information provides a unique identifier of a message from a WAMP communication perspective, a user 
perspective and a component task perspective. As such it provides a complete foundation for storing task
information in a database. 

LIEStudio components store tasks in the database in order to store results, as a persistent store for to follow
the progress of a longer living task and to rerun a task.

WAMP Workflow Schema
====================

All communication is always wrapped in a "task" ensuring that task meta data is send along with every request.
Within the task is the request itself targeting a WAMP method. The WAMP method always accepts the task as only
argument (without optional keyword arguments). The task contains all data required by the method. This way
of working allows the WAMP method to safely extend it's capabilities without breaking the API. Intelligence is
with the WAMP service.

The data required by the WAMP method is always and "input" and a "config" object. Input can be "inline" not having
any dependencies or link to another WAMP method as provider. In the latter case the input itself forms the new
task call to the WAMP method serving as input provider.
A config object holds additional metadata (optionally) required by the WAMP method to perform the action. It is
best practice for the component hosting the WAMP service to define defaults for the configuration. When using the
config API, the default configuration can be overruled by the custom configuration provided. 

The combination of task metadata, input and config data in one object allows the targeted WAMP method to store
an atomic task object in the database using the unique task ID. As the same task ID is also returned with the 
methods (intermediate) results, the sender knows where to find the data or how to ask for it.
This also allows a workflow to be repeated in the same manner as all data and settings are known at every stage
of the workflow.

A task datamodel construct for a PLANTS docking run could look like this:

.. code-block:: python

    {
        '_taskMeta': {
            '_realm':,         
            '_package_name':,  
            '_authid':,        
            '_class_name':,    
            '_authrole':,      
            '_authmethod':,    
            '_session':,       
            '_status':,        
            '_init_time':,     
            '_update_time':,   
            '_status_message':,
            '_id':            
        },
        '_dataType': 'wampMethod',
        '_wampUrl': 'liestudio.docking.run',
        '_inputDict': {
                'ligandFile': {
                '_dataType': 'wampMethod',
                '_wampUrl': 'liestudio.structure.get',
                '_configDict': {
                    'structure_id': '100203',
                    'format': 'mol2'
                },
            },
            'proteinFile': {
                '_dataType': 'inlineSource',
                '_data: '<structure inline pdb>'
            }
        },
        '_outputDict': {
            '<results>'
        },
        '_configDict': {
            'method', 'plants',
            'workdir': '/home/workdir',
            'bindingsite_center': [0,0,0]
        }
    }

Message type
------------

LIEStudio tasks communicate data in a number of predefined types indicated by the ``_dataType`` tag:

* **wampMethod**: the data to be send or retrieved is handeled by a different WAMP method and thus can
  be regarded as a new task. A data construct if type wampMethod is required to have the
  ``_wampUrl`` tag with the fully qualified WAMP URL of the method to be called and either a
  ``_inputDict`` or ``_configDict`` depending on the method specifications. Other tags are optional.
* **inlineSource**: the data is send inline using the ``_data`` tag. 
* **fileSource**: the data is located in a file at the ``_url`` location with optional tags ``_fileType``.

``_inputDict``:
WAMP methods may accept an arbitrary number of input values. These can be regarded as the arguments 
of a Python function. Keyword arguments are stored in the ``_configDict``.
Each input is of a certain data type as described above.

``_outputDict``:
Equal to the ``_inputDict`` in capabilities.

``_configDict``:
Keyword values to the method.lie\_componentbase\.tests package
=================================

.. automodule:: lie_componentbase.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_componentbase.tests.module

Submodules
----------

lie\_componentbase\.tests\.wamp\_api\_test module
-------------------------------------------------

.. automodule:: lie_componentbase.tests.wamp_api_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_componentbase\.tests\.module package
=========================================

.. automodule:: lie_componentbase.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_componentbase\.tests\.module\.system\_test module
------------------------------------------------------

.. automodule:: lie_componentbase.tests.module.system_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_componentbase\.config package
==================================

.. automodule:: lie_componentbase.config
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_componentbase\.config\.config\_format module
-------------------------------------------------

.. automodule:: lie_componentbase.config.config_format
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.config\.config\_handler module
--------------------------------------------------

.. automodule:: lie_componentbase.config.config_handler
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.config\.config\_io module
---------------------------------------------

.. automodule:: lie_componentbase.config.config_io
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.config\.config\_orm\_handler module
-------------------------------------------------------

.. automodule:: lie_componentbase.config.config_orm_handler
    :members:
    :undoc-members:
    :show-inheritance:


lie_componentbase
=================

.. toctree::
   :maxdepth: 4

lie\_componentbase package
==========================

.. automodule:: lie_componentbase
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_componentbase.config
    lie_componentbase.tests

Submodules
----------

lie\_componentbase\.application\_session module
-----------------------------------------------

.. automodule:: lie_componentbase.application_session
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.component\_manager module
---------------------------------------------

.. automodule:: lie_componentbase.component_manager
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.logger module
---------------------------------

.. automodule:: lie_componentbase.logger
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.runner module
---------------------------------

.. automodule:: lie_componentbase.runner
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.util module
-------------------------------

.. automodule:: lie_componentbase.util
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.wamp\_schema module
---------------------------------------

.. automodule:: lie_componentbase.wamp_schema
    :members:
    :undoc-members:
    :show-inheritance:

lie\_componentbase\.wamp\_taskmeta module
-----------------------------------------

.. automodule:: lie_componentbase.wamp_taskmeta
    :members:
    :undoc-members:
    :show-inheritance:


lie\_auth\.tests\.wamp package
==============================

.. automodule:: mdstudio_auth.tests.wamp
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_auth\.tests\.wamp\.api\_test module
----------------------------------------

.. automodule:: mdstudio_auth.tests.wamp.api_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_auth\.tests\.module package
================================

.. automodule:: mdstudio_auth.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_auth\.tests\.module\.user\_test module
-------------------------------------------

.. automodule:: mdstudio_auth.tests.module.user_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_auth\.oauth package
========================

.. automodule:: mdstudio_auth.oauth
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_auth\.oauth\.client module
-------------------------------

.. automodule:: mdstudio_auth.oauth.client
    :members:
    :undoc-members:
    :show-inheritance:

lie\_auth\.oauth\.request\_validator module
-------------------------------------------

.. automodule:: mdstudio_auth.oauth.request_validator
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_auth
========

.. toctree::
   :maxdepth: 4

   mdstudio_auth
lie\_auth\.tests package
========================

.. automodule:: mdstudio_auth.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_auth.tests.module
    mdstudio_auth.tests.wamp

lie\_auth package
=================

.. automodule:: mdstudio_auth
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_auth.oauth
    mdstudio_auth.tests

Submodules
----------

lie\_auth\.authorizer module
----------------------------

.. automodule:: mdstudio_auth.authorizer
    :members:
    :undoc-members:
    :show-inheritance:

lie\_auth\.password\_retrieval module
-------------------------------------

.. automodule:: mdstudio_auth.password_retrieval
    :members:
    :undoc-members:
    :show-inheritance:

lie\_auth\.sendmail module
--------------------------

.. automodule:: mdstudio_auth.sendmail
    :members:
    :undoc-members:
    :show-inheritance:

lie\_auth\.util module
----------------------

.. automodule:: mdstudio_auth.util
    :members:
    :undoc-members:
    :show-inheritance:

lie\_auth\.wamp\_services module
--------------------------------

.. automodule:: mdstudio_auth.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie_scriptrunner
================

.. toctree::
   :maxdepth: 4

   lie_scriptrunner
lie\_scriptrunner package
=========================

.. automodule:: lie_scriptrunner
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_scriptrunner\.wamp\_services module
----------------------------------------

.. automodule:: lie_scriptrunner.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie\_md package
===============

.. automodule:: mdstudio_gromacs
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_md\.\_\_main\_\_ module
----------------------------

.. automodule:: mdstudio_gromacs.__main__
    :members:
    :undoc-members:
    :show-inheritance:

lie\_md\.wamp\_services module
------------------------------

.. automodule:: mdstudio_gromacs.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_gromacs
======

.. toctree::
   :maxdepth: 4

   mdstudio_gromacs
lie\_plants\_docking\.tests\.wamp package
=========================================

.. automodule:: lie_plants_docking.tests.wamp
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_plants\_docking\.tests\.wamp\.api\_test module
---------------------------------------------------

.. automodule:: lie_plants_docking.tests.wamp.api_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_plants\_docking\.tests\.module package
===========================================

.. automodule:: lie_plants_docking.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_plants\_docking\.tests\.module\.docking\_test module
---------------------------------------------------------

.. automodule:: lie_plants_docking.tests.module.docking_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_plants\_docking package
============================

.. automodule:: lie_plants_docking
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_plants_docking.tests

Submodules
----------

lie\_plants\_docking\.clustering module
---------------------------------------

.. automodule:: lie_plants_docking.clustering
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.docking\_base module
------------------------------------------

.. automodule:: lie_plants_docking.docking_base
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.docking\_settings module
----------------------------------------------

.. automodule:: lie_plants_docking.docking_settings
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.plants\_conf module
-----------------------------------------

.. automodule:: lie_plants_docking.plants_conf
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.plants\_docking module
--------------------------------------------

.. automodule:: lie_plants_docking.plants_docking
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.utils module
----------------------------------

.. automodule:: lie_plants_docking.utils
    :members:
    :undoc-members:
    :show-inheritance:

lie\_plants\_docking\.wamp\_services module
-------------------------------------------

.. automodule:: lie_plants_docking.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie_plants_docking
==================

.. toctree::
   :maxdepth: 4

   lie_plants_docking
lie\_plants\_docking\.tests package
===================================

.. automodule:: lie_plants_docking.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_plants_docking.tests.module
    lie_plants_docking.tests.wamp

lie_authorizer
==============

.. toctree::
   :maxdepth: 4

lie\_authorizer package
=======================

.. automodule:: lie_authorizer
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_authorizer\.wamp\_services module
--------------------------------------

.. automodule:: lie_authorizer.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio\_structures package
=======================

.. automodule:: mdstudio_structures
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

mdstudio\_structures\.wamp\_services module
--------------------------------------

.. automodule:: mdstudio_structures.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_structures
==============

.. toctree::
   :maxdepth: 4

   mdstudio_structures
lie\_haddock package
====================

.. automodule:: mdstudio_haddock
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_haddock\.HADDOCKCmdLineApp module
--------------------------------------

.. automodule:: mdstudio_haddock.HADDOCKCmdLineApp
    :members:
    :undoc-members:
    :show-inheritance:

lie\_haddock\.HADDOCKInterface module
-------------------------------------

.. automodule:: mdstudio_haddock.HADDOCKInterface
    :members:
    :undoc-members:
    :show-inheritance:

lie\_haddock\.run\_haddock module
---------------------------------

.. automodule:: mdstudio_haddock.run_haddock
    :members:
    :undoc-members:
    :show-inheritance:

lie\_haddock\.wamp\_services module
-----------------------------------

.. automodule:: mdstudio_haddock.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_haddock
===========

.. toctree::
   :maxdepth: 4

   mdstudio_haddock
lie\_schema package
===================

.. automodule:: lie_schema
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_schema\.wamp\_services module
----------------------------------

.. automodule:: lie_schema.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie_schema
==========

.. toctree::
   :maxdepth: 4

   lie_schema
lie\_logger\.tests package
==========================

.. automodule:: lie_logger.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_logger.tests.module

lie_logger
==========

.. toctree::
   :maxdepth: 4

   lie_logger
lie\_logger\.tests\.module package
==================================

.. automodule:: lie_logger.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_logger\.tests\.module\.logger\_test module
-----------------------------------------------

.. automodule:: lie_logger.tests.module.logger_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_logger package
===================

.. automodule:: lie_logger
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_logger\.\_\_main\_\_ module
--------------------------------

.. automodule:: lie_logger.__main__
    :members:
    :undoc-members:
    :show-inheritance:

lie\_logger\.wamp\_services module
----------------------------------

.. automodule:: lie_logger.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_amber
=========

.. toctree::
   :maxdepth: 4

   mdstudio_amber
lie\_amber package
==================

.. automodule:: mdstudio_amber
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_amber\.ambertools module
-----------------------------

.. automodule:: mdstudio_amber.ambertools
    :members:
    :undoc-members:
    :show-inheritance:

lie\_amber\.cli\_runner module
------------------------------

.. automodule:: mdstudio_amber.cli_runner
    :members:
    :undoc-members:
    :show-inheritance:

lie\_amber\.settings module
---------------------------

.. automodule:: mdstudio_amber.settings
    :members:
    :undoc-members:
    :show-inheritance:

lie\_amber\.wamp\_services module
---------------------------------

.. automodule:: mdstudio_amber.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie\_config\.tests\.module package
==================================

.. automodule:: lie_config.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_config\.tests\.module\.config\_test module
-----------------------------------------------

.. automodule:: lie_config.tests.module.config_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_config\.tests\.module\.decorator\_class module
---------------------------------------------------

.. automodule:: lie_config.tests.module.decorator_class
    :members:
    :undoc-members:
    :show-inheritance:

lie\_config\.tests\.module\.decorator\_test module
--------------------------------------------------

.. automodule:: lie_config.tests.module.decorator_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_config\.tests\.module\.io\_test module
-------------------------------------------

.. automodule:: lie_config.tests.module.io_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_config\.tests\.module\.orm\_test module
--------------------------------------------

.. automodule:: lie_config.tests.module.orm_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_config package
===================

.. automodule:: lie_config
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_config.tests

Submodules
----------

lie\_config\.\_\_main\_\_ module
--------------------------------

.. automodule:: lie_config.__main__
    :members:
    :undoc-members:
    :show-inheritance:

lie\_config\.wamp\_services module
----------------------------------

.. automodule:: lie_config.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie\_config\.tests\.wamp package
================================

.. automodule:: lie_config.tests.wamp
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_config\.tests\.wamp\.api module
------------------------------------

.. automodule:: lie_config.tests.wamp.api
    :members:
    :undoc-members:
    :show-inheritance:


lie\_config\.tests package
==========================

.. automodule:: lie_config.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_config.tests.module
    lie_config.tests.wamp

lie_config
==========

.. toctree::
   :maxdepth: 4

   lie_config
lie\_atb package
================

.. automodule:: mdstudio_atb
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_atb\.\_\_main\_\_ module
-----------------------------

.. automodule:: mdstudio_atb.__main__
    :members:
    :undoc-members:
    :show-inheritance:

lie\_atb\.atb\_api\_py2 module
------------------------------

.. automodule:: mdstudio_atb.atb_api_py2
    :members:
    :undoc-members:
    :show-inheritance:

lie\_atb\.settings module
-------------------------

.. automodule:: mdstudio_atb.settings
    :members:
    :undoc-members:
    :show-inheritance:

lie\_atb\.wamp\_services module
-------------------------------

.. automodule:: mdstudio_atb.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_atb
=======

.. toctree::
   :maxdepth: 4

   mdstudio_atb
lie\_user\.tests\.module package
================================

.. automodule:: mdstudio_user.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_user\.tests\.module\.user\_test module
-------------------------------------------

.. automodule:: mdstudio_user.tests.module.user_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_user package
=================

.. automodule:: mdstudio_user
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_user.tests

Submodules
----------

lie\_user\.password\_retrieval module
-------------------------------------

.. automodule:: mdstudio_user.password_retrieval
    :members:
    :undoc-members:
    :show-inheritance:

lie\_user\.sendmail module
--------------------------

.. automodule:: mdstudio_user.sendmail
    :members:
    :undoc-members:
    :show-inheritance:

lie\_user\.util module
----------------------

.. automodule:: mdstudio_user.util
    :members:
    :undoc-members:
    :show-inheritance:

lie\_user\.wamp\_services module
--------------------------------

.. automodule:: mdstudio_user.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


lie\_user\.tests\.wamp package
==============================

.. automodule:: mdstudio_user.tests.wamp
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_user\.tests\.wamp\.api\_test module
----------------------------------------

.. automodule:: mdstudio_user.tests.wamp.api_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_user\.tests package
========================

.. automodule:: mdstudio_user.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_user.tests.module
    mdstudio_user.tests.wamp

mdstudio_user
========

.. toctree::
   :maxdepth: 4

lie\_graph\.tests\.module package
=================================

.. automodule:: lie_graph.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_graph\.tests\.module\.graph\_test module
---------------------------------------------

.. automodule:: lie_graph.tests.module.graph_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.tests\.module\.graphalgorithms\_test module
-------------------------------------------------------

.. automodule:: lie_graph.tests.module.graphalgorithms_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.tests\.module\.graphaxis\_test module
-------------------------------------------------

.. automodule:: lie_graph.tests.module.graphaxis_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.tests\.module\.graphdict\_test module
-------------------------------------------------

.. automodule:: lie_graph.tests.module.graphdict_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.tests\.module\.graphmathop\_test module
---------------------------------------------------

.. automodule:: lie_graph.tests.module.graphmathop_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.tests\.module\.graphorm\_test module
------------------------------------------------

.. automodule:: lie_graph.tests.module.graphorm_test
    :members:
    :undoc-members:
    :show-inheritance:


lie\_graph\.tests package
=========================

.. automodule:: lie_graph.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_graph.tests.module

lie_graph
=========

.. toctree::
   :maxdepth: 4

   lie_graph
lie\_graph\.io package
======================

.. automodule:: lie_graph.io
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_graph\.io\.io\_dict\_parser module
---------------------------------------

.. automodule:: lie_graph.io.io_dict_parser
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.io\.io\_helpers module
----------------------------------

.. automodule:: lie_graph.io.io_helpers
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.io\.io\_json\_format module
---------------------------------------

.. automodule:: lie_graph.io.io_json_format
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.io\.io\_pgf\_format module
--------------------------------------

.. automodule:: lie_graph.io.io_pgf_format
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.io\.io\_tgf\_format module
--------------------------------------

.. automodule:: lie_graph.io.io_tgf_format
    :members:
    :undoc-members:
    :show-inheritance:


lie\_graph package
==================

.. automodule:: lie_graph
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_graph.io
    lie_graph.tests

Submodules
----------

lie\_graph\.config\_graph module
--------------------------------

.. automodule:: lie_graph.config_graph
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph module
------------------------

.. automodule:: lie_graph.graph
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_algorithms module
------------------------------------

.. automodule:: lie_graph.graph_algorithms
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_axis\_class module
-------------------------------------

.. automodule:: lie_graph.graph_axis_class
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_axis\_methods module
---------------------------------------

.. automodule:: lie_graph.graph_axis_methods
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_dict module
------------------------------

.. automodule:: lie_graph.graph_dict
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_helpers module
---------------------------------

.. automodule:: lie_graph.graph_helpers
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_io module
----------------------------

.. automodule:: lie_graph.graph_io
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_math\_operations module
------------------------------------------

.. automodule:: lie_graph.graph_math_operations
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_mixin module
-------------------------------

.. automodule:: lie_graph.graph_mixin
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_orm module
-----------------------------

.. automodule:: lie_graph.graph_orm
    :members:
    :undoc-members:
    :show-inheritance:

lie\_graph\.graph\_query module
-------------------------------

.. automodule:: lie_graph.graph_query
    :members:
    :undoc-members:
    :show-inheritance:


lie\_db\.tests\.module package
==============================

.. automodule:: mdstudio_db.tests.module
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_db\.tests\.module\.db\_test module
---------------------------------------

.. automodule:: mdstudio_db.tests.module.db_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_db\.tests\.module\.test\_db\_methods module
------------------------------------------------

.. automodule:: mdstudio_db.tests.module.test_db_methods
    :members:
    :undoc-members:
    :show-inheritance:


lie\_db package
===============

.. automodule:: mdstudio_db
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_db.tests

Submodules
----------

lie\_db\.cache\_dict module
---------------------------

.. automodule:: mdstudio_db.cache_dict
    :members:
    :undoc-members:
    :show-inheritance:

lie\_db\.db\_methods module
---------------------------

.. automodule:: mdstudio_db.db_methods
    :members:
    :undoc-members:
    :show-inheritance:

lie\_db\.wamp\_services module
------------------------------

.. automodule:: mdstudio_db.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_db
======

.. toctree::
   :maxdepth: 4

   mdstudio_db
lie\_db\.tests package
======================

.. automodule:: mdstudio_db.tests
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    mdstudio_db.tests.module

lie\_workflow package
=====================

.. automodule:: mdstudio_workflow
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    lie_workflow.tests

Submodules
----------

lie\_workflow\.common module
----------------------------

.. automodule:: lie_workflow.common
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.task\_metadata module
------------------------------------

.. automodule:: lie_workflow.task_metadata
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.task\_specs module
---------------------------------

.. automodule:: lie_workflow.task_specs
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.wamp\_services module
------------------------------------

.. automodule:: lie_workflow.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.workflow module
------------------------------

.. automodule:: lie_workflow.workflow
    :members:
    :undoc-members:
    :show-inheritance:


lie\_workflow\.tests package
============================

.. automodule:: lie_workflow.tests
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_workflow\.tests\.dummy\_task\_runners module
-------------------------------------------------

.. automodule:: lie_workflow.tests.dummy_task_runners
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.tests\.module\_branched\_test module
---------------------------------------------------

.. automodule:: lie_workflow.tests.module_branched_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.tests\.module\_breakpoint\_test module
-----------------------------------------------------

.. automodule:: lie_workflow.tests.module_breakpoint_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.tests\.module\_io\_test module
---------------------------------------------

.. automodule:: lie_workflow.tests.module_io_test
    :members:
    :undoc-members:
    :show-inheritance:

lie\_workflow\.tests\.module\_linear\_test module
-------------------------------------------------

.. automodule:: lie_workflow.tests.module_linear_test
    :members:
    :undoc-members:
    :show-inheritance:


lie_workflow
============

.. toctree::
   :maxdepth: 4

   lie_workflow
lie\_cli package
================

.. automodule:: mdstudio_cli
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

lie\_cli\.\_\_main\_\_ module
-----------------------------

.. automodule:: mdstudio_cli.__main__
    :members:
    :undoc-members:
    :show-inheritance:

lie\_cli\.wamp\_services module
-------------------------------

.. automodule:: mdstudio_cli.wamp_services
    :members:
    :undoc-members:
    :show-inheritance:


mdstudio_cli
=======

.. toctree::
   :maxdepth: 4

   mdstudio_cli
