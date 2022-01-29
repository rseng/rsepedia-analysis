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