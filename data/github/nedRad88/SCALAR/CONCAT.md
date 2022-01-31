# SCALAR - Streaming ChALlenge plAtfoRm
![Logo](./logo.PNG)


SCALAR is the first of its kind, platform to organize Machine Learning competitions on Data Streams.
It was used to organize a real-time ML competition on IEEE Big Data Cup Challenges 2019.

## Features

**Data stream mining competitions** - With SCALAR you can organize the real-time competitions on Data Streams. 
It is inspired by [Kaggle](https://www.kaggle.com/), platform for offline machine learning competitions.



**Simple, user friendly interface** - SCALAR has a WEB application that allows you to easily browse and 
subscribe to the competitions. Its simple and intuitive design, let's you to easily upload the datasets, and create the 
competition. 

![Dialog competition](./add_competition.png)

**Live results and leaderboard** - Since the competition is real-time, the results are also updated in real time. 
During the competition, you can follow performance of your model in the WEB application. Leaderboard will show how 
do you compare to other users, and live chart shows the comparison with baseline model.

**Secure, bi-directional streaming communication** - We use a combination of `gRPC` and `Protobuf` to provide secure, 
low latency bi-directional streaming communication between server and users.

![Architecture](./Architecture.png)

**Freedom to choose a programming language** - SCALAR lets users to choose their preferred environment. The only 
requirement is to be able to communicate through `gRPC` and `Protobuf`, which is supported for many programming 
languages: Python, Java, C++, GO... Additionally, SCALAR provides support for R. Apart from that, users can choose 
their setup, environment and additional resources to train better models.

**** 

## Getting Started

The project is done in Python and organized in Docker containers.
Each service is a separate Docker container.

### Prerequisites

To run the platform locally, Docker is needed:

[Install Docker](https://docs.docker.com/get-docker/)

Also, Docker compose should be installed:

[Install Docker compose](https://docs.docker.com/compose/install/)

### Running

Running is done using Docker-compose.

### Quick setup (to test the functionality of the platform)

- Run ```setup.py``` and follow the instructions to setup the environment. The script will set up the 
time zone and create the docker network for all containers.

 - Once the ```setup.py``` finished successfully, the platform can be run by:
```
docker-compose up
```

### In depth setup (If you want to use it for organizing competitions and enable the registration):

 - Download the code locally and then adjust the [config.json](provider/my_application/config.json) and [docker-compose.yml](./docker-compose.yml) files. More details in [config-ReadMe.txt](provider/config-ReadMe.txt) and in [docker-compose-ReadMe.txt](./docker-compose-ReadMe.txt).

 - Set up an email account which will be used to send the registration confirmation message and authentication token.
For that, you will need to set up your email account to allow the access of less secure apps.
For a quick start, update only email information in [config.json](provider/my_application/config.json).

 - In [docker-compose.yml](./docker-compose.yml) update only the local paths to mount a persistent volumes, following the [docker-compose-ReadMe.txt](./docker-compose-ReadMe.txt).

 - Run [setup.py](./setup.py) and follow the instructions to setup the environment. The script will set up the 
time zone and create the docker network for all containers.
```python3 setup.py```

 - Once the ```setup.py``` finished successfully, the platform can be run by:
```
docker-compose up
```

This command will pull all necessary containers and run them.
When all services are up, web application will be available on [localhost:80](http://localhost:80)

To log in to the platform, you can use default credentials: `admin:admin`

For the test purposes, the test competition will be created automatically.

It will be scheduled to start ```5 minutes``` after starting the platform.

Once you log in, you will be able to see the competition under ```Competitions/Coming``` tab.

In order to subscribe to the competition, click on it and then click the ```Subscribe``` button.
![Subscribe to competition](./subscribe_competition.png)

### Setting up the client

Navigate to [example_data](./example_data) directory, and run:

```python3 client_setup.py```

Once the necessary packages have been installed, go to [Python client directory](./example_data/user_code/Code/Python)
and edit the [client.py](./example_data/user_code/Code/Python/client.py) file. Copy the ```Competition code``` and
```Secret key``` from a competition page and add it in ```client.py``` as shown in the figure below:

![Client setup](./client_setup.png)

Once the competition has started, run ```client.py```, and you should be able to see how the messages and predictions are exchanged.
Then you will be able to see the live chart and leaderboard on the competition page. (You will have to refresh the page to get new measures.)

To get to know around the platform use the the [Quick Start Guide](./SCALAR_Quick_Start_Guide.pdf). 
To create and participate in the competition use the provided [examples](./example_data).
* **Nedeljko Radulovic**
* **Dihia Boulegane**
* **Albert Bifet**
* **Nenad Stojanovic**


## Acknowledgments

Open source Docker containers were used:
* [MongoDB](https://hub.docker.com/_/mongo)
* [Spark Docker container by GettyImages](https://hub.docker.com/r/gettyimages/spark)
* [MySQL](https://hub.docker.com/_/mysql)
* [Kafka by Wurstmeister](https://hub.docker.com/r/wurstmeister/kafka)
* [Zookeeper by Wurstmeister](https://hub.docker.com/r/wurstmeister/zookeeper)

## References

* [Boulegane, Dihia, et al. "Real-Time Machine Learning Competition on Data Streams at the IEEE Big Data 2019." 2019 IEEE International Conference on Big Data (Big Data). IEEE, 2019.](https://ieeexplore.ieee.org/abstract/document/9006357?casa_token=f0mJeR8-WfYAAAAA:yEt_Mix9dumrPpo64uPBbI0XI4Kvfim4Pkg5xNVVzXqK4AGToX0XcJPKgETkE1hs86Pcc0u5xYc)
* [IEEE Big Data Cup 2019](https://bigmine.github.io/real-time-ML-competition/index.html)


<!-- 
Thank you for contributing with a PR!

Please fill the description of change(s) and/or indicate if it fixes an open issue (optional).

To ease the merge process please review the attached checklist.
-->

Changes proposed in this pull request:

* 
* 
* Fixes # [Optional]

---

### Checklist
#### Implementation
- [ ] Implementation is correct (it performs its intended function).
- [ ] Code is consistent with the framework.
- [ ] Code is properly documented.
- [ ] PR description covers ALL the changes performed.
- [ ] Files changed (update, add, delete) are in the PR's scope (no extra files are included).

#### Tests
- [ ] New functionality is tested.
- [ ] ALL tests pass with no errors.
### CONTRIBUTING

We are happy to have you reading this file as we welcome contributions from community.
Here you will find the information on how to start contributing to `SCALAR`.

When contributing to this repository, please first discuss the change you wish to make via issue, 
email, or any other method with the owners of this repository before making a change.

## Contribuiton to a GitHub Project

If you are not familiar with GitHub and how to contribute to a github project, check the following links as a good first step:

 - [The Beginner's guide to contributing to a GitHub project](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/)

 - [Make a clean pull request](https://github.com/MarcDiethelm/contributing/blob/master/README.md)
 
If you have an idea what can be improved or added to the system, or you found a bug, report an [issue](https://github.com/nedRad88/SCALAR/issues).
We will respond promptly to your proposal.

## Code style and documentation

 - Python code shall comply with [PEP 8](https://www.python.org/dev/peps/pep-0008/) 
 - Documentation shall be in docstring.
 - The names for new modules and methods should be intuitive.
 - If possible, provide a way to test the new module.

## Submitting updates

Send your GitHub pull request to [SCALAR](https://github.com/nedRad88/SCALAR/pulls) with a clear statement of what you have done.
Ensure the PR description clearly describes the problem and solution. Include the relevant issue number if applicable.
Test the functionality of the module before sending the PR. 
Please follow our code conventions and write clear messages for your commits.



---
title: '`SCALAR` - A Platform for Real-time Machine Learning Competitions on Data Streams'
bibliography: references.bib
tags:
  - Python
  - Stream Data Mining
  - Real-time Machine Learning
  - Information Flow Processing
authors:
  - name: Nedeljko Radulovic
    affiliation: 1
  - name: Dihia Boulegane
    affiliation: "1, 2"
  - name: Albert Bifet
    affiliation: "1, 3"
affiliations:
 - name: LTCI, Télécom Paris, IP-Paris, Paris, France
   index: 1
 - name: Orange Labs, Grenoble, France
   index: 2
 - name: University of Waikato, New Zealand
   index: 3
date: 26 August 2020
nocite: '@*'
---

# Summary

`SCALAR` is a new platform for running real-time machine learning competitions on data streams. Following the intent of Kaggle, which serves as a platform for organizing machine learning competitions adapted for batch learning, we propose `SCALAR` as a novel platform explicitly designed for stream learning in real-time. `SCALAR` supports both classification and regression problems in the data streaming setting. It has been developed in `Python`, using state of the art open-source solutions: `Apache Kafka`, `Apache Spark`, `gRPC`, `Protobuf`, and `Docker`. 

# Statement of need 

Existing machine learning platforms for research and competitions, such as [`Kaggle`](https://www.kaggle.com/) and [`Alibaba Tianchi`](https://tianchi.aliyun.com/) provide static data sets that users leverage to do batch training and scoring. These popular platforms, especially `Kaggle`, enable users to propose better solutions for a wide range of machine learning problems and applications, pushing forward the whole research community.

While these platforms enable collaborative model building, they do not provide real-time settings and data streaming instead of static data sets. Inspired by the great reception and adoption that these platforms have had, `SCALAR` was built for real-time environments by supporting data streaming to fill this gap, helping users test and improve their methodologies in a real-time machine learning scenario. 

On `SCALAR`, data are continuously released in batches during defined time intervals. Predictions for each current batch sent before designated deadline are evaluated in real-time, and the results are shown on a live leaderboard.

`SCALAR` can be applied to provide real-time machine learning support in multiple scenarios, one of the more notable and recent ones was its role as the core organizing and coordination platform for the first Real-time Machine Learning Competition on Data Streams [@boulegane2019real] as part of the [IEEE Big Data 2019 Cup Challenges](http://bigdataieee.org/BigData2019/BigDataCupChallenges.html).


# Streaming learning setting


Most of the existing platforms for data science competitions are tailored to offline learning where a static dataset is made available to the participants before the competition starts. This dataset is divided into training and test sets. The training set is used to build and train the model, which is then scored on the test set. 

In online learning, data arrive in a stream of records (instances) and the model needs to be trained incrementally as new instances are released. Since the data arrive at high speed, predictions have to be issued within a short time. Considering this specific setup of the data stream mining scenario (Figure \ref{fig:stream_mining}), every model should use a limited 
amount of time and memory, process one instance at a time and inspect it only once. The model must be able to predict at any time [@bifet2010moa].

![Stream data mining scenario\label{fig:stream_mining}](stream_mining.png)

The model is continuosly updated using the labeled instances and then evaluated on new unlabeled instances. This scenario represents the prequential evaluation scenario or "test-then-train" setting, to make it possible we first assume that the records in the data stream arrive in batches and that each batch can be of size `1` or more. Then, we define two different kinds of batches:

* Initial batch: This batch is used to build and train the initial model. It is aimed to avoid the cold start problem and as such, contains both features and targets. The initial batch usually has a large number of instances.

* Regular batch: The regular batch contains new test instances while providing training instances of previous batches that are used to evaluate the quality of the model up to the current time.

![Initial and regular batches in the data stream\label{fig:online_learning}](online_learning.jpg)

# Architecture

To support all the aforementioned requirements for stream data mining and to be able to organize competitions in such a scenario, a specifically dedicated platform was needed. To the best of our knowledge, there doesn't exist such a platform. Thus we propose `SCALAR`. Figure \ref{architecture} highlights the architecture of the platform that includes a web application, and a streaming engine following the fundamentals of Information Flow Processing (IFP)[@Cugola:12].

![Architecture of the platform\label{fig:architecture}](Architecture.png)

Layers and components' descriptions:

* Data sources: `SCALAR` enables competitions by providing data in the form of a `CSV` file, which is then used to recreate a continuous stream following predefined competition settings such as the time interval between two releases, and the number of instances at a time.
    
* Data layer: represents the persistence node in the system where different kinds of information are stored. `MongoDB` is used for unstructured data (e.g., user’s predictions, records from the stream, evaluation metrics, etc.) and `MySQL` for structured data (competition information, user information).
    
* Streaming Engine: this layer is responsible for handling the streams. From the `CSV` files, `Kafka` recreates multiple streams to be sent to multiple users. `SCALAR` provides a secure and independent bi-directional streaming communication between the user and the streaming engine. This allows each user to consume training and test instances and submit the respective predictions according to the competition predefined data format. The resource server is the `API` that handles all authenticated requests from the client application, whereas the authorization server is in charge of users' identification. The Online evaluation engine handles both the stream of instances, and the streams of participants' predictions in order to compute and update the evaluation metrics online before storing them into the dedicated database.
    
* Application layer: consists of two parts the web application, and client application. The web application represents a user-friendly interface that allows participants to register, subscribe to competitions, and follow leaderboard and model performance online. The client application is provided to accommodate participants' code to solve the machine learning problem at hand. It delivers a secure communication channel to receive the stream of training and test instances and send their respective predictions to the server.

# Acknowledgements

We would like to thank Nenad Stojanovic for his contribution to the project and to acknowledge support for this project 
from the EDF (Electricité de France) and OpenML.


# References

Getting started with SCALAR
***************************

# SCALAR - Streaming ChALlenge plAtfoRm
![Logo](./_static/images/logo.png)


SCALAR is the first of its kind, platform to organize Machine Learning competitions on Data Streams.
It was used to organize a real-time ML competition on IEEE Big Data Cup Challenges 2019.

## Features

**Data stream mining competitions** - With SCALAR you can organize the real-time competitions on Data Streams. 
It is inspired by [Kaggle](https://www.kaggle.com/), platform for offline machine learning competitions.



**Simple, user friendly interface** - SCALAR has a WEB application that allows you to easily browse and 
subscribe to the competitions. Its simple and intuitive design, let's you to easily upload the datasets, and create the 
competition. 

![Dialog competition](./add_competition.png)

**Live results and leaderboard** - Since the competition is real-time, the results are also updated in real time. 
During the competition, you can follow performance of your model in the WEB application. Leaderboard will show how 
do you compare to other users, and live chart shows the comparison with baseline model.

**Secure, bi-directional streaming communication** - We use a combination of `gRPC` and `Protobuf` to provide secure, 
low latency bi-directional streaming communication between server and users.

![Architecture](./Architecture.png)

**Freedom to choose a programming language** - SCALAR lets users to choose their preferred environment. The only 
requirement is to be able to communicate through `gRPC` and `Protobuf`, which is supported for many programming 
languages: Python, Java, C++, GO... Additionally, SCALAR provides support for R. Apart from that, users can choose 
their setup, environment and additional resources to train better models.

**** 

## Getting Started

The project is done in Python and organized in Docker containers.
Each service is a separate Docker container.

### Prerequisites

To run the platform locally, Docker is needed:

[Install Docker](https://docs.docker.com/get-docker/)

Also, Docker compose should be installed:

[Install Docker compose](https://docs.docker.com/compose/install/)

### Running

Running is done using Docker-compose.

 - Download the code locally and then adjust the [config.json](provider/my_application/config.json) and [docker-compose.yml](./docker-compose.yml) files. More details in [config-ReadMe.txt](provider/config-ReadMe.txt) and in [docker-compose-ReadMe.txt](./docker-compose-ReadMe.txt).

 - Set up an email account which will be used to send the registration confirmation message and authentication token.
For that, you will need to set up your email account to allow the access of less secure apps.
For a quick start, update only email information in [config.json](provider/my_application/config.json).

 - In [docker-compose.yml](./docker-compose.yml) update only the local paths to mount a persistent volumes, following the [docker-compose-ReadMe.txt](./docker-compose-ReadMe.txt).

 - To run the SCALAR application using Docker-compose, first create the Docker bridge network on your local machine:
```
docker network create --driver bridge provider_network --subnet=172.22.0.0/16 --ip-range=172.22.0.0/24 --gateway=172.22.0.1

```
You can choose the IP ranges according to your preferences.

 - Once the [config.json](provider/my_application/config.json) and [docker-compose.yml](./docker-compose.yml) have been set up and Docker network has been created,
  the platform can be run by:

```
docker-compose up
```

This command will pull all necessary containers and run them.
When all services are up, web application will be available on [localhost:80](http://localhost:80)

To log in to the platform, you can use default credentials: `admin:admin`
Default `Test` datastream, for creating the test competition is available under `Datastreams` tab.


To get to know around the platform use the the [Quick Start Guide](./SCALAR_Quick_Start_Guide.pdf). 
To create and participate in the competition use the provided [examples](./example_data).
* **Nedeljko Radulovic**
* **Dihia Boulegane**
* **Albert Bifet**
* **Nenad Stojanovic**


## Acknowledgments

Open source Docker containers were used:
* [MongoDB](https://hub.docker.com/_/mongo)
* [Spark Docker container by GettyImages](https://hub.docker.com/r/gettyimages/spark)
* [MySQL](https://hub.docker.com/_/mysql)
* [Kafka by Wurstmeister](https://hub.docker.com/r/wurstmeister/kafka)
* [Zookeeper by Wurstmeister](https://hub.docker.com/r/wurstmeister/zookeeper)

## References

* [Boulegane, Dihia, et al. "Real-Time Machine Learning Competition on Data Streams at the IEEE Big Data 2019." 2019 IEEE International Conference on Big Data (Big Data). IEEE, 2019.](https://ieeexplore.ieee.org/abstract/document/9006357?casa_token=f0mJeR8-WfYAAAAA:yEt_Mix9dumrPpo64uPBbI0XI4Kvfim4Pkg5xNVVzXqK4AGToX0XcJPKgETkE1hs86Pcc0u5xYc)
* [IEEE Big Data Cup 2019](https://bigmine.github.io/real-time-ML-competition/index.html)


angular-random-string
=====================

Angular service that generates a random alphanumeric string certain length.

## install

````bash
	$ bower install --save angular-random-string 
````

##usage

1. Load script

```html
	<script src="/path/to/angular-random-string.js"></script>
```

2. Add to your module

```javascript
	var app = angular.module('MyApp', ['angularRandomString']);
```

3. Inject to controller/service/etc.

```javascript
	app.controller('myController', function(randomString){
		var str = randomString(); // -> random alphanumeric string 10 char length
		var str32 = randomString(32); // -> random alphanumeric string 32 char length
	});
```<a name="0.2.5"></a>

# 0.2.5 (2015-02-23)
- build v0.2.5
- bug: revert UMD support due to breaking changes  (#288, #289, #290)
- bug: fix extend (PR #286)
- chore: fix typos in CHANGELOG

<a name="0.2.4"></a>

# 0.2.4 (2015-02-18)
- build v0.2.4
- Fixed jshint isuses
- added UMD support #273
- fixed broken tests
- updated bower and npm dependencies
- added .editorconfig file
- updated LICENSE #268

<a name="0.2.3"></a>

# 0.2.3 (2015-10-11)
- build v0.2.3
- Fixed jshint issues
- Updated mixed-up dates in change log

<a name="0.2.2"></a>

# 0.2.2 (2015-05-29)
- build v0.2.2
- fix(localStorage): parsing safety #230

<a name="0.2.1"></a>

# 0.2.1 (2015-05-18)

## Breaking Changes
- build v0.2.1
- refac(common): remove deprecated code
- fix(localStorage): load/save objects #225
- Fix for bug introduced in 0.2.0 boolean values that not in objects are not converted properly

<a name="0.2.0"></a>

# 0.2.0 (2015-05-10)

## Breaking Changes
- build v0.2.0
- fromJson was replaced by JSON.parse with reviver fn
- Allow multiple keys for `.remove()`
- Fixed wrong angular dependence version.
- docs(README): Video Tutorial
- Update Documentation
- Set individual expiry for cookies
- docs(README.md): get started
- docs(README.md): gitter badge
- Added Gitter badge
- refa(src): remove duplicated stringification of value while storing
- style(src): indentation
- fixed issue noted in issue #193
- Changes to clearAll function - accept RegExp
- add --save argument to install with bower
- docs(README.md): cookie.clearAll hash/title
- docs(README.md): expand cookie.clearAll
- Update README.md
- fix(LICENSE): update copyright year
- fix(README.md): add \n just for aesthetic
- docs(demo): better example and clearAll
- Update README.md
- fix(README.md): add badges
- Improved documentation for setStorageCookieDomain.
- Fix a minor typo in a comment
- docs(REAMME.md): version

<a name="0.1.1"></a>
# 0.1.1 (2014-10-06)


## Breaking Changes
- update your `index.html` file to reference angular-local-storage at its new
  path inside the `dist` directory `/angular-local-storage/dist/angular-local-storage.js`
angular-local-storage
=====================
An Angular module that gives you access to the browsers local storage

[![NPM version][npm-image]][npm-url]
[![Build status][travis-image]][travis-url]
[![Test coverage][coveralls-image]][coveralls-url]
[![Dependency Status][david-image]][david-url]
[![License][license-image]][license-url]
[![Downloads][downloads-image]][downloads-url]

##Table of contents:
- [![Gitter][gitter-image]][gitter-url]
- [Get Started](#get-started)
- [Video Tutorial](https://www.youtube.com/watch?v=I4iB0kOSmx8)
- [Development](#development)
- [Configuration](#configuration)
 - [setPrefix](#setprefix)
 - [setStorageType](#setstoragetype)
 - [setStorageCookie](#setstoragecookie)
 - [setStorageCookieDomain](#setstoragecookiedomain)
 - [setNotify](#setnotify)
 - [Example](#configuration-example)
- [API Documentation](#api-documentation)
 - [isSupported](#issupported)
 - [getStorageType](#getstoragetype)
 - [set](#set)
 - [get](#get)
 - [keys](#keys)
 - [remove](#remove)
 - [clearAll](#clearall)
 - [bind](#bind)
 - [deriveKey](#derivekey)
 - [length](#length)
 - [cookie](#cookie)
    - [isSupported](#cookieissupported)
    - [set](#cookieset)
    - [get](#cookieget)
    - [remove](#cookieremove)
    - [clearAll](#cookieclearall)

##Get Started
**(1)** You can install angular-local-storage using 2 different ways:<br/>
**Git:**
clone & build [this](https://github.com/grevory/angular-local-storage.git) repository<br/>
**Bower:**
```bash
$ bower install angular-local-storage --save
```
**npm:**
```bash
$ npm install angular-local-storage
```
**(2)** Include `angular-local-storage.js` (or `angular-local-storage.min.js`) from the [dist](https://github.com/grevory/angular-local-storage/tree/master/dist) directory in your `index.html`, after including Angular itself.

**(3)** Add `'LocalStorageModule'` to your main module's list of dependencies.

When you're done, your setup should look similar to the following:

```html
<!doctype html>
<html ng-app="myApp">
<head>

</head>
<body>
    ...
    <script src="//ajax.googleapis.com/ajax/libs/angularjs/1.1.5/angular.min.js"></script>
    <script src="bower_components/js/angular-local-storage.min.js"></script>
    ...
    <script>
        var myApp = angular.module('myApp', ['LocalStorageModule']);

    </script>
    ...
</body>
</html>
```
##Configuration
###setPrefix
You could set a prefix to avoid overwriting any local storage variables from the rest of your app<br/>
**Default prefix:** `ls.<your-key>`
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setPrefix('yourAppName');
});
```
###setStorageType
You could change web storage type to localStorage or sessionStorage<br/>
**Default storage:** `localStorage`
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setStorageType('sessionStorage');
});
```
###setStorageCookie
Set cookie options (usually in case of fallback)<br/>
**expiry:** number of days before cookies expire (0 = does not expire). **default:** `30`<br/>
**path:** the web path the cookie represents. **default:** `'/'`
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setStorageCookie(45, '<path>');
});
```
###setStorageCookieDomain
Set the cookie domain, since this runs inside a the `config()` block, only providers and constants can be injected.  As a result, `$location` service can't be used here, use a hardcoded string or `window.location`.<br/>
**No default value**
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setStorageCookieDomain('<domain>');
});
```

For local testing (when you are testing on localhost) set the domain to an empty string ''. Setting the domain to 'localhost' will not work on all browsers (eg. Chrome) since some browsers only allow you to set domain cookies for registry controlled domains, i.e. something ending in .com or so, but not IPs **or intranet hostnames** like localhost. </br>

###setNotify
Send signals for each of the following actions:<br/>
**setItem** , default: `true`<br/>
**removeItem** , default: `false`
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setNotify(true, true);
});
```
###Configuration Example
Using all together
```js
myApp.config(function (localStorageServiceProvider) {
  localStorageServiceProvider
    .setPrefix('myApp')
    .setStorageType('sessionStorage')
    .setNotify(true, true)
});
```
##API Documentation
##isSupported
Checks if the browser support the current storage type(e.g: `localStorage`, `sessionStorage`).  
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  if(localStorageService.isSupported) {
    //...
  }
  //...
});
```
###getStorageType
**Returns:** `String`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  var storageType = localStorageService.getStorageType(); //e.g localStorage
  //...
});
```
###set
Directly adds a value to local storage.<br/>
If local storage is not supported, use cookies instead.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function submit(key, val) {
   return localStorageService.set(key, val);
  }
  //...
});
```
###get
Directly get a value from local storage.<br/>
If local storage is not supported, use cookies instead.<br/>
**Returns:** `value from local storage`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function getItem(key) {
   return localStorageService.get(key);
  }
  //...
});
```
###keys
Return array of keys for local storage, ignore keys that not owned.<br/>
**Returns:** `value from local storage`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  var lsKeys = localStorageService.keys();
  //...
});
```
###remove
Remove an item(s) from local storage by key.<br/>
If local storage is not supported, use cookies instead.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function removeItem(key) {
   return localStorageService.remove(key);
  }
  //...
  function removeItems(key1, key2, key3) {
   return localStorageService.remove(key1, key2, key3);
  }
});
```
###clearAll
Remove all data for this app from local storage.<br/>
If local storage is not supported, use cookies instead.<br/>
**Note:** Optionally takes a regular expression string and removes matching.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function clearNumbers(key) {
   return localStorageService.clearAll(/^\d+$/);
  }
  //...
  function clearAll() {
   return localStorageService.clearAll();
  }
});
```
###bind
Bind $scope key to localStorageService.  
**Usage:** `localStorageService.bind(scope, property, value[optional], key[optional])`  
***key:*** The corresponding key used in local storage  
**Returns:** deregistration function for this listener.
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  localStorageService.set('property', 'oldValue');
  $scope.unbind = localStorageService.bind($scope, 'property');

  //Test Changes
  $scope.update = function(val) {
    $scope.property = val;
    $timeout(function() {
      alert("localStorage value: " + localStorageService.get('property'));
    });
  }
  //...
});
```
```html
<div ng-controller="MainCtrl">
  <p>{{property}}</p>
  <input type="text" ng-model="lsValue"/>
  <button ng-click="update(lsValue)">update</button>
  <button ng-click="unbind()">unbind</button>
</div>
```

###deriveKey
Return the derive key  
**Returns** `String`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  localStorageService.set('property', 'oldValue');
  //Test Result
  console.log(localStorageService.deriveKey('property')); // ls.property
  //...
});
```
###length
Return localStorageService.length, ignore keys that not owned.  
**Returns** `Number`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  var lsLength = localStorageService.length(); // e.g: 7
  //...
});
```
##Cookie
Deal with browser's cookies directly.
##cookie.isSupported
Checks if cookies are enabled in the browser.  
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  if(localStorageService.cookie.isSupported) {
    //...
  }
  //...
});
```
###cookie.set
Directly adds a value to cookies.<br/>
**Note:** Typically used as a fallback if local storage is not supported.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function submit(key, val) {
   return localStorageService.cookie.set(key, val);
  }
  //...
});
```
**Cookie Expiry** Pass a third argument to specify number of days to expiry
```js
    localStorageService.cookie.set(key,val,10)
```
sets a cookie that expires in 10 days.
###cookie.get
Directly get a value from a cookie.<br/>
**Returns:** `value from local storage`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function getItem(key) {
   return localStorageService.cookie.get(key);
  }
  //...
});
```
###cookie.remove
Remove directly value from a cookie.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function removeItem(key) {
   return localStorageService.cookie.remove(key);
  }
  //...
});
```
###cookie.clearAll
Remove all data for this app from cookie.<br/>
**Returns:** `Boolean`
```js
myApp.controller('MainCtrl', function($scope, localStorageService) {
  //...
  function clearAll() {
   return localStorageService.cookie.clearAll();
  }
});
```

Check out the full demo at http://gregpike.net/demos/angular-local-storage/demo.html

##Development:
* Don't forget about tests.
* If you planning add some feature please create issue before.

Clone the project:
```sh
$ git clone https://github.com/<your-repo>/angular-local-storage.git
$ npm install
$ bower install
```
Run the tests:
```sh
$ grunt test
```
**Deploy:**<br/>
Run the build task, update version before(bower,package)
```sh
$ grunt dist
$ git tag 0.*.*
$ git push origin master --tags
```

[npm-image]: https://img.shields.io/npm/v/angular-local-storage.svg?style=flat-square
[npm-url]: https://npmjs.org/package/angular-local-storage
[travis-image]: https://img.shields.io/travis/grevory/angular-local-storage.svg?style=flat-square
[travis-url]: https://travis-ci.org/grevory/angular-local-storage
[coveralls-image]: https://img.shields.io/coveralls/grevory/angular-local-storage.svg?style=flat-square
[coveralls-url]: https://coveralls.io/r/grevory/angular-local-storage
[david-image]: http://img.shields.io/david/grevory/angular-local-storage.svg?style=flat-square
[david-url]: https://david-dm.org/grevory/angular-local-storage
[license-image]: http://img.shields.io/npm/l/angular-local-storage.svg?style=flat-square
[license-url]: LICENSE
[downloads-image]: http://img.shields.io/npm/dm/angular-local-storage.svg?style=flat-square
[downloads-url]: https://npmjs.org/package/angular-local-storage
[gitter-image]: https://badges.gitter.im/Join%20Chat.svg
[gitter-url]: https://gitter.im/grevory/angular-local-storage?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
# angular-socket-io
Bower Component for using AngularJS with [Socket.IO](http://socket.io/),
based on [this](http://briantford.com/blog/angular-socket-io.html).


## Install

1. `bower install angular-socket-io` or [download the zip](https://github.com/btford/angular-socket-io/archive/master.zip).
2. Made sure the Socket.IO client lib is loaded. It's often served at `/socket.io/socket.io.js`.
3. Include the `socket.js` script provided by this component into your app.
4. Add `btford.socket-io` as a module dependency to your app.


## Usage

This module exposes a `socketFactory`, which is an API for instantiating
sockets that are integrated with Angular's digest cycle.


### Making a Socket Instance

```
// in the top-level module of the app
angular.module('myApp', [
  'btford.socket-io',
  'myApp.MyCtrl'
]).
factory('mySocket', function (socketFactory)) {
  return socketFactory();
});
```

With that, you can inject your `mySocket` service into controllers and
other serivices within your application!

### Using Your Socket Instance

Building on the example above:

```
// in the top-level module of the app
angular.module('myApp', [
  'btford.socket-io',
  'myApp.MyCtrl'
]).
factory('mySocket', function (socketFactory)) {
  return socketFactory();
});
```


## API

For the most part, this component works exactly like you would expect.
The only API addition is `socket.forward`, which makes it easier to add/remove listeners in a way that works with [AngularJS's scope](http://docs.angularjs.org/api/ng.$rootScope.Scope).

### `socket.on` / `socket.addListener`
Takes an event name and callback.
Works just like the method of the same name from Socket.IO.

### `socket.removeListener`
Takes an event name and callback.
Works just like the method of the same name from Socket.IO.

### `socket.emit`
Sends a message to the server.
Optionally takes a callback.

Works just like the method of the same name from Socket.IO.

### `socket.forward`

`socket.forward` allows you to forward the events received by Socket.IO's socket to AngularJS's event system.
You can then listen to the event with `$scope.$on`.
By default, socket-forwarded events are namespaced with `socket:`.

The first argument is a string or array of strings listing the event names to be forwarded.
The second argument is optional, and is the scope on which the events are to be broadcast.
If an argument is not provided, it defaults to `$rootScope`.
As a reminder, broadcasted events are propagated down to descendant scopes.

#### Examples

An easy way to make socket error events available across your app:

```javascript
// in the top-level module of the app
angular.module('myApp', [
  'btford.socket-io',
  'myApp.MyCtrl'
]).
factory('mySocket', function (socketFactory)) {
  var mySocket = socketFactory();
  mySocket.forward('error');
  return mySocket;
});

// in one of your controllers
angular.module('myApp.MyCtrl', []).
  controller('MyCtrl', function ($scope) {
    $scope.$on('socket:error', function (ev, data) {

    });
  });
```

Avoid duplicating event handlers when a user navigates back and forth between routes:

```javascript
angular.module('myMod', ['btford.socket-io']).
  controller('MyCtrl', function ($scope, socket) {
    socket.forward('someEvent', $scope);
    scope.$on('socket:someEvent', function (ev, data) {
      $scope.theData = data;
    });
  });
```


### `socketFactory({ ioSocket: }}`

This option allows you to provide the `socket` service with a `Socket.IO socket` object to be used internally.
This is useful if you want to connect on a different path, or need to hold a reference to the `Socket.IO socket` object for use elsewhere.

```javascript
angular.module('myApp', [
  'btford.socket-io'
]).
factory('mySocket', function (socketFactory)) {
  var myIoSocket = io.connect('/some/path');

  mySocket = socketFactory({
    ioSocket: myIoSocket
  });

  return mySocket;
});
```

### `socketFactory({ scope: })`

This option allows you to set the scope on which `$broadcast` is forwarded to when using the `forward` method.
It defaults to `$rootScope`.


### `socketFactory({ prefix: })`

The default prefix is `socket:`.


#### Example

To remove the prefix:

```javascript
angular.module('myApp', [
  'btford.socket-io'
]).
config(function (socketProvider) {
  socketProvider.prefix('');
});
```



## Migrating from 0.2 to 0.3

`angular-socket-io` version `0.3` changes X to make fewer assumptions
about the lifecycle of the socket. Previously, the assumption was that your
application has a single socket created at config time. While this holds
for most apps I've seen, there's no reason you shouldn't be able to
lazily create sockets, or have multiple connections.

In `0.2`, `angular-socket-io` exposed a `socket` service. In `0.3`, it
instead exposes a `socketFactory` service which returns socket instances.
Thus, getting the old API is as simple as making your own `socket` service
with `socketFactory`. The examples below demonstrate how to do this.

### Simple Example

In most cases, adding the following to your app should suffice:

```javascript
// ...
factory('socket', function (socketFactory)) {
  return socketFactory();
});
// ...
```

### Example with Configuration

Before:

```javascript
angular.module('myApp', [
  'btford.socket-io'
]).
config(function (socketProvider)) {
  socketProvider.prefix('foo~');
  socketProvider.ioSocket(io.connect('/some/path'));
}).
controller('MyCtrl', function (socket) {
  socket.on('foo~bar', function () {
    $scope.bar = true;
  });
});
```

After:

```javascript
angular.module('myApp', [
  'btford.socket-io'
]).
factory('socket', function (socketFactory)) {
  return socketFactory({
    prefix: 'foo~',
    ioSocket: io.connect('/some/path')
  });
}).
controller('MyCtrl', function (socket) {
  socket.on('foo~bar', function () {
    $scope.bar = true;
  });
});
```


## License
MIT
my\_application.repositories package
====================================

Submodules
----------

my\_application.repositories.BaselineToMongo module
---------------------------------------------------

.. automodule:: my_application.repositories.BaselineToMongo
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.repositories.CompetitionRepository module
---------------------------------------------------------

.. automodule:: my_application.repositories.CompetitionRepository
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.repositories.KafkaToMongo module
------------------------------------------------

.. automodule:: my_application.repositories.KafkaToMongo
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: my_application.repositories
   :members:
   :undoc-members:
   :show-inheritance:
API Reference
**************

my\_application package
=======================

Subpackages
-----------

.. toctree::
   :maxdepth: 2

   my_application.repositories

Submodules
----------

my\_application.auth module
---------------------------

.. automodule:: my_application.auth
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.consumer module
-------------------------------

.. automodule:: my_application.consumer
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.producer module
-------------------------------

.. automodule:: my_application.producer
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.repository module
---------------------------------

.. automodule:: my_application.repository
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.sparkEvaluator module
-------------------------------------

.. automodule:: my_application.sparkEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.sparkToMongo module
-----------------------------------

.. automodule:: my_application.sparkToMongo
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.stream\_server module
-------------------------------------

.. automodule:: my_application.stream_server
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.subscription\_auth module
-----------------------------------------

.. automodule:: my_application.subscription_auth
   :members:
   :undoc-members:
   :show-inheritance:

my\_application.views module
----------------------------

.. automodule:: my_application.views
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: my_application
   :members:
   :undoc-members:
   :show-inheritance:
.. SCALAR documentation master file, created by
   sphinx-quickstart on Wed Oct 14 18:39:55 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SCALAR's documentation!
==================================

.. image:: ./_static/images/logo.png
   :width: 600px
   :alt: scalar-logo

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Getting started <get_started.md>
   API Reference <my_application.rst>
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
