Change Log
==========

This document records all notable changes to Boatswain.

This project adheres to `Semantic Versioning <http://semver.org/>`_.

`Unreleased`
------------
* Removed python 2 and 3.4 tests
* Python 2 support will no longer be actively maintained

`1.0.4`_
--------
* Updated dynamic string message to variable

`1.0.3`_
--------
* Added a non-zero exit code when a build fails
* Added a build summary show succesful and non-succesfully built images
* Added a keep-building (-k) command line argument
* Building now stops at the first failed image by default

`1.0.2`_
--------
* Correct citation.cff file

`1.0.1`_
--------
* Fixed some packaging things

`1.0.0`_
--------

* Fixed help text of push command
* Fixed extraction of image id from docker response in some cases
* Windows compatibility
* Appveyor windows tests are passing


`0.7.0`_ (2017-04-03)
---------------------

* Added a 'before' and 'command' key to the build definition.
  This is a list of commands that need to be staged into the context directory.
* Default verbosity only shows 1 progress bar for all images
* Changed progress indication to full white block

`0.6.0`_ (2017-03-09)
---------------------

* Added the tree command which will print the tree of the boatswain file
* Added quiet and extra verbose modes

`0.5.1`_ (2017-02-10)

* Fixed issue with printing unicode text from the docker stream

`0.5.0`_ (2017-02-10)
---------------------

* Implemented push command
* Build will now greedily try to build images instead of throwing an
  exception at the first error.
* Added error messages to failing builds
* Standardized return values (e.g. always a list)
* Refactored to reduce code duplication in boatswain class

`0.4.0`_ (2017-02-09)
---------------------

* Progress timer now increases every second
* Improved error reporting (No longer uses an exception)

`0.3.0`_ (2017-02-08)
---------------------

* Added a whole bunch of tests
* Added the clean command
* Changed file layout from recursive to using from

`0.2.0`_ (2017-02-06)
---------------------

* Added dry-run option
* Added ability to build only one image

`0.1.0`_ (2017-02-02)
---------------------

* Initial release


.. _0.1.0: https://github.com/nlesc-sherlock/boatswain/commit/f8b85edd3ed9f21c04fa846eae1af7abed8d0d77
.. _0.2.0: https://github.com/nlesc-sherlock/boatswain/compare/f8b85ed...0.2.0
.. _0.3.0: https://github.com/nlesc-sherlock/boatswain/compare/0.2.0...0.3.0
.. _0.4.0: https://github.com/nlesc-sherlock/boatswain/compare/0.3.0...0.2.0
.. _0.5.0: https://github.com/nlesc-sherlock/boatswain/compare/0.4.0...0.5.0
.. _0.5.1: https://github.com/nlesc-sherlock/boatswain/compare/0.5.0...0.5.1
.. _0.6.0: https://github.com/nlesc-sherlock/boatswain/compare/0.5.1...0.6.0
.. _0.7.0: https://github.com/nlesc-sherlock/boatswain/compare/0.6.0...0.7.0
.. _1.0.0: https://github.com/nlesc-sherlock/boatswain/compare/0.7.0...1.0.0
.. _1.0.1: https://github.com/nlesc-sherlock/boatswain/compare/1.0.0...1.0.1
.. _1.0.2: https://github.com/nlesc-sherlock/boatswain/compare/1.0.1...1.0.2
.. _1.0.3: https://github.com/nlesc-sherlock/boatswain/compare/1.0.2...1.0.3
.. _1.0.4: https://github.com/nlesc-sherlock/boatswain/compare/1.0.3...1.0.4
.. image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/nlesc/boatswain
.. image:: https://img.shields.io/github/last-commit/NLeSC/boatswain
   :target: https://github.com/nlesc/boatswain
.. image:: https://img.shields.io/github/license/NLeSC/boatswain
   :target: https://github.com/nlesc/boatswain
.. image:: https://img.shields.io/badge/rsd-boatswain-00a3e3.svg?labelColor=gray&color=00a3e3
   :target: https://research-software.nl/software/boatswain
.. image:: https://img.shields.io/pypi/v/boatswain
   :target: https://pypi.org/project/boatswain/
.. image:: https://zenodo.org/badge/80722427.svg
   :target: https://zenodo.org/badge/latestdoi/80722427
.. image:: https://bestpractices.coreinfrastructure.org/projects/3758/badge
    :target: https://bestpractices.coreinfrastructure.org/projects/3758
.. image:: https://travis-ci.org/NLeSC/boatswain.svg?branch=master
    :target: https://travis-ci.org/NLeSC/boatswain
.. image:: https://ci.appveyor.com/api/projects/status/5n7uj8ownch05e34/branch/master?svg=true
    :target: https://ci.appveyor.com/project/NLeSC/boatswain/branch/master
.. image:: https://api.codacy.com/project/badge/Grade/3271fedeb94b4d8ba24273d9cba0c852
    :target: https://www.codacy.com/app/NLeSC/boatswain?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/boatswain&amp;utm_campaign=Badge_Grade

Boatswain
=========
Boatswain is a simple build system for docker images.

It is especially useful when you have multiple docker images that
depend on each other.


Installation
============

Boatswain is a simple python script you can install with pip

::

    $ pip install boatswain


Usage
=====
Create a file called boatswain.yml for your project with the following
syntax, which is heavily based on docker-compose.

.. code-block:: yaml

    version: 1.0                    # The version of the boatswain yaml
    organisation: boatswain         # Your dockerhub organisation
    images:
        image1:pytest:              # the key will be used to tag the image
            context: docker/image1  # The path of the dockerfile
        image2:pytest:
            context: docker/image2
            from: image1:pytest     # This image depends on the other image
        image3:pytest:
            context: docker/image3
            from: image2:pytest
        image4:pytest:
            context: docker/image4
            tag: image12:pytest     # This image will be tagged with this

Building
--------

You can build the images defined in the boatswain file using the following
command. This invokes the docker build process for each image in the order
so the dependencies are built before the dependent image.

::

    $ boatswain build

Cleaning
--------

You can clean the images that are defined in the boatswain file.

::

    $ boatswain clean

Pushing
-------

You can push the built images to dockerhub using the `push` command.
It will be pushed to `organisation/imagetag` on dockerhub.

::

    $ boatswain push

Extra Options
=============
-h
    Display the options

-q, --quiet
    Don't display any output

-k, --keep_building
    Keep building images even if an error occurs

--dryrun
    Do not actually execute the command, just go through the motions

-b <boatswain_file>, --boatswain_file <boatswain_file>
    Override the default boatswain file (boatswain.yml)

-f, --force
    Force building the images, even if they already exist
    (only for build)

Debugging your build
====================
When your build does not go the way you expected boatswain
can display some more information by using:

-v
    Verbose mode, displays a build progress for each image

-vv
    Very verbose mode, displays the output of the docker build process

--debug
    Debug mode, displays debug information of boatswain
    as well as the output of the docker build process
.. image:: https://travis-ci.org/nlesc-sherlock/boatswain.svg?branch=master
    :target: https://travis-ci.org/nlesc-sherlock/boatswain


Boatswain
=========
Boatswain is a simple build system for docker images.

It is especially usefull when you have multiple docker images that
depend on each other.


Installation
============

Boatswain is a simple python script you can install with pip

::

    $ pip install boatswain


Usage
=====
Create a file called boatswain.yml for your project with the following
syntax, which is heavily based on docker-compose.

.. code-block:: yaml

    version: 1.0                    # The version of the boatswain yaml
    organisation: boatswain         # Your dockerhub organisation
    images:
        image1:pytest:              # the key will be used to tag the image
            context: docker/image1  # The path of the dockerfile
        image2:pytest:
            context: docker/image2
            from: image1:pytest     # This image depends on the other image
        image3:pytest:
            context: docker/image3
            from: image2:pytest
        image4:pytest:
            context: docker/image4
            tag: image12:pytest     # This image will be tagged with this
