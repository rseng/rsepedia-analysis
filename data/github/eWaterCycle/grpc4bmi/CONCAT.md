# grpc4bmi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1462641.svg)](https://doi.org/10.5281/zenodo.1462641)
[![CI](https://github.com/eWaterCycle/grpc4bmi/workflows/CI/badge.svg)](https://github.com/eWaterCycle/grpc4bmi/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/grpc4bmi/badge/?version=latest)](https://grpc4bmi.readthedocs.io/en/latest/?badge=latest)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=grpc4bmi&metric=alert_status)](https://sonarcloud.io/dashboard?id=grpc4bmi)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=grpc4bmi&metric=coverage)](https://sonarcloud.io/dashboard?id=grpc4bmi)

## Purpose

This software allows you to wrap your Basic Model Interface (BMI) implementation ([https://github.com/csdms/bmi](https://github.com/csdms/bmi)) in a server process and communicate with it via the included Python client. The communication is serialized to protocol buffers by gRPC ([https://grpc.io/](https://grpc.io/)) and occurs over network ports.

## Installation

Optionally, create your virtual environment and activate it, Then, run

```bash
pip install grpc4bmi
```

on the client (Python) side. If your server model is implemented in Python, do the same in the server environment (e.g. docker container). If the model is implemented in R, run instead

```bash
pip install grpc4bmi[R]
```

in the server environment. For bleeding edge version from GitHub use

```bash
pip install git+https://github.com/eWaterCycle/grpc4bmi.git#egg=grpc4bmi
```

Finally if the model is implemented in C or C++, clone this git repo and run

```bash
make
make install
```

in the cpp folder.

## Usage

### Model written in Python

For inspiration look at the example in the test directory. To start a server process that allows calls to your BMI implementation, type

```bash
run-bmi-server --name <PACKAGE>.<MODULE>.<CLASS> --port <PORT> --path <PATH>
```

where ```<PACKAGE>, <MODULE>``` are the Python package and module containing your implementation, ```<CLASS>``` is your
bmi model class name, ```<PORT>``` is any available port on the host system, and optionally ```<PATH>``` denotes an
additional path that should be added to the system path to make your implementation work. The name option above is
optional, and if not provided the script will look at the environment variables ```BMI_PACKAGE```, ```BMI_MODULE``` and
```BMI_CLASS```. Similarly, the port can be defined by the environment variable ```BMI_PORT```.
This software assumes that your implementation constructor has no parameters.

### Model written in C/C++ (beta)

Create an executable along the lines of cpp/run-bmi-server.cc. You can copy the file and replace the function

```C++
Bmi* create_model_instance()
{
    /* Return your new BMI instance pointer here... */
}
```

with the instantiation of your model BMI. The model needs to implement the csdms BMI for C, but you may also implement our more object-oriented C++ interface [BmiCppExtension](https://github.com/eWaterCycle/grpc4bmi/blob/master/cpp/bmi_cpp_extension.h).

### Model written in R

The grpc4bmi Python package can also run BMI models written in R if the model is a subclass of [AbstractBmi](https://github.com/eWaterCycle/bmi-r/blob/master/R/abstract-bmi.R#L9)
See [https://github.com/eWaterCycle/bmi-r](https://github.com/eWaterCycle/bmi-r) for instruction on R and Docker.

Run the R model a server with

```bash
run-bmi-server --lang R [--path <R file with BMI model>] --name [<PACKAGE>::]<CLASS> --port <PORT>
```

For example with [WALRUS](https://github.com/eWaterCycle/grpc4bmi-examples/tree/master/walrus) use

```bash
run-bmi-server --lang R --path ~/git/eWaterCycle/grpc4bmi-examples/walrus/walrus-bmi.r --name WalrusBmi --port 55555
```

### The client side

The client side has only a Python implementation. The default BMI client assumes a running server process on a given port.

```python
from grpc4bmi.bmi_grpc_client import BmiClient
import grpc
mymodel = BmiClient(grpc.insecure_channel("localhost:<PORT>"))
print mymodel.get_component_name()
mymodel.initialize(<FILEPATH>)
...further BMI calls...
```

The package contains also client implementation that own the server process, either as a Python subprocess or a docker
image or a singularity image running the ```run-bmi-server``` script. For instance

```python
from grpc4bmi.bmi_client_subproc import BmiClientSubProcess
mymodel = BmiClientSubProcess(<PACKAGE>.<MODULE>.<CLASS>)
```

will automatically launch the server in a sub-process and

```python
from grpc4bmi.bmi_client_subproc import BmiClientDocker
mymodel = BmiClientDocker(<IMAGE>,<PORT>)
```

will launch a docker container, assuming that a gRPC BMI server will start and exposes the port ```<PORT>```.

```python
from grpc4bmi.bmi_client_singularity import BmiClientSingularity
mymodel = BmiClientSingularity(<IMAGE>,<PORT>)
```

will launch a singularity container, assuming that a gRPC BMI server will start and exposes the port ```<PORT>```.

For more documentation see [https://grpc4bmi.readthedocs.io/](https://grpc4bmi.readthedocs.io/).

## Development: generating the gRPC code

When developers change the proto-file, it is necessary to install gRPC tools Python packages in your Python environment:

```bash
pip install -r requirements.txt
pip install -e .
# For R integration also install the R extras with
pip install -e .[R]
```

and install the C++ runtime and `protoc` command as described in <https://github.com/google/protobuf/blob/master/src/README.md>.
After this, simply executing the `proto_gen.sh` script should do the job.

## Future work

More language bindings are underway.
.. _pip-install:

Installing requirements
=======================

To use `grpc4bmi`_ the Python package should be installed.

If your model uses some virtual environment with installed dependencies (e.g. Anaconda or virtualenv), activate this environment before installing grpc4bmi.

Install the grpc4bmi python package with pip:

.. code-block:: sh

    $ pip install grpc4bmi

This will install the latest release of the package; for the most recent github revision type instead

.. code-block:: sh

    $ pip install git+https://github.com/eWaterCycle/grpc4bmi.git#egg=grpc4bmi

.. _grpc4bmi: https://pypi.org/project/grpc4bmi/
.. _usage:

Using the client
================

We assume that service is always dedicated to a single client, addressing a BMI model with multiple users at the same time results in undefined behavior.

.. _python-grpc4bmi-client:

Python BMI Client
.................

For a given running BMI service process connected to networking port ``<PORT>``, we can start communicating with this server by instantiating the :class:`grpc4bmi.bmi_grpc_client.BmiClient` python class:

.. code-block:: python

    import grpc
    from grpc4bmi.bmi_grpc_client import BmiClient

    mymodel = BmiClient(grpc.insecure_channel("localhost:<PORT>"))


For the example model launched in :ref:`python-example`, the component name can be retrieved following the usual BMI syntax,

.. code-block:: python

    print(mymodel.get_component_name())
    Hello world


Python Subprocess
.................

This python class launches a BMI server upon creation,

.. code-block:: python

    from grpc4bmi.bmi_client_subproc import BmiClientSubProcess

    model = BmiClientSubProcess(<PACKAGE>.<MODULE>.<CLASS>)


The code above will execute ``run-bmi-server`` in a python subprocess and automatically listen to the appropriate port. Note that this requires your client to run in the same python environment as your model.

:ref:`running-python` Python server explains the roles of ``<PACKAGE>``, ``<MODULE>`` and ``<CLASS>``.

Polyglot CLI
------------

Once you have started a GRPC server you can test it by connecting to it using the `Polyglot - a universal grpc command line client`_.

Polyglot requires Java and the `polglot.yar` file can be downloaded at https://github.com/dinowernli/polyglot/releases

The following commands expects a GRPC server running on localhost on port 55555.

To get the component name use

.. code-block:: sh

    echo '{}' | java -jar polyglot.jar call --endpoint=localhost:55555 --full_method=bmi.BmiService/getComponentName


.. _Polyglot - a universal grpc command line client: https://github.com/grpc-ecosystem/polyglot
grpc4bmi package
================

.. automodule:: grpc4bmi
   :members:
   :undoc-members:
   :show-inheritance:

Submodules
----------

grpc4bmi.bmi\_client\_docker module
-----------------------------------

.. automodule:: grpc4bmi.bmi_client_docker
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_client\_singularity module
----------------------------------------

.. automodule:: grpc4bmi.bmi_client_singularity
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_client\_subproc module
------------------------------------

.. automodule:: grpc4bmi.bmi_client_subproc
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_grpc\_client module
---------------------------------

.. automodule:: grpc4bmi.bmi_grpc_client
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_grpc\_legacy\_server module
-----------------------------------------

.. automodule:: grpc4bmi.bmi_grpc_legacy_server
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_grpc\_server module
---------------------------------

.. automodule:: grpc4bmi.bmi_grpc_server
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.bmi\_r\_model module
-----------------------------

.. automodule:: grpc4bmi.bmi_r_model
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.reserve module
-----------------------

.. automodule:: grpc4bmi.reserve
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.run\_server module
---------------------------

.. automodule:: grpc4bmi.run_server
   :members:
   :undoc-members:
   :show-inheritance:

grpc4bmi.utils module
---------------------

.. automodule:: grpc4bmi.utils
   :members:
   :undoc-members:
   :show-inheritance:

Command line tools
==================

run-bmi-server
--------------

.. argparse::
   :module: grpc4bmi.run_server
   :func: build_parser
   :prog: run-bmi-server
.. grpc4bmi documentation master file, created by
   sphinx-quickstart on Thu Mar 28 10:52:55 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grpc4bmi's documentation!
====================================

.. include:: introduction.rst

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing
   server/index
   usage
   container/building
   container/usage
   cli
   python_api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Introduction
============

The Basic Modeling Interface (BMI, see https://github.com/csdms/bmi) is a multi-language library interface tailored to earth system models. This software allows you to wrap a BMI implementation in a server process and communicate with it via the included python client. The communication is serialized to protocol buffers by GRPC (https://grpc.io/) and occurs over network ports. On the server side, we support BMI implementations in :ref:`python <pythonservice>`, R or C/C++. Fortran models need to be linked against the C-version of the BMI. On the client side, we expose the standard python BMI (https://github.com/csdms/bmi-python/).

This setup enables you to wrap your BMI-enabled model in a Docker (https://www.docker.com/) or Singularity (https://singularity.lbl.gov/) container  and communicate with it from a python process on the host machine.

.. image:: _static/design-overview.svg
Python API
==========

.. toctree::
   :maxdepth: 4

   grpc4bmi
R
=

Grpc4bmi allows you to wrap a Hydrological model written in the `R language`_ into a GRPC server.

.. _R language: https://www.r-project.org/


Installing Requirements
-----------------------

The `bmi-r`_ package can be installed using the following `devtools`_ command

.. code-block:: R

    devtools::install_github("eWaterCycle/bmi-r")


Creating
--------

A model must implement the basic model interface (bmi).

This can be done by sub-classing the AbstractBmi class found in the `bmi-r`_ R package.


A model (in the example called `mymodel`) can than be given a basic model interface with something like

.. code-block:: R

    library(R6)
    library(bmi)
    library(mymodel)

    MyModelBmi <- R6Class(
        inherit = AbstractBmi,
        public = list(
            getComponentName = function() return('mymodel'),
            bmi_initialize = function(config_file) {
                # TODO Construct & initialize mymodel model
            },
            update = function() {
                # TODO evolve mymodel model to next time step
            },
            # TODO implement all other bmi functions
        )
    )


For an example of a BMI interface of the `Wageningen Lowland Runoff Simulator (WALRUS)`_ see `walrus-bmi.r`_

.. _bmi-r: https://github.com/eWaterCycle/bmi-r
.. _devtools: https://devtools.r-lib.org/
.. _Wageningen Lowland Runoff Simulator (WALRUS): https://github.com/ClaudiaBrauer/WALRUS
.. _walrus-bmi.r: https://github.com/eWaterCycle/grpc4bmi-examples/blob/master/walrus/walrus-bmi.r

Running
-------

Once the model has an BMI interface it can be run as a GRPC server by installing the `grpc4bmi[R]` Python package with

.. code-block:: sh

    pip install grpc4bmi[R]

The server can be started with

.. code-block:: sh

    run-bmi-server --lang R [--path <R file with BMI model>] --name [<PACKAGE>::]<CLASS> --port <PORT>

For the WALRUS model the command is

.. code-block:: sh

    run-bmi-server --lang R --path ~/git/eWaterCycle/grpc4bmi-examples/walrus/walrus-bmi.r --name WalrusBmi --port 55555

The Python grpc4bmi :ref:`usage` can then be used to connect to the server.
Note that the ``--port`` and ``--path`` arguments also can be specified as the respective environment variables ``BMI_PORT`` and ``BMI_PATH``.
.. _pythonservice:

Python
======

If you have a BMI-compliant model written in python, grpc4bmi provides a quick way to set up a BMI service.

Installing Requirements
-----------------------

The grpc4bmi Python package should be :ref:`installed <pip-install>`.


Creating
--------

To obtain a python BMI for your model, install the `Python bmi package (bmipy) <https://pypi.org/project/bmipy/>`_ and implement the :class:`bmipy.Bmi` abstract base class for your model. For exposing this model as a GRPC service, it is necessary to have a constructor without arguments: all initialization state will be presented to the model via the configuration file in the ``initialize`` method.

.. _running-python:

Running
-------

The installation of the grpc4bmi package installs the ``run-bmi-server`` command. You can run your model as a service by typing

.. code-block:: sh

    $ run-bmi-server --name <PACKAGE>.<MODULE>.<CLASS>

where ``<PACKAGE>``, ``<MODULE>`` are the python package and module containing your python BMI model, which should contain a python class ``<CLASS>`` that implements Bmi. The script assumes that this class does not take any constructor arguments. Upon running, the server will report which networking port it has decided to use on the terminal. This port will later be needed by BMI clients to communicate with your service.
The port can also be specified by adding the option ``--port <PORT>`` or pre-define the environment variable ``BMI_PORT`` (the latter takes precedence over the former).
An extra system path can be specified by adding the option ``--path <PATH>`` or pre-define the environment variable ``BMI_PATH``.


.. _python-example:

Example
-------

As an example, suppose we have a package

.. code-block:: sh

    $ mypackage
    $ - __init__.py
    $ - mymodule.py

and inside the ``mymodule.py`` the bmi implementation

.. code-block:: python

    from bmi import Bmi

    class MyBmi(Bmi):
        def __init__(self):
            ...
        def get_component_name(self):
            return "Hello world"

Then we launch this toy model as a service by executing

.. code-block:: sh

    $ run-bmi-server --name mypackage.mymodule.MyBmi

This will report the chosen port number in the standard output stream. It can be used to connect to the service via the BMI :ref:`grpc python client <python-grpc4bmi-client>`.Creating a BMI server
=====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   python
   R
   Cpp
C/C++/Fortran
=============

.. _install_cpp:

Installing Requirements
-----------------------
For native programming languages it is necessary to install and compile the C++ bindings of GRPC and protobuf on your system:

.. code-block:: sh

    git clone -b $(curl -L https://grpc.io/release) --depth=1 https://github.com/grpc/grpc
    cd grpc
    git submodule update --init --recursive
    wget -q -O cmake-linux.sh https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5-Linux-x86_64.sh
    sudo sh cmake-linux.sh -- --skip-license --prefix=/usr/local
    rm cmake-linux.sh
    mkdir cmake/build && cd cmake/build
    /usr/local/bin/cmake ../.. -DgRPC_INSTALL=ON -DgRPC_SSL_PROVIDER=package -DgRPC_BUILD_TESTS=OFF -DBUILD_SHARED_LIBS=ON
    sudo make -j4 install
    sudo ldconfig

You will also need to compile grpc4bmi

.. code-block:: sh

    git clone --depth=1 https://github.com/eWaterCycle/grpc4bmi.git
    cd grpc4bmi && git submodule update --init
    cd cpp
    mkdir -p build && cd build && cmake .. && sudo make install


Creating
--------

The grpc4bmi package comes with a C++ abstract base class that contains the BMI functions. The `header file <https://github.com/eWaterCycle/grpc4bmi/blob/master/cpp/bmi_class.h>`_ will
be copied to your system include path upon the installation steps above. Write an implementation of the ``Bmi`` class using your model time step code and data structures. You don't have to worry about global variables in your model code: with grpc4bmi every model instance runs in its own memory space. For the same reason, the ``get_value_ptr`` and ``set_value_ptr`` methods can be safely ignored, they are never called through the grpc process bridge.

Running
-------

Since native language lack reflection, it is necessary to make your own ``run_bmi_server`` program. We provide a function ``run_bmi_server(Bmi*, int*, char*)`` in the ``bmi_grpc_server.h`` header that can be called with your model instance (see the example below). To compile your server binary, it is necessary to link against grpc4bmi and protobuf libraries.
The program will accept a single optional argument which is the port the server will run on. The port can also be specified using the BMI_PORT environment variable. The default port is 50051.

.. _example_cpp:

Example
-------

To create a BMI to your model, write a header file in which you declare the overridden functions of the base class ``Bmi`` in the included file ``bmi_class.h``.

my_bmi_model.h:

.. code-block:: cpp

    #include <bmi_class.h>

    class MyBmiModel: public Bmi
    {
        public:
            MyBmiModel();
            int initialize(const char* config_file) override;
            ...
            int get_component_name(char* name) const override;
    };

Write your implementation of the basic modeling interface in the corresponding source file

my_bmi_model.cc:

.. code-block:: cpp

    #include <my_bmi_model.h>
    #include <cstring>

    MyBmiModel::MyBmiModel(){}
    int MyBmiModel::initialize(const char* config_file)
    {
        /* ...initialize the model from config_file... */
        return BMI_SUCCESS;
    }
    ...
    int MyBmiModel::get_component_name(char* name) const
    {
        strcpy(name, "Hello world");
        return BMI_SUCCESS;
    }

Now the BMI server can be simply be implemented as

run_my_bmi_model.cc:

.. code-block:: cpp

    #include "bmi_grpc_server.h"
    #include "my_bmi_model.h"

    int main(int argc, char* argv[])
    {
        Bmi* model = new HypeBmi();
        run_bmi_server(model, argc, argv);
        delete model;
        return 0;
    }

This binary will need to be linked against grpc4bmi and the protobuf libraries:

.. code-block:: sh

    g++ -o my_bmi_server run_my_bmi_model.o my_bmi_model.o `pkg-config --libs protobuf grpc++ grpc` -Wl,--no-as-needed -lgrpc++_reflection -ldl -lgrpc4bmi



Fortran
.......

In case you have a Fortran model, we advice to write the corresponding functions in Fortran first and export them to the implementation, e.g.

my_bmi_model.f90:

.. code-block:: fortran

    subroutine get_component_name(name) bind(c, name="get_component_name_f")
        use, intrinsic ::iso_c_binding
        implicit none
        character(kind=c_char), intent(out) :: name(*)
        name(1:11)="Hello world"
        name(12)=c_null_char

Now it is possible to call this function from the BMI C implementation as follows,

my_bmi_model.cc:

.. code-block:: cpp

    extern "C" void get_component_name_f(char*)
    int MyBmiModel::get_component_name(char* name) const
    {
        get_component_name_f(name);
        return BMI_SUCCESS;
    }
Using the container clients
===========================

.. _docker_client:

Docker
------

Grpc4bmi can run containers with `Docker engine`_.

Use the :class:`grpc4bmi.bmi_client_docker.BmiClientDocker` class to start a Docker container and get a client to interact with the model running inside the container.



For example the PCR-GLOBWB model can be started in a Docker container with

.. code-block:: python

    model = BmiClientDocker(image='ewatercycle/pcrg-grpc4bmi:latest', image_port=55555,
                            input_dir="./input",
                            output_dir="./output")
    # Interact with model
    model.initialize('config.cfg')

    # Stop container
    del model

.. _Docker engine: https://docs.docker.com/

Singularity
-----------

Grpc4bmi can run containers on `Singularity`_.

The Docker images build :ref:`previously <building-docker-image>` can be either run directly or converted to singularity image file and run.

To run a Docker image directly use `docker://<docker image name>` as singularity image name.

To convert a Docker image to a singularity image file use

.. code-block:: sh

    singularity build  docker://<docker image name> <singularity image filename>


Use the :class:`grpc4bmi.bmi_client_singularity.BmiClientSingularity` class to start a Singularity container and get a client to interact with the model running inside the container.

.. code-block:: python

    from grpc4bmi.bmi_client_singularity import BmiClientSingularity
    image = '<docker image name of grpc4bmi server of a bmi model>'
    client = BmiClientSingularity(image, input_dir='<directory with models input data files>')

For example for the wflow Docker image the commands would be the following

.. code-block:: python

    from grpc4bmi.bmi_client_singularity import BmiClientSingularity
    image = 'docker://ewatercycle/wflow-grpc4bmi:latest'
    client = BmiClientSingularity(image, input_dir='wflow_rhine_sbm', output_dir='wflow_output')

.. _Singularity: https://www.sylabs.io/guides/latest/user-guide/
.. _building-docker-image:

Building a docker image
=======================

The biggest advantage of using grpc4bmi is that you can embed the model code in a container like a `Docker`_ image. The grpc bridge allows you to address it from the host machine with the python BMI.

To establish this, install your BMI model and grpc4bmi inside the container, and let ``run-bmi-server`` act as the entry point of the docker image.


Python
------

The docker file for the model container simply contains the installation instructions of grpc4bmi and the BMI-enabled model itself, and as entrypoint the ``run-bmi-server`` command. For the :ref:`python example <python-example>` the Docker file will read

.. code-block:: Dockerfile

    FROM ubuntu:bionic
    MAINTAINER your name <your email address>

    # Install grpc4bmi
    RUN pip install git+https://github.com/eWaterCycle/grpc4bmi.git#egg=grpc4bmi

    # Install here your BMI model:
    RUN git clone <MODEL-URL> /opt/mymodeldir

    # Run bmi server
    ENTRYPOINT ["run-bmi-server", "--name", "mypackage.mymodule.MyBmi", "--path", "/opt/mymodeldir"]

    # Expose the magic grpc4bmi port
    EXPOSE 55555

The port 55555 is the internal port in the Docker container that the model communicates over. It is the default port for ``run_bmi_server`` and also the default port that all clients listen to.

R
-

The Docker image can be made by writing a `Dockerfile` file like

.. code-block:: Dockerfile

    FROM r-base
    LABEL maintainer="Your name <your email address>"

    RUN apt update && apt install -t unstable -y python3-dev python3-pip git && \
      pip3 install git+https://github.com/eWaterCycle/grpc4bmi.git#egg=grpc4bmi[R]

    RUN install.r remotes && installGithub.r eWaterCycle/bmi-r
    RUN install.r <R mymodel library from CRAN>

    # Copy BMI interface of model into Docker image
    RUN mkdir /opt/
    COPY mymodel-bmi.r /opt/

    # Config file and forcing file will be mounted at /data
    RUN mkdir /data
    WORKDIR /data
    VOLUME /data

    ENV BMI_PORT=55555

    CMD ["run-bmi-server", "--lang", "R", "--path", "/opt/mymodel-bmi.r", "--name", "mymodel"]

    EXPOSE 55555


The WALRUS model has a `Dockerfile`_  file which can be used as an example.

.. _Dockerfile: https://github.com/eWaterCycle/grpc4bmi-examples/blob/master/walrus/Dockerfile

C/C++/Fortran
-------------

For native languages you need to compile you BMI model inside the container let your bmi server runner binary program act as the entry point. The protobuf, grpc and grpc4bmi libraries need to be installed in your docker image, which means that the :ref:`installation instructions <install_cpp>` must be adopted in your Docker file. Then, include the installation of the model itself and the bmi run binary that you have written (as described :ref:`here <example_cpp>`). Finally the entry point in the docker file should be the launch of this binary and you should expose port 55555. For the C++ example :ref:`C++ example <example_cpp>`

.. code-block:: Dockerfile

    # ...download, compile and install grpc and grpc4bmi...
    # ...download, compile and install my_bmi_model...
    # Run bmi server
    ENTRYPOINT ["my_bmi_server"]

    # Expose the magic grpc4bmi port
    EXPOSE 55555

Building and Publishing
-----------------------

The Docker image can be build with

.. code-block:: sh

    docker build -t <image name> .

The Docker image can be published at `Docker Hub`_ by creating a repository and pushing it with

.. code-block:: sh

   docker push <image name>

The example WALRUS model is published at https://cloud.docker.com/u/ewatercycle/repository/docker/ewatercycle/walrus-grpc4bmi.

The Docker image can then be started with the grpc4bmi :ref:`docker client <docker_client>`.

.. _Docker: https://docs.docker.com/
.. _Docker Hub: https://hub.docker.com/
