# MagnumPI

Potential Integrator - A package for the calculation of scattering angles and cross sections.

The web application is hosted on [https://magnumpi.gitlab.io/magnumpi](https://magnumpi.gitlab.io/magnumpi).

For installation instructions see [INSTALL.md](INSTALL.md).

## Directory layout

* [magnumpi/](magnumpi/), C++ header files
* [src/](src/), C++ source code
* [external/](external/), C++ dependencies
* [app/](app/), C++ source code for command line tools
* [testsuite/](testsuite/), C++ tests
* [python/](python/), Python bindings
* [pymagnumpi/](pymagnumpi/), Python package and Python command line tools
* [web/assembly/](web/assembly/), JavaScript npm package called `@magnumpi/wasm` with WebAssembly version of C++ library compiled with [emscripten](https://emscripten.org/)
* [web/magnumpi-app/](web/magnumpi-app/), [Single Page web application](https://en.wikipedia.org/wiki/Single-page_application) using `@magnumpi/wasm` npm package
* [web/ws/python-connexion/](web/ws/python-connexion/), Web service written in Python
* [web/ws/cpp-pistache-server/](web/ws/cpp-pistache-server/), Web service written in C++
* [web/cgi/](web/cgi/), C++ CGI scripts

During CI job the following files get generated:

* [doc/notes/notes.pdf](https://gitlab.com/magnumpi/magnumpi/-/jobs/artifacts/master/file/build/doc/notes/notes.pdf?job=cpp)
# Quick compilation instructions

## I. Requirements

To compile the code you need

* [cmake](https://cmake.org/)
* C++ compiler like [GCC](https://gcc.gnu.org/)
* (optional) [pybind11](https://pybind11.readthedocs.io/) for Python bindings
* (optional) [emscripten](https://emscripten.org/) for WebAssembly module
* (optional) pdflatex, doxygen and/or tex4ht for doc generation

## II. After cloning, create a build directory, cmake and make

```shell
cd /path/to/magnumpi
mkdir build
cd build
cmake ..
make
```

## III. Prepare the environment and run one of the example scripts:

```shell
. set-magnumpi-local
app/calc_scat ../input/scat_hs.json
```

NOTES:

* Step II can be repeated if you want to support multiple builds.
  As an example, a debug and doc build can be set up as follows:

  ```shell
  cd /path/to/magnumpi
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -Dwith-doc=1 ..
  make
  ```

* After updating, repeat the 'make' step in your build directories.

* The steps above assume that you run the code in the build directory.
  If you plan to install the code, add a --prexix option to configure
  and type 'make install' after make. In the installation directory,
  use bin/set-magnumpi instead of set-magnumpi-local to prepare the
  environment.

* Various bits of documentation can be built in the doc/ subdirectory
  (in the build directory). To that end, install pdflatex, doxygen
  and/or tex4ht, reconfigure and do a make inside doc/. The LaTeX
  compilation may silently fail if not all its required sub-packages are
  installed. In order to get the LateX errors, go to the LaTeX directory
  (e.g. doc/paper1) and type 'make LATEX_INTERACTION_MODE='. Without
  that environment setting, 'make' will instruct LaTeX to ignore errors
  and bail out. tex4ht will create xml/html and xml/mathml versions of
  the documents. When doxygen is available, source code documentation
  will be extracted from annotations in the code. The documentation will
  appear in doc/doxygen.
# Web assembly of MagnumPI

The JavaScript npm package of the MagnumPI C++ library.
MagnumPI can be used to calculate scattering angles and several types of cross-sections between 2 particles.

This directory contains a [emscripten bindings](https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html) file (`src/wa.cpp`) to make C++ functions available to JS using a [WebAssembly](https://webassembly.org/) file.

## Installation

Install using

```shell
npm install @magnumpi/wasm
```

Or for [yarn](https://yarnpkg.com/) users

```shell
yarn add @magnumpi/wasm
```

## Run

The library can be run in [NodeJS](https://nodejs.org/) and [Web Worker](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API) environments. See below for a number of examples.

### Using CommonJS require and promises

```javascript
const magnumpicalc = require('@magnumpi/wasm')();
// Wait for wasm to be loaded
magnumpicalc.then((magnumpi) => {
    console.log(magnumpi.eV);

    // Convert eV to Joule
    const joule = magnumpi.in_si('{"value": 1e+5, "unit": "eV"}', magnumpi.dimensions.energy);
    console.log(joule);
});
```

### Using dynamic import and await

```shell
node --experimental-repl-await
```

```js
// Using dynamic import and ES modules
const magnumpicalc = await import('@magnumpi/wasm');
const magnumpi = await magnumpicalc.default();
magnumpi.eV
```

### Using ES module file

Create test file called `ev.mjs` with

```js
import magnumpicalc from '@magnumpi/wasm';
const magnumpi = await magnumpicalc();
magnumpi.eV
```

Run with

```shell
node ev.mjs
```

### Use in web application web worker

See [../magnumpi-app/](https://gitlab.com/magnumpi/magnumpi/-/tree/master/web/magnumpi-app/README.md) how to use the library in a web worker of a web application.

### Example calculation

In the [src/magnumpi.test.js](https://gitlab.com/magnumpi/magnumpi/-/tree/master/web/assembly/src/magnumpi.test.js) file there is a test that calculates the cross section of a FiniteWell potential, the test shows how the available methods can be connected together.

## Development

### Requirements

Requires

* emscripten SDK see https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#sdk-download-and-install
* make

The `emcc` executable of the emscripten SDK should be in your PATH.

### Dependencies

Install development dependencies with

```bash
npm install
```

### Build

To build `dist/` folder with wasm file and js wrappers use

```bash
npm run build
```

To clean `dist/` folder use

```shell
npm run clean
```

## Tests

There are also unit tests written in [jest](https://jestjs.io/) test framework which can be run with

```shell
npm run test
```
[OpenAPI](https://www.openapis.org/) specification of MagnumPI web service.

For each calc_* executable maps config file to a request and output files to a response.

# Generated server stubs

First install
```bash
npm install @openapitools/openapi-generator-cli 
```

## C++

```
npx openapi-generator generate -i openapi.yml -g cpp-pistache-server -o cpp-pistache-server \
--invoker-package magnumpi.ws --api-package magnumpi.ws.api \
--model-package magnumpi.ws.model -additional-properties=helpersPackage=magnumpi.ws.helpers
```

## Python

Instead of using generator, make use of https://connexion.readthedocs.io/ directly.
Also generator does not handle potential oneof

```
# npx openapi-generator generate -i openapi.yml -g python-flask -o python-flask
cookiecutter https://github.com/nlesc/python-template.git 
```

## PHP

```
npx openapi-generator generate -i openapi.yml -g php-lumen -o php-lumen
```

## NodeJS

Nodejs server generator is currently broken
```
npx openapi-generator generate -i openapi.yml -g nodejs-server -o nodejs-server
```
# MagnumPI web service

Web service which exposes endpoints to run MagnumPI calculations.

* OpenAPI specification at /openapi.json
* OpenAPI UI Console at /ui
* Quick calculations run directly.
* Slower calculations run using [celery job queue](https://docs.celeryproject.org/en/stable/).

## Install & run

```shell
# Install magnumpi Python package in current Python environment by following steps in ../../pymagnumpi/DEVELOP.md
# Install other deps
pip install .
# Star Task queue
docker run -d -p 6379:6379 redis
# Start task worker
export MAGNUMPI_INPUTDATA_DIR=$PWD/../../..
celery -A magnumpi_webservice.tasks worker
# Start webservice in another terminal
gunicorn -w 4 -b 0.0.0.0:8888 magnumpi_webservice.serve:app
```

Will run web service on http://localhost:8888, with the swagger ui on http://localhost:8888/ui

Run with docker-compose from root of repo with

```shell
docker-compose -f web/ws/python-connexion/docker-compose.yml up
```

Which will run web service on http://localhost:8888

## Development

### Install dependencies

```shell
pip install -e .[dev]
```

### Run webservice in debug mode

```shell
python -m magnumpi_webservice.serve
```

### Tests

Test can be run with

```shell
pytest
```

### Lint

Lint with pycodestyle

```shell
pycodestyle src/magnumpi_webservice
```
# REST API Server for Potential Integrator

## Overview
This API Server was generated by the [OpenAPI Generator](https://openapi-generator.tech) project.
It uses the [Pistache](https://github.com/oktal/pistache) Framework.

## Files organization
The Pistache C++ REST server generator creates three folders:
- `api`: This folder contains the handlers for each method specified in the OpenAPI definition. Every handler extracts
the path and body parameters (if any) from the requests and tries to parse and possibly validate them.
Once this step is completed, the main API class calls the corresponding abstract method that should be implemented
by the developer (a basic implementation is provided under the `impl` folder)
- `impl`: As written above, the implementation folder contains, for each API, the corresponding implementation class,
which extends the main API class and implements the abstract methods.
Every method receives the path and body parameters as constant reference variables and a reference to the response
object, that should be filled with the right response and sent at the end of the method with the command:
response.send(returnCode, responseBody, [mimeType])
- `model`: This folder contains the corresponding class for every object schema found in the OpenAPI specification.

The main folder contains also a file with a main that can be used to start the server.
Of course, is you should customize this file based on your needs

## Installation
First of all, you need to download and install the libraries listed [here](#libraries-required).

Once the libraries are installed, in order to compile and run the server please follow the steps below:
```bash
mkdir build
cd build
cmake ..
make
```

Once compiled run the server:

```bash
cd build
./api-server
```

Test with:
```
curl -d '@i../../../../input/scat_hs_colonna.json' -H "Content-Type: application/json" -X POST  http://localhost:8080/scattering_angle
```

## Libraries required
- [pistache](http://pistache.io/quickstart)
- [JSON for Modern C++](https://github.com/nlohmann/json/#integration): Please download the `json.hpp` file and
put it under the model folder

## Namespaces
magnumpi.ws.api
magnumpi.ws.model
# Calculators as CGI scripts

## Hello world

Compile

```sh
g++ -I../../pi_support -std=c++14 src/hello.cpp -o cgi-bin/hello
```

Test Directly

```sh
echo '{"name":"me"}' | ./cgi-bin/hello
Content-type: application/json
{"message":"Hello World!","name":"me"}
```

```sh
docker build -t magnumpi:cgi .
docker run -p 8888:80 magnumpi:cgi
curl -X POST -H "Content-Type: application/json" -d '{"name":"me"}' http://localhost:8888/cgi-bin/hello
```

## Scattering angle

Compile

```sh
g++ -std=c++14 -I../../pi_support -I../.. -o cgi-bin/calc_scat_colonna ../../src/beast_post.cpp ../../src/json_nlohmann.cpp ../../src/units.cpp ../../src/data_set.cpp ../../src/spline.cpp ../../src/potential.cpp ../../src/scat_colonna.cpp ../../src/scat_colonna_dcs.cpp src/calc_scat_colonna.cpp
```

Test directly

```sh
cat ../../input/He_He+/input/gerade_scat_inline.json | cgi-bin/calc_scat_colonna
```

Against Docker container:

```sh
curl -X POST -H "Content-Type: application/json" -d @../../input/He_He+/input/gerade_scat_inline.json http://localhost:8888/cgi-bin/calc_scat_colonna
```
# MagnumPI web application

![screenshot](magnumpi-webapp-screencast.gif "Screenshot")

## Available Scripts

To use the LXCatJSONOnline potential type the LXCat web service must be running on a known endpoint, default is `http://localhost:3000`.

### Dependencies

Installs dependencies of application and build wasm.

```shell
cd ..
# From web/
# Install all dependepencies
npm install
# Build web assembly
npm run build --workspace=assembly
cd -
```

By using [npm packages](https://docs.npmjs.com/cli/v7/using-npm/workspaces) the assembly/ dir is included in the node_modules of magnumpi-app.

### Development mode

```shell
npm run dev
```

Runs the app in the development mode.
Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

The page will reload if you make edits.
You will also see any lint errors in the console.

### Tests

To run unit tests do

```shell
npm run test:unit
```

### Linting & formatting

The code should be formatted using [prettier](https://prettier.io/) and have no [eslint](https://eslint.org/) errors.

Automatic formatting can be done with

```shell
npm run format
```

Linting can be done with

```shell
npm run lint
```

See prettier and eslint websites for how to configure editor to perform automatic linting and formatting.

### Build

The build command requires bash and perl to fix up the generated files.

To make a production build use

```shell
npm run build
```

Builds the app for production to the `dist/` folder.
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.
Your app is ready to be deployed!

For local deployment use

```bash
npm run serve
```

Open [http://localhost:5000](http://localhost:5000) to view it in the browser.

## Docker

The web application can be build & run using Docker

To build from repo root execute following

```bash
cd ../..
docker build -t passingxsams/magnumpi:webapp -f web/magnumpi-app/Dockerfile .
```

To run do

```bash
docker run -p 5000:5000 passingxsams/magnumpi:webapp
```

Will run web application on [http://localhost:5000](http://localhost:5000)

## Credits

This project was bootstrapped with [Vite](https://vitejs.dev/).
The code in this directory aims to create tables chi(b) for a given energy E,
using the Viehland algorithm. At present is uses a mix of Viehland and Colonna
calculations, and the setting 'use Colonna' setting is not always respected.
The code also uses a b-grid as if the result is going to be used for a cross
section calculation (the b-points are the positions of a Clemshaw-Curtis
integration). This may be unexpected for a user who specifies a b-range in the
iput file.
# To build Python with cmake and conda

## Install deps in conda env

(Boost should not be installed in OS. Otherwise import problems)

```shell
. ~/miniconda39/bin/activate
conda install -n base -c conda-forge mamba
mamba create -y -n magnumpi -c conda-forge pybind11 pytest cmake compilers boost-cpp
conda activate magnumpi
```

## Compile

```shell
# In root
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -Dwith-doc=0  ..
make -j 4
```

## Run Python scripts

```shell
. set-magnumpi-local
python3 python/app/get_cs.py ../input/cs.json
```
# Scripts that use the Python Bindings for MagnumPI

This directory contains scripts that do calculations using
the python bindings for the MagnumPI code.

NOTES

 * Source set-magnumpi-local (which reside in the top-level build directory),
   that will set relevant environment variables. It will also append the
   python/ subdirectory of the build directory to the PYTHONPATH environment
   variable.

 * The Python3 script `get_scat.py` requires one input file name,
   for example `input/scat_rscp_colonna.json`. It calculates
   the scattering angle chi(b,E) for user-specified impact parameter
   (b) and energy (E) ranges. It calls the C++ function
   `magnumpi::get_scattering_angle`. The resulting JSON object
   is written to the standard output stream (by Pyton code).

 * The Python3 script `get_pot.py` requires one input file name,
   for example `input/HH.json`. It prints a table of values of the
   potential V(r) and its first and second derivatives for various
   values of the inter-molecular separation r.

 * The Python3 script `get_cs.py` requires one input file name,
   for example `input/cs.in`. It prints the requested cross sections
   `Q_l(eps)` for user-selected values `l` in `[0,l_max)` as a function
   of a set of user-specified energy values.

As an example, from the build directory the `get_scat.py` script can be
invoked as (everything on one line):

```
    python3 python/app/get_scat.py ../input/scat_rscp_colonna.json
```

# Python Bindings for MagnumPI

This directory contains python bindings for some of the interfaces
of MagnumPI and some associated native Python code.

NOTES

 * the Python bindings in backend.cpp are compiled into a
   shared library, which dynamically links against the
   MagnumPI library. The bindings use the pybind11 library.

 * For autotools build, the shared library that contains the bindings is
   produced in the `.libs/` subdirectory of the ideas/jvdijk/python build
   directory.
# Develop

Documentation for developers of the pymagnumpi package.

## Non-Python requirements

* C++ compiler
* make
* cmake
* [pybind11](https://pybind11.readthedocs.io/)
* boost

### Install with conda

The requirements can be installed using [miniconda](https://docs.conda.io/en/latest/miniconda.html) with

```shell
conda install -n base -c conda-forge mamba
mamba create -y -n magnumpi -c conda-forge pybind11 cmake compilers boost-cpp make
conda activate magnumpi
```

(Boost should not be installed in OS. Otherwise import problems)

## Compile shared library

```shell
# In root
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Dwith-doc=0  ..
make -j 4
cd pymagnumpi
make install
# installs *.so into /pymagnumpi/src/magnumpi
```

## Install Python requirements

```shell
# In pymagnumpi/
pip install -e .[dev]
```

## Tests

The tests in `tests/` folder can be run with

```shell
# In pymagnumpi/
pytest
```

To get test coverage (obviously will only cover Python code, not C++ code)

```shell
pytest --cov magnumpi --cov-branch --cov-report=html --cov-report=term --cov-report=xml
```

### Lint

Lint with pycodestyle

```shell
pycodestyle src/magnumpi
```

## Docs

Docs for CLI, API including _libmagnumpi.so can be build with

```shell
# In pymagnumpi/
cd docs
make html
xdg-open _build/html/index.html
```

TODO host docs on [readthedocs.io](https://readthedocs.io)

## Build

To build binary wheel for current OS.
Needs [shared library](#compile-shared-library) to have been compiled.

```shell
# In pymagnumpi/
python -m build --wheel .
```

TODO make for multiple OS / python versions, maybe use https://github.com/pypa/cibuildwheel or multiple conda environments.

## Upload to PyPI

```shell
# In pymagnumpi/
twine upload dist/*
```
# MagnumPI Python package

Python package that makes C++ MagnumPI library available in Python.

## Install

```shell
pip install magnumpi
```

## Command line tool

The `magnumpi` command line tool can be used to perform calculations. See below for some examples.

The example input files can be found in [https://gitlab.com/magnumpi/magnumpi](https://gitlab.com/magnumpi/magnumpi). The following commands are run in root of that repo.

The example input files include env vars that are expanded in C++, and those env vars needs to be set

```shell
export MAGNUMPI_INPUTDATA_DIR=$PWD
```

### Calculate cross section

```shell
magnumpi cs input/cs.json
```

It will output a JSON formatted string with the cross section at the requested energy values.

```json
{
  "title": "sigma_l(eps)",
  "unit": "m^2",
  "column_axis": {
    "title": "l",
    "labels": [
      "1",
      "2",
      "3"
    ]
  },
  "row_axis": {
    "title": "epsilon",
    "data": {
      "unit": "eV",
      "values": [
        9.999999999999961e-05,
...
        99999.9999999996
      ]
    }
  },
  "columns": [
    [
      8.092877898598035e-18,
      5.906052452597951e-18,
      9.479830033419527e-18
    ],
...
  ]
}
```

### Calculate scattering angle

```shell
magnumpi scat input/scat_rscp_colonna.json
```

It will output a JSON formatted string with the scattering angle at the requested impact and energy ranges.

```json
{
  "column_axis": {
    "data": {
      "unit": "eV",
      "values": [
        9.999999999999961e-05,
...
        99999.9999999996
      ]
    },
    "title": "E"
  },
  "columns": [
    [
      3.141592653589793,
...
      0.0002028565858154252
    ]
  ],
  "row_axis": {
    "data": {
      "unit": "m",
      "values": [
        0.0,
        1e-11,
...
        3e-10
      ]
    },
    "title": "b"
  },
  "title": "chi(b,E)",
  "unit": "rad"
}
```

### Calculate derivates of potential

```shell
magnumpi pot --rmin 2e-10 --rmax 3e-10 --rdel 1e-11 input/HH.json
#distance [m],V [V],dV/dr [V/m],d2V/dr2 [V/m^2]
2e-10,-1.0541933465153312e-19,4.225791224186558e-09,-104.47263840488723
2.1e-10,-6.843852916251855e-20,3.1738317744082248e-09,-102.61821539091953
2.2e-10,-4.160587069680212e-20,2.2191271972441033e-09,-86.91630626111781
2.2999999999999998e-10,-2.342850334973455e-20,1.450262373088819e-09,-66.65514992303468
2.4e-10,-1.1925072469313695e-20,8.828186881170004e-10,-47.2578142103886
2.5e-10,-5.176180186826646e-21,4.937041134467853e-10,-31.22042921959424
2.6000000000000003e-10,-1.5831173447527182e-21,2.4494767170692424e-10,-19.19034553797456
2.7000000000000005e-10,6.012825711190221e-23,9.760433642444005e-11,-10.835896548407579
2.8000000000000007e-10,5.952064614391767e-22,1.8400813738867518e-11,-5.428755491601337
2.900000000000001e-10,5.69477929651543e-22,-1.8149521987992936e-11,-2.1789629148688308
```

## Programmatic access

```python
from magnumpi import get_scattering_angles
reqeust = ...
response = get_scattering_angles(request)
print(response)
```

## Develop

For development of this package see [DEVELOP.md](DEVELOP.md) document.
# MagnumPI Input Files

This directory contains four types of input files:

 * Files with names of the form `*_int.json` describe interactions and contain
   a single element `interaction`. An `interaction` element contains two
   elements, a `potential` and a `reduced_mass`. The latter is needed for
   the calculation of collision integrals.

 * Files with names of the form `*_pot.json` describe potentialsand contain a
   single element `potential`. These suffice for the calculation of all other
   quantities: scattering angles, cross sections, reduced cross sections and
   reduced collision integrals.

 * Configuration files: these have names `*_scat.json`, `*_cs.json`,
   `*_cs_star.json`, `*_omega.json` and `*_omega_star.json`. These contain the
   numerical parameters that are needed to do a calculation of scattering
   angles, cross sections, reduced cross sections, collision integrals and
   reduced collision integrals. Typically a single configuration file can be
   used for multiple studies by combining it with different potential or
   interaction files.

 * Files with extension `_lut.dat`: these are old-style unannotated lookup
   tables that are referred to by some potential files.

Potential and interaction files do not necessarily describe a particular
real-world case: a Lennard Jones interaction file may use dummy parameter
values such as 1 Angstrom, 1 eV and 1 amu, and are used mostly for testing.
When such file describes a real interaction between two particles of a
particular type, this is made clear by the name of the file. As an example,
the file `heplus_he_ungerade_int.json` describes the potential for ungerade
interaction of a helium atom and a helium ion. In all cases, the files
contain annotation elements that describe the origin and purpose of the
file contents.

############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

1. use the search functionality `here <https://github.com//magnumpi_webservice/issues>`__ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

1. use the search functionality `here <https://github.com//magnumpi_webservice/issues>`__ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the MagnumPI webservice repository on GitHub;
1. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
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
reported by contacting the project team at s.verhoeven@esciencecenter.nl. All
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
Command line interface
======================

.. argparse::
   :module: magnumpi.cli
   :func: make_parser
   :prog: magnumpi
.. MagnumPI documentation master file, created by
   sphinx-quickstart on Fri Jul 30 11:10:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MagnumPI's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   cli
   API reference <api/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
