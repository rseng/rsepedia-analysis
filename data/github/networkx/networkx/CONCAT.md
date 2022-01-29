# pip requirements files

## Index

- [`default.txt`](default.txt)
  Default requirements
- [`extra.txt`](extra.txt)
  Optional requirements that may require extra steps to install
- [`example.txt`](example.txt)
  Requirements for gallery examples
- [`test.txt`](test.txt)
  Requirements for running test suite
- [`doc.txt`](doc.txt)
  Requirements for building the documentation (see `../doc/`)
- [`developer.txt`](developer.txt)
  Requirements for developers
- [`release.txt`](release.txt)
  Requirements for making releases

## Examples

### Installing requirements

```bash
$ pip install -U -r requirements/default.txt
```

### Running the tests

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/test.txt
```
<!--
Please run black to format your code.
See https://networkx.org/documentation/latest/developer/contribute.html for details.
-->
---
name: Bug report
about: 'Please describe the problem you have encountered'
---

<!-- If you have a general question about NetworkX, please use the discussions tab to create a new discussion -->

<!--- Provide a general summary of the issue in the Title above -->


### Current Behavior
<!--- Tell us what happens instead of the expected behavior -->

### Expected Behavior
<!--- Tell us what should happen -->

### Steps to Reproduce
<!--- Provide a minimal example that reproduces the bug -->

### Environment
<!--- Please provide details about your local environment -->
Python version:
NetworkX version:


### Additional context
<!--- Add any other context about the problem here, screenshots, etc. -->
# Building docs

We use Sphinx for generating the API and reference documentation.

Pre-built versions can be found at

    https://networkx.org/

for both the stable and the latest (i.e., development) releases.

## Instructions

After installing NetworkX and its dependencies, install the Python
packages needed to build the documentation by entering::

    pip install -r requirements/doc.txt

in the root directory.

To build the HTML documentation, enter::

    make html

in the ``doc/`` directory.  This will generate a ``build/html`` subdirectory
containing the built documentation.

To build the PDF documentation, enter::

    make latexpdf

You will need to have LaTeX installed for this.
