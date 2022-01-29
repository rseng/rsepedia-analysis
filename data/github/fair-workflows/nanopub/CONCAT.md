# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [1.2.10] - 2021-09-01
### Changed
* Use latest yatiml version instead of pinned version

## [1.2.9] - 2021-09-01
### Added
* Include LICENSE file in python setup

## [1.2.8] - 2021-09-01

### Added
* Also publish sdist when publishing to pypi

## [1.2.7] - 2021-06-25

### Changed
* Prevent `setup_nanopub_profile` from ever overwriting key pair

### Fixed
* Pin `click` at version 7.1.2, as versions >8 break `setup_nanopub_profile`

## [1.2.6] - 2021-04-30

### Added
* Search result dicts now contain nanopublication label too, if provided by the grlc endpoint.

## [1.2.5] - 2021-03-05

### Fixed
* Fix bug that overwrites optional pubinfo and prov in `from_assertion()` calls.

## [1.2.4] - 2021-03-04

### Fixed
* Fix bug where user rdf was being mutated during publishing.

## [1.2.3] - 2021-02-05

### Added
* Added new `publication_attributed_to` argument to `Publication.from_assertion()`. Allows the `pubinfo` attribution to be manually set if desired.

## [1.2.2] - 2021-01-29

### Fixed
* Fix FileAlreadyExists bug in `setup_nanopub_profile`

## [1.2.1] - 2021-01-22

### Changed
* Rename `setup_profile` to `setup_nanopub_profile` to avoid ambiguity/clashes with other tools

### Fixed
* Make `nanopub` package compatible Windows operating system
* Added UTF-8 related flags to nanopub-java (in java call) to fix issues with certain characters on certain java builds
* Make regex in orcid validation accept ids ending with 'X'

## [1.2.0] - 2020-12-23

### Added
* Added Zenodo badge to README
* Pagination of results for search methods of `NanopubClient`

### Changed
* `nanopub-java` dependency is installed upon installation instead of upon runtime.
* search methods of `NanopubClient` return iterator instead of list


## [1.1.0] - 2020-12-17

### Added
* `.zenodo.json` for linking to zenodo
* `pubkey` option to methods of `NanopubClient` that allows searching for publications 
    signed with the given pubkey. For these methods:
    - `find_nanopubs_with_text`
    - `find_nanopubs_with_pattern`
    - `find_things`
* `filter_retracted` option to methods of `NanopubClient` that allows searching for publications 
    that are note retracted. For these methods:
    - `find_nanopubs_with_text`
    - `find_nanopubs_with_pattern`
    - `find_things`
* `NanopubClient.find_retractions_of` method to search retractions of a given nanopublication.
* `Publication.signed_with_public_key` property: the public key that the publication was signed with.
* `Publication.is_test_publication` property: denoting whether this is a publicaion on the test server.

### Changed
* Improved error message by pointing to documentation instead of Readme upon ProfileErrors

### Fixed
* Catch FileNotFoundError when profile.yml does not exist, raise ProfileError with useful messageinstead.
* Fixed broken link to documentation in README.md

## [1.0.0] - 2020-12-08

NB: All changes before [1.0.0] are collapsed in here (even though there were multiple pre-releases)
### Added
- `nanopub.client` module with the NanopubClient class that implements:
  * Searching (being a client with a direct (but incomplete) mapping to the nanopub server grlc endpoint):
    * `find_nanopubs_with_text` method
    * `find_nanopubs_with_pattern` method
    * `find_things` method
  * Fetching:
    * `fetch` method to fetch a nanopublication
  * Publishing:
    * Publish a statement using `claim` method
    * Publish a `Publication` object with `publish` method
  * Retracting:
    * Publish a retraction of an existing nanopublication created by this user (i.e. signed with same RSA key)
  
  * Test server functionality
    * Client can optionally be set to publish to (and fetch from) the nanopub test servers only.

- `nanopub.publication` module
  * `Publication` class to represent a nanopublication. 
  Includes `from_assertion` class method to construct a Publication object
  from an assertion graph
  * `replace_in_rdf` helper method to replace values in RDF
- `nanopub.java_wrapper` module, provides an interface to the nanopub-java tool for
  signing and publishing nanopublications.
- `nanopub.profile` module, getters and setters for the nanopub user profile
- `nanopub.setup_profile`, interactive command-line client to setup user profile
- `nanopub.namespaces`, often-used RDF namespaces
- `examples/`, holds a few notebooks that serve as examples of using the library
- User documentation
![Build Status](https://github.com/fair-workflows/nanopub/workflows/Python%20application/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/nanopub/badge/?version=latest)](https://nanopub.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/fair-workflows/nanopub/badge.svg?branch=main)](https://coveralls.io/github/fair-workflows/nanopub?branch=main)
[![PyPI version](https://badge.fury.io/py/nanopub.svg)](https://badge.fury.io/py/nanopub)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4491/badge)](https://bestpractices.coreinfrastructure.org/projects/4491)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![DOI](https://zenodo.org/badge/302247101.svg)](https://zenodo.org/badge/latestdoi/302247101)

# nanopub
The ```nanopub``` library provides a high-level, user-friendly python interface for searching, publishing and retracting nanopublications.

Nanopublications are a formalized and machine-readable way of communicating
the smallest possible units of publishable information. See [the documentation](https://nanopub.readthedocs.io/en/latest/getting-started/what-are-nanopubs.html)
for more information.

# Documentation

Checkout the [user documentation](https://nanopub.readthedocs.io/)

# Setup
Install using pip:
```
pip install nanopub
```

To publish to the nanopub server you need to setup your profile. This allows the nanopub server to identify you. Run 
the following command in the terminal:
```
setup_nanopub_profile
```
This will ask you a few questions, then it will use that information to add and store RSA keys to sign your nanopublications with, (optionally) publish a nanopublication with your name and ORCID iD to declare that you are using using these RSA keys, and store your ORCID iD to automatically add as author to the
provenance of any nanopublication you will publish using this library.

## Quick Start


### Publishing nanopublications
```python

from nanopub import Publication, NanopubClient
from rdflib import Graph, URIRef, RDF, FOAF

# Create the client, that allows searching, fetching and publishing nanopubs
client = NanopubClient()

# Either quickly publish a statement to the server
client.claim('All cats are gray')

# Or: 1. construct a desired assertion (a graph of RDF triples)
my_assertion = Graph()
my_assertion.add( (URIRef('www.example.org/timbernerslee'), RDF.type, FOAF.Person) )

# 2. Make a Publication object with this assertion
publication = Publication.from_assertion(assertion_rdf=my_assertion)

# 3. Publish the Publication object. The URI at which it is published is returned.
publication_info = client.publish(publication)
print(publication_info)
```


### Searching for nanopublications
```python
from nanopub import NanopubClient

# Search for all nanopublications containing the text 'fair'
results = client.find_nanopubs_with_text('fair')
print(results)
```

### Fetching nanopublications and inspecting them
```python
# Fetch the nanopublication at the specified URI
publication = client.fetch('http://purl.org/np/RApJG4fwj0szOMBMiYGmYvd5MCtRle6VbwkMJUb1SxxDM')

# Print the RDF contents of the nanopublication
print(publication)

# Iterate through all triples in the assertion graph
for s, p, o in publication.assertion:
    print(s, p, o)

```
                                         
## Dependencies
The ```nanopub``` library currently uses the [```nanopub-java```](https://github.com/Nanopublication/nanopub-java) tool for signing and publishing new nanopublications. This is automatically installed by the library.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

## Our Standards

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

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at r.richardson@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/fair-workflows/nanopub/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/fair-workflows/nanopub/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. install the dependencies required to run `nanopub` in development:

    ```bash
    pip install -r requirements.txt
    pip install -r requirements-dev.txt
    ```

1. make sure the existing tests still work by running ``pytest``. Note that any pull requests to the nanopub repository on github will automatically trigger running of the test suite;
1. check that the code is in accordance with the PEP8 style guide, by running ``flake8 . --count --show-source --statistics``,
configuration is in `tox.ini`.
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the nanopub repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

## You want to update the nanopub-java dependency
Run `bin/nanopub-java --download` and move the generated 
nanopub-*-jar-with-dependencies.jar file to `nanopub/lib/` folder.
# Searching the nanopub server
The `NanopubClient` provides methods for searching the nanopub server. It provides
an (uncomplete) mapping to the [nanopub server grlc endpoint](http://grlc.nanopubs.lod.labs.vu.nl/api/local/local).

## Text search
Search for all nanopublications containing some text using 
`NanopubClient.find_nanopubs_with_text()`
```python
from nanopub import NanopubClient
client = NanopubClient()
results = client.find_nanopubs_with_text('fair')
```

## Triple pattern search
Search for nanopublications whose assertions contain triples that match a specific pattern.
```python
from nanopub import NanopubClient
client = NanopubClient()
# Search for nanopublications whose assertions contain triples that are ```rdf:Statement```s.
results = client.find_nanopubs_with_pattern(
                pred='http://www.w3.org/1999/02/22-rdf-syntax-ns#type',
                obj='http://www.w3.org/1999/02/22-rdf-syntax-ns#Statement')
```

## Search on introduced concept
Search for any nanopublications that introduce a concept of the given type, that contain 
text with the given search term.
```python
from nanopub import NanopubClient
client = NanopubClient()
# Search for nanopublications that introduce a concept that is a ```p-plan:Step```.
results = client.find_things('http://purl.org/net/p-plan#Step')
```

## Interpreting search results
Each search method returns a generator of dicts depicting matching nanopublications.

Each dict has the following key-value pairs:
* `date`: The date and time the nanopublication was created.
* `description`: A description of the nanopublication that was parsed from the nanopublication RDF.
* `np`: The URI of the matching nanopublication.

Example results (from `NanopubClient.find_nanopubs_with_text('fair')`):
```python
>>> print(list(results))
[{'date': '2020-05-01T08:05:25.575Z',
  'description': 'The primary objective of the VODAN Implementation Network is '
                 'to showcase the creation and deployment of FAIR data related '
                 'to COVID-19',
  'np': 'http://purl.org/np/RAdDKjIGPt_2mE9oJtB3YQX6wGGdCC8ZWpkxEIoHsxOjE'},
 {'date': '2020-05-14T09:34:53.554Z',
  'description': 'FAIR IN community',
  'np': 'http://purl.org/np/RAPE0A-NrIZDeX3pvFJr0uHshocfXuUj8n_J3BkY0sMuU'}]
```

## Returning retracted publications in search
By default nanopublications that have a valid retraction do not show up in search results.
A valid retraction is a retraction that is signed with the same public key as 
the nanopublication that it retracts. 
You can toggle this behavior with the `filter_retracted` parameter,
here is an example with `NanopubClient.find_nanopubs_with_text`:
```python
from nanopub import NanopubClient
client = NanopubClient()
# Search for nanopublications containing the text fair, also returning retracted publications.
results = client.find_nanopubs_with_text('fair', filter_retracted=False)
```

## Filtering search results for a particular publication key
You can filter search results to publications that are signed with
a specific publication key (effectively filtering on publications from a single author).
You use the `pubkey` argument for that. 
Here is an example with `NanopubClient.find_nanopubs_with_text`:
```python
from nanopub import NanopubClient, profile
# Search for nanopublications containing the text 'test',
# filtering on publications signed with my publication key.
client = NanopubClient(use_test_server=True)
my_public_key = profile.get_public_key()
results = client.find_nanopubs_with_text('test', pubkey=my_public_key)
```
# Fetching nanopublications
You can fetch nanopublications from the nanopub server using
`NanopubClient.fetch()`. The resulting object is a `Publication`
object that you can use to inspect the nanopublication.
```python
from nanopub import NanopubClient

# Fetch the nanopublication at the specified URI
client = NanopubClient()
publication = client.fetch('http://purl.org/np/RApJG4fwj0szOMBMiYGmYvd5MCtRle6VbwkMJUb1SxxDM')

# Print the RDF contents of the nanopublication
print(publication)

# Iterate through all triples in the assertion graph
for s, p, o in publication.assertion:
    print(s, p, o)

# Iterate through the publication info
for s, p, o in publication.pubinfo:
    print(s, p, o)

# Iterate through the provenance graph
for s, p, o in publication.provenance:
    print(s,p,o)

# See the concept that is introduced by this nanopublication (if any)
print(publication.introduces_concept)
```
# The nanopub test server
Throughout this documentation we make use of the 
[nanopub test server](http://test-server.nanopubs.lod.labs.vu.nl/)
by setting `use_test_server=True` when instantiating `NanopubClient`:
```python
>>> from nanopub import NanopubClient
>>> client = NanopubClient(use_test_server=True)
```
This will search and fetch from, and publish to the [nanopub test server](http://test-server.nanopubs.lod.labs.vu.nl/).

When learning about nanopub using the testserver is a good idea, because:
* You are free to experiment with publishing without polluting the production server.
* You can draft a publication and know exactly what it will look like on the nanopub server without polluting the production server.
* When searching (and to a lesser extent fetching) you are not putting an unnecessary load on the production server.

## Test purl URIs do not point to the test server
There is one caveat when using the test server that can be confusing:
The purl URI (for example: [http://purl.org/np/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8](http://server.nanopubs.lod.labs.vu.nl/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8))
points to the [nanopub production server](http://server.nanopubs.lod.labs.vu.nl/) 
resulting in a 404 page not found error.

A manual workaround is:
1. Open [http://purl.org/np/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8](http://purl.org/np/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8)
 in your browser
2. Notice that the URL changed to [http://server.nanopubs.lod.labs.vu.nl/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8](http://server.nanopubs.lod.labs.vu.nl/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8).
3. Replace 'server' with 'test-server': [http://test-server.nanopubs.lod.labs.vu.nl/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8](http://test-server.nanopubs.lod.labs.vu.nl/RA71u9tYPd7ZQifE_6hXjqVim6pkweuvjoi-8ehvLvzg8).

> **NB**: `NanopubClient.fetch()` does this for you if `use_test_server=True`.
# Setup instructions

## Install nanopub library
Install using pip:
```
pip install nanopub
```
## Nanopub-java dependency
The ```nanopub``` library currently uses the [nanopub-java](https://github.com/Nanopublication/nanopub-java)
tool for signing and publishing new nanopublications. This is automatically installed by the library.

### Java
If you want to publish nanopublications you need to have the java runtime environment installed,
this might already be installed on your system. You can check this for unix:
```
java --version
```
Or [follow these instructions for windows](http://kb.mit.edu/confluence/pages/viewpage.action?pageId=6750761)

### Installing java
If java is not installed [follow these instructions](https://www.java.com/en/download/help/index_installing.html)

## Setup for users new to python
We recommend using [anaconda](https://www.anaconda.com/products/individual) 
to install python and manage python dependencies

## Setup your profile

To publish to the nanopub server you need to setup your profile (note that you can use
fetch and search functionality without a profile). This allows the nanopub server to identify you.

Run the following interactive command:
```
setup_nanopub_profile
```
This will setup the following:

### Stored profile
A local version of the profile will be stored in the
nanopub user config dir (by default `HOMEDIR/.nanopub/profile.yml`)

### RSA keys
It will add and store RSA keys to sign your nanopublications. By
default they are stored under `HOMEDIR/.nanopub/id_rsa` and `HOMEDIR/.nanopub/id_rsa.pub`.

### ORCID iD
This includes your [ORCID iD](https://orcid.org/) (i.e. https://orcid.org/0000-0000-0000-0000).
If you don't have an ORCID iD yet, you need to [register](https://orcid.org/register). We use
the ORCID iD to automatically add as author to the provenance of any nanopublication you will publish
using this library.

### Introductory nanopublication
We encourage you to make use of `setup_nanopub_profile`'s option 
to publish your profile to the nanopub servers. This links your ORCID iD
to your RSA key, thereby making all your publications linkable to you.
Here is an [example introductory nanopublicaiton](http://purl.org/np/RAy1CYBfBYFd_TFI8Z_jr3taf6fB9u-grqsKyLzTmMvQI).

The link to this nanopublication is also stored in your profile.
# What are nanopublications?
Nanopublications are a formalized and machine-readable way of communicating the smallest possible units
of publishable information. This could be, for example, the outcome of a scientific study or a claim
made by a particular scientist.

Nanopublications are searchable, citable, and contain authorship and attribution
information. The aim is to encourage individual scientific results to be released in a traceable and
interoperable format. As such, nanopublications are an effective [FAIR](https://www.go-fair.org/fair-principles/)
means of communicating scientific claims and results. Read more about them at [http://nanopub.org/](http://nanopub.org/).

## Different elements of a nanopublication
_From [nanopub.org](http://nanopub.org/wordpress/?page_id=65) documentation (2020/12/02)_

![Schematic representation of a nanopub](../img/nanopub.png "Schematic representation of a nanopub")

As can be seen in this image, a nanopublication has three basic elements:

1. Assertion: The assertion is the main content of a nanopublication
in the form of an small atomic unit of information
2. Provenance: This part describes how the assertion above came to be.
This can include the scientific methods that were used to generate the assertion,
for example a reference to the kind of study that was performed and its parameters.
3. Publication Info:  This part contains metadata about the nanopublication as a whole,
such as when and by whom it was created and the license terms for its reuse.
# Using the nanopublication's namespace
In a nanopublication you often want to refer to a concept that is not
defined somewhere on the WWW.
In that case it makes sense to make use of the namespace of the nanopublication itself, 
see for example this assertion that uses `nanopub-uri#timbernerslee` to refer
to the concept Tim Berner's Lee.
```
@prefix sub: <http://purl.org/np/RA_j6TPcnoQJ_XkISjugTgaRsFGLhpbZCC3mE7fXs0REI#> .

sub:assertion {
    sub:timbernerslee a <http://xmlns.com/foaf/0.1/Person> .
}
```
## Using blank nodes 
But how do you make use of the nanopublication's namespace if you do not have
access to the published nanopublication URI yet? We solve that by making use of
blank nodes.

Upon publication, any blank nodes in the rdf graph are replaced with the nanopub's URI, with the blank node name as a
fragment. For example, if the blank node is called 'timbernerslee', that would result in a URI composed of the
nanopub's (base) URI, followed by #timbernslee. We can thus use blank nodes to refer to new concepts, making use of the namespace of the 
to-be-published URI.

An example:

```python
>>> import rdflib
>>> from nanopub import Publication, NanopubClient
>>> 
>>> my_assertion = rdflib.Graph()
>>> 
>>> # We want to introduce a new concept in our publication: Tim Berners Lee
>>> tim = rdflib.BNode('timbernerslee')
>>> 
>>> # We assert that he is a person
>>> my_assertion.add((tim, rdflib.RDF.type, rdflib.FOAF.Person) )
>>> 
>>> # And create a publication object for this assertion
>>> publication = Publication.from_assertion(assertion_rdf=my_assertion)
>>> 
>>> # Let's publish this to the test server
>>> client = NanopubClient(use_test_server=True)
>>> client.publish(publication)
Published to http://purl.org/np/RAdaZsPRcY5usXFKwSBfz9g-HOu-Bo1XmmhQc4g7uESgU
```
View the full nanopublication [here](http://purl.org/np/RAdaZsPRcY5usXFKwSBfz9g-HOu-Bo1XmmhQc4g7uESgU).

As you can see in the assertion, the 'timbernerslee' blank node is replaced with 
a uri in the nanopublication's namespace:
```
@prefix sub: <http://purl.org/np/RAdaZsPRcY5usXFKwSBfz9g-HOu-Bo1XmmhQc4g7uESgU#> .

sub:assertion {
    sub:timbernerslee a <http://xmlns.com/foaf/0.1/Person> .
}
```

## Introducing a concept
You can optionally specify that the Publication introduces a 
particular concept using blank nodes. 
The pubinfo graph will note that this nanopub npx:introduces the concept.
The concept should be a blank node (rdflib.term.BNode), 
and is converted to a URI derived from the nanopub's URI 
with a fragment (#) made from the blank node's name.

An example:
```python
>>> import rdflib
>>> from nanopub import Publication, NanopubClient
>>> 
>>> my_assertion = rdflib.Graph()
>>> 
>>> # We want to introduce a new concept in our publication: Tim Berners Lee
>>> tim = rdflib.BNode('timbernerslee')
>>> 
>>> # We assert that he is a person
>>> my_assertion.add((tim, rdflib.RDF.type, rdflib.FOAF.Person) )
>>> 
>>> # We can create a publication introducing this new concept
>>> publication = Publication.from_assertion(assertion_rdf=my_assertion,
>>>                                          introduces_concept=tim)
>>> 
>>> # Let's publish this to the test server
>>> client = NanopubClient(use_test_server=True)
>>> client.publish(publication)
Published to http://purl.org/np/RAq9gFEgxlOyG9SSDZ5DmBbyGet2z6pkrdWXIVYa6U6qI
Published concept to http://purl.org/np/RAq9gFEgxlOyG9SSDZ5DmBbyGet2z6pkrdWXIVYa6U6qI#timbernerslee
```
Note that `NanopubClient.publish()` now also prints the published concept URI.

View the full nanopublication [here](http://purl.org/np/RAq9gFEgxlOyG9SSDZ5DmBbyGet2z6pkrdWXIVYa6U6qI).

The publication info of the nanopublication denotes that this nanopublication introduces the 'timbernerslee' concept:
```
@prefix npx: <http://purl.org/nanopub/x/> .
@prefix sub: <http://purl.org/np/RAq9gFEgxlOyG9SSDZ5DmBbyGet2z6pkrdWXIVYa6U6qI#> .
@prefix this: <http://purl.org/np/RAq9gFEgxlOyG9SSDZ5DmBbyGet2z6pkrdWXIVYa6U6qI> .

sub:pubInfo {
   this: npx:introduces sub:timbernerslee .
}
```
# Retracting a nanopublication
A nanopublication is persistent, you can never edit nor delete it.
You can however retract a nanopublication.
This is done by publishing a new nanopublication that states that you
retract the original publication. You can use `NanopubClient.retract()`:
```python
>>> from nanopub import NanopubClient
>>> client = NanopubClient(use_test_server=True)
>>> client.retract('http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ')
Published to http://purl.org/np/RAv75Xhhz5jv--Nnu9RDqIGy2xHr74REGC4vtOSxrwX4c
```
View the full retraction nanopublication [here](http://purl.org/np/RAv75Xhhz5jv--Nnu9RDqIGy2xHr74REGC4vtOSxrwX4c).

The assertion states that the researcher (denoted by the ORCID iD from your profile)
retracts the provided nanopublication:
```
@prefix npx: <http://purl.org/nanopub/x/> .
@prefix sub: <http://purl.org/np/RAv75Xhhz5jv--Nnu9RDqIGy2xHr74REGC4vtOSxrwX4c#> .

sub:assertion {
    <https://orcid.org/0000-0000-0000-0000> npx:retracts <http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ> .
}
```
By default nanopublications that have a valid retraction do not show up in search results.
A valid retraction is a retraction that is signed with the same public key as 
the nanopublication that it retracts.

## Retracting a nanopublication that is not yours
By default we do not retract nanopublications that are not yours (i.e. signed with another public key). 
If you try to do this it will trigger an AssertionError:
```python
>>> from nanopub import NanopubClient
>>> client = NanopubClient(use_test_server=True)
>>> not_my_nanopub_uri = 'http://purl.org/np/RAr6rs7o8Sr5OGCs0127ah37DYUvgiWzjOuCvV-OSusAk'
>>> client.retract(not_my_nanopub_uri)
---------------------------------------------------------------------------
AssertionError                            Traceback (most recent call last)
<ipython-input-30-7141d9e82fbc> in <module>
      1 not_my_nanopub_uri = 'http://purl.org/np/RAr6rs7o8Sr5OGCs0127ah37DYUvgiWzjOuCvV-OSusAk'
----> 2 client.retract(not_my_nanopub_uri)

~/projects/fair-workflows/nanopub/nanopub/client.py in retract(self, uri, force)
    265         """
    266         if not force:
--> 267             self._check_public_keys_match(uri)
    268         assertion_rdf = rdflib.Graph()
    269         orcid_id = profile.get_orcid_id()

~/projects/fair-workflows/nanopub/nanopub/client.py in _check_public_keys_match(self, uri)
    245                               f'this one: {their_public_key}')
    246             if their_public_key != profile.get_public_key():
--> 247                 raise AssertionError('The public key in your profile does not match the public key'
    248                                      'that the publication that you want to retract is signed '
    249                                      'with. Use force=True to force retraction anyway.')

AssertionError: The public key in your profile does not match the public keythat the publication that you want to retract is signed with. Use force=True to force retraction anyway.
```
We can use force=True to override this behavior:
```python
client.retract(not_my_nanopub_uri, force=True)
```

## Find retractions of a given nanopublication
You can find out whether a given publication is retracted 
and what the nanopublications are that retract it using `NanopubClient.find_retractions_of`:
```python
>>> from nanopub import NanopubClient
>>> client = NanopubClient(use_test_server=True)
>>> # This URI has 1 retraction:
>>> client.find_retractions_of('http://purl.org/np/RAirauh-vy5f7UJEMTm08C5bh5pnWD-abb-qk3fPYWCzc')
['http://purl.org/np/RADjlGIB8Vqt7NbG1kqzw-4aIV_k7nyIRirMhPKEYVSlc']
>>> # This URI has no retractions
>>> client.find_retractions_of('http://purl.org/np/RAeMfoa6I05zoUmK6sRypCIy3wIpTgS8gkum7vdfOamn8')
[]
```
# Publishing nanopublications
The `nanopub` library provides an intuitive API that makes publishing nanopublications much easier. 
The rationale is that you often do not want to worry about the details of composing 
the RDF that is often the same in each nanopublication. Instead you should focus on the 
content of your nanopublication: the assertion.

## Prerequisits for publishing
Before you can publish you should [setup your profile](../getting-started/setup)

## Quickly publishing nanopublications using `claim`
You can quickly publish a nanopublicaiton with a single simple statement using the `claim` method:
```python
>>> from nanopub import NanopubClient

>>> # Create the client (we use use_test_server=True to point to the test server)
>>> client = NanopubClient(use_test_server=True)

>>> # Publish a simple statement to the server
>>> client.claim('All cats are gray')
Published to http://purl.org/np/RA47eJP2UBJCWuJ324c6Qw0OwtCb8wCrprwSk39am7xck
```
View the resulting nanopublication [here](http://purl.org/np/RA47eJP2UBJCWuJ324c6Qw0OwtCb8wCrprwSk39am7xck).

The generated RDF makes use of the Hypotheses and Claims Ontology ([HYCL](http://purl.org/petapico/o/hycl))

This is the assertion part of the nanopublication, denoting the statement:
```
@prefix hycl: <http://purl.org/petapico/o/hycl#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix sub: <http://purl.org/np/RA47eJP2UBJCWuJ324c6Qw0OwtCb8wCrprwSk39am7xck#> .

sub:assertion {
    sub:mystatement a hycl:Statement ;
        rdfs:label "All cats are gray" .
}
```

The provenance part of the nanopublication denotes that the ORCID iD from the profile claimed the
statement:
```
@prefix hycl: <http://purl.org/petapico/o/hycl#> .
@prefix prov: <http://www.w3.org/ns/prov#> .
@prefix sub: <http://purl.org/np/RA47eJP2UBJCWuJ324c6Qw0OwtCb8wCrprwSk39am7xck#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

sub:provenance {
    sub:assertion prov:generatedAtTime "2020-12-01T10:39:09.427920"^^xsd:dateTime ;
        prov:wasAttributedTo <https://orcid.org/0000-0000-0000-0000> .

    <https://orcid.org/0000-0000-0000-0000> hycl:claims sub:mystatement .
}
```
 
## A simple recipe for publishing RDF triples
You can use `Publication` objects to easily publish nanopublications with your assertion 
(think of the assertion as the content of your nanopublication).

This is a 3-step recipe that works for most cases:
 1) Construct a desired assertion using [`rdflib`](https://rdflib.readthedocs.io/en/stable/).
 2) Make a `Publication` object using the assertion, making use of `Publication.from_assertion()`.
 3) Publish the `Publication` object using `NanopubClient.publish()`.
 
Here is a minimal example:
```python
>>> import rdflib
>>> from nanopub import NanopubClient, Publication
>>> 
>>> # Create the client (we use use_test_server=True to point to the test server)
>>> client = NanopubClient(use_test_server=True)
>>> 
>>> # 1. construct a desired assertion (a graph of RDF triples) using rdflib
>>> my_assertion = rdflib.Graph()
>>> my_assertion.add((rdflib.URIRef('www.example.org/timbernerslee'),
>>>                   rdflib.RDF.type,
>>>                   rdflib.FOAF.Person))
>>> 
>>> # 2. Make a Publication object with this assertion
>>> publication = Publication.from_assertion(assertion_rdf=my_assertion)
>>> 
>>> # 3. Publish the Publication object.
>>> publication_info = client.publish(publication)
Published to http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ
```
View the resulting nanopublication [here](http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ).

This is the resulting assertion part of the nanopublication:
```
@prefix sub: <http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ#> .

sub:assertion {
    <https://www.example.org/timbernerslee> a <http://xmlns.com/foaf/0.1/Person> .
}
```

The library automatically adds relevant RDF triples for the provenance part of the nanopublication:
```
@prefix prov: <http://www.w3.org/ns/prov#> .
@prefix sub: <http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

sub:provenance {
    sub:assertion prov:generatedAtTime "2020-12-01T10:44:32.367084"^^xsd:dateTime .
}
```
as well as for the publication info part of the nanopublication:
```
@prefix npx: <http://purl.org/nanopub/x/> .
@prefix prov: <http://www.w3.org/ns/prov#> .
@prefix sub: <http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ#> .
@prefix this: <http://purl.org/np/RAfk_zBYDerxd6ipfv8fAcQHEzgZcVylMTEkiLlMzsgwQ> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

sub:pubInfo {
    sub:sig npx:hasAlgorithm "RSA" ;
        npx:hasPublicKey "MIGfMA0GCSqGSIb3DQEBAQUAA4GNADCBiQKBgQCmso7vmRO/Cp4Pt0RkJJkV5qfc1WFYU/jMtkdxxb5+lfIVXNV97XQnM1Tj4fkb/W6jkP6fHl8mj8Q7hl7VgUnQ6I+B7cMGpxW9Z8Br+JNx8DPMMt08VCH5+JMENPRKl91r7rF/YPWCAgL9eqXSixCNMNAj5RBmMTQoPuRkpgmt1wIDAQAB" ;
        npx:hasSignature "aPZMJ3Md6X1PHYvXJiNoRUni9+1oS9faCfiPRRCrj4K/uZPN0J/znjxGuCUxoZRJ4b4RfSxmHFGRKfCFusJX+7Y3xuxYx4GYHzYhBciK7T5pO02V4w6sdwHLKd5E+Wcl0PTr2t3lEjq6yzY98wEXlZLAbaRDBJvzpg5xORifQDw=" ;
        npx:hasSignatureTarget this: .

    this: prov:generatedAtTime "2020-12-01T10:44:32.367084"^^xsd:dateTime ;
        prov:wasAttributedTo <https://orcid.org/0000-0000-0000-0000> .
}
```
# Setting publication info and provenance
Here we show how you can control the publication info and provenance parts 
of the nanopublication.

## Specifying where the nanopublication is derived from
You can specify that the nanopub's assertion is derived from another URI (such as an existing nanopublication):
```python
import rdflib
from nanopub import Publication

my_assertion = rdflib.Graph()
my_assertion.add((rdflib.term.BNode('timbernserslee'), rdflib.RDF.type, rdflib.FOAF.Person))


publication = Publication.from_assertion(
    assertion_rdf=my_assertion,
    derived_from=rdflib.URIRef('http://www.example.org/another-nanopublication'))
```
Note that ```derived_from``` may also be passed a list of URIs.

The provenance part of the publication will denote:
```
@prefix sub: <http://purl.org/nanopub/temp/mynanopub#> .
@prefix prov: <http://www.w3.org/ns/prov#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

:provenance {
    sub:assertion prov:wasDerivedFrom <http://www.example.org/another-nanopublication> .
}
```

## Attributing the assertion to someone
You can attribute the assertion to someone by specifying the `assertion_attributed_to` argument:
```python
import rdflib
from nanopub import Publication

my_assertion = rdflib.Graph()
my_assertion.add((rdflib.term.BNode('timbernserslee'), rdflib.RDF.type, rdflib.FOAF.Person))

publication = Publication.from_assertion(
    assertion_rdf=my_assertion,
    assertion_attributed_to=rdflib.URIRef('https://orcid.org/0000-0000-0000-0000'))
```

The provenance part of the publication will denote:
```
@prefix : <http://purl.org/nanopub/temp/mynanopub#> .
@prefix prov: <http://www.w3.org/ns/prov#> .

:provenance {
    :assertion prov:wasAttributedTo <https://orcid.org/0000-0000-0000-0000> .
}
```
Note: Often the assertion should be attributed to yourself.
Instead of passing your ORCID iD to `assertion_attributed_to`,
you can easily tell nanopub to attribute the assertion to
the ORCID iD in your profile by setting `attribute_assertion_to_profile=True`.

## Specifying custom provenance triples
You can add your own triples to the provenance graph of the nanopublication
by passing them in an `rdflib.Graph` object to the `provenance_rdf` argument:
```python
import rdflib
from nanopub import namespaces, Publication

my_assertion = rdflib.Graph()
my_assertion.add((rdflib.term.BNode('timbernserslee'), rdflib.RDF.type, rdflib.FOAF.Person))

provenance_rdf = rdflib.Graph()
provenance_rdf = provenance_rdf.add((rdflib.term.BNode('timbernserslee'),
                                     namespaces.PROV.actedOnBehalfOf,
                                     rdflib.term.BNode('markzuckerberg')))
publication = Publication.from_assertion(assertion_rdf=my_assertion,
                                         provenance_rdf=provenance_rdf)
```

## Specifying custom publication info triples
You can add your own triples to the publication info graph of the nanopublication
by passing them in an `rdflib.Graph` object to the `pubinfo_rdf` argument:
```python
import rdflib
from nanopub import namespaces, Publication

my_assertion = rdflib.Graph()
my_assertion.add((rdflib.term.BNode('timbernserslee'), rdflib.RDF.type, rdflib.FOAF.Person))

pubinfo_rdf = rdflib.Graph()
pubinfo_rdf = pubinfo_rdf.add((rdflib.term.BNode('activity'),
                               rdflib.RDF.type,
                               namespaces.PROV.Activity))
publication = Publication.from_assertion(assertion_rdf=my_assertion,
                                         pubinfo_rdf=pubinfo_rdf)
```
