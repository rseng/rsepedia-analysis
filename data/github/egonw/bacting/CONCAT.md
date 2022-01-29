# Contributions

There are various ways to contribute to Bacting, and this does not only involve development, but can also be
checking of an existing Bioclipse scripts can already be run with Bacting (the Bioclipse 2.6 API is
[only partially ported](https://github.com/egonw/bacting/projects/2)).

## User experience

We love to hear from you about the Bioclipse scripts you have written in the past. We like to learn
from this to see what functionality needs to be reported. You can use [this issue template](https://github.com/egonw/bacting/issues/new?assignees=&labels=enhancement&template=feature-request--bioclipse-api-method.md&title=)
for that.

## Code development

Developing Bioclipse is not trivial. Developing Bacting is easier but still has a learning curve.
It is recommended to do all code development in branches, allowing others to suggest improvements
and generally do peer review. If you fork the repository, it is suggested to keep your `master`
branch identical to the upstream version.

### Porting Bioclipse managers

There are still many unported managers. Some are easy to port, some are hard. Some managers focus on
the Bioclipse graphical user interface, which in Bacting is not available. Porting such a manager will
be a lot harder. Other managers cannot be ported yet because the used libraries are not yet available
from Maven Central. However, there is plenty of work that can still be done.

Once you decided it is needed to port a manager, one easy way
to start a new manager is to just copy/paste the full folder of a simple manager, like that of OPSIN
found [here](https://github.com/egonw/bacting/tree/master/managers-cheminfo/net.bioclipse.managers.opsin).
The copied files can then be updated with a text editor to have unique names.

With the folder structure in place, you can start copy/pasting manager content from the Bioclipse project.
Make sure to always:

* keep copyright and license information for all copied Bioclipse source code
* ideally, separate copy/pasting the code and making any updates
* the API should match that of the manager interface, not of the Bioclipse implementation

The latter means you should be aware that Bacting does not have:

* a `IProgressMonitor` and the Bacting API should exclude those method parameters
* use `String file` instead of `IFile file`

#### Testing

Testing for ported managers is welcome buit not required. Bioclipse itself has a test suite.

### Developing a new manager

Bacting is not limited to existing managers and new managers are welcome, wrapping around a yet unused
bioinformatics of cheminformatics Java library. It is recommended to use the
[Issue tracker](https://github.com/egonw/bacting/issues) first to discuss the plans and to make sure
the Java library is indeed not used yet.

New managers can be developed in this repository, but do not have to be. The section
about starting with the folder structure of a simple manager in "Porting Bioclipse managers"
can be used here too.

#### Testing

Testing for new managers is recommended. 
# Bacting

[![License](https://img.shields.io/badge/License-EPL%201.0-red.svg)](https://opensource.org/licenses/EPL-1.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2638709.svg)](https://doi.org/10.5281/zenodo.2638709)
[![build](https://github.com/egonw/bacting/workflows/build/badge.svg)](https://github.com/egonw/bacting/actions?query=workflow%3Abuild)
[![Maven Central](https://img.shields.io/maven-central/v/io.github.egonw.bacting/bacting.svg?label=Maven%20Central)](https://search.maven.org/search?q=g:%22io.github.egonw.bacting%22%20AND%20a:%22bacting%22)
[![codecov](https://codecov.io/gh/egonw/bacting/branch/master/graph/badge.svg?token=E1NGWVWL04)](https://codecov.io/gh/egonw/bacting)

Bacting := acting as the Bioclipse TNG (The Next Generation)

Bacting is an open-source platform for chemo- and bioinformatics based on [Bioclipse](https://scholia.toolforge.org/topic/Q1769726)
that defines a number of common domain objects and wraps common functionality, providing a toolkit independent, scriptable solution to
handle data from the life sciences. Like Bioclipse, Bacting is written in the Java language, making use in Java-derived
languages like [Groovy](https://en.wikipedia.org/wiki/Apache_Groovy) easy, but also accessible to Python. Deposition of the Bacting package on
[Maven Central](https://search.maven.org/search?q=g:%22io.github.egonw.bacting%22%20AND%20a:%22bacting%22) allows it
to be easily used in Groovy scripts with `@Grab` instructions.

## How to cite

If you use this software, please cite the article in JOSS:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02558/status.svg)](https://doi.org/10.21105/joss.02558)

# Install

For the below use cases, Bacting is actually installed on demand. In Groovy this is done with
`@Grab` and in Python with `from pybacting import cdk` (see [pybacting](https://github.com/cthoyt/pybacting))
or `scyjava.config`. This section explains how Bacting can be installed from
the source code.

## From the source code

First, you need a working [Maven installation](https://www.google.nl/search?q=install+maven) and the code is tested with
Java 11, and 14, and can be installed with:

```shell
mvn clean install -Dgpg.skip -Dmaven.javadoc.skip=true
```

## Making releases

Before making a release, update the version number in this `README.md` and in `CITATION.cff`.

Releases are created by the release manager and requires permission to submit the release to Maven Central
(using an approved Sonatype ([oss.sonatype.org](http://oss.sonatype.org/)) account).
If these requirements are fulfilled then the following two commands to the job:

```shell
export MAVEN_OPTS="--add-opens=java.base/java.util=ALL-UNNAMED --add-opens=java.base/java.lang.reflect=ALL-UNNAMED --add-opens=java.base/java.text=ALL-UNNAMED --add-opens=java.desktop/java.awt.font=ALL-UNNAMED"
mvn versions:set -DnewVersion=0.0.31
mvn deploy -P release
```

### Making snapshots

```shell
mvn versions:set -DnewVersion=0.0.31-SNAPSHOT
mvn deploy
```

### Updating the JavaDoc

After the release is avaiable on Maven Central, the [JavaDoc](https://egonw.github.io/bacting-api/)
needs to be updated. The JavaDoc is generated with the below command, and the results are stored
in the `target/site/apidocs/` folder:

```shell
mvn clean javadoc:javadoc javadoc:aggregate
```

That created content needs to be copied into the `docs/` folder of
[this git repository](https://github.com/egonw/bacting-api/).

# Usage

## Groovy

It can be used in [Groovy](https://en.wikipedia.org/wiki/Apache_Groovy) by including the
Bacting managers you need. The following example tells Groovy to download the `CDKManager`
and instantiate it for the given workspace location (as it if was running in Bioclipse
itself), and then converts a [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
string to a Bioclipse `IMolecule` data object:

```groovy
@Grab(group='io.github.egonw.bacting', module='managers-cdk', version='0.0.31')

workspaceRoot = "."
def cdk = new net.bioclipse.managers.CDKManager(workspaceRoot);

println cdk.fromSMILES("COC")
```

## Python

Bacting can also be used in Python 3.7 with [pybacting](https://github.com/cthoyt/pybacting) and
[scyjava](https://github.com/scijava/scyjava).

### Pybacting

Pybacting can be installed with `pip install pybacting` (or `pip3`, depending on your platform).
The above code example looks like:

```python
from pybacting import cdk

print(cdk.fromSMILES("COC"))
```

Pybacting uses a specific Bacting version, so check the [website](https://github.com/cthoyt/pybacting)
to see which Bacting version you are using.

### Scyjava

Scyjava can be installed with `pip install scyjava` (or `pip3`, depending on your platform).
The code example looks like:

```python
from scyjava import config, jimport
config.add_endpoints('io.github.egonw.bacting:managers-cdk:0.0.31')

workspaceRoot = "."
cdkClass = jimport("net.bioclipse.managers.CDKManager")
cdk = cdkClass(workspaceRoot)

print(cdk.fromSMILES("COC"))
```

# Code examples

Full code examples can be found in the following sources:

* Open Notebooks for Wikidata, including:
  * [script that compares a SMILES string with Wikidata, and creates QuickStatements for missing information](https://github.com/egonw/ons-wikidata/blob/master/Wikidata/createWDitemsFromSMILES.groovy)
  * [script that reads melting points from an Excel spreadsheet to enter into Wikidata](https://github.com/egonw/ons-wikidata/blob/master/MeltingPoints/createQuickStatements.groovy)
* Open Notebooks for WikiPathays:
  * [script to recognize IUPAC names in WikiPathways](https://github.com/egonw/ons-wikipathways/blob/master/WikiPathways/getLabelsWithIUPACNames.groovy)
* Some examples from [A lot of Bioclipse Scripting Language examples](https://bioclipse.github.io/bioclipse.scripting/):
  * [FullPathWikiPathways.groovy](https://bioclipse.github.io/bioclipse.scripting/code/FullPathWikiPathways.code.html)
  * [XMLIsWellFormed.groovy](https://bioclipse.github.io/bioclipse.scripting/code/XMLIsWellFormed.code.html)
  * [XMLListNamespaces.groovy](https://bioclipse.github.io/bioclipse.scripting/code/XMLListNamespaces.code.html)
  * [PerceiveCDKAtomTypes.groovy](https://bioclipse.github.io/bioclipse.scripting/code/PerceiveCDKAtomTypes.code.html)
  * [InChIKeyGenerate.groovy](https://bioclipse.github.io/bioclipse.scripting/code/InChIKeyGenerate.code.html)
  * [ParseIUPACName.groovy](https://bioclipse.github.io/bioclipse.scripting/code/ParseIUPACName.code.html)
  * [ChemSpiderResolve.groovy](https://bioclipse.github.io/bioclipse.scripting/code/ChemSpiderResolve.code.html)

## API Coverage

For the time being, the coverage of the original API is *incomplete*.
Particularly, manager functionality around graphical UX
in the original Bioclipse may never be implemented. Each Bacting release will implement more APIs and
the release notes will mention which managers and which methods have been added.
An overview of the supports APIs can be found in [this overview](https://github.com/egonw/bacting/projects/2).

For a description of the API, I refer to the book
[A lot of Bioclipse Scripting Language examples](https://bioclipse.github.io/bioclipse.scripting/) that
Jonathan and I compiled. However, a [JavaDoc API](https://egonw.github.io/bacting-api/) is also available.

All Bacting scripts will be backwards compatible with Bioclipse. If you want to install Bioclipse
and see its wonderful UX in actions, [download Bioclipse 2.6.2 here](https://sourceforge.net/projects/bioclipse/files/bioclipse2/bioclipse2.6.2/).

## Using SNAPSHOT versions

You may need to occassionally delete the
modules cached by Groovy, by doing something like, to remove earlier SNAPSHOT versions:

```shell
\rm -Rf ~/.groovy/grapes/io.github.egonw.bacting/
```

# Copyright and authors

Code in this repository contains mostly code that originated from Bioclipse
and the headers of the individual source code files describe who contributed to that code of that class, but unfortunately this code
ownership is not always clear. I refer to the various [Bioclipse code repositories](https://github.com/bioclipse)
for the git history for detailed information.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

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
reported by contacting the project team at egon.willighagen+bacting@gmail.com. All
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

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Some porting basics

Porting [Bioclipse](https://github.com/bioclipse) functionality is not trivial, but once you picked up the pattern, not too hard. Some key differences:

- Bacting does not have `IFile`, but uses `String` instead (matching the API of Bioclipse, instead of the implementations)
- Bacting does not have `IProgressMonitor` (just remove it)
- code should be copied as exact as possible
- copyright ownership must be kept; copying is non-trvial so adding yourself as copyright owner is acceptable
- the test method must be copied too, and [code coverage](https://app.codecov.io/gh/egonw/bacting) should be taken into account
---
title: 'Bacting: a next generation, command line version of Bioclipse'
tags:
  - bioinformatics
  - cheminformatics
  - Bioclipse
authors:
  - name: Egon Willighagen
    orcid: 0000-0001-7542-0286
    affiliation: 1
affiliations:
 - name: Dept of Bioinformatics - BiGCaT, NUTRIM, Maastricht University
   index: 1
date: 14 June 2021
bibliography: paper.bib
---

# Summary

Bioclipse ([https://github.com/bioclipse](https://github.com/bioclipse)) was originally
developed as an interactive user interface (UI) based on Eclipse for research in the fields
of biology and chemistry [@bioclipse1]. It was later extended with scripting
functionality and scripts could be written in JavaScript, Python, and Groovy [@bioclipse2].
An innovative aspect of the second Bioclipse version was that Bioclipse plugins could inject
domain specific functionality into the scripting language. This was done using
OSGi ([https://www.osgi.org/](https://www.osgi.org/))
and the Spring Framework ([https://spring.io/](https://spring.io/)),
making so-called *managers* accessible in scripts. However, there have not been any
recent Bioclipse releases. Bacting is a next generation,
command line version of Bioclipse, that is more easily updated, built, released, and used. A subset
of the original functionality is available, and some managers have already been updated to
use more recent versions of dependencies.

# Statement of Need

While Bioclipse has served our research for many years, a number of limitations has made
this increasingly hard. For example, the dependency of Bioclipse on the Eclipse UI requires the scripts to be run inside
a running Bioclipse application. This makes repeatedly running of a script needlessly hard
and use in continuous integration systems or use on computing platforms impossible. A second
problem was that the build and release system of Bioclipse was complex, making it hard for
others to repeat creating new releases. This is reflected in the lack of recent releases
and complicates the process for external developers wishing to make patches.

These needs triggered a next generation design of Bioclipse: 1. the managers providing
the domain-specific functionality would need to be usable on the command line; 2. building
the Bioclipse managers should be possible on the command line, ideally with continuous build
systems; 3. Bacting should be easy to install and reuse.

# Implementation

To keep the option open to backport new functionality to Bioclipse, the API is copied as
precisely as possible. However, there are some differences. For example, there is only
a single manager class, and no longer interfaces for both the scripting language or for
running Bioclipse user interface. This means that the *IFile* to *String*
translations in the API do not exist in Bioclipse. Furthermore, there are currently
no progress monitors. That said, the source code implementing the method is otherwise
identical and easily translated back to the original Bioclipse source code.

This is done by separating the Bioclipse code from the Bacting manager implementations.
The latter is mainly described in this paper, and the former is found as much more
stable code in the GitHub [https://github.com/egonw/bacting-bioclipse](https://github.com/egonw/bacting-bioclipse)
repository. The code is identical to the original Bioclipse code, but mavenized in
this repository, allowing deployment on Maven Central.

The Bacting manager are found in the [https://github.com/egonw/bacting](https://github.com/egonw/bacting)
repository and while the managers in this repository share most of the code with the original
Bioclipse implementations, they are still considered new implementations and therefore
are tested using JUnit. A second important difference is that Bioclipse documentation was
found on the manager interfaces, but in Bacting the JavaDoc is found is found in the
implementations of the managers. A final difference is how the managers are used:
because they are not injected into the scripting language, each manager needs to be
created manually, which requires one extra line of code for each manager.

## Continuous integration and releases

Bacting is hosted on GitHub and takes advantage of the integrations with Zenodo for automatic
archiving of releases (see [https://github.com/egonw/bacting/releases](https://github.com/egonw/bacting/releases))
and with GitHub Actions for continuous integration (see [https://github.com/egonw/bacting/actions](https://github.com/egonw/bacting/actions)).
Maven is used as a build system and automatically downloads the dependencies when compiling the source code.
GitHub Actions compiles the source code regularly with Java 8, 11, and 14. During the process the JUnit 5 unit
tests are run and the compilation aborted when there are testing failures. The extend to which
the tests execute code in the managers is tested with JaCoCo ([https://www.jacoco.org/jacoco/](https://www.jacoco.org/jacoco/))
and reported online with Codecov at ([https://codecov.io/gh/egonw/bacting](https://codecov.io/gh/egonw/bacting)).

Releases are made at irregular intervals, but often triggered by downstream uses that need additional
Bioclipse functionality to be ported. Releases are created
with the `mvn release:prepare` and `mvn release:perform` process that tags the commit, updates the
version numbers, and uploads the release to Maven Central. Second, a changelog is written for the
GitHub releases page, which triggers the archiving on Zenodo (see
[https://doi.org/10.5281/zenodo.2638709](https://doi.org/10.5281/zenodo.2638709)). Finally, at that
moment the JavaDoc is also generated and uploaded to another GitHub repository
(see [https://github.com/egonw/bacting-api](https://github.com/egonw/bacting-api))
making it available online with GitHub pages at [https://egonw.github.io/bacting-api/](https://egonw.github.io/bacting-api/).

## Updated dependencies of managers

The *cdk* manager wrapping Chemistry Development Kit functionality was updated to
version 2.3, released in 2017 [@Mayfield2019; @Willighagen2017]. The *opsin* manager was
updated to use OPSIN version 2.5.0, released in 2020 [@Lowe2011]. The *bridgedb*
manager was updated to BridgeDb version 2.3.10, released in 2020 [@Brenninkmeijer2020; @vanIersel2010].

# Ported Functionality

Bioclipse has a long list of managers and so far only a subset has been ported, which is briefly described in this table:

| Bacting Manager      | Functionality                                                                        |
| -------------------- | ------------------------------------------------------------------------------------ |
| bioclipse            | Bioclipse manager with common functionality                                          |
| ui                   | Bioclipse manager with user interface functionality                                  |
| report               | Manager that provides an API to create HTML reports                                  |
| cdk                  | Chemistry Development Kit for cheminformatics functionality [@Willighagen2017]       |
| inchi                | Methods for generating and validating InChIs and InChIKeys [@Spjuth2013]             |
| pubchem              | Methods to interact with the PubChem databases                                       |
| chemspider           | Methods to interact with the Chemspider databases                                    |
| rdf                  | Resource Description Framework (RDF) functionality, using Apache Jena                |
| opsin                | Access to the OPSIN library for parsing IUPAC names [@Lowe2011]                      |
| bridgedb             | Access to the BridgeDb library for identifier mapping [@vanIersel2010]               |
| biojava              | Access to the Biojava library for sequence functionality [@Holland2008]              |

The functionality of the Bioclipse managers is partly documented in the 
[A lot of Bioclipse Scripting Language examples](https://bioclipse.github.io/bioclipse.scripting/) booklet,
of which several scripts are available as Bacting examples. For example, the
[FullPathWikiPathways.groovy](https://bioclipse.github.io/bioclipse.scripting/code/FullPathWikiPathways.code.html)
page from this booklet shows both the Bioclipse version of the script as well as the Bacting version.

## Grabbing Bacting from Groovy

Use of Bacting in the Groovy language takes advantage of the fact that it is available from Maven Central,
allowing `@Grab` to be used to dynamically download the code as in this example for the *cdk* manager:

```groovy
@Grab(
  group='io.github.egonw.bacting',
  module='managers-cdk', version='0.0.15'
)

def cdk = new net.bioclipse.managers.CDKManager(".");

println cdk.fromSMILES("COC")
```

Similarly, Bacting can be used in Python using [scyjava](https://pypi.org/project/scyjava/).

# Use cases

Bioclipse scripts have been in use in our group in various research lines to automate repetitive work.
Various scripts have now been ported to Bacting and several are now available as open notebook science
repositories at [https://github.com/egonw/ons-wikidata](https://github.com/egonw/ons-wikidata),
[https://github.com/egonw/ons-chebi](https://github.com/egonw/ons-chebi), and
[https://github.com/egonw/ons-wikipathways](https://github.com/egonw/ons-wikipathways). The scripts in these repositories are
used to populate Wikidata with chemical structures to support the Scholia
project [@Nielsen2017; @Willighagen2018], the WikiPathways project [@Slenter2018], and
feed additional metabolite identifiers into Wikidata for creation of BridgeDb identifier mapping databases
in an implementation study of the ELIXIR Metabolomics Community [@Willighagen2020; @vanRijswijk2017].
Furthermore, Bacting is used to populate Wikidata with OECD Testing Guidelines in the [NanoCommons](https://www.nanocommons.eu/)
project and extend the eNanoMapper ontology (see [https://github.com/egonw/ons-wikidata/blob/master/OECD/convertToOWL.groovy](https://github.com/egonw/ons-wikidata/blob/master/OECD/convertToOWL.groovy)) [@Hastings2015],
to generate RDF for a public data set in the [NanoSolveIT](https://nanosolveit.eu/) project (see [https://github.com/NanoSolveIT/10.1021-acsnano.8b07562](https://github.com/NanoSolveIT/10.1021-acsnano.8b07562)) [@Afantitis2020], to create a booklet with data about the SARS-CoV-2
and related coronavirusses (see [https://github.com/egonw/SARS-CoV-2-Queries](https://github.com/egonw/SARS-CoV-2-Queries)), and to support 
Various of these use cases are ongoing and are not yet published, which is planned.

# Acknowledgements

We acknowledge the contributions of the Bioclipse developers which have been
ported here into the Bacting software.

# References
# Patch description

A concise description of what this patch does. Where appropriate, mention issue reports, etc.

Changes proposed in this pull request:

 - 
 - 
 - 

# Pull request checks

When reusing Bioclipse code, please ensure:

 - [ ] does all code have proper copyright and license info?
 - [ ] all reused code has an OSI-approved license?
 - [ ] all new code has an Eclipse Public License (1.0 or later and GPL exception)?

@egonw
---
name: 'Bug report: unexpected output of Bacting API method'
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

## Bacting API method with unexpected output

A short description of the unexpected output. Please include the code you used, and the output in Bacting and what you expected or got with Bioclipse below:

### Expected Output

```
The expected output
```

### Actual Output

```
The actual output with Bacting
```

## Additional context

Add any other context or screenshots about the feature request here.
---
name: 'Feature request: Bioclipse API method'
about: Suggest a Bioclipse API method to be ported
title: ''
labels: enhancement
assignees: ''

---

## Bioclipse API method to be ported

A short description of the Bioclipse API method you like to see ported. Ideally, provide as much context as possible, such as journal articles describing or using the manager, links to source code, etc.

## Example use

A short code example where the API method is used.

## Additional context

Add any other context or screenshots about the feature request here.
