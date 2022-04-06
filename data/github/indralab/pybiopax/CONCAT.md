PyBioPAX: A python implementation of the BioPAX object model
------------------------------------------------------------
[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![DOI](https://zenodo.org/badge/261255657.svg)](https://zenodo.org/badge/latestdoi/261255657)
[![Build](https://github.com/indralab/pybiopax/workflows/Tests/badge.svg)](https://github.com/indralab/pybiopax/actions)
[![Documentation](https://readthedocs.org/projects/pybiopax/badge/?version=latest)](https://pybiopax.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/pybiopax.svg)](https://badge.fury.io/py/pybiopax)
[![Python 3](https://img.shields.io/pypi/pyversions/pybiopax.svg)](https://www.python.org/downloads/release/python-357/)

PyBioPAX implements the BioPAX level 3 object model (http://www.biopax.org/release/biopax-level3-documentation.pdf) as a set of
Python classes. It exposes API functions to read OWL files into this
object model, and to dump OWL files from this object model.
This allows for the processing and creation of BioPAX models natively in
Python.

Installation
------------
PyBioPAX can be installed from PyPI as a package:

```bash
$ pip install pybiopax
```

Usage
-----
Reading an OWL file into a BioPaxModel object:

```python
import pybiopax
model = pybiopax.model_from_owl_file('test.owl')
```


Writing a BioPaxModel into an OWL file:

```python
import pybiopax
pybiopax.model_to_owl_file(model, 'test.owl')
```

Querying Pathway Commons to get a BioPaxModel object:

```python
import pybiopax
model = pybiopax.model_from_pc_query('pathsfromto', ['MAP2K1'], ['MAPK1'])
```

Working with the elements of the Python object model:

```python
import pybiopax
model = pybiopax.model_from_pc_query('pathsfromto', ['MAP2K1'], ['MAPK1'])

# Each BioPaxModel instance has an objects attribute which is a dict
# whose keys are object URIs (strings) and values are BioPaxObject instances.
assert isinstance(model.objects, dict)
assert all(isinstance(obj, pybiopax.biopax.BioPaxObject)
           for obj in model.objects.values())

# Let's look at a specific object
bcr = model.objects['BiochemicalReaction_4f689747397d98089c551022a3ae2d88']

# This is a BiochemicalReaction which has a left and a right side. All list/set
# types per the BioPAX specification are represented as lists in the Python
# object model
# Both left and right consist of a single protein
left = bcr.left[0]
assert isinstance(left, pybiopax.biopax.Protein)
assert left.display_name == 'ERK1-2'
right = bcr.right[0]
assert isinstance(right, pybiopax.biopax.Protein)
assert right.display_name == 'ERK1-2-active'
```

We can also use the `pybiopax.paths` module to construct iterators over
objects based on a string specification from a given starting point.
Continuing from the block of code above, we take the BiochemicalReaction
`bcr` and link to reactants on its left hand side, then linking to their
entity references, and finally linking back to all the physical entities
that those are references of.

```python
from pybiopax.paths import find_objects

erks = find_objects(bcr, 'left/entity_reference/entity_reference_of')
```

Contribution and support
------------------------
To contribute to the code, please submit a pull request after
reading the [contribution guidelines](https://github.com/indralab/pybiopax/blob/master/CONTRIBUTING.md).
To report bugs and issues, or ask questions related to PyBioPAX, please
submit an [issue](https://github.com/indralab/pybiopax/issues).

Funding
-------
Development of this software was supported by the Defense Advanced Research
Projects Agency under award W911NF-15-1-0544 and the National Cancer Institute
under award U54-CA225088.
Git fork/PR workflow
--------------------
This repository uses the [forking model](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow)
for collaboration. In this model,
each developer forks the main (indralab/pybiopax) repository, pushes code
only to branches in their own fork, and then submits pull requests to
`indralab`. After cloning your own fork of `pybiopax`, you should add `indralab`
as a remote to be able to track the latest changes by doing

```
git remote add indralab git@github.com:indralab/pybiopax.git
```

When a PR is submitted from a branch, any further changes can be pushed
to that same branch even after the PR has been opened; those changes
are automatically appended to the PR. Please always check the box on Github
allowing the maintainers of the repo to make changes to the PR.

In addition, as a convention, we only merge PRs whose branches are rebased
on top of the latest indralab/master. This means that instead of merging
indralab/master into your own branch to resolve conflicts, you should always
rebase on top of indralab/master and force push your branches if
needed (you can do this even if a PR from that branch is already open).
Consistent with this convention, in general, you should not use `git pull`
to update your local fork. Rather, use `git fetch --all`,
`git merge --ff-only`, `git rebase` or `git reset --hard` as needed to
get to the desired state. PRs are always merged using a separate merge commit,
ensuring that merges clearly correspond to the inclusion of a specific
feature or the fixing of a specific issue. In some cases, for instance, if
a branch includes many trial and error commits, the maintainers may squash
some commits before merging.

Pull requests
-------------
Always submit PRs via the indralab repository.
Give your PR a concise and clear title that describes without excessive detail
what the PR does. You should give more details in the description, pointing
out the important changes made and any additional remarks that are relevant.
If the PR fixes any issues, you can add "Fixes #xxx" to the text of the PR,
which, when merged, will also automatically close the issue.
The branch itself should have a short but recognizable name related to the
feature it adds or fixes rather than a generic name (e.g. patch, fix).

Commit messages
---------------
The commit message should typically consist of a single line describing what
the commit does. A good commit is one for which a clear and concise commit
message is necessary and sufficient - irrespective of how much code the commit
changes. A good set of guidelines can be found
[here](https://chris.beams.io/posts/git-commit/).

Code style
----------
Please follow [PEP8 guidelines](https://www.python.org/dev/peps/pep-0008/)
when implementing new code. If modifying existing
code, we ask that you do not mix extensive stylistic changes with meaningful
code changes. Any stylistic changes should be submitted in a separate PR.

The most important stylistic requirements are:
- use 4 spaces for indentation instead of tabs
- wrap lines to max 80 characters
- name variables and functions all lowercase with underscore as a separator
(e.g. `some_variable`)
- name classes with starting letters capitalized and no separator
(e.g. `SomeClass`)

In addition, functions or classes that are not meant to be part of the API
of the given module (for instance helper functions that a user wouldn't
directly call) should be prefixed with an underscore. These then won't show
up and clutter the auto-generated API documentation.

Python version compatibility
----------------------------
PyBioPAX is compatible with Python versions 3.6 and later. Therefore, features
only available in Python versions later than 3.6 should not be used.

Documentation
-------------
All API functions (i.e. functions that a user can call) and classes need to be
documented via docstrings. Please follow the
[NumPy documentation style](https://numpydoc.readthedocs.io/en/latest/format.html)
when adding docstrings to functions and classes.

The docstring
- is surrounded by triple double-quotes,
- starts with a one line short summary on the same line as the starting quotes,
- after the short summary and an empty line, can contain an arbitrary length
extended summary,
- lists all arguments, their types and descriptions in a Parameters block
- lists all returned values, their types and descriptions in a Returns block

To verify that the documentation build is working, go into the `doc` folder
and run `make html`. Warnings and errors indicate any issues during the build
process. The resulting HTML documentation can be opened with a browser from
`doc/_build/html/index.html` and inspected to verify that it looks as
intended.

Testing
-------
PyBioPAX is tested using the `pytests` package and we recommend running tests
locally with `tox` as `tox -e py`.

All new functionalities added should also be tested. Similarly, fixed bugs
should have regression tests added. Normally, any test file with `test` in its
name and any functions/classes that have `test/Test` in their names in these
files will be automatically discovered and tested. Tests should generally be
included in `pybiopax/tests`, and new tests should be placed in the appropriate
existing file, if possible. Otherwise, add a new file using the
`test_<module_name>.py` naming convention. Where possible, tests should be
short and focused. If the newly added test requires special dependencies or
other preliminary setup, the `.github/workflows/tests.yml` configuration for
GitHub Actions needs to be updated to make the test work. PRs will not be
merged unless tests are passing.

Logging
-------
Instead of using `print` for printing information to stdout, use the `logging`
module to first create an approproately named logger,
as for instance `logger = logging.getLogger(__name__)` and then use the
appropriate level of logging (typically debug, info, warning or error) with
the logger object to print messages. The configuration of the logging format
is uniform across PyBioPAX without further configuration needed for each
individual logger instance. In addition, by using `__name__` to instantiate
the logger, the hierarchy of logger objects across the module is maintained
making it easier to control logging at various levels.

New dependencies
----------------
When adding new features, using built-in Python libraries or
packages that are already standard dependencies of PyBioPAX is strongly
preferred.  In case a new dependency needs to be used, that dependency needs to
be
- added to the install list or one of the extras list in setup.py
- added to the installation instructions in the documentation if any special
  instructions are needed for setup
- added to doc/requirements.txt as an installed dependency
- added to .github/workflows/tests.yml unless installed via setup.py
This folder contains source material for a Journal of Open Source Software
submission:
- paper.md: the main manuscript
- paper.bib: BibTex citations for the manuscript
- build.sh: a script to compile the manuscript into a pdf

Submission via: https://joss.theoj.org/papers/new
---
title: 'PyBioPAX: biological pathway exchange in Python'

tags:
- Python
- systems biology
- biological pathways
- networks

authors:
- name: Benjamin M. Gyori
  orcid: 0000-0001-9439-5346
  affiliation: 1
- name: Charles Tapley Hoyt
  orcid: 0000-0003-4423-4370
  affiliation: 1

affiliations:
- name: Laboratory of Systems Pharmacology, Harvard Medical School
  index: 1
  ror: 03vek6s52

date: 7 December 2021
bibliography: paper.bib
repository: indralab/pybiopax
---

# Statement of need

Understanding the complex molecular processes governing how cells respond to
external stimuli crucially relies on prior knowledge about signaling,
regulatory, and metabolic pathways. Standardized representations are
necessary to exchange such pathway knowledge and allow interoperability between
tools. BioPAX [@demir2010biopax] is a widely used pathway exchange format that
is formally defined in the [BioPAX Language Specification](http://www.biopax.org/release/biopax-level3-documentation.pdf).
BioPAX is serialized into the Web Ontology Language (OWL)
format, typically as RDF/XML. Software support for parsing, serializing, and
finding patterns in BioPAX models is implemented in the Paxtools Java package
[@demir2013using]. However, interacting with Paxtools is difficult from Python,
and requires running a Java Virtual Machine via cross-language frameworks such
as pyjnius [@pyjnius]. Therefore, there is a need for native Python software support
for BioPAX to facilitate integration with widely used systems biology tools
(e.g., PySB [@lopez2013pysb], Tellurium [@medley2018tellurium], PyBEL
[@hoyt2017pybel]), and pathway analysis workflows more generally.

## State of the field

Support for the BioPAX language is implemented in the Paxtools Java package
[@demir2013using] and a wrapper extension around it called PaxtoolsR enabling
its usage from an R environment [@luna2016paxtoolsr]. A graphical tool for the
visualization of BioPAX called ChiBE [@babur2010chibe] is also available as a
Java package. There also exist dedicated analysis packages for pathway
enrichment such as the BioPAX-Parser Java package [@agapito2020biopax] and
several tools solving the conversion of BioPAX representations into modeling
formalisms such as the BioASF Java package [@haydarlou2016bioasf]. Overall,
however, there is no Python library implementing the BioPAX object model, and
enabling the manipulation and analysis of BioPAX models.

# Summary

We present PyBioPAX, a Python software package to process and manipulate BioPAX
models. PyBioPAX implements the BioPAX Level 3 object model as a set of Python
classes, and implements a BioPAX OWL processor to deserialize BioPAX content
from OWL files or strings into these objects. Once a BioPAX model and all its
linked elements are deserialized into Python objects, they can be traversed and
modified in memory. PyBioPAX supports serialization of BioPAX models into
OWL/XML files compatible with other tools in the BioPAX ecosystem.

PyBioPAX implements the BioPAX OWL semantics where object attributes can be
subtyped (e.g., "display name" is a subtype of "name") using Python property
attributes and getter/setter functions. It also supports exposing
"inverse links" between objects; for example, a BioPAX Xref object, which
represents a cross-reference, exposes a list of `xref_of` links back to the
objects of which it is a cross-reference. Again, the coherence of these links at
the level of a BioPAX model is guaranteed through the use of Python property
attributes. The inverse links contribute to the efficient traversal of BioPAX
models by allowing to link from e.g., one participant of a reaction to the
reaction itself and its other participants. To facilitate model traversal,
PyBioPAX provides a module to iterate over linked objects that satisfy a path
constraint string specification from a given starting object.

PyBioPAX also provides a client to the Pathway Commons web
service [@rodchenkov2020pathway] that makes three different graph query types
available: paths-from-to, paths-between, and neighborhood to extract subsets of
knowledge aggregated from structured sources in Pathway Commons
(e.g., Reactome [@jassal2020reactome]) as BioPAX models. PyBioPAX further provides
web service clients for processing BioPAX content from other pathway databases
including NetPath [@kandasamy2010netpath], and multiple members of the
BioCyc database collection [@karp2019biocyc].

# Case studies

In the following case studies, we demonstrate the role of PyBioPAX in
qualitative and quantitative analyses driven by BioPAX models.

## Traversing Pathway Commons

We demonstrate using PyBioPAX to process the Pathway Commons
version 12 (PC12) "detailed" model [BioPAX OWL file](https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.Detailed.BIOPAX.owl.gz),
to traverse it, and then to extract several biologically motivated motifs
corresponding to the following questions:

1. Which controllers of the catalyses of biochemical reactions require a
   co-factor?
2. Which controllers of the catalyses of biochemical reactions are in a
   phosphorylated state?
3. Which biochemical reactions constitute a simple phosphorylation event?
4. Which complexes contain a protein bound to one or more small molecules?
5. What are all the features (e.g., post-translational modifications, fragments)
   of a given protein?

Our implementations of these queries in the
corresponding [Jupyter notebook](https://nbviewer.org/github/indralab/pybiopax/blob/0.1.0/notebooks/Explore%20Pathway%20Commons.ipynb)
identified nearly 4M objects in PC12, 83 controllers that need co-factors, 1,283
controllers that are in a phosphorylated state, 15,332 simple phosphorylation
reactions, 13,338 proteins bound to a single small molecule, and 184 proteins
bound to two more small molecules.

Additionally, PyBioPAX enabled us to write queries to find superlative
entities.  For instance, we found that the protein with the most modifications
was NOTCH1, with 38 modifications. We further found that the RNA transcript of
KTN1 had the most interactions (947), and AR had the most interactions of any
protein (106).

## Gene set enrichment on Reactome pathways

Expert-curated pathways have been used as a means of dimensionality reduction
and interpretation of transcriptomics data. However, most prior methods are limited
to using pre-defined pathway lists (e.g., [@emon2020] only includes KEGG
pathways). Here, we demonstrate using PyBioPAX to implement a similar workflow
that is generally applicable to any pathway definition originating from BioPAX
content, represented as PyBioPAX models.

First, we obtained all human pathways as PyBioPAX models through PyBioPAX's API
for the Reactome web service. We then traversed each model to identify
physical entities representing proteins, aggregate their cross-references,
and ultimately construct a list of HGNC gene identifiers for each pathway.
Second, we collected curated transcriptomics experiments from the
CREEDS database [@wang2016] that list the
differentially expressed (DE) genes resulting from select drug perturbations,
gene knockouts, gene overexpressions, and diseases.

Finally, we used Fisher's exact test in an all-by-all
comparison of the lists of DE genes for each perturbation experiment against
the lists of genes whose proteins are present in each pathway. From this
matrix we identified anti-correlations between drug perturbation experiments and
gene perturbation experiments via the Pearson correlation coefficient. For
example, this highlighted a strong relationship between estradiol and GPER1,
suggesting GPER1 activation as a mechanism of action for estradiol.

The corresponding Jupyter notebook can be
found [here](https://nbviewer.org/github/indralab/pybiopax/blob/0.1.0/notebooks/Pathway%20Anticorrelation%20Analysis.ipynb).

# Availability and usage

PyBioPAX is available as a package on [PyPI](https://pypi.org/project/pybiopax)
with the source code available
at [https://github.com/indralab/pybiopax](https://github.com/indralab/pybiopax)
and documentation available
at [https://pybiopax.readthedocs.io/](https://pybiopax.readthedocs.io). The
repository also contains an interactive Jupyter notebook tutorial and notebooks
for the two case studies described above.

In addition to our case studies, PyBioPAX has been integrated into INDRA
[@indra] and serves as the primary entry point for processing BioPAX content
into INDRA Statements through the traversal of a BioPAX model. It has also been
used in [@weber-etal-2021-extend] to process BioPAX content from Reactome into a
node-edge graph used to train a machine-learning model used to improve natural
language processing.

# Acknowledgements

Development of this software was supported by the Defense Advanced Research
Projects Agency under award W911NF-15-1-0544 and the National Cancer Institute
under award U54-CA225088.

# References
.. pybiopax documentation master file, created by
   sphinx-quickstart on Thu May 28 14:40:01 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyBioPAX documentation
======================

.. toctree::
   :maxdepth: 3

   modules/api
   modules/biopax
   modules/pc_client
   modules/paths
   modules/xml_util


* :ref:`genindex`
* :ref:`modindex`
Path finding
============

.. automodule:: pybiopax.paths
    :members:
    :show-inheritance:
BioPAX object model
-------------------

.. automodule:: pybiopax.biopax
    :members:
    :show-inheritance:

.. automodule:: pybiopax.biopax.model
    :members:
    :show-inheritance:

.. automodule:: pybiopax.biopax.base
    :members:
    :show-inheritance:

.. automodule:: pybiopax.biopax.interaction
    :members:
    :show-inheritance:

.. automodule:: pybiopax.biopax.physical_entity
    :members:
    :show-inheritance:

.. automodule:: pybiopax.biopax.util
    :members:
    :show-inheritance:
Pathway Commons client
======================

.. automodule:: pybiopax.pc_client
    :members:
    :show-inheritance:
PyBioPAX API
============

.. automodule:: pybiopax.api
    :members:
    :show-inheritance:
XML/OWL processing utilities
============================

.. automodule:: pybiopax.xml_util
    :members:
    :show-inheritance:
