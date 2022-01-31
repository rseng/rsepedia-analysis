* Textnets version:
* Python version:
* Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

There are many ways you can contribute.

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/jboynyc/textnets/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Textnets could always use more documentation, whether as part of the
official Textnets docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/jboynyc/textnets/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up **textnets** for local development.

1. Fork the ``textnets`` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/textnets.git

3. Install your local copy into a virtual environment. This is how you set up
   your fork for local development::

    $ cd textnets/
    $ poetry install
    $ poetry install -E doc

   If you use `nix <https://nixos.org/nix>`__, you can also invoke
   ``nix-shell`` in the repository to quickly create a development environment
   with the included ``shell.nix`` file.

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, format and lint your changes and run the
   unit tests::

    $ make lint
    $ make test

6. Commit your changes::

    $ git add .
    $ git commit -m "Your detailed description of your changes."

7. Push you changes and submit a pull request through GitHub::

    $ git push origin name-of-your-bugfix-or-feature

   Alternately, if you'd rather avoid using GitHub, email a patch to the
   maintainer. See https://git-send-email.io/ for instructions.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add
   illustrative code examples to the tutorial.
3. The pull request should work for Python 3.7, 3.8, and 3.9. Check
   https://github.com/jboynyc/textnets/actions/workflows/ci.yml to make sure
   that the tests pass.

Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ poetry version patch # possible: major / minor / patch
$ git commit -a -m "Bump to version $(poetry version -s)"
$ git tag -a v$(poetry version -s)
$ make push

Tagged releases are automatically published to PyPI.
=======
Credits
=======

Author & Maintainer
-------------------

`John D. Boy <https://www.jboy.space>`__, Leiden University

Contributors
------------

- Jeremy Foote (`jdfoote <https://github.com/jdfoote>`__) contributed to the documentation.
- Your name here?
=======
History
=======

0.7.0 (2021-11-12)
------------------
* Adds abilitiy to save and load an instance of `Corpus`, `Textnet` and
  `params` to and from file using `Corpus.save`, `load_corpus`, `Textnet.save`,
  `load_textnet`, `params.save` and `params.load`. The same file can be used
  for all three kinds of objects, so all relevant data for a project can be
  saved in one file.
* Some further optimization of backbone extraction.
* Adds bipartite centrality measures (HITS, CoHITS and BiRank) and a bipartite
  clustering coefficient.
* Improved testing and type hints.
* Expanded documentation with advanced topics, including the new save/load
  feature and interacting with other libraries for network analysis and machine
  learning. Docs now use the PyData theme.
* Improvements to visualization. When plotting, nodes and edges can now be
  scaled by any attribute.
* Breaking change: Term weighing now happens in the ``corpus`` submodule, so
  the ``sublinear`` argument has to be passed to the methods for term
  extraction (``tokenized``, ``noun_phrases`` and ``ngrams``). This change will
  make it easier to add additional term extraction and weighing options.
* Adds ``tn.init_seed()`` utility to quickly initialize pseudorandom number
  generator.
* Adds Python 3.10 compatibility.
* Updates dependencies, including ``igraph`` with some relevant upstream
  changes contributed by yours truly, as well as spaCy.

0.6.0 (2021-10-14)
------------------
* Adds `params` as a container for global parameters. This makes it possible to
  fix the random seed and to change the resolution parameter for the community
  detection algorithm, among others. If the parameter ``autodownload`` is set
  to true, **textnets** will attempt to download all required spaCy language
  models automatically.
* Added HTML representation for the root module that displays versions of key
  dependencies.
* Added back string representations of `Corpus` and `TextnetBase`-derived
  classes.
* Adds a `Corpus.from_dict` method.
* `Corpus` now exposes the ``lang`` attribute, so the corpus language can be
  set after initialization of a class instance.
* The bipartite layout optionally used by `Textnet.plot` is now horizontal, so
  node types are arranged in columns rather than rows. That way node labels are
  less likely to overlap.
* Adds ``label_nodes`` argument to the `Textnet.plot` method to label both types
  of nodes. Defaults to ``False``.
* Adds ``node_opacity`` and ``edge_opacity`` arguments for `Textnet.plot`.
* Makes polygons marking clusters more visually appealing by adding opacity.
* Probably fixes `a bug <https://github.com/jboynyc/textnets/issues/30>`_ that
  would occasionally result in an exception being raised during plotting
  (``IndexError: color index too large``).
* When initializing an instance of the `Textnet` class, you can now optionally
  pass the argument ``connected=True``, in which case only the largest
  component of the underlying network is kept. When creating a one-mode
  projection using `Textnet.project`, a ``connected`` argument can also be
  passed.
* Adds `TextnetBase.save_graph` to save the underlying graph (for instance, for
  further processing in Gephi).
* Improved and extended documentation and docstrings.
* Update dependencies.

0.5.4 (2021-09-24)
------------------
* Fix the cross-platform build and deploy pipeline.
* Create binary packages for conda-forge.
* Otherwise, no substantive change from previous release.

0.5.3 (2021-09-24)
------------------
* Adds Catalan, Macedonian and Russian language models.
* Significantly speeds up backbone extraction by implementing the disparity
  filter integrand in Cython. (If the compiled extension cannot be loaded for
  some reason, it falls back on an interpreted function.)
* `PyPI <http://pypi.org/project/textnets>`_ *should* now receive binary wheels
  for Mac, Windows and Linux (via GitHub Actions) to ease installation on each
  platform.
* Improved type annotations.
* Update several dependencies.

0.5.2 (2021-08-24)
------------------
* Improve the handling of edge cases when initializing the `Corpus` and
  `Textnet` classes, such as empty data being provided.
* Added ability to run the tutorial in the documentation interactively using
  `thebe <https://thebelab.readthedocs.io/>`_.
* Update to spacy 3.1 and bump other dependencies.

0.5.1 (2021-07-06)
------------------
* Adds `Corpus.ngrams` method as alternative to `Corpus.noun_phrases`. This is
  useful when working in languages that do not have noun chunks, such as
  Chinese.
* Fixes a bug in `Corpus.from_files`.
* Introduces HTML representations of core classes for nicer integration in
  Jupyter notebooks.
* Updates several dependencies.

0.5.0 (2021-06-28)
------------------
* Migrate continuous integration testing from Travis to GitHub Actions.
* Continuous integration tests now run for MacOS and Windows too.
* Update to Spacy 3 and bump other dependency versions.
* Improvements to documentation.
* Handle dependencies and build project using Poetry (PEP 517 and 518).
* Remove deprecated command-line interface.

0.4.11 (2020-11-09)
-------------------
* Python 3.9 compatibility!
* Updated documentation with conda-forge installation option.
* Bump versions for numerous dependencies.

0.4.10 (2020-09-14)
-------------------
* Add ``cairocffi`` dependency and update installation docs.
* Bump ``leidenalg`` dependency to version 0.8.1.

0.4.9 (2020-07-15)
------------------
* Add ``color_clusters`` option to `Textnet` plotting methods. This colors
  nodes according to their partition using a bespoke color palette.

0.4.8 (2020-07-10)
------------------
* The `Corpus` class now handles missing data (#13).
* Support for more corpus languages. If no statistical language model is
  available, `Corpus` tries to use a basic ("blank") model.
* Improved documentation around dependencies and language support.
* Added tests.

0.4.7 (2020-07-01)
------------------
* No substantive change from previous release.

0.4.6 (2020-07-01)
------------------
* Bump spacy dependency to version 2.3 because it includes several new language
  models.

0.4.5 (2020-06-29)
------------------
* `Textnet.plot` and `ProjectedTextnet.plot` now accept arguments to selectively
  suppress node or edge labels. ``node_label_filter`` and ``edge_label_filter``
  take a function that is mapped to the iterator of nodes and edges. Only nodes
  or edges for which the function returns ``True`` are displayed in the plot.
* `Corpus` now has a useful string representation.
* Documentation updates, particularly to show the label filter functionality.

0.4.4 (2020-06-19)
------------------

* Methods to report centrality measures in `TextnetBase` now return
  `pandas.Series` objects. This has some nice benefits, like seeing node labels
  alongside centrality measures and being able to call ``.hist()`` on them to
  visualize the distribution.
* Scaling of nodes by centrality in plots should bring out differences more
  clearly now.
* Improved and expanded tutorial. Among other things, it now uses short codes
  to specify language models.

0.4.3 (2020-06-17)
------------------

* Python 3.7 compatibility is here.
* New ``circular_layout`` option for `Textnet.plot`. This is based on "`Tidier
  Drawings <https://www.reingold.co/graph-drawing.shtml>`_" and looks very nice
  for some bipartite graphs.
* String representation of `Textnet` instances now gives helpful information.
* Updated documentation to note changed Python version requirement.

0.4.2 (2020-06-16)
------------------

* `ProjectedTextnet.plot` now takes an argument, ``alpha``, that allows for
  pruning the graph in order to visualize its "backbone." This is useful when
  working with hairball graphs, which is common when creating textnets. Right
  now, it uses Serrano et al.'s disparity filter. That means that edges with an
  alpha value greater than the one specified are discarded, so lower values
  mean more extreme pruning.
* Language models can now be specified using a short ISO language code.
* Bipartite networks can now be plotted using a layered layout (by Kozo
  Sugiyama). Simply pass ``sugiyama_layout=True`` to `Textnet.plot`.
* Incremental improvements to documentation.

0.4.1 (2020-06-12)
------------------

* Documented `TextnetBase` methods to output lists of nodes ranked by various
  centrality measures: `top_betweenness` and several more.
* Added `top_cluster_nodes` to output list of top nodes per cluster found via
  community detection. This is useful when trying to interpret such clusters as
  themes/topics (in the projected word-to-word graph) or as groupings (in the
  document-to-document graph).
* Small additions to documentation.

0.4.0 (2020-06-11)
------------------

Lots of changes, some of them breaking, but overall just providing nicer
abstractions over the underlying pandas and igraph stuff.

* Introduced `TextnetBase` and `ProjectedTextnet` classes, and made `Textnet` a
  descendant of the former.
* Improved code modularity to make it easier to add features.
* `Corpus` is now based on a Series rather than a DataFrame.
* Added methods for creating an instance of `Corpus`: `from_df`, `from_csv`,
  `from_sql`.
* Expanded and improved documentation.
* Added bibliography to documentation using a Sphinx bibtex plugin.
* A first contributor!

0.3.6 (2020-06-03)
------------------

* Small change to *finally* get automatic deployments to PyPI to work.

0.3.5 (2020-06-03)
------------------

* Overall improvements to documentation.
* Added ``label_edges`` argument to `Textnet.plot`.

0.3.4 (2020-06-02)
------------------

* Integrated self-contained example that can be downloaded as Jupyter notebook
  into tutorial.
* Still trying to get automatic deployments to PyPI working.

0.3.3 (2020-06-02)
------------------

* More documentation.
* Attempt to get automatic deployments to PyPI working.

0.3.2 (2020-06-02)
------------------

* Set up continuous integration with Travis CI.
* Set up pyup.io dependency safety checks.
* Expanded documentation.
* A logo!

0.3.2 (2020-05-31)
------------------

* Further improvements to documentation.

0.3.1 (2020-05-31)
------------------

* Improvements to documentation.

0.3.0 (2020-05-31)
------------------

* First release on PyPI.
=====================================
Textnets: text analysis with networks
=====================================

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/jboynyc/textnets-binder/trunk?filepath=Tutorial.ipynb
   :alt: Launch on Binder

.. image:: https://github.com/jboynyc/textnets/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/jboynyc/textnets/actions/workflows/ci.yml
   :alt: CI status

.. image:: https://readthedocs.org/projects/textnets/badge/?version=stable
   :target: https://textnets.readthedocs.io/en/stable/?badge=stable
   :alt: Documentation Status

.. image:: https://anaconda.org/conda-forge/textnets/badges/installer/conda.svg
   :target: https://anaconda.org/conda-forge/textnets
   :alt: Install with conda

.. image:: https://joss.theoj.org/papers/10.21105/joss.02594/status.svg
   :target: https://doi.org/10.21105/joss.02594
   :alt: Published in Journal of Open Source Software

**textnets** represents collections of texts as networks of documents and
words. This provides novel possibilities for the visualization and analysis of
texts.

.. figure:: https://textnets.readthedocs.io/en/dev/_static/impeachment-statements.svg
   :alt: Bipartite network graph

   Network of U.S. Senators and words used in their official statements
   following the acquittal vote in the 2020 Senate impeachment trial (`source
   <https://www.jboy.space/blog/enemies-foreign-and-partisan.html>`_).

The ideas underlying **textnets** are presented in this paper:

  Christopher A. Bail, "`Combining natural language processing and network
  analysis to examine how advocacy organizations stimulate conversation on social
  media`__," *Proceedings of the National Academy of Sciences of the United States
  of America* 113, no. 42 (2016), 11823–11828, doi:10.1073/pnas.1607151113.

__ https://doi.org/10.1073/pnas.1607151113

Initially begun as a Python implementation of `Chris Bail's textnets package
for R`_, **textnets** now comprises unique features for term extraction and
weighing, visualization, and analysis.

.. _`Chris Bail's textnets package for R`: https://github.com/cbail/textnets/

**textnets** is free software under the terms of the GNU General Public License
v3.

Features
--------

**textnets** builds on `spaCy`_, a state-of-the-art library for
natural-language processing, and `igraph`_ for network analysis. It uses the
`Leiden algorithm`_ for community detection, which is able to perform community
detection on the bipartite (word–group) network.

.. _`igraph`: http://igraph.org/python/
.. _`Leiden algorithm`: https://doi.org/10.1038/s41598-019-41695-z
.. _`spaCy`: https://spacy.io/

**textnets** seamlessly integrates with Python's excellent `scientific stack`_.
That means that you can use **textnets** to analyze and visualize your data in
Jupyter notebooks!

.. _`scientific stack`: https://numfocus.org/

**textnets** is easily installable using the ``conda`` and ``pip`` package
managers. It requires Python 3.7 or higher.

Read `the documentation <https://textnets.readthedocs.io>`_ to learn more about
the package's features.

Citation
--------

Using **textnets** in a scholarly publication? Please cite this paper:

.. code-block:: bibtex

   @article{Boy2020,
     author   = {John D. Boy},
     title    = {textnets},
     subtitle = {A {P}ython Package for Text Analysis with Networks},
     journal  = {Journal of Open Source Software},
     volume   = {5},
     number   = {54},
     pages    = {2594},
     year     = {2020},
     doi      = {10.21105/joss.02594},
   }

Learn More
----------

==================  =============================================
**Documentation**   https://textnets.readthedocs.io/
**Repository**      https://github.com/jboynyc/textnets
**Issues & Ideas**  https://github.com/jboynyc/textnets/issues
**Conda-Forge**     https://anaconda.org/conda-forge/textnets
**PyPI**            https://pypi.org/project/textnets/
**DOI**             `10.21105/joss.02594 <https://doi.org/10.21105/joss.02594>`_
**Archive**         `10.5281/zenodo.3866676 <https://doi.org/10.5281/zenodo.3866676>`_
==================  =============================================

.. image:: https://textnets.readthedocs.io/en/dev/_static/textnets-logo.svg
   :alt: textnets logo
   :target: https://textnets.readthedocs.io
   :align: center
   :width: 140
.. include:: ../CONTRIBUTING.rst
.. highlight:: python

===============
Advanced Topics
===============

Saving and loading your project
-------------------------------

In this example, we define a ``project_file`` to store the configuration
parameters, corpus, and textnet. If the file exists, they are loaded from file;
else they are created and saved to file.

.. code:: python

   from pathlib import Path
   import textnets as tn


   working_dir = Path(".")
   project_file = working_dir / "my_project.db"

   if project_file.exists():
       tn.params.load(project_file)
       corpus = tn.load_corpus(project_file)
       net = tn.load_textnet(project_file)
   else:
       my_params = {"seed": 42, "autodownload": True}
       tn.params.update(my_params)
       corpus = tn.Corpus(tn.examples.digitalisierung, lang="de")
       net = tn.Textnet(corpus.noun_phrases(normalize=True))
       tn.params.save(project_file)
       corpus.save(project_file)
       net.save(project_file)

This code would only require the corpus and textnet to be created once.
Subsequent runs of the script could skip ahead to visualization or analysis.
This saves time, but also helps ensure the reproducibility of results.

Using alternate community detection algorithms
----------------------------------------------

By default, **textnets** will use the Leiden algorithm to find communities in
bipartite and projected networks. You can, however, also use other algorithms.

(These examples assume that you have already created a bipartite `Textnet`
called ``net``.)

Implemented in igraph
~~~~~~~~~~~~~~~~~~~~~

When plotting a textnet, you can supply the arguments ``show_clusters`` or
``color_clusters``. These accept a boolean value, but you can also pass a
`VertexClustering <igraph.VertexClustering>`, which is the data structure used
by ``igraph``.

If you want to use Blondel et al.'s multilevel algorithm to color the nodes of
a projected textnet, you can do so as follows:

.. code:: python

   terms = net.project(node_type="term")

   # initialize the random seed before running community detection
   tn.init_seed()
   part = terms.graph.community_multilevel(weights="weight")

   print("Modularity: ", terms.graph.modularity(part, weights="weight"))

   terms.plot(label_nodes=True, color_clusters=part)

Alternately, we can also overwrite the textnet's ``clusters`` property::

   terms.clusters = part

To return to the default (clusters detected by the Leiden algorithm), simply delete the clusters property::

   del terms.clusters

Implemented in leidenalg
~~~~~~~~~~~~~~~~~~~~~~~~

The ``leidenalg`` package is installed as a dependency of **textnets**. It can
produce a variety of different partition types, and in some cases, you may want
to use a different one that the default. In this example, ``leidenalg`` is
instructed to use the method of "asymptotic surprise" to determine the graph
partition.

.. code:: python

   import leidenalg as la

   terms.clusters = la.find_partition(terms.graph,
                                      partition_type=la.SurpriseVertexPartition, 
                                      weights="weight",
                                      n_iterations=-1,
                                      seed=tn.params["seed"])

After setting the clusters like this, you can plot the network as before. You
can also output a list of nodes per partition::

   terms.top_cluster_nodes()

Implemented in cdlib
~~~~~~~~~~~~~~~~~~~~

The Community Discovery Library (`cdlib <https://cdlib.readthedocs.io/>`__)
implements a wide range of algorithms for community detection that aren't
available in ``igraph``. Some of them are also able to perform community
detection on the bipartite network.

In order to run this example, you first have to install ``cdlib``.

.. code:: python

   from cdlib.algorithms import infomap_bipartite, paris

The first example applies the Infomap community detection algorithm to the
bipartite network::

   # initialize the random seed before running community detection
   tn.init_seed()
   bi_node_community_map = infomap_bipartite(net.graph.to_networkx()).to_node_community_map()

   # overwrite clusters detected by Leiden algorithm
   net.clusters = bi_node_community_map
   print("Modularity: ", net.modularity)

   net.plot(label_nodes=True, show_clusters=True)

This example applies the Paris hierarchical clustering algorithm to the
projected network::

   docs = net.project(node_type="doc")

   # initialize the random seed before running community detection
   tn.init_seed()
   docs_node_community_map = paris(docs.graph.to_networkx()).to_node_community_map()

   # overwrite clusters detected by Leiden algorithm
   docs.clusters = docs_node_community_map
   print("Modularity: ", docs.modularity)

   docs.plot(label_nodes=True, color_clusters=True)

Implemented in karateclub
~~~~~~~~~~~~~~~~~~~~~~~~~

`Karate Club <https://karateclub.readthedocs.io/>`__ is a library of
machine-learning methods to apply to networks. Among other things, it also
implements community detection algorithms. Here's an example for using
community detection from ``karateclub`` with **textnets**.

This example requires you to first have installed ``karateclub``.

.. code:: python

   from karateclub import SCD

   cd = SCD(seed=tn.params["seed"])
   cd.fit(net.graph.to_networkx())

   net.clusters = list(cd.get_memberships().values())
   print("Modularity: ", net.modularity)

   np.plot(color_clusters=True, label_nodes=True)

Additional measures for centrality analysis
-------------------------------------------

The `tutorial` provides examples of using BiRank, betweenness, closeness and
(weighted and unweighted) degree to analyze a textnet. The `NetworkX
<https://networkx.org>`__ library implements a large variety of other
centrality measures that may also prove helpful that aren't available in
``igraph``, the library that ``textnets`` builds on, including additional
centrality measures for bipartite networks.

This example requires ``networkx`` to be installed.

.. code:: python

   import networkx as nx

   bi_btwn = nx.algorithms.bipartite.betweenness_centrality(net.graph.to_networkx())
   net.nodes["btwn"] = list(bi_btwn.values())
   docs.plot(scale_nodes_by="btwn")

   katz_centrality = nx.katz_centrality(docs.graph.to_networkx(), weight="weight")
   docs.nodes["katz"] = list(katz_centrality.values())
   docs.plot(scale_nodes_by="katz")

Alternative methods of term extraction and weighing
---------------------------------------------------

By default, **textnets** leverages spaCy language models to break up your
corpus when you call `noun_phrases`, `ngrams` or `tokenized`, and it uses
*tf-idf* term weights. There are many alternative ways of extracting terms and
weighing them, and by defining a simple function, you can use them with
**textnets**.

This example uses `YAKE! <http://yake.inesctec.pt/>`__, the popular library for
keyword extraction, to extract keywords from a corpus and weighs them according
to their significance.

This example requires ``yake`` to be installed.

.. code:: python

   import textnets as tn
   from yake import KeywordExtractor

   def yake(
      corpus: tn.Corpus,
      lang: str="en",
      ngram_size: int=3,
      top: int=50,
      window: int=2
   ) -> tn.corpus.TidyText:
      """Use YAKE keyword extraction to break up corpus."""
      kw = KeywordExtractor(
               lan=lang,
               n=ngram_size,
               top=top,
               windowsSize=window
           )
      tt = []
      for label, doc in corpus.documents.iteritems():
          for term, sig in kw.extract_keywords(doc):
              tt.append({"label": label, "term": term, "term_weight": 1-sig, "n": 1})
      return tn.corpus.TidyText(tt).set_index("label")

The result of calling ``yake`` on an instance of ``Corpus`` can be passed to
``Textnet``.
API Reference
=============

.. autosummary::
   :toctree: _api
   :template: custom-module-template.rst
   :recursive:

   textnets
.. include:: ../HISTORY.rst
=============
Bibliography
=============

.. bibliography:: refs.bib
   :all:
.. include:: ../AUTHORS.rst
.. highlight:: shell

============
Installation
============

**textnets** is included in `conda-forge`_ and the `Python Package Index`_, so
it can either be installed using `conda`_ or `pip`_.

.. _conda-forge: https://anaconda.org/conda-forge/textnets/
.. _`Python Package Index`: https://pypi.org/project/textnets/
.. _conda: https://conda.io/
.. _pip: https://pip.pypa.io

.. note::

   Please note that **textnets** requires Python 3.7 or newer to run.

Using conda
-----------

This is the preferred method for most users. The `Anaconda Python
distribution`_ is an easy way to get up and running with Python, especially if
you are on a Mac or Windows system.

.. _Anaconda Python distribution: https://www.anaconda.com/products/individual

Once it is installed you can use its package manager ``conda`` to install
**textnets**::

   $ conda install -c conda-forge textnets

This tells conda to install **textnets** from the conda-forge channel.

If you don't know how to enter this command, you can use the Anaconda Navigator
instead. It provides a graphical interface that allows you to install new
packages.

.. admonition:: Installing **textnets** in Anaconda Navigator

   1. Go to the **Environments** tab.
   2. Click the **Channels** button.
   3. Click the **Add** button.
   4. Enter the channel URL https://conda.anaconda.org/conda-forge/
   5. Hit your keyboard's **Enter** key.
   6. Click the **Update channels** button.
   7. Now you can install **textnets** in a new environment. (Make sure the
      package filter on the **Environments** tab is set to "all.")

Using pip
---------

Alternately, if you already have Python installed, you can use its package
manger to install **textnets**. In a `virtual environment`_, run::

   $ python -m pip install textnets

.. _`virtual environment`: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

Plotting
--------

.. sidebar::

    In rare cases you may have to `install CFFI`_ separately for plotting to
    work.

.. _install CFFI: https://cffi.readthedocs.io/en/latest/installation.html

**textnets** installs the `Cairo`_ graphics library as a dependency for
plotting. If you are using a Mac without Anaconda, you will probably have to
install Cairo separately using the `Homebrew`_ package manager.

.. _Cairo: https://www.cairographics.org/
.. _Homebrew: https://formulae.brew.sh/formula/cairo

Language Support
----------------

**textnets** can try to download the `language models`_ you need "on the fly"
if you set the ``autodownload`` parameter to ``True``. (It is off by default
because language models are frequently many hundreds of megabytes in size and
probably shouldn't be downloaded on a metered connection.)

>>> import textnets as tn
>>> tn.params["autodownload"] = True

You can also install the models manually by issuing a command like::

   $ python -m spacy download en_core_web_sm

After updating **textnets** you may also need to update the language models.
Run the following command to check::

   $ python -m spacy validate

.. _`language models`: https://spacy.io/usage/models#download

If there are no language models available for your corpus language, there may
be some `basic support <https://spacy.io/usage/models#languages>`_. Even in
that case, some languages (including Korean, Vietnamese, Thai, Russian, and
Ukrainian) require additional installs for tokenization support. Consult the
spaCy documentation for details.
========
Tutorial
========

This tutorial walks you through all the steps required to analyze and visualize
your data using **textnets**. The tutorial first presents a self-contained
example before addressing miscellaneous other issues related to using
**textnets**.

Example
-------

.. tip::

   Download this example as a Jupyter notebook so you can follow along:
   :jupyter-download:notebook:`tutorial`.

   You can also make this tutorial "live" so you can adjust the example code
   and re-run it.

   .. thebe-button:: Do it live!

To use **textnets** in a project, you typically start with the following import:

.. jupyter-execute::

   import textnets as tn

You can set a fixed seed to ensure that results are reproducible across runs of
your script (see :cite:t:`Sandve2013`):

.. jupyter-execute::

   tn.params["seed"] = 42

Construct the corpus from the example data:

.. jupyter-execute::

   corpus = tn.Corpus(tn.examples.moon_landing)

What is this `moon_landing` example all about? (Hint: Click on the output below
to see what's in our corpus.)

.. jupyter-execute::

   corpus

.. note::

   Hat tip to Chris Bail for this example data!

Next, we create the textnet:

.. jupyter-execute::

   t = tn.Textnet(corpus.tokenized(), min_docs=1)

We're using `tokenized` with all defaults, so **textnets** is removing stop
words, applying stemming, and removing punctuation marks, numbers, URLs and the
like. However, we're overriding the default setting for ``min_docs``, opting to
keep even words that appear in only one document (that is, a single newspaper
headline).

Let's take a look:

.. jupyter-execute::

   t.plot(label_nodes=True,
          show_clusters=True)

The ``show_clusters`` options marks the partitions found by the Leiden
community detection algorithm (see :doc:`here <la:multiplex>`). It identifies
document--term groups that appear to form part of the same theme in the texts.

You may be wondering: why is the moon drifting off by itself in the network
plot? That's because the word moon appears exactly once in each document, so
its *tf-idf* value for each document is 0.

Let's visualize the same thing again, but this time scale the nodes according
to their BiRank (a centrality measure for bipartite networks) and the edges
according to their weights.

.. jupyter-execute::

   t.plot(label_nodes=True,
          show_clusters=True,
          scale_nodes_by="birank",
          scale_edges_by="weight")

We can also visualize the projected networks.

First, the network of newspapers:

.. jupyter-execute::

    papers = t.project(node_type="doc")
    papers.plot(label_nodes=True)

As before in the bipartite network, we can see the *Houston Chronicle*,
*Chicago Tribune* and *Los Angeles Times* cluster more closely together.

Next, the term network:

.. jupyter-execute::

   words = t.project(node_type="term")
   words.plot(label_nodes=True,
              show_clusters=True)

Aside from visualization, we can also analyze our corpus using network metrics.
For instance, documents with high betweenness centrality (or "cultural
betweenness"; :cite:t:`Bail2016`) might link together themes, thereby stimulating
exchange across symbolic divides.

.. jupyter-execute::

   papers.top_betweenness()

As we can see, the *Los Angeles Times* is a cultural bridge linking the
headline themes of the East Coast newspapers to the others.

.. jupyter-execute::

   words.top_betweenness()

It's because the *Times* uses the word "walk" in its headline, linking the "One
Small Step" cluster to the "Man on Moon" cluster.

We can produce the term network plot again, this time scaling nodes according
to their betweenness centrality, and pruning edges from the network using
"backbone extraction" :cite:p`Serrano2009`.

We can also use ``color_clusters`` (instead of ``show_clusters``) to color
nodes according to their partition.

And we can filter node labels, labeling only those nodes that have a
betweenness centrality score above the median. This is particularly useful in
high-order networks where labeling every single node would cause too much
visual clutter.

.. jupyter-execute::

   words.plot(label_nodes=True,
              scale_nodes_by="betweenness",
              color_clusters=True,
              alpha=0.5,
              edge_width=[10*w for w in words.edges["weight"]],
              edge_opacity=0.4,
              node_label_filter=lambda n: n.betweenness() > words.betweenness.median())

Wrangling Text & Mangling Data
------------------------------

How to go from this admittedly contrived example to working with your own data?
The following snippets are meant to help you get started. The first thing is to
get your data in the right shape.

A textnet is built from a collection—or *corpus*—of texts, so we use the
`Corpus` class to get our data ready. Each of the following snippets assumes
that you have imported `Corpus` and `Textnet` like in the preceding example.

From a Dictionary
~~~~~~~~~~~~~~~~~

You may already have your texts in a Python data structure, such as a
dictionary mapping document labels (keys) to documents (values). In that case,
you can use the `from_dict <Corpus.from_dict>` method to construct your
`Corpus`.

.. code:: python

   data = {f"Documento {label+1}": doc for label, doc in enumerate(docs)}
   corpus = tn.Corpus.from_dict(data, lang="it")

You can specify which `language model <https://spacy.io/models>`__ you would
like to use using the ``lang`` argument. The default is English, but you don’t
have to be monolingual to use **textnets**. (Languages in `LANGS` are fully
supported since we can use spacy's statistical language models. Other languages
are only partially supported, so `noun_phrases` will likely not function.)

From Pandas
~~~~~~~~~~~

`Corpus` can read documents directly from pandas' `Series <pd:pandas.Series>`
or `DataFrame <pd:pandas.DataFrame>`; mangling your data into the appropriate
format should only take :doc:`one or two easy steps
<pd:getting_started/intro_tutorials/10_text_data>`. The important thing is to
have the texts in one column, and the document labels as the index.

.. code:: python

   corpus = tn.Corpus(series, lang="nl")
   # or alternately:
   corpus = tn.Corpus.from_df(df, doc_col="tekst", lang="nl")

If you do not specify ``doc_col``, **textnets** assumes that the first column
containing strings is the one you meant.

From a database or CSV file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also use `Corpus` to load your documents from a database or
comma-separated value file using `from_sql` and `from_csv` respectively.

.. code:: python

   import sqlite3

   with sqlite3.connect("documents.db") as conn:
       articles = tn.Corpus.from_sql("SELECT title, text FROM articles", conn)

As before, you do can specify a ``doc_col`` to specify which column contains
your texts. You can also specify a ``label_col`` containing document labels. By
default, `from_sql` uses the first column as the ``label_col`` and the first
column after that containing strings as the ``doc_col``.

.. code:: python

   blog = tn.Corpus.from_csv("blog-posts.csv",
                             label_col="slug",
                             doc_col="summary"
                             sep=";")

Both `from_sql` and `from_csv` accept additional keyword arguments that are
passed to `pandas.read_sql` and `pandas.read_csv` respectively.

From Files
~~~~~~~~~~

Perhaps you have each document you want to include in your textnet stored on
disk in a separate text file. For such cases, `Corpus` comes with a utility,
`from_files`. You can simply pass a path to it using a `globbing
<https://en.wikipedia.org/wiki/Glob_(programming)>`__ pattern:

.. code:: python

   corpus = tn.Corpus.from_files("/path/to/texts/*.txt")

You can also pass it a list of paths:

.. code:: python

   corpus = tn.Corpus.from_files(["kohl.txt", "schroeder.txt", "merkel.txt"],
                                 doc_labels=["Kohl", "Schröder", "Merkel"],
                                 lang="de")

You can optionally pass explicit labels for your documents using the argument
``doc_labels``. Without this, labels are inferred from file names by stripping
off the file suffix.

Break It Up
~~~~~~~~~~~

The textnet is built from chunks of texts. `Corpus` offers three methods for
breaking your texts into chunks: `tokenized`, `ngrams`, and `noun_phrases`. The
first breaks your texts up into individual words, the second into n-grams of
desired size, while the third looks for `noun phrases
<https://en.wikipedia.org/wiki/Noun_phrase>`__ such as “my husband,” “our prime
minister,” or “the virus.”

.. code:: python

   trigrams = corpus.ngrams(3)

.. code:: python

   np = corpus.noun_phrases(remove=["Lilongwe", "Mzuzu", "Blantyre"])

.. warning::
   For large corpora, some of these operations can be computationally intense.
   Use your friendly neighborhood HPC cluster or be prepared for your laptop to
   get hot.

An optional boolean argument, ``sublinear``, can be passed to `tokenized`,
`ngrams`, and `noun_phrases`. It decides whether to use sublinear (logarithmic)
scaling when calculating *tf-idf* term weights. The default is ``True``,
because sublinear scaling is considered good practice in the information
retrieval literature :cite:p:`Manning2008`, but there may be good reason to
turn it off.

Calling these methods results in another data frame, which we can feed to
`Textnet` to make our textnet.

Make Connections
----------------

A textnet is a `bipartite network
<https://en.wikipedia.org/wiki/Bipartite_graph>`__  of *terms* (words or
phrases) and *documents* (which often represent the people or groups who
authored them). We create the textnet from the processed corpus using the
`Textnet` class.

.. code:: python

   t = tn.Textnet(np)

`Textnet` takes a few optional arguments. The most important one is
``min_docs``. It determines how many documents a term must appear in to be
included in the textnet. A term that appears only in a single document creates
no link, so the default value is 2. However, this can lead to a very noisy
graph, and usually only terms that appear in a significant proportion of
documents really indicate latent topics, so it is common to pass a higher
value.

``connected`` is a boolean argument that decides whether only the largest
connected component of the resulting network should be kept. It defaults to
``False``.

``doc_attrs`` allows setting additional attributes for documents that become
node attributes in the resulting network. For instance, if texts represent
views of members of different parties, we can set a party attribute.

.. code:: python

   t = tn.Textnet(corpus.tokenized(), doc_attr=df[["party"]].to_dict())

Seeing Results
--------------

You are now ready to see the first results. `Textnet` comes with a utility
method, `plot <Textnet.plot>`, which allows you to quickly visualize the
bipartite network.

For bipartite network, it can be helpful to use a layout option, such as
``bipartite_layout``, ``circular_layout``, or ``sugiyama_layout``, which help
to spatially separate the two node types.

You may want terms that are used in more documents to appear bigger in the
plot. In that case, use the ``scale_nodes_by`` argument with the value
``degree``. Other useful options include ``label_term_nodes``,
``label_doc_nodes``, and ``label_edges``. These are all boolean options, so
simply pass the value ``True`` to enable them.

Finally, enabling ``show_clusters`` will draw polygons around detected groups
of nodes with a community structure.

Projecting
----------

Depending on your research question, you may be interested either in how terms
or documents are connected. You can project the bipartite network into a
single-mode network of either kind.

.. code:: python

   groups = t.project(node_type="doc")
   print(groups.summary)

The resulting network only contains nodes of the chosen type (``doc`` or
``term``). Edge weights are calculated, and node attributes are maintained.

Like the bipartite network, the projected textnet also has a `plot
<ProjectedTextnet.plot>` method. This takes an optional argument, ``alpha``,
which can help "de-clutter" the resulting visualization by removing edges. The
value for this argument is a significance value, and only edges with a
significance value at or below the chosen value are kept. What remains in the
pruned network is called the "backbone" in the network science literature.
Commonly chosen values for ``alpha`` are in the range between 0.2 and 0.6 (with
lower values resulting in more aggressive pruning).

In visualizations of the projected network, you may want to scale nodes
according to centrality. Pass the argument ``scale_nodes_by`` with a value of
"betweenness," "closeness," "degree," "strength," or "eigenvector_centrality."

Label nodes using the boolean argument ``label_nodes``. As above,
``show_clusters`` will mark groups of nodes with a community structure.

Analysis
--------

The tutorial above gives some examples of using centrality measures to analyze
your corpus. Aside from `top_betweenness`, the package also provides the
methods `top_closeness`, `top_degree` (for unweighted degree), `top_strength`
(for weighted degree), and `top_ev` (for eigenvector centrality). In the
bipartite network, you can use `top_birank`, `top_hits` and `top_cohits` to see
nodes ranked by variations of a bipartite centrality measure :cite:p:`He2017`.
By default, they each output the ten top nodes for each centrality measure.

In addition, you can use `top_cluster_nodes <TextnetBase.top_cluster_nodes>` to
help interpret the community structure of your textnet. Clusters can either be
interpreted as latent themes (in the word network) or as groupings of documents
using similar words or phrases (in the document network).

Saving
------

You can save both the network that underlies a textnet as well as
visualizations. Assuming you want to save the projected term network, called
``words``, that we created above, you can do so as follows:

.. code:: python

   words.save_graph("term_network.gml")

This will create a file in the current directory in Graph Modeling Language
(GML) format. This can then be opened by Pajek, yEd, Gephi and other programs.
Consult the docs for ``save_graph`` for a list of supported formats.

If instead you want to save a plot of a network, the easiest thing is to pass
the ``target`` keyword to the `Textnet.plot` method.

.. code:: python

   words.plot(label_nodes=True, color_clusters=True, target="term_network.svg")

Supported file formats include PNG, EPS and SVG.
.. include:: ../README.rst

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   tutorial
   advanced

   reference

   bibliography

* :ref:`genindex`

.. toctree::
   :maxdepth: 2
   :caption: Project Information

   contributing
   history
   authors
   license
{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module attributes

   .. autosummary::
      :toctree:
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
      :nosignatures:
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :template: custom-class-template.rst
      :nosignatures:
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:
   :special-members: __call__, __add__, __mul__

   {% block methods %}
   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :nosignatures:
   {% for item in methods %}
      {%- if not item.startswith('_') %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
