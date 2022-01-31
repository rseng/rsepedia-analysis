
---
title: 'AuDoLab: Automatic document labelling and classification  for extremely unbalanced data'
tags:
  - Python
  - One-class SVM
  - Unsupervised Document Classification
  - One-class Document Classification
  - LDA Topic Modelling
  - Out-of-domain Training Data
authors:
  - name: Arne Tillmann
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Anton Thielmann
    affiliation: 1
  - name: Gillian Kant
    affiliation: 1
  - name: Christoph Weisser
    affiliation: "1, 2"
  - name: Benjamin Säfken
    affiliation: "1, 2"
  - name: Alexander Silbersdorff
    affiliation: "1, 2"
  - name: Thomas Kneib
    affiliation: "1, 2"
affiliations:
 - name: Georg-August-Universität Göttingen, Göttingen, Germany
   index: 1
 - name: Campus-Institut Data Science (CIDAS), Göttingen, Germany
   index: 2
date: 24 April 2021
bibliography: paper.bib
---

# Summary

AuDoLab provides a novel approach to one-class document classification for heavily imbalanced datasets, even if labelled training data is not available.
Our package enables the user to create specific out-of-domain training data to classify a heavily underrepresented target class
in a document dataset using a recently developed integration of Web Scraping, Latent Dirichlet Allocation Topic Modelling and One-class Support Vector Machines [@Thielmann]. AuDoLab can achieve high quality results even on highly specific classification problems without the need to invest in the time and cost intensive
labelling of training documents by humans. Hence, AuDoLab has a broad range of scientific research or business real world applications. In the following, a few potential use cases will be briefly discussed that should illustrate the broad range of applications in various domains. For example AuDoLab could be used to identify emails with very specific topics such as fraud or money laundering that might have an extremely low prevalence. Similarly, AuDoLab could be used in the medical field to classify medical documents that are concerned with very specific topics such as heart attacks or dental problems. Furthermore, AuDoLab may be used to identify legal documents with very specific topics such as machine learning. Note that, the only limiting factor to the broad range of use cases, is the availability of out-of-domain training data, that can be generated via Web Scraping from IEEEXplore [@IEEE], ArXiv or PubMed. Given that a broad range of training documents can be obtained from these websites AuDoLab has a correspondingly broad range of applications. The following section provides an overview of AuDoLab. AuDoLab can be installed conveniently via pip. A detailed description of the package and installation and can be found in the packages repository or on the documentation website.^[https://AuDoLab.readthedocs.io]

# Statement of need

Unsupervised document classification is mainly performed to gain insight into the underlying topics of large text corpora.
In this process, documents covering highly underrepresented topics have only a minor impact on the algorithm's topic definitions. As a result, underrepresented topics can sometimes be "overlooked" and documents are assigned topic prevalences that do not reflect the underlying content.
Thus, labeling underrepresented topics in large text corpora is often done manually and can therefore be very labour-intensive and time-consuming.
AuDoLab enables the user to tackle this problem and perform unsupervised one-class document classification for heavily underrepresented document classes.
This leverages the results of one-class document classification using One-class Support Vector Machines (SVM) [@Scholkopf; @Manevitz] and extends them to the use case of severely imbalanced datasets.
This adaptation and extension is achieved by implementing a multi-level classification rule as visualised in the graph below.

![Classification Procedure.\label{fig:test2}](figures/tree.PNG){ width=100% }

The first part of the package allows the user to scrape training documents (scientific papers), ideally covering the target topic in which the user is interested, from IEEEXplore [@IEEE], ArXiv or PubMed. The user can search for multiple search terms and specify an individual search query and, in the case of IEEEXplore, scrape additional information such as the author names or the number of citations. Thus, an individually labelled (e.g., via author-keywords) training data set is created. Through the integration of pre-labelled out-of-domain training data, the problem of the heavily underrepresented target class can be circumvented, as large enough training corpora can be automatically generated.
Subsequently, the text data is preprocessed for the classification part. The text preprocessing includes common Natural Language Processing (NLP) text preprocessing techniques such as stopword removal and lemmatization.  As  document  representations  the  term  frequency-inverse  document  frequency (tf-idf) representations are chosen. The tf-idf scores are computed on a joint corpus from the web-scraped out-of-domain training data and the target text data.

The second and main part of the classification rule lies in the training of the one-class SVM [@Scholkopf]. As a training corpus, only the out-of-domain training data is used.  By adjusting hyperparameters, the user can create a strict or relaxed classification rule, that reflects the users belief about the prevalence of the target class inside the target data set and the quality of the scraped out-of-domain training data. The last part of the classification rule enables the user to control the classifiers results with the help of Latent Dirichlet Allocation (LDA) topic models [@Blei] (and e.g., wordclouds). Additionally, the user can generate interactive plots depicting the identified topics during the LDA topic modelling [@ldavis].

The second step can be repeated, depending on the users perceived quality of the classification results.

## Comparison with existing tools

At the moment no Python Package with a comparable functionality of AuDoLab is available, since AuDoLab is based on a novel and recently published classification prodcedure [@Thielmann].
Thereby, AuDoLab uses and integrates in particular a combination of Web Scraping, Topic Modelling and One-class Classifcation for which various individual packages are available. Details on the statistical methodology can be found in [@Thielmann]. An application of the methodology on a data set of patent data can found in [@Thielmann2021]. For Topic Modelling available packages are the LDA algorithm as implemented in the package Gensim [@rehurek_lrec] or the package TTLocVis [@Kant2020] for short and sparse text. Visual representations of the topics can be implemented with LDAvis [@ldavis]. The One-class SVM classification package is availabe in Scikit-learn [@scikit-learn]. Alternative further research could explore Deep Learning Algorithms as well [@Saefken2020; @Saefken2021].

# Acknowledgements

We thank the Campus-Institut Data Science (CIDAS), Göttingen, Germany for funding this project.

# References
.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/ArneTillmann/AuDoLab/issues.

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

AuDoLab could always use more documentation, whether as part of the
official AuDoLab docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/ArneTillmann/AuDoLab/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `AuDoLab` for local development.

1. Fork the `AuDoLab` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/AuDoLab.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv AuDoLab
    $ cd AuDoLab/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 AuDoLab tests
    $ python setup.py test or pytest
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.5, 3.6, 3.7 and 3.8, and for PyPy. Check
   https://travis-ci.com/ArneTillmann/AuDoLab/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ pytest tests.test_AuDoLab


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ bump2version patch # possible: major / minor / patch
$ git push
$ git push --tags

Travis will then deploy to PyPI if tests pass.
=======
Credits
=======

Development Lead
----------------

* Arne Tillmann <arne.tillmann.vellmar@gmail.com>

Contributors
------------
* Anton Thielmann <anton.thielmann@stud.uni-goettingen.de>
* Christoph Weisser <c.weisser@stud.uni-goettingen.de>
* Benjamin Säfken <benjamin.saefken@uni-goettingen.de>
* Alexander Silbersdorff <asilbersdorff@uni-goettingen.de>
* Thomas Kneib <Thomas.Kneib@wiwi.uni-goettingen.de>

Special thanks
--------------



This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
=======
History
=======

0.0.1 (2021-04-07)
------------------

* First release on PyPI.
=======
AuDoLab
=======

.. image:: https://img.shields.io/pypi/v/AuDoLab.svg
        :target: https://pypi.python.org/pypi/AuDoLab

.. image:: https://api.travis-ci.com/ArneTillmann/AuDoLab.svg?branch=main&status=passed
        :target: https://travis-ci.com/ArneTillmann/AuDoLab

.. image:: https://readthedocs.org/projects/audolab/badge/?version=latest
        :target: https://audolab.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://joss.theoj.org/papers/10.21105/joss.03719/status.svg
        :target: https://doi.org/10.21105/joss.03719

With AuDoLab you can perform Latend Direchlet Allocation on highly imbalanced datasets.

=======
Summary
=======

AuDoLab provides a novel approach to one-class document classification for heavily imbalanced datasets, even if labelled training data is not available. Our package enables the user to create specific out-of-domain training data to classify a heavily underrepresented target class in a document dataset using a recently developed integration of Web Scraping, Latent Dirichlet Allocation Topic Modelling and One-class Support Vector Machines. AuDoLab can achieve high quality results even on higly specific classification problems without the need to invest in the time and cost intensive labelling of training documents by humans. Hence, AuDoLab has a broad range of scientific research or business applications.

Unsupervised document classification is mainly performed to gain insight into the underlying topics of large text corpora. In this process, documents covering highly underrepresented topics have only a minor impact on the algorithm's topic definitions. As a result, underrepresented topics can sometimes be "overlooked" and documents are assigned topic prevalences that do not reflect the underlying content. Thus, labeling underrepresented topics in large text corpora is often done manually and can therefore be very labour-intensive and time-consuming. AuDoLab enables the user to tackle this problem and perform unsupervised one-class document classification for heavily underrepresented document classes.

============
Installation
============


Stable release
--------------

To install AuDoLab, run this command in your terminal (bash, PowerShell, etc.), given that you have python 3 and pip installed :

.. code-block:: console

    $ pip install AuDoLab

This is the preferred method to install AuDoLab, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for AuDoLab can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/ArneTillmann/AuDoLab

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/ArneTillmann/AuDoLab/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/ArneTillmann/AuDoLab
.. _tarball: https://github.com/ArneTillmann/AuDoLab/tarball/master

=====
Usage
=====
Before the actuall usage you want to download the stopwords for nltk by running::

    import nltk
    nltk.download('stopwords')

inside a python console.
To use AuDoLab in a project::

    from AuDoLab import AuDoLab

Then you want to create an instance of the AuDoLab class

    audo = AuDoLab.AuDoLab()

In this example we used publicly available data from the nltk package::

    from nltk.corpus import reuters
    import numpy as np
    import pandas as pd

    data = []

    for fileid in reuters.fileids():
        tag, filename = fileid.split("/")
        data.append(
            (filename,
             ", ".join(
                 reuters.categories(fileid)),
                reuters.raw(fileid)))

    data = pd.DataFrame(data, columns=["filename", "categories", "text"])

Then you want to scrape abstracts, e.g. from IEEE with the abstract scraper::

    scraped_documents = audo.get_ieee("https://ieeexplore.ieee.org/search
                                       /searchresult.jsp?newsearch=true&
                                       queryText=cotton&highlight=true&
                                       returnFacets=ALL&returnType=SEARCH&
                                       matchPubs=true&rowsPerPage=100&
                                       pageNumber=1\",
                                       pages=1)

The data as well as the scraped papers need to be preprocessed before use in the
classifier::

    preprocessed_target = audo.text_cleaning(data=data, column="text")

    preprocessed_paper = audo.text_cleaning(
        data=scraped_documents, column="abstract")

    target_tfidf, training_tfidf = audo.tf_idf(
        data=preprocessed_target,
        papers=preprocessed_paper,
        data_column="lemma",
        papers_column="lemma",
        features=100000,
    )

Afterwards we can train and use the classifiers and choose the desired
one::

    o_svm_result = audo.one_class_svm(
        training=training_tfidf,
        predicting=target_tfidf,
        nus=np.round(np.arange(0.001, 0.5, 0.01), 7),
        quality_train=0.9,
        min_pred=0.001,
        max_pred=0.05,
    )

    result = audo.choose_classifier(preprocessed_target, o_svm_result, 0)

And finally you can estimate the topics of the data::

    lda_target = audo.lda_modeling(data=result, num_topics=5)

    audo.lda_visualize_topics(type="pyldavis")

* Free software: GNU General Public License v3
* Documentation: https://AuDoLab.readthedocs.io.
AuDoLab.subclasses.one\_class\_svm
==================================

.. automodule:: AuDoLab.subclasses.one_class_svm
   :members:
   :undoc-members:
   :show-inheritance:
AuDoLab.subclasses.abstractscraper_pubmed
=========================================

.. automodule:: AuDoLab.subclasses.abstractscraper_pubmed
   :members:
   :undoc-members:
   :show-inheritance:
.. include:: ../CONTRIBUTING.rst
.. include:: ../HISTORY.rst
.. include:: ../AUTHORS.rst
AuDoLab.subclasses.abstractscraper
==================================

.. automodule:: AuDoLab.subclasses.abstractscraper
   :members:
   :undoc-members:
   :show-inheritance:
AuDoLab.subclasses.abstractscraper_arxiv
========================================

.. automodule:: AuDoLab.subclasses.abstractscraper_arxiv
   :members:
   :undoc-members:
   :show-inheritance:
AuDoLab.subclasses.preprocessing
================================

.. automodule:: AuDoLab.subclasses.preprocessing
   :members:
   :undoc-members:
   :show-inheritance:
.. include:: ../README.rst
AuDoLab.subclasses.tf\_idf
==========================

.. automodule:: AuDoLab.subclasses.tf_idf
   :members:
   :undoc-members:
   :show-inheritance:
Classes and Methods
===================

.. toctree::
   :maxdepth: 4


AuDoLab.AuDoLab main module
---------------------------

.. automodule:: AuDoLab.AuDoLab
 :members:
 :undoc-members:
 :show-inheritance:

Subpackages
-----------

.. toctree::
  :maxdepth: 4


  AuDoLab.abstractscraper
  AuDoLab.abstractscraper_arxiv
  AuDoLab.abstractscraper_pubmed
  AuDoLab.preprocessing
  AuDoLab.tf_idf
  AuDoLab.one_class_svm
  AuDoLab.lda
Welcome to AuDoLab's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   readme
   modules
   contributing
   authors
   history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
AuDoLab.subclasses.lda
======================

.. automodule:: AuDoLab.subclasses.lda
   :members:
   :undoc-members:
   :show-inheritance:
