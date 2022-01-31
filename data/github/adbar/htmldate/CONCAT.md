## Changelog


### 1.0.0
- faster and more accurate encoding detection
- simplified code base
- include support for Python 3.10 and dropped support for Python 3.5


### 0.9.1
- improved generic date parsing (thanks @RadhiFadlillah)
- specific support for French and Indonesian (thanks @RadhiFadlillah)
- additional evaluation for English news sites (kudos to @coreydockser & @rahulbot)
- bugs fixed


### 0.9.0
- improved exhaustive search
- simplified code
- bug fixes
- removed support for Python 3.4

### 0.8.1
- bugfixes

### 0.8.0
- `dateparser` and `regex` modules fully integrated
- patterns added for coverage
- smarter HTML doc loading

### 0.7.3
- dependencies updated and reduced: switch from `requests` to bare `urllib3`, make `chardet` standard and `cchardet` optional
- fixes: downloads, `OverflowError` in extraction

### 0.7.2
- compatibility with Python 3.9
- better speed and accuracy

### 0.7.1
- technical release: package requirements and docs wording

### 0.7.0
- code base and performance improved
- minimum date available as option
- support for Turkish patterns and CMS idiosyncrasies (thanks @evolutionoftheuniverse)

### 0.6.3
- more efficient code
- additional evaluation data

### 0.6.2
- performance and documentation improved

### 0.6.1
- code base restructured
- bugs fixed and further tests
- restored retro-compatibility with Python 3.4

### 0.6.0
- reduced number of packages dependencies
- introduced and tested optional dependencies
- more detailed documentation on readthedocs

### 0.5.6
- tests on Windows
- compataibility and code linting

### 0.5.5
- tests on Linux & MacOS
- bugs removed

### 0.5.4
- manually set maximum date
- better precision
- temporarily dropped support for Python 3.4

### 0.5.3
- coverage extension

### 0.5.2
- small bugs and coverage issues removed
- streamlined utils
- documentation added

### 0.5.1
- bugs corrected and cleaner code
- more errors caught and better test coverage

### 0.5.0
- significant speed-up after code profiling
- better support of free text detection (DE/EN)

### 0.4.1
- fixed lxml dependency
- reordered XPath-expressions

### 0.4.0
- refined and combined XPath-expressions
- better extraction of dates in free text
- better coverage and consistency issues solved

### 0.3.1
- improved consistency and further tests

### 0.3.0
- improvements in markup analysis along with more tests
- higher resolution for free text detection (e.g. DD/MM/YY)
- download mode (serial on command-line)

### 0.2.2
- better code consistency
- tested for Python2 and 3 with tox and coverage stats

### 0.2.1
- refined date comparisons
- debug and logging options
- more tests and test files
- extensive search can be disabled

### 0.2.0
- refined targeting of HTML structure
- better extraction logic for plain text cases
- further tests

### 0.1.2
- better extraction
- logging
- further tests
- settings

### 0.1.1
- tests functions (tox and pytest)
- retro-compatibility (python2)
- minor improvements

### 0.1.0
- minimum viable package
## How to contribute

Thank you for considering contributing to htmldate!

Here are some important resources:

  * [List of currently open issues](https://github.com/adbar/htmldate/issues) (no pretention to exhaustivity!)
  * [How to Contribute to Open Source](https://opensource.guide/how-to-contribute/)

There are many ways to contribute, you could:

  * Improve the documentation
  * Find bugs and submit bug reports
  * Submit feature requests
  * Write tutorials or blog posts
  * Write code

## Submitting changes

Please send a [GitHub Pull Request to htmldate](https://github.com/adbar/htmldate/pull/new/master) with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)).

**Working on your first Pull Request?** You can learn how from this series: [How to Contribute to an Open Source Project on GitHub](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github)

For further questions you can contact me on [GitHub issues](), [Twitter](https://twitter.com/adbarbaresi) or by [E-Mail](http://adrien.barbaresi.eu/contact.html)


Thanks,

Adrien---
title: 'htmldate: A Python package to extract publication dates from web pages'
tags:
  - Python
  - metadata extraction
  - date parsing
  - web scraping
  - natural language processing
authors:
  - name: Adrien Barbaresi
    orcid: 0000-0002-8079-8694
    affiliation: 1
affiliations:
 - name: Berlin-Brandenburg Academy of Sciences
   index: 1
date: 29 July 2020
bibliography: paper.bib
---



# Introduction


### Rationale


Metadata extraction is part of data mining and knowledge extraction. Being able to better qualify content allows for insights based on descriptive or typological information (e.g., content type, authors, categories), better bandwidth control (e.g., by knowing when webpages have been updated), or optimization of indexing (e.g., caches, language-based heuristics). It is useful for applications including database management, business intelligence, or data visualization. This particular effort is part of a methodological approach to derive information from web documents in order to build text databases for research, chiefly linguistics and natural language processing. Dates are critical components since they are relevant both from a philological standpoint and in the context of information technology.

Although text is ubiquitous on the Web, extracting information from web pages can prove to be difficult. Web documents come in different shapes and sizes mostly because of the wide variety of genres, platforms, and content management systems, and not least because of greatly diverse publication goals. In most cases, immediately accessible data on retrieved webpages do not carry substantial or accurate information: neither the URL nor the server response provide a reliable way to date a web document, that is to find out when it has been published or possibly modified. In that case it is necessary to fully parse the document or apply robust scraping patterns on it. Improving extraction methods for web collections can hopefully allow for combining both the quantity resulting from broad web crawling and the quality obtained by accurately extracting text and metadata and by rejecting documents which do not match certain criteria.


### Research context


Fellow colleagues are working on a lexicographic information platform [@GeykenEtAl:2017] at the language center of the Berlin-Brandenburg Academy of Sciences ([dwds.de](https://www.dwds.de/)). The platform hosts and provides access to a series of metadata-enhanced web corpora [@Barbaresi:2016]. Information on publication and modification dates is crucial to be able to make sense of linguistic data, that is, in the case of lexicography to determine precisely when a given word was used for the first time and how its use evolves through time.

Large "offline" web text collections are now standard among the research community in linguistics and natural language processing. The construction of such text corpora notably involves "crawling, downloading, 'cleaning' and de-duplicating the data, then linguistically annotating it and loading it into a corpus query tool" [@Kilgarriff:2007]. Web crawling [@Olston:2010] involves a significant number of design decisions and turning points in data processing, without which data and applications turn into a "Wild West" [@JoGebru:2020]. Researchers face a lack of information regarding the content, whose adequacy, focus, and quality are the object of a post hoc evaluation [@Baroni:2009]. Comparably, web corpora (i.e., document collections) usually lack metadata gathered with or obtained from documents. Between opportunistic and restrained data collection [@Barbaresi:2015], a significant challenge lies in the ability to extract and pre-process web data to meet scientific expectations with respect to corpus quality.


# Functionality


``htmldate`` finds original and updated publication dates of web pages using heuristics on HTML code and linguistic patterns. It operates both within Python and from the command-line. URLs, HTML files, or HTML trees are given as input, and the library outputs a date string in the desired format or ``None`` as the output is thouroughly verified in terms of plausibility and adequateness.

The package features a combination of tree traversal and text-based extraction, and the following methods are used to date HTML documents:

1. Markup in header: common patterns are used to identify relevant elements (e.g., ``link`` and ``meta`` elements) including Open Graph protocol attributes and a large number of content management systems idiosyncrasies
1. HTML code: The whole document is then searched for structural markers: ``abbr`` and ``time`` elements as well as a series of attributes (e.g. ``postmetadata``)
1. Bare HTML content: A series of heuristics is run on text and markup:

    - in ``fast`` mode the HTML page is cleaned and precise patterns are targeted
    - in ``extensive`` mode all potential dates are collected and a disambiguation algorithm determines the best one

Finally, a date is returned if a valid cue could be found in the document, corresponding to either the last update or the original publishing statement (the default), which allows for switching between original and updated dates. The output string defaults to ISO 8601 YMD format.

``htmldate`` is compatible with all recent versions of Python (currently 3.4 to 3.9). It is designed to be computationally efficient and used in production on millions of documents. All the steps needed from web page download to HTML parsing, scraping, and text analysis are handled, including batch processing. It is distributed under the GNU General Public License v3.0. Markup-based extraction is multilingual by nature, and text-based refinements for better coverage currently support German, English and Turkish.


# State of the art

Diverse extraction and scraping techniques are routinely used on web document collections by companies and research institutions alike. Content extraction mostly draws on Document Object Model (DOM) examination, that is, on considering a given HTML document as a tree structure whose nodes represent parts of the document to be operated on. Less thorough and not necessarily faster alternatives use superficial search patterns such as regular expressions in order to capture desirable excerpts.


## Alternatives

There are comparable software solutions in Python. The following date extraction packages are open-source and work out-of-the-box:

- ``articleDateExtractor`` detects, extracts, and normalizes the publication date of an online article or blog post [@articleDateExtractor],
- ``date_guesser`` extracts publication dates from a web pages along with an accuracy measure which is not tested here [@dateguesser],
- ``goose3`` can extract information for embedded content [@goose3],
- ``htmldate`` is the software package described here; it is designed to extract original and updated publication dates of web pages [@Barbaresi:2019],
- ``newspaper`` is mostly geared towards newspaper texts [@newspaper],
- ``news-please`` is a news crawler that extracts structured information [@HamborgEtAl:2017],


Two alternative packages are not tested here but that also could be used:

- ``datefinder`` [@datefinder] features pattern-based date extraction for texts written in English,
- if dates are nowhere to be found, using ``CarbonDate`` [@carbondate] can be an option, however this is computationally expensive.


## Benchmark

#### Test set

The experiments below are run on a collection of documents that are either typical for Internet articles (news outlets, blogs, including smaller ones) or non-standard and thus harder to process. They were selected from large collections of web pages in German. For the sake of completeness, a few documents in other languages were added (English, European languages, Chinese, and Arabic).

#### Evaluation

The evaluation script is available in the project repository: ``tests/comparison.py``. The tests can be reproduced by cloning the repository, installing all necessary packages and running the evaluation script with the data provided in the ``tests`` directory.

Only documents with dates that are clearly able to be determined are considered for this benchmark. A given day is taken as unit of reference, meaning that results are converted to ``%Y-%m-%d`` format if necessary in order to make them comparable.

#### Time

The execution time (best of 3 tests) cannot be easily compared in all cases as some solutions perform a whole series of operations which are irrelevant to this task.

#### Errors

``goose3``'s output is not always meaningful and/or in a standardized format, so these cases were discarded. news-please seems to have trouble with some encodings (e.g., in Chinese), in which case it leads to an exception.


## Results

The results in Table 1 show that date extraction is not a completely solved task but one for which extractors have to resort to heuristics and guesses. The figures documenting recall and accuracy capture the real-world performance of the tools as the absence of a date output impacts the result.


| Python Package | Precision | Recall | Accuracy | F-Score | Time |
| --- | --- | --- | --- | --- | --- |
| newspaper 0.2.8 	| 0.888 	| 0.407 	| 0.387 	| 0.558 	| 81.6 |
| goose3 3.1.6 		| 0.887 	| 0.441 	| 0.418 	| 0.589 	| 15.5 |
| date_guesser 2.1.4 	| 0.809 	| 0.553 	| 0.489 	| 0.657 	| 40.0 |
| news-please 1.5.3  	| 0.823 	| 0.660 	| 0.578 	| 0.732 	| 69.6 |
| articleDateExtractor 0.20 | 0.817 	| 0.635 	| 0.556 	| 0.714		| 6.8 |
| htmldate 0.7.0 *(fast)* | **0.903** 	| 0.907 	| 0.827 	| 0.905 	| **2.4** |
| htmldate[all] 0.7.0 *(extensive)* | 0.889 	| **1.000** 	| **0.889** 	| **0.941** 	| 3.8 |

: 225 web pages containing identifiable dates (as of 2020-07-29)


Precision describes if the dates given as output are correct: ``newspaper`` and ``goose3`` fare well precision-wise but they fail to extract dates in a large majority of cases (poor recall). The difference in accuracy between ``date_guesser`` and ``newspaper`` is consistent with tests described on the website of the former.

It turns out that ``htmldate`` performs better than the other solutions overall. It is also noticeably faster than the strictly comparable packages (``articleDateExtractor`` and ``date_guesser``). Despite being measured on a sample, the higher accuracy and faster processing time are highly significant. Especially for smaller news outlets, websites, and blogs, as well as pages written in languages other than English (in this case mostly but not exclusively German), ``htmldate`` greatly extends date extraction coverage without sacrificing precision.



#### Note on the different versions:

- ``htmldate[all]`` means that additional components are added for performance and coverage. They can be installed with ``pip/pip3/pipenv htmldate[all]`` and result in differences with respect to accuracy (due to further linguistic analysis) and potentially speed (faster date parsing).
- The fast mode does not output as many dates (lower recall) but its guesses are more often correct (better precision).


# Acknowledgements

This work has been supported by the ZDL research project (*Zentrum für digitale Lexikographie der deutschen Sprache*, [zdl.org](https://www.zdl.org/)). Thanks to Yannick Kozmus (evaluation), user evolutionoftheuniverse (patterns for Turkish) and further [contributors](https://github.com/adbar/htmldate/graphs/contributors) for testing and working on the package. Thanks to Daniel S. Katz, Geoff Bacon and Maarten van Gompel for reviewing this JOSS submission.

The following Python modules have been of great help: ``lxml``, ``ciso8601``, and ``dateparser``. A few patterns are derived from ``python-goose``, ``metascraper``, ``newspaper`` and ``articleDateExtractor``; this package extends their coverage and robustness significantly.


# References


htmldate: find the publication date of web pages
================================================

.. image:: https://img.shields.io/pypi/v/htmldate.svg
    :target: https://pypi.python.org/pypi/htmldate
    :alt: Python package

.. image:: https://img.shields.io/pypi/pyversions/htmldate.svg
    :target: https://pypi.python.org/pypi/htmldate
    :alt: Python versions

.. image:: https://readthedocs.org/projects/htmldate/badge/?version=latest
    :target: https://htmldate.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/codecov/c/github/adbar/htmldate.svg
    :target: https://codecov.io/gh/adbar/htmldate
    :alt: Code Coverage

.. image:: https://static.pepy.tech/badge/htmldate/month
    :target: https://pepy.tech/project/htmldate
    :alt: Downloads

|

:Code:           https://github.com/adbar/htmldate
:Documentation:  https://htmldate.readthedocs.io
:Issue tracker:  https://github.com/adbar/htmldate/issues

|

Find original and updated publication dates of any web page. From the command-line or within Python, all the steps needed from web page download to HTML parsing, scraping, and text analysis are included.

In a nutshell
-------------

|

.. image:: docs/htmldate-demo.gif
    :alt: Demo as GIF image
    :align: center
    :width: 80%
    :target: https://htmldate.readthedocs.org/

|

With Python:

.. code-block:: python

    >>> from htmldate import find_date
    >>> find_date('http://blog.python.org/2016/12/python-360-is-now-available.html')
    '2016-12-23'
    >>> find_date('https://netzpolitik.org/2016/die-cider-connection-abmahnungen-gegen-nutzer-von-creative-commons-bildern/', original_date=True)
    '2016-06-23'

On the command-line:

.. code-block:: bash

    $ htmldate -u http://blog.python.org/2016/12/python-360-is-now-available.html
    '2016-12-23'


Features
--------


-  Compatible with all recent versions of Python (see above)
-  Multilingual, robust and efficient (used in production on millions of documents)
-  URLs, HTML files, or HTML trees are given as input (includes batch processing)
-  Output as string in any date format (defaults to `ISO 8601 YMD <https://en.wikipedia.org/wiki/ISO_8601>`_)
-  Detection of both original and updated dates


*htmldate* finds original and updated publication dates of web pages using heuristics on HTML code and linguistic patterns. It provides following ways to date a HTML document:

1. **Markup in header**: Common patterns are used to identify relevant elements (e.g. ``link`` and ``meta`` elements) including `Open Graph protocol <http://ogp.me/>`_ attributes and a large number of CMS idiosyncrasies
2. **HTML code**: The whole document is then searched for structural markers: ``abbr`` and ``time`` elements as well as a series of attributes (e.g. ``postmetadata``)
3. **Bare HTML content**: A series of heuristics is run on text and markup:

  - in ``fast`` mode the HTML page is cleaned and precise patterns are targeted
  - in ``extensive`` mode all potential dates are collected and a disambiguation algorithm determines the best one


Performance
-----------

=============================== ========= ========= ========= ========= =======
500 web pages containing identifiable dates (as of 2021-09-24)
-------------------------------------------------------------------------------
Python Package                  Precision Recall    Accuracy  F-Score   Time
=============================== ========= ========= ========= ========= =======
articleDateExtractor 0.20       0.769     0.691     0.572     0.728     3.3x
date_guesser 2.1.4              0.738     0.544     0.456     0.626     20x
goose3 3.1.9                    0.821     0.453     0.412     0.584     8.2x
htmldate[all] 0.9.1 (fast)      **0.839** 0.906     0.772     0.871     **1x**
htmldate[all] 0.9.1 (extensive) 0.825     **0.990** **0.818** **0.900** 1.7x
newspaper3k 0.2.8               0.729     0.630     0.510     0.675     8.4x
news-please 1.5.21              0.769     0.691     0.572     0.728     30x
=============================== ========= ========= ========= ========= =======

For complete results and explanations see the `evaluation page <https://htmldate.readthedocs.io/en/latest/evaluation.html>`_.


Installation
------------

This Python package is tested on Linux, macOS and Windows systems, it is compatible with Python 3.6 upwards. It is available on the package repository `PyPI <https://pypi.org/>`_ and can notably be installed with ``pip`` (``pip3`` where applicable): ``pip install htmldate`` and optionally ``pip install htmldate[speed]``.


Documentation
-------------

For more details on installation, Python & CLI usage, **please refer to the documentation**: `htmldate.readthedocs.io <https://htmldate.readthedocs.io/>`_


License
-------

*htmldate* is distributed under the `GNU General Public License v3.0 <https://github.com/adbar/htmldate/blob/master/LICENSE>`_. If you wish to redistribute this library but feel bounded by the license conditions please try interacting `at arms length <https://www.gnu.org/licenses/gpl-faq.html#GPLInProprietarySystem>`_, `multi-licensing <https://en.wikipedia.org/wiki/Multi-licensing>`_ with `compatible licenses <https://en.wikipedia.org/wiki/GNU_General_Public_License#Compatibility_and_multi-licensing>`_, or `contacting me <https://github.com/adbar/htmldate#author>`_.

See also `GPL and free software licensing: What's in it for business? <https://www.techrepublic.com/blog/cio-insights/gpl-and-free-software-licensing-whats-in-it-for-business/>`_


Author
------

This effort is part of methods to derive information from web documents in order to build `text databases for research <https://www.dwds.de/d/k-web>`_ (chiefly linguistic analysis and natural language processing). Extracting and pre-processing web texts to the exacting standards of scientific research presents a substantial challenge for those who conduct such research. There are web pages for which neither the URL nor the server response provide a reliable way to find out when a document was published or modified. For more information:

.. image:: https://joss.theoj.org/papers/10.21105/joss.02439/status.svg
   :target: https://doi.org/10.21105/joss.02439
   :alt: JOSS article

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3459599.svg
   :target: https://doi.org/10.5281/zenodo.3459599
   :alt: Zenodo archive

.. code-block:: shell

    @article{barbaresi-2020-htmldate,
      title = {{htmldate: A Python package to extract publication dates from web pages}},
      author = "Barbaresi, Adrien",
      journal = "Journal of Open Source Software",
      volume = 5,
      number = 51,
      pages = 2439,
      url = {https://doi.org/10.21105/joss.02439},
      publisher = {The Open Journal},
      year = 2020,
    }

-  Barbaresi, A. "`htmldate: A Python package to extract publication dates from web pages <https://doi.org/10.21105/joss.02439>`_", Journal of Open Source Software, 5(51), 2439, 2020. DOI: 10.21105/joss.02439
-  Barbaresi, A. "`Generic Web Content Extraction with Open-Source Software <https://hal.archives-ouvertes.fr/hal-02447264/document>`_", Proceedings of KONVENS 2019, Kaleidoscope Abstracts, 2019.
-  Barbaresi, A. "`Efficient construction of metadata-enhanced web corpora <https://hal.archives-ouvertes.fr/hal-01371704v2/document>`_", Proceedings of the `10th Web as Corpus Workshop (WAC-X) <https://www.sigwac.org.uk/wiki/WAC-X>`_, 2016.

You can contact me via my `contact page <https://adrien.barbaresi.eu/>`_ or `GitHub <https://github.com/adbar>`_.


Contributing
------------

`Contributions <https://github.com/adbar/htmldate/blob/master/CONTRIBUTING.md>`_ are welcome!

Feel free to file issues on the `dedicated page <https://github.com/adbar/htmldate/issues>`_. Thanks to the `contributors <https://github.com/adbar/htmldate/graphs/contributors>`_ who submitted features and bugfixes!

Kudos to the following software libraries:

-  `ciso8601 <https://github.com/closeio/ciso8601>`_, `lxml <http://lxml.de/>`_, `dateparser <https://github.com/scrapinghub/dateparser>`_
-  A few patterns are derived from the `python-goose <https://github.com/grangier/python-goose>`_, `metascraper <https://github.com/ianstormtaylor/metascraper>`_, `newspaper <https://github.com/codelucas/newspaper>`_ and `articleDateExtractor <https://github.com/Webhose/article-date-extractor>`_ libraries. This module extends their coverage and robustness significantly.
Evaluation
==========

Sources
-------

Date-annotated HTML pages
^^^^^^^^^^^^^^^^^^^^^^^^^

- BBAW collection (multilingual): Adrien Barbaresi, Lukas Kozmus.
- Additional English news pages: `Data Culture Group <https://dataculturegroup.org>`_ at Northeastern University.


Reproducing the evaluation
--------------------------

1. Install the packages specified in ``eval-requirements.txt``
2. Run the script ``comparison.py``
Evaluation
==========


Although text is ubiquitous on the Web, extracting information from web pages can prove to be difficult. In most cases, immediately accessible data on retrieved webpages do not carry substantial or accurate information: neither the URL nor the server response provide a reliable way to date a web document, that is find when it was written or modified. Content extraction mostly draws on Document Object Model (DOM) examination, that is on considering a given HTML document as a tree structure whose nodes represent parts of the document to be operated on. Less thorough and not necessarily faster alternatives use superficial search patterns such as regular expressions in order to capture desirable excerpts.


Alternatives
------------

There are comparable software solutions in Python, the following date extraction packages are open-source and work out-of-the-box:

- `articleDateExtractor <https://github.com/Webhose/article-date-extractor>`_ detects, extracts and normalizes the publication date of an online article or blog post,
- `date_guesser <https://github.com/mitmedialab/date_guesser>`_ extracts publication dates from a web pages along with an accuracy measure (not used here),
- `goose3 <https://github.com/goose3/goose3>`_ can extract information for embedded content,
- `htmldate <https://github.com/adbar/htmldate>`_ is the software package described here, it is designed to extract original and updated publication dates of web pages,
- `newspaper <https://github.com/codelucas/newspaper>`_ is mostly geared towards newspaper texts,
- `news-please <https://github.com/fhamborg/news-please>`_ is a news crawler that extracts structured information.

Two alternative packages are not tested here but could be used in addition:

- `datefinder <https://github.com/akoumjian/datefinder>`_ features pattern-based date extraction for texts written in English,
- if the date is nowhere to be found `carbon dating <https://github.com/oduwsdl/CarbonDate>`_ the web page can be an option, however this is computationally expensive.


Description
-----------

**Test set**: the experiments below are run on a collection of documents which are either typical for Internet articles (news outlets, blogs, including smaller ones) or non-standard and thus harder to process. They were selected from `large collections of web pages in German <https://www.dwds.de/d/k-web>`_. For the sake of completeness a few documents in other languages were added (mostly in English and French but also in other European languages, Chinese, Japanese and Arabic).

**Evaluation**: only documents with dates that are clearly to be determined are considered for this benchmark. A given day is taken as unit of reference, meaning that results are converted to ``%Y-%m-%d`` format if necessary in order to make them comparable. The evaluation script is available on the project repository: `tests/comparison.py <https://github.com/adbar/htmldate/blob/master/tests/comparison.py>`_. To reproduce the tests just clone the repository, install all necessary packages and run the evaluation script with the data provided in the *tests* directory.

**Time**: the execution time (best of 3 tests) cannot be easily compared in all cases as some solutions perform a whole series of operations which are irrelevant to this task.

**Errors:** *goose3*'s output isn't always meaningful and/or in a standardized format, these cases were discarded. *news-please* seems to have trouble with some encodings (e.g. in Chinese), in which case it leads to an exception.


Results
-------

The results below show that **date extraction is not a completely solved task** but one for which extractors have to resort to heuristics and guesses. The figures documenting recall and accuracy capture the real-world performance of the tools as the absence of a date output impacts the result.


=============================== ========= ========= ========= ========= =======
500 web pages containing identifiable dates (as of 2021-09-24)
-------------------------------------------------------------------------------
Python Package                  Precision Recall    Accuracy  F-Score   Time
=============================== ========= ========= ========= ========= =======
articleDateExtractor 0.20       0.769     0.691     0.572     0.728     3.3x
date_guesser 2.1.4              0.738     0.544     0.456     0.626     20x
goose3 3.1.9                    0.821     0.453     0.412     0.584     8.2x
htmldate[all] 0.9.1 (fast)      **0.839** 0.906     0.772     0.871     **1x**
htmldate[all] 0.9.1 (extensive) 0.825     **0.990** **0.818** **0.900** 1.7x
newspaper3k 0.2.8               0.729     0.630     0.510     0.675     8.4x
news-please 1.5.21              0.769     0.691     0.572     0.728     30x
=============================== ========= ========= ========= ========= =======


Additional data for new pages in English collected by the `Data Culture Group <https://dataculturegroup.org>`_ at Northeastern University.

Precision describes if the dates given as output are correct: *goose3* fares well precision-wise but it fails to extract dates in a large majority of cases (poor recall). The difference in accuracy between *date_guesser* and *newspaper* is consistent with tests described on the `website of the former <https://github.com/mitmedialab/date_guesser>`_.

It turns out that *htmldate* performs better than the other solutions overall. It is also noticeably faster than the strictly comparable packages (*articleDateExtractor* and most certainly *date_guesser*). Despite being measured on a sample, **the higher accuracy and faster processing time are highly significant**. Especially for smaller news outlets, websites and blogs, as well as pages written in languages other than English (in this case mostly but not exclusively German), *htmldate* greatly extends date extraction coverage without sacrificing precision.


Note on the different versions:

- *htmldate[all]* means that additional components are added for performance and coverage, which results in differences with respect to accuracy (due to further linguistic analysis) and potentially speed (faster date parsing). They can be installed with ``pip/pip3/pipenv htmldate[all]``.
- The fast mode does not output as many dates (lower recall) but its guesses are more often correct (better precision).


Older Results
-------------


=============================== ========= ========= ========= ========= =======
225 web pages containing identifiable dates (as of 2020-07-29)
-------------------------------------------------------------------------------
Python Package                  Precision Recall    Accuracy  F-Score   Time
=============================== ========= ========= ========= ========= =======
articleDateExtractor 0.20       0.817     0.635     0.556     0.714     6.8
date_guesser 2.1.4              0.809     0.553     0.489     0.657     40.0
goose3 3.1.6                    0.887     0.441     0.418     0.589     15.5
htmldate 0.7.0 (fast)           **0.903** 0.907     0.827     0.905     **2.4**
htmldate[all] 0.7.0 (extensive) 0.889     **1.000** **0.889** **0.941** 3.8
newspaper 0.2.8                 0.888     0.407     0.387     0.558     81.6
news-please 1.5.3               0.823     0.660     0.578     0.732     69.6
=============================== ========= ========= ========= ========= =======


=============================== ========= ========= ========= ========= =======
225 web pages containing identifiable dates (as of 2020-11-03)
-------------------------------------------------------------------------------
Python Package                  Precision Recall    Accuracy  F-Score   Time
=============================== ========= ========= ========= ========= =======
articleDateExtractor 0.20       0.817     0.635     0.556     0.714     3.5x
date_guesser 2.1.4              0.809     0.553     0.489     0.657     21x
goose3 3.1.6                    0.887     0.441     0.418     0.589     7.7x
htmldate[all] 0.7.2 (fast)      **0.899** 0.917     0.831     0.908     **1x**
htmldate[all] 0.7.2 (extensive) 0.893     **1.000** **0.893** **0.944** 1.6x
newspaper3k 0.2.8               0.888     0.407     0.387     0.558     40x
news-please 1.5.13              0.823     0.660     0.578     0.732     31x
=============================== ========= ========= ========= ========= =======
Core functions
==============

.. contents:: **Contents**
    :backlinks: none


Handling date extraction
------------------------

.. autofunction:: htmldate.core.find_date

.. autofunction:: htmldate.core.examine_header

.. autofunction:: htmldate.core.search_page


Useful internal functions
-------------------------

.. autofunction:: htmldate.extractors.try_ymd_date

.. autofunction:: htmldate.extractors.custom_parse

.. autofunction:: htmldate.extractors.regex_parse

.. autofunction:: htmldate.extractors.extract_url_date

.. autofunction:: htmldate.extractors.extract_partial_url_date

.. autofunction:: htmldate.extractors.external_date_parser


Helpers
-------

.. autofunction:: htmldate.extractors.convert_date

.. autofunction:: htmldate.extractors.date_validator

.. autofunction:: htmldate.utils.load_html

.. autofunction:: htmldate.utils.fetch_url
Options
=======

.. contents:: **Contents**
    :backlinks: none


Configuration
-------------


Input format
~~~~~~~~~~~~

The module expects strings as shown above, it is also possible to use already parsed HTML (i.e. a LXML tree object):

.. code-block:: python

    >>> from htmldate import find_date
    >>> from lxml import html
    >>> mytree = html.fromstring('<html><body><span class="entry-date">July 12th, 2016</span></body></html>')
    >>> find_date(mytree)
    '2016-07-12'

An external module can be used for download, as described in versions anterior to 0.3. This example uses the legacy mode with `requests <http://docs.python-requests.org/>`_ as external module.

.. code-block:: python

    >>> from htmldate.core import find_date
    # using requests
    >>> import requests
    >>> r = requests.get('https://creativecommons.org/about/')
    >>> find_date(r.text)
    '2017-11-28' # may have changed since
    # using htmldate's own fetch_url function
    >>> from htmldate.utils import fetch_url
    >>> htmldoc = fetch_url('https://blog.wikimedia.org/2018/06/28/interactive-maps-now-in-your-language/')
    >>> find_date(htmldoc)
    '2018-06-28'
    # or simply
    >>> find_date('https://blog.wikimedia.org/2018/06/28/interactive-maps-now-in-your-language/') # URL detected
    '2018-06-28'


Date format
~~~~~~~~~~~

The output format of the dates found can be set in a format known to Python's ``datetime`` module, the default being ``%Y-%m-%d``:

.. code-block:: python

    >>> find_date('https://www.gnu.org/licenses/gpl-3.0.en.html', outputformat='%d %B %Y')
    '18 November 2016' # may have changed since


.. autofunction:: htmldate.validators.output_format_validator


Original date
~~~~~~~~~~~~~

Although the time delta between the original publication and the "last modified" statement is usually a matter of hours or days at most, it can be useful in some contexts to prioritize the original publication date during extraction:

.. code-block:: python

    >>> find_date('https://netzpolitik.org/2016/die-cider-connection-abmahnungen-gegen-nutzer-von-creative-commons-bildern/') # default setting
    '2019-06-24'
    >>> find_date('https://netzpolitik.org/2016/die-cider-connection-abmahnungen-gegen-nutzer-von-creative-commons-bildern/', original_date=True) # modified behavior
    '2016-06-23'


Settings
--------

See ``settings.py`` file:

.. automodule:: htmldate.settings
   :members:
   :show-inheritance:
   :undoc-members:

The module can then be re-compiled locally to apply changes to the settings.
htmldate: find the publication date of web pages
================================================

.. image:: https://img.shields.io/pypi/v/htmldate.svg
    :target: https://pypi.python.org/pypi/htmldate
    :alt: Python package

.. image:: https://img.shields.io/pypi/pyversions/htmldate.svg
    :target: https://pypi.python.org/pypi/htmldate
    :alt: Python versions

.. image:: https://img.shields.io/codecov/c/github/adbar/htmldate.svg
    :target: https://codecov.io/gh/adbar/htmldate
    :alt: Code Coverage

.. image:: https://static.pepy.tech/badge/htmldate/month
    :target: https://pepy.tech/project/htmldate
    :alt: Downloads

|

:Code:           https://github.com/adbar/htmldate
:Documentation:  https://htmldate.readthedocs.io
:Issue tracker:  https://github.com/adbar/htmldate/issues

|

.. image:: htmldate-demo.gif
    :alt: Demo as GIF image
    :align: center
    :width: 80%
    :target: https://htmldate.readthedocs.org/

|

Find original and updated publication dates of any web page. From the command-line or within Python, all the steps needed from web page download to HTML parsing, scraping, and text analysis are included.

In a nutshell, with Python:

.. code-block:: python

    >>> from htmldate import find_date
    >>> find_date('http://blog.python.org/2016/12/python-360-is-now-available.html')
    '2016-12-23'
    >>> find_date('https://netzpolitik.org/2016/die-cider-connection-abmahnungen-gegen-nutzer-von-creative-commons-bildern/', original_date=True)
    '2016-06-23'

On the command-line:

.. code-block:: bash

    $ htmldate -u http://blog.python.org/2016/12/python-360-is-now-available.html
    '2016-12-23'

|

.. contents:: **Contents**
    :backlinks: none

|

Features
--------


-  Compatible with all recent versions of Python (see above)
-  Multilingual, robust and efficient (used in production on millions of documents)
-  URLs, HTML files, or HTML trees are given as input (includes batch processing)
-  Output as string in any date format (defaults to `ISO 8601 YMD <https://en.wikipedia.org/wiki/ISO_8601>`_)
-  Detection of both original and updated dates


*htmldate* finds original and updated publication dates of web pages using heuristics on HTML code and linguistic patterns. It provides following ways to date a HTML document:

1. **Markup in header**: Common patterns are used to identify relevant elements (e.g. ``link`` and ``meta`` elements) including `Open Graph protocol <http://ogp.me/>`_ attributes and a large number of CMS idiosyncrasies
2. **HTML code**: The whole document is then searched for structural markers: ``abbr`` and ``time`` elements as well as a series of attributes (e.g. ``postmetadata``)
3. **Bare HTML content**: A series of heuristics is run on text and markup:

  - in ``fast`` mode the HTML page is cleaned and precise patterns are targeted
  - in ``extensive`` mode all potential dates are collected and a disambiguation algorithm determines the best one

The output is thouroughly verified in terms of plausibility and adequateness and the library outputs a date string, corresponding to either the last update or the original publishing statement (the default), in the desired format (defaults to `ISO 8601 YMD format <https://en.wikipedia.org/wiki/ISO_8601>`_).

Markup-based extraction is multilingual by nature, text-based refinements for better coverage currently support German, English and Turkish.


Installation
------------

This Python package is tested on Linux, macOS and Windows systems, it is compatible with Python 3.6 upwards. It is available on the package repository `PyPI <https://pypi.org/>`_ and can notably be installed with ``pip`` or ``pipenv``:

.. code-block:: bash

    $ pip install htmldate # pip3 install on systems where both Python 2 and 3 are installed
    $ pip install --upgrade htmldate # to make sure you have the latest version
    $ pip install git+https://github.com/adbar/htmldate.git # latest available code (see build status above)

Additional libraries can be installed to enhance efficiency: ``cchardet`` and ``ciso8601`` (for speed). They may not work on all platforms and have thus been singled out although installation is recommended:

.. code-block:: bash

    $ pip install htmldate[speed] # install with additional functionality

You can also install or update the packages separately, *htmldate* will detect which ones are present on your system and opt for the best available combination.

*For infos on dependency management of Python packages see* `this discussion thread <https://stackoverflow.com/questions/41573587/what-is-the-difference-between-venv-pyvenv-pyenv-virtualenv-virtualenvwrappe>`_.


With Python
-----------

All the functions of the module are currently bundled in *htmldate*.

In case the web page features easily readable metadata in the header, the extraction is straightforward. A more advanced analysis of the document structure is sometimes needed:

.. code-block:: python

    >>> from htmldate import find_date
    >>> find_date('http://blog.python.org/2016/12/python-360-is-now-available.html')
    '# DEBUG analyzing: <h2 class="date-header"><span>Friday, December 23, 2016</span></h2>'
    '# DEBUG result: 2016-12-23'
    '2016-12-23'

``htmldate`` can resort to a guess based on a complete screening of the document (``extensive_search`` parameter) which can be deactivated:

.. code-block:: python

    >>> find_date('https://creativecommons.org/about/')
    '2017-08-11' # has been updated since
    >>> find_date('https://creativecommons.org/about/', extensive_search=False)
    >>>

Already parsed HTML (that is a LXML tree object):

.. code-block:: python

    # simple HTML document as string
    >>> htmldoc = '<html><body><span class="entry-date">July 12th, 2016</span></body></html>'
    >>> find_date(htmldoc)
    '2016-07-12'
    # parsed LXML tree
    >>> from lxml import html
    >>> mytree = html.fromstring('<html><body><span class="entry-date">July 12th, 2016</span></body></html>')
    >>> find_date(mytree)
    '2016-07-12'

Change the output to a format known to Python's ``datetime`` module, the default being ``%Y-%m-%d``:

.. code-block:: python

    >>> find_date('https://www.gnu.org/licenses/gpl-3.0.en.html', outputformat='%d %B %Y')
    '18 November 2016'  # may have changed since
    >>> find_date('http://blog.python.org/2016/12/python-360-is-now-available.html', outputformat='%Y-%m-%dT%H:%M:%S%z')
    '2016-12-23T05:11:00-0500'

Although the time delta between original publication and "last modified" info is usually a matter of hours or days, it can be useful to prioritize the **original publication date**:

.. code-block:: python

    >>> find_date('https://netzpolitik.org/2016/die-cider-connection-abmahnungen-gegen-nutzer-von-creative-commons-bildern/', original_date=True)  # modified behavior
    '2016-06-23'

For more information see `options page <options.html>`_.


On the command-line
-------------------

A command-line interface is included:

.. code-block:: bash

    $ htmldate -u http://blog.python.org/2016/12/python-360-is-now-available.html
    '2016-12-23'
    $ wget -qO- "http://blog.python.org/2016/12/python-360-is-now-available.html" | htmldate
    '2016-12-23'

For usage instructions see ``htmldate -h``:

.. code-block:: bash

    $ htmldate --help
    htmldate [-h] [-f] [-i INPUTFILE] [--original] [-min MINDATE] [-max MAXDATE] [-u URL] [-v] [--version]
    optional arguments:
        -h, --help            show this help message and exit
        -f, --fast            fast mode: disable extensive search
        -i INPUTFILE, --inputfile INPUTFILE
                              name of input file for batch processing (similar to wget -i)
        --original            original date prioritized
        -min MINDATE, --mindate MINDATE
                              earliest acceptable date (YYYY-MM-DD)
        -max MAXDATE, --maxdate MAXDATE
                              latest acceptable date (YYYY-MM-DD)
        -u URL, --URL URL     custom URL download
        -v, --verbose         increase output verbosity
        --version             show version information and exit


The batch mode ``-i`` takes one URL per line as input and returns one result per line in tab-separated format:

.. code-block:: bash

    $ htmldate --fast -i list-of-urls.txt


License
-------

*htmldate* is distributed under the `GNU General Public License v3.0 <https://github.com/adbar/htmldate/blob/master/LICENSE>`_. If you wish to redistribute this library but feel bounded by the license conditions please try interacting `at arms length <https://www.gnu.org/licenses/gpl-faq.html#GPLInProprietarySystem>`_, `multi-licensing <https://en.wikipedia.org/wiki/Multi-licensing>`_ with `compatible licenses <https://en.wikipedia.org/wiki/GNU_General_Public_License#Compatibility_and_multi-licensing>`_, or `contacting me <https://github.com/adbar/htmldate#author>`_.

See also `GPL and free software licensing: What's in it for business? <https://www.techrepublic.com/blog/cio-insights/gpl-and-free-software-licensing-whats-in-it-for-business/>`_


Author
------

This effort is part of methods to derive information from web documents in order to build `text databases for research <https://www.dwds.de/d/k-web>`_ (chiefly linguistic analysis and natural language processing). Extracting and pre-processing web texts to the exacting standards of scientific research presents a substantial challenge for those who conduct such research. There are web pages for which neither the URL nor the server response provide a reliable way to find out when a document was published or modified. For more information:

.. image:: https://joss.theoj.org/papers/10.21105/joss.02439/status.svg
   :target: https://doi.org/10.21105/joss.02439
   :alt: JOSS article

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3459599.svg
   :target: https://doi.org/10.5281/zenodo.3459599
   :alt: Zenodo archive

.. code-block:: shell

    @article{barbaresi-2020-htmldate,
      title = {{htmldate: A Python package to extract publication dates from web pages}},
      author = "Barbaresi, Adrien",
      journal = "Journal of Open Source Software",
      volume = 5,
      number = 51,
      pages = 2439,
      url = {https://doi.org/10.21105/joss.02439},
      publisher = {The Open Journal},
      year = 2020,
    }

-  Barbaresi, A. "`htmldate: A Python package to extract publication dates from web pages <https://doi.org/10.21105/joss.02439>`_", Journal of Open Source Software, 5(51), 2439, 2020. DOI: 10.21105/joss.02439
-  Barbaresi, A. "`Generic Web Content Extraction with Open-Source Software <https://hal.archives-ouvertes.fr/hal-02447264/document>`_", Proceedings of KONVENS 2019, Kaleidoscope Abstracts, 2019.
-  Barbaresi, A. "`Efficient construction of metadata-enhanced web corpora <https://hal.archives-ouvertes.fr/hal-01371704v2/document>`_", Proceedings of the `10th Web as Corpus Workshop (WAC-X) <https://www.sigwac.org.uk/wiki/WAC-X>`_, 2016.

You can contact me via my `contact page <https://adrien.barbaresi.eu/>`_ or `GitHub <https://github.com/adbar>`_.


Contributing
------------

`Contributions <https://github.com/adbar/htmldate/blob/master/CONTRIBUTING.md>`_ are welcome!

Feel free to file issues on the `dedicated page <https://github.com/adbar/htmldate/issues>`_. Thanks to the `contributors <https://github.com/adbar/htmldate/graphs/contributors>`_ who submitted features and bugfixes!

Kudos to the following software libraries:

-  `ciso8601 <https://github.com/closeio/ciso8601>`_, `lxml <http://lxml.de/>`_, `dateparser <https://github.com/scrapinghub/dateparser>`_
-  A few patterns are derived from the `python-goose <https://github.com/grangier/python-goose>`_, `metascraper <https://github.com/ianstormtaylor/metascraper>`_, `newspaper <https://github.com/codelucas/newspaper>`_ and `articleDateExtractor <https://github.com/Webhose/article-date-extractor>`_ libraries. This module extends their coverage and robustness significantly.


Going further
-------------

Known caveats
~~~~~~~~~~~~~

The granularity may not always match the desired output format. If only information about the year could be found and the chosen date format requires to output a month and a day, the result is 'padded' to be located at the middle of the year, in that case the 1st of January.

Besides, there are pages for which no date can be found, ever:

.. code-block:: python

    >>> r = requests.get('https://example.com')
    >>> htmldate.find_date(r.text)
    >>>

If the date is nowhere to be found, it might be worth considering `carbon dating <https://github.com/oduwsdl/CarbonDate>`_ the web page, however this is computationally expensive. In addition, `datefinder <https://github.com/akoumjian/datefinder>`_ features pattern-based date extraction for texts written in English.

Tests
~~~~~

A series of webpages triggering different structural and content patterns is included for testing purposes:

.. code-block:: bash

    $ pytest tests/unit_tests.py


.. toctree::
   :maxdepth: 2

   corefunctions
   evaluation
   options


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
