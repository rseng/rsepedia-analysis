---
title: 'TRUNAJOD: A text complexity library to enhance natural language processing'
tags:
  - Python
  - natural language processing
  - machine learning
  - text complexity
  - text coherence
authors:
  - name: Diego A. Palma
    orcid: 0000-0003-1540-7164
    affiliation: 1
  - name: Christian Soto
    affiliation: 1
  - name: Mónica Veliz
    affiliation: 1
  - name: Bruno Karelovic
    affiliation: 1
  - name: Bernardo Riffo
    affiliation: 1
affiliations:
 - name: Universidad de Concepción
   index: 1
date: 3 March 2020
bibliography: paper.bib
---

# Summary

We present `TRUNAJOD`, a text complexity analysis tool that includes a wide variety of linguistics measurements that can be extracted from texts as an approximation for readability, coherence, and cohesion. The features that `TRUNAJOD` can extract from the text are based on the literature and can be separated into the following categories: discourse markers, emotions, entity grid-based measurements, givenness, lexical-semantic norms, semantic measures, surface proxies, etc. In this first version of `TRUNAJOD`, we mainly support the Spanish language, but several features support any language that has proper natural language processing POS tagging and dependency parsing capabilities. Finally, we show how TRUNAJOD could be used in applied research.

# Statement of need

`TRUNAJOD` aims to address three challenges:

1. A standardized API for text complexity measurements
2. An open-source code, so any researcher in the linguistics field could contribute to it
3. Easy-to-build applications and tools that rely upon text complexity assessment

Other tools aim to make it easy for the public to get coherence and cohesion metrics. One such tool is TAACO [@crossley2019tool], which is written in Python and can be freely downloaded. A problem with TAACO is that it is a desktop application, which encloses the code. This makes it impossible to contribute modifications or new features, as it is a closed system. Moreover, it does not implement other relevant features to assess cohesion and coherence of discourse, for example, entity grid-based features. One open-source project with this purpose is the `Cohere` framework [@smith2016cohere], which is written in a mixture of Java and Python. However, it does not seem to be actively maintained and it does not implement other measurements that could be used by other researchers. On the other hand, most of the tools only support English languages and do not provide support for a plethora of metrics available in a comprehensible API. `TRUNAJOD` aims to be different, in the sense that we do not present a closed system, but rather, an open-source project, trying to follow the best Python development patterns. Furthermore, we rely on `spaCy`, enabling us to support not only one language but multiple languages for coherence and cohesion tasks, which enables `TRUNAJOD` to improve its performance when `spaCy` does, promoting collaboration.

Moreover, `TRUNAJOD` not only implements state-of-the-art measurements for text complexity assessment, but also bundles new sets of predictors for this task. In this sense, `TRUNAJOD` ́s contributions are:

* Fixing periphrasis of texts, because many NLP tools have issues dealing with periphrasis. In this release, this only applies to Spanish.
* Adding heuristics for measurements based on clause count. TRUNAJOD provides a new algorithm for clause segmentation.
* TRUNAJOD provides several approximations to narrativity in these new clause segmentation-dependent indices.

Text complexity assessment is a natural language processing task that can be applied to multiple problems, such as automatic summarization, automatic essay scoring, automatic summary evaluation, intelligent tutoring systems, and so on. Text complexity is usually related to the readability of a text, which is dependent on several of its intrinsic properties, mainly cohesion and coherence.

Automatic coherence evaluation is an open problem, and there have been several studies addressing it. On one side, text coherence assessment has been related to how sentences connect either semantically, or by co-referencing noun phrases. In the semantic view of coherence, Latent Semantic Analysis (LSA) [@foltz1998measurement] has been widely used because of its simplicity. In essence, sentences are represented as vectors, and the coherence of the text is computed as the average sentence similarity, using similarity vector measurements (such as cosine distance). This approach has drawbacks, such as that sentence ordering does not matter. To solve this issue, other methodologies have been proposed based on discourse theory, in particular the centering theory. One such approach is entity grids [@barzilay2008modeling] and entity graphs [@guinaudeau2013graph] that treat coherence as to how are entities take different roles between sentences and how are they connected in the text. `TRUNAJOD` implements all these models, and thus `TRUNAJOD` can compute coherence based on sentence similarities using word vectors. `TRUNAJOD` also provides an API for dealing with entity grids and entity graphs, to extract such measurements.

A downside of the previous approaches for text complexity is that they only capture either CORPUS-based semantics or relationships between entities. Additionally, entity grids rely heavily on the dependency parser at hand and the co-reference resolution used because an entity might be mentioned in several ways across a text. The problem with this is that these measurements might be noisy depending on the use case, and simpler measurements would fit better in such cases. In `TRUNAJOD`, we compute several surface proxies that have been used by several state-of-the-art text assessment tools [@mcnamara2014automated] [@page1994computer]. Such surface proxies try to approximate intrinsic properties of the text such as narrativity, connectiveness, givenness, cohesion, and coherence. `TRUNAJOD` includes classical measurements such as word count, sentence count, pronoun-noun ratio, type-token ratios, frequency index, etc. Moreover, `TRUNAJOD` comes bundled with heuristics to compute clause count-based metrics, such as subordination and clause length, among others. To achieve this, `TRUNAJOD` adds periphrasis tags to the text to heuristically segment clauses.

One drawback of using surface proxies (shallow measurements) is that they do not capture all the properties of the text, and just rely on approximations that are captured from raw text. In some use cases, this is not desirable (e.g., automatic essay scoring, intelligent tutoring systems), and these measurements should be complemented with other measurements that are desirable in those use cases, such as lexical-semantic norms and word emotions. Lexical semantic norms are norms for words that are related to the psychological degree of activation of a word in the reader [@guasch2016spanish]. Some examples of such variables are concreteness, imageability, familiarity, arousal, valence. These variables might be used for reading comprehension tools (e.g., in reading comprehension assessments, it is desirable that feedback is concrete). Moreover, emotions could be used in such cases, and even in opinion mining [@sidorov2012empirical]. `TRUNAJOD` comes bundled with both types of features, and thus, average lexical-semantic norms and emotions could be extracted from the text. 

`TRUNAJOD` architectural design is shown in Figure 1.

![`TRUNAJOD` architectural design.](imgs/figure1.png){ width=50% }

Basically, `TRUNAJOD` API takes as input a spaCy Doc and `TRUNAJOD` models (lemmatizer, synonym map, lexical-semantic norm map, etc.) and then it can compute supported text complexity metrics. It is worth noting that `TRUNAJOD` has available downloadable models from its GitHub repository, but currently only Spanish models are available. Nevertheless, it should be straightforward adding models for different languages.

# Acknowledgements

This research was supported by FONDEF (Chile) under Grant IT17I0051 "Desarrollo de una herramienta computacional para 
evaluación automática de textos del Sistema escolar chileno" ("Development of a computational tool for automatic 
assessment of Chilean school texts")

# References
# Changelog


## v0.1.1

* Add `spaCy` as a required dependency (#23).
* Fix pypi `long_description` (#23).
* Migrate from `travis-ci` to `Github` actions (#22).
* Add Type Hints to TTR module (#18)
* Add Type Hints to discourse Markers (#16)
<p align="center">
<img style="width: 30%; height: 30%" src="https://raw.githubusercontent.com/dpalmasan/TRUNAJOD2.0/master/imgs/trunajod_logo.png">
</p>

# TRUNAJOD: A text complexity library for text analysis built on spaCy

<p align="center">
<a href="https://github.com/dpalmasan/TRUNAJOD2.0/actions"><img alt="Actions Status" src="https://github.com/dpalmasan/TRUNAJOD2.0/workflows/Test/badge.svg"></a>
<a href="https://trunajod20.readthedocs.io/en/stable/?badge=stable"><img alt="Documentation Status" src="https://readthedocs.org/projects/trunajod20/badge/?version=stable"></a>
<img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/trunajod">
<a href="https://github.com/dpalmasan/TRUNAJOD2.0/blob/master/LICENSE"><img alt="License: MIT" src="https://img.shields.io/github/license/dpalmasan/TRUNAJOD2.0"></a>
<a href="https://pypi.org/project/TRUNAJOD/"><img alt="PyPI" src="https://img.shields.io/pypi/v/TRUNAJOD"></a>
<a href="https://pepy.tech/project/trunajod"><img alt="Downloads" src="https://static.pepy.tech/badge/TRUNAJOD"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
<a href="https://spacy.io"><img alt="Built with spaCy" src="https://img.shields.io/badge/built%20with-spaCy-09a3d5.svg"></a>
<a href="https://doi.org/10.21105/joss.03153"><img alt="JOSS paper" src="https://joss.theoj.org/papers/10.21105/joss.03153/status.svg"></a>
<a href="https://doi.org/10.5281/zenodo.4707403"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4707403.svg" alt="DOI"></a>
</p>

``TRUNAJOD`` is a Python library for text complexity analysis build on the 
high-performance [spaCy](https://spacy.io/) library. With all the basic NLP capabilities provided by
spaCy (dependency parsing, POS tagging, tokenizing), ``TRUNAJOD`` focuses on extracting
measurements from texts that might be interesting for different applications and use cases.
While most of the indices could be computed for different languages, currently we mostly support 
Spanish. We are happy if you contribute with indices implemented for your language!

## Features

* Utilities for text processing such as lemmatization, POS checkings.
* Semantic measurements from text such as average coherence between sentences and average synonym overlap.
* Giveness measurements such as pronoun density and pronoun noun ratio.
* Built-in emotion lexicon to compute emotion calculations based on words in the text.
* Lexico-semantic norm dataset to compute lexico-semantic variables from text.
* Type token ratio (TTR) based metrics, and tunnable TTR metrics.
* A built-in syllabizer (currently only for spanish).
* Discourse markers based measurements to obtain measures of connectivity inside the text.
* Plenty of surface proxies of text readability that can be computed directly from text.
* Measurements of parse tree similarity as an approximation to syntactic complexity.
* Parse tree correction to add periphrasis and heuristics for clause count, all based on linguistics experience.
* Entity Grid and entity graphs model implementation as a measure of coherence.
* An easy to use and user-friendly API.

## Installation

`TRUNAJOD` can be installed by running `pip install trunajod`. It requires Python 3.6.2+ to run.

## Getting Started

Using this package has some other pre-requisites. It assumes that you already have your model set up on spacy. If not, please first install or download a model (for Spanish users, a spanish model). Then you can get started with the following code snippet.

You can download pre-build `TRUNAJOD` models from the repo, under the `models` directory.

Below is a small snippet of code that can help you in getting started with this lib. Don´t forget to take a look at the [documentation](https://trunajod20.readthedocs.io/en/latest).

The example below assumes you have the `es_core_news_sm` spaCy Spanish model installed. You can install the model running: `python -m spacy download es_core_news_sm`. For other models, please check [spaCy docs](https://spacy.io/usage/models).


```python
from TRUNAJOD import surface_proxies
from TRUNAJOD.entity_grid import EntityGrid
from TRUNAJOD.lexico_semantic_norms import LexicoSemanticNorm
import pickle
import spacy
import tarfile


class ModelLoader(object):
    """Class to load model."""
    def __init__(self, model_file):
        tar = tarfile.open(model_file, "r:gz")
        self.crea_frequency = {}
        self.infinitive_map = {}
        self.lemmatizer = {}
        self.spanish_lexicosemantic_norms = {}
        self.stopwords = {}
        self.wordnet_noun_synsets = {}
        self.wordnet_verb_synsets = {}

        for member in tar.getmembers():
            f = tar.extractfile(member)
            if "crea_frequency" in member.name:
                self.crea_frequency = pickle.loads(f.read())
            if "infinitive_map" in member.name:
                self.infinitive_map = pickle.loads(f.read())
            if "lemmatizer" in member.name:
                self.lemmatizer = pickle.loads(f.read())
            if "spanish_lexicosemantic_norms" in member.name:
                self.spanish_lexicosemantic_norms = pickle.loads(f.read())
            if "stopwords" in member.name:
                self.stopwords = pickle.loads(f.read())
            if "wordnet_noun_synsets" in member.name:
                self.wordnet_noun_synsets = pickle.loads(f.read())
            if "wordnet_verb_synsets" in member.name:
                self.wordnet_verb_synsets = pickle.loads(f.read())


# Load TRUNAJOD models
model = ModelLoader("trunajod_models_v0.1.tar.gz")

# Load spaCy model
nlp = spacy.load("es_core_news_sm", disable=["ner", "textcat"])

example_text = (
    "El espectáculo del cielo nocturno cautiva la mirada y suscita preguntas"
    "sobre el universo, su origen y su funcionamiento. No es sorprendente que "
    "todas las civilizaciones y culturas hayan formado sus propias "
    "cosmologías. Unas relatan, por ejemplo, que el universo ha"
    "sido siempre tal como es, con ciclos que inmutablemente se repiten; "
    "otras explican que este universo ha tenido un principio, "
    "que ha aparecido por obra creadora de una divinidad."
)

doc = nlp(example_text)

# Lexico-semantic norms
lexico_semantic_norms = LexicoSemanticNorm(
    doc,
    model.spanish_lexicosemantic_norms,
    model.lemmatizer
)

# Frequency index
freq_index = surface_proxies.frequency_index(doc, model.crea_frequency)

# Clause count (heurístically)
clause_count = surface_proxies.clause_count(doc, model.infinitive_map)

# Compute Entity Grid
egrid = EntityGrid(doc)

print("Concreteness: {}".format(lexico_semantic_norms.get_concreteness()))
print("Frequency Index: {}".format(freq_index))
print("Clause count: {}".format(clause_count))
print("Entity grid:")
print(egrid.get_egrid())
```

This should output:

```
Concreteness: 1.95
Frequency Index: -0.7684649336888104
Clause count: 10
Entity grid:
{'ESPECTÁCULO': ['S', '-', '-'], 'CIELO': ['X', '-', '-'], 'MIRADA': ['O', '-', '-'], 'UNIVERSO': ['O', '-', 'S'], 'ORIGEN': ['X', '-', '-'], 'FUNCIONAMIENTO': ['X', '-', '-'], 'CIVILIZACIONES': ['-', 'S', '-'], 'CULTURAS': ['-', 'X', '-'], 'COSMOLOGÍAS': ['-', 'O', '-'], 'EJEMPLO': ['-', '-', 'X'], 'TAL': ['-', '-', 'X'], 'CICLOS': ['-', '-', 'X'], 'QUE': ['-', '-', 'S'], 'SE': ['-', '-', 'O'], 'OTRAS': ['-', '-', 'S'], 'PRINCIPIO': ['-', '-', 'O'], 'OBRA': ['-', '-', 'X'], 'DIVINIDAD': ['-', '-', 'X']}
```

## A real world example

`TRUNAJOD` lib was used to make `TRUNAJOD` web app, which is an application to assess text complexity and to check the adquacy of a text to a particular school level. To achieve this, several `TRUNAJOD` indices were analyzed for multiple Chilean school system texts (from textbooks), and latent features were created. Here is a snippet:

```python
"""Example of TRUNAJOD usage."""
import glob

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import spacy
import textract  # To read .docx files
import TRUNAJOD.givenness
import TRUNAJOD.ttr
from TRUNAJOD import surface_proxies
from TRUNAJOD.syllabizer import Syllabizer

plt.rcParams["figure.figsize"] = (11, 4)
plt.rcParams["figure.dpi"] = 200


nlp = spacy.load("es_core_news_sm", disable=["ner", "textcat"])

features = {
    "lexical_diversity_mltd": [],
    "lexical_density": [],
    "pos_dissimilarity": [],
    "connection_words_ratio": [],
    "grade": [],
}
for filename in glob.glob("corpus/*/*.docx"):
    text = textract.process(filename).decode("utf8")
    doc = nlp(text)
    features["lexical_diversity_mltd"].append(
        TRUNAJOD.ttr.lexical_diversity_mtld(doc)
    )
    features["lexical_density"].append(surface_proxies.lexical_density(doc))
    features["pos_dissimilarity"].append(
        surface_proxies.pos_dissimilarity(doc)
    )
    features["connection_words_ratio"].append(
        surface_proxies.connection_words_ratio(doc)
    )

    # In our case corpus was organized as:
    # corpus/5B/5_2_55.docx where the folder that
    # contained the doc, contained the school level, in
    # this example 5th grade
    features["grade"].append(filename.split("/")[1][0])

df = pd.DataFrame(features)


fig, axes = plt.subplots(2, 2)

sns.boxplot(x="grade", y="lexical_diversity_mltd", data=df, ax=axes[0, 0])
sns.boxplot(x="grade", y="lexical_density", data=df, ax=axes[0, 1])
sns.boxplot(x="grade", y="pos_dissimilarity", data=df, ax=axes[1, 0])
sns.boxplot(x="grade", y="connection_words_ratio", data=df, ax=axes[1, 1])
```

Which yields:

<img width="600" height="480" src="https://raw.githubusercontent.com/dpalmasan/TRUNAJOD2.0/master/imgs/figure2.png">

### _TRUNAJOD_ web app example

`TRUNAJOD` web app backend was built using `TRUNAJOD` lib. A demo video is shown below (it is in Spanish):

[![TRUNAJOD demo](https://img.youtube.com/vi/wl3ImqEVjeQ/0.jpg)](https://www.youtube.com/watch?v=wl3ImqEVjeQ)

## Contributing to _TRUNAJOD_

Bug reports and fixes are always welcome! Feel free to file issues, or ask for a feature request. We use `Github` issue tracker for this. If you'd like to contribute, feel free to submit a pull request. For more questions you can contact me at `dipalma (at) udec (dot) cl`.

More details can be found in
[CONTRIBUTING](https://github.com/dpalmasan/TRUNAJOD2.0/blob/master/CONTRIBUTING.md).


## References

If you find anything of this useful, feel free to cite the following papers, from which a lot of this python library was made for (I am also in the process of submitting this lib to an open software journal):

1. [Palma, D., & Atkinson, J. (2018). Coherence-based automatic essay assessment. IEEE Intelligent Systems, 33(5), 26-36.](https://ieeexplore.ieee.org/abstract/document/8506398/)
2. [Palma, D., Soto, C., Veliz, M., Riffo, B., & Gutiérrez, A. (2019, August). A Data-Driven Methodology to Assess Text Complexity Based on Syntactic and Semantic Measurements. In International Conference on Human Interaction and Emerging Technologies (pp. 509-515). Springer, Cham.](https://link.springer.com/chapter/10.1007/978-3-030-25629-6_79)

```bib
@article{Palma2021,
  doi = {10.21105/joss.03153},
  url = {https://doi.org/10.21105/joss.03153},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {60},
  pages = {3153},
  author = {Diego A. Palma and Christian Soto and Mónica Veliz and Bruno Karelovic and Bernardo Riffo},
  title = {TRUNAJOD: A text complexity library to enhance natural language processing},
  journal = {Journal of Open Source Software}
}

@article{palma2018coherence,
  title={Coherence-based automatic essay assessment},
  author={Palma, Diego and Atkinson, John},
  journal={IEEE Intelligent Systems},
  volume={33},
  number={5},
  pages={26--36},
  year={2018},
  publisher={IEEE}
}

@inproceedings{palma2019data,
  title={A Data-Driven Methodology to Assess Text Complexity Based on Syntactic and Semantic Measurements},
  author={Palma, Diego and Soto, Christian and Veliz, M{\'o}nica and Riffo, Bernardo and Guti{\'e}rrez, Antonio},
  booktitle={International Conference on Human Interaction and Emerging Technologies},
  pages={509--515},
  year={2019},
  organization={Springer}
}
```
# Contributing to TRUNAJOD

:+1::tada: Thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to `TRUNAJOD`. These are mostly guidelines, however there are some special steps you need to be aware of when creating a pull request so we end up with a standardized code. Please feel free to propose changes to this document in a PR.

#### Table Of Contents

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  * [Git Commit Messages](#git-commit-messages)
  * [Python Styleguide](#python-styleguide)
  * [Documentation Styleguide](#documentation-styleguide)

## How Can I Contribute?

### Reporting Bugs

If you encounter a bug or something in the functionality does not seem right, feel free to file an issue in the Github issue tracker. Ideally, provide an isolated example on how to reproduce the issue. The minimal information to help us assess the bug/issue:

* Clear description of the bug
* Current behavior and expected behavior
* Ideally a minimal example on how to reproduce it. We know that it is not always straightforward doing this, so this is optional.

### Suggesting Enhancements

Feel free to suggest any enhancement. Note that if you'd like a certain new index to be implemented, please provide in the description a reference paper or article, so we understand better what the enhancement is about and how can we make sure it is correctly implemented. 

We accept any kind of enhancements, so if you feel that usability, API, documentation can be improved, feel free to submit an issue.

### Your First Code Contribution

Unsure where to begin contributing to `TRUNAJOD`? You can start by looking through these `good-first-issue` and `help-wanted` issues:

* [Good first issue][good-first-issue] - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues][help-wanted] - issues which should be a bit more involved than `good first issue` issues.

### Pull Requests

Feel free to submit a pull request if you have a proposal for a feature or bugfix. However, make sure that the changes description are clear and you provide unit tests for the new code. Moreover, make sure that all the worflows in the CI (Github actions) pass, as we'd prefer not having regressions on the test cases. Finally, make sure to add citations if you add an algorithm. This will involve updating the `.bib` file for the updated module and making sure that you add the citation in the docstring of the new function. For example:

```python
def d_estimate(
    doc: Doc, min_range: int = 35, max_range: int = 50, trials: int = 5
) -> float:
    r"""Compute D measurement for lexical diversity.

    The measurement is based in :cite:`richards2000measuring`. We pick ``n``
    numbers of tokens, varying ``N`` from ``min_range`` up to ``max_range``.
    For each ``n`` we do the following:

    # Continues...
    """
```

## Styleguides

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line

### Python Styleguide

We recommend using `pre-commit` hooks, which are already set up in `.pre-commit-config.yaml`. In particular, we are using the following checks:

* [black](https://github.com/psf/black)
* [flake8](https://pypi.org/project/flake8/)
* [reorder_python_imports](https://github.com/asottile/reorder_python_imports)
* [pydocstyle](https://pypi.org/project/pydocstyle/)

Therefore, before committing, we recommend installing [pre-commit](https://pre-commit.com/). To have everything set up, please follow these steps:

* `pip install pre-commit`
* `pre-commit install`

Now, before committing, `pre-commit` hooks will be run, and you can ensure the checks are passed before the commit is created. Finally, make sure that you use typehints for the new functions you define, as we are planning to add `mypy` checks to the workflows.

### Documentation Styleguide

We use restructured text format, this is important as the documentation is hosted in [Read the Docs](https://readthedocs.org/). In particular, `TRUNAJOD` docs can be found in [https://trunajod20.readthedocs.io/en/latest/](https://trunajod20.readthedocs.io/en/latest/).

You can test the docs locally by running `sphinx-build -a -b html -W docs/ docs/_build/`. You might need to run `pip install -r docs/requirements.txt` in order to properly generate the documentation. Finally, you can check the docs by opening `docs/_build/index.html` file in your favorite browser.

[good-first-issue]:https://github.com/dpalmasan/TRUNAJOD2.0/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22
[help-wanted]:https://github.com/dpalmasan/TRUNAJOD2.0/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
.. mdinclude:: ../README.md

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api_reference/root
.. _ref-api-reference-givenness:

Givenness
=========

.. automodule:: TRUNAJOD.givenness
    :members:

.. bibliography:: givenness_ref.bib.. _ref-api-reference-entity_grid:

Entity Grids
============

.. automodule:: TRUNAJOD.entity_grid
    :members:

.. bibliography:: entity_grid.bib.. _ref-api-reference-lexico_semantic_norms:

Lexico-Semantic Norms
=====================

.. automodule:: TRUNAJOD.lexico_semantic_norms
    :members:

.. bibliography:: lexico_ref.bib
.. _ref-api-reference-discourse_markers:

Discourse Markers
=================

.. automodule:: TRUNAJOD.discourse_markers
    :members:

.. bibliography:: discourse_markers_ref.bib.. _ref-api-reference-utils:

Utils
=====

.. automodule:: TRUNAJOD.utils
    :members:.. _api-reference:

*************
API Reference
*************

.. toctree::
   :maxdepth: 2

   discourse_markers
   emotions
   entity_grid
   givenness
   lexico_semantic_norms
   semantic_measures
   surface_proxies
   syllabizer
   ttr
   utils
   .. _ref-api-reference-emotions:

Emotions
========

.. automodule:: TRUNAJOD.emotions
    :members:

.. bibliography:: emotions_ref.bib.. _ref-api-reference-ttr:

Type Token Ratios
=================

.. automodule:: TRUNAJOD.ttr
    :members:

.. bibliography:: ttr.bib.. _ref-api-reference-surface_proxies:

Surface Proxies
================

.. automodule:: TRUNAJOD.surface_proxies
    :members:.. _ref-api-reference-syllabizer:

Syllabizer
==========

.. automodule:: TRUNAJOD.syllabizer
    :members:.. _ref-api-reference-semantic_measurements:

Semantic Measures
=================

.. automodule:: TRUNAJOD.semantic_measures
    :members:


.. bibliography:: semantic_ref.bib
