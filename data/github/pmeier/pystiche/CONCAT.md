---
title: 'pystiche: A Framework for Neural Style Transfer'
tags:
  - Python
  - Neural Style Transfer
  - framework
  - PyTorch
authors:
  - name: Philip Meier
    orcid: 0000-0002-5184-1622
    affiliation: 1
  - name: Volker Lohweg
    orcid: 0000-0002-3325-7887
    affiliation: 1
affiliations:
 - name: >
     inIT –- Institute Industrial IT, 
     Technische Hochschule Ostwestfalen-Lippe (TH-OWL)
   index: 1
   
date: 15 October 2020
bibliography: docs/source/references.bib
---

# Summary

The seminal work of Gatys, Ecker, and Bethge gave birth to the field of 
_Neural Style Transfer_ (NST) in 2016 [@GEB2016]. The general idea behind an NST can be 
conveyed with only three images and two symbols:

![](docs/source/graphics/banner/banner.jpg)

In words: an NST describes the merger between the content and artistic style of two 
arbitrary images. This idea is nothing new in the field of 
_Non-photorealistic rendering_ (NPR) [@GG2001]. What distinguishes NST from traditional 
NPR approaches is its generality: an NST only needs a single arbitrary content and 
style image as input and thus "makes -- for the first time -- a generalized style 
transfer practicable" [@SID2017].

Besides peripheral tasks, an NST at its core is the definition of an optimization 
criterion called _perceptual loss_, which estimates the perceptual quality of the 
stylized image. Usually the perceptual loss comprises a deep neural network that needs 
to supply encodings of images from various depths.

`pystiche` is a framework for NST written in Python and built upon the _Deep Learning_ 
(DL) framework PyTorch [@PGM+2019]. It provides modular and efficient implementations 
for commonly used perceptual losses [@MV2015; @GEB2016; @LW2016] as well as neural net 
architectures [@SZ2014; @KSH2012]. This enables users to mix current state-of-the-art 
techniques with new ideas with ease. 

Due to its vivid nature, the field of NST gained a lot of traction in the short time 
after its emergence [@JYF+2019]. While many new techniques or augmentations have been 
developed, the field lacks standardization, which is especially evident in the 
reference implementations of the authors. `pystiche` aims to fill this gap.

# Statement of Need

Currently, unlike DL, there exist no library or framework for implementing NST. Thus, 
authors of new NST techniques either implement everything from scratch or base their 
implementation upon existing ones of other authors. Both ways have their downsides: 
while the former dampens innovation due to the lengthy implementation of reusable 
parts, with the latter the author inherits the technical debt due to the rapid 
development pace of DL hard- and software. In order to overcome this, `pystiche` 
pursues similar goals as DL frameworks:

1. **Accessibility**
   Starting off with NST can be quite overwhelming due to the sheer amount of 
   techniques one has to know and be able to deploy. `pystiche` aims to provide an 
   easy-to-use interface that reduces the necessary prior knowledge about NST and DL 
   to a minimum.
2. **Reproducibility**
   Implementing NST from scratch is not only inconvenient but also error-prone. 
   `pystiche` aims to provide reusable tools that let developers focus on their ideas 
   rather than worrying about everything around it.

`pystiche`s core audience are researchers, but its easy-to-use user interface 
opens up the field of NST for recreational use by laypersons.

# Acknowledgements

This contribution is part of the project _Fused Security Features_, which is funded by 
the _Ministry for Culture and Science of North Rhine-Westphalia_ (MKW NRW) under the 
Grant ID `005-1703-0013`. The authors thank Julian Bültemeier for extensive internal 
testing.

# References
# Scripts

This folder contains some scripts

## `generate_requirements_rtd.py`

The documentation, built by [Read the Docs (RTD)](https://readthedocs.org/), first 
installs `docs/requirements-rtd.txt` before installing `pystiche`. This has two 
reasons:

1. Installing PyTorch distributionswith CUDA support exceeds their memory limit. Thus, 
   we need to make sure to install it with CPU support only. 
2. The additional dependencies to build the documentation only live in `tox.ini`. Thus, 
   we need to extract them.

This script automatically populates `docs/requirements-rtd.txt`. It requires

- `pyyaml` and
- `light-the-torch>=0.2`

to be installed.

## `perform_model_optimization`

Trainings script for the model-optimization example. Requires the root directory of the 
dataset used for training as positional argument.
Contributing
============

First and foremost: Thank you for your interest in ``pystiche`` s development! We
appreciate all contributions be it code or something else.

Guide lines
-----------

``pystiche`` uses the `GitHub workflow <https://guides.github.com/introduction/flow/>`_
. Below is small guide how to make your first contribution.

.. note::

  The following guide assumes that `git <https://git-scm.com/>`_,
  `python <https://www.python.org/>`_, and `pip <https://pypi.org/project/pip/>`_ ,
  are available on your system. If that is not the case, follow the official
  installation instructions.

.. note::

  ``pystiche`` officially supports Python ``3.6``, ``3.7``, and ``3.8``. To ensure
  backwards compatibility, the development should happen on the minimum Python
  version, i. e. ``3.6``.

1. Fork ``pystiche`` on GitHub

  Navigate to `pmeier/pystiche <https://github.com/pmeier/pystiche>`_ on GitHub and
  click the **Fork** button in the top right corner.

2. Clone your fork to your local file system

  Use ``git clone`` to get a local copy of ``pystiche`` s repository that you can work
  on:

  .. code-block:: sh

    $ PYSTICHE_ROOT="pystiche"
    $ git clone "https://github.com/pmeier/pystiche.git" $PYSTICHE_ROOT

3. Setup your development environment

  .. code-block:: sh

    $ cd $PYSTICHE_ROOT
    $ virtualenv .venv --prompt="(pystiche) "
    $ source .venv/bin/activate
    $ pip install -r requirements-dev.txt
    $ pre-commit install

  .. note::

    While ``pystiche`` s development requirements are fairly lightweight, it is still
    recommended to install them in a virtual environment rather than system wide. If you
    do not have ``virtualenv`` installed, you can do so by running
    ``pip install --user virtualenv``.

4. Create a branch for local development

  Use ``git checkout`` to create local branch with a descriptive name:

  .. code-block:: sh

    $ PYSTICHE_BRANCH="my-awesome-feature-or-bug-fix"
    $ git checkout -b $PYSTICHE_BRANCH

  Now make your changes. Happy Coding!

5. Use ``tox`` to run various checks

  .. code-block:: sh

    $ tox

  .. note::

    Running ``tox`` is equivalent to running

    .. code-block:: sh

      $ tox -e lint-style
      $ tox -e lint-typing
      $ tox -e tests-integration
      $ tox -e tests-galleries
      $ tox -e tests-docs

    You can find details what the individual commands do below of this guide.

6. Commit and push your changes

  If all checks are passing you can commit your changes an push them to your fork:

  .. code-block:: sh

    $ git add .
    $ git commit -m "Descriptive message of the changes made"
    $ git push -u origin $PYSTICHE_BRANCH

  .. note::

    For larger changes, it is good practice to split them in multiple small commits
    rather than one large one. If you do that, make sure to run the test suite before
    every commit. Furthermore, use ``git push`` without any parameters for consecutive
    commits.

7. Open a Pull request (PR)

  1. Navigate to `pmeier/pystiche/pulls <https://github.com/pmeier/pystiche/pulls>`_ on
     GitHub and click on the green button "New pull request".
  2. Click on "compare across forks" below the "Compare changes" headline.
  3. Select your fork for "head repository" and your branch for "compare" in the
     drop-down menus.
  4. Click the the green button "Create pull request".

  .. note::

    If the time between the branch being pushed and the PR being opened is not too
    long, GitHub will offer you a yellow box after step 1. If you click the button,
    you can skip steps 2. and 3.

.. note::

  Steps 1. to 3. only have to performed once. If you want to continue contributing,
  make sure to branch from the current ``master`` branch. You can use ``git pull``

  .. code-block:: sh

    $ git checkout master
    $ git pull origin
    $ git checkout -b "my-second-awesome-feature-or-bug-fix"

  If you forgot to do that or if since the creation of your branch many commits have
  been made to the ``master`` branch, simply rebase your branch on top of it.

  .. code-block:: sh

    $ git checkout master
    $ git pull origin
    $ git checkout "my-second-awesome-feature-or-bug-fix"
    $ git rebase master

If you are contributing bug-fixes or
documentation improvements, you can open a
`pull request (PR) <https://github.com/pmeier/pystiche/pulls>`_ without further
discussion. If on the other hand you are planning to contribute new features, please
open an `issue <https://github.com/pmeier/pystiche/issues>`_ and discuss the feature
with us first.

Every PR is subjected to multiple automatic checks (continuous integration, CI) as well
as a manual code review that it has to pass before it can be merged. The automatic
checks are performed by `tox <https://tox.readthedocs.io/en/latest/>`_. You can find
details and instructions how to run the checks locally below.

Code format and linting
-----------------------

``pystiche`` uses `isort <https://timothycrosley.github.io/isort/>`_ to sort the
imports, `black <https://black.readthedocs.io/en/stable/>`_ to format the code, and
`flake8 <https://flake8.pycqa.org/en/latest/>`_ to enforce
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ compliance. To format and check the
code style, run

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -e lint-style

.. note::

  Amongst others, ``isort``, ``black``, and ``flake8`` are run by
  `pre-commit <https://pre-commit.com/>`_ before every commit.

Furthermore, ``pystiche_papers`` is
`PEP561 <https://www.python.org/dev/peps/pep-0561/>`_ compliant and checks the type
annotations with `mypy <http://mypy-lang.org/>`_. To check the static typing, run

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -e lint-typing

For convenience, you can run all lint checks with

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -f lint


Test suite
----------

``pystiche`` uses `pytest <https://docs.pytest.org/en/stable/>`_ to run the test suite.
You can run it locally with

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox

.. note::

  ``pystiche_papers`` adds the following custom options with the
  corresponding ``@pytest.mark.*`` decorators:
  - ``--skip-large-download``: ``@pytest.mark.large_download``
  - ``--skip-slow``: ``@pytest.mark.slow``
  - ``--run-flaky``: ``@pytest.mark.flaky``

  Options prefixed with ``--skip`` are run by default and skipped if the option is
  given. Options prefixed with ``--run`` are skipped by default and run if the option
  is given.

  These options are passed through ``tox`` if given after a ``--`` flag. For example,
  the CI invokes the test suite with

  .. code-block:: sh

    cd $PYSTICHE_ROOT
    source .venv/bin/activate
    tox -- --skip-large-download


Documentation
-------------

To build the html documentation locally, run

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -e docs-html

To build the latex (PDF) documentation locally, run

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -e docs-latex

To build both, run

.. code-block:: sh

  cd $PYSTICHE_ROOT
  source .venv/bin/activate
  tox -f docs

.. note::

  Building the documentation triggers a
  `sphinx gallery <https://sphinx-gallery.github.io/stable/index.html>`_ build by
  default for the example galleries. This which will take some time to complete. To get 
  around this, ``pystiche`` offers two environment variables:

  - ``PYSTICHE_PLOT_GALLERY``: If ``False``, the code inside the galleries is not
    executed. See the
    `official sphinx-gallery documentation <https://sphinx-gallery.github.io/stable/configuration.html#without-execution>`_
    for details. Defaults to ``True``.
  - ``PYSTICHE_DOWNLOAD_GALLERY``: If ``True``, downloads pre-built
    galleries and uses them instead of rebuilding. For the ``master`` the galleries are
    at most six hours old. Defaults to ``False``.

  Both environment variables are evaluated with :func:`~distutils.util.strtobool`.
.. start-badges

.. list-table::
    :stub-columns: 1

    * - package
      - |license| |status|
    * - citation
      - |pyopensci| |joss|
    * - code
      - |black| |mypy| |lint|
    * - tests
      - |tests| |coverage|
    * - docs
      - |docs| |rtd|

.. end-badges

|logo|

``pystiche``
============

``pystiche`` (pronounced
`/ˈpaɪˈstiʃ/ <http://ipa-reader.xyz/?text=%CB%88pa%C9%AA%CB%88sti%CA%83>`_ ) is a
framework for
`Neural Style Transfer (NST) <https://github.com/ycjing/Neural-Style-Transfer-Papers>`_
built upon `PyTorch <https://pytorch.org>`_. The name of the project is a pun on
*pastiche* `meaning <https://en.wikipedia.org/wiki/Pastiche>`_:

    A pastiche is a work of visual art [...] that imitates the style or character of
    the work of one or more other artists. Unlike parody, pastiche celebrates, rather
    than mocks, the work it imitates.

.. image:: docs/source/graphics/banner/banner.jpg
    :alt: pystiche banner

``pystiche`` has similar goals as Deep Learning (DL) frameworks such as PyTorch:

1. **Accessibility**
    Starting off with NST can be quite overwhelming due to the sheer amount of
    techniques one has to know and be able to deploy. ``pystiche`` aims to provide an
    easy-to-use interface that reduces the necessary prior knowledge about NST and DL
    to a minimum.
2. **Reproducibility**
    Implementing NST from scratch is not only inconvenient but also error-prone.
    ``pystiche`` aims to provide reusable tools that let developers focus on their
    ideas rather than worrying about bugs in everything around it.


Installation
============

``pystiche`` is a proper Python package and can be installed with ``pip``. The latest
release can be installed with

.. code-block:: sh

  pip install pystiche

Usage
=====

``pystiche`` makes it easy to define the optimization criterion for an NST task fully
compatible with PyTorch. For example, the banner above was generated with the following
``criterion``:

.. code-block:: python

  from pystiche import enc, loss

  mle = enc.vgg19_multi_layer_encoder()

  perceptual_loss = loss.PerceptualLoss(
      content_loss=loss.FeatureReconstructionLoss(
          mle.extract_encoder("relu4_2")
      ),
      style_loss=loss.MultiLayerEncodingLoss(
          mle,
          ("relu1_1", "relu2_1", "relu3_1", "relu4_1", "relu5_1"),
          lambda encoder, layer_weight: ops.GramOLoss(
              encoder, score_weight=layer_weight
          ),
          score_weight=1e3,
      ),
  )

For the full example, head over to the example
`NST with pystiche <https://pystiche.readthedocs.io/en/latest/galleries/examples/beginner/example_nst_with_pystiche.html>`_.

Documentation
=============

For

- `detailed installation instructions <https://pystiche.readthedocs.io/en/latest/getting_started/installation.html>`_,
- a `gallery of usage examples <https://pystiche.readthedocs.io/en/latest/galleries/examples/index.html>`_,
- the `API reference <https://pystiche.readthedocs.io/en/latest/api/index.html>`_,
- the `contributing guidelines <https://pystiche.readthedocs.io/en/latest/getting_started/contributing.html>`_,

or anything else, head over to the `documentation <https://pystiche.readthedocs.io/en/latest/>`_.

Citation
========

If you use this software, please cite it as

.. code-block:: bibtex

  @Article{ML2020,
    author  = {Meier, Philip and Lohweg, Volker},
    journal = {Journal of Open Source Software {JOSS}},
    title   = {pystiche: A Framework for Neural Style Transfer},
    year    = {2020},
    doi     = {10.21105/joss.02761},
  }

.. |logo|
  image:: logo.svg
    :target: https://pystiche.org
    :alt: pystiche logo

.. |license|
  image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://opensource.org/licenses/BSD-3-Clause
    :alt: License

.. |status|
  image:: https://www.repostatus.org/badges/latest/active.svg
    :target: https://www.repostatus.org/#active
    :alt: Project Status: Active

.. |pyopensci|
  image:: https://tinyurl.com/y22nb8up
    :target: https://github.com/pyOpenSci/software-review/issues/25
    :alt: pyOpenSci

.. |joss|
  image:: https://joss.theoj.org/papers/10.21105/joss.02761/status.svg
    :target: https://doi.org/10.21105/joss.02761
    :alt: JOSS

.. |black|
  image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. |mypy|
  image:: http://www.mypy-lang.org/static/mypy_badge.svg
    :target: http://mypy-lang.org/
    :alt: mypy

.. |lint|
  image:: https://github.com/pmeier/pystiche/workflows/lint/badge.svg
    :target: https://github.com/pmeier/pystiche/actions?query=workflow%3Alint+branch%3Amaster
    :alt: Lint status via GitHub Actions

.. |tests|
  image:: https://github.com/pmeier/pystiche/workflows/tests/badge.svg
    :target: https://github.com/pmeier/pystiche/actions?query=workflow%3Atests+branch%3Amaster
    :alt: Test status via GitHub Actions

.. |coverage|
  image:: https://codecov.io/gh/pmeier/pystiche/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/pmeier/pystiche
    :alt: Test coverage

.. |docs|
  image:: https://github.com/pmeier/pystiche/workflows/docs/badge.svg
    :target: https://github.com/pmeier/pystiche/actions?query=workflow%3Adocs+branch%3Amaster
    :alt: Docs status via GitHub Actions

.. |rtd|
  image:: https://img.shields.io/readthedocs/pystiche?label=latest&logo=read%20the%20docs
    :target: https://pystiche.readthedocs.io/en/latest/?badge=latest
    :alt: Latest documentation hosted on Read the Docs
.. _usage_examples:

``pystiche`` usage examples
===========================

.. note::

  Although a GPU is not a requirement, it is strongly advised to run these examples
  with one. If you don't have access to a GPU, the execution time of the examples might
  increase by multiple orders of magnitude. The total running time provided at the end
  of each example is measured using a GPU.
Advanced
--------
Beginner
--------
Literature Reference
====================

.. bibliography:: references.bib
.. image:: ../../logo.svg
    :alt: pystiche logo

Welcome to ``pystiche`` 's documentation!
=========================================

``pystiche`` (pronounced
`/ˈpaɪˈstiʃ/ <http://ipa-reader.xyz/?text=%CB%88pa%C9%AA%CB%88sti%CA%83>`_ ) is a
framework for
`Neural Style Transfer (NST) <https://github.com/ycjing/Neural-Style-Transfer-Papers>`_
built upon `PyTorch <https://pytorch.org>`_. The name of the project is a pun on
*pastiche* `meaning <https://en.wikipedia.org/wiki/Pastiche>`_:

    A pastiche is a work of visual art [...] that imitates the style or character of
    the work of one or more other artists. Unlike parody, pastiche celebrates, rather
    than mocks, the work it imitates.

``pystiche`` has similar goals as Deep Learning (DL) frameworks such as PyTorch:

1. **Accessibility**
    Starting off with NST can be quite overwhelming due to the sheer amount of
    techniques one has to know and be able to deploy. ``pystiche`` aims to provide an
    easy-to-use interface that reduces the necessary prior knowledge about NST and DL
    to a minimum.
2. **Reproducibility**
    Implementing NST from scratch is not only inconvenient but also error-prone.
    ``pystiche`` aims to provide reusable tools that let developers focus on their
    ideas rather than worrying about bugs in everything around it.

.. toctree::
  :maxdepth: 1
  :caption: Contents

  Getting Started <getting_started/index>
  Gist <gist/index>
  Usage Examples <galleries/examples/index>
  CLI <cli/index>
  Package Reference <api/index>
  Literature Reference <literature>
Perceptual loss
===============

The identification of content and style are core elements of a Neural Style Transfer
(NST). The agreement of the content and style of two images is measured with the
``content_loss`` and ``style_loss``, respectively.


Operators
---------

In ``pystiche`` these losses are implemented :class:`~pystiche.loss.Loss` s.
:class:`~pystiche.loss.Loss` s are differentiated between two  types:
:class:`~pystiche.loss.RegularizationLoss` and
:class:`~pystiche.loss.ComparisonLoss`. A
:class:`~pystiche.loss.RegularizationLoss` works without any context while a
:class:`~pystiche.loss.ComparisonLoss` compares two images. Furthermore,
``pystiche`` differentiates between two different domains an
:class:`~pystiche.loss.Loss` can work on: :class:`~pystiche.ops.op.PixelOperator`
and :class:`~pystiche.ops.EncodingOperator` . A :class:`~pystiche.ops.PixelOperator`
operates directly on the ``input_image`` while an
:class:`~pystiche.ops.EncodingOperator` encodes it first.

In total ``pystiche`` supports four archetypes:

+-------------------------------------------------------+-----------------------------------------------------------------------+
| :class:`~pystiche.loss.Loss`                          | Builtin examples                                                      |
+=======================================================+=======================================================================+
| :class:`~pystiche.ops.PixelRegularizationOperator`    | - :class:`~pystiche.loss.TotalVariationLoss` :cite:`MV2015`           |
+-------------------------------------------------------+-----------------------------------------------------------------------+
| :class:`~pystiche.ops.EncodingRegularizationOperator` |                                                                       |
+-------------------------------------------------------+-----------------------------------------------------------------------+
| :class:`~pystiche.ops.PixelComparisonOperator`        |                                                                       |
+-------------------------------------------------------+-----------------------------------------------------------------------+
| :class:`~pystiche.ops.EncodingComparisonOperator`     | - :class:`~pystiche.loss.FeatureReconstructionLoss` :cite:`MV2015`    |
|                                                       | - :class:`~pystiche.loss.GramLoss` :cite:`GEB2016`                    |
|                                                       | - :class:`~pystiche.loss.MRFLoss` :cite:`LW2016`                      |
+-------------------------------------------------------+-----------------------------------------------------------------------+

Multi-layer encoder
-------------------

One of the main improvements of NST compared to traditional approaches is that the
agreement is not measured in the pixel or a handcrafted feature space, but rather in
the learned feature space of a Convolutional Neural Network called ``encoder``.
Especially variants of the ``style_loss`` depend upon encodings, i. e. feature maps,
from various layers of the encoder.

``pystiche`` offers a
:class:`~pystiche.enc.MultiLayerEncoder` that enables to extract all required encodings
after a single forward pass. If the same operator should be applied to different layers
of a :class:`~pystiche.enc.MultiLayerEncoder`, a
:class:`~pystiche.loss.MultiLayerEncodingLoss` can be used.


Perceptual loss
---------------

The :class:`~pystiche.loss.PerceptualLoss` combines all :class:`~pystiche.ops.Operator`
s in a single measure acting as joint optimization criterion. How the optimization is
performed will be detailed in the next section.
Optimization
============

The merging of the identified content and style with a Neural Style Transfer (NST) is
posed as an optimization problem. The optimization is performed on the basis of a
:class:`~pystiche.loss.PerceptualLoss`. A distinction is made between two different
approaches.

Image optimization
------------------

In its basic form, an NST optimizes the pixels of the ``input_image`` directly. That
means they are iteratively adapted to reduce the perceptual loss. This
process is called *image optimization* and can be performed in ``pystiche`` with a
:func:`~pystiche.optim.image_optimization` .

Model optimization
------------------

While the image optimization approach yields the highest quality results, the
computation is quite expensive and usually takes multiple minutes to complete for a
single image. *Model optimization* on the other hand trains a model called
``transformer`` to perform the stylization. The training is performed with the same
perceptual loss as before, but now the ``transformer`` weights are used as optimization
parameters. The training is even more time consuming but afterwards the stylization is
performed in a single forward pass of the ``input_image`` through the ``transformer``.
The quality however, while still high, is lower than for image optimisation approaches
since the ``transformer`` cannot finetune the ``output_image``. In ``pystiche`` a model
optimization can be performed with a
:func:`~pystiche.optim.model_optimization` .

.. note::
  Due to the differences in execution time image and model optimization approaches are
  often dubbed *slow* and *fast* respectively.
.. _gist:

Gist
====

From a high viewpoint, Neural Style Transfer (NST) can be described with only three
images and two symbols:

.. image:: ../graphics/banner/banner.jpg
    :alt: pystiche banner

Not only the quality of the results but also the underlying steps are comparable to the
work of human artisans or craftsmen. Such a manual style transfer can be roughly
divided into three steps:

1. The content or motif of an image needs to be identified. That means one has to
   identify which parts of the image are essential and on the other hand which details
   can be discarded.
2. The style of an image, such as color, shapes, brush strokes, needs to be identified.
   Usually that means one has to intensively study of the works of the original artist.
3. The identified content and style have to be merged together. This can be the most
   difficult step, since it usually requires a lot of skill to match the style of
   another artist.

In principle an NST performs the same steps, albeit fully automatically. This is
nothing new in the field of computational style transfers. What makes NST stand out is
its generality: NST only needs a single arbitrary content and style image as input and
thus "makes -- for the first time -- a generalized style transfer practicable."
:cite:`SID2017`.

The following sections provide the gist of how these three steps are performed with
``pystiche`` as part of an NST . Afterwards head over to the
:ref:`usage examples <usage_examples>` to see ``pystiche`` in action.

.. toctree::
  :maxdepth: 2

  Perceptual loss <loss>
  Optimization <optim>
``pystiche``
============

.. automodule:: pystiche

.. autofunction:: home

Objects
-------

.. autoclass:: ComplexObject
  :members:
    _properties,
    extra_properties,
    properties,
    _named_children,
    extra_named_children,
    named_children,
.. autoclass:: LossDict
  :members:
    __setitem__,
    aggregate,
    total,
    backward,
    item,
    __mul__
.. autoclass:: Module
  :members: torch_repr

Math
----

.. autofunction:: nonnegsqrt
.. autofunction:: gram_matrix
.. autofunction:: cosine_similarity
``pystiche.meta``
=================

.. automodule:: pystiche.meta

.. autofunction:: tensor_meta
.. autofunction:: is_scalar_tensor
.. autofunction:: is_conv_module
.. autofunction:: conv_module_meta
.. autofunction:: is_pool_module
.. autofunction:: pool_module_meta
``pystiche.misc``
=================

.. automodule:: pystiche.misc

.. autofunction:: get_input_image
.. autofunction:: get_device
.. autofunction:: reduce
.. autofunction:: build_complex_obj_repr
``pystiche.loss``
=================

.. automodule:: pystiche.loss

.. autoclass:: Loss
    :members:
        set_input_guide,
        forward,
    :undoc-members:
    :show-inheritance:

.. autoclass:: RegularizationLoss
    :members:
        input_enc_to_repr,
        calculate_score,
    :undoc-members:
    :show-inheritance:

.. autoclass:: ComparisonLoss
    :members:
        set_target_image,
        input_enc_to_repr,
        target_enc_to_repr,
        calculate_score,
    :undoc-members:
    :show-inheritance:


Container
---------

.. autoclass:: LossContainer
    :members:
        set_input_guide,
        set_target_image,

.. autoclass:: MultiLayerEncodingLoss

.. autoclass:: MultiRegionLoss
    :members:
        set_regional_input_guide,
        set_regional_target_image,

.. autoclass:: PerceptualLoss
    :members:
        regional_content_guide,
        set_content_guide,
        regional_style_image,
        set_style_image,


Regularization
--------------

.. autoclass:: TotalVariationLoss
    :show-inheritance:


Comparison
----------

.. autoclass:: FeatureReconstructionLoss
    :show-inheritance:
.. autoclass:: GramLoss
    :show-inheritance:
.. autoclass:: MRFLoss
    :members: scale_and_rotate_transforms
    :show-inheritance:
``pystiche.pyramid``
====================

.. automodule:: pystiche.pyramid


Image pyramid
-------------

.. autoclass:: ImagePyramid
.. autoclass:: OctaveImagePyramid
  :show-inheritance:


Pyramid level
-------------

.. autoclass:: PyramidLevel
  :members: resize_image, resize_guide
``pystiche.enc``
================

.. automodule:: pystiche.enc

.. autoclass:: Encoder
  :members:
    forward,
    propagate_guide
.. autoclass:: SequentialEncoder

.. autoclass:: MultiLayerEncoder
  :members:
    __contains__,
    verify,
    register_layer,
    __call__,
    forward,
    clear_cache,
    encode,
    propagate_guide,
    trim,
    extract_encoder

.. autoclass:: SingleLayerEncoder
  :members:
    forward,
    propagate_guide

Models
------

.. autoclass:: ModelMultiLayerEncoder
  :members:
    state_dict_url,
    collect_modules,
    load_state_dict,
    load_state_dict_from_url

VGG
^^^

.. autoclass:: VGGMultiLayerEncoder
.. autofunction:: vgg11_multi_layer_encoder
.. autofunction:: vgg11_bn_multi_layer_encoder
.. autofunction:: vgg13_multi_layer_encoder
.. autofunction:: vgg13_bn_multi_layer_encoder
.. autofunction:: vgg16_multi_layer_encoder
.. autofunction:: vgg16_bn_multi_layer_encoder
.. autofunction:: vgg19_multi_layer_encoder
.. autofunction:: vgg19_bn_multi_layer_encoder

AlexNet
^^^^^^^

.. autoclass:: AlexNetMultiLayerEncoder
.. autofunction:: alexnet_multi_layer_encoder
Package reference
=================

.. toctree::
  :maxdepth: 1

  pystiche
  pystiche.data
  pystiche.demo
  pystiche.enc
  pystiche.image
  pystiche.loss
  pystiche.loss.functional
  pystiche.meta
  pystiche.misc
  pystiche.optim
  pystiche.pyramid
``pystiche.image``
==================

.. automodule:: pystiche.image


Utilities
---------

.. autofunction:: calculate_aspect_ratio
.. autofunction:: image_to_edge_size
.. autofunction:: edge_to_image_size
.. autofunction:: extract_batch_size
.. autofunction:: extract_num_channels
.. autofunction:: extract_image_size
.. autofunction:: extract_edge_size
.. autofunction:: extract_aspect_ratio


I/O
---

.. autofunction:: read_image
.. autofunction:: write_image
.. autofunction:: show_image


Guides
------

.. autofunction:: verify_guides
.. autofunction:: read_guides
.. autofunction:: write_guides
.. autofunction:: guides_to_segmentation
.. autofunction:: segmentation_to_guides
``pystiche.loss.functional``
============================

.. automodule:: pystiche.loss.functional
  :members:
  :undoc-members:
``pystiche.optim``
==================

.. automodule:: pystiche.optim

.. autofunction:: default_image_optimizer
.. autofunction:: image_optimization
.. autofunction:: pyramid_image_optimization

.. autofunction:: default_model_optimizer
.. autofunction:: model_optimization
.. autofunction:: multi_epoch_model_optimization
``pystiche.data``
=================

.. automodule:: pystiche.data

.. autoclass:: LocalImage
  :members: read
.. autoclass:: LocalImageCollection
  :members: read

.. autoclass:: DownloadableImage
  :members: generate_file, download, read
.. autoclass:: DownloadableImageCollection
  :members: download, read
.. include:: ../../../CONTRIBUTING.rst
.. _installation:

Installation
============

The latest **stable** version can be installed with

.. code-block:: sh

  pip install pystiche

The latest **potentially unstable** version can be installed with

.. code-block::

  pip install git+https://github.com/pmeier/pystiche@master


Installation of PyTorch
-----------------------

``pystiche`` is built upon `PyTorch <https://pytorch.org>`_ and depends on
``torch`` and ``torchvision``. By default, a ``pip install`` of ``pystiche`` tries to
install the PyTorch distributions precompiled for the latest CUDA release. If you use
another version or don't have a CUDA-capable GPU, we encourage you to try
`light-the-torch <https://github.com/pmeier/light-the-torch>`_ for a convenient
installation:

.. code-block:: sh

  pip install light-the-torch
  ltt install pystiche

Otherwise, please follow the
`official installation instructions of PyTorch <https://pytorch.org/get-started/>`_ for
your setup before you install ``pystiche``.

.. note::

  While ``pystiche`` is designed to be fully functional without a GPU, most tasks
  require significantly more time to perform on a CPU.
Getting started
===============

.. toctree::
  :maxdepth: 1

  Installation <installation>
  Contributing <contributing>
