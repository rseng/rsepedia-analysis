# libfmp

This repository contains the Python package libfmp. This package goes hand in hand with the FMP Notebooks, a collection of educational material for teaching and learning Fundamentals of Music Processing (FMP) with a particular focus on the audio domain. For detailed explanations and example appliciations of the libfmp-functions we refer to the FMP Notebooks:

https://audiolabs-erlangen.de/FMP

The FMP notebooks also contain a dedicated notebook for libfmp:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/B/B_libfmp.html

There is also an API documentation for libfmp:

https://meinardmueller.github.io/libfmp

If you use the package libfmp, please consider the following references.

## References

Meinard Müller and Frank Zalkow. [libfmp: A Python Package for Fundamentals of Music Processing.](https://joss.theoj.org/papers/10.21105/joss.03326) Journal of Open Source Software (JOSS), 6(63), 2021.

Meinard Müller and Frank Zalkow. [FMP Notebooks: Educational Material for Teaching and Learning Fundamentals of Music Processing.](https://archives.ismir.net/ismir2019/paper/000069.pdf) Proceedings of the International Conference on Music Information Retrieval (ISMIR), pp. 573&ndash;580, Delft, The Netherlands, 2019.

Meinard Müller. [Fundamentals of Music Processing &ndash; Using Python and Jupyter Notebooks.](http://www.music-processing.de/) Springer Verlag, 2nd edition, 2021.

Meinard Müller. [An Educational Guide Through the FMP Notebooks for Teaching and Learning Fundamentals of Music Processing.](https://www.mdpi.com/2624-6120/2/2/18) Signals, 2(2): 245&ndash;285, 2021.

## Statement of Need

The libfmp package bundles core concepts from the music information retrieval (MIR) field in the form of well-documented and easy-to-use Python functions. It is designed to aid students with the transition from being learners (e.g., studying the FMP notebooks) to becoming researchers by providing proper software support for building and experimenting with complex MIR pipelines. Going beyond and complementing existing Python packages (such as librosa), the libfmp package contains (previously unpublished) reference implementations of MIR algorithms from the literature and new Python implementations of previously published MATLAB toolboxes. The functionality of libfmp addresses diverse MIR tasks such as tuning estimation, music structure analysis, audio thumbnailing, chord recognition, tempo estimation, beat and local pulse tracking, fragment-level music retrieval, and audio decomposition.

## Installing

With Python >= 3.6, you can install libfmp using the Python package manager pip:

```
pip install libfmp
```

## Contributing

The libfmp-package has been developed in the context of the FMP notebooks. Being an integral part, all libfmp-functions need to manually synchronized with text passages, explanations, and the code in the FMP notebooks. Of course, we are happy for suggestions and contributions. However, to facilitate the synchronization, we would be grateful for either directly contacting us via email (meinard.mueller@audiolabs-erlangen.de) or for creating [an issue](https://github.com/meinardmueller/libfmp/issues) in our GitHub repository. Please do not submit a pull request without prior consultation with us.

If you want to report an issue with libfmp or seek support, please use the same communication channels (email or GitHub issue).

## Tests

The functions of libmfp are also covered in the [FMP notebooks](https://audiolabs-erlangen.de/FMP). There, you find several test cases for the functions, showing typical input-output behaviors. Beyond these tests, the FMP notebooks offer extensive explanations of these functions. Thus, we consider FMP as a replacement for conventional unit tests.

Furthermore, we provide a small script that tests one function of each subpackage from libfmp. Rather than covering the full functionality of libfmp, it only verifies the correct import structure within the libfmp package.

There are two options for executing the test script. The first is just to run the script, which results in no output if there are no errors.

```
python test_examples.py
```

The second option is to use [pytest](https://pytest.org), which results in a more instructive output. pytest is available when installing libfmp with the extra requirements for testing.

```
pip install 'libfmp[tests]'
pytest test_examples.py
```

## Acknowledgements

The main authors of libfmp, Meinard Müller and Frank Zalkow, are associated with the International Audio Laboratories Erlangen, which are a joint institution of the Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU) and Fraunhofer Institute for Integrated Circuits IIS. We thank the German Research Foundation (DFG) for various research grants that allow us for conducting fundamental research in music processing. Furthermore, we thank the various people who have contributed to libfmp with code and suggestions. In particular, we want to thank (in alphabetic order) Stefan Balke, Michael Krause, Patricio Lopez-Serrano, Julian Reck, Sebastian Rosenzweig, Angel Villar-Corrales, Christof Weiß, and Tim Zunner.
# Contributing to libfmp

The libfmp-package goes hand in hand with the FMP notebooks. In particular, we need to manually synchronize all libfmp-functions with text passages, explanations, and code in the FMP notebooks. Of course, we are happy for suggestions and contributions. However, to facilitate the synchronization, we would be grateful for either directly contacting us via email (meinard.mueller@audiolabs-erlangen.de) or for creating [an issue](https://github.com/meinardmueller/libfmp/issues) in our Github repository. Please do not submit a pull request without prior consultation with us.

By contributing, you agree that your contributions will be licensed under the MIT License.
Tempo and Beat Tracking (libfmp.c6)
===================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C6/C6.html

.. automodule:: libfmp.c6
    :members:
    :undoc-members:
.. automodule:: libfmp.c6.c6s1_onset_detection
    :members:
    :undoc-members:
.. automodule:: libfmp.c6.c6s1_peak_picking
    :members:
    :undoc-members:
.. automodule:: libfmp.c6.c6s2_tempo_analysis
    :members:
    :undoc-members:
.. automodule:: libfmp.c6.c6s3_adaptive_windowing
    :members:
    :undoc-members:
.. automodule:: libfmp.c6.c6s3_beat_tracking
    :members:
    :undoc-members:
Music Structure Analysis (libfmp.c4)
====================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C4/C4.html

.. automodule:: libfmp.c4
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s1_annotation
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s2_ssm
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s2_synthetic_ssm
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s2_threshold
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s3_thumbnail
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s4_novelty_kernel
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s4_structure_feature
    :members:
    :undoc-members:
.. automodule:: libfmp.c4.c4s5_evaluation
    :members:
    :undoc-members:
Fourier Analysis of Signals (libfmp.c2)
=======================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C2/C2.html

.. automodule:: libfmp.c2
    :members:
    :undoc-members:
.. automodule:: libfmp.c2.c2_complex
    :members:
    :undoc-members:
.. automodule:: libfmp.c2.c2_digitization
    :members:
    :undoc-members:
.. automodule:: libfmp.c2.c2_fourier
    :members:
    :undoc-members:
.. automodule:: libfmp.c2.c2_interference
    :members:
    :undoc-members:
.. automodule:: libfmp.c2.c2_interpolation
    :members:
    :undoc-members:
Chord Recognition (libfmp.c5)
=============================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C5/C5.html

.. automodule:: libfmp.c5
    :members:
    :undoc-members:
.. automodule:: libfmp.c5.c5s1_basic_theory_harmony
    :members:
    :undoc-members:
.. automodule:: libfmp.c5.c5s2_chord_rec_template
    :members:
    :undoc-members:
.. automodule:: libfmp.c5.c5s3_chord_rec_hmm
    :members:
    :undoc-members:
Content-Based Audio Retrieval (libfmp.c7)
=========================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C7/C7.html

.. automodule:: libfmp.c7
    :members:
    :undoc-members:
.. automodule:: libfmp.c7.c7s1_audio_id
    :members:
    :undoc-members:
.. automodule:: libfmp.c7.c7s2_audio_matching
    :members:
    :undoc-members:
.. automodule:: libfmp.c7.c7s3_version_id
    :members:
    :undoc-members:
Index
=====
Musically Informed Audio Decomposition (libfmp.c8)
==================================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C8/C8.html

.. automodule:: libfmp.c8
    :members:
    :undoc-members:
.. automodule:: libfmp.c8.c8s1_hps
    :members:
    :undoc-members:
.. automodule:: libfmp.c8.c8s2_f0
    :members:
    :undoc-members:
.. automodule:: libfmp.c8.c8s2_salience
    :members:
    :undoc-members:
.. automodule:: libfmp.c8.c8s3_nmf
    :members:
    :undoc-members:
Basics (libfmp.b)
=================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/B/B.html

.. automodule:: libfmp.b
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_annotation
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_audio
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_layout
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_plot
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_sonification
    :members:
    :undoc-members:
.. automodule:: libfmp.b.b_test_module
    :members:
    :undoc-members:
Music Synchronization (libfmp.c3)
=================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C3/C3.html

.. automodule:: libfmp.c3
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s1_audio_feature
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s1_post_processing
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s1_transposition_tuning
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s2_dtw
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s2_dtw_plot
    :members:
    :undoc-members:
.. automodule:: libfmp.c3.c3s3_tempo_curve
    :members:
    :undoc-members:
Getting Started
===============

You can install libfmp using the Python package manager pip:

.. code-block:: bash

    pip install libfmp

Beyond the API documentation of this webpage, you find extensive explanations of libfmp's functionality in the FMP Notebooks:

https://www.audiolabs-erlangen.de/FMP

In particular, there are dedicated notebooks on how to get started with FMP and on libfmp.

https://www.audiolabs-erlangen.de/resources/MIR/FMP/B/B_GetStarted.html
https://www.audiolabs-erlangen.de/resources/MIR/FMP/B/B_libfmp.html
Libfmp API Documentation
========================

This webpage contains the API documentation for the Python package libfmp.
This package goes hand in hand with the FMP Notebooks, a collection of educational material for teaching and learning Fundamentals of Music Processing (FMP) with a particular focus on the audio domain.
For detailed explanations and example applications of the libfmp-functions, we refer to the FMP Notebooks:

http://audiolabs-erlangen.de/FMP

The source code for the package libfmp is hosted at GitHub:

https://github.com/meinardmueller/libfmp

If you use the package libfmp, please consider the following references.

.. [#] Meinard Müller and Frank Zalkow. `libfmp: A Python Package for Fundamentals of Music Processing. <https://joss.theoj.org/papers/10.21105/joss.03326>`_ Journal of Open Source Software (JOSS), 6(63), 2021.

.. [#] Meinard Müller and Frank Zalkow. `FMP Notebooks: Educational Material for Teaching and Learning Fundamentals of Music Processing. <https://archives.ismir.net/ismir2019/paper/000069.pdf>`_ Proceedings of the International Conference on Music Information Retrieval (ISMIR), pp. 573–580, Delft, The Netherlands, 2019.

.. [#] Meinard Müller. `Fundamentals of Music Processing – Using Python and Jupyter Notebooks. <http://www.music-processing.de/>`_ Springer Verlag, 2nd edition, 2021.

.. [#] Meinard Müller. `An Educational Guide Through the FMP Notebooks for Teaching and Learning Fundamentals of Music Processing. <https://www.mdpi.com/2624-6120/2/2/18>`_ Signals, 2(2): 245–285, 2021.

.. toctree::
    :hidden:

    getting_started


.. toctree::
    :caption: API Documentation
    :maxdepth: 1
    :hidden:

    index_b
    index_c1
    index_c2
    index_c3
    index_c4
    index_c5
    index_c6
    index_c7
    index_c8

.. toctree::
    :caption: Reference
    :maxdepth: 1
    :hidden:

    genindex
    py-modindex
Module Index
============
Music Representations (libfmp.c1)
=================================

The `FMP notebooks <https://www.audiolabs-erlangen.de/FMP>`_ provide detailed textbook-like explanations of central techniques and algorithms implemented in the libfmp.
The part of FMP related to this module is available at the following URL:

https://www.audiolabs-erlangen.de/resources/MIR/FMP/C1/C1.html

.. automodule:: libfmp.c1
    :members:
    :undoc-members:
.. automodule:: libfmp.c1.c1s1_sheet_music
    :members:
    :undoc-members:
.. automodule:: libfmp.c1.c1s2_symbolic_rep
    :members:
    :undoc-members:
.. automodule:: libfmp.c1.c1s3_audio_rep
    :members:
    :undoc-members:
