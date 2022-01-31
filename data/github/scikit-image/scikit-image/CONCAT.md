# scikit-image: Image processing in Python
[![Image.sc forum](https://img.shields.io/badge/dynamic/json.svg?label=forum&url=https%3A%2F%2Fforum.image.sc%2Ftags%2Fscikit-image.json&query=%24.topic_list.tags.0.topic_count&colorB=brightgreen&suffix=%20topics&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAABPklEQVR42m3SyyqFURTA8Y2BER0TDyExZ+aSPIKUlPIITFzKeQWXwhBlQrmFgUzMMFLKZeguBu5y+//17dP3nc5vuPdee6299gohUYYaDGOyyACq4JmQVoFujOMR77hNfOAGM+hBOQqB9TjHD36xhAa04RCuuXeKOvwHVWIKL9jCK2bRiV284QgL8MwEjAneeo9VNOEaBhzALGtoRy02cIcWhE34jj5YxgW+E5Z4iTPkMYpPLCNY3hdOYEfNbKYdmNngZ1jyEzw7h7AIb3fRTQ95OAZ6yQpGYHMMtOTgouktYwxuXsHgWLLl+4x++Kx1FJrjLTagA77bTPvYgw1rRqY56e+w7GNYsqX6JfPwi7aR+Y5SA+BXtKIRfkfJAYgj14tpOF6+I46c4/cAM3UhM3JxyKsxiOIhH0IO6SH/A1Kb1WBeUjbkAAAAAElFTkSuQmCC)](https://forum.image.sc/tags/scikit-image)
[![Stackoverflow](https://img.shields.io/badge/stackoverflow-Ask%20questions-blue.svg)](https://stackoverflow.com/questions/tagged/scikit-image)
[![project chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://skimage.zulipchat.com)
[![codecov.io](https://codecov.io/github/scikit-image/scikit-image/coverage.svg?branch=main)](https://codecov.io/github/scikit-image/scikit-image?branch=main)

- **Website (including documentation):** [https://scikit-image.org/](https://scikit-image.org)
- **User forum:** [https://forum.image.sc/tag/scikit-image](https://forum.image.sc/tag/scikit-image)
- **Developer forum:** [https://discuss.scientific-python.org/c/contributor/skimage](https://discuss.scientific-python.org/c/contributor/skimage)
- **Source:** [https://github.com/scikit-image/scikit-image](https://github.com/scikit-image/scikit-image)
- **Benchmarks:** [https://pandas.pydata.org/speed/scikit-image/](https://pandas.pydata.org/speed/scikit-image/)

## Installation from binaries

- **Debian/Ubuntu:** ``sudo apt-get install python-skimage``
- **OSX:** ``pip install scikit-image``
- **Anaconda:** ``conda install -c conda-forge scikit-image``

Also see [installing ``scikit-image``](INSTALL.rst).

## Installation from source

Install dependencies using:

```
pip install -r requirements.txt
```

Then, install scikit-image using:

```
$ pip install .
```

If you plan to develop the package, you may run it directly from source:

```
$ pip install -e .  # Do this once to add package to Python path
```

Every time you modify Cython files, also run:

```
$ python setup.py build_ext -i  # Build binary extensions
```

## License (Modified BSD)

Copyright (C) 2011, the scikit-image team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of skimage nor the names of its contributors may be
    used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

## Citation

If you find this project useful, please cite:

> Stéfan van der Walt, Johannes L. Schönberger, Juan Nunez-Iglesias,
> François Boulogne, Joshua D. Warner, Neil Yager, Emmanuelle
> Gouillart, Tony Yu, and the scikit-image contributors.
> *scikit-image: Image processing in Python*. PeerJ 2:e453 (2014)
> https://doi.org/10.7717/peerj.453
[scikit-image Code of Conduct](doc/source/conduct/code_of_conduct.md)
# pip requirements files

## Index

- [default.txt](default.txt)
  Default requirements
- [docs.txt](docs.txt)
  Documentation requirements
- [optional.txt](optional.txt)
  Optional requirements. All of these are installable without a compiler through pypi.
- [extras.txt](extras.txt)
  Optional requirements that require a compiler to install.
- [test.txt](test.txt)
  Requirements for running test suite
- [build.txt](build.txt)
  Requirements for building from the source repository

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

## Justification for blocked versions

* Cython 0.28.2 was empircally found to fail tests while other patch-releases 0.28.x do not
* Cython 0.29.0 erroneously sets the `__path__` to `None`. See https://github.com/cython/cython/issues/2662
* Cython 0.29.18 fails due to a bad definition of M_PI. See https://github.com/cython/cython/issues/3622
* matplotlib 3.0.0 is not used because of a bug that collapses 3D axes (see https://github.com/scikit-image/scikit-image/pull/3474 and https://github.com/matplotlib/matplotlib/issues/12239).
* pillow 7.1.0 fails on png files, See https://github.com/scikit-image/scikit-image/issues/4548
* pillow 7.1.1 fails due to https://github.com/python-pillow/Pillow/issues/4518
* pillow 8.3.0 broke array coercion when dtype was specified https://github.com/python-pillow/Pillow/pull/5572
* imread 0.7.2 fails due to build failure https://github.com/luispedro/imread/issues/36
* sphinx-gallery 0.8.0 is banned due to bug introduced on binder: https://github.com/scikit-image/scikit-image/pull/4959#issuecomment-687653537
## Description

<!--
(Note: for guidance on how to use `scikit-image`, please post instead on https://forum.image.sc/tag/scikit-image)
-->

## Way to reproduce
```python
# Place the full code we need to recreate your issue here
# upload all necessary images to github too!
```


## Version information
```python
# Paste the output of the following python commands
from __future__ import print_function
import sys; print(sys.version)
import platform; print(platform.platform())
import skimage; print(f'scikit-image version: {skimage.__version__}')
import numpy; print(f'numpy version: {numpy.__version__}')
```

```python
# your output here

```
## Description

<!-- If this is a bug-fix or enhancement, state the issue # it closes -->
<!-- If this is a new feature, reference what paper it implements. -->


## Checklist

<!-- It's fine to submit PRs which are a work in progress! -->
<!-- But before they are merged, all PRs should provide: -->
- [Docstrings for all functions](https://github.com/numpy/numpy/blob/master/doc/example.py)
- Gallery example in `./doc/examples` (new features only)
- Benchmark in `./benchmarks`, if your changes aren't covered by an
  existing benchmark
- Unit tests
- Clean style in [the spirit of PEP8](https://www.python.org/dev/peps/pep-0008/)
- Descriptive commit messages (see below)

<!-- For detailed information on these and other aspects see -->
<!-- the scikit-image contribution guidelines. -->
<!-- https://scikit-image.org/docs/dev/contribute.html -->

## For reviewers

<!-- Don't remove the checklist below. -->
- Check that the PR title is short, concise, and will make sense 1 year
  later.
- Check that new functions are imported in corresponding `__init__.py`.
- Check that new features, API changes, and deprecations are mentioned in
  `doc/release/release_dev.rst`.
- There is a bot to help automate backporting a PR to an older branch. For
  example, to backport to v0.19.x after merging, add the following in a PR
  comment: `@meeseeksdev backport to v0.19.x`
- To run benchmarks on a PR, add the `run-benchmark` label. To rerun, the label
  can be removed and then added again. The benchmark output can be checked in
  the "Actions" tab.
# Benchmark CI

<!-- Author: @jaimergp -->
<!-- Last updated: 2021.07.06 -->
<!-- Describes the work done as part of https://github.com/scikit-image/scikit-image/pull/5424 -->

## How it works

The `asv` suite can be run for any PR on GitHub Actions (check workflow `.github/workflows/benchmarks.yml`) by adding a `run-benchmark` label to said PR. This will trigger a job that will run the benchmarking suite for the current PR head (merged commit) against the PR base (usually `main`).

We use `asv continuous` to run the job, which runs a relative performance measurement. This means that there's no state to be saved and that regressions are only caught in terms of performance ratio (absolute numbers are available but they are not useful since we do not use stable hardware over time). `asv continuous` will:

* Compile `scikit-image` for _both_ commits. We use `ccache` to speed up the process, and `mamba` is used to create the build environments.
* Run the benchmark suite for both commits, _twice_  (since `processes=2` by default).
* Generate a report table with performance ratios:
    * `ratio=1.0` -> performance didn't change.
    * `ratio<1.0` -> PR made it slower.
    * `ratio>1.0` -> PR made it faster.

Due to the sensitivity of the test, we cannot guarantee that false positives are not produced. In practice, values between `(0.7, 1.5)` are to be considered part of the measurement noise. When in doubt, running the benchmark suite one more time will provide more information about the test being a false positive or not.

## Running the benchmarks on GitHub Actions

1. On a PR, add the label `run-benchmark`.
2. The CI job will be started. Checks will appear in the usual dashboard panel above the comment box.
3. If more commits are added, the label checks will be grouped with the last commit checks _before_ you added the label.
4. Alternatively, you can always go to the `Actions` tab in the repo and [filter for `workflow:Benchmark`](https://github.com/scikit-image/scikit-image/actions?query=workflow%3ABenchmark). Your username will be assigned to the `actor` field, so you can also filter the results with that if you need it.

## The artifacts

The CI job will also generate an artifact. This is the `.asv/results` directory compressed in a zip file. Its contents include:

* `fv-xxxxx-xx/`. A directory for the machine that ran the suite. It contains three files:
    * `<baseline>.json`, `<contender>.json`: the benchmark results for each commit, with stats.
    * `machine.json`: details about the hardware.
* `benchmarks.json`: metadata about the current benchmark suite.
* `benchmarks.log`: the CI logs for this run.
* This README.

## Re-running the analysis

Although the CI logs should be enough to get an idea of what happened (check the table at the end), one can use `asv` to run the analysis routines again.

1. Uncompress the artifact contents in the repo, under `.asv/results`. This is, you should see `.asv/results/benchmarks.log`, not `.asv/results/something_else/benchmarks.log`. Write down the machine directory name for later.
2. Run `asv show` to see your available results. You will see something like this:

```
$> asv show

Commits with results:

Machine    : Jaimes-MBP
Environment: conda-py3.9-cython-numpy1.20-scipy

    00875e67

Machine    : fv-az95-499
Environment: conda-py3.7-cython-numpy1.17-pooch-scipy

    8db28f02
    3a305096
```

3. We are interested in the commits for `fv-az95-499` (the CI machine for this run). We can compare them with `asv compare` and some extra options. `--sort ratio` will show largest ratios first, instead of alphabetical order. `--split` will produce three tables: improved, worsened, no changes. `--factor 1.5` tells `asv` to only complain if deviations are above a 1.5 ratio. `-m` is used to indicate the machine ID (use the one you wrote down in step 1). Finally, specify your commit hashes: baseline first, then contender!

```
$> asv compare --sort ratio --split --factor 1.5 -m fv-az95-499 8db28f02 3a305096

Benchmarks that have stayed the same:

       before           after         ratio
     [8db28f02]       [3a305096]
     <ci-benchmark-check~9^2>
              n/a              n/a      n/a  benchmark_restoration.RollingBall.time_rollingball_ndim
      1.23±0.04ms       1.37±0.1ms     1.12  benchmark_transform_warp.WarpSuite.time_to_float64(<class 'numpy.float64'>, 128, 3)
       5.07±0.1μs       5.59±0.4μs     1.10  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (192, 192, 192), (192, 192, 192))
      1.23±0.02ms       1.33±0.1ms     1.08  benchmark_transform_warp.WarpSuite.time_same_type(<class 'numpy.float32'>, 128, 3)
       9.45±0.2ms       10.1±0.5ms     1.07  benchmark_rank.Rank3DSuite.time_3d_filters('majority', (32, 32, 32))
       23.0±0.9ms         24.6±1ms     1.07  benchmark_interpolation.InterpolationResize.time_resize((80, 80, 80), 0, 'symmetric', <class 'numpy.float64'>, True)
         38.7±1ms         41.1±1ms     1.06  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (2048, 2048), (192, 192, 192))
       4.97±0.2μs       5.24±0.2μs     1.05  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (2048, 2048), (2048, 2048))
       4.21±0.2ms       4.42±0.3ms     1.05  benchmark_rank.Rank3DSuite.time_3d_filters('gradient', (32, 32, 32))

...
```

If you want more details on a specific test, you can use `asv show`. Use `-b pattern` to filter which tests to show, and then specify a commit hash to inspect:

```
$> asv show -b time_to_float64 8db28f02

Commit: 8db28f02 <ci-benchmark-check~9^2>

benchmark_transform_warp.WarpSuite.time_to_float64 [fv-az95-499/conda-py3.7-cython-numpy1.17-pooch-scipy]
  ok
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
  --                                                                N / order
  --------------- --------------------------------------------------------------------------------------------------------------
      dtype_in       128 / 0     128 / 1      128 / 3     1024 / 0    1024 / 1    1024 / 3    4096 / 0    4096 / 1    4096 / 3
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
    numpy.uint8    2.56±0.09ms   523±30μs   1.28±0.05ms   130±3ms     28.7±2ms    81.9±3ms   2.42±0.01s   659±5ms    1.48±0.01s
    numpy.uint16   2.48±0.03ms   530±10μs   1.28±0.02ms   130±1ms    30.4±0.7ms   81.1±2ms    2.44±0s     653±3ms    1.47±0.02s
   numpy.float32    2.59±0.1ms   518±20μs   1.27±0.01ms   127±3ms     26.6±1ms    74.8±2ms   2.50±0.01s   546±10ms   1.33±0.02s
   numpy.float64   2.48±0.04ms   513±50μs   1.23±0.04ms   134±3ms     30.7±2ms    85.4±2ms   2.55±0.01s   632±4ms    1.45±0.01s
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
  started: 2021-07-06 06:14:36, duration: 1.99m
```

## Other details

### Skipping slow or demanding tests

To minimize the time required to run the full suite, we trimmed the parameter matrix in some cases and, in others, directly skipped tests that ran for too long or require too much memory. Unlike `pytest`, `asv` does not have a notion of marks. However, you can `raise NotImplementedError` in the setup step to skip a test. In that vein, a new private function is defined at `benchmarks.__init__`: `_skip_slow`. This will check if the `ASV_SKIP_SLOW` environment variable has been defined. If set to `1`, it will raise `NotImplementedError` and skip the test. To implement this behavior in other tests, you can add the following attribute:

```python
from . import _skip_slow  # this function is defined in benchmarks.__init__

def time_something_slow():
    pass

time_something.setup = _skip_slow
```
This directory contains meta-files related to the Lewiner marching cubes
algorithm. These are not used by the algorithm, but can be convenient
for development/maintenance:
    
* MarchingCubes.cpp - the original algorithm, this is ported to Cython
* LookupTable.h - the original LUTs, these are ported to Python
* createluts.py - scrip to generate Python luts from the .h file
* visual_test.py - script to compare visual results of marchingcubes algorithms
# Version 0.16

- The following functions are deprecated and will be removed in 0.18:
  ``skimage.measure.compare_mse``,
  ``skimage.measure.compare_nrmse``,
  ``skimage.measure.compare_pnsr``,
  ``skimage.measure.compare_ssim``
  Their functionality still exists, but under the new ``skimage.metrics``
  submodule under different names.
- Additionally, three new functions have been added to ``skimage.metrics``:
  ``skimage.metrics.variation_of_information``
  ``skimage.metrics.adapted_rand_error``
  ``skimage.metrics.contingency_table``
- A new example of plotting these evaluation metrics has been added to the docs.

# Version 0.15

- ``skimage.feature.canny`` now uses a more accurate Gaussian filter
  internally; output values will be different from 0.14.
- ``skimage.filters.threshold_niblack`` and
  ``skimage.filters.threshold_sauvola``
  now accept a tuple as ``window_size`` besides integers.

# Version 0.14

- ``skimage.filters.gaussian_filter`` has been removed. Use
  ``skimage.filters.gaussian`` instead.
- ``skimage.filters.gabor_filter`` has been removed. Use
  ``skimage.filters.gabor`` instead.
- The old syntax support for ``skimage.transform.integrate`` has been removed.
- The ``normalise`` parameter of ``skimage.feature.hog`` was removed due to
  incorrect behavior: it only applied a square root instead of a true
  normalization. If you wish to duplicate the old behavior, set
  ``transform_sqrt=True``.
- ``skimage.measure.structural_similarity`` has been removed. Use
  ``skimage.measure.compare_ssim`` instead.
- In ``skimage.measure.compare_ssim``, the `dynamic_range` has been removed in
  favor of '`data_range`.
- In ``skimage.restoration.denoise_bilateral``, the `sigma_range` kwarg has
  been removed in favor of `sigma_color`.
- ``skimage.measure.marching_cubes`` has been removed in favor of
  ``skimage.measure.marching_cubes_lewiner``.
- ``ntiles_*`` parameters have been removed from
  ``skimage.exposure.equalize_adapthist``. Use ``kernel_size`` instead.
- ``skimage.restoration.nl_means_denoising`` has been removed in
  favor of ``skimage.restoration.denoise_nl_means``.
- ``skimage.measure.LineModel`` has been removed in favor of
  ``skimage.measure.LineModelND``.
- In ``skimage.feature.hog`` visualise has been changed to visualize.
- `freeimage` plugin of ``skimage.io`` has been removed.

# Version 0.13

- `skimage.filter` has been removed. Use `skimage.filters` instead.
- `skimage.filters.canny` has been removed.
  `canny` is available only from `skimage.feature` now.
- Deprecated filters `hsobel`, `vsobel`, `hscharr`, `vscharr`, `hprewitt`,
  `vprewitt`, `roberts_positive_diagonal`, `roberts_negative_diagonal` have
  been removed from `skimage.filters.edges`.
- The `sigma` parameter of `skimage.filters.gaussian` and the `selem` parameter
  of `skimage.filters.median` have been made optional, with default
  values.
- The `clip_negative` parameter of `skimage.util.dtype_limits` is now set
  to `None` by default, equivalent to `True`, the former value. In version
  0.15, will be set to `False`.
- The `circle` parameter of `skimage.transform.radon` and `skimage.transform.iradon`
  are now set to `None` by default, equivalent to `False`, the former value. In version
  0.15, will be set to `True`.
- Parameters ``ntiles_x``, ``ntiles_y`` have been removed from
  ``skimage.exposure.equalize_adapthist``.
- The ``freeimage`` io plugin is no longer supported, and will use ``imageio``
  instead.  We will completely remove the ``freeimage`` plugin in Version 0.14.

# Version 0.12

- ``equalize_adapthist`` now takes a ``kernel_size`` keyword argument, replacing
  the ``ntiles_*`` arguments.
- The functions ``blob_dog``, ``blob_log`` and ``blob_doh`` now return float
  arrays instead of integer arrays.
- ``transform.integrate`` now takes lists of tuples instead of integers
  to define the window over which to integrate.
- `reverse_map` parameter in `skimage.transform.warp` has been removed.
- `enforce_connectivity` in `skimage.segmentation.slic` defaults to ``True``.
- `skimage.measure.fit.BaseModel._params`,
  `skimage.transform.ProjectiveTransform._matrix`,
  `skimage.transform.PolynomialTransform._params`,
  `skimage.transform.PiecewiseAffineTransform.affines_*` attributes
  have been removed.
- `skimage.filters.denoise_*` have moved to `skimage.restoration.denoise_*`.
- `skimage.data.lena` has been removed.

# Version 0.11

- The ``skimage.filter`` subpackage has been renamed to ``skimage.filters``.
- Some edge detectors returned values greater than 1--their results are now
  appropriately scaled with a factor of ``sqrt(2)``.

# Version 0.10

- Removed ``skimage.io.video`` functionality due to broken gstreamer bindings

# Version 0.9

- No longer wrap ``imread`` output in an ``Image`` class
- Change default value of `sigma` parameter in ``skimage.segmentation.slic``
  to 0
- ``hough_circle`` now returns a stack of arrays that are the same size as the
  input image. Set the ``full_output`` flag to True for the old behavior.
- The following functions were deprecated over two releases:
  `skimage.filter.denoise_tv_chambolle`,
  `skimage.morphology.is_local_maximum`, `skimage.transform.hough`,
  `skimage.transform.probabilistic_hough`,`skimage.transform.hough_peaks`.
  Their functionality still exists, but under different names.

Version 0.4
-----------
- Switch mask and radius arguments for ``median_filter``

Version 0.3
-----------
- Remove ``as_grey``, ``dtype`` keyword from ImageCollection
- Remove ``dtype`` from imread
- Generalise ImageCollection to accept a load_func
User Guide
==========

```{toctree}
---
maxdepth: 2
---

user_guide/getting_started
user_guide/numpy_images
user_guide/data_types
user_guide/plugins
user_guide/video
user_guide/visualization
user_guide/transforming_image_data
user_guide/geometrical_transform
user_guide/tutorials
user_guide/getting_help
```
(core_dev)=

Core Developer Guide
====================

Welcome, new core developer!  The core team appreciate the quality of
your work, and enjoy working with you; we have therefore invited you
to join us.  Thank you for your numerous contributions to the project
so far.

This document offers guidelines for your new role.  First and
foremost, you should familiarize yourself with the project's
{doc}`mission, vision, and values <values>`.  When in
doubt, always refer back here.

As a core team member, you gain the responsibility of shepherding
other contributors through the review process; here are some
guidelines.

All Contributors Are Treated The Same
-------------------------------------

You now have the ability to push changes directly to the main
branch, but should never do so; instead, continue making pull requests
as before and in accordance with the 
{doc}`general contributor guide <contribute>`.

As a core contributor, you gain the ability to merge or approve
other contributors' pull requests.  Much like nuclear launch keys, it
is a shared power: you must merge *only after* another core has
approved the pull request, *and* after you yourself have carefully
reviewed it.  (See {ref}`sec:reviewing` and especially
{ref}`sec:understand` below.) To ensure a clean git history,
use GitHub's [Squash and Merge][gh_sqmrg]
feature to merge, unless you have a good reason not to do so.

[gh_sqmrg]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/merging-a-pull-request#merging-a-pull-request-on-github

(sec:reviewing)=
Reviewing
---------

### How to Conduct A Good Review

*Always* be kind to contributors. Nearly all of `scikit-image` is
volunteer work, for which we are tremendously grateful. Provide
constructive criticism on ideas and implementations, and remind
yourself of how it felt when your own work was being evaluated as a
novice.

`scikit-image` strongly values mentorship in code review.  New users
often need more handholding, having little to no git
experience. Repeat yourself liberally, and, if you don’t recognize a
contributor, point them to our development guide, or other GitHub
workflow tutorials around the web. Do not assume that they know how
GitHub works (e.g., many don't realize that adding a commit
automatically updates a pull request). Gentle, polite, kind
encouragement can make the difference between a new core developer and
an abandoned pull request.

When reviewing, focus on the following:

1. **API:** The API is what users see when they first use
   `scikit-image`. APIs are difficult to change once released, so
   should be simple, [functional][wiki_functional] (i.e. not
   carry state), consistent with other parts of the library, and
   should avoid modifying input variables.  Please familiarize
   yourself with the project's [deprecation policy][dep_pol]

2. **Documentation:** Any new feature should have a gallery
   example, that not only illustrates but explains it.

3. **The algorithm:** You should understand the code being modified or
   added before approving it.  (See {ref}`sec:understand`
   below.) Implementations should do what they claim,
   and be simple, readable, and efficient.

4. **Tests:** All contributions to the library *must* be tested, and
   each added line of code should be covered by at least one test. Good
   tests not only execute the code, but explores corner cases.  It is tempting
   not to review tests, but please do so.

[wiki_functional]: https://en.wikipedia.org/wiki/Functional_programming
[dep_pol]: https://scikit-image.org/docs/dev/contribute.html#deprecation-cycle

Other changes may be *nitpicky*: spelling mistakes, formatting,
etc. Do not ask contributors to make these changes, and instead
make the changes by [pushing to their branch][gh_push]
or using GitHub’s [suggestion][gh_suggest] [feature][gh_feedback].
(The latter is preferred because it gives the contributor a choice in
whether to accept the changes.)

[gh_push]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/committing-changes-to-a-pull-request-branch-created-from-a-fork
[gh_suggest]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/commenting-on-a-pull-request
[gh_feedback]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/incorporating-feedback-in-your-pull-request

Our default merge policy is to squash all PR commits into a single
commit. Users who wish to bring the latest changes from ``main``
into their branch should be advised to merge, not to rebase.  Even
when merge conflicts arise, don’t ask for a rebase unless you know
that a contributor is experienced with git. Instead, rebase the branch
yourself, force-push to their branch, and advise the contributor on
how to force-pull.  If the contributor is no longer active, you may
take over their branch by submitting a new pull request and closing
the original. In doing so, ensure you communicate that you are not
throwing the contributor's work away!

Please add a note to a pull request after you push new changes; GitHub
does not send out notifications for these.

(sec:understand)=
### Merge Only Changes You Understand

*Long-term maintainability* is an important concern.  Code doesn't
merely have to *work*, but should be *understood* by multiple core
developers.  Changes will have to be made in the future, and the
original contributor may have moved on.

Therefore, *do not merge a code change unless you understand it*. Ask
for help freely: we have a long history of consulting community
members, or even external developers, for added insight where needed,
and see this as a great learning opportunity.

While we collectively "own" any patches (and bugs!) that become part
of the code base, you are vouching for changes you merge.  Please take
that responsibility seriously.

In practice, if you are the second core developer reviewing and approving a
given pull request, you typically merge it (again, using GitHub's Squash and
Merge feature) in the wake of your approval. What are the exceptions to this
process? If the pull request has been particularly controversial or the
subject of much debate (e.g., involving API changes), then you would want to
wait a few days before merging. This waiting time gives others a chance to
speak up in case they are not fine with the current state of the pull request.
Another exceptional situation is one where the first approving review happened
a long time ago and many changes have taken place in the meantime.

When squashing commits GitHub concatenates all commit messages.
Please edit the resulting message so that it gives a concise, tidy
overview of changes. For example, you may want to grab the
description from the PR itself, and delete lines such as "pep8 fix",
"apply review comments", etc. Please retain all Co-authored-by
entries.

Closing issues and pull requests
--------------------------------

Sometimes, an issue must be closed that was not fully resolved. This can be
for a number of reasons:

- the person behind the original post has not responded to calls for
  clarification, and none of the core developers have been able to reproduce
  their issue;
- fixing the issue is difficult, and it is deemed too niche a use case to
  devote sustained effort or prioritize over other issues; or
- the use case or feature request is something that core developers feel
  does not belong in scikit-image,

among others. Similarly, pull requests sometimes need to be closed without
merging, because:

- the pull request implements a niche feature that we consider not worth the
  added maintenance burden;
- the pull request implements a useful feature, but requires significant
  effort to bring up to scikit-image's standards, and the original
  contributor has moved on, and no other developer can be found to make the
  necessary changes; or
- the pull request makes changes that do not align with our values, such as
  increasing the code complexity of a function significantly to implement a
  marginal speedup,

among others.

All these may be valid reasons for closing, but we must be wary not to alienate
contributors by closing an issue or pull request without an explanation. When
closing, your message should:

- explain clearly how the decision was made to close. This is particularly
  important when the decision was made in a community meeting, which does not
  have as visible a record as the comments thread on the issue itself;
- thank the contributor(s) for their work; and
- provide a clear path for the contributor or anyone else to appeal the
  decision.

These points help ensure that all contributors feel welcome and empowered to
keep contributing, regardless of the outcome of past contributions.

Further resources
-----------------

As a core member, you should be familiar with community and developer
resources such as:

-  Our {doc}`contributor guide <contribute>`
-  Our [community guidelines](https://scikit-image.org/community_guidelines.html)
-  [PEP8](https://www.python.org/dev/peps/pep-0008/) for Python style
-  [PEP257](https://www.python.org/dev/peps/pep-0257/) and the
   [NumPy documentation guide][numpydoc]
   for docstrings. (NumPy docstrings are a superset of PEP257. You
   should read both.)
-  The scikit-image [tag on StackOverflow][so_tag]
-  The scikit-image [tag on forum.image.sc](https://forum.image.sc/tags/scikit-image)
-  Our [developer forum][ml]
-  Our [chat room](https://skimage.zulipchat.com/)

[numpydoc]: https://docs.scipy.org/doc/numpy/docs/howto_document.html
[so_tag]: https://stackoverflow.com/questions/tagged/scikit-image
[ml]: https://discuss.scientific-python.org/c/contributor/skimage

You are not required to monitor all of the social resources.

Inviting New Core Members
-------------------------

Any core member may nominate other contributors to join the core team.
Nominations happen on a private email list,
<skimage-core@python.org>. As of this writing, there is no hard-and-fast
rule about who can be nominated; at a minimum, they should have: been
part of the project for at least six months, contributed
significant changes of their own, contributed to the discussion and
review of others' work, and collaborated in a way befitting our
community values.

Contribute To This Guide!
-------------------------

This guide reflects the experience of the current core developers.  We
may well have missed things that, by now, have become second
nature—things that you, as a new team member, will spot more easily.
Please ask the other core developers if you have any questions, and
submit a pull request with insights gained.

Conclusion
----------

We are excited to have you on board!  We look forward to your
contributions to the code base and the community.  Thank you in
advance!
Glossary
========

Work in progress

```{glossary}

array 
    Numerical array, provided by the {class}`numpy.ndarray` object. In
    ``scikit-image``, images are NumPy arrays with dimensions that
    correspond to spatial dimensions of the image, and color channels for
    color images. See {ref}`numpy`.

channel
    Typically used to refer to a single color channel in a color image. RGBA
    images have an additional alpha (transparency) channel. Functions use a
    ``channel_axis`` argument to specify which axis of an array corresponds
    to channels. Images without channels are indicated via
    ``channel_axis=None``. Aside from the functions in ``skimage.color``, most
    functions with a ``channel_axis`` argument just apply the same operation
    across each channel. In this case, the "channels" do not strictly need to
    represent color or alpha information, but may be any generic batch
    dimension over which to operate.

circle
    The perimeter of a {term}`disk`.

contour
    Curve along which a 2-D image has a constant value. The interior
    (resp. exterior) of the contour has values greater (resp. smaller)
    than the contour value.

contrast
    Differences of intensity or color in an image, which make objects
    distinguishable. Several functions to manipulate the contrast of an
    image are available in {mod}`skimage.exposure`. See {ref}`exposure`.

disk
    A filled-in {term}`circle`.

float
    Representation of real numbers, for example as {obj}`np.float32` or
    {obj}`np.float64`. See {ref}`data_types`. Some operations on images
    need a float datatype (such as multiplying image values with
    exponential prefactors in {func}`filters.gaussian`), so that 
    images of integer type are often converted to float type internally. Also
    see {term}`int` values.

float values
    See {term}`float`.

histogram
    For an image, histogram of intensity values, where the range of
    intensity values is divided into bins and the histogram counts how
    many pixel values fall in each bin. See
    {func}`exposure.histogram`.

int
    Representation of integer numbers, which can be signed or not, and
    encoded on one, two, four or eight bytes according to the maximum value
    which needs to be represented. In ``scikit-image``, the most common
    integer types are {obj}`np.int64` (for large integer values) and
    {obj}`np.uint8` (for small integer values, typically images of labels
    with less than 255 labels). See {ref}`data_types`. 

int values
    See {term}`int`.

iso-valued contour
    See {term}`contour`.

labels
    An image of labels is of integer type, where pixels with the same
    integer value belong to the same object. For example, the result of a
    segmentation is an image of labels. {func}`measure.label` labels
    connected components of a binary image and returns an image of
    labels. Labels are usually contiguous integers, and
    {func}`segmentation.relabel_sequential` can be used to relabel
    arbitrary labels to sequential (contiguous) ones.

label image
    See {term}`labels`.

pixel
    Smallest element of an image. An image is a grid of pixels, and the
    intensity of each pixel is variable. A pixel can have a single
    intensity value in grayscale images, or several channels for color
    images. In ``scikit-image``, pixels are the individual elements of
    ``numpy arrays`` (see {ref}`numpy`).
    Also see {term}`voxel`.

segmentation
    Partitioning an image into multiple objects (segments), for
    example an object of interest and its background. The output of a
    segmentation is typically an image of {term}`labels`, where
    the pixels of different objects have been attributed different
    integer labels. Several segmentation algorithms are available in
    {mod}`skimage.segmentation`.

voxel
    {term}`pixel` (smallest element of an image) of a
    three-dimensional image.
```
# Our mission

scikit-image aims to be the reference library for scientific image analysis in
Python. We accomplish this by:

- being **easy to use and install**. We are careful in taking on new
  dependencies, and sometimes cull existing ones, or make them optional. All
  functions in our API have thorough docstrings clarifying expected inputs and
  outputs.
- providing a **consistent API**. Conceptually identical arguments have the
  same name and position in a function signature.
- **ensuring correctness**. Test coverage is close to 100% and code is reviewed by
  at least two core developers before being included in the library.
- **caring for users’ data**. We have a [functional API][functional] and don't modify
  input arrays unless explicitly directed to do so.
- promoting **education in image processing**, with extensive pedagogical
  documentation.

(sec:values)=
Our values
----------

- We are inclusive. We continue to welcome and mentor newcomers who are
  making their first contribution.
- We are community-driven. Decisions about the API and features are driven by
  our users' requirements, not by the whims of the core team. (See
  {ref}`governance`.)
- We serve scientific applications primarily, over “consumer” image editing in
  the vein of Photoshop or GIMP. This often means prioritizing n-dimensional
  data support, and rejecting implementations of “flashy” filters that have
  little scientific value.
- We value simple, readable implementations over getting every last ounce of
  performance. Readable code that is easy to understand, for newcomers and
  maintainers alike, makes it easier to contribute new code as well as prevent
  bugs. This means that we will prefer a 20% slowdown if it reduces lines of
  code two-fold, for example.
- We value education and documentation. All functions should have NumPy-style
  [docstrings][numpydoc], preferably with examples, as well as gallery
  examples that showcase how that function is used in a scientific application.
  Core developers take an active role in finishing documentation examples.
- We don't do magic. We use NumPy arrays instead of fancy façade objects
  [^np], and we prefer to educate users rather than make decisions on their
  behalf.  This does not preclude [sensible defaults][defaults].

This document
-------------

Much in the same way that the [Zen of Python][zen] and PEP8 guide style and
implementation details in most Python code, this guide is meant to guide
decisions about the future of scikit-image, be it in terms of code style,
whether to accept new functionality, or whether to take on new dependencies,
among other things.

References
----------

To find out more about the history of this document, please read the following:

- [Original blog post][blog]
- [The GitHub issue][issue]
- [The image.sc forum post][forum]
- [The SKIP GitHub pull request][skip_pr]

% Links
% -----

[functional]: https://en.wikipedia.org/wiki/Functional_programming
[blog]: https://ilovesymposia.com/2018/07/13/the-road-to-scikit-image-1-0/
[numpydoc]: https://docs.scipy.org/doc/numpy/docs/howto_document.html
[defaults]: https://forum.image.sc/t/request-for-comment-road-to-scikit-image-1-0/20099/4
[zen]: https://www.python.org/dev/peps/pep-0020/
[issue]: https://github.com/scikit-image/scikit-image/issues/3263
[forum]: https://forum.image.sc/t/request-for-comment-road-to-scikit-image-1-0/20099
[skip_pr]: https://github.com/scikit-image/scikit-image/pull/3585
[cc0]: https://creativecommons.org/publicdomain/zero/1.0/
[ccby]: https://dancohen.org/2013/11/26/cc0-by/

Copyright
---------

This document is dedicated to the public domain with the Creative Commons CC0
[license][cc0]. Attribution to this source is encouraged where appropriate, as per
[CC0+BY][ccby].

[^np]: The use of NumPy arrays was the most supported of the statement's
       components, together with the points about inclusivity, mentorship, and
       documentation. We had +1s from Mark Harfouche, Royi Avital, and Greg Lee,
       among others.
scikit-image: image processing in Python
========================================

[scikit-image](https://scikit-image.org) builds on
[scipy.ndimage](https://docs.scipy.org/doc/scipy/reference/ndimage.html) to
provide a versatile set of image processing routines in
[Python](https://www.python.org).

This library is developed by its community, and
{doc}`contributions <contribute>` are most welcome!

Read about our {doc}`mission, vision, and values <values>` and how we
{doc}`govern <skips/1-governance>` the project.

Major proposals to the project are documented in {doc}`SKIPs <skips/index>`.

Homepage
--------
<https://scikit-image.org>

User Forum
----------
<https://forum.image.sc/tags/scikit-image>

Developer Forum
---------------
<https://discuss.scientific-python.org/c/contributor/skimage>

Source, bugs, and patches
-------------------------
<https://github.com/scikit-image/scikit-image>

Contact
-------
Stefan van der Walt `<stefanv at berkeley.edu>`
scikit-image Code of Conduct
============================

scikit-image adopts the [SciPy code of conduct][scipy_coc].

[scipy_coc]: https://github.com/scipy/scipy/blob/master/doc/source/dev/conduct/code_of_conduct.rst

Introduction
------------

This code of conduct applies to all spaces managed by the scikit-image project,
including all public and private forums, issue trackers, wikis, blogs,
Twitter, and any other communication channel used by our community.  The scikit-image
project does not organise in-person events, however events related to our
community should have a code of conduct similar in spirit to this one.

This code of conduct should be honored by everyone who participates in
the scikit-image community formally or informally, or claims any affiliation with the
project, in any project-related activities and especially when representing the
project, in any role.

This code is not exhaustive or complete. It serves to distill our common
understanding of a collaborative, shared environment and goals. Please try to
follow this code in spirit as much as in letter, to create a friendly and
productive environment that enriches the surrounding community.


Specific Guidelines
-------------------

We strive to:

1. Be open. We invite anyone to participate in our community. We prefer to use
   public methods of communication for project-related messages, unless
   discussing something sensitive. This applies to messages for help or
   project-related support, too; not only is a public support request much more
   likely to result in an answer to a question, it also ensures that any
   inadvertent mistakes in answering are more easily detected and corrected.

2. Be empathetic, welcoming, friendly, and patient. We work together to resolve
   conflict, and assume good intentions. We may all experience some frustration
   from time to time, but we do not allow frustration to turn into a personal
   attack. A community where people feel uncomfortable or threatened is not a
   productive one.

3. Be collaborative. Our work will be used by other people, and in turn we will
   depend on the work of others. When we make something for the benefit of the
   project, we are willing to explain to others how it works, so that they can
   build on the work to make it even better. Any decision we make will affect
   users and colleagues, and we take those consequences seriously when making
   decisions.

4. Be inquisitive. Nobody knows everything! Asking questions early avoids many
   problems later, so we encourage questions, although we may direct them to
   the appropriate forum. We will try hard to be responsive and helpful.

5. Be careful in the words that we choose.  We are careful and respectful in
   our communication and we take responsibility for our own speech. Be kind to
   others. Do not insult or put down other participants.  We will not accept
   harassment or other exclusionary behaviour, such as:

   - Violent threats or language directed against another person.
   - Sexist, racist, or otherwise discriminatory jokes and language.
   - Posting sexually explicit or violent material.
   - Posting (or threatening to post) other people's personally identifying information ("doxing").
   - Sharing private content, such as emails sent privately or non-publicly,
     or unlogged forums such as IRC channel history, without the sender's consent.
   - Personal insults, especially those using racist or sexist terms.
   - Unwelcome sexual attention.
   - Excessive profanity. Please avoid swearwords; people differ greatly in their sensitivity to swearing.
   - Repeated harassment of others. In general, if someone asks you to stop, then stop.
   - Advocating for, or encouraging, any of the above behaviour.


Diversity Statement
-------------------

The scikit-image project welcomes and encourages participation by everyone. We are
committed to being a community that everyone enjoys being part of. Although
we may not always be able to accommodate each individual's preferences, we try
our best to treat everyone kindly.

No matter how you identify yourself or how others perceive you: we welcome you.
Though no list can hope to be comprehensive, we explicitly honour diversity in:
age, culture, ethnicity, genotype, gender identity or expression, language,
national origin, neurotype, phenotype, political beliefs, profession, race,
religion, sexual orientation, socioeconomic status, subculture and technical
ability, to the extent that these do not conflict with this code of conduct.


Though we welcome people fluent in all languages, scikit-image development is
conducted in English.

Standards for behaviour in the scikit-image community are detailed in the Code of
Conduct above. Participants in our community should uphold these standards
in all their interactions and help others to do so as well (see next section).


Reporting Guidelines
--------------------

We know that it is painfully common for internet communication to start at or
devolve into obvious and flagrant abuse.  We also recognize that sometimes
people may have a bad day, or be unaware of some of the guidelines in this Code
of Conduct. Please keep this in mind when deciding on how to respond to a
breach of this Code.

For clearly intentional breaches, report those to the Code of Conduct committee
(see below). For possibly unintentional breaches, you may reply to the person
and point out this code of conduct (either in public or in private, whatever is
most appropriate). If you would prefer not to do that, please feel free to
report to the Code of Conduct Committee directly, or ask the Committee for
advice, in confidence.

You can report issues to the scikit-image Code of Conduct committee, at
<skimage-conduct@groups.io>. Currently, the committee consists of:

- Emmanuelle Gouillart <emmanuelle.gouillard@normalesup.org> (chair)
- Carol Willing <willingc@gmail.com>
- François Boulogne <devel@sciunto.org>
- Egor Panfilov <egor.v.panfilov@gmail.com>

If your report involves any members of the committee, or if they feel they have
a conflict of interest in handling it, then they will recuse themselves from
considering your report. Alternatively, if for any reason you feel
uncomfortable making a report to the committee, then you can also contact:

- Any member of the {ref}`scikit-image Steering Council <governance>`, or
- Senior [NumFOCUS staff](https://numfocus.org/code-of-conduct#persons-responsible): 
  <conduct@numfocus.org>.


Incident reporting resolution & Code of Conduct enforcement
-----------------------------------------------------------

*This section summarizes the most important points, more details can be found
in the* {ref}`reporting manual <CoC_reporting_manual>`.

We will investigate and respond to all complaints. The scikit-image Code of Conduct
Committee and the scikit-image Steering Committee (if involved) will protect the
identity of the reporter, and treat the content of complaints as confidential
(unless the reporter agrees otherwise).

In case of severe and obvious breaches, e.g. personal threat or violent, sexist
or racist language, we will immediately disconnect the originator from scikit-image
communication channels; please see the manual for details.

In cases not involving clear severe and obvious breaches of this code of
conduct, the process for acting on any received code of conduct violation
report will be:

1. acknowledge report is received
2. reasonable discussion/feedback
3. mediation (if feedback didn't help, and only if both reporter and reportee agree to this)
4. enforcement via transparent decision (see {ref}`CoC_resolutions`) by the
   Code of Conduct Committee

The committee will respond to any report as soon as possible, and at most
within 72 hours.


Endnotes
--------

We are thankful to the groups behind the following documents, from which we
drew content and inspiration:

- [The SciPy project](https://www.scipy.org/)
- [The Apache Foundation Code of Conduct](https://www.apache.org/foundation/policies/conduct.html)
- [The Contributor Covenant](https://www.contributor-covenant.org/version/1/4/code-of-conduct)
- [Jupyter Code of Conduct](https://github.com/jupyter/governance/tree/master/conduct)
- [Open Source Guides - Code of Conduct](https://opensource.guide/code-of-conduct/)
(skip_list)=
# scikit-image proposals (SKIPS)

SKIPs document any major changes or proposals to the scikit-image project.

## Sections

```{toctree}
---
maxdepth: 1
---
0-skip-process
1-governance
2-values
3-transition-to-v1
```

```{toctree}
---
hidden:
---
template
```
.. _howto_contribute:

How to contribute to scikit-image
=================================

Developing Open Source is great fun! Join us on the `scikit-image
developer forum <https://discuss.scientific-python.org/c/contributor/skimage>`_ and tell us
which of the following challenges you'd like to solve.

* Mentoring is available for those new to scientific programming in Python.
* If you're looking for something to implement or to fix, you can browse the
  `open issues on GitHub <https://github.com/scikit-image/scikit-image/issues?q=is%3Aopen>`__.
* The technical detail of the `development process`_ is summed up below.
  Refer to the :doc:`gitwash <gitwash/index>` for a step-by-step tutorial.

.. contents::
   :local:

Development process
-------------------

Here's the long and short of it:

1. If you are a first-time contributor:

   * Go to `https://github.com/scikit-image/scikit-image
     <https://github.com/scikit-image/scikit-image>`_ and click the
     "fork" button to create your own copy of the project.

   * Clone the project to your local computer::

      git clone https://github.com/your-username/scikit-image.git

   * Change the directory::

      cd scikit-image

   * Add the upstream repository::

      git remote add upstream https://github.com/scikit-image/scikit-image.git

   * Now, you have remote repositories named:

     - ``upstream``, which refers to the ``scikit-image`` repository
     - ``origin``, which refers to your personal fork

.. note::

    Although our code is hosted on `github
    <https://github.com/scikit-image/>`_, our dataset is stored on `gitlab
    <https://gitlab.com/scikit-image/data>`_ and fetched with `pooch
    <https://github.com/fatiando/pooch>`_. New data must be submitted on
    gitlab. Once merged, the data registry ``skimage/data/_registry.py``
    in the main codebase on github must be updated.

2. Develop your contribution:

   * Pull the latest changes from upstream::

      git checkout main
      git pull upstream main

   * Create a branch for the feature you want to work on. Use a sensible name,
     such as 'transform-speedups'::

      git checkout -b transform-speedups

   * Commit locally as you progress (with ``git add`` and ``git commit``).
     Please write `good commit messages
     <https://vxlabs.com/software-development-handbook/#good-commit-messages>`_.

3. To submit your contribution:

   * Push your changes back to your fork on GitHub::

      git push origin transform-speedups

   * Enter your GitHub username and password (repeat contributors or advanced
     users can remove this step by `connecting to GitHub with SSH
     <https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_).

   * Go to GitHub. The new branch will show up with a green "pull request"
     button -- click it.

   * If you want, post on the `developer forum
     <https://discuss.scientific-python.org/c/contributor/skimage>`_ to explain your changes or
     to ask for review.

For a more detailed discussion, read these :doc:`detailed documents
<gitwash/index>` on how to use Git with ``scikit-image`` (:ref:`using-git`).

4. Review process:

   * Reviewers (the other developers and interested community members) will
     write inline and/or general comments on your pull request (PR) to help
     you improve its implementation, documentation, and style.  Every single
     developer working on the project has their code reviewed, and we've come
     to see it as a friendly conversation from which we all learn and the
     overall code quality benefits.  Therefore, please don't let the review
     discourage you from contributing: its only aim is to improve the quality
     of the project, not to criticize (we are, after all, very grateful for the
     time you're donating!).

   * To update your pull request, make your changes on your local repository
     and commit. As soon as those changes are pushed up (to the same branch as
     before) the pull request will update automatically.

   * Continuous integration (CI) services are triggered after each pull request
     submission to build the package, run unit tests, measure code coverage,
     and check the coding style (PEP8) of your branch. The tests must pass
     before your PR can be merged. If CI fails, you can find out why by
     clicking on the "failed" icon (red cross) and inspecting the build and
     test logs.

   * A pull request must be approved by two core team members before merging.

5. Document changes

   If your change introduces any API modifications, please update
   ``doc/source/api_changes.txt``.

   If your change introduces a deprecation, add a reminder to ``TODO.txt``
   for the team to remove the deprecated functionality in the future.

.. note::

   To reviewers: if it is not obvious from the PR description, add a short
   explanation of what a branch did to the merge message and, if closing a
   bug, also add "Closes #123" where 123 is the issue number.


Divergence between ``upstream main`` and your feature branch
------------------------------------------------------------

If GitHub indicates that the branch of your PR can no longer
be merged automatically, merge the main branch into yours::

   git fetch upstream main
   git merge upstream/main

If any conflicts occur, they need to be fixed before continuing.  See
which files are in conflict using::

   git status

Which displays a message like::

   Unmerged paths:
     (use "git add <file>..." to mark resolution)

     both modified:   file_with_conflict.txt

Inside the conflicted file, you'll find sections like these::

   <<<<<<< HEAD
   The way the text looks in your branch
   =======
   The way the text looks in the main branch
   >>>>>>> main

Choose one version of the text that should be kept, and delete the
rest::

   The way the text looks in your branch

Now, add the fixed file::

   git add file_with_conflict.txt

Once you've fixed all merge conflicts, do::

   git commit

.. note::

   Advanced Git users are encouraged to `rebase instead of merge
   <https://scikit-image.org/docs/dev/gitwash/development_workflow.html#rebasing-on-trunk>`__,
   but we squash and merge most PRs either way.

Build environment setup
-----------------------

Please refer to :ref:`installing-scikit-image` for development installation
instructions.

Guidelines
----------

* All code should have tests (see `test coverage`_ below for more details).
* All code should be documented, to the same
  `standard <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_ as NumPy and SciPy.
* For new functionality, always add an example to the gallery (see
  :ref:`Sphinx-Gallery<sphinx_gallery>` below for more details).
* No changes are ever merged without review and approval by two core team members.
  There are two exceptions to this rule. First, pull requests which affect
  only the documentation require review and approval by only one core team
  member in most cases. If the maintainer feels the changes are large or
  likely to be controversial, two reviews should still be encouraged. The
  second case is that of minor fixes which restore CI to a working state,
  because these should be merged fairly quickly. Reach out on the
  `developer forum <https://discuss.scientific-python.org/c/contributor/skimage>`_ if
  you get no response to your pull request.
  **Never merge your own pull request.**

Stylistic Guidelines
--------------------

* Set up your editor to remove trailing whitespace.  Follow `PEP08
  <https://www.python.org/dev/peps/pep-0008/>`__.  Check code with pyflakes / flake8.

* Use numpy data types instead of strings (``np.uint8`` instead of
  ``"uint8"``).

* Use the following import conventions::

   import numpy as np
   import matplotlib.pyplot as plt
   from scipy import ndimage as ndi

   # only in Cython code
   cimport numpy as cnp
   cnp.import_array()

* When documenting array parameters, use ``image : (M, N) ndarray``
  and then refer to ``M`` and ``N`` in the docstring, if necessary.

* Refer to array dimensions as (plane), row, column, not as x, y, z. See
  :ref:`Coordinate conventions <numpy-images-coordinate-conventions>`
  in the user guide for more information.

* Functions should support all input image dtypes.  Use utility functions such
  as ``img_as_float`` to help convert to an appropriate type.  The output
  format can be whatever is most efficient.  This allows us to string together
  several functions into a pipeline, e.g.::

   hough(canny(my_image))

* Use ``Py_ssize_t`` as data type for all indexing, shape and size variables
  in C/C++ and Cython code.

* Use relative module imports, i.e. ``from .._shared import xyz`` rather than
  ``from skimage._shared import xyz``.

* Wrap Cython code in a pure Python function, which defines the API. This
  improves compatibility with code introspection tools, which are often not
  aware of Cython code.

* For Cython functions, release the GIL whenever possible, using
  ``with nogil:``.


Testing
-------

See the testing section of the Installation guide.

Test coverage
-------------

Tests for a module should ideally cover all code in that module,
i.e., statement coverage should be at 100%.

To measure the test coverage, install
`pytest-cov <https://pytest-cov.readthedocs.io/en/latest/>`__
(using ``pip install pytest-cov``) and then run::

  $ make coverage

This will print a report with one line for each file in `skimage`,
detailing the test coverage::

  Name                                             Stmts   Exec  Cover   Missing
  ------------------------------------------------------------------------------
  skimage/color/colorconv                             77     77   100%
  skimage/filter/__init__                              1      1   100%
  ...


Building docs
-------------

To build docs, run ``make`` from the ``doc`` directory. ``make help`` lists
all targets. For example, to build the HTML documentation, you can run:

.. code:: sh

    make html

Then, all the HTML files will be generated in ``scikit-image/doc/build/html/``.
To rebuild a full clean documentation, run:

.. code:: sh

    make clean
    make html

Requirements
~~~~~~~~~~~~

`Sphinx <http://www.sphinx-doc.org/en/stable/>`_,
`Sphinx-Gallery <https://sphinx-gallery.github.io>`_,
and LaTeX are needed to build the documentation.

**Sphinx:**

Sphinx and other python packages needed to build the documentation
can be installed using: ``scikit-image/requirements/docs.txt`` file.

.. code:: sh

    pip install -r requirements/docs.txt

.. _sphinx_gallery:

**Sphinx-Gallery:**

The above install command includes the installation of
`Sphinx-Gallery <https://sphinx-gallery.github.io>`_, which we use to create
the :ref:`examples_gallery`.
Refer to the Sphinx-Gallery documentation for complete instructions on syntax and usage.

If you are contributing an example to the gallery or editing an existing one,
build the docs (see above) and open a web browser to check how your edits
render at ``scikit-image/doc/build/html/auto_examples/``: navigate to the file
you have added or changed.

When adding an example, visit also
``scikit-image/doc/build/html/auto_examples/index.html`` to check how the new
thumbnail renders on the gallery's homepage. To change the thumbnail image,
please refer to `this section
<https://sphinx-gallery.github.io/stable/configuration.html#choosing-thumbnail>`_
of the Sphinx-Gallery docs.

Note that gallery examples should have a maximum figure width of 8 inches.

**LaTeX Ubuntu:**

.. code:: sh

    sudo apt-get install -qq texlive texlive-latex-extra dvipng

**LaTeX Mac:**

Install the full `MacTex <https://www.tug.org/mactex/>`__ installation or
install the smaller
`BasicTex <https://www.tug.org/mactex/morepackages.html>`__ and add *ucs*
and *dvipng* packages:

.. code:: sh

    sudo tlmgr install ucs dvipng

Fixing Warnings
~~~~~~~~~~~~~~~

-  "citation not found: R###" There is probably an underscore after a
   reference in the first line of a docstring (e.g. [1]\_). Use this
   method to find the source file: $ cd doc/build; grep -rin R####

-  "Duplicate citation R###, other instance in..."" There is probably a
   [2] without a [1] in one of the docstrings

-  Make sure to use pre-sphinxification paths to images (not the
   \_images directory)

Deprecation cycle
-----------------

If the behavior of the library has to be changed, a deprecation cycle must be
followed to warn users.

- a deprecation cycle is *not* necessary when:

    * adding a new function, or
    * adding a new keyword argument to the *end* of a function signature, or
    * fixing what was buggy behavior

- a deprecation cycle is necessary for *any breaking API change*, meaning a
    change where the function, invoked with the same arguments, would return a
    different result after the change. This includes:

    * changing the order of arguments or keyword arguments, or
    * adding arguments or keyword arguments to a function, or
    * changing a function's name or submodule, or
    * changing the default value of a function's arguments.

Usually, our policy is to put in place a deprecation cycle over two releases.

For the sake of illustration, we consider the modification of a default value in
a function signature. In version N (therefore, next release will be N+1), we
have

.. code-block:: python

    def a_function(image, rescale=True):
        out = do_something(image, rescale=rescale)
        return out

that has to be changed to

.. code-block:: python

    def a_function(image, rescale=None):
        if rescale is None:
            warn('The default value of rescale will change '
                 'to `False` in version N+3.', stacklevel=2)
            rescale = True
        out = do_something(image, rescale=rescale)
        return out

and in version N+3

.. code-block:: python

    def a_function(image, rescale=False):
        out = do_something(image, rescale=rescale)
        return out

Here is the process for a 2-release deprecation cycle:

- In the signature, set default to `None`, and modify the docstring to specify
  that it's `True`.
- In the function, _if_ rescale is set to `None`, set to `True` and warn that the
  default will change to `False` in version N+3.
- In ``doc/release/release_dev.rst``, under deprecations, add "In
  `a_function`, the `rescale` argument will default to `False` in N+3."
- In ``TODO.txt``, create an item in the section related to version N+3 and write
  "change rescale default to False in a_function".

Note that the 2-release deprecation cycle is not a strict rule and in some
cases, the developers can agree on a different procedure upon justification
(like when we can't detect the change, or it involves moving or deleting an
entire function for example).

Scikit-image uses warnings to highlight changes in its API so that users may
update their code accordingly. The ``stacklevel`` argument sets the location in
the callstack where the warnings will point. In most cases, it is appropriate
to set the ``stacklevel`` to ``2``.  When warnings originate from helper
routines internal to the scikit-image library, it is may be more appropriate to
set the ``stacklevel`` to ``3``. For more information, see the documentation of
the `warn <https://docs.python.org/3/library/warnings.html#warnings.warn>`__
function in the Python standard library.

To test if your warning is being emitted correctly, try calling the function
from an IPython console. It should point you to the console input itself
instead of being emitted by the files in the scikit-image library.

* **Good**: ``ipython:1: UserWarning: ...``
* **Bad**: ``scikit-image/skimage/measure/_structural_similarity.py:155: UserWarning:``

Bugs
----

Please `report bugs on GitHub <https://github.com/scikit-image/scikit-image/issues>`_.

Benchmarks
----------

While not mandatory for most pull requests, we ask that performance related
PRs include a benchmark in order to clearly depict the use-case that is being
optimized for. A historical view of our snapshots can be found on
at the following `website <https://pandas.pydata.org/speed/scikit-image/>`_.

In this section we will review how to setup the benchmarks,
and three commands ``asv dev``, ``asv run`` and ``asv continuous``.

Prerequisites
~~~~~~~~~~~~~
Begin by installing `airspeed velocity <https://asv.readthedocs.io/en/stable/>`_
in your development environment. Prior to installation, be sure to activate your
development environment, then if using ``venv`` you may install the requirement with::

  source skimage-dev/bin/activate
  pip install asv

If you are using conda, then the command::

  conda activate skimage-dev
  conda install asv

is more appropriate. Once installed, it is useful to run the command::

  asv machine

To let airspeed velocity know more information about your machine.

Writing a benchmark
~~~~~~~~~~~~~~~~~~~
To write  benchmark, add a file in the ``benchmarks`` directory which contains a
a class with one ``setup`` method and at least one method prefixed with ``time_``.

The ``time_`` method should only contain code you wish to benchmark.
Therefore it is useful to move everything that prepares the benchmark scenario
into the ``setup`` method. This function is called before calling a ``time_``
method and its execution time is not factored into the benchmarks.

Take for example the ``TransformSuite`` benchmark:

.. code-block:: python

  import numpy as np
  from skimage import transform

  class TransformSuite:
      """Benchmark for transform routines in scikit-image."""

      def setup(self):
          self.image = np.zeros((2000, 2000))
          idx = np.arange(500, 1500)
          self.image[idx[::-1], idx] = 255
          self.image[idx, idx] = 255

      def time_hough_line(self):
          result1, result2, result3 = transform.hough_line(self.image)

Here, the creation of the image is completed in the ``setup`` method, and not
included in the reported time of the benchmark.

It is also possible to benchmark features such as peak memory usage. To learn
more about the features of `asv`, please refer to the official
`airpseed velocity documentation <https://asv.readthedocs.io/en/latest/writing_benchmarks.html>`_.

Also, the benchmark files need to be importable when benchmarking old versions
of scikit-image. So if anything from scikit-image is imported at the top level,
it should be done as:

.. code-block:: python

    try:
        from skimage import metrics
    except ImportError:
        pass

The benchmarks themselves don't need any guarding against missing features,
only the top-level imports.

To allow tests of newer functions to be marked as "n/a" (not available)
rather than "failed" for older versions, the setup method itself can raise a
NotImplemented error.  See the following example for the registration module:

.. code-block:: python

    try:
        from skimage import registration
    except ImportError:
        raise NotImplementedError("registration module not available")

Testing the benchmarks locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prior to running the true benchmark, it is often worthwhile to test that the
code is free of typos. To do so, you may use the command::

  asv dev -b TransformSuite

Where the ``TransformSuite`` above will be run once in your current environment
to test that everything is in order.

Running your benchmark
~~~~~~~~~~~~~~~~~~~~~~

The command above is fast, but doesn't test the performance of the code
adequately. To do that you may want to run the benchmark in your current
environment to see the performance of your change as you are developing new
features. The command ``asv run -E existing`` will specify that you wish to run
the benchmark in your existing environment. This will save a significant amount
of time since building scikit-image can be a time consuming task::

  asv run -E existing -b TransformSuite

Comparing results to main
~~~~~~~~~~~~~~~~~~~~~~~~~

Often, the goal of a PR is to compare the results of the modifications in terms
speed to a snapshot of the code that is in the main branch of the
``scikit-image`` repository. The command ``asv continuous`` is of help here::

  asv continuous main -b TransformSuite

This call will build out the environments specified in the ``asv.conf.json``
file and compare the performance of the benchmark between your current commit
and the code in the main branch.

The output may look something like::

  $ asv continuous main -b TransformSuite
  · Creating environments
  · Discovering benchmarks
  ·· Uninstalling from conda-py3.7-cython-numpy1.15-scipy
  ·· Installing 544c0fe3 <benchmark_docs> into conda-py3.7-cython-numpy1.15-scipy.
  · Running 4 total benchmarks (2 commits * 2 environments * 1 benchmarks)
  [  0.00%] · For scikit-image commit 37c764cb <benchmark_docs~1> (round 1/2):
  [...]
  [100.00%] ··· ...ansform.TransformSuite.time_hough_line           33.2±2ms

  BENCHMARKS NOT SIGNIFICANTLY CHANGED.

In this case, the differences between HEAD and main are not significant
enough for airspeed velocity to report.

It is also possible to get a comparison of results for two specific revisions
for which benchmark results have previously been run via the `asv compare`
command::

    asv compare v0.14.5 v0.17.2

Finally, one can also run ASV benchmarks only for a specific commit hash or
release tag by appending ``^!`` to the commit or tag name. For example to run
the skimage.filter module benchmarks on release v0.17.2::

    asv run -b Filter v0.17.2^!
    .. _installing-scikit-image:

Installing scikit-image
==============================================================================

How you should install ``scikit-image`` depends on your needs and skills:

- Simplest solution:
  `scientific Python distribution <#scientific-python-distributions>`_.

- If you can install Python packages and work in virtual environments:

  - `pip <#install-via-pip>`_

  - `conda <#install-via-conda>`_

- Easy solution but with pitfalls: `system package manager <#system-package-managers>`_ (yum, apt, ...).

- `You're looking to contribute to scikit-image <#installing-scikit-image-for-contributors>`_.

Supported platforms
------------------------------------------------------------------------------

- Windows 64-bit on x86 processors
- Mac OS X on x86 processors
- Linux 64-bit on x86 processors

For information on other platforms, see `other platforms <#other-platforms>`_.

Version check
------------------------------------------------------------------------------

To see whether ``scikit-image`` is already installed or to check if an install has
worked, run the following in a Python shell or Jupyter notebook:

.. code-block:: python

  import skimage
  print(skimage.__version__)

or, from the command line:

.. code-block:: sh

   python -c "import skimage; print(skimage.__version__)"

(Try ``python3`` if ``python`` is unsuccessful.)

You'll see the version number if ``scikit-image`` is installed and
an error message otherwise.

Scientific Python distributions
------------------------------------------------------------------------------

In a single install these give you Python,
``scikit-image`` and libraries it depends on, and other useful scientific
packages. They install into an isolated environment, so they won't conflict
with any existing installed programs.

Drawbacks are that the install can be large and you may not get
the most recent ``scikit-image``.

We recommend one of these distributions:

- `Anaconda <https://www.anaconda.com/distribution/>`_
- `Python(x,y) <https://python-xy.github.io/>`_
- `WinPython <https://winpython.github.io/>`_

When using the ``scikit-image``
documentation, make sure it's for the version you've installed (see
`Version check <#version-check>`_ above).


Installation via pip and conda
------------------------------------------------------------------------------

These install only ``scikit-image`` and its dependencies; pip has an option to
include related packages.

.. _install-via-pip:

pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prerequisites to a pip install: You're able to use your system's command line to
install packages and are using a
`virtual environment
<https://towardsdatascience.com/virtual-environments-104c62d48c54?gi=2532aa12906#ee81>`_
(any of
`several
<https://stackoverflow.com/questions/41573587/what-is-the-difference-between-venv-pyvenv-pyenv-virtualenv-virtualenvwrappe>`_\
).

While it is possible to use pip without a virtual environment, it is not advised: 
virtual environments create a clean Python environment that does not interfere 
with any existing system installation, can be easily removed, and contain only
the package versions your application needs. They help avoid a common
challenge known as 
`dependency hell <https://en.wikipedia.org/wiki/Dependency_hell>`_.

To install the current ``scikit-image`` you'll need at least Python 3.6. If
your Python is older, pip will find the most recent compatible version.

.. code-block:: sh

  # Update pip
  python -m pip install -U pip
  # Install scikit-image
  python -m pip install -U scikit-image

To include a selection of other scientific Python packages that expand
``scikit-image``'s capabilities to include, e.g., parallel processing, you
can install the package ``scikit-image[optional]``:

.. code-block:: sh

    python -m pip install -U scikit-image[optional]

.. warning::

    Please do not use the command ``sudo`` and ``pip`` together as ``pip`` may
    overwrite critical system libraries which may require you to reinstall your
    operating system.

.. _install-via-conda:

conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Miniconda is a bare-essentials version of the Anaconda package; you'll need to
install packages like ``scikit-image`` yourself. Like Anaconda, it installs
Python and provides virtual environments.

- `conda documentation <https://docs.conda.io>`_
- `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
- `conda-forge <https://conda-forge.org>`_, a conda channel maintained
  with the latest ``scikit-image`` package

Once you have your conda environment set up, you can install ``scikit-image``
with the command:

.. code-block:: sh

    conda install scikit-image

System package managers
------------------------------------------------------------------------------

Using a package manager (``yum``, ``apt-get``, etc.) to install ``scikit-image``
or other Python packages is not your best option:

- You're likely to get an older version.

- You'll probably want to make updates and add new packages outside of
  the package manager, leaving you with the same kind of
  dependency conflicts you see when using pip without a virtual environment.

- There's an added risk because operating systems use Python, so if you
  make system-wide Python changes (installing as root or using sudo),
  you can break the operating system.


Downloading all demo datasets
------------------------------------------------------------------------------

Some of the data used in our examples is hosted online and is not installed
by default by the procedures explained above. Data are downloaded once, at the
first call, but this requires an internet connection. If you prefer downloading
all the demo datasets to be able to work offline, you can run this command:

.. code-block:: sh

    python -c 'from skimage.data import download_all; download_all()'

or call ``download_all()`` in your favourite interactive Python environment
(IPython, Jupyter notebook, ...).

Other platforms
------------------------------------------------------------------------------

We still support Windows 32-bit on x86 processors but urge switching
to Windows 64-bit.

Unsupported platforms include:

1. Linux on 32-bit x86 processors.
2. Linux on 32-bit on ARM processors (Raspberry Pi running Raspbian):

   - While we do not officially support this distribution, we point users to
     `piwheels <https://wwww.piwheels.org>`_
     and their
     `scikit-image's specific page <https://www.piwheels.org/project/scikit-image/>`_.

   - You may need to install additional system dependencies listed for
     `imagecodecs <https://www.piwheels.org/project/imagecodecs/>`_.
     See
     `issue 4721 <https://github.com/scikit-image/scikit-image/issues/4721>`_.

3. Linux on 64-bit ARM processors (Nvidia Jetson):

   - Follow the conversation on
     `issue 4705 <https://github.com/scikit-image/scikit-image/issues/4705>`_.

Although these platforms lack official support, many of the core
developers have experience with them and can help with questions.

If you want to install on an unsupported platform, try
`building from source <#building-from-source>`_.

Tell us which other platforms you'd like to see ``scikit-image`` on!
We are very interested in how ``scikit-image`` gets
`used <https://github.com/scikit-image/scikit-image/issues/4375>`_.

If you'd like to package ``scikit-image`` for an as-yet-unsupported platform,
`reach out on GitHub <https://github.com/scikit-image/scikit-image/issues>`_.


Additional help
------------------------------------------------------------------------------

If you still have questions, reach out through

- our `user forum <https://forum.image.sc/tags/scikit-image>`_
- our `developer forum <https://discuss.scientific-python.org/c/contributor/skimage>`_
- our `chat channel <https://skimage.zulipchat.com/>`_
- `Stack Overflow <https://stackoverflow.com/questions/tagged/scikit-image>`_


To suggest a change in these instructions,
`please open an issue on GitHub <https://github.com/scikit-image/scikit-image/issues/new>`_.


Installing scikit-image for contributors
========================================

We are assuming that you have a default Python environment already configured on
your computer and that you intend to install ``scikit-image`` inside of it.

We also make a few more assumptions about your system:

- You have a C compiler set up.
- You have a C++ compiler set up.
- You are running a version of Python compatible with our system as listed
  in our `setup.py file <https://github.com/scikit-image/scikit-image/blob/main/setup.py#L212>`_.
- You've cloned the git repository into a directory called ``scikit-image``.
  You have set up the `upstream` remote to point to our repository and `origin`
  to point to your fork.


This directory contains the following files:

.. code-block::

    scikit-image
    ├── asv.conf.json
    ├── azure-pipelines.yml
    ├── benchmarks
    ├── CODE_OF_CONDUCT.md
    ├── CONTRIBUTING.rst
    ├── CONTRIBUTORS.txt
    ├── doc
    ├── INSTALL.rst
    ├── LICENSE.txt
    ├── Makefile
    ├── MANIFEST.in
    ├── README.md
    ├── RELEASE.txt
    ├── requirements
    ├── requirements.txt
    ├── setup.cfg
    ├── setup.py
    ├── skimage
    ├── TODO.txt
    ├── tools

All commands below are assumed to be running from the ``scikit-image``
directory containing the files above.


Build environment setup
------------------------------------------------------------------------------

Once you've cloned your fork of the scikit-image repository,
you should set up a Python development environment tailored for scikit-image.
You may choose the environment manager of your choice.
Here we provide instructions for two popular environment managers:
``venv`` (pip based) and ``conda`` (Anaconda or Miniconda).

venv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using ``venv``, you may find the following bash commands useful:

.. code-block:: sh

  # Create a virtualenv named ``skimage-dev``
  python -m venv skimage-dev
  # Activate it. On Linux and MacOS:
  source skimage-dev/bin/activate
  # Make sure that pip is up to date
  pip install --upgrade pip
  # Install all development and runtime dependencies of scikit-image
  pip install -r <(cat requirements/*.txt)
  # Build and install scikit-image from source
  pip install -e . -vv
  # Test your installation
  pytest skimage

On Windows, please use ``skimage-dev\Scripts\activate`` on the activation step.

conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using conda for development, we
recommend adding the conda-forge channel for the most up-to-date version
of many dependencies.
Some dependencies we use (for testing and documentation) are not available
from the default Anaconda channel. Please follow the official
`conda-forge installation instructions <https://conda-forge.org/#about>`_
before you get started.

.. code-block:: sh

  # Create a conda environment named ``skimage-dev``
  conda create --name skimage-dev
  # Activate it
  conda activate skimage-dev
  # Install major development and runtime dependencies of scikit-image
  conda install `for i in requirements/{default,build,test}.txt; do echo -n " --file $i "; done`
  # Install scikit-image from source
  pip install -e . -vv
  # Test your installation
  pytest skimage

Updating the installation
------------------------------------------------------------------------------

When updating your installation, it is often necessary to recompile submodules
that have changed. Do so with the following commands:

.. code-block:: sh

    # Grab the latest source
    git checkout main
    git pull upstream main
    # Update the installation
    pip install -e . -vv

Testing
-------

``scikit-image`` has an extensive test suite that ensures correct
execution on your system.  The test suite must pass before a pull
request can be merged, and tests should be added to cover any
modifications to the code base.

We use the `pytest <https://docs.pytest.org/en/latest/>`__
testing framework, with tests located in the various
``skimage/submodule/tests`` folders.

Our testing requirements are listed below:

.. include:: ../../requirements/test.txt
   :literal:


Run all tests using:

.. code-block:: sh

    pytest skimage

Or the tests for a specific submodule:

.. code-block:: sh

    pytest skimage/morphology

Or tests from a specific file:

.. code-block:: sh

    pytest skimage/morphology/tests/test_grey.py

Or a single test within that file:

.. code-block:: sh

    pytest skimage/morphology/tests/test_grey.py::test_3d_fallback_black_tophat

Use ``--doctest-modules`` to run doctests. For example, run all tests and all
doctests using:

.. code-block:: sh

    pytest --doctest-modules skimage

Warnings during testing phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scikit-image tries to catch all warnings in its development builds to ensure
that crucial warnings from dependencies are not missed. This might cause
certain tests to fail if you are building scikit-image with versions of
dependencies that were not tested at the time of the release. To disable
failures on warnings, export the environment variable
``SKIMAGE_TEST_STRICT_WARNINGS`` with a value of `0` or `False` and run the
tests:

.. code-block:: sh

   export SKIMAGE_TEST_STRICT_WARNINGS=False
   pytest --pyargs skimage

Platform-specific notes
------------------------------------------------------------------------------

**Windows**

If you experience the error ``Error:unable to find vcvarsall.bat`` it means
that your computer does not have recommended compilers for Python. You can
either download and install Windows compilers from `here`_  or use
`MinGW compilers`_ . If using `MinGW`, make sure to correctly configure
``distutils`` by modifying (or create, if not existing) the configuration file
``distutils.cfg`` (located for example at
``C:\Python26\Lib\distutils\distutils.cfg``) to contain::

  [build]
   compiler=mingw32

A run-through of the compilation process for Windows is included in
our `setup of Azure Pipelines`_ (a continuous integration service).

.. _setup of Azure Pipelines: https://github.com/scikit-image/scikit-image/blob/main/azure-pipelines.yml
.. _here: https://wiki.python.org/moin/WindowsCompilers#Microsoft_Visual_C.2B-.2B-_14.0_standalone:_Visual_C.2B-.2B-_Build_Tools_2015_.28x86.2C_x64.2C_ARM.29
.. _MinGW compilers: http://www.mingw.org/wiki/howto_install_the_mingw_gcc_compiler_suite

**Debian and Ubuntu**

Install suitable compilers:

.. code-block:: sh

  sudo apt-get install build-essential


Full requirements list
----------------------
**Build Requirements**

.. include:: ../../requirements/build.txt
   :literal:

**Runtime Requirements**

.. include:: ../../requirements/default.txt
   :literal:

**Test Requirements**

.. include:: ../../requirements/test.txt
   :literal:

**Documentation Requirements**

.. include:: ../../requirements/docs.txt
   :literal:

**Optional Requirements**

You can use ``scikit-image`` with the basic requirements listed above, but some
functionality is only available with the following installed:

* `SimpleITK <http://www.simpleitk.org/>`__
    Optional I/O plugin providing a wide variety of `formats <https://itk.org/Wiki/ITK_File_Formats>`__.
    including specialized formats using in medical imaging.

* `Astropy <https://www.astropy.org>`__
    Provides FITS I/O capability.

* `PyAMG <https://pyamg.org/>`__
    The ``pyamg`` module is used for the fast ``cg_mg`` mode of random
    walker segmentation.

* `Dask <https://dask.org/>`__
    The ``dask`` module is used to speed up certain functions.


.. include:: ../../requirements/optional.txt
  :literal:


**Extra Requirements**

These requirements have been included as a convenience, but are not widely
installable through PyPI on our supported platforms. As such, we keep them in
a separate list for more advanced members of our community to install.

* `imread <https://pythonhosted.org/imread/>`__
    Optional I/O plugin providing most standard `formats <https://pythonhosted.org//imread/formats.html>`__.

.. include:: ../../requirements/extras.txt
  :literal:

Help with contributor installation
------------------------------------------------------------------------------

See `Additional help <#additional-help>`_ above.
.. |nbsp| unicode:: 0xA0
   :trim:

.. contents:: |nbsp|

Release notes
=============

.. include:: ../release/_release_notes_for_docs.rst
.. toctree:: 
   :hidden: 
 
   gitwash/index

.. include:: ../../CONTRIBUTING.rst
****************************
scikit-image |version| docs
****************************

scikit-image is an image processing toolbox for SciPy_. View the
`release notes <release_notes.html>`_ for |version|.

.. _SciPy: https://www.scipy.org


Sections
========

.. toctree::
   :hidden:

   overview
   api/api
   api_changes
   install
   user_guide
   glossary
   contribute
   license
   auto_examples/index
   release_notes
   core_developer
   skips/index
   values
   conduct/code_of_conduct

.. list-table::
   :class: contentstable

   * - `Overview <overview.html>`_

       Introduction to scikit-image.

     - `API Reference <api/api.html>`_ (`changes <api_changes.html>`_)

       Documentation for the functions included in scikit-image.

   * - `Mission Statement <values.html>`_

       Our mission, vision, and values.

     - `Governance <skips/1-governance.html>`_

       How decisions are made in scikit-image.

   * - `Installation <install.html>`_

       How to install scikit-image.

     - `User Guide <user_guide.html>`_

       Usage guidelines.

   * - `Glossary <glossary.html>`_

       Definitions of common terms.

     - `Contribute <contribute.html>`_

       Take part in development.  Core developers, please `read your guide
       <core_developer.html>`_.


   * - `Gallery <auto_examples/index.html>`_

       Introductory generic and domain-specific examples.

     - `License Info <license.html>`_

       Conditions on the use and redistribution of this package.

   * - `SKIPs <skips/index.html>`_

       **s**\ ci\ **k**\ it-\ **i**\ mage **p**\ roposals, documents describing major changes to the
       library.

     - `Code of Conduct <conduct/code_of_conduct.html>`_

       Community interaction guidelines.


Indices
=======

.. list-table::
   :class: contentstable

   * - :ref:`search`

       Search this documentation.

     - :ref:`genindex`

       All functions, classes, terms.
.. include:: ../../INSTALL.rst
:orphan:

.. _CoC_reporting_manual:

scikit-image Code of Conduct - How to follow up on a report
-----------------------------------------------------------

This is the manual followed by scikit-image's Code of Conduct Committee. It's used
when we respond to an issue to make sure we're consistent and fair.

Enforcing the Code of Conduct impacts our community today and for the future.
It's an action that we do not take lightly. When reviewing enforcement
measures, the Code of Conduct Committee will keep the following values and
guidelines in mind:

* Act in a personal manner rather than impersonal.  The Committee can engage
  the parties to understand the situation, while respecting the privacy and any
  necessary confidentiality of reporters.  However, sometimes it is necessary
  to communicate with one or more individuals directly: the Committee's goal is
  to improve the health of our community rather than only produce a formal
  decision.

* Emphasize empathy for individuals rather than judging behavior, avoiding
  binary labels of "good" and "bad/evil". Overt, clear-cut aggression and
  harassment exists and we will be address that firmly.  But many scenarios
  that can prove challenging to resolve are those where normal disagreements
  devolve into unhelpful or harmful behavior from multiple parties.
  Understanding the full context and finding a path that re-engages all is
  hard, but ultimately the most productive for our community.

* We understand that email is a difficult medium and can be isolating.
  Receiving criticism over email, without personal contact, can be
  particularly painful.  This makes it especially important to keep an
  atmosphere of open-minded respect of the views of others.  It also means
  that we must be transparent in our actions, and that we will do everything
  in our power to make sure that all our members are treated fairly and with
  sympathy.

* Discrimination can be subtle and it can be unconscious. It can show itself
  as unfairness and hostility in otherwise ordinary interactions.  We know
  that this does occur, and we will take care to look out for it.  We would
  very much like to hear from you if you feel you have been treated unfairly,
  and we will use these procedures to make sure that your complaint is heard
  and addressed.

* Help increase engagement in good discussion practice: try to identify where
  discussion may have broken down and provide actionable information, pointers
  and resources that can lead to positive change on these points.

* Be mindful of the needs of new members: provide them with explicit support
  and consideration, with the aim of increasing participation from
  underrepresented groups in particular.

* Individuals come from different cultural backgrounds and native languages.
  Try to identify any honest misunderstandings caused by a non-native speaker
  and help them understand the issue and what they can change to avoid causing
  offence.  Complex discussion in a foreign language can be very intimidating,
  and we want to grow our diversity also across nationalities and cultures.

*Mediation*: voluntary, informal mediation is a tool at our disposal.  In
contexts such as when two or more parties have all escalated to the point of
inappropriate behavior (something sadly common in human conflict), it may be
useful to facilitate a mediation process. This is only an example: the
Committee can consider mediation in any case, mindful that the process is meant
to be strictly voluntary and no party can be pressured to participate. If the
Committee suggests mediation, it should:

* Find a candidate who can serve as a mediator.
* Obtain the agreement of the reporter(s). The reporter(s) have complete
  freedom to decline the mediation idea, or to propose an alternate mediator.
* Obtain the agreement of the reported person(s).
* Settle on the mediator: while parties can propose a different mediator than
  the suggested candidate, only if common agreement is reached on all terms can
  the process move forward.
* Establish a timeline for mediation to complete, ideally within two weeks.

The mediator will engage with all the parties and seek a resolution that is
satisfactory to all.  Upon completion, the mediator will provide a report
(vetted by all parties to the process) to the Committee, with recommendations
on further steps.  The Committee will then evaluate these results (whether
satisfactory resolution was achieved or not) and decide on any additional
action deemed necessary.


How the committee will respond to reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the committee (or a committee member) receives a report, they will first
determine whether the report is about a clear and severe breach (as defined
below).  If so, immediate action needs to be taken in addition to the regular
report handling process.

Clear and severe breach actions
+++++++++++++++++++++++++++++++

We know that it is painfully common for internet communication to start at or
devolve into obvious and flagrant abuse.  We will deal quickly with clear and
severe breaches like personal threats, violent, sexist or racist language.

When a member of the Code of Conduct committee becomes aware of a clear and
severe breach, they will do the following:

* Immediately disconnect the originator from all scikit-image communication channels.
* Reply to the reporter that their report has been received and that the
  originator has been disconnected.
* In every case, the moderator should make a reasonable effort to contact the
  originator, and tell them specifically how their language or actions
  qualify as a "clear and severe breach".  The moderator should also say
  that, if the originator believes this is unfair or they want to be
  reconnected to scikit-image, they have the right to ask for a review, as below, by
  the Code of Conduct Committee.
  The moderator should copy this explanation to the Code of Conduct Committee.
* The Code of Conduct Committee will formally review and sign off on all cases
  where this mechanism has been applied to make sure it is not being used to
  control ordinary heated disagreement.

Report handling
+++++++++++++++

When a report is sent to the committee they will immediately reply to the
reporter to confirm receipt. This reply must be sent within 72 hours, and the
group should strive to respond much quicker than that.

If a report doesn't contain enough information, the committee will obtain all
relevant data before acting. The committee is empowered to act on the Steering
Council’s behalf in contacting any individuals involved to get a more complete
account of events.

The committee will then review the incident and determine, to the best of their
ability:

* What happened.
* Whether this event constitutes a Code of Conduct violation.
* Who are the responsible party(ies).
* Whether this is an ongoing situation, and there is a threat to anyone's
  physical safety.

This information will be collected in writing, and whenever possible the
group's deliberations will be recorded and retained (i.e. chat transcripts,
email discussions, recorded conference calls, summaries of voice conversations,
etc).

It is important to retain an archive of all activities of this committee to
ensure consistency in behavior and provide institutional memory for the
project.  To assist in this, the default channel of discussion for this
committee will be a private mailing list accessible to current and future
members of the committee as well as members of the Steering Council upon
justified request. If the Committee finds the need to use off-list
communications (e.g. phone calls for early/rapid response), it should in all
cases summarize these back to the list so there's a good record of the process.

The Code of Conduct Committee should aim to have a resolution agreed upon within
two weeks. In the event that a resolution can't be determined in that time, the
committee will respond to the reporter(s) with an update and projected timeline
for resolution.


.. _CoC_resolutions:

Resolutions
~~~~~~~~~~~

The committee must agree on a resolution by consensus. If the group cannot reach
consensus and deadlocks for over a week, the group will turn the matter over to
the Steering Council for resolution.


Possible responses may include:

* Taking no further action

  - if we determine no violations have occurred.
  - if the matter has been resolved publicly while the committee was considering responses.

* Coordinating voluntary mediation: if all involved parties agree, the
  Committee may facilitate a mediation process as detailed above.
* Remind publicly, and point out that some behavior/actions/language have been
  judged inappropriate and why in the current context, or can but hurtful to
  some people, requesting the community to self-adjust.
* A private reprimand from the committee to the individual(s) involved. In this
  case, the group chair will deliver that reprimand to the individual(s) over
  email, cc'ing the group.
* A public reprimand. In this case, the committee chair will deliver that
  reprimand in the same venue that the violation occurred, within the limits of
  practicality. E.g., for a discussion violation it will be on the
  originating forum, but for a chat room discussion where the
  person/context may be gone, another means may be better suited.
  The group may also publish the message elsewhere for documentation purposes.
* A request for a public or private apology, assuming the reporter agrees to
  this idea: they may at their discretion refuse further contact with the
  violator. The chair will deliver this request. The committee may, if it
  chooses, attach "strings" to this request: for example, the group may ask a
  violator to apologize in order to retain their membership on the forum.
* A "mutually agreed upon hiatus" where the committee asks the individual to
  temporarily refrain from community participation. If the individual chooses
  not to take a temporary break voluntarily, the committee may issue a
  "mandatory cooling off period".
* A permanent or temporary ban from some or all scikit-image spaces (forums,
  gitter.im, etc.). The group will maintain records of all such bans so that
  they may be reviewed in the future or otherwise maintained.

Once a resolution is agreed upon, but before it is enacted, the committee will
contact the original reporter and any other affected parties and explain the
proposed resolution. The committee will ask if this resolution is acceptable,
and must note feedback for the record.

Finally, the committee will make a report to the scikit-image Steering Council (as
well as the scikit-image core team in the event of an ongoing resolution, such as a
ban).

The committee will never publicly discuss the issue; all public statements will
be made by the chair of the Code of Conduct Committee or by the scikit-image Steering
Council.


Conflicts of Interest
~~~~~~~~~~~~~~~~~~~~~

In the event of any conflict of interest, a committee member must immediately
notify the other members, and recuse themselves if necessary.
.. _set-up-fork:

==================
 Set up your fork
==================

First you follow the instructions for :ref:`forking`.

Overview
========

::

   git clone git@github.com:your-user-name/scikit-image.git
   cd scikit-image
   git remote add upstream https://github.com/scikit-image/scikit-image.git

In detail
=========

Clone your fork
---------------

#. Clone your fork to the local computer with ``git clone
   git@github.com:your-user-name/scikit-image.git``
#. Investigate.  Change directory to your new repo: ``cd scikit-image``. Then
   ``git branch -a`` to show you all branches.  You'll get something
   like::

      * main
      remotes/origin/main

   This tells you that you are currently on the ``main`` branch, and
   that you also have a ``remote`` connection to ``origin/main``.
   What remote repository is ``remote/origin``? Try ``git remote -v`` to
   see the URLs for the remote.  They will point to your github fork.

   Now you want to connect to the upstream `scikit-image github`_ repository, so
   you can merge in changes from trunk.

.. _linking-to-upstream:

Linking your repository to the upstream repo
--------------------------------------------

::

   cd scikit-image
   git remote add upstream https://github.com/scikit-image/scikit-image.git

``upstream`` here is just the arbitrary name we're using to refer to the
main `scikit-image`_ repository at `scikit-image github`_.

Note that we've used ``https://`` for the URL rather than ``git@``.  The
``https://`` URL is read only.  This means we that we can't accidentally
(or deliberately) write to the upstream repo, and we are only going to
use it to merge into our own code.

Just for your own satisfaction, show yourself that you now have a new
'remote', with ``git remote -v show``, giving you something like::

   upstream	https://github.com/scikit-image/scikit-image.git (fetch)
   upstream	https://github.com/scikit-image/scikit-image.git (push)
   origin	git@github.com:your-user-name/scikit-image.git (fetch)
   origin	git@github.com:your-user-name/scikit-image.git (push)

.. include:: links.inc

.. _development-workflow:

####################
Development workflow
####################

You already have your own forked copy of the `scikit-image`_ repository, by
following :ref:`forking`. You have :ref:`set-up-fork`. You have configured
git by following :ref:`configure-git`.  Now you are ready for some real work.

Workflow summary
================

In what follows we'll refer to the upstream scikit-image ``main`` branch, as
"trunk".

* Don't use your ``main`` branch for anything.  Consider deleting it.
* When you are starting a new set of changes, fetch any changes from trunk,
  and start a new *feature branch* from that.
* Make a new branch for each separable set of changes |emdash| "one task, one
  branch" (`ipython git workflow`_).
* Name your branch for the purpose of the changes - e.g.
  ``bugfix-for-issue-14`` or ``refactor-database-code``.
* If you can possibly avoid it, avoid merging trunk or any other branches into
  your feature branch while you are working.
* If you do find yourself merging from trunk, consider :ref:`rebase-on-trunk`
* Ask on the `scikit-image developer forum`_ if you get stuck.
* Ask for code review!

This way of working helps to keep work well organized, with readable history.
This in turn makes it easier for project maintainers (that might be you) to see
what you've done, and why you did it.

See `linux git workflow`_ and `ipython git workflow`_ for some explanation.

Consider deleting your main branch
====================================

It may sound strange, but deleting your own ``main`` branch can help reduce
confusion about which branch you are on.  See `deleting master on github`_ for
details (replacing `master` with `main`).

.. _update-mirror-trunk:

Update the mirror of trunk
==========================

First make sure you have done :ref:`linking-to-upstream`.

From time to time you should fetch the upstream (trunk) changes from github::

   git fetch upstream

This will pull down any commits you don't have, and set the remote branches to
point to the right commit.  For example, 'trunk' is the branch referred to by
(remote/branchname) ``upstream/main`` - and if there have been commits since
you last checked, ``upstream/main`` will change after you do the fetch.

.. _make-feature-branch:

Make a new feature branch
=========================

When you are ready to make some changes to the code, you should start a new
branch.  Branches that are for a collection of related edits are often called
'feature branches'.

Making an new branch for each set of related changes will make it easier for
someone reviewing your branch to see what you are doing.

Choose an informative name for the branch to remind yourself and the rest of us
what the changes in the branch are for.  For example ``add-ability-to-fly``, or
``buxfix-for-issue-42``.

::

    # Update the mirror of trunk
    git fetch upstream
    # Make new feature branch starting at current trunk
    git branch my-new-feature upstream/main
    git checkout my-new-feature

Generally, you will want to keep your feature branches on your public github_
fork of `scikit-image`_.  To do this, you `git push`_ this new branch up to your
github repo.  Generally (if you followed the instructions in these pages, and by
default), git will have a link to your github repo, called ``origin``.  You push
up to your own repo on github with::

   git push origin my-new-feature

In git >= 1.7 you can ensure that the link is correctly set by using the
``--set-upstream`` option::

   git push --set-upstream origin my-new-feature

From now on git will know that ``my-new-feature`` is related to the
``my-new-feature`` branch in the github repo.

.. _edit-flow:

The editing workflow
====================

Overview
--------

::

   # hack hack
   git add my_new_file
   git commit -am 'NF - some message'
   git push

In more detail
--------------

#. Make some changes
#. See which files have changed with ``git status`` (see `git status`_).
   You'll see a listing like this one::

     # On branch ny-new-feature
     # Changed but not updated:
     #   (use "git add <file>..." to update what will be committed)
     #   (use "git checkout -- <file>..." to discard changes in working directory)
     #
     #	modified:   README
     #
     # Untracked files:
     #   (use "git add <file>..." to include in what will be committed)
     #
     #	INSTALL
     no changes added to commit (use "git add" and/or "git commit -a")

#. Check what the actual changes are with ``git diff`` (`git diff`_).
#. Add any new files to version control ``git add new_file_name`` (see
   `git add`_).
#. To commit all modified files into the local copy of your repo,, do
   ``git commit -am 'A commit message'``.  Note the ``-am`` options to
   ``commit``. The ``m`` flag just signals that you're going to type a
   message on the command line.  The ``a`` flag |emdash| you can just take on
   faith |emdash| or see `why the -a flag?`_ |emdash| and the helpful use-case
   description in the `tangled working copy problem`_. The `git commit`_ manual
   page might also be useful.
#. To push the changes up to your forked repo on github, do a ``git
   push`` (see `git push`_).

Ask for your changes to be reviewed or merged
=============================================

When you are ready to ask for someone to review your code and consider a merge:

#. Go to the URL of your forked repo, say
   ``https://github.com/your-user-name/scikit-image``.
#. Use the 'Switch Branches' dropdown menu near the top left of the page to
   select the branch with your changes:

   .. image:: branch_dropdown.png

#. Click on the 'Pull request' button:

   .. image:: pull_button.png

   Enter a title for the set of changes, and some explanation of what you've
   done.  Say if there is anything you'd like particular attention for - like a
   complicated change or some code you are not happy with.

   If you don't think your request is ready to be merged, just say so in your
   pull request message.  This is still a good way of getting some preliminary
   code review.

Some other things you might want to do
======================================

Delete a branch on github
-------------------------

::

   git checkout main
   # delete branch locally
   git branch -D my-unwanted-branch
   # delete branch on github
   git push origin :my-unwanted-branch

(Note the colon ``:`` before ``test-branch``.  See also:
https://help.github.com/en/github/using-git/managing-remote-repositories

Several people sharing a single repository
------------------------------------------

If you want to work on some stuff with other people, where you are all
committing into the same repository, or even the same branch, then just
share it via github.

First fork scikit-image into your account, as from :ref:`forking`.

Then, go to your forked repository github page, say
``https://github.com/your-user-name/scikit-image``

Click on the 'Admin' button, and add anyone else to the repo as a
collaborator:

   .. image:: pull_button.png

Now all those people can do::

    git clone git@githhub.com:your-user-name/scikit-image.git

Remember that links starting with ``git@`` use the ssh protocol.

Your collaborators can then commit directly into that repo with the
usual::

     git commit -am 'ENH - much better code'
     git push origin main # pushes directly into your repo

Explore your repository
-----------------------

To see a graphical representation of the repository branches and
commits::

   gitk --all

To see a linear list of commits for this branch::

   git log

You can also look at the `network graph visualizer`_ for your github
repo.

Finally the :ref:`fancy-log` ``lg`` alias will give you a reasonable text-based
graph of the repository.

.. _rebase-on-trunk:

Rebasing on trunk
-----------------

Let's say you thought of some work you'd like to do. You
:ref:`update-mirror-trunk` and :ref:`make-feature-branch` called
``cool-feature``. At this stage trunk is at some commit, let's call it E. Now
you make some new commits on your ``cool-feature`` branch, let's call them A, B,
C.  Maybe your changes take a while, or you come back to them after a while.  In
the meantime, trunk has progressed from commit E to commit (say) G::

          A---B---C cool-feature
         /
    D---E---F---G trunk

At this stage you consider merging trunk into your feature branch, and you
remember that this here page sternly advises you not to do that, because the
history will get messy. Most of the time you can just ask for a review, and not
worry that trunk has got a little ahead.  But sometimes, the changes in trunk
might affect your changes, and you need to harmonize them.  In this situation
you may prefer to do a rebase.

rebase takes your changes (A, B, C) and replays them as if they had been made to
the current state of ``trunk``.  In other words, in this case, it takes the
changes represented by A, B, C and replays them on top of G. After the rebase,
your history will look like this::

                  A'--B'--C' cool-feature
                 /
    D---E---F---G trunk

See `rebase without tears`_ for more detail.

To do a rebase on trunk::

    # Update the mirror of trunk
    git fetch upstream
    # go to the feature branch
    git checkout cool-feature
    # make a backup in case you mess up
    git branch tmp cool-feature
    # rebase cool-feature onto trunk
    git rebase --onto upstream/main upstream/main cool-feature

In this situation, where you are already on branch ``cool-feature``, the last
command can be written more succinctly as::

    git rebase upstream/main

When all looks good you can delete your backup branch::

   git branch -D tmp

If it doesn't look good you may need to have a look at
:ref:`recovering-from-mess-up`.

If you have made changes to files that have also changed in trunk, this may
generate merge conflicts that you need to resolve - see the `git rebase`_ man
page for some instructions at the end of the "Description" section. There is
some related help on merging in the git user manual - see `resolving a merge`_.

.. _recovering-from-mess-up:

Recovering from mess-ups
------------------------

Sometimes, you mess up merges or rebases. Luckily, in git it is
relatively straightforward to recover from such mistakes.

If you mess up during a rebase::

   git rebase --abort

If you notice you messed up after the rebase::

   # reset branch back to the saved point
   git reset --hard tmp

If you forgot to make a backup branch::

   # look at the reflog of the branch
   git reflog show cool-feature

   8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
   278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
   26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj
   ...

   # reset the branch to where it was before the botched rebase
   git reset --hard cool-feature@{2}

.. _rewriting-commit-history:

Rewriting commit history
------------------------

.. note::

   Do this only for your own feature branches.

There's an embarrassing typo in a commit you made? Or perhaps the you
made several false starts you would like the posterity not to see.

This can be done via *interactive rebasing*.

Suppose that the commit history looks like this::

    git log --oneline
    eadc391 Fix some remaining bugs
    a815645 Modify it so that it works
    2dec1ac Fix a few bugs + disable
    13d7934 First implementation
    6ad92e5 * masked is now an instance of a new object, MaskedConstant
    29001ed Add pre-nep for a copule of structured_array_extensions.
    ...

and ``6ad92e5`` is the last commit in the ``cool-feature`` branch. Suppose we
want to make the following changes:

* Rewrite the commit message for ``13d7934`` to something more sensible.
* Combine the commits ``2dec1ac``, ``a815645``, ``eadc391`` into a single one.

We do as follows::

    # make a backup of the current state
    git branch tmp HEAD
    # interactive rebase
    git rebase -i 6ad92e5

This will open an editor with the following text in it::

    pick 13d7934 First implementation
    pick 2dec1ac Fix a few bugs + disable
    pick a815645 Modify it so that it works
    pick eadc391 Fix some remaining bugs

    # Rebase 6ad92e5..eadc391 onto 6ad92e5
    #
    # Commands:
    #  p, pick = use commit
    #  r, reword = use commit, but edit the commit message
    #  e, edit = use commit, but stop for amending
    #  s, squash = use commit, but meld into previous commit
    #  f, fixup = like "squash", but discard this commit's log message
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    # However, if you remove everything, the rebase will be aborted.
    #

To achieve what we want, we will make the following changes to it::

    r 13d7934 First implementation
    pick 2dec1ac Fix a few bugs + disable
    f a815645 Modify it so that it works
    f eadc391 Fix some remaining bugs

This means that (i) we want to edit the commit message for
``13d7934``, and (ii) collapse the last three commits into one. Now we
save and quit the editor.

Git will then immediately bring up an editor for editing the commit
message. After revising it, we get the output::

    [detached HEAD 721fc64] FOO: First implementation
     2 files changed, 199 insertions(+), 66 deletions(-)
    [detached HEAD 0f22701] Fix a few bugs + disable
     1 files changed, 79 insertions(+), 61 deletions(-)
    Successfully rebased and updated refs/heads/my-feature-branch.

and the history looks now like this::

     0f22701 Fix a few bugs + disable
     721fc64 ENH: Sophisticated feature
     6ad92e5 * masked is now an instance of a new object, MaskedConstant

If it went wrong, recovery is again possible as explained :ref:`above
<recovering-from-mess-up>`.

.. include:: links.inc
.. _forking:

======================================================
Making your own copy (fork) of scikit-image
======================================================

You need to do this only once.  The instructions here are very similar
to the instructions at https://help.github.com/en/github/getting-started-with-github/fork-a-repo |emdash|
please see that page for more detail.  We're repeating some of it here just to give
the specifics for the `scikit-image`_ project, and to suggest some default names.

Set up and configure a github account
=====================================

If you don't have a github account, go to the github page, and make one.

You then need to configure your account to allow write access |emdash| see
the ``Generating SSH keys`` help on `github help`_.

Create your own forked copy of `scikit-image`_
======================================================

#. Log into your github account.
#. Go to the `scikit-image`_ github home at `scikit-image github`_.
#. Click on the *fork* button:

   .. image:: forking_button.png

   Now, after a short pause and some 'Hardcore forking action', you
   should find yourself at the home page for your own forked copy of `scikit-image`_.

.. include:: links.inc

.. _configure-git:

===============
 Configure git
===============

.. _git-config-basic:

Overview
========

Your personal git configurations are saved in the ``.gitconfig`` file in
your home directory.

Here is an example ``.gitconfig`` file::

  [user]
          name = Your Name
          email = you@yourdomain.example.com

  [alias]
          ci = commit -a
          co = checkout
          st = status
          stat = status
          br = branch
          wdiff = diff --color-words

  [core]
          editor = vim

  [merge]
          summary = true

You can edit this file directly or you can use the ``git config --global``
command::

  git config --global user.name "Your Name"
  git config --global user.email you@yourdomain.example.com
  git config --global alias.ci "commit -a"
  git config --global alias.co checkout
  git config --global alias.st "status -a"
  git config --global alias.stat "status -a"
  git config --global alias.br branch
  git config --global alias.wdiff "diff --color-words"
  git config --global core.editor vim
  git config --global merge.summary true

To set up on another computer, you can copy your ``~/.gitconfig`` file,
or run the commands above.

In detail
=========

user.name and user.email
------------------------

It is good practice to tell git_ who you are, for labeling any changes
you make to the code.  The simplest way to do this is from the command
line::

  git config --global user.name "Your Name"
  git config --global user.email you@yourdomain.example.com

This will write the settings into your git configuration file,  which
should now contain a user section with your name and email::

  [user]
        name = Your Name
        email = you@yourdomain.example.com

Of course you'll need to replace ``Your Name`` and ``you@yourdomain.example.com``
with your actual name and email address.

Aliases
-------

You might well benefit from some aliases to common commands.

For example, you might well want to be able to shorten ``git checkout``
to ``git co``.  Or you may want to alias ``git diff --color-words``
(which gives a nicely formatted output of the diff) to ``git wdiff``

The following ``git config --global`` commands::

  git config --global alias.ci "commit -a"
  git config --global alias.co checkout
  git config --global alias.st "status -a"
  git config --global alias.stat "status -a"
  git config --global alias.br branch
  git config --global alias.wdiff "diff --color-words"

will create an ``alias`` section in your ``.gitconfig`` file with contents
like this::

  [alias]
          ci = commit -a
          co = checkout
          st = status -a
          stat = status -a
          br = branch
          wdiff = diff --color-words

Editor
------

You may also want to make sure that your editor of choice is used ::

  git config --global core.editor vim

Merging
-------

To enforce summaries when doing merges (``~/.gitconfig`` file again)::

   [merge]
      log = true

Or from the command line::

  git config --global merge.log true

.. _fancy-log:

Fancy log output
----------------

This is a very nice alias to get a fancy log output; it should go in the
``alias`` section of your ``.gitconfig`` file::

    lg = log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)[%an]%Creset' --abbrev-commit --date=relative

You use the alias with::

    git lg

and it gives graph / text output something like this (but with color!)::

    * 6d8e1ee - (HEAD, origin/my-fancy-feature, my-fancy-feature) NF - a fancy file (45 minutes ago) [Matthew Brett]
    *   d304a73 - (origin/placeholder, placeholder) Merge pull request #48 from hhuuggoo/master (2 weeks ago) [Jonathan Terhorst]
    |\  
    | * 4aff2a8 - fixed bug 35, and added a test in test_bugfixes (2 weeks ago) [Hugo]
    |/  
    * a7ff2e5 - Added notes on discussion/proposal made during Data Array Summit. (2 weeks ago) [Corran Webster]
    * 68f6752 - Initial implimentation of AxisIndexer - uses 'index_by' which needs to be changed to a call on an Axes object - this is all very sketchy right now. (2 weeks ago) [Corr
    *   376adbd - Merge pull request #46 from terhorst/master (2 weeks ago) [Jonathan Terhorst]
    |\  
    | * b605216 - updated joshu example to current api (3 weeks ago) [Jonathan Terhorst]
    | * 2e991e8 - add testing for outer ufunc (3 weeks ago) [Jonathan Terhorst]
    | * 7beda5a - prevent axis from throwing an exception if testing equality with non-axis object (3 weeks ago) [Jonathan Terhorst]
    | * 65af65e - convert unit testing code to assertions (3 weeks ago) [Jonathan Terhorst]
    | *   956fbab - Merge remote-tracking branch 'upstream/master' (3 weeks ago) [Jonathan Terhorst]
    | |\  
    | |/

Thanks to Yury V. Zaytsev for posting it.

.. include:: links.inc
.. _install-git:

=============
 Install git
=============

Overview
========

================ =============
Debian / Ubuntu  ``sudo apt-get install git``
Fedora           ``sudo yum install git-core``
Windows          Download and install msysGit_
OS X             Use the git-osx-installer_
================ =============

In detail
=========

See the git page for the most recent information.

Have a look at the github install help pages available from `github help`_

There are good instructions here:
https://book.git-scm.com/book/en/v2/Getting-Started-Installing-Git

.. include:: links.inc
.. _maintainer-workflow:

###################
Maintainer workflow
###################

This page is for maintainers |emdash| those of us who merge our own or other
peoples' changes into the upstream repository.

Being as how you're a maintainer, you are completely on top of the basic stuff
in :ref:`development-workflow`.

The instructions in :ref:`linking-to-upstream` add a remote that has read-only
access to the upstream repo.  Being a maintainer, you've got read-write access.

It's good to have your upstream remote have a scary name, to remind you that
it's a read-write remote::

    git remote add upstream-rw git@github.com:scikit-image/scikit-image.git
    git fetch upstream-rw

*******************
Integrating changes
*******************

Let's say you have some changes that need to go into trunk
(``upstream-rw/main``).

The changes are in some branch that you are currently on.  For example, you are
looking at someone's changes like this::

    git remote add someone https://github.com/someone/scikit-image.git
    git fetch someone
    git branch cool-feature --track someone/cool-feature
    git checkout cool-feature

So now you are on the branch with the changes to be incorporated upstream.  The
rest of this section assumes you are on this branch.

A few commits
=============

If there are only a few commits, consider rebasing to upstream::

    # Fetch upstream changes
    git fetch upstream-rw
    # rebase
    git rebase upstream-rw/main

Remember that, if you do a rebase, and push that, you'll have to close any
github pull requests manually, because github will not be able to detect the
changes have already been merged.

A long series of commits
========================

If there are a longer series of related commits, consider a merge instead::

    git fetch upstream-rw
    git merge --no-ff upstream-rw/main

The merge will be detected by github, and should close any related pull requests
automatically.

Note the ``--no-ff`` above.  This forces git to make a merge commit, rather than
doing a fast-forward, so that these set of commits branch off trunk then rejoin
the main history with a merge, rather than appearing to have been made directly
on top of trunk.

Check the history
=================

Now, in either case, you should check that the history is sensible and you have
the right commits::

    git log --oneline --graph
    git log -p upstream-rw/main..

The first line above just shows the history in a compact way, with a text
representation of the history graph. The second line shows the log of commits
excluding those that can be reached from trunk (``upstream-rw/main``), and
including those that can be reached from current HEAD (implied with the ``..``
at the end). So, it shows the commits unique to this branch compared to trunk.
The ``-p`` option shows the diff for these commits in patch form.

Push to trunk
=============

::

    git push upstream-rw my-new-feature:main

This pushes the ``my-new-feature`` branch in this repository to the ``main``
branch in the ``upstream-rw`` repository.

.. include:: links.inc
.. _git-development:

=====================
 Git for development
=====================

Contents:

.. toctree::
   :maxdepth: 2

   forking_hell
   set_up_fork
   configure_git
   development_workflow
   maintainer_workflow
.. _using-git:

Working with *scikit-image* source code
================================================

Contents:

.. toctree::
   :maxdepth: 2

   git_intro
   git_install
   following_latest
   patching
   git_development
   git_resources


================
 Making a patch
================

You've discovered a bug or something else you want to change
in `scikit-image`_ .. |emdash| excellent!

You've worked out a way to fix it |emdash| even better!

You want to tell us about it |emdash| best of all!

The easiest way is to make a *patch* or set of patches.  Here
we explain how.  Making a patch is the simplest and quickest,
but if you're going to be doing anything more than simple
quick things, please consider following the
:ref:`git-development` model instead.

.. _making-patches:

Making patches
==============

Overview
--------

::

   # tell git who you are
   git config --global user.email you@yourdomain.example.com
   git config --global user.name "Your Name Comes Here"
   # get the repository if you don't have it
   git clone https://github.com/scikit-image/scikit-image.git
   # make a branch for your patching
   cd scikit-image
   git branch the-fix-im-thinking-of
   git checkout the-fix-im-thinking-of
   # hack, hack, hack
   # Tell git about any new files you've made
   git add somewhere/tests/test_my_bug.py
   # commit work in progress as you go
   git commit -am 'BF - added tests for Funny bug'
   # hack hack, hack
   git commit -am 'BF - added fix for Funny bug'
   # make the patch files
   git format-patch -M -C main

Then, send the generated patch files to the `scikit-image
developer forum`_ |emdash| where we will thank you warmly.

In detail
---------

#. Tell git who you are so it can label the commits you've
   made::

      git config --global user.email you@yourdomain.example.com
      git config --global user.name "Your Name Comes Here"

#. If you don't already have one, clone a copy of the
   `scikit-image`_ repository::

      git clone https://github.com/scikit-image/scikit-image.git
      cd scikit-image

#. Make a 'feature branch'.  This will be where you work on
   your bug fix.  It's nice and safe and leaves you with
   access to an unmodified copy of the code in the main
   branch::

      git branch the-fix-im-thinking-of
      git checkout the-fix-im-thinking-of

#. Do some edits, and commit them as you go::

      # hack, hack, hack
      # Tell git about any new files you've made
      git add somewhere/tests/test_my_bug.py
      # commit work in progress as you go
      git commit -am 'BF - added tests for Funny bug'
      # hack hack, hack
      git commit -am 'BF - added fix for Funny bug'

   Note the ``-am`` options to ``commit``. The ``m`` flag just
   signals that you're going to type a message on the command
   line.  The ``a`` flag |emdash| you can just take on faith |emdash|
   or see `why the -a flag?`_.

#. When you have finished, check you have committed all your
   changes::

      git status

#. Finally, make your commits into patches.  You want all the
   commits since you branched from the ``main`` branch::

      git format-patch -M -C main

   You will now have several files named for the commits::

      0001-BF-added-tests-for-Funny-bug.patch
      0002-BF-added-fix-for-Funny-bug.patch

   Send these files to the `scikit-image developer forum`_.

When you are done, to switch back to the main copy of the
code, just return to the ``main`` branch::

   git checkout main

Moving from patching to development
===================================

If you find you have done some patches, and you have one or
more feature branches, you will probably want to switch to
development mode.  You can do this with the repository you
have.

Fork the `scikit-image`_ repository on github |emdash| :ref:`forking`.
Then::

   # checkout and refresh main branch from main repo
   git checkout main
   git pull origin main
   # rename pointer to main repository to 'upstream'
   git remote rename origin upstream
   # point your repo to default read / write to your fork on github
   git remote add origin git@github.com:your-user-name/scikit-image.git
   # push up any branches you've made and want to keep
   git push origin the-fix-im-thinking-of

Then you can, if you want, follow the
:ref:`development-workflow`.

.. include:: links.inc
.. _git-resources:

=============
git resources
=============

Tutorials and summaries
=======================

* `github help`_ has an excellent series of how-to guides.
* `learn.github`_ has an excellent series of tutorials
* The `pro git book`_ is a good in-depth book on git. 
* A `git cheat sheet`_ is a page giving summaries of common commands.
* The `git user manual`_ 
* The `git tutorial`_
* The `git community book`_
* `git ready`_ |emdash| a nice series of tutorials
* `git casts`_ |emdash| video snippets giving git how-tos.
* `git magic`_ |emdash| extended introduction with intermediate detail
* The `git parable`_ is an easy read explaining the concepts behind git.
* `git foundation`_ expands on the `git parable`_.
* Fernando Perez' git page |emdash| `Fernando's git page`_ |emdash| many
  links and tips
* A good but technical page on `git concepts`_
* `git svn crash course`_: git for those of us used to subversion_

Advanced git workflow
=====================

There are many ways of working with git; here is a posts on the
rules of thumb that other projects have come up with:

* Linus Torvalds on `linux git workflow`_ .  Summary; use the git tools
  to make the history of your edits as clean as possible; merge from
  upstream edits as little as possible in branches where you are doing
  active development.

Manual pages online
===================

You can get these on your own machine with (e.g) ``git help push`` or
(same thing) ``git push --help``, but, for convenience, here are the
online manual pages for some common commands:

* `git add`_
* `git branch`_
* `git checkout`_
* `git clone`_
* `git commit`_
* `git config`_
* `git diff`_
* `git log`_
* `git pull`_
* `git push`_
* `git remote`_
* `git status`_

.. include:: links.inc
.. _following-latest:

=============================
 Following the latest source
=============================

These are the instructions if you just want to follow the latest
*scikit-image* source, but you don't need to do any development for now.

The steps are:

* :ref:`install-git`
* get local copy of the `scikit-image github`_ git repository
* update local copy from time to time

Get the local copy of the code
==============================

From the command line::

   git clone https://github.com/scikit-image/scikit-image.git

You now have a copy of the code tree in the new ``scikit-image`` directory.

Updating the code
=================

From time to time you may want to pull down the latest code.  Do this with::

   cd scikit-image
   git pull

The tree in ``scikit-image`` will now have the latest changes from the initial
repository.

.. include:: links.inc
==============
 Introduction
==============

These pages describe a git_ and github_ workflow for the `scikit-image`_
project.

There are several different workflows here, for different ways of
working with *scikit-image*.

This is not a comprehensive git reference, it's just a workflow for our
own project.  It's tailored to the github hosting service. You may well
find better or quicker ways of getting stuff done with git, but these
should get you started.

For general resources for learning git, see :ref:`git-resources`.

.. include:: links.inc
.. _numpy:

==================================
A crash course on NumPy for images
==================================

Images in ``scikit-image`` are represented by NumPy ndarrays. Hence, many 
common operations can be achieved using standard NumPy methods for 
manipulating arrays::

    >>> from skimage import data
    >>> camera = data.camera()
    >>> type(camera)
    <type 'numpy.ndarray'>

Retrieving the geometry of the image and the number of pixels::

    >>> camera.shape
    (512, 512)
    >>> camera.size
    262144

Retrieving statistical information about image intensity values::

    >>> camera.min(), camera.max()
    (0, 255)
    >>> camera.mean()
    118.31400299072266

NumPy arrays representing images can be of different integer or float
numerical types. See :ref:`data_types` for more information about these
types and how ``scikit-image`` treats them.


NumPy indexing
--------------

NumPy indexing can be used both for looking at the pixel values and to
modify them::

    >>> # Get the value of the pixel at the 10th row and 20th column
    >>> camera[10, 20]
    153
    >>> # Set to black the pixel at the 3rd row and 10th column
    >>> camera[3, 10] = 0

Be careful! In NumPy indexing, the first dimension (``camera.shape[0]``)
corresponds to rows, while the second (``camera.shape[1]``) corresponds
to columns, with the origin (``camera[0, 0]``) at the top-left corner.
This matches matrix/linear algebra notation, but is in contrast to
Cartesian (x, y) coordinates. See `Coordinate conventions`_ below for
more details.

Beyond individual pixels, it is possible to access/modify values of
whole sets of pixels using the different indexing capabilities of NumPy.

Slicing::

    >>> # Set the first ten lines to "black" (0)
    >>> camera[:10] = 0

Masking (indexing with masks of booleans)::

    >>> mask = camera < 87
    >>> # Set to "white" (255) the pixels where mask is True
    >>> camera[mask] = 255

Fancy indexing (indexing with sets of indices)::

    >>> inds_r = np.arange(len(camera))
    >>> inds_c = 4 * inds_r % len(camera)
    >>> camera[inds_r, inds_c] = 0

Masks are very useful when you need to select a set of pixels on which
to perform the manipulations. The mask can be any boolean array
of the same shape as the image (or a shape broadcastable to the image shape).
This can be used to define a region of interest, for example, a disk::

    >>> nrows, ncols = camera.shape
    >>> row, col = np.ogrid[:nrows, :ncols]
    >>> cnt_row, cnt_col = nrows / 2, ncols / 2
    >>> outer_disk_mask = ((row - cnt_row)**2 + (col - cnt_col)**2 >
    ...                    (nrows / 2)**2)
    >>> camera[outer_disk_mask] = 0

.. image:: ../auto_examples/numpy_operations/images/sphx_glr_plot_camera_numpy_001.png
    :width: 45%
    :target: ../auto_examples/numpy_operations/plot_camera_numpy.html

Boolean operations from NumPy can be used to define even more complex masks::

    >>> lower_half = row > cnt_row
    >>> lower_half_disk = np.logical_and(lower_half, outer_disk_mask)
    >>> camera = data.camera()
    >>> camera[lower_half_disk] = 0


Color images
------------

All of the above remains true for color images. A color image is a
NumPy array with an additional trailing dimension for the channels::

    >>> cat = data.chelsea()
    >>> type(cat)
    <type 'numpy.ndarray'>
    >>> cat.shape
    (300, 451, 3)

This shows that ``cat`` is a 300-by-451 pixel image with three channels
(red, green, and blue). As before, we can get and set the pixel values::

    >>> cat[10, 20]
    array([151, 129, 115], dtype=uint8)
    >>> # Set the pixel at (50th row, 60th column) to "black"
    >>> cat[50, 60] = 0
    >>> # set the pixel at (50th row, 61st column) to "green"
    >>> cat[50, 61] = [0, 255, 0]  # [red, green, blue]

We can also use 2D boolean masks for 2D multichannel images, as we did with
the grayscale image above:

.. plot::

    Using a 2D mask on a 2D color image

    >>> from skimage import data
    >>> cat = data.chelsea()
    >>> reddish = cat[:, :, 0] > 160
    >>> cat[reddish] = [0, 255, 0]
    >>> plt.imshow(cat)

The example color images included in ``skimage.data`` have channels stored
along the last axis, although other software may follow different conventions.
The scikit-image library functions supporting color images have a
``channel_axis`` argument that can be used to specify which axis of an array
corresponds to channels.

.. _numpy-images-coordinate-conventions:

Coordinate conventions
----------------------

Because ``scikit-image`` represents images using NumPy arrays, the
coordinate conventions must match. Two-dimensional (2D) grayscale images
(such as `camera` above) are indexed by rows and columns (abbreviated to
either ``(row, col)`` or ``(r, c)``), with the lowest element ``(0, 0)``
at the top-left corner. In various parts of the library, you will
also see ``rr`` and ``cc`` refer to lists of row and column
coordinates. We distinguish this convention from ``(x, y)``, which commonly
denote standard Cartesian coordinates, where ``x`` is the horizontal coordinate,
``y`` - the vertical one, and the origin is at the bottom left
(Matplotlib axes, for example, use this convention).

In the case of multichannel images, any dimension (array axis) can be used for
color channels, and is denoted by ``channel`` or ``ch``. Prior to scikit-image
0.19, this channel dimension was always last, but in the current release the
channel dimension can be specified by a ``channel_axis`` argument. Functions
that require multichannel data default to ``channel_axis=-1``. Otherwise,
functions default to ``channel_axis=None``, indicating that no axis is
assumed to correspond to channels.

Finally, for volumetric (3D) images, such as videos, magnetic resonance imaging
(MRI) scans, confocal microscopy, etc., we refer to the leading dimension
as ``plane``, abbreviated as ``pln`` or ``p``.

These conventions are summarized below:

.. table:: *Dimension name and order conventions in scikit-image*

  =========================   =============================
  Image type                  Coordinates
  =========================   =============================
  2D grayscale                (row, col)
  2D multichannel (eg. RGB)   (row, col, ch)
  3D grayscale                (pln, row, col)
  3D multichannel             (pln, row, col, ch)
  =========================   =============================

Note that the position of ``ch`` is controlled by the ``channel_axis``
argument.

|

Many functions in ``scikit-image`` can operate on 3D images directly::

    >>> rng = np.random.default_rng()
    >>> im3d = rng.random((100, 1000, 1000))
    >>> from skimage import morphology
    >>> from scipy import ndimage as ndi
    >>> seeds = ndi.label(im3d < 0.1)[0]
    >>> ws = morphology.watershed(im3d, seeds)

In many cases, however, the third spatial dimension has lower resolution
than the other two. Some ``scikit-image`` functions provide a ``spacing``
keyword argument to help handle this kind of data::

    >>> from skimage import segmentation
    >>> slics = segmentation.slic(im3d, spacing=[5, 1, 1], channel_axis=None)

Other times, the processing must be done plane-wise. When planes are stacked
along the leading dimension (in agreement with our convention), the following
syntax can be used::

    >>> from skimage import filters
    >>> edges = np.empty_like(im3d)
    >>> for pln, image in enumerate(im3d):
    ...     # Iterate over the leading dimension 
    ...     edges[pln] = filters.sobel(image)


Notes on the order of array dimensions
--------------------------------------

Although the labeling of the axes might seem arbitrary, it can have a
significant effect on the speed of operations. This is because modern
processors never retrieve just one item from memory, but rather a whole
chunk of adjacent items (an operation called prefetching). Therefore,
processing of elements that are next to each other in memory is faster
than processing them when they are scattered, even if the number of operations
is the same::

    >>> def in_order_multiply(arr, scalar):
    ...     for plane in list(range(arr.shape[0])):
    ...         arr[plane, :, :] *= scalar
    ...
    >>> def out_of_order_multiply(arr, scalar):
    ...     for plane in list(range(arr.shape[2])):
    ...         arr[:, :, plane] *= scalar
    ...
    >>> import time
    >>> rng = np.random.default_rng()
    >>> im3d = rng.random((100, 1024, 1024))
    >>> t0 = time.time(); x = in_order_multiply(im3d, 5); t1 = time.time()
    >>> print("%.2f seconds" % (t1 - t0))  # doctest: +SKIP
    0.14 seconds
    >>> s0 = time.time(); x = out_of_order_multiply(im3d, 5); s1 = time.time()
    >>> print("%.2f seconds" % (s1 - s0))  # doctest: +SKIP
    1.18 seconds
    >>> print("Speedup: %.1fx" % ((s1 - s0) / (t1 - t0)))  # doctest: +SKIP
    Speedup: 8.6x


When the last/rightmost dimension becomes even larger the speedup is
even more dramatic. It is worth thinking about *data locality* when
developing algorithms. In particular, ``scikit-image`` uses C-contiguous
arrays by default.
When using nested loops, the last/rightmost dimension of the array
should be in the innermost loop of the computation. In the example
above, the ``*=`` numpy operator iterates over all remaining dimensions.


A note on the time dimension
----------------------------

Although ``scikit-image`` does not currently provide functions to
work specifically with time-varying 3D data, its compatibility with
NumPy arrays allows us to work quite naturally with a 5D array of the
shape (t, pln, row, col, ch)::

    >>> for timepoint in image5d:  # doctest: +SKIP
    ...     # Each timepoint is a 3D multichannel image
    ...     do_something_with(timepoint)

We can then supplement the above table as follows:

.. table:: *Addendum to dimension names and orders in scikit-image*

  ========================   =========================================
  Image type                 coordinates
  ========================   =========================================
  2D color video             (t, row, col, ch)
  3D color video             (t, pln, row, col, ch)
  ========================   =========================================
=================================
Getting help on using ``skimage``
=================================

Besides the user guide, there exist other opportunities to get help on
using ``skimage``.

Examples gallery
----------------

The :ref:`examples_gallery` gallery provides graphical examples of
typical image processing tasks. By a quick glance at the different
thumbnails, the user may find an example close to a typical use case of
interest. Each graphical example page displays an introductory paragraph,
a figure, and the source code that generated the figure. Downloading the
Python source code enables one to modify quickly the example into a case
closer to one's image processing applications.

Users are warmly encouraged to report on their use of ``skimage`` on the
:ref:`mailing_list`, in order to propose more examples in the future.
Contributing examples to the gallery can be done on github (see
:doc:`../contribute`).

Search field
------------

The ``quick search`` field located in the navigation bar of the html
documentation can be used to search for specific keywords (segmentation,
rescaling, denoising, etc.).

API Discovery
-------------

NumPy provides a ``lookfor`` function to search API functions. 
By default ``lookfor`` will search the NumPy API.
NumPy lookfor example:
```np.lookfor('eigenvector') ```

But it can be used to search in modules, by passing in the module
name as a string:

``` np.lookfor('boundaries', 'skimage') ```

or the module itself.
```
> import skimage
> np.lookfor('boundaries', skimage)
```

Docstrings
----------

Docstrings of ``skimage`` functions are formatted using `Numpy's
documentation standard
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_,
starting with a ``Parameters`` section for the arguments and a
``Returns`` section for the objects returned by the function. Also, most
functions include one or more examples.


.. _mailing_list:

Mailing-list
------------

The scikit-image mailing-list is scikit-image@python.org (users
should `join
<https://mail.python.org/mailman3/lists/scikit-image.python.org/>`_ before posting). This
mailing-list is shared by users and developers, and it is the right
place to ask any question about ``skimage``, or in general, image
processing using Python.  Posting snippets of code with minimal examples
ensures to get more relevant and focused answers.

We would love to hear from how you use ``skimage`` for your work on the
mailing-list!
============================================
Geometrical transformations of images
============================================

Cropping, resizing and rescaling images
---------------------------------------

.. currentmodule:: skimage.transform

Images being NumPy arrays (as described in the :ref:`numpy` section), cropping
an image can be done with simple slicing operations. Below we crop a 100x100
square corresponding to the top-left corner of the astronaut image. Note that
this operation is done for all color channels (the color dimension is the last,
third dimension):: 

   >>> from skimage import data
   >>> img = data.astronaut()
   >>> top_left = img[:100, :100]

In order to change the shape of the image, :mod:`skimage.color` provides several
functions described in :ref:`sphx_glr_auto_examples_transform_plot_rescale.py`
.

.. literalinclude:: ../../examples/transform/plot_rescale.py
    :language: python
    :start-after: import matplotlib.pyplot as plt
    :end-before: fig, axes


.. image:: ../auto_examples/transform/images/sphx_glr_plot_rescale_001.png
   :target: ../auto_examples/transform/plot_rescale.html
   :align: center
   :width: 80%

Projective transforms (homographies)
------------------------------------

`Homographies <https://en.wikipedia.org/wiki/Homography>`_
are transformations of a Euclidean space that preserve the alignment of points.
Specific cases of homographies correspond to the conservation of more
properties, such as parallelism (affine transformation), shape (similar
transformation) or distances (Euclidean transformation). The different types
of homographies available in scikit-image are presented in
:ref:`sphx_glr_auto_examples_transform_plot_transform_types.py`.

Projective transformations can either be created using the explicit
parameters (e.g. scale, shear, rotation and translation)::

   from skimage import data
   from skimage import transform
   from skimage import img_as_float
   
   tform = transform.EuclideanTransform(
      rotation=np.pi / 12.,
      translation = (100, -20)
      )

or the full transformation matrix::

   from skimage import data
   from skimage import transform
   from skimage import img_as_float
   
   matrix = np.array([[np.cos(np.pi/12), -np.sin(np.pi/12), 100],
                      [np.sin(np.pi/12), np.cos(np.pi/12), -20],
                      [0, 0, 1]])
   tform = transform.EuclideanTransform(matrix)

The transformation matrix of a transform is available as its ``tform.params``
attribute. Transformations can be composed by multiplying matrices with the
``@`` matrix multiplication operator.

Transformation matrices use
`Homogeneous coordinates <https://en.wikipedia.org/wiki/Homogeneous_coordinates>`_,
which are the extension of Cartesian coordinates used in Euclidean geometry to
the more general projective geometry. In particular, points at infinity can be
represented with finite coordinates.

Transformations can be applied to images using :func:`skimage.transform.warp`::

   img = img_as_float(data.chelsea())
   tf_img = transform.warp(img, tform.inverse)

.. image:: ../auto_examples/transform/images/sphx_glr_plot_transform_types_001.png
   :target: ../auto_examples/transform/plot_transform_types.html
   :align: center
   :width: 80%

The different transformations in :mod:`skimage.transform` have a ``estimate``
method in order to estimate the parameters of the transformation from two sets
of points (the source and the destination), as explained in the
:ref:`sphx_glr_auto_examples_transform_plot_geometric.py` tutorial::

   text = data.text()

   src = np.array([[0, 0], [0, 50], [300, 50], [300, 0]])
   dst = np.array([[155, 15], [65, 40], [260, 130], [360, 95]])

   tform3 = transform.ProjectiveTransform()
   tform3.estimate(src, dst)
   warped = transform.warp(text, tform3, output_shape=(50, 300))


.. image:: ../auto_examples/transform/images/sphx_glr_plot_geometric_002.png
   :target: ../auto_examples/transform/plot_geometric.html
   :align: center
   :width: 80%


The ``estimate`` method uses least-squares optimization to minimize the distance
between source and optimization.
Source and destination points can be determined manually, or using the
different methods for feature detection available in :mod:`skimage.feature`, 
such as

 * :ref:`sphx_glr_auto_examples_features_detection_plot_corner.py`,  
 * :ref:`sphx_glr_auto_examples_features_detection_plot_orb.py`,
 * :ref:`sphx_glr_auto_examples_features_detection_plot_brief.py`,
 * etc.

and matching points using :func:`skimage.feature.match_descriptors` before
estimating transformation parameters. However, spurious matches are often made,
and it is advisable to use the RANSAC algorithm (instead of simple
least-squares optimization) to improve the robustness to outliers, as explained
in :ref:`sphx_glr_auto_examples_transform_plot_matching.py`.

.. image:: ../auto_examples/transform/images/sphx_glr_plot_matching_001.png
   :target: ../auto_examples/transform/plot_matching.html
   :align: center
   :width: 80%

Examples showing applications of transformation estimation are

 * stereo matching
   :ref:`sphx_glr_auto_examples_transform_plot_fundamental_matrix.py` and
 * image rectification :ref:`sphx_glr_auto_examples_transform_plot_geometric.py`  

The ``estimate`` method is point-based, that is, it uses only a set of points
from the source and destination images. For estimating translations (shifts),
it is also possible to use a *full-field* method using all pixels, based on
Fourier-space cross-correlation. This method is implemented by
:func:`skimage.registration.register_translation` and explained in the
:ref:`sphx_glr_auto_examples_registration_plot_register_translation.py`
tutorial.

.. image:: ../auto_examples/registration/images/sphx_glr_plot_register_translation_001.png
   :target: ../auto_examples/registration/plot_register_translation.html
   :align: center
   :width: 80%


The
:ref:`sphx_glr_auto_examples_registration_plot_register_rotation.py` tutorial
explains a variant of this full-field method for estimating a rotation, by
using first a log-polar transformation.
Data visualization
------------------

Data visualization takes an important place in image processing. Data can be
a single 2D grayscale image or a more complex one with multidimensional aspects: 3D
in space, timelapse, multiple channels.

Therefore, the visualization strategy will depend on the data complexity and
a range of tools external to scikit-image can be used for this purpose.
Historically, scikit-image provided viewer tools but powerful packages
are now available and must be preferred.


Matplotlib
^^^^^^^^^^

`Matplotlib <https://matplotlib.org/>`__ is a library able to generate static
plots, which includes image visualization.

Plotly
^^^^^^

`Plotly <https://dash.plotly.com/>`__ is a plotting library relying on web
technologies with interaction capabilities.

Mayavi
^^^^^^

`Mayavi <https://docs.enthought.com/mayavi/mayavi/>`__ can be used to visualize
3D images.

Napari
^^^^^^

`Napari <https://napari.org/>`__ is a multi-dimensional image viewer. It’s
designed for browsing, annotating, and analyzing large multi-dimensional images.
Image Segmentation
------------------

Image segmentation is the task of labeling the pixels of objects of
interest in an image.

In this tutorial, we will see how to segment objects from a background.
We use the ``coins`` image from ``skimage.data``. This image shows
several coins outlined against a darker background. The segmentation of
the coins cannot be done directly from the histogram of gray values,
because the background shares enough gray levels with the coins that a
thresholding segmentation is not sufficient.

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_001.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

::

    >>> from skimage import data
    >>> from skimage.exposure import histogram
    >>> coins = data.coins()
    >>> hist, hist_centers = histogram(coins)

Simply thresholding the image leads either to missing significant parts
of the coins, or to merging parts of the background with the
coins. This is due to the inhomogeneous lighting of the image.

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_002.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

A first idea is to take advantage of the local contrast, that is, to
use the gradients rather than the gray values.

Edge-based segmentation
~~~~~~~~~~~~~~~~~~~~~~~

Let us first try to detect edges that enclose the coins. For edge
detection, we use the `Canny detector
<https://en.wikipedia.org/wiki/Canny_edge_detector>`_ of ``skimage.feature.canny``

::

    >>> from skimage.feature import canny
    >>> edges = canny(coins/255.)

As the background is very smooth, almost all edges are found at the
boundary of the coins, or inside the coins.

::

    >>> from scipy import ndimage as ndi
    >>> fill_coins = ndi.binary_fill_holes(edges)

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_003.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

Now that we have contours that delineate the outer boundary of the coins,
we fill the inner part of the coins using the
``ndi.binary_fill_holes`` function, which uses mathematical morphology
to fill the holes.

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_004.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

Most coins are well segmented out of the background. Small objects from
the background can be easily removed using the ``ndi.label``
function to remove objects smaller than a small threshold.

::

    >>> label_objects, nb_labels = ndi.label(fill_coins)
    >>> sizes = np.bincount(label_objects.ravel())
    >>> mask_sizes = sizes > 20
    >>> mask_sizes[0] = 0
    >>> coins_cleaned = mask_sizes[label_objects]

However, the segmentation is not very satisfying, since one of the coins
has not been segmented correctly at all. The reason is that the contour
that we got from the Canny detector was not completely closed, therefore
the filling function did not fill the inner part of the coin.

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_005.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

Therefore, this segmentation method is not very robust: if we miss a
single pixel of the contour of the object, we will not be able to fill
it. Of course, we could try to dilate the contours in order to
close them. However, it is preferable to try a more robust method.

Region-based segmentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Let us first determine markers of the coins and the background. These
markers are pixels that we can label unambiguously as either object or
background. Here, the markers are found at the two extreme parts of the
histogram of gray values:

::

    >>> markers = np.zeros_like(coins)
    >>> markers[coins < 30] = 1
    >>> markers[coins > 150] = 2

We will use these markers in a watershed segmentation. The name watershed
comes from an analogy with hydrology. The `watershed transform
<https://en.wikipedia.org/wiki/Watershed_%28image_processing%29>`_ floods
an image of elevation starting from markers, in order to determine the catchment
basins of these markers. Watershed lines separate these catchment basins,
and correspond to the desired segmentation.

The choice of the elevation map is critical for good segmentation.
Here, the amplitude of the gradient provides a good elevation map. We
use the Sobel operator for computing the amplitude of the gradient::

    >>> from skimage.filters import sobel
    >>> elevation_map = sobel(coins)

From the 3-D surface plot shown below, we see that high barriers effectively
separate the coins from the background.

.. image:: data/elevation_map.jpg
    :align: center

and here is the corresponding 2-D plot:

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_006.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

The next step is to find markers of the background and the coins based on the
extreme parts of the histogram of gray values::

    >>> markers = np.zeros_like(coins)
    >>> markers[coins < 30] = 1
    >>> markers[coins > 150] = 2

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_007.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

Let us now compute the watershed transform::

    >>> from skimage.segmentation import watershed
    >>> segmentation = watershed(elevation_map, markers)

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_008.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

With this method, the result is satisfying for all coins. Even if the
markers for the background were not well distributed, the barriers in the
elevation map were high enough for these markers to flood the entire
background.

We remove a few small holes with mathematical morphology::

    >>> segmentation = ndi.binary_fill_holes(segmentation - 1)

We can now label all the coins one by one using ``ndi.label``::

    >>> labeled_coins, _ = ndi.label(segmentation)

.. image:: ../auto_examples/applications/images/sphx_glr_plot_coins_segmentation_009.png
   :target: ../auto_examples/applications/plot_coins_segmentation.html
   :align: center

========================
How to parallelize loops
========================

In image processing, we frequently apply the same algorithm
on a large batch of images. In this paragraph, we propose to
use `joblib <https://joblib.readthedocs.io>`_ to parallelize
loops. Here is an example of such repetitive tasks:

.. code-block:: python

    from skimage import data, color, util
    from skimage.restoration import denoise_tv_chambolle
    from skimage.feature import hog

    def task(image):
        """
        Apply some functions and return an image.
        """
        image = denoise_tv_chambolle(image[0][0], weight=0.1, channel_axis=-1)
        fd, hog_image = hog(color.rgb2gray(image), orientations=8,
                            pixels_per_cell=(16, 16), cells_per_block=(1, 1),
                            visualize=True)
        return hog_image


    # Prepare images
    hubble = data.hubble_deep_field()
    width = 10
    pics = util.view_as_windows(hubble, (width, hubble.shape[1], hubble.shape[2]), step=width)

To call the function ``task`` on each element of the list ``pics``, it is
usual to write a for loop. To measure the execution time of this loop, you can
use ipython and measure the execution time with ``%timeit``.

.. code-block:: python

    def classic_loop():
        for image in pics:
            task(image)


    %timeit classic_loop()

Another equivalent way to code this loop is to use a comprehension list which has the same efficiency.

.. code-block:: python

    def comprehension_loop():
        [task(image) for image in pics]

    %timeit comprehension_loop()

``joblib`` is a library providing an easy way to parallelize for loops once we have a comprehension list.
The number of jobs can be specified.

.. code-block:: python

    from joblib import Parallel, delayed
    def joblib_loop():
        Parallel(n_jobs=4)(delayed(task)(i) for i in pics)

    %timeit joblib_loop()
I/O Plugin Infrastructure
-------------------------
A plugin consists of two files, the source and the descriptor ``.ini``.  Let's
say we'd like to provide a plugin for ``imshow`` using ``matplotlib``.  We'll
call our plugin ``mpl``::

  skimage/io/_plugins/mpl.py
  skimage/io/_plugins/mpl.ini

The name of the ``.py`` and ``.ini`` files must correspond.  Inside the
``.ini`` file, we give the plugin meta-data::

  [mpl] <-- name of the plugin, may be anything
  description = Matplotlib image I/O plugin
  provides = imshow <-- a comma-separated list, one or more of
                        imshow, imsave, imread, _app_show

The "provides"-line lists all the functions provided by the plugin.  Since our
plugin provides ``imshow``, we have to define it inside ``mpl.py``::

  # This is mpl.py

  import matplotlib.pyplot as plt

  def imshow(img):
      plt.imshow(img)

Note that, by default, ``imshow`` is non-blocking, so a special function
``_app_show`` must be provided to block the GUI.  We can modify our plugin to
provide it as follows::

  [mpl]
  provides = imshow, _app_show

::

  # This is mpl.py

  import matplotlib.pyplot as plt

  def imshow(img):
      plt.imshow(img)

  def _app_show():
      plt.show()

Any plugin in the ``_plugins`` directory is automatically examined by
``skimage.io`` upon import.  You may list all the plugins on your
system::

  >>> import skimage.io as io
  >>> io.find_available_plugins()
  {'gtk': ['imshow'],
   'matplotlib': ['imshow', 'imread', 'imread_collection'],
   'pil': ['imread', 'imsave', 'imread_collection'],
   'test': ['imsave', 'imshow', 'imread', 'imread_collection'],}

or only those already loaded::

  >>> io.find_available_plugins(loaded=True)
  {'matplotlib': ['imshow', 'imread', 'imread_collection'],
   'pil': ['imread', 'imsave', 'imread_collection']}

A plugin is loaded using the ``use_plugin`` command::

  >>> import skimage.io as io
  >>> io.use_plugin('pil') # Use all capabilities provided by PIL

or

::

  >>> io.use_plugin('pil', 'imread') # Use only the imread capability of PIL

Note that, if more than one plugin provides certain functionality, the
last plugin loaded is used.

To query a plugin's capabilities, use ``plugin_info``::

  >>> io.plugin_info('pil')
  >>>
  {'description': 'Image reading via the Python Imaging Library',
   'provides': 'imread, imsave'}

Tutorials
=========

.. toctree::
   :maxdepth: 1

   tutorial_segmentation
   tutorial_parallelization
.. _data_types:

===================================
Image data types and what they mean
===================================

In ``skimage``, images are simply numpy_ arrays, which support a variety of
data types [1]_, *i.e.* "dtypes". To avoid distorting image intensities (see
`Rescaling intensity values`_), we assume that images use the following dtype
ranges:

=========  =================================
Data type  Range
=========  =================================
uint8      0 to 255
uint16     0 to 65535
uint32     0 to 2\ :sup:`32` - 1
float      -1 to 1 or 0 to 1
int8       -128 to 127
int16      -32768 to 32767
int32      -2\ :sup:`31` to 2\ :sup:`31` - 1
=========  =================================

Note that float images should be restricted to the range -1 to 1 even though
the data type itself can exceed this range; all integer dtypes, on the other
hand, have pixel intensities that can span the entire data type range. With a
few exceptions, *64-bit (u)int images are not supported*.

Functions in ``skimage`` are designed so that they accept any of these dtypes,
but, for efficiency, *may return an image of a different dtype* (see `Output
types`_). If you need a particular dtype, ``skimage`` provides utility
functions that convert dtypes and properly rescale image intensities (see
`Input types`_). You should **never use** ``astype`` on an image, because it
violates these assumptions about the dtype range::

   >>> from skimage.util import img_as_float
   >>> image = np.arange(0, 50, 10, dtype=np.uint8)
   >>> print(image.astype(float)) # These float values are out of range.
   [  0.  10.  20.  30.  40.]
   >>> print(img_as_float(image))
   [ 0.          0.03921569  0.07843137  0.11764706  0.15686275]


Input types
===========

Although we aim to preserve the data range and type of input images, functions
may support only a subset of these data-types. In such
a case, the input will be converted to the required type (if possible), and
a warning message printed to the log if a memory copy is needed. Type
requirements should be noted in the docstrings.

The following utility functions in the main package are available to developers
and users:

=============  =================================
Function name  Description
=============  =================================
img_as_float   Convert to floating point (integer types become 64-bit floats)
img_as_ubyte   Convert to 8-bit uint.
img_as_uint    Convert to 16-bit uint.
img_as_int     Convert to 16-bit int.
=============  =================================

These functions convert images to the desired dtype and *properly rescale their
values*::

   >>> from skimage.util import img_as_ubyte
   >>> image = np.array([0, 0.5, 1], dtype=float)
   >>> img_as_ubyte(image)
   array([  0, 128, 255], dtype=uint8)

Be careful! These conversions can result in a loss of precision, since 8 bits
cannot hold the same amount of information as 64 bits::

   >>> image = np.array([0, 0.5, 0.503, 1], dtype=float)
   >>> image_as_ubyte(image)
   array([  0, 128, 128, 255], dtype=uint8)

Note that ``img_as_float`` will preserve the precision of floating point types
and does not automatically rescale the range of floating point inputs.

Additionally, some functions take a ``preserve_range`` argument where a range
conversion is convenient but not necessary. For example, interpolation in
``transform.warp`` requires an image of type float, which should have a range
in [0, 1]. So, by default, input images will be rescaled to this range.
However, in some cases, the image values represent physical measurements, such
as temperature or rainfall values, that the user does not want rescaled.
With ``preserve_range=True``, the original range of the data will be
preserved, even though the output is a float image. Users must then ensure
this non-standard image is properly processed by downstream functions, which
may expect an image in [0, 1]. In general, unless a function has a
``preserve_range=False`` keyword argument, floating point inputs will not
be automatically rescaled.


    >>> from skimage import data
    >>> from skimage.transform import rescale
    >>> image = data.coins()
    >>> image.dtype, image.min(), image.max(), image.shape
    (dtype('uint8'), 1, 252, (303, 384))
    >>> rescaled = rescale(image, 0.5)
    >>> (rescaled.dtype, np.round(rescaled.min(), 4),
    ...  np.round(rescaled.max(), 4), rescaled.shape)
    (dtype('float64'), 0.0147, 0.9456, (152, 192))
    >>> rescaled = rescale(image, 0.5, preserve_range=True)
    >>> (rescaled.dtype, np.round(rescaled.min()),
    ...  np.round(rescaled.max()), rescaled.shape
    (dtype('float64'), 4.0, 241.0, (152, 192))


Output types
============

The output type of a function is determined by the function author and is
documented for the benefit of the user.  While this requires the user to
explicitly convert the output to whichever format is needed, it ensures that no
unnecessary data copies take place.

A user that requires a specific type of output (e.g., for display purposes),
may write::

   >>> from skimage.util import img_as_uint
   >>> out = img_as_uint(sobel(image))
   >>> plt.imshow(out)


Working with OpenCV
===================

It is possible that you may need to use an image created using ``skimage`` with
OpenCV_ or vice versa. OpenCV image data can be accessed (without copying) in
NumPy (and, thus, in scikit-image).
OpenCV uses BGR (instead of scikit-image's RGB) for color images, and its
dtype is uint8 by default (See `Image data types and what they mean`_). BGR stands
for Blue Green Red.

Converting BGR to RGB or vice versa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The color images in ``skimage`` and OpenCV have 3 dimensions: width, height and
color. RGB and BGR use the same color space, except the order of colors is reversed.

Note that in ``scikit-image`` we usually refer to ``rows`` and ``columns`` instead
of width and height (see :ref:`numpy-images-coordinate-conventions`).

For an image with colors along the last axis, the following instruction
effectively reverses the order of the colors, leaving the rows and columns
unaffected.

    >>> image = image[:, :, ::-1]

Using an image from OpenCV with ``skimage``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If cv_image is an array of unsigned bytes, ``skimage`` will understand it by
default. If you prefer working with floating point images, :func:`img_as_float`
can be used to convert the image::

    >>> from skimage.util import img_as_float
    >>> image = img_as_float(any_opencv_image)

Using an image from ``skimage`` with OpenCV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reverse can be achieved with :func:`img_as_ubyte`::

    >>> from skimage.util import img_as_ubyte
    >>> cv_image = img_as_ubyte(any_skimage_image)


Image processing pipeline
=========================

This dtype behavior allows you to string together any ``skimage`` function
without worrying about the image dtype.  On the other hand, if you want to use
a custom function that requires a particular dtype, you should call one of the
dtype conversion functions (here, ``func1`` and ``func2`` are ``skimage``
functions)::

   >>> from skimage.util import img_as_float
   >>> image = img_as_float(func1(func2(image)))
   >>> processed_image = custom_func(image)

Better yet, you can convert the image internally and use a simplified
processing pipeline::

   >>> def custom_func(image):
   ...     image = img_as_float(image)
   ...     # do something
   ...
   >>> processed_image = custom_func(func1(func2(image)))


Rescaling intensity values
==========================

When possible, functions should avoid blindly stretching image intensities
(e.g. rescaling a float image so that the min and max intensities are
0 and 1), since this can heavily distort an image. For example, if you're
looking for bright markers in dark images, there may be an image where no
markers are present; stretching its input intensity to span the full range
would make background noise look like markers.

Sometimes, however, you have images that should span the entire intensity
range but do not. For example, some cameras store images with 10-, 12-, or
14-bit depth per pixel. If these images are stored in an array with dtype
uint16, then the image won't extend over the full intensity range, and thus,
would appear dimmer than it should. To correct for this, you can use the
``rescale_intensity`` function to rescale the image so that it uses the full
dtype range::

   >>> from skimage import exposure
   >>> image = exposure.rescale_intensity(img10bit, in_range=(0, 2**10 - 1))

Here, the ``in_range`` argument is set to the maximum range for a 10-bit image.
By default, ``rescale_intensity`` stretches the values of ``in_range`` to match
the range of the dtype. ``rescale_intensity`` also accepts strings as inputs
to ``in_range`` and ``out_range``, so the example above could also be written
as::

   >>> image = exposure.rescale_intensity(img10bit, in_range='uint10')


Note about negative values
==========================

People very often represent images in signed dtypes, even though they only
manipulate the positive values of the image (e.g., using only 0-127 in an int8
image). For this reason, conversion functions *only spread the positive values*
of a signed dtype over the entire range of an unsigned dtype. In other words,
negative values are clipped to 0 when converting from signed to unsigned
dtypes. (Negative values are preserved when converting between signed dtypes.)
To prevent this clipping behavior, you should rescale your image beforehand::

   >>> image = exposure.rescale_intensity(img_int32, out_range=(0, 2**31 - 1))
   >>> img_uint8 = img_as_ubyte(image)

This behavior is symmetric: The values in an unsigned dtype are spread over
just the positive range of a signed dtype.


References
==========

.. _numpy: https://docs.scipy.org/doc/numpy/user/
.. [1] https://docs.scipy.org/doc/numpy/user/basics.types.html
.. _OpenCV: https://opencv.org/
Getting started
---------------

``scikit-image`` is an image processing Python package that works with
:mod:`numpy` arrays. The package is imported as ``skimage``: ::

    >>> import skimage

Most functions of ``skimage`` are found within submodules: ::

    >>> from skimage import data
    >>> camera = data.camera()

A list of submodules and functions is found on the `API reference
<https://scikit-image.org/docs/stable/api/api.html>`_ webpage.

Within scikit-image, images are represented as NumPy arrays, for
example 2-D arrays for grayscale 2-D images ::

    >>> type(camera)
    <type 'numpy.ndarray'>
    >>> # An image with 512 rows and 512 columns
    >>> camera.shape
    (512, 512)

The :mod:`skimage.data` submodule provides a set of functions returning
example images, that can be used to get started quickly on using
scikit-image's functions: ::

    >>> coins = data.coins()
    >>> from skimage import filters
    >>> threshold_value = filters.threshold_otsu(coins)
    >>> threshold_value
    107

Of course, it is also possible to load your own images as NumPy arrays
from image files, using :func:`skimage.io.imread`: ::

    >>> import os
    >>> filename = os.path.join(skimage.data_dir, 'moon.png')
    >>> from skimage import io
    >>> moon = io.imread(filename)

Use `natsort <https://pypi.org/project/natsort/>`_ to load multiple images ::

    >>> import os
    >>> from natsort import natsorted, ns
    >>> from skimage import io
    >>> list_files = os.listdir('.')
    >>> list_files
    ['01.png', '010.png', '0101.png', '0190.png', '02.png']
    >>> list_files = natsorted(list_files)
    >>> list_files
    ['01.png', '02.png', '010.png', '0101.png', '0190.png']
    >>> image_list = []
    >>> for filename in list_files:
    ...   image_list.append(io.imread(filename))
Handling Video Files
--------------------

Sometimes it is necessary to read a sequence of images from a standard video
file, such as .avi and .mov files.

In a scientific context, it is usually better to avoid these formats in favor
of a simple directory of images or a multi-dimensional TIF. Video formats are
more difficult to read piecemeal, typically do not support random frame access
or research-minded meta data, and use lossy compression if not carefully
configured. But video files are in widespread use, and they are easy to share,
so it is convenient to be equipped to read and write them when necessary.

Tools for reading video files vary in their ease of installation and use, their
disk and memory usage, and their cross-platform compatibility.  This is a
practical guide.

A Workaround: Convert the Video to an Image Sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a one-off solution, the simplest, surest route is to convert the video to a
collection of sequentially-numbered image files, often called an image
sequence. Then the images files can be read into an `ImageCollection` by
`skimage.io.imread_collection`. Converting the video to frames can be done
easily in `ImageJ <https://imagej.nih.gov/ij/>`__, a cross-platform, GUI-based
program from the bio-imaging community, or `FFmpeg <https://www.ffmpeg.org/>`__, a
powerful command-line utility for manipulating video files.

In FFmpeg, the following command generates an image file from each frame in a
video. The files are numbered with five digits, padded on the left with zeros.

.. code-block:: bash

   ffmpeg -i "video.mov" -f image2 "video-frame%05d.png"

More information is available in an `FFmpeg tutorial on image sequences 
<https://en.wikibooks.org/wiki/FFMPEG_An_Intermediate_Guide/image_sequence#Making_an_Image_Sequence_from_a_video>`__.

Generating an image sequence has disadvantages: they can be large and unwieldy,
and generating them can take some time. It is generally preferable to work
directly with the original video file. For a more direct solution, we need to
execute FFmpeg or LibAV from Python to read frames from the video.
FFmpeg and LibAV are two large open-source
projects that decode video from the sprawling variety of formats used in the
wild. There are several ways to use them from Python. Each, unfortunately,
has some disadvantages.


PyAV
^^^^

`PyAV <http://docs.mikeboers.com/pyav/develop/>`__ uses FFmpeg's (or LibAV's) libraries
to read image data directly from the video file. It invokes them using Cython
bindings, so it is very fast.

.. code-block:: python

   import av
   v = av.open('path/to/video.mov')

PyAV's API reflects the way frames are stored in a video file.

.. code-block:: python

   for packet in container.demux():
       for frame in packet.decode():
           if frame.type == 'video':
               img = frame.to_image()  # PIL/Pillow image
               arr = np.asarray(img)  # numpy array
               # Do something!

Adding Random Access to PyAV
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Video` class in `PIMS <https://github.com/soft-matter/pims>`__
invokes PyAV and adds additional functionality to solve a common
problem in scientific applications, accessing a video by frame
number. Video file formats are designed to be searched in an
approximate way, by time, and they do not support an efficient means
of seeking a specific frame number. PIMS adds this missing
functionality by decoding (but not reading) the entire video at and
producing an internal table of contents that supports indexing by
frame.

.. code-block:: python

   import pims
   v = pims.Video('path/to/video.mov')
   v[-1]  # a 2D numpy array representing the last frame

MoviePy
^^^^^^^

`Moviepy <https://zulko.github.io/moviepy>`__ invokes FFmpeg through a
subprocess, pipes the decoded video from FFmpeg
into RAM, and reads it out. This approach is straightforward, but it can be
brittle, and it's not workable for large videos that exceed available RAM.
It works on all platforms if FFmpeg is installed.

Since it does not link to FFmpeg's underlying libraries, it is easier to
install but about `half as fast <https://gist.github.com/mikeboers/6843684>`__.

.. code-block:: python

    from moviepy.editor import VideoFileClip
    myclip = VideoFileClip("some_video.avi")

Imageio
^^^^^^^^

`Imageio <http://imageio.github.io/>`_ takes the same approach as MoviePy. It
supports a wide range of other image file formats as well.

.. code-block:: python

    import imageio
    filename = '/tmp/file.mp4'
    vid = imageio.get_reader(filename,  'ffmpeg')

    for num, image in vid.iter_data():
        print(image.mean())

    metadata = vid.get_meta_data()

OpenCV
^^^^^^

Finally, another solution is the `VideoReader
<https://docs.opencv.org/2.4/modules/highgui/doc/reading_and_writing_images_and_video.html#videocapture-open>`__
class in OpenCV, which has bindings to FFmpeg. If you need OpenCV for other reasons,
then this may be the best approach.
============================================
Image adjustment: transforming image content
============================================

Color manipulation
------------------

.. currentmodule:: skimage.color

Most functions for manipulating color channels are found in the submodule
:mod:`skimage.color`.

Conversion between color models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Color images can be represented using different `color spaces
<https://en.wikipedia.org/wiki/Color_space>`_. One of the most common
color spaces is the `RGB space
<https://en.wikipedia.org/wiki/RGB_color_model>`_, where an image has red,
green and blue channels. However, other color models are widely used,
such as the `HSV color model
<https://en.wikipedia.org/wiki/HSL_and_HSV>`_, where hue, saturation and
value are independent channels, or the `CMYK model
<https://en.wikipedia.org/wiki/CMYK_color_model>`_ used for printing.

:mod:`skimage.color` provides utility functions to convert images
to and from different color spaces. Integer-type arrays can be
transformed to floating-point type by the conversion operation::

    >>> # bright saturated red
    >>> red_pixel_rgb = np.array([[[255, 0, 0]]], dtype=np.uint8)
    >>> color.rgb2hsv(red_pixel_rgb)
    array([[[ 0.,  1.,  1.]]])
    >>> # darker saturated blue
    >>> dark_blue_pixel_rgb = np.array([[[0, 0, 100]]], dtype=np.uint8)
    >>> color.rgb2hsv(dark_blue_pixel_rgb)
    array([[[ 0.66666667,  1.        ,  0.39215686]]])
    >>> # less saturated pink
    >>> pink_pixel_rgb = np.array([[[255, 100, 255]]], dtype=np.uint8)
    >>> color.rgb2hsv(pink_pixel_rgb)
    array([[[ 0.83333333,  0.60784314,  1.        ]]])



Conversion from RGBA to RGB - Removing alpha channel through alpha blending
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Converting an RGBA image to an RGB image by alpha blending it with a
background is realized with :func:`rgba2rgb` ::

    >>> from skimage.color import rgba2rgb
    >>> from skimage import data
    >>> img_rgba = data.logo()
    >>> img_rgb = rgba2rgb(img_rgba)

Conversion between color and gray values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Converting an RGB image to a grayscale image is realized with
:func:`rgb2gray` ::

    >>> from skimage.color import rgb2gray
    >>> from skimage import data
    >>> img = data.astronaut()
    >>> img_gray = rgb2gray(img)

:func:`rgb2gray` uses a non-uniform weighting of color channels, because of the
different sensitivity of the human eye to different colors. Therefore,
such a weighting ensures `luminance preservation
<https://en.wikipedia.org/wiki/Grayscale#Converting_color_to_grayscale>`_
from RGB to grayscale::

    >>> red_pixel = np.array([[[255, 0, 0]]], dtype=np.uint8)
    >>> color.rgb2gray(red_pixel)
    array([[ 0.2125]])
    >>> green_pixel = np.array([[[0, 255, 0]]], dtype=np.uint8)
    >>> color.rgb2gray(green_pixel)
    array([[ 0.7154]])


Converting a grayscale image to RGB with :func:`gray2rgb` simply
duplicates the gray values over the three color channels.

Image inversion
~~~~~~~~~~~~~~~

An inverted image is also called complementary image. For binary images, True values
become False and conversely. For grayscale images, pixel values are replaced by the
difference of the maximum value of the data type and the actual value. For RGB
images, the same operation is done for each channel. This operation can be achieved
with :py:func:`skimage.util.invert`::

    >>> from skimage import util
    >>> img = data.camera()
    >>> inverted_img = util.invert(img)

Painting images with labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`label2rgb` can be used to superimpose colors on a grayscale image
using an array of labels to encode the regions to be represented with the
same color.


.. image:: ../auto_examples/segmentation/images/sphx_glr_plot_join_segmentations_001.png
   :target: ../auto_examples/segmentation/plot_join_segmentations.html
   :align: center
   :width: 80%



.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_color_exposure_plot_tinting_grayscale_images.py`
  * :ref:`sphx_glr_auto_examples_segmentation_plot_join_segmentations.py`
  * :ref:`sphx_glr_auto_examples_segmentation_plot_rag_mean_color.py`


.. _exposure:

Contrast and exposure
---------------------

.. currentmodule:: skimage.exposure

Image pixels can take values determined by the ``dtype`` of the image
(see :ref:`data_types`), such as 0 to 255 for ``uint8`` images or ``[0,
1]`` for floating-point images. However, most images either have a
narrower range of values (because of poor contrast), or have most pixel
values concentrated in a subrange of the accessible values.
:mod:`skimage.exposure` provides functions that spread the intensity
values over a larger range.

A first class of methods compute a nonlinear function of the intensity,
that is independent of the pixel values of a specific image. Such methods
are often used for correcting a known non-linearity of sensors, or
receptors such as the human eye. A well-known example is `Gamma
correction <https://en.wikipedia.org/wiki/Gamma_correction>`_, implemented
in :func:`adjust_gamma`.

Other methods re-distribute pixel values according to the *histogram* of
the image. The histogram of pixel values is computed with
:func:`skimage.exposure.histogram`::

    >>> image = np.array([[1, 3], [1, 1]])
    >>> exposure.histogram(image)
    (array([3, 0, 1]), array([1, 2, 3]))

:func:`histogram` returns the number of pixels for each value bin, and
the centers of the bins. The behavior of :func:`histogram` is therefore
slightly different from the one of :func:`numpy.histogram`, which returns
the boundaries of the bins.

The simplest contrast enhancement :func:`rescale_intensity` consists in
stretching pixel values to the whole allowed range, using a linear
transformation::

    >>> from skimage import exposure
    >>> text = data.text()
    >>> text.min(), text.max()
    (10, 197)
    >>> better_contrast = exposure.rescale_intensity(text)
    >>> better_contrast.min(), better_contrast.max()
    (0, 255)

Even if an image uses the whole value range, sometimes there is very
little weight at the ends of the value range. In such a case, clipping
pixel values using percentiles of the image improves the contrast (at the
expense of some loss of information, because some pixels are saturated by
this operation)::

    >>> moon = data.moon()
    >>> v_min, v_max = np.percentile(moon, (0.2, 99.8))
    >>> v_min, v_max
    (10.0, 186.0)
    >>> better_contrast = exposure.rescale_intensity(
    ...                                     moon, in_range=(v_min, v_max))

The function :func:`equalize_hist` maps the cumulative distribution
function (cdf) of pixel values onto a linear cdf, ensuring that all parts
of the value range are equally represented in the image. As a result,
details are enhanced in large regions with poor contrast. As a further
refinement, histogram equalization can be performed in subregions of the
image with :func:`equalize_adapthist`, in order to correct for exposure
gradients across the image. See the example
:ref:`sphx_glr_auto_examples_color_exposure_plot_equalize.py`.

.. image:: ../auto_examples/color_exposure/images/sphx_glr_plot_equalize_001.png
   :target: ../auto_examples/color_exposure/plot_equalize.html
   :align: center
   :width: 90%


.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_color_exposure_plot_equalize.py`


.. _skip_3_transition_v1:

==========================================
SKIP 3 — Transitioning to scikit-image 1.0
==========================================

:Author: Juan Nunez-Iglesias <juan.nunez-iglesias@monash.edu>
:Status: Final
:Type: Standards Track
:Created: 2021-07-15
:Resolved: 2021-09-13
:Resolution: Rejected
:Version effective: None

Abstract
--------

scikit-image is preparing to release version 1.0. This is potentially an
opportunity to clean up the API, including backwards incompatible changes. Some
of these changes involve changing return values without changing function
signatures, which can ordinarily only be done by adding an otherwise useless
keyword argument (such as ``new_return_style=True``) whose default value
changes over several releases. The result is *still* a backwards incompatible
change, but made over a longer time period.

Despite being in beta and in a 0.x series of releases, scikit-image is used
extremely broadly, and any backwards incompatible changes are likely to be
disruptive. This SKIP proposes a process to ensure that the community is aware
of upcoming changes, and can adapt their libraries *or* their declared
scikit-image version dependencies accordingly.

Motivation and Scope
--------------------

scikit-image has grown organically over the past 12 years, with functionality
being added by a broad community of contributors from different backgrounds.
This has resulted in various parts of the API being inconsistent: for example,
``skimage.transform.warp`` inverts the order of coordinates, so that a
translation of (45, 32) actually moves the values in a NumPy array by 32 along
the 0th axis, and 45 along the 1st, *but only in 2D*.

Additionally, as our user base has grown, it has become apparent that certain
early API choices turned out to be more confusing than helpful. For example,
scikit-image will automatically convert images to various data types,
*rescaling them in the process*. A uint8 image in the range [0, 255] will
automatically be converted to a float64 image in [0, 1]. This might initially
seem reasonable, but, for consistency, uint16 images in [0, 65535] are rescaled
to [0, 1] floats, and uint16 images with 12-bit range in [0, 4095], which are
common in microscopy, are rescaled to [0, 0.0625]. These silent conversions
have resulted in much user confusion.

Changing this convention would require adding a ``preserve_range=`` keyword
argument to almost *all* scikit-image functions, whose default value would
change from False to True over 4 versions. Eventually, the change would be
backwards-incompatible, no matter how gentle we made the deprecation curve.

Given the accumulation of potential API changes that have turned out to be too
burdensome and noisy to fix with a standard deprecation cycle, principally
because they involve changing function outputs for the same inputs, it makes
sense to make all those changes in a transition to version 1.0 -- semantic
versioning, which we use, explicitly allows breaking API changes on major
version updates [6]_. However, we must acknowledge that (1) an enormous number
of projects depend on scikit-image and would thus be affected by backwards
incompatible changes, and (2) it is not yet common practice in the scientific
Python community to put upper version bounds on dependencies, so it is very
unlikely that anyone used ``scikit-image<1.*`` in their dependency list (though
this is slowly changing [5]_).

Given the above, we need to come up with a way to notify all our users that
this change is coming, while also allowing them to silence any warnings once
they have been noted.

Detailed description
--------------------

It is beyond the scope of this document to list all of the proposed API changes
for scikit-image 1.0, many of which have yet to be decided upon. Indeed, the
scope and ambition of the 1.0 transition could grow if this SKIP is accepted.
The SKIP instead proposes a mechanism for warning users about upcoming breaking
changes. A meta-issue tracking the proposed changes can be found on GitHub,
scikit-image/scikit-image#5439 [7]_. Some examples are briefly included below
for illustrative purposes:

- Stop rescaling input arrays when the dtype must be coerced to float.
- Stop swapping coordinate axis order in different contexts, such as drawing or
  warping.
- Allow automatic return of non-NumPy types, so long as they are coercible to
  NumPy with ``numpy.asarray``.
- Harmonizing similar parameters in different functions to have the same name;
  for example, we currently have ``random_seed``, ``random_state``, ``seed``,
  or ``sample_seed`` in different functions, all to mean the same thing.
- Changing ``measure.regionprops`` to return a dictionary instead of a list.
- Combine functions that have the same purpose, such as ``watershed``,
  ``slic``, or ``felzenschwalb``, into a common namespace. This would make it
  easier for new users to find out which functions they should try out for a
  specific task.

The question is, how do we make this transition while causing as little
disruption as possible?

This document proposes releasing 0.19 as the final 0.x series release, then
immediately releasing a nearly identical 0.20 release that warns users about
breaking changes in 1.0, thus giving them an opportunity to pin their
scikit-image dependency to 0.19.x. The warning would also point users to a
transition guide to prepare their code for 1.0. See `Implementation` for
details.

This approach ensures that all users get ample warning, and a chance to ensure
that their scripts and libraries will continue to work after 1.0 is released.
Users who don't have the time or inclination to make the transition will be
able to pin their dependencies correctly. Those who prefer to be on the cutting
edge will also be able to plan around the 1.0 release and update their code
correctly, in sync with scikit-image.

Related Work
------------

`pandas` released 1.0.0 in January 2020, including many backwards-incompatible
API changes [3]_. `scipy` released version 1.0 in 2017, but, given its stage of
maturity and position at the base of the scientific Python ecosystem, opted not
to make major breaking changes [4]_. However, SciPy has adopted a policy of
adding upper-bounds on dependencies [5]_, acknowledging that the ecosystem as a
whole makes backwards incompatible changes on a 2 version deprecation cycle.

Implementation
--------------

The details of the proposal are as follows:

- scikit-image 0.19 will be the final *true* 0.x release. It contains some new
  features, bug fixes, and several API changes following on from deprecations
  in 0.17.
- shortly after 0.19, we release 0.20, which is identical except that it emits
  a warning at import time. The warning reads something like the following:
  "scikit-image 1.0 will be released later this year and will contain breaking
  changes. To ensure your code keeps running, please install
  ``scikit-image<=0.19.*``. To silence this warning but still depend on
  scikit-image after 1.0 is released, install ``scikit-image!=0.20.*``." The
  warning also contains a link for further details, and instructions for
  managing the dependency in both conda and pip environments.
- After 0.20, we make all the API changes we need, without deprecation cycles.
  Importantly, for every API change, we add a line to a "scikit-image 1.0
  transition guide" in the documentation, which maps every changed
  functionality in the library from its old form to its new form. These changes
  are tracked on a GitHub issue [7]_ and in the 1.0 milestone [8]_.
- Once the transition has happened in the repository, we release 1.0.0a0, an
  alpha release which contains a global warning pointing to the transition
  guide, as well as all of the new functionality. We also release 0.21, which
  contains the same warning but is functionally identical to 0.19. This gives
  authors who chose to pin to ``scikit-image!=0.20.*`` a chance to make the
  migration to 1.0.
- After at least one month, we release 1.0.
- We continue to maintain a 0.19.x branch with bug fixes for a year, in order
  to give users time to transition to the new API.

Backward compatibility
----------------------

This proposal breaks backwards compatibility in numerous places in the library.

Alternatives
------------

New package naming
..................

Instead of breaking compatibility in the ``scikit-image`` package, we could
leave that package at 0.19, and release a *new* package, e.g.
``scikit-image1``, which starts at 1.0 and imports as ``skimage1``. This would
obviate the need for users to pin their scikit-image version — users depending
on skimage 0.x would be able to use that library "in perpetuity."

Ultimately, the core developers felt that this approach could unnecessarily
fragment the community, between those that continue using 0.19 and those that
shift to 1.0. Ultimately, the transition of downstream code to 1.0 would be
equally painful as the proposed approach, but the pressure to make the switch
would be decreased, as everyone installing ``scikit-image`` would still get the
old version.

Continuous deprecation over multiple versions
.............................................

This transition could occur gradually over many versions. For example, for
functions automatically converting and rescaling float inputs, we could add a
``preserve_range`` keyword argument that would initially default to False, but
the default value of False would be deprecated, with users getting a warning to
switch to True. After the switch, we could (optionally) deprecate the
argument, arriving, after a further two releases, at the same place:
scikit-image no longer rescales data automatically, there are no
unnecessary keyword arguments lingering all over the API.

Of course, this kind of operation would have to be done simultaneously over all
of the above proposed changes.

Ultimately, the core team felt that this approach generates more work for both
the scikit-image developers and the developers of downstream libraries, for
dubious benefit: ultimately, later versions of scikit-image will still be
incompatible with prior versions, although over a longer time scale.

Not making the proposed API changes
...................................

Another possibility is to reject backwards incompatible API changes outright,
except in extreme cases. The core team feels that this is essentially
equivalent to pinning the library at 0.19.

Discussion
----------

In early July 2021, the core team held a series of meetings to discuss this
approach. The minutes of this meeting are in the scikit-image meeting notes
repository [9]_.

Ongoing discussion will happen on the user forum [10]_, the
developer forum [11]_, and GitHub discussion [7]_. Specific links to relevant
posts will be added to this document before acceptance.

Resolution
----------

This SKIP was discussed most extensively in a thread on the mailing list in
July 2021 [12]_. In the end, many and core developers felt that this plan
posed too big a risk of either changing code behavior silently or eroding
goodwill in the community, or both. Matthew Brett wrote [13]_:

    I'm afraid I wasn't completely sure whether the 1.0 option would
    result in breaking what I call the Konrad Hinsen rule for scientific
    software:

    """
    Under (virtually) no circumstances should new versions of a scientific
    package silently give substantially different results for the same
    function / method call from a previous version of the package.
    """

Matthew further wrote [14]_ that if we *don't* break the Hinsen rule, but
instead break users' unpinned scripts, we will lose a lot of goodwill from the
community:

    If you make all these break (if they are lucky) or give completely
    wrong results, it's hard to imagine you aren't going to cause
    significant damage to the rest-of-iceberg body of users who are not on
    the mailing list.

Riadh Fezzani, one of our core developers, felt strongly that SemVer [6]_ was
sufficient to protect users [15]_:

    In scikit-image, we adopted the semantic versioning as it
    is largely adopted in the engineering community. This convention manages
    API breaking and that's what we are doing by releasing v1.0

Even taking this view, though, it cannot address the issue of external
scikit-image "documentation", such as a decade's worth of accumulated
StackOverflow answers, that would be made obsolete by a breaking 1.0 release,
as pointed out by Josh Warner [16]_:

    It's also worth considering that there is a substantial corpus of
    scikit-image teaching material out there. The majority we do not control,
    so cannot be updated or edited. The first hits on YouTube for tutorials
    are not the most recent, but older ones with lots of views.

Nor can it address the issue of *gradually* migrating a code base from the old
API to the new API, as pointed out by Tom Caswell [17]_:

    Put another way, you do not want to put a graduate student in the position
    of saying "I _want_ to use the new API, but I have 10k LoC of inherited
    code using the old API .....".

Ultimately, all these concerns add up to a compelling case to rejecting the
SKIP. Juan Nunez-Iglesias wrote on the mailing list [18]_:

    My proposal going forward is to reject SKIP-3 and create a SKIP-4 proposing
    the skimage2 package.

The SKIP is therefore rejected.

References and Footnotes
------------------------

All SKIPs should be declared as dedicated to the public domain with the CC0
license [1]_, as in `Copyright`, below, with attribution encouraged with CC0+BY
[2]_.

.. [1] CC0 1.0 Universal (CC0 1.0) Public Domain Dedication,
   https://creativecommons.org/publicdomain/zero/1.0/
.. [2] https://dancohen.org/2013/11/26/cc0-by/
.. [3] https://pandas.pydata.org/pandas-docs/stable/whatsnew/v1.0.0.html#backwards-incompatible-api-changes
.. [4] https://docs.scipy.org/doc/scipy/reference/release.1.0.0.html
.. [5] https://github.com/scipy/scipy/pull/12862
.. [6] https://semver.org/
.. [7] https://github.com/scikit-image/scikit-image/discussions/5651
.. [8] https://github.com/scikit-image/scikit-image/milestones/1.0
.. [9] https://github.com/scikit-image/meeting-notes/blob/main/2021/july-api-meetings.md
.. [10] https://forum.image.sc/tag/scikit-image
.. [11] https://discuss.scientific-python.org/c/contributor/skimage
.. [12] https://mail.python.org/archives/list/scikit-image@python.org/thread/DSV6PEYVJ4RZRUWWV5SBNF7FFRERTSCF/
.. [13] https://mail.python.org/archives/list/scikit-image@python.org/message/UYARUQM5LBWXIAWBAPNHIQIDRKUUDTEK/
.. [14] https://mail.python.org/archives/list/scikit-image@python.org/message/63ZGG7DY5SWVM62XASHMCPFAG6KPJCMT/
.. [15] https://mail.python.org/archives/list/scikit-image@python.org/message/HXI7YVCN6IFF5TL54JBP5QRUDHKTTYRR/
.. [16] https://mail.python.org/archives/list/scikit-image@python.org/message/HRZGMOJLD2WDIO3JXQV3PRWKIUOVOF7P/
.. [17] https://mail.python.org/archives/list/scikit-image@python.org/message/GFXBQYKDACDCH7BGNEGOU7LKHR2LPFX6/
.. [18] https://mail.python.org/archives/list/scikit-image@python.org/message/5J4W63BXFQTT4GHPTZFH3AM4QHAXOW5R/

Copyright
---------

This document is dedicated to the public domain with the Creative Commons CC0
license [1]_. Attribution to this source is encouraged where appropriate, as per
CC0+BY [2]_.
.. _governance:

====================================================
SKIP 1 — scikit-image governance and decision-making
====================================================

:Author: Juan Nunez-Iglesias <juan.nunez-iglesias@monash.edu>
:Author: Stéfan van der Walt <stefanv@berkeley.edu>
:Author: Josh Warner
:Author: François Boulogne
:Author: Emmanuelle Gouillart
:Author: Mark Harfouche
:Author: Lars Grüter
:Author: Egor Panfilov
:Status: Final
:Type: Process
:Created: 2019-07-02
:Resolved: 2019-09-25
:Resolution: https://github.com/scikit-image/scikit-image/pull/4182
:skimage-Version: 0.16

Abstract
========

The purpose of this document is to formalize the governance process used by the
scikit-image project, to clarify how decisions are made and how the various
elements of our community interact.

This is a consensus-based community project. Anyone with an interest in the
project can join the community, contribute to the project design, and
participate in the decision making process. This document describes how that
participation takes place, how to find consensus, and how deadlocks are
resolved.

Roles And Responsibilities
==========================

The Community
-------------
The scikit-image community consists of anyone using or working with the project
in any way.

Contributors
------------
A community member can become a contributor by interacting directly with the
project in concrete ways, such as:

- proposing a change to the code via a GitHub pull request;
- reporting issues on our
  `GitHub issues page <https://github.com/scikit-image/scikit-image/issues>`_;
- proposing a change to the documentation,
  `website <https://github.com/scikit-image/scikit-image-web>`_, or
  `tutorials <https://github.com/scikit-image/skimage-tutorials>`_ via a
  GitHub pull request;
- discussing the design of the library, website, or tutorials on the
  `developer forum <https://discuss.scientific-python.org/c/contributor/skimage>`_,
  in the
  `project chat room <https://skimage.zulipchat.com/>`_, or in existing issues and pull
  requests; or
- reviewing
  `open pull requests <https://github.com/scikit-image/scikit-image/pulls>`_,

among other possibilities. Any community member can become a contributor, and
all are encouraged to do so. By contributing to the project, community members
can directly help to shape its future.

Contributors are encouraged to read the
:doc:`contributing guide <../contribute>`.

Core developers
---------------
Core developers are community members that have demonstrated continued
commitment to the project through ongoing contributions. They
have shown they can be trusted to maintain scikit-image with care. Becoming a
core developer allows contributors to merge approved pull requests, cast votes
for and against merging a pull-request, and be involved in deciding major
changes to the API, and thereby more easily carry on with their project related
activities. Core developers appear as organization members on the scikit-image
`GitHub organization <https://github.com/orgs/scikit-image/people>`_. Core
developers are expected to review code contributions while adhering to the
:ref:`core developer guide <core_dev>`.

New core developers can be nominated by any existing core developer.
Discussion about new core developer nominations is one of the few activities
that takes place on the project's private management list. The decision to
invite a new core developer must be made by “lazy consensus”, meaning unanimous
agreement by all responding existing core developers. Invitation must take
place at least one week after initial nomination, to allow existing members
time to voice any objections.

Steering Council
----------------
The Steering Council (SC) members are core developers who have additional
responsibilities to ensure the smooth running of the project. SC members are
expected to participate in strategic planning, approve changes to the
governance model, and make decisions about funding granted to the project
itself. (Funding to community members is theirs to pursue and manage.) The
purpose of the SC is to ensure smooth progress from the big-picture
perspective. Changes that impact the full project require analysis informed by
long experience with both the project and the larger ecosystem. When the core
developer community (including the SC members) fails to reach such a consensus
in a reasonable timeframe, the SC is the entity that resolves the issue.

The steering council is fixed in size to five members. This may be expanded in
the future. The initial steering council of scikit-image consists of Stéfan
van der Walt, Juan Nunez-Iglesias, Emmanuelle Gouillart, Josh Warner, and
Zachary Pincus. The SC membership is revisited every January. SC members who do
not actively engage with the SC duties are expected to resign. New members are
added by nomination by a core developer. Nominees should have demonstrated
long-term, continued commitment to the project and its :doc:`values <../values>`. A
nomination will result in discussion that cannot take more than a month and
then admission to the SC by consensus.

The scikit-image steering council may be contacted at
`skimage-steering@groups.io <mailto:skimage-steering@groups.io>`__.

Decision Making Process
=======================

Decisions about the future of the project are made through discussion with all
members of the community. All non-sensitive project management discussion takes
place on the project
`developer forum <https://discuss.scientific-python.org/c/contributor/skimage>`_
and the `issue tracker <https://github.com/scikit-image/scikit-image/issues>`_.
Occasionally, sensitive discussion may occur on a private list.

Decisions should be made in accordance with the :doc:`mission, vision and
values <../values>` of the scikit-image project.

Scikit-image uses a “consensus seeking” process for making decisions. The group
tries to find a resolution that has no open objections among core developers.
Core developers are expected to distinguish between fundamental objections to a
proposal and minor perceived flaws that they can live with, and not hold up the
decision-making process for the latter.  If no option can be found without
objections, the decision is escalated to the SC, which will itself use
consensus seeking to come to a resolution. In the unlikely event that there is
still a deadlock, the proposal will move forward if it has the support of a
simple majority of the SC. Any vote must be backed by a :ref:`scikit-image
proposal (SKIP) <skip>`.

Decisions (in addition to adding core developers and SC membership as above)
are made according to the following rules:

- **Minor documentation changes**, such as typo fixes, or addition / correction of a
  sentence (but no change of the scikit-image.org landing page or the “about”
  page), require approval by a core developer *and* no disagreement or requested
  changes by a core developer on the issue or pull request page (lazy
  consensus). Core developers are expected to give “reasonable time” to others
  to give their opinion on the pull request if they’re not confident others
  would agree.

- **Code changes and major documentation changes** require agreement by *two*
  core developers *and* no disagreement or requested changes by a core developer
  on the issue or pull-request page (lazy consensus).

- **Changes to the API principles** require a :ref:`skip` and follow the
  decision-making process outlined above.

- **Changes to this governance model or our mission, vision, and values**
  require a :ref:`skip` and follow the decision-making process outlined above,
  *unless* there is unanimous agreement from core developers on the change.

If an objection is raised on a lazy consensus, the proposer can appeal to the
community and core developers and the change can be approved or rejected by
escalating to the SC, and if necessary, a SKIP (see below).

.. _skip:

Improvement proposals (SKIPs)
=============================

For all votes, a formal proposal must have been made public and discussed
before the vote. After discussion has taken place, the key advocate of the
proposal must create a consolidated document summarizing the discussion, called
a SKIP, on which the core team votes. The lifetime of a SKIP detailed in
:ref:`skip0`.

A list of all existing SKIPs is available :ref:`here <skip_list>`.

Copyright
=========

This document is based on the `scikit-learn governance document
<https://scikit-learn.org/stable/governance.html>`_ and is placed in the public
domain.
.. _skip0:

============================
SKIP 0 — Purpose and Process
============================

:Author: Jarrod Millman <millman@berkeley.edu>
:Author: Juan Nunez-Iglesias <juan.nunez-iglesias@monash.edu>
:Author: Stéfan van der Walt <stefanv@berkeley.edu>
:Status: Active
:Type: Process
:Created: 2019-07-30


What is a SKIP?
---------------

SKIP stands for **s**\ ci\ **k**\ it-\ **i**\ mage **p**\ roposal. A SKIP is a design document providing
information to the community, or describing a new feature for
scikit-image or its processes or environment. The SKIP should provide a
rationale for the proposed change as well as a concise technical
specification, if applicable.

We intend SKIPs to be the primary mechanisms for proposing major new
features, for collecting community input on an issue, and for
documenting the design decisions that have gone into scikit-image. The SKIP
author is responsible for building consensus within the community and
documenting dissenting opinions.

Because the SKIPs are maintained as text files in a versioned
repository, their revision history is the historical record of the
feature proposal [1]_.


Types
^^^^^

There are three kinds of SKIPs:

1. A **Standards Track** SKIP describes a new feature or implementation
   for scikit-image.

2. An **Informational** SKIP describes a scikit-image design issue, or provides
   general guidelines or information to the Python community, but does not
   propose a new feature. Informational SKIPs do not necessarily represent a
   scikit-image community consensus or recommendation, so users and
   implementers are free to ignore Informational SKIPs.

3. A **Process** SKIP describes a process surrounding scikit-image, or
   proposes a change to (or an event in) a process. Process SKIPs are
   like Standards Track SKIPs but apply to areas other than the scikit-image
   library itself. They may propose an implementation, but not to
   scikit-image's codebase; they require community consensus. Examples include
   procedures, guidelines, changes to the decision-making process, and changes
   to the tools or environment used in scikit-image development. Any meta-SKIP
   is also considered a Process SKIP.


SKIP Workflow
-------------

The SKIP process begins with a new idea for scikit-image. A SKIP should contain
a single key proposal or new idea. Small enhancements or patches often don't
need a SKIP and can be injected into the scikit-image development workflow
with a pull request to the scikit-image `repo`_. The more focused the
SKIP, the more likely it is to be accepted.

Each SKIP must have a champion---someone who writes the SKIP using the style
and format described below, shepherds the discussions in the appropriate
forums, and attempts to build community consensus around the idea.  The SKIP
champion (a.k.a. Author) should first attempt to ascertain whether the idea is
suitable for a SKIP. Posting to the scikit-image `developer forum`_ is the best
way to do this.

The proposal should be submitted as a draft SKIP via a `GitHub pull
request`_ to the ``doc/source/skips`` directory with the name
``skip-<n>.rst`` where ``<n>`` is an appropriately assigned number (e.g.,
``skip-35.rst``). The draft must use the :ref:`skip_template` file.

Once the PR is in place, the SKIP should be announced on the developer
forum for discussion (comments on the PR itself should be restricted to
minor editorial and technical fixes).

At the earliest convenience, the PR should be merged (regardless of whether it
is accepted during discussion). A SKIP that outlines a coherent argument and
that is considered reasonably complete should be merged optimistically,
regardless of whether it is accepted during discussion. Additional PRs may be
made by the author to update or expand the SKIP, or by maintainers to set its
status, discussion URL, etc.

Standards Track SKIPs consist of two parts, a design document and a
reference implementation. It is generally recommended that at least a
prototype implementation be co-developed with the SKIP, as ideas that sound
good in principle sometimes turn out to be impractical. Often it makes sense
for the prototype implementation to be made available as PR to the scikit-image
repo, as long as it is properly marked as WIP (work in progress).


Review and Resolution
^^^^^^^^^^^^^^^^^^^^^

SKIPs are discussed on the developer forum. The possible paths of the
status of SKIPs are as follows:

.. image:: _static/skip-flowchart.png

All SKIPs should be created with the ``Draft`` status.

Eventually, after discussion, there may be a consensus that the SKIP
should be accepted – see the next section for details. At this point
the status becomes ``Accepted``.

Once a SKIP has been ``Accepted``, the reference implementation must be
completed. When the reference implementation is complete and incorporated
into the main source code repository, the status will be changed to ``Final``.

To allow gathering of additional design and interface feedback before
committing to long term stability for a language feature or standard library
API, a SKIP may also be marked as "Provisional". This is short for
"Provisionally Accepted", and indicates that the proposal has been accepted for
inclusion in the reference implementation, but additional user feedback is
needed before the full design can be considered "Final". Unlike regular
accepted SKIPs, provisionally accepted SKIPs may still be Rejected or Withdrawn
even after the related changes have been included in a release.

Wherever possible, it is considered preferable to reduce the scope of a
proposal to avoid the need to rely on the "Provisional" status (e.g. by
deferring some features to later SKIPs), as this status can lead to version
compatibility challenges in the wider ecosystem.

A SKIP can also be assigned status ``Deferred``. The SKIP author or a
core developer can assign the SKIP this status when no progress is being made
on the SKIP.

A SKIP can also be ``Rejected``. Perhaps after all is said and done it
was not a good idea. It is still important to have a record of this
fact. The ``Withdrawn`` status is similar---it means that the SKIP author
themselves has decided that the SKIP is actually a bad idea, or has
accepted that a competing proposal is a better alternative.

When a SKIP is ``Accepted``, ``Rejected``, or ``Withdrawn``, the SKIP should be
updated accordingly. In addition to updating the status field, at the very
least the ``Resolution`` header should be added with a link to the relevant
post on the discussion forum.

SKIPs can also be ``Superseded`` by a different SKIP, rendering the
original obsolete. The ``Replaced-By`` and ``Replaces`` headers
should be added to the original and new SKIPs respectively.

Process SKIPs may also have a status of ``Active`` if they are never
meant to be completed, e.g. SKIP 0 (this SKIP).


How a SKIP becomes Accepted
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A SKIP is ``Accepted`` by consensus of all interested contributors. We
need a concrete way to tell whether consensus has been reached. When
you think a SKIP is ready to accept, start a topic on the
developer forum with a subject like:

  Proposal to accept SKIP #<number>: <title>

In the body of your email, you should:

* link to the latest version of the SKIP,

* briefly describe any major points of contention and how they were
  resolved,

* include a sentence like: "If there are no substantive objections
  within 7 days from this email, then the SKIP will be accepted; see
  SKIP 0 for more details."

For an equivalent example in the NumPy library, see: https://mail.python.org/pipermail/numpy-discussion/2018-June/078345.html

After you send the email, you should make sure to link to the email
thread from the ``Discussion`` section of the SKIP, so that people can
find it later.

Generally the SKIP author will be the one to send this email, but
anyone can do it – the important thing is to make sure that everyone
knows when a SKIP is on the verge of acceptance, and give them a final
chance to respond. If there's some special reason to extend this final
comment period beyond 7 days, then that's fine, just say so in the
email. You shouldn't do less than 7 days, because sometimes people are
travelling or similar and need some time to respond.

In general, the goal is to make sure that the community has consensus,
not provide a rigid policy for people to try to game. When in doubt,
err on the side of asking for more feedback and looking for
opportunities to compromise.

If the final comment period passes without any substantive objections,
then the SKIP can officially be marked ``Accepted``. You should send a
followup email notifying the list (celebratory emoji optional but
encouraged 🎉✨), and then update the SKIP by setting its ``:Status:``
to ``Accepted``, and its ``:Resolution:`` header to a link to your
followup email.

If there *are* substantive objections, then the SKIP remains in
``Draft`` state, discussion continues as normal, and it can be
proposed for acceptance again later once the objections are resolved.

In unusual cases, when no consensus can be reached between core developers, the
`scikit-image Steering Council`_ may be asked to decide whether a controversial
SKIP is ``Accepted``.


Maintenance
^^^^^^^^^^^

In general, Standards track SKIPs are no longer modified after they have
reached the Final state, as the code and project documentation are considered
the ultimate reference for the implemented feature. They may, however, be
updated under special circumstances.

Process SKIPs may be updated over time to reflect changes
to development practices and other details. The precise process followed in
these cases will depend on the nature and purpose of the SKIP being updated.


Format and Template
-------------------

SKIPs are UTF-8 encoded text files using the reStructuredText_ format.  Please
see the :ref:`skip_template` file and the reStructuredTextPrimer_ for more
information.  We use Sphinx_ to convert SKIPs to HTML for viewing on the web
[2]_.


Header Preamble
^^^^^^^^^^^^^^^

Each SKIP must begin with a header preamble.  The headers
must appear in the following order.  Headers marked with ``*`` are
optional.  All other headers are required. ::

    :Author: <list of authors' real names and optionally, email addresses>
    :Status: <Draft | Active | Accepted | Deferred | Rejected |
             Withdrawn | Final | Superseded>
    :Type: <Standards Track | Process>
    :Created: <date created on, in dd-mmm-yyyy format>
  * :Requires: <skip numbers>
  * :skimage-Version: <version number>
  * :Replaces: <skip number>
  * :Replaced-By: <skip number>
  * :Resolution: <url>

The Author header lists the names, and optionally the email addresses
of all the authors of the SKIP.  The format of the Author header
value must be

    Random J. User <address@dom.ain>

if the email address is included, and just

    Random J. User

if the address is not given.  If there are multiple authors, each should be on
a separate line.


Discussion
----------

- https://github.com/scikit-image/scikit-image/pull/3585


References and Footnotes
------------------------

.. [1] This historical record is available by the normal git commands
   for retrieving older revisions, and can also be browsed on
   `GitHub <https://github.com/scikit-image/scikit-image/tree/main/doc/source/skips>`_.

.. [2] The URL for viewing SKIPs on the web is
   https://scikit-image.org/docs/stable/skips/

.. _repo: https://github.com/scikit-image/scikit-image

.. _developer forum: https://discuss.scientific-python.org/c/contributor/skimage

.. _issue tracker: https://github.com/scikit-image/scikit-image/issues

.. _scikit-image Steering Council:
   https://scikit-image.org/docs/stable/skips/1-governance.html

.. _`GitHub pull request`: https://github.com/scikit-image/scikit-image/pulls

.. _reStructuredText: https://docutils.sourceforge.io/rst.html

.. _reStructuredTextPrimer: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

.. _Sphinx: http://www.sphinx-doc.org/en/stable/


Copyright
---------

This document has been placed in the public domain.
=======================================
SKIP 2 — scikit-image mission statement
=======================================

:Author: Juan Nunez-Iglesias <juan.nunez-iglesias@monash.edu>
:Author: Stéfan van der Walt <stefanv@berkeley.edu>
:Author: Josh Warner
:Author: François Boulogne
:Author: Emmanuelle Gouillart
:Author: Mark Harfouche
:Author: Lars Grüter
:Author: Egor Panfilov
:Author: Gregory Lee
:Status: Active
:Type: Process
:Created: 2018-12-08
:Resolved:
:Resolution:
:skimage-Version: 0.16

Abstract
--------

scikit-image should adopt the document below as its mission statement. This
statement will feature prominently in the scikit-image home page and readme,
as well as the contributor and core developer guides. Decisions about the API
and the future of the library would be referenced against this document. (See
:ref:`governance`.)

In July 2018, I (Juan) published a blog post that broadly outlined what I would
want from a roadmap for scikit-image [1]_, but requested comments from the
community before it would be finalized. I consider that we have collected
comments for long enough and can move forward with formal adoption. Most
comments were positive, so below I'll just summarize the “negative” comments
under “rejected ideas”.

Detailed description
--------------------

(Or: What problem does this proposal solve?)

Over the past few years, scikit-image has been slightly “adrift”, with new and
old contributors coming in, adding what small bits they need at the time, and
disappearing again for a while. This is *fine* and we want to encourage more of
it, but it also lacks direction. Additionally, without a concerted roadmap to
concentrate effort, many of these contributions just fall by the wayside, as it
is difficult for new contributors to live up to our stringent (and largely
unwritten) standards of code.

Implementation
--------------

Our mission
***********

scikit-image aims to be the reference library for scientific image analysis in
Python. We accomplish this by:

- being **easy to use and install**. We are careful in taking on new
  dependencies, and sometimes cull existing ones, or make them optional. All
  functions in our API have thorough docstrings clarifying expected inputs and
  outputs.
- providing a **consistent API**. Conceptually identical arguments have the
  same name and position in a function signature.
- **ensuring correctness**. Test coverage is close to 100% and code is reviewed by
  at least two core developers before being included in the library.
- **caring for users’ data**. We have a functional [2]_ API and don't modify
  input arrays unless explicitly directed to do so.
- promoting **education in image processing**, with extensive pedagogical
  documentation.

Our values
**********

- We are inclusive. We continue to welcome and mentor newcomers who are
  making their first contribution.
- We are community-driven. Decisions about the API and features are driven by
  our users' requirements, not by the whims of the core team. (See
  :ref:`governance`.)
- We serve scientific applications primarily, over “consumer” image editing in
  the vein of Photoshop or GIMP. This often means prioritizing n-dimensional
  data support, and rejecting implementations of “flashy” filters that have
  little scientific value.
- We value simple, readable implementations over getting every last ounce of
  performance. Readable code that is easy to understand, for newcomers and
  maintainers alike, makes it easier to contribute new code as well as prevent
  bugs. This means that we will prefer a 20% slowdown if it reduces lines of
  code two-fold, for example.
- We value education and documentation. All functions should have NumPy-style
  docstrings [3]_, preferably with examples, as well as gallery
  examples that showcase how that function is used in a scientific application.
  Core developers take an active role in finishing documentation examples.
- We don't do magic. We use NumPy arrays instead of fancy façade objects
  [#np]_, and we prefer to educate users rather than make decisions on their
  behalf.  This does not preclude sensible defaults. [4]_

This document
*************

Much in the same way that the Zen of Python [5]_ and PEP8 guide style and
implementation details in most Python code, this guide is meant to guide
decisions about the future of scikit-image, be it in terms of code style,
whether to accept new functionality, or whether to take on new dependencies,
among other things.

References
**********

To find out more about the history of this document, please read the following:

- Original blog post [1]_
- The GitHub issue [6]_
- The image.sc forum post [7]_
- The SKIP GitHub pull request [8]_

.. rubric:: Footnotes

.. [#np] The use of NumPy arrays was the most supported of the statement's
   components, together with the points about inclusivity, mentorship, and
   documentation. We had +1s from Mark Harfouche, Royi Avital, and Greg Lee,
   among others.

Backward compatibility
----------------------

This SKIP formalizes what had been the unwritten culture of scikit-image, so it
does not raise any backward compatibility concerns.

Alternatives
------------

Two topics in the original discussion were ultimately rejected, detailed below:

Handling metadata
*****************

In my original post, I suggested that scikit-image should have some form of
metadata handling before 1.0. Among others, Mark Harfouche, Curtis Rueden, and
Dan Allan all advised that (a) maybe scikit-image doesn't *need* to handle
metadata, and can instead focus on being a robust lower-level library that
another like XArray can use to include metadata handling, and (b) anyway,
metadata support can be added later without breaking the 1.0 API. I think these
are very good points and furthermore metadata handling is super hard and I
don't mind keeping this off our plate for the moment.

Magical thinking
****************

Philipp Hanslovsky suggested [9]_ that, regarding "doing magic", it is
advisable in some contexts, and a good solution is to provide a magic layer
built on top of the non-magical one. I agree with this assessment, but, until
1.0, scikit-image should remain the non-magic layer.

Discussion
----------

See References below.

References
----------

.. [1] https://ilovesymposia.com/2018/07/13/the-road-to-scikit-image-1-0/
.. [2] https://en.wikipedia.org/wiki/Functional_programming
.. [3] https://docs.scipy.org/doc/numpy/docs/howto_document.html
.. [4] https://forum.image.sc/t/request-for-comment-road-to-scikit-image-1-0/20099/4
.. [5] https://www.python.org/dev/peps/pep-0020/
.. [6] https://github.com/scikit-image/scikit-image/issues/3263
.. [7] https://forum.image.sc/t/request-for-comment-road-to-scikit-image-1-0/20099
.. [8] https://github.com/scikit-image/scikit-image/pull/3585
.. [9] https://forum.image.sc/t/request-for-comment-road-to-scikit-image-1-0/20099/3
.. [10] CC0 1.0 Universal (CC0 1.0) Public Domain Dedication,
   https://creativecommons.org/publicdomain/zero/1.0/
.. [11] https://dancohen.org/2013/11/26/cc0-by/

Copyright
---------

This document is dedicated to the public domain with the Creative Commons CC0
license [10]_. Attribution to this source is encouraged where appropriate, as per
CC0+BY [11]_.
.. _skip_template:

==================================
SKIP X — Template and Instructions
==================================

:Author: <list of authors' real names and, optionally, email addresses>
:Status: <Draft | Active | Accepted | Deferred | Rejected | Withdrawn |
          Final | Superseded>
:Type: <Standards Track | Process>
:Created: <date created on, in yyyy-mm-dd format>
:Resolved: <date resolved, in yyyy-mm-dd format>
:Resolution: <url> (required for Accepted | Rejected | Withdrawn)
:Version effective: <version-number> (for accepted SKIPs)

Abstract
--------

The abstract should be a short description of what the SKIP will achieve.


Motivation and Scope
--------------------

This section describes the need for the proposed change. It should describe the
existing problem, who it affects, what it is trying to solve, and why. This
section should explicitly address the scope of and key requirements for the
proposed change.


Detailed description
--------------------

This section should provide a detailed description of the proposed change. It
should include examples of how the new functionality would be used, intended
use-cases, and pseudocode illustrating its use.

Related Work
------------

This section should list relevant and/or similar technologies, possibly in
other libraries. It does not need to be comprehensive, just list the major
examples of prior and relevant art.

Implementation
--------------

This section lists the major steps required to implement the SKIP. Where
possible, it should be noted where one step is dependent on another, and which
steps may be optionally omitted. Where it makes sense, each step should
include a link related pull requests as the implementation progresses.

Any pull requests or developmt branches containing work on this SKIP should
be linked to from here. (A SKIP does not need to be implemented in a single
pull request if it makes sense to implement it in discrete phases).


Backward compatibility
----------------------

This section describes the ways in which the SKIP breaks backward
compatibility.


Alternatives
------------

If there were any alternative solutions to solving the same problem, they
should be discussed here, along with a justification for the chosen
approach.


Discussion
----------

This section may just be a bullet list including links to any discussions
regarding the SKIP, but could also contain additional comments about that
discussion:

- This includes links to discussion forum threads or relevant GitHub discussions.


References and Footnotes
------------------------
All SKIPs should be declared as dedicated to the public domain with the CC0
license [1]_, as in `Copyright`, below, with attribution encouraged with CC0+BY
[2]_.

.. [1] CC0 1.0 Universal (CC0 1.0) Public Domain Dedication,
   https://creativecommons.org/publicdomain/zero/1.0/
.. [2] https://dancohen.org/2013/11/26/cc0-by/


Copyright
---------

This document is dedicated to the public domain with the Creative Commons CC0
license [1]_. Attribution to this source is encouraged where appropriate, as per
CC0+BY [2]_.
Announcement: scikits-image 0.6
===============================

We're happy to announce the 6th version of scikits-image!

Scikits-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website

  http://skimage.org

New Features
------------
- Packaged in Debian as ``python-skimage``
- Template matching
- Fast user-defined image warping
- Adaptive thresholding
- Structural similarity index
- Polygon, circle and ellipse drawing
- Peak detection
- Region properties
- TiffFile I/O plugin

... along with some bug fixes and performance tweaks.

Contributors to this release
----------------------------
- Vincent Albufera
- David Cournapeau
- Christoph Gohlke
- Emmanuelle Gouillart
- Pieter Holtzhausen
- Zachary Pincus
- Johannes Schönberger
- Tom (tangofoxtrotmike)
- James Turner
- Stefan van der Walt
- Tony S Yu
Announcement: scikits.image 0.3
===============================

After a brief (!) absence, we're back with a new and shiny version of
scikits.image, the image processing toolbox for SciPy.

This release runs under all major operating systems where
Python (>=2.6 or 3.x), NumPy and SciPy can be installed.

For more information, visit our website

  http://scikits-image.org

or the examples gallery at

   http://scikits-image.org/docs/0.3/auto_examples/

New Features
------------
- Shortest paths
- Total variation denoising
- Hough and probabilistic Hough transforms
- Radon transform with reconstruction
- Histogram of gradients
- Morphology, including watershed, connected components
- Faster homography transformations (rotations, zoom, etc.)
- Image dtype conversion routines
- Line drawing
- Better image collection handling
- Constant time median filter
- Edge detection (canny, sobel, etc.)
- IO: Freeimage, FITS, Qt and other image loaders; video support.
- SIFT feature loader
- Example data-sets

... as well as many bug fixes and minor updates.

Contributors for this release
-----------------------------
Martin Bergtholdt
Luis Pedro Coelho
Chris Colbert
Damian Eads
Dan Farmer
Emmanuelle Gouillart
Brian Holt
Pieter Holtzhausen
Thouis (Ray) Jones
Lee Kamentsky
Almar Klein
Kyle Mandli
Andreas Mueller
Neil Muller
Zachary Pincus
James Turner
Stefan van der Walt
Gael Varoquaux
Tony Yu
scikit-image 0.18.3
===================

We're happy to announce the release of scikit-image v0.18.3!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This is a small bugfix release for compatibility with Pooch 1.5 and SciPy 1.7.

Bug fixes
---------
- Only import from Pooch's public API. This resolves an import failure with
  Pooch 1.5.0. (#5531, backport of #5529)
- Do not use deprecated ``scipy.linalg.pinv2`` in ``random_walker`` when
  using the multigrid solver. (#5531, backport of #5437)

3 authors added to this release [alphabetical by first name or login]
---------------------------------------------------------------------
David Manthey
Gregory Lee
Mark Harfouche

3 reviewers added to this release [alphabetical by first name or login]
-----------------------------------------------------------------------
Gregory Lee
Juan Nunez-Iglesias
Mark Harfouche


scikit-image 0.18.2
===================

We're happy to announce the release of scikit-image v0.18.2!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This release mostly serves to add wheels for the aarch64 architecture; it also
fixes a couple of minor bugs.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org

Bug fixes
---------
- allow either SyntaxError or OSError for truncated JPG (#5315, #5334)
- Fix sphinx: role already being registered (#5319, #5335)

Development process
-------------------
- Update pyproject.toml to ensure pypy compatibility and aarch compatibility (#5326, #5328)
- Build aarch64 wheels (#5197, #5210)
- See if latest Ubuntu image fixes QEMU CPU detection issue (#5227, #5233)
- Rename `master` to `main` throughout (#5243, #5295)
- Fixup test for INSTALL_FROM_SDIST (#5283, #5296)
- Remove unecessary manual installation of packages from before_install (#5298)
- Use manylinux2010 for python 3.9+ (#5303, #5310)
- add numpy version specification on aarch for cpython 3.8 (#5374, #5375)

7 authors added to this release [alphabetical by first name or login]
---------------------------------------------------------------------
- François Boulogne
- Janakarajan Natarajan
- Juan Nunez-Iglesias
- John Lee
- Mark Harfouche
- MeeseeksMachine
- Stéfan van der Walt


9 reviewers added to this release [alphabetical by first name or login]
-----------------------------------------------------------------------
- Alexandre de Siqueira
- Gregory R. Lee
- Juan Nunez-Iglesias
- Marianne Corvellec
- Mark Harfouche
- Matti Picus
- Matthias Bussonnier
- Riadh Fezzani
- Stéfan van der Walt


scikit-image 0.18.1
===================

This is a bug fix release and contains the following two bug fixes:

- Fix indexing error for labelling in large (>2GB) arrays (#5143, #5151)
- Only use retry_if_failed with recent pooch (#5148)

See below for the new features and API changes in 0.18.0.

Announcement: scikit-image 0.18.0
=================================

We're happy to announce the release of scikit-image v0.18.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This release of scikit-image drops support for Python 3.6 in accordance with
the `NEP-29 Python and Numpy version support community standard
<https://numpy.org/neps/nep-0029-deprecation_policy.html>`_: Python 3.7 or
newer is required to run this version.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org


New Features
------------

- Add the iterative Lucas-Kanade (iLK) optical flow method (#4161)
- Add Feret diameter in region properties (#4379, #4820)
- Add functions to compute Euler number and Crofton perimeter estimation (#4380)
- Add a function to compute the Hausdorff distance (#4382)
- Added 3D support for many filters in ``skimage.filters.rank``.
- An experimental implementation of trainable pixel segmentation, aiming for
  compatibility with the scikit-learn API, has been added to
  ``skimage.future``. Try it out! (#4739)
- Add a new function ``segmentation.expand_labels`` to dilate labels while
  preventing overlap (#4795)
- It is now possible to pass extra measurement functions to
  ``measure.regionprops`` and ``regionprops_table`` (#4810)
- Add rolling ball algorithm for background subtraction (#4851)
- New sample images have been added in the ``data`` subpackage: ``data.eagle``
  (#4922), ``data.human_mitosis`` (#4939), ``data.cells3d`` (#4951), and
  ``data.vortex`` (#5041). Also note that the image for ``data.camera`` has
  been changed due to copyright issues (#4913).
- ``skimage.feature.structure_tensor`` now supports 3D (and nD) images as input
  (#5002)
- Many thresholding methods can now receive a precomputed histogram as input,
  resulting in significant speedups if multiple methods are tried on the same
  image, or if a fast histogram method is used. (#5006)
- ``measure.regionprops`` now supports multichannel intensity images (#5037)

Documentation
-------------

- Add an example to the flood fill tutorial (#4619)
- Docstring enhancements for marching cubes and find_contours (#4641)
- A new tutorial presenting a cell biology example has been added to the
  gallery (#4648). Special thanks to Pierre Poulain and Fred Bernard
  (Université de Paris and Institut Jacques Monod) for scientific review of
  this example!
- Improve register rotation example with notes and references (#4723)
- Add versionadded for new scalar type support for "scale" param in
  ``transform.AffineTransform`` (#4733)
- New tutorial on `visualizing 3D data <https://scikit-image.org/docs/dev/auto_examples/applications/plot_3d_image_processing.html>`_ (#4850)
- Add example for 3D adaptive histogram equalization (AHE) (#4658)
- Automatic formatting of docstrings for improved consistency (#4849)
- Improved docstring for ``rgb2lab`` (#4839) and ``marching_cubes`` (#4846)
- Improved docstring for ``measure.marching_cubes``, mentioning how to decimate a
  mesh using mayavi (#4846)
- Document how to contribute a gallery example. (#4857)
- Fix and improve entropy example (#4904)
- expand the benchmarking section of the developer docs (#4905)
- Improved docstring for ``util.random_noise`` (#5001)
- Improved docstrings for ``morphology.h_maxima`` and ``morphology.h_minima``
  (#4929).
- Improved docstring for ``util.img_as_int`` (#4888).
- A new example demonstrates interactive exploration of regionprops results
  using the PyData stack (pandas, seaborn) at
  <https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_regionprops.html>`_
  (#5010).
- Documentation has been added to explain
  `how to download example datasets <https://scikit-image.org/docs/dev/install.html#downloading-all-demo-datasets>`_
  which are not installed with scikit-image (#4984). Similarly, the contributor
  guide has been updated to mention how to host new datasets in a gitlab
  repository (#4892).
- The `benchmarking section of the developer documentation <https://scikit-image.org/docs/dev/contribute.html#benchmarks>`_
  has been expanded (#4905).
- Added links to the image.sc forum in example pages (#5094, #5096)
- Added missing datasets to gallery examples (#5116, #5118)
- Add farid filters in __all__, to populate the documentation (#5128, #5129)
- Proofread gallery example for rank filters. (#5126, #5136)

Improvements
------------

- float32 support for SLIC (#4683), ORB (#4684, #4697), BRIEF (#4685),
  ``pyramid_gaussian`` (#4696), Richardson-Lucy deconvolution (#4880)
- In ``skimage.restoration.richardson_lucy``, computations are now done in
  single-precision when the input image is single-precision. This can give a
  substantial performance improvement when working with single precision data.
- Richardson-Lucy deconvolution now has a ``filter_epsilon`` keyword argument
  to avoid division by very small numbers (#4823)
- Add default level parameter (max-min) / 2 in ``measure.find_contours`` (#4862)
- The performance of the SLIC superpixels algorithm
  (``skimage.segmentation.slice``) was improved for the case where a mask
  is supplied by the user (#4903). The specific superpixels produced by
  masked SLIC will not be identical to those produced by prior releases.
- ``exposure.adjust_gamma`` has been accelerated for ``uint8`` images by using
  a look-up table (LUT) (#4966).
- ``measure.label`` has been accelerated for boolean input images, by using
  ``scipy.ndimage``'s implementation for this case (#4945).
- ``util.apply_parallel`` now works with multichannel data (#4927).
- ``skimage.feature.peak_local_max`` supports now any Minkowski distance.
- We now use sparse cross-correlation to accelerate local thresholding
  functions (#4912)
- ``morphology.convex_hull_image`` now uses much less memory by checking hull
  inequalities in sequence (#5020)
- Polygon rasterization is more precise and will no longer potentially exclude
  input vertices. (#5029)
- Add data optional requirements to allow pip install scikit-image[data]
  (#5105, #5111)
- OpenMP support in MSVC (#4924, #5111)
- Restandardize handling of Multi-Image files (#2815, #5132)
- Consistent zoom boundary behavior across SciPy versions (#5131, #5133)

API Changes
-----------

- ``skimage.restoration.richardson_lucy`` returns a single-precision output
  when the input is single-precision. Prior to this release, double-precision
  was always used. (#4880)
- The default value of ``threshold_rel`` in ``skimage.feature.corner`` has
  changed from 0.1 to None, which corresponds to letting
  ``skimage.feature.peak_local_max`` decide on the default. This is currently
  equivalent to ``threshold_rel=0``.
- In ``measure.label``, the deprecated ``neighbors`` parameter has been
  removed. (#4942)
- The image returned by ``data.camera`` has changed because of copyright
  issues (#4913).

Bug fixes
---------

- A bug in ``label2rgb`` has been fixed when the input image had np.uint8
  dtype (#4661)
- Fixed incorrect implementation of ``skimage.color.separate_stains`` (#4725)
- Many bug fixes have been made in ``peak_local_max`` (#2592, #4756, #4760,
  #5047)
- Fix bug in ``random_walker`` when input labels have negative values (#4771)
- PSF flipping is now correct for Richardson-Lucy deconvolution work in >2D (#4823)
- Fix equalize_adapthist (CLAHE) for clip value 1.0 (#4828)
- For the RANSAC algorithm, improved the case where all data points are
  outliers, which was previously raising an error
  (#4844)
- An error-causing bug has been corrected for the ``bg_color`` parameter in
  ``label2rgb`` when its value is a string (#4840)
- A normalization bug was fixed in ``metrics.variation_of_information``
  (#4875)
- Euler characteristic property of ``skimage.measure.regionprops`` was erroneous
  for 3D objects, since it did not take tunnels into account. A new implementation
  based on integral geometry fixes this bug (#4380).
- In ``skimage.morphology.selem.rectangle`` the ``height`` argument
  controlled the width and the ``width`` argument controlled the height.
  They have been replaced with ``nrow`` and ``ncol``. (#4906)
- ``skimage.segmentation.flood_fill`` and ``skimage.segmentation.flood``
  now consistently handle negative values for ``seed_point``.
- Segmentation faults in ``segmentation.flood`` have been fixed (#4948, #4972)
- A segfault in ``draw.polygon`` for the case of 0-d input has been fixed
  (#4943).
- In ``registration.phase_cross_correlation``, a ``ValueError`` is raised when
  NaNs are found in the computation (as a result of NaNs in input images).
  Before this fix, an incorrect value could be returned where the input images
  had NaNs (#4886).
- Fix edge filters not respecting padding mode (#4907)
- Use v{} for version tags with pooch (#5104, #5110)
- Fix compilation error in XCode 12 (#5107, #5111)

Deprecations
------------

- The ``indices`` argument in ``skimage.feature.peak_local_max`` has been
  deprecated. Indices will always be returned. (#4752)
- In ``skimage.feature.structure_tensor``, an ``order`` argument has been
  introduced which will default to 'rc' starting in version 0.20. (#4841)
- ``skimage.feature.structure_tensor_eigvals`` has been deprecated and will be
  removed in version 0.20. Use ``skimage.feature.structure_tensor_eigenvalues``
  instead.
- The ``skimage.viewer`` subpackage and the ``skivi`` script have been
  deprecated and will be removed in version 0.20. For interactive visualization
  we recommend using dedicated tools such as `napari <https://napari.org>`_ or
  `plotly <https://plotly.com>`_. In a similar vein, the ``qt`` and ``skivi``
  plugins of ``skimage.io`` have been deprecated
  and will be removed in version 0.20. (#4941, #4954)
- In ``skimage.morphology.selem.rectangle`` the arguments ``width`` and
  ``height`` have been deprecated. Use ``nrow`` and ``ncol`` instead.
- The explicit setting ``threshold_rel=0` was removed from the Examples of the
  following docstrings: ``skimage.feature.BRIEF``,
  ``skimage.feature.corner_harris``, ``skimage.feature.corner_shi_tomasi``,
  ``skimage.feature.corner_foerstner``, ``skimage.feature.corner_fast``,
  ``skimage.feature.corner_subpix``, ``skimage.feature.corner_peaks``,
  ``skimage.feature.corner_orientations``, and
  ``skimage.feature._detect_octave``.
- In ``skimage.restoration._denoise``, the warning regarding
  ``rescale_sigma=None`` was removed.
- In ``skimage.restoration._cycle_spin``, the ``# doctest: +SKIP`` was removed.

Development process
-------------------

- Fix #3327: Add functionality for benchmark coverage (#3329)
- Release process notes have been improved. (#4228)
- ``pyproject.toml`` has been added to the sdist.
- Build and deploy dev/master documentation using GitHub Actions (#4852)
- Website now deploys itself (#4870)
- build doc on circle ci and link artifact (#4881)
- Benchmarks can now run on older scikit-image commits (#4891)
- Website analytics are tracked using plausible.io and can be visualized on
  https://plausible.io/scikit-image.org (#4893)
- Artifacts for the documentation build are now found in each pull request
  (#4881).
- Documentation source files can now be written in Markdown in addition to
  ReST, thanks to ``myst`` (#4863).
- update trove classifiers and tests for Python 3.9 + fix pytest config (#5052)
- fix Azure Pipelines, pytest config, and trove classifiers for Python 3.8 (#5054)
- Moved our testing from Travis to GitHub Actions (#5074)
- We now build our wheels on GitHub Actions on the main repo using
  cibuildwheel. Many thanks to the matplotlib and scikit-learn developers for
  paving the way for us! (#5080)
- Disable Travis-CI builds (#5099, #5111)
- Improvements to CircleCI build: no parallelization and caching) (#5097, #5119)

Other Pull Requests
-------------------

- Manage iradon input and output data type (#4298)
- random walker: Display a warning when the probability is outsite [0,1] for a given tol (#4631)
- MAINT: remove unused cython file (#4633)
- Forget legacy data dir (#4662)
- Setup longdesc markdown and switch to 0.18dev (#4663)
- Optional pooch dependency (#4666)
- Adding new default values to functions on doc/examples/segmentation/plot_ncut (#4676)
- Reintroduced convert with a strong deprecation warning (#4681)
- In release notes, better describe skimage's relationship to ecosystem (#4689)
- Perform some todo tasks for 0.18 (#4690)
- Perform todo tasks for 0.17! (#4691)
- suppressing warnings from gallery examples (#4692)
- release notes for 0.17.2 (#4702)
- Fix gallery example mentioning deprecated argument (#4706)
- Specify the encoding of files opened in the setup phase (#4713)
- Remove duplicate fused type definition (#4724)
- Blacklist cython version 0.29.18 (#4730)
- Fix CI failures related to conversion of np.floating to dtype (#4731)
- Fix Ci failures related to array ragged input numpy deprecation (#4735)
- Unwrap decorators before resolving link to source (sphinx.ext.linkcode) (#4740)
- Fix plotting error in j-invariant denoising tutorial (#4744)
- Highlight all source lines with HTML doc "source" links (sphinx.ext.linkcode) (#4746)
- Turn checklist boxes into bullet points inside the pull request template (#4747)
- Deprecate (min_distance < 1) and (footprint.size < 2) in peak_local_max (#4753)
- forbid dask 2.17.0 to fix CI (#4758)
- try to fix ci which is broken because of pyqt5 last version (#4788)
- Remove unused variable in j invariant docs (#4792)
- include all md files in manifest.in (#4793)
- Remove additional "::" to make plot directive work. (#4798)
- Use optipng to compress images/thumbnails in our gallery (#4800)
- Fix runtime warning in blob.py (#4803)
- Add TODO task for sphinx-gallery>=0.9.0 to remove enforced thumbnail_size (#4804)
- Change SSIM code example to use real MSE (#4807)
- Let biomed example load image data with Pooch. (#4809)
- Tweak threshold_otsu error checking - closes #4811 (#4812)
- Ensure assert messages from Cython rank filters are informative (#4815)
- Simplify equivalent_diameter function (#4819)
- DOC: update subpackage descriptions (#4825)
- style: be explicit when stacking arrays (#4826)
- MAINT: import Iterable from collections.abc (Python 3.9 compatibility) (#4834)
- Silence several warnings in the test suite (#4837)
- Silence a few RuntimeWarnings in the test suite (#4838)
- handle color string mapping correctly (#4840)
- DOC: Autoformat docstrings in ``io.*.py`` (#4845)
- Update min req for pillow due to CVE-2020-10379 and co. (#4861)
- DOC: First pass at format conversion, rst -> myst (#4863)
- Fixed typo in comment (#4867)
- Alternative wording for install guide PR #4750 (#4871)
- DOC: Clarify condition on unique vertices returned by marching cubes (#4872)
- Remove unmaintained wiki page link in contributor guidelines (#4873)
- new matomo config (#4879)
- Fix Incorrect documentation for skimage.util.img_as_int Issue (#4888)
- Minor edit for proper doc rendering (#4897)
- Changelog back-log (#4898)
- minor refactoring in phase_cross_correlation (#4901)
- Fix draw.circle/disk deprecation message, fixes #4884 (#4908)
- Add versionchanged tag for new opt param in measure.find_contours() (#4909)
- Declare build dependencies (#4920)
- Replace words with racial connotations (#4921)
- Fixes to apply_parallel for functions working with multichannel data (#4927)
- Improve description of h_maxima and h_minima functions (#4928) (#4929)
- CI: Skip doc build for PYTHONOPTIMIZE=2 (#4930)
- MAINT: Remove custom fused type in skimage/morphology/_max_tree.pyx (#4931)
- MAINT: remove numpydoc option, issue fixed in numpydoc 1.0 (#4932)
- modify development version string to allow use with NumpyVersion (#4947)
- CI: Add verbose option to avoid travis timeout for OSX install script  (#4956)
- Fix CI: ban sphinx-gallery 0.8.0 (#4960)
- Alias for data.chelsea: data.cat() (#4962)
- Fix typo. (#4963)
- CI: Use Travis wait improved to avoid timeout for OSX builds (#4965)
- Small enhancement in "Contour finding" example: Removed unused variable n (#4967)
- MAINT: remove unused imports (#4968)
- MAINT: Remove conditional import on networkx (#4970)
- forbid latest version of pyqt (#4973)
- Remove warnings/explicit settings on feature, restoration (#4974)
- Docstring improvements for label and regionprops_label (#4983)
- try to fix timeout problem with circleci (#4986)
- improve Euler number example (#4989)
- [website] Standardize Documentation index page. (#4990)
- Proofread INSTALL file. (#4991)
- Catch leftover typos in INSTALL file. (#4992)
- Let tifffile.imread handle additional keyword arguments (#4997)
- Update docstring for random_noise function (#5001)
- Update sphinx mapping for sklearn and numpy (#5003)
- Update docstring slic superpixels (#5014)
- Bump numpy versions to match scipy (kinda) (#5016)
- Fix usage of numpy.pad for old versions of numpy (#5017)
- [MRG] Update documentation to new data.camera() (#5018)
- bumped plotly requirement for docs (#5021)
- Fix IndexError when calling hough_line_peaks with too few angles (#5024)
- Code simplification after latest numpy bump (#5027)
- Fixes broken link to CODE_OF_CONDUCT.md (#5030)
- Specify whether core dev should merge right after second approving review. (#5040)
- Update pytest configuration to include ``test_`` functions (#5044)
- MAINT Build fix for pyodide (#5059)
- reduce OSX build time so that Travis is happy (#5067)
- DOC: document the normalized kernel in prewitt_h, prewitt_v (#5076)
- Some minor tweaks to CI (#5079)
- removed usage of numpy's private functions from util.arraycrop (#5081)
- peak_local_max: remove deprecated `indices` argument from examples (#5082)
- Replace np.bool, np.float, and np.int with bool, float, and int (#5103, #5108)
- change plausible script to track outbound links (#5115, #5123)
- Remove Python 3.6 support (#5117, #5125)
- Optimize ensure_spacing (#5062, #5135)


52 authors added to this release [alphabetical by first name or login]
----------------------------------------------------------------------

A warm thank you to all contributors who added to this release. A fraction of contributors were first-time contributors to open source and a much larger fraction first-time contributors to scikit-image. It's a great feeling for maintainers to welcome new contributors, and the diversity of scikit-image contributors is surely a big strength of the package.

- Abhishek Arya
- Abhishek Patil
- Alexandre de Siqueira
- Ben Nathanson
- Cameron Blocker
- Chris Roat
- Christoph Gohlke
- Clement Ng
- Corey Harris
- David McMahon
- David Mellert
- Devi Sandeep
- Egor Panfilov
- Emmanuelle Gouillart
- François Boulogne
- Genevieve Buckley
- Gregory R. Lee
- Harry Kwon
- iofall (cedarfall)
- Jan Funke
- Juan Nunez-Iglesias
- Julian Gilbey
- Julien Jerphanion
- kalpana
- kolibril13 (kolibril13)
- Kushaan Gupta
- Lars Grüter
- Marianne Corvellec
- Mark Harfouche
- Marvin Albert
- Matthias Bussonnier
- Max Frei
- Nathan
- neeraj3029 (neeraj3029)
- Nick
- notmatthancock (matt)
- OGordon100 (OGordon100)
- Owen Solberg
- Riadh Fezzani
- Robert Haase
- Roman Yurchak
- Ronak Sharma
- Ross Barnowski
- Ruby Werman
- ryanlu41 (ryanlu41)
- Sebastian Wallkötter
- Shyam Saladi
- Stefan van der Walt
- Terence Honles
- Volker Hilsenstein
- Wendy Mak
- Yogendra Sharma

41 reviewers added to this release [alphabetical by first name or login]
------------------------------------------------------------------------

- Abhishek Arya
- Abhishek Patil
- Alexandre de Siqueira
- Ben Nathanson
- Chris Roat
- Clement Ng
- Corey Harris
- Cris Luengo
- David Mellert
- Egor Panfilov
- Emmanuelle Gouillart
- François Boulogne
- Gregory R. Lee
- Harry Kwon
- Jan Funke
- Juan Nunez-Iglesias
- Julien Jerphanion
- kalpana
- Kushaan Gupta
- Lars Grüter
- Marianne Corvellec
- Mark Harfouche
- Marvin Albert
- neeraj3029
- Nick
- OGordon100
- Riadh Fezzani
- Robert Haase
- Ross Barnowski
- Ruby Werman
- ryanlu41
- Scott Trinkle
- Sebastian Wallkötter
- Stanley_Wang
- Stefan van der Walt
- Steven Brown
- Stuart Mumford
- Terence Honles
- Volker Hilsenstein
- Wendy Mak
Announcement: scikits-image 0.4
===============================

We're happy to announce the 0.4 release of scikits-image, an image processing
toolbox for SciPy.

Please visit our examples gallery to see what we've been up to:

   http://scikits-image.org/docs/0.4/auto_examples/

Note that, in this release, we renamed the module from ``scikits.image`` to
``skimage``, to work around name space conflicts with other scikits (similarly,
the machine learning scikit is now imported as ``sklearn``).

A big shout-out also to everyone currently at SciPy India; have fun, and
remember to join the scikits-image sprint!

This release runs under all major operating systems where Python (>=2.6 or
3.x), NumPy and SciPy can be installed.

For more information, visit our website

  http://scikits-image.org

New Features
------------
- Module rename from ``scikits.image`` to ``skimage``
- Contour finding
- Grey-level co-occurrence matrices
- Skeletonization and medial axis transform
- Convex hull images
- New test data sets
- GDAL I/O plugin

... as well as some bug fixes.

Contributors to this release
----------------------------
* Andreas Mueller
* Christopher Gohlke
* Emmanuelle Gouillart
* Neil Yager
* Nelle Varoquaux
* Riaan van den Dool
* Stefan van der Walt
* Thouis (Ray) Jones
* Tony S Yu
* Zachary Pincus
Announcement: scikits-image 0.5
===============================

We're happy to announce the 0.5 release of scikits-image, our image processing
toolbox for SciPy.

For more information, please visit our website

  http://scikits-image.org

New Features
------------
- Consistent intensity rescaling and improved range conversion.
- Random walker segmentation.
- Harris corner detection.
- Otsu thresholding.
- Block views, window views and montage.
- Plugin for Christoph Gohlke's "tifffile".
- Peak detection.
- Improved FreeImage wrappers and meta-data reading.
- 8-neighbor and background labelling.

... along with updates to the documentation and website, and a number of bug
fixes.

Contributors to this release
----------------------------
* Andreas Mueller
* Brian Holt
* Christoph Gohlke
* Emmanuelle Gouillart
* Michael Aye
* Nelle Varoquaux
* Nicolas Pinto
* Nicolas Poilvert
* Pieter Holtzhausen
* Stefan van der Walt
* Tony S Yu
* Warren Weckesser
* Zachary Pincus
Announcement: scikit-image 0.X.0
================================

We're happy to announce the release of scikit-image v0.X.0!

scikit-image is an image processing library for the scientific Python
ecosystem that includes algorithms for segmentation, geometric
transformations, feature detection, registration, color space
manipulation, analysis, filtering, morphology, and more.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org


New Features
------------



Improvements
------------



API Changes
-----------



Bugfixes
--------



Deprecations
------------



Contributors to this release
----------------------------
Announcement: scikit-image 0.14.4
=================================

We're happy to announce the release of scikit-image v0.14.4!

As a reminder, 0.14.x is the final version of scikit-image with support for
Python 2.7, and will receive critical bug fixes until Jan 1, 2020. If you
are using Python 3.5 or later, you should upgrade to scikit-image 0.15.x.

This is a bugfix release, and contains the following changes from v0.14.3:

Bug Fixes
---------
- Fix float32 support in denoise_bilateral and denoise_tv_bregman (#3937)
- Fixup test for RANSAC: don't pick duplicate samples #3901 (#3916)

Other Pull Requests
-------------------
- Backport PR #3943 on branch v0.14.x (Update the joblib link in tutorial_parallelization.rst) (#3944)

3 authors added to this release [alphabetical by first name or login]
---------------------------------------------------------------------
- Alexandre de Siqueira
- Juan Nunez-Iglesias
- Mark Harfouche


Announcement: scikit-image 0.14.3
=================================

As a reminder, 0.14.x is the final version of scikit-image with support for
Python 2.7, and will receive critical bug fixes until Jan 1, 2020. If you
are using Python 3.5 or later, you should upgrade to scikit-image 0.15.x.

This is a bugfix release, and contains the following changes from v0.14.2:

API Changes
-----------
-  morphology.local_maxima now returns a boolean array instead of uint8 (#3749,
   #3752)


Bug Fixes
---------
- _marching_cubes_lewiner_cy: mark char as signed (#3587, #3678)
- Fix potential use of NULL pointer (#3696)
- pypi: explicitly exclude Python 3.1, 3.2, and 3.3 (#3726)
- Reduce default tolerance in threshold_li (#3622) (#3781)
- Denoising functions now accept float32 images (#3449) (#3486) (#3880)

Other Pull Requests
-------------------

- BLD: pin cython's language_level (#3716)
- Build tools: Upgrade xcode to 9.4 on v0.14.x branch (#3724)
- Get rid of the requirements-parser dependency (#3534, #3727)
- Add small galleries in the API (#2940, #3728)
- Correctly ignore release notes auto-generated for docs (#3656, #3737)
- Fix qt 5.12 pinning for 0.14.x branch. (#3744, #3753)
- Minor fixes to documentation and testing infrastructure - backports #3870 and #3869 (#3881)
- Set astropy minimum requirement to 1.2 to help the CIs. (#3767, #3770)
- Avoid NumPy warning while stacking arrays. (#3768, #3771)
- Fix human readable error message on a bad build. (#3223, #3790)
- Unify LICENSE files for easier interpretation (#3791, #3792)
- Documentation formatting / compilation fixes - Backport of #3838 to v0.14.x (#3885)
- Fix build by using latest ``wheel`` package (scikit-image/scikit-image-wheels#10)


12 authors added to this release [alphabetical by first name]
-------------------------------------------------------------
- Andrew Murray
- Christoph Gohlke
- Egor Panfilov
- François Boulogne
- Johannes Schönberger
- Juan Nunez-Iglesias
- Lars Grueter
- Mark Harfouche
- Matthew Bowden
- Nehal J Wani
- Nelle Varoquaux
- Stefan van der Walt
- Thomas A Caswell

... and, as always, a special mention to Matthias Bussonnier's Meeseeks Box,
which remains invaluable for our backports.

4 committers added to this release [alphabetical by first name or login]
------------------------------------------------------------------------
- Josh Warner
- Juan Nunez-Iglesias
- Mark Harfouche
- Stefan van der Walt

5 reviewers added to this release [alphabetical by first name or login]
-----------------------------------------------------------------------
- Egor Panfilov
- François Boulogne
- Juan Nunez-Iglesias
- Mark Harfouche
- Stefan van der Walt


Announcement: scikit-image 0.14.2
=================================

This release handles an incompatibility between scikit-image and NumPy
1.16.0, released on January 13th 2019.

It contains the following changes from 0.14.1:

API changes
-----------
- ``skimage.measure.regionprops`` no longer removes singleton dimensions from
  label images (#3284). To recover the old behavior, replace
  ``regionprops(label_image)`` calls with
  ``regionprops(np.squeeze(label_image))``

Bug fixes
---------
- Address deprecation of NumPy ``_validate_lengths`` (backport of #3556)
- Correctly handle the maximum number of lines in Hough transforms
  (backport of #3514)
- Correctly implement early stopping criterion for rank kernel noise
  filter (backport of #3503)
- Fix ``skimage.measure.regionprops`` for 1x1 inputs (backport of #3284)

Enhancements
------------
- Rewrite of ``local_maxima`` with flood-fill (backport of #3022, #3447)

Build Process & Testing
-----------------------
- Dedicate a ``--pre`` build in appveyor (backport of #3222)
- Avoid Travis-CI failure regarding ``skimage.lookfor`` (backport of #3477)
- Stop using the ``pytest.fixtures`` decorator (#3558)
- Filter out DeprecationPendingWarning for matrix subclass (#3637)
- Fix matplotlib test warnings and circular import (#3632)

Contributors & Reviewers
------------------------
- François Boulogne
- Emmanuelle Gouillart
- Lars Grüter
- Mark Harfouche
- Juan Nunez-Iglesias
- Egor Panfilov
- Stefan van der Walt


Announcement: scikit-image 0.14.1
=================================

We're happy to announce the release of scikit-image v0.14.1!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This is our first release under our Long Term Support for 0.14 policy. As a
reminder, 0.14 is the last release to support Python 2.7, but it will be
updated with bug fixes and popular features until January 1st, 2020.

This release contains the following changes from 0.14.0:


Bug fixes
---------
- ``skimage.color.adapt_rgb`` was applying input functions to the wrong axis
  (#3097)
- ``CollectionViewer`` now indexes correctly (it had been broken by an update
  to NumPy indexing) (#3288)
- Handle deprecated indexing-by-list and NumPy ``matrix`` from NumPy 1.15
  (#3238, #3242, #3292)
- Fix incorrect inertia tensor calculation (#3303) (Special thanks to JP Cornil
  for reporting this bug and for their patient help with this fix)
- Fix missing comma in ``__all__`` listing of ``moments_coord_central``, so it
  and ``moments_normalized`` can now be correctly imported from the ``measure``
  namespace (#3374)
- Fix background color in ``label2rgb(..., kind='avg')`` (#3280)

Enhancements
------------
- "Reflect" mode in transforms now works fine when an image dimension has size
  1 (#3174)
- ``img_as_float`` now allows single-precision (32-bit) float arrays to pass
  through unmodified, rather than being up-converted to 64-bit (#3110, #3052,
  #3391)
- Speed up rgb2gray computation (#3187)
- The scikit-image viewer now works with different PyQt versions (#3157)
- The ``cycle_spin`` function for enhanced denoising works single-threaded
  when dask is not installed now (#3218)
- scikit-image's ``io`` module will no longer inadvertently set the matplotlib
  backend when imported (#3243)
- Fix deprecated ``get`` keyword from dask in favor of ``scheduler`` (#3366)
- Add missing ``cval`` parameter to threshold_local (#3370)


API changes
-----------
- Remove deprecated ``dynamic_range`` in ``measure.compare_psnr`` (#3313)

Documentation
-------------
- Improve the documentation on data locality (#3127)
- Improve the documentation on dealing with video (#3176)
- Update broken link for Canny filter documentation (#3276)
- Fix incorrect documentation for the ``center`` parameter of
  ``skimage.transform.rotate`` (#3341)
- Fix incorrect formatting of docstring in ``measure.profile_line`` (#3236)

Build process / development
---------------------------
- Ensure Cython is 0.23.4 or newer (#3171)
- Suppress warnings during testing (#3143)
- Fix skimage.test (#3152)
- Don't upload artifacts to AppVeyor (there is no way to delete them) (#3315)
- Remove ``import *`` from the scikit-image package root (#3265)
- Allow named non-core contributors to issue MeeseeksDev commands (#3357,
  #3358)
- Add testing in Python 3.7 (#3359)
- Add license file to the binary distribution (#3322)
- ``lookfor`` is no longer defined in ``__init__.py`` but rather imported to it
  (#3162)
- Add ``pyproject.toml`` to ensure Cython is present before building (#3295)
- Add explicit Python version Trove classifiers for PyPI (#3417)
- Ignore known test failures in 32-bit releases, allowing 32-bit wheel builds
  (#3434)
- Ignore failure to raise floating point warnings on certain ARM platforms
  (#3337)
- Fix tests to be compatible with PyWavelets 1.0 (#3406)

Credits
-------
Made with commits from (alphabetical by last name):

- François Boulogne
- Genevieve Buckley
- Sean Budd
- Matthias Bussonnier
- Sarkis Dallakian
- Christoph Deil
- François-Michel De Rainville
- Emmanuelle Gouillart
- Yaroslav Halchenko
- Mark Harfouche
- Jonathan Helmus
- Gregory Lee
- @Legodev
- Matt McCormick
- Juan Nunez-Iglesias
- Egor Panfilov
- Jesse Pangburn
- Johannes Schönberger
- Stefan van der Walt

Reviewed by (alphabetical by last name):

- François Boulogne
- Emmanuelle Gouillart
- Mark Harfouche
- Juan Nunez-Iglesias
- Egor Panfilov
- Stéfan van der Walt
- Josh Warner

And with the special support of [MeeseeksDev](https://github.com/MeeseeksBox),
created by Matthias Bussonnier


Announcement: scikit-image 0.14.0
=================================

We're happy to announce the release of scikit-image v0.14.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This is the last major release with official support for Python 2.7. Future
releases will be developed using Python 3-only syntax.

However, 0.14 is a long-term support (LTS) release and will receive bug fixes
and backported features deemed important (by community demand) until January
1st 2020 (end of maintenance for Python 2.7; see PEP 373 for details).

For more information, examples, and documentation, please visit our website:

http://scikit-image.org


New Features
------------
- Lookfor function to search across the library: ``skimage.lookfor``. (#2713)
- nD support for ``skimage.transform.rescale``, ``skimage.transform.resize``,
  and ``skimage.transform.pyramid_*`` transforms. (#1522)
- Chan-Vese segmentation algorithm. (#1957)
- Manual segmentation with matplotlib for fast data annotation:
  ``skimage.future.manual_polygon_segmentation``,
  ``skimage.future.manual_lasso_segmentation``. (#2584)
- Hysteresis thresholding:
  ``skimage.filters.apply_hysteresis_threshold``. (#2665)
- Segmentation with morphological snakes:
  ``skimage.segmentation.morphological_chan_vese`` (2D),
  ``skimage.segmentation.morphological_geodesic_active_contour`` (2D and 3D). (#2791)
- nD support for image moments: ``skimage.measure.moments_central``,
  ``skimage.measure.moments_central``, ``skimage.measure.moments_normalized``,
  ``skimage.measure.moments_hu``. This change leads to 3D/nD compatibility for
  many regionprops. (#2603)
- Image moments from coordinate input: ``skimage.measure.moments_coords``,
  ``skimage.measure.moments_coords_central``. (#2859)
- Added 3D support to ``blob_dog`` and ``blob_log``. (#2854)
- Inertia tensor and its eigenvalues can now be computed outside of
  regionprops; available in ``skimage.measure.inertia_tensor``. (#2603)
- Cycle-spinning function for approximating shift-invariance by averaging
  results from a series of spatial shifts:
  ``skimage.restoration.cycle_spin``. (#2647)
- Haar-like feature: ``skimage.feature.haar_like_feature``,
  ``skimage.feature.haar_like_feature_coord``,
  ``skimage.feature.draw_haar_like_feature``. (#2848)
- Data generation with random_shapes function:
  ``skimage.draw.random_shapes``. (#2773)
- Subset of LFW (Labeled Faces in the Wild) database:
  ``skimage.data.cbcl_face_database``. (#2905)
- Fully reworked montage function (now with a better padding behavior):
  ``skimage.util.montage``. (#2626)
- YDbDr colorspace conversion routines: ``skimage.color.rgb2ydbdr``,
  ``skimage.color.ydbdr2rgb``. (#3018)


Improvements
------------
- ``VisuShrink`` method for ``skimage.restoration.denoise_wavelet``. (#2470)
- New ``max_ratio`` parameter for ``skimage.feature.match_descriptors``. (#2472)
- ``skimage.transform.resize`` and ``skimage.transform.rescale`` have a new
  ``anti_aliasing`` option to avoid aliasing artifacts when down-sampling
  images. (#2802)
- Support for multichannel images for ``skimage.feature.hog``. (#2870)
- Non-local means denoising (``skimage.restoration.denoise_nl_means``) has
  a new optional parameter, ``sigma``, that can be used to specify the noise
  standard deviation. This enables noise-robust patch distance estimation. (#2890)
- Mixed dtypes support for ``skimage.measure.compare_ssim``,
  ``skimage.measure.compare_psnr``, etc. (#2893)
- New ``alignment`` parameter in ``skimage.feature.plot_matches``. (#2955)
- New ``seed`` parameter in ``skimage.transform.probabilistic_hough_line``. (#2960)
- Various performance improvements. (#2821, #2878, #2967, #3035, #3056, #3100)


Bugfixes
--------
- Fixed ``skimage.measure.regionprops.bbox_area`` returning incorrect value. (#2837)
- Changed gradient and L2-Hys norm computation in ``skimage.feature.hog``
  to closely follow the paper. (#2864)
- Fixed ``skimage.color.convert_colorspace`` not working for YCbCr, YPbPr. (#2780)
- Fixed incorrect composition of projective transformation with inverse transformation. (#2826)
- Fixed bug in random walker appearing when seed pixels are isolated inside pruned zones. (#2946)
- Fixed ``rescale`` not working properly with different rescale factors in multichannel case. (#2959)
- Fixed float and integer dtype support in ``skimage.util.invert``. (#3030)
- Fixed ``skimage.measure.find_contours`` raising StopIteration on Python 3.7. (#3038)
- Fixed platform-specific issues appearing in Windows and/or 32-bit environments. (#2867, #3033)


API Changes
-----------
- ``skimage.util.montage.`` namespace has been removed, and
  ``skimage.util.montage.montage2d`` function is now available as
  ``skimage.util.montage2d``.
- ``skimage.morphology.binary_erosion`` now uses ``True`` as border
  value, and is now consistent with ``skimage.morphology.erosion``.


Deprecations
------------
- ``freeimage`` plugin has been removed from ``skimage.io``.
- ``skimage.util.montage2d`` is deprecated and will be removed in 0.15.
  Use ``skimage.util.montage`` function instead.
- ``skimage.novice`` is deprecated and will be removed in 0.16.
- ``skimage.transform.resize`` and ``skimage.transform.rescale`` have a new
  ``anti_aliasing`` option that avoids aliasing artifacts when down-sampling
  images. This option will be enabled by default in 0.15.
- ``regionprops`` will use row-column coordinates in 0.16. You can start
  using them now with ``regionprops(..., coordinates='rc')``. You can silence
  warning messages, and retain the old behavior, with
  ``regionprops(..., coordinates='xy')``. However, that option will go away
  in 0.16 and result in an error. This change has a number of consequences.
  Specifically, the "orientation" region property will measure the
  anticlockwise angle from a *vertical* line, i.e. from the vector (1, 0) in
  row-column coordinates.
- ``skimage.morphology.remove_small_holes`` ``min_size`` argument is deprecated
  and will be removed in 0.16. Use ``area_threshold`` instead.


Contributors to this release
----------------------------

- Alvin
- Norman Barker
- Brad Bazemore
- Leonid Bloch
- Benedikt Boecking
- Jirka Borovec
- François Boulogne
- Larry Bradley
- Robert Bradshaw
- Matthew Brett
- Floris van Breugel
- Alex Chum
- Yannick Copin
- Nethanel Elzas
- Kira Evans
- Christoph Gohlke
- GGoussar
- Jens Glaser
- Peter Goldsborough
- Emmanuelle Gouillart
- Ben Hadfield
- Mark Harfouche
- Scott Heatwole
- Gregory R. Lee
- Guillaume Lemaitre
- Theodore Lindsay
- Kevin Mader
- Jarrod Millman
- Vinicius Monego
- Pradyumna Narayana
- Juan Nunez-Iglesias
- Kesavan PS
- Egor Panfilov
- Oleksandr Pavlyk
- Justin Pinkney
- Robert Pollak
- Jonathan Reich
- Émile Robitaille
- Rose Zhao
- Alex Rothberg
- Arka Sadhu
- Max Schambach
- Johannes Schönberger
- Sourav Singh
- Kesavan Subburam
- Matt Swain
- Saurav R. Tuladhar
- Nelle Varoquaux
- Viraj
- David Volgyes
- Stefan van der Walt
- Thomas Walter
- Scott Warchal
- Josh Warner
- Nicholas Weir
- Sera Yang
- Chiang, Yi-Yo
- corrado9999
- ed1d1a8d
- eepaillard
- leaprovenzano
- mikigom
- mrastgoo
- mutterer
- pmneila
- timhok
- zhongzyd


We'd also like to thank all the people who contributed their time to perform the reviews:

- Leonid Bloch
- Jirka Borovec
- François Boulogne
- Matthew Brett
- Thomas A Caswell
- Kira Evans
- Peter Goldsborough
- Emmanuelle Gouillart
- Almar Klein
- Gregory R. Lee
- Joan Massich
- Juan Nunez-Iglesias
- Faraz Oloumi
- Daniil Pakhomov
- Egor Panfilov
- Dan Schult
- Johannes Schönberger
- Steven Silvester
- Alexandre de Siqueira
- Nelle Varoquaux
- Stefan van der Walt
- Josh Warner
- Eric Wieser


Full list of changes
--------------------
This release is the result of 14 months of work.
It contains the following 186 merged pull requests by 67 committers:

- n-dimensional rescale, resize, and pyramid transforms (#1522)
- Segmentation: Implementation of a simple Chan-Vese Algorithm (#1957)
- JPEG quality argument in imsave (#2063)
- improve geometric models fitting (line, circle) using LSM (#2433)
- Improve input parameter handling in `_sift_read` (#2452)
- Remove broken test in `_shared/tests/test_interpolation.py` (#2454)
- [MRG] Pytest migration (#2468)
- Add VisuShrink method for `denoise_wavelet` (#2470)
- Ratio test for descriptor matching (#2472)
- Make HOG visualization use midpoints of orientation bins (#2525)
- DOC: Add example for rescaling/resizing/downscaling (#2560)
- Gallery random walker: Rescale image range to -1, 1 (#2575)
- Update conditional requirement for PySide (#2578)
- Add configuration file for `pep8_speaks` (#2579)
- Manual segmentation tool with matplotlib (#2584)
- Website updates (documentation build) (#2585)
- Update the release process notes (#2593)
- Defer matplotlib imports (#2596)
- Spelling: replaces colour by color (#2598)
- Add nD support to image moments computation (#2603)
- Set xlim and ylim in rescale gallery example (#2606)
- Reduce runtime of local_maxima gallery example (#2608)
- MAINT _shared.testing now contains pytest's useful functions (#2614)
- error message misspelled, integral to integer (#2615)
- Respect standard notations for images in functions arguments (#2617)
- MAINT: remove unused argument in private inpainting function (#2618)
- MAINT: some minor edits on Chan Vese segmentation (#2619)
- Fix UserWarning: Unknown section Example (#2620)
- Eliminate some TODOs for 0.14 (#2621)
- Clean up and fix bug in ssim tests (#2622)
- Add padding_width to montage2d and add montage_rgb (#2626)
- Add tests covering erroneous input to morphology.watershed (#2631)
- Fix name of code coverage tool (#2638)
- MAINT: Remove undefined attributes in skimage.filters (#2643)
- Improve the support for 1D images in `color.gray2rgb`  (#2645)
- ENH: add cycle spinning routine (#2647)
- as_gray replaces as_grey in imread() and load() (#2652)
- Fix AppVeyor pytest execution (#2658)
- More TODOs for 0.14 (#2659)
- pin sphinx to <1.6 (#2662)
- MAINT: use relative imports instead of absolute ones (#2664)
- Add hysteresis thresholding function (#2665)
- Improve hysteresis docstring (#2669)
- Add helper functions img_as_float32 and img_as_float64 (#2673)
- Remove unnecessary assignment in pxd file. (#2683)
- Unused var and function call in documentation example (#2684)
- Make `imshow_collection` to plot images on a grid of convenient aspect ratio (#2689)
- Fix typo in Chan-Vese docstrings (#2692)
- Fix data type error with marching_cubes_lewiner(allow_degenerate=False) (#2694)
- Add handling for uniform arrays when finding local extrema. (#2699)
- Avoid unnecessary copies in skimage.morphology.label (#2701)
- Deprecate `visualise` in favor of `visualize` in `skimage.feature.hog` (#2705)
- Remove alpha channel when saving to jpg format (#2706)
- Tweak in-place installation instructions (#2712)
- Add `skimage.lookfor` function (#2713)
- Speedup image dtype conversion by switching to `asarray` (#2715)
- MAINT reorganizing CI-related scripts (#2718)
- added rect function to draw module (#2719)
- Remove duplicate parameter in `skimage.io.imread` docstring (#2725)
- Add support for 1D arrays for grey erosion (#2727)
- Build with Xcode 9 beta 3, MacOS 10.12 (#2730)
- Travis docs one platform (#2732)
- Install documentation build requirements on Travis-CI (#2737)
- Add reference papers for `restoration.inpaint_biharmonic` (#2738)
- Completely remove `freeimage` plugin from `skimage.io` (#2744)
- Implementation and test fix for shannon_entropy calculation. (#2749)
- Minor cleanup (#2750)
- Add notes on testing to CONTRIBUTING (#2751)
- Update OSX install script (#2752)
- fix bug in horizontal seam_carve and seam_carve test. issue :#2545 (#2754)
- Recommend merging instead of rebasing, to lower contribution barrier (#2757)
- updated second link, first link still has paywall (#2768)
- DOC: set_color docstring, in-place said explicitly (#2771)
- Add module for generating random, labeled shapes (#2773)
- Ignore known failures (#2774)
- Update testdoc (#2775)
- Remove bento support (#2776)
- AppVeyor supports dot-file-style (#2779)
- Fix bug in `color.convert_colorspace` for YCbCr, YPbPr (#2780)
- Reorganizing requirements (#2781)
- WIP: Deal with long running command on travis (#2782)
- Deprecate the novice module (#2742) (#2784)
- Document mentioning deprecations in the release notes (#2785)
- [WIP] FIX Swirl center coordinates are reversed (#2790)
- Implementation of the Morphological Snakes (#2791)
- Merge TASKS.txt with CONTRIBUTING.txt (#2800)
- Add Gaussian filter-based antialiasing to resize (#2802)
- Add morphological snakes to release notes (#2803)
- Return empty array if hough_line_peaks detects nothing (#2805)
- Add W503 to pep8speaks ignore. (#2816)
- Slice PIL palette correctly using extreme image value. (#2818)
- Move INSTALL to top-level (#2819)
- Make simple watershed fast again (#2821)
- The gallery now points to the stable docs (#2822)
- Adapt AppVeyor to use Python.org dist, and remove install script (#2823)
- Remove pytest yield (#2824)
- Bug fix in projective transformation composition with inverse transformation (#2826)
- FIX: add estimate_sigma to __all__ in restoration module (#2829)
- Switch from LaTeX to MathJax in doc build (#2832)
- Docstring fixes for better formula formatting (#2834)
- Fix regionprops.bbox_area bug (#2837)
- MAINT: add Python 3.6 to appveyor, small edits (#2840)
- Allow convex area calculation in 3D for regionprops (#2847)
- [MRG] DOC fix documentation build (#2851)
- Change default args from list to tuple in `feature.draw_multiblock_lbp` (#2852)
- Add 3D support to `blob_dog` and `blob_log` (#2854)
- Update compare_nrmse docstring (#2855)
- Fix link order in example (#2858)
- Add Computation of Image Moments to Coordinates (#2859)
- Revert gradient formula, modify the deprecation warning, and fix L2-Hys norm in `skimage.feature.hog` (#2864)
- OverflowError: Python int too large to convert to C long on win-amd64-py2.7 (#2867)
- Fix `skimage.measure.centroid` and add test coverage (#2869)
- Add multichannel support to `feature.hog` (#2870)
- Remove scipy version check in `active_contour` (#2871)
- Update DOI reference in `measure.compare_ssim` (#2872)
- Fix randomness and expected ranges for RGB in `test_random_shapes`. (#2877)
- Nl means fixes for large datasets (#2878)
- Make `test_random_shapes` use internally shipped testing tools (#2879)
- DOC: Update docstring for is_low_constrast to match function signature (#2883)
- Update URL in RAG docstring (#2885)
- Fix spelling typo in NL means docstring (#2887)
- noise-robust patch distance estimation for non-local means (#2890)
- Allow mixed dtypes in compare_ssim, compare_psnr, etc. (#2893)
- EHN add Haar-like feature (#2896)
- Add CBCL face database subset to `skimage.data` (#2897)
- EXA example for haar like features (#2898)
- Install documentation dependencies on all builds (#2900)
- Improve LineModelND doc strings (#2903)
- Add a subset of LFW dataset to `skimage.data` (#2905)
- Update default parameter values in the docstring of `skimage.restoration.unsupervised_wiener` (#2906)
- Revert "Add CBCL face database subset to `skimage.data`" (#2907)
- remove unused parameter 'n_segments' in `_enforce_label_connectivity_cython()` (#2908)
- Update six version to make pytest_cov work (#2909)
- Fix typos in `draw._random_shapes._generate_triangle_mask` docstring (#2914)
- do not assume 3 channels during non-local means denoising (#2922)
- add missing cdef in _integral_image_3d (non-local means) (#2923)
- Replace `morphology.remove_small_holes` argument `min_size` with `area_threshold` (#2924)
- Ensure warning to provide bool array is warranted (#2930)
- Remove copyright notice with permission of the author (Thomas Lewiner) (#2932)
- Fix link to Windows binaries in README. (#2934)
- Handle NumPy 1.14 API changes (#2935)
- Specify `gradient` parameter docstring in `compare_ssim` (#2937)
- Fixed broken link on LBP documentation (#2941)
- Corrected bug related to border value of morphology.binary_erosion (#2945)
- Correct bug in random walker when seed pixels are isolated inside pruned zones (#2946)
- Fix Cython compilation warnings in NL Means and Watershed (#2947)
- Add `alignment` parameter to `feature.plot_matches` (#2955)
- Raise warning when attempting to save boolean image (#2957)
- Allow different rescale factors in multichannel warp (#2959)
- Add seed parameter to probabilistic_hough_line (#2960)
- Minor style fixes for #2946 (#2961)
- Build on fewer AppVeyor platforms to avoid timeout (#2962)
- Watershed segmentation: make usable for large arrays (#2967)
- Mark data_range as being a float (#2971)
- Use correct NumPy version comparison in pytest configuration (#2975)
- Handle matplotlib 2.2 pre-release deprecations (#2977)
- Bugfix LineModelND.residuals does not use the optional parameter `params` (#2979)
- Return empty list on flat images with hough_ellipse #2820 (#2996)
- Add release notes for 0.13.1 (#2999)
- MAINT: PIL removed saving RGBA images as jpeg files (#3004)
- Ensure stdev is always nonnegative in _mean_std (#3008)
- Add citation information to README (#3013)
- Add YDbDr colorspace conversion routines (#3018)
- Minor style and documentation updates for #2859 (#3023)
- `draw.random_shapes` API improvements (#3029)
- Type dependent inversion (#3030)
- Fix ValueError: Buffer dtype mismatch, expected 'int64_t' but got 'int' on win_amd64 (#3033)
- Replace pow function calls in Cython modules to fix performance issues on Windows (#3035)
- Add __pycache__ and .cache to .gitignore. (#3037)
- Fix RuntimeError: generator raised StopIteration on Python 3.7 (#3038)
- Fix invert tests (#3039)
- Fix examples not displaying figures (#3040)
- Correct reference for the coins sample image (#3042)
- Switch to basis numpy int dtypes in dtype_range (#3050)
- speedup img_as_float by making division multiplication and avoiding unnecessary allocation (#3056)
- For sparse CG solver, provide atol=0 keyword for SciPy >= 1.1 (#3063)
- Update dependencies and deprecations to fix Travis builds (#3072)
- Sanitizing marching_cubes_lewiner spacing input argument (#3074)
- Allow convex_hull_image on empty images (#3076)
- v0.13.x: Backport NumPy 1.14 compatibility (#3085)
- Force Appveyor to fail on failed tests (#3093)
- Add `threshold_local` to `filters` module namespace (#3096)
- Replace grey by gray where no deprecation is needed (#3098)
- Optimize _probabilistic_hough_line function (#3100)
- Rebuild docs upon deploy to ensure Javascript is generated (#3104)
- Fix random gallery script generation (#3106)
Announcement: scikit-image 0.10.0
=================================

We're happy to announce the release of scikit-image v0.10.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

   http://scikit-image.org


New Features
------------

In addition to many bug fixes, (speed) improvements and new examples, the 118
pull requests (1112 commits) merged for this release include the following new
or improved features (PR number in brackets):

- BRIEF, ORB and CENSURE features (#834)
- Blob Detection (#903)
- Phase unwrapping (#644)
- Wiener deconvolution (#800)
- IPython notebooks in examples gallery (#1000)
- Luv colorspace conversion (#798)
- Viewer overlays (#810)
- ISODATA thresholding (#859)
- A new rank filter for summation (#844)
- Faster MCP with anisotropy support (#854)
- N-d peak finding (`peak_local_max`) (#906)
- `imread_collection` support for all plugins (#862)
- Enforce SLIC superpixels connectivity and add SLIC-zero (#857, #864)
- Correct mesh orientation (for use in external visualizers) in
  marching cubes algorithm (#882)
- Loading from URL and alpha support in novice module (#916, #946)
- Equality for regionprops (#956)


API changes
-----------

The following backward-incompatible API changes were made between 0.9 and 0.10:

- Removed deprecated functions in `skimage.filter.rank.*`
- Removed deprecated parameter `epsilon` of `skimage.viewer.LineProfile`
- Removed backwards-compatability of `skimage.measure.regionprops`
- Removed {`ratio`, `sigma`} deprecation warnings of `skimage.segmentation.slic`
  and also remove explicit `sigma` parameter from doc-string example
- Changed default mode of random_walker segmentation to 'cg_mg' > 'cg' > 'bf',
  depending on which optional dependencies are available.
- Removed deprecated `out` parameter of `skimage.morphology.binary_*`
- Removed deprecated parameter `depth` in `skimage.segmentation.random_walker`
- Removed deprecated logger function in `skimage/__init__.py`
- Removed deprecated function `filter.median_filter`
- Removed deprecated `skimage.color.is_gray` and `skimage.color.is_rgb`
  functions
- Removed deprecated `skimage.segmentation.visualize_boundaries`
- Removed deprecated `skimage.morphology.greyscale_*`
- Removed deprecated `skimage.exposure.equalize`


Contributors to this release
----------------------------

This release was made possible by the collaborative efforts of many
contributors, both new and old.  They are listed in alphabetical order by
surname:

- Raphael Ackermann
- Ankit Agrawal
- Maximilian Albert
- Pietro Berkes
- Vighnesh Birodkar
- François Boulogne
- Olivier Debeir
- Christoph Deil
- Jaidev Deshpande
- Jaime Frio
- Jostein Bø Fløystad
- Neeraj Gangwar
- Christoph Gohlke
- Michael Hansen
- Almar Klein
- Jeremy Metz
- Juan Nunez-Iglesias
- François Orieux
- Guillem Palou
- Rishabh Raj
- Thomas Robitaille
- Michal Romaniuk
- Johannes L. Schönberger
- Steven Silvester
- Julian Taylor
- Gregor Thalhammer
- Matthew Trentacoste
- Siva Prasad Varma
- Guillem Palou Visa
- Stefan van der Walt
- Joshua Warner
- Tony S Yu
- radioxoma
Announcement: scikits-image 0.8.0
=================================

We're happy to announce the 8th version of scikit-image!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

    http://scikit-image.org


New Features
------------

- New rank filter package with many new functions and a very fast underlying
  local histogram algorithm, especially for large structuring elements
  `skimage.filter.rank.*`
- New function for small object removal
  `skimage.morphology.remove_small_objects`
- New circular hough transformation `skimage.transform.hough_circle`
- New function to draw circle perimeter `skimage.draw.circle_perimeter` and
  ellipse perimeter `skimage.draw.ellipse_perimeter`
- New dense DAISY feature descriptor `skimage.feature.daisy`
- New bilateral filter `skimage.filter.denoise_bilateral`
- New faster TV denoising filter based on split-Bregman algorithm
  `skimage.filter.denoise_tv_bregman`
- New linear hough peak detection `skimage.transform.hough_peaks`
- New Scharr edge detection `skimage.filter.scharr`
- New geometric image scaling as convenience function
  `skimage.transform.rescale`
- New theme for documentation and website
- Faster median filter through vectorization `skimage.filter.median_filter`
- Grayscale images supported for SLIC segmentation
- Unified peak detection with more options `skimage.feature.peak_local_max`
- `imread` can read images via URL and knows more formats `skimage.io.imread`

Additionally, this release adds lots of bug fixes, new examples, and
performance enhancements.


Contributors to this release
----------------------------

This release was only possible due to the efforts of many contributors, both
new and old.

- Adam Ginsburg
- Anders Boesen Lindbo Larsen
- Andreas Mueller
- Christoph Gohlke
- Christos Psaltis
- Colin Lea
- François Boulogne
- Jan Margeta
- Johannes Schönberger
- Josh Warner (Mac)
- Juan Nunez-Iglesias
- Luis Pedro Coelho
- Marianne Corvellec
- Matt McCormick
- Nicolas Pinto
- Olivier Debeir
- Paul Ivanov
- Sergey Karayev
- Stefan van der Walt
- Steven Silvester
- Thouis (Ray) Jones
- Tony S Yu
Announcement: scikit-image 0.17.2
=================================

We're happy to announce the release of scikit-image v0.17.2, which is a bug-fix
release.

Bug fixes
---------

- We made pooch an optional dependency, since it has been added as required
  dependency by mistake (#4666), and we fixed a bug about the path used for pooch
  to download data (#4662)
- The support of float 32 images was corrected for slic segmentation,
  ORB and BRIEF feature detectors (#4683, #4684, #4685, #4696, #4697)
- We removed deprecated arguments (#4691)
   * ``mask``, ``shift_x``, and ``shift_y`` from ``skimage.filters.median``
   * ``beta1`` and ``beta2`` from ``skimage.filters.frangi``
   * ``beta1`` and ``beta2`` from ``skimage.filters.hessian``
   * ``dtype`` from ``skimage.io.imread``
   * ``img`` from skimage.morphology.skeletonize_3d.
- Gallery examples were updated to suppress warnings and take into account new
  default values in some functions (#4692 and #4676)



6 authors added to this release [alphabetical by first name or login]
---------------------------------------------------------------------
- Alexandre de Siqueira
- Emmanuelle Gouillart
- François Boulogne
- Juan Nunez-Iglesias
- Mark Harfouche
- Riadh Fezzani



Announcement: scikit-image 0.17.1
=================================

We're happy to announce the release of scikit-image v0.17.1!


scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.


For more information, examples, and documentation, please visit our website:

https://scikit-image.org

Many thanks to the 54 authors who contributed the amazing number of 213 merged
pull requests! scikit-image is a community-based project and we are happy that
this number includes first-time contributors to scikit-image.

Special thanks for the release to the Cython team, who helped us make our code
compatible with their coming Cython 3.0 release. 

New Features
------------

- Hyperparameter calibration of denoising algorithms with
  `restoration.calibrate_denoiser` (#3824), with corresponding
  gallery example and tutorial.
- `measure.profile_line` has a new `reduce_func` parameter to accept a
  reduction operation to be computed on pixel values along the profile (#4206)
- nD windows for reducing spectral leakage when computing the FFT of
  n-dimensional images, with `filters.window` (#4252) (with new gallery example)
- Add Minkowski distance metric support to corner_peak (#4218)
- `util.map_array` was introduced to map a set of pixel values to another one
  (for example to map region labels to the size of regions in an image of
  labels) #4612 and #4646
- Masked marching cubes (#3829)
- The SLIC superpixel algorithm now accepts a mask to exclude some parts of the
  image and force the superpixel boundaries to follow the boundary of the mask
  (#3850)
- Pooch -- on the fly download of datasets from github: we introduced the
  possibility to include larger datasets in the `data` submodule, thanks to the
  `pooch` library. `data.download_all` fetches all datasets. (#3945)
- Starting with this version, our gallery examples now have links to run the
  example notebook on a binder instance. (#4543)

New doc tutorials and gallery examples have been added to the use of regionprops_table (#4348)
geometrical transformations (#4385), and the registration of rotation and
scaling with no shared center (#4515). A new section on registration has been
added to the gallery (#4575).

Improvements
------------

- scikit-image aims at being fully compatible with 3D arrays, and when possible
  with nD arrays. nD support has been added to color conversion functions
  (#4418), to the CLAHE `exposure.equalize_adapthist` algorithm (#4598) 
  and to the Sobel, Scharr, and Prewitt filters (#4347).
- Multichannel support for denoise_tv_bregman (#4446)
- The memory footprint of `segmentation.relabel_sequential` has been reduced in
  the case of labels much larger than the number of labels (#4612)
- Random ellipses are now possible in `draw.random_shapes` (#4493)
- Add border conditions to ridge filters (#4396)
- `segmentation.random_walker` new Jacobi preconditioned conjugate gradient mode
  (#4359) and minor corrections #4630
- Warn when rescaling with NaN in exposure.intensity_range (#4265)

We have also improved the consistency of several functions regarding the way
they handle data types

- Make dtype consistent in filters.rank functions (#4289)
- Fix colorconv float32 to double cast (#4296)
- Prevent radon from upcasting float32 arrays to double (#4297)
- Manage iradon_sart input and output data type (#4300)

API Changes
-----------

- When used with floating point inputs, ``denoise_wavelet`` no longer rescales
  the range of the data or clips the output to the range [0, 1] or [-1, 1].
  For non-float inputs, rescaling and clipping still occurs as in prior
  releases (although with a bugfix related to the scaling of ``sigma``).
- For 2D input, edge filters (Sobel, Scharr, Prewitt, Roberts, and Farid)
  no longer set the boundary pixels to 0 when a mask is not supplied. This was
  changed because the boundary mode for `scipy.ndimage.convolve` is now
  ``'reflect'``, which allows meaningful values at the borders for these
  filters. To retain the old behavior, pass
  ``mask=np.ones(image.shape, dtype=bool)`` (#4347)
- When ``out_range`` is a range of numbers and not a dtype in
  :func:`skimage.exposure.rescale_intensity`, the output data type will always
  be float (#4585)
- The values returned by :func:`skimage.exposure.equalize_adapthist` will be
  slightly different from previous versions due to different rounding behavior
  (#4585)
- Move masked_register_translation from feature to registration (#4503)
- Move register_translation from skimage.feature to skimage.registration (#4502)
- Move watershed from morphology to segmentation (#4443)
- Rename draw.circle() to draw.disk() (#4428)
- The forward and backward maps returned by :func:`skimage.segmentation.relabel_sequential`
  are no longer NumPy arrays, but more memory-efficient `ArrayMap` objects that behave
  the same way for mapping. See the ``relabel_sequential`` documentation for more details.
  To get NumPy arrays back, cast it as a NumPy array: ``np.asarray(forward_map)`` (#4612)


Bugfixes
--------

- ``denoise_wavelet``: For user-supplied `sigma`, if the input image gets
  rescaled via ``img_as_float``, the same scaling will be applied to `sigma` to
  preserve the relative scale of the noise estimate. To restore the old,
  behaviour, the user can manually specify ``rescale_sigma=False``.
- Fix Frangi artefacts around the image (#4343)
- Fix Negative eigenvalue in inertia_tensor_eigvals due to floating point precision (#4589)
- Fix morphology.flood for F-ordered images (#4556)
- Fix h_maxima/minima strange behaviors on floating point image input (#4496)
- Fix peak_local_max coordinates ordering (#4501)
- Sort naturally peaks coordinates of same amplitude in peak_local_max (#4582)
- Fix denoise_nl_means data type management (#4322)
- Update rescale_intensity to prevent under/overflow and produce proper output dtype (#4585)

(other small bug fixes are part of the list of other pull requests at the end)

Deprecations
------------
The minimal supported Python version by this release is 3.6.

- Parameter ``inplace`` in skimage.morphology.flood_fill has been deprecated
  in favor of ``in_place`` and will be removed in version scikit-image 0.19.0
  (#4250).
- ``skimage.segmentation.circle_level_set`` has been deprecated and will be
  removed in 0.19. Use ``skimage.segmentation.disk_level_set`` instead.
- ``skimage.draw.circle`` has been deprecated and will be removed in 0.19.
  Use ``skimage.draw.disk`` instead.
- Deprecate filter argument in iradon due to clash with python keyword (#4158)
- Deprecate marching_cubes_classic (#4287)
- Change label2rgb default background value from -1 to 0 (#4614)
- Deprecate rgb2grey and grey2rgb (#4420)
- Complete deprecation of circle in morphsnakes (#4467)
- Deprecate non RGB image conversion in rgb2gray (#4838, #4439), and deprecate
  non gray scale image conversion in gray2rgb (#4440)

The list of other pull requests is given at the end of this document, after the
list of authors and reviewers.

54 authors added to this release [alphabetical by first name or login]
----------------------------------------------------------------------

- aadideshpande (aadideshpande)
- Alexandre de Siqueira
- Asaf Kali
- Cedric
- D-Bhatta (D-Bhatta)
- Danielle
- Davis Bennett
- Dhiren Serai
- Dylan Cutler
- Egor Panfilov
- Emmanuelle Gouillart
- Eoghan O'Connell
- Eric Jelli
- Eric Perlman
- erjel (erjel)
- Evan Widloski
- François Boulogne
- Gregory R. Lee
- Hazen Babcock
- Jan Eglinger
- Joshua Batson
- Juan Nunez-Iglesias
- Justin Terry
- kalvdans (kalvdans)
- Karthikeyan Singaravelan
- Lars Grüter
- Leengit (Leengit)
- leGIT-bot (leGIT-bot)
- LGiki
- Marianne Corvellec
- Mark Harfouche
- Marvin Albert
- mellertd (Dave Mellert)
- Miguel de la Varga
- Mostafa Alaa
- Mojdeh Rastgoo (mrastgoo)
- notmatthancock (matt)
- Ole Streicher
- Riadh Fezzani
- robroooh (robroooh)
- SamirNasibli
- schneefux (schneefux)
- Scott Sievert
- Stefan van der Walt
- Talley Lambert
- Tim Head (betatim)
- Thomas A Caswell
- Timothy Sweetser
- Tony Tung
- Uwe Schmidt
- VolkerH (VolkerH)
- Xiaoyu Wu
- Yuanqin Lu
- Zaccharie Ramzi
- Zhōu Bówēi 周伯威


35 reviewers added to this release [alphabetical by first name or login]
------------------------------------------------------------------------
- Alexandre de Siqueira
- Asaf Kali
- D-Bhatta
- Egor Panfilov
- Emmanuelle Gouillart
- Eoghan O'Connell
- erjel
- François Boulogne
- Gregory R. Lee
- Hazen Babcock
- Jacob Quinn Shenker
- Jirka Borovec
- Josh Warner
- Joshua Batson
- Juan Nunez-Iglesias
- Justin Terry
- Lars Grüter
- Leengit
- leGIT-bot
- Marianne Corvellec
- Mark Harfouche
- Marvin Albert
- mellertd
- Miguel de la Varga
- Riadh Fezzani
- robroooh
- SamirNasibli
- Stefan van der Walt
- Timothy Sweetser
- Tony Tung
- Uwe Schmidt
- VolkerH
- Xiaoyu Wu
- Zhōu Bówēi 周伯威


Other Pull Requests
*******************
- [WIP] DOC changing the doc in plot_glcm (#2789)
- Document tophat in the gallery (#3609)
- More informative error message on boolean images for regionprops  (#4156)
- Refactor/fix threshold_multiotsu (#4178)
- Sort the generated API documentation alphabetically (#4208)
- Fix the random Linux build fails in travis CI (#4227)
- Initialize starting vector for `scipy.sparse.linalg.eigsh` to ensure reproducibility in graph_cut (#4251)
- Add histogram matching test (#4254)
- MAINT: use SciPy's implementation of convolution method (#4267)
- Improve CSS for SKIP rendering (#4271)
- Add toggle for prompts in docstring examples next to copybutton (#4273)
- Tight layout for glcm example in gallery (#4285)
- Forward port 0.16.2 release notes (#4290)
- Fix typo in `hog` docstring (#4302)
- pyramid functions take preserve_range kwarg (#4310)
- Create test and fix types (#4311)
- Deprecate numpy.pad wrapping (#4313)
- Clarify merge policy in core contributor guide (#4315)
- Regionprops is empty bug (#4316)
- Add check to avoid import craching (#4319)
- Fix typo in `simple_metrics` docstring (#4323)
- Make peak_local_max exclude_border independent and anisotropic (#4325)
- Fix blob_log/blob_dog and their corresponding tests (#4327)
- Add section on closing issues to core dev guide (#4328)
- Use gaussian filter output array if provided (#4329)
- Move cython pinning forward (#4330)
- Add python 3.8 to the build matrix (#4331)
- Avoid importing mathematical functions from scipy as told ;) (#4332)
- Add dtype keyword argument to block reduce and small documentation changes (#4334)
- Add explicit use of 32-bit int in fast_exp (#4338)
- Fix single precision cast to double in slic (#4339)
- Change `measure.block_reduce` to accept explicit `func_kwargs` kwd (#4341)
- Fix equalize_adapthist border artifacts (#4349)
- Make hough_circle_peaks respect min_xdistance, min_ydistance (#4350)
- Deprecate CONTRIBUTORS.txt and replace by git shortlog command (#4351)
- Add warning on pillow version if reading a MPO image (#4354)
- Minor documentation improvement in `measure.block_reduce` (#4355)
- Add example to highlight regionprops_table (#4356)
- Remove code that tries to avoid upgrading large dependencies from setup.py (#4362)
- Fix float32 promotion in cubic interpolation (#4363)
- Update to the new way of generating Sphinx search box (#4367)
- clarify register_translation example description (#4368)
- Bump scipy minimum version to 1.0.1 (#4372)
- Fixup OSX Builds by skipping building with numpy 1.18.0 (#4376)
- Bump pywavelets to 0.5.2 (#4377)
- mini-galleries for classes as well in API doc (#4381)
- gallery: Fix typo + reduce the angle to a reasonable value (#4386)
- setup: read long description from README (#4392)
- Do not depend on test execution order for success (#4393)
- _adapthist module refactoring and memory use reduction (#4395)
- Documentation fixes for transform (rescale, warp_polar) (#4401)
- DOC: specify the meaning of m in ransac formula (#4404)
- Updating link to values in core developer guide (#4405)
- Fix subtract_mean underflow correction (#4409)
- Fix hanging documentation build in Azure (#4411)
- Fix warnings regarding invalid escape sequences. (#4414)
- Fix the URLs in skimage.transform.pyramids (#4415)
- Fix profile_line interpolation errors (#4416)
- MAINT: replace circle_level_set by disk_level_set (#4421)
- Add stacklevel=2 to deprecation warnings in skimage.measure.marching_cubes (#4422)
- Deprecate rank.tophat and rank.bottomhat (#4423)
- Add gray2rgba and deprecate RGBA support in gray2rgb (#4424)
- ISSUE_TEMPLATE: add note about image.sc forum (#4429)
- Fix the link in skips.1-governance (#4432)
- Fix the dead link in skimage.feature.canny (#4433)
- Fix use_quantiles behavior in canny (#4437)
- Remove redundant checks for threshold values in Canny (#4441)
- Difference of Gaussians function (#4445)
- Fix test for denoise_tv_bregman accepting float32 and float64 as inputs (#4448)
- Standardize colon usage in docstrings (#4449)
- Bump numpy version to 1.15.1 (#4452)
- Set minimum tifffile version to fix numpy incompatibility (#4453)
- Cleanup warnings regarding denoise_wavelet (#4456)
- Address FutureWarning from numpy in subdtype check in reginoprops (#4457)
- Skip warnings in doctests for warning module (#4458)
- Skip doctests for deprecated functions rank.tophat rank.bottomhat since they emit warnings (#4459)
- Skip morphology.watershed doctest since it was moved and emits a warning (#4460)
- Use rgba2rgb directly where rgb kind is inferred (#4461)
- Cleanup corner peaks warnings (#4463)
- Fix edgecase bugs in segmentation.relabel_sequential (#4465)
- Fix deltaE cmc close colors bug (#4469)
- Fix bool array warping (#4470)
- Fix bool array profile_line (#4471)
- Fix values link in governance (#4472)
- Improving example on filters (#4479)
- reduce runtime of non-local means tests (#4480)
- Add sponsor button (#4481)
- reduced the duration of the longest tests (#4487)
- tiny improvements to haar feature examples (#4490)
- Add min version to sphinx-gallery >= 0.3.1 to work with py3.8 (#4498)
- Fix KeyError in find_contours (#4505)
- Fix bool array save with imageio plugin (#4512)
- Fixing order of elements in docstrings of skimage/color/colorconv (#4518)
- Fix exposure_adapthist return when clip_limit == 1 (#4519)
- Adding info on venv activation on Windows (#4521)
- Fix similarity transform scale (#4524)
- Added explanation in the example of `segmentation/plot_label.py` to make the background transparent (#4527)
- Add example code for generating structuring elements. (#4528)
- Block imread version 0.7.2 due to build failure (#4529)
- Maint: edits to suppress some warnings (unused imports, blank lines) (#4530)
- MNT: remove duplicate nogil specification (#4546)
- Block pillow 7.1.0, see #4548 (#4551)
- Fix binder requirements (#4555)
- Do not enforce pil plugin in skimage.data (#4560)
- Remove "backport to 0.14" in github template (#4561)
- Fix inconsistency in docstring (filters.median) (#4562)
- Disable key check for texlive in travis-mac as a temporary workaround (#4565)
- Bump Pywavelets min requirement to 1.1.1 (#4568)
- Strip backslash in sphinx 3.0.0 (#4569)
- Remove binary specification from match_descriptors docstring (#4571)
- Remove importing skimage.transform as tf (#4576)
- Add note to remove option in doc config when numpydoc will be patched (#4578)
- update task in TODO.txt (#4579)
- Rename convert to _convert, as it is a private function (#4590)
- Do not overwrite data module in plot_skeleton.py (#4591)
- [CI fix] add import_array in cython files where numpy is cimport-ed (#4592)
- Recommend cnp.import_array in contribution guide (#4593)
- Add example of natsort usage in documentation (#4599)
- Fix broken and permanently moved links (#4600)
- Fix typo in cython import_array (#4602)
- Update min required sphinx version for sphinx-copybutton (#4604)
- Clarify error message when montaging multichannel nD images and multichannel=False (#4607)
- Fix register_translation warning message (#4609)
- Add notes on deprecation warnings in marching_cube_* and gray2rgb (#4610)
- Improve loading speed of our gallery by reducing the thumbnail size (#4613)
- Fixed wrong behaviour of `exposure.rescale_intensity` for constant input. (#4615)
- Change math formatting in the docstrings (#4617)
- Add .mypy_cache to the gitignore (#4620)
- typo fixes for register rotation gallery example (#4623)
- Userguide: add a visualization chapter (#4627)
- Fix deprecation warnings due to invalid escape sequences.  (#4628)
- add docstring examples for moments_hu and centroid (#4632)
- Update pooch registry with new file location (#4635)
- Misleading "ValueError: Input array has to be either 3- or 4-dimensional" in montage (#4638)
- Fix broken link (#4639)
- AffineTransform: Allow a single value for 'scale' to apply to both sx & sy (#4642)
- Fix CI - cython 3.0a4 (#4643)
- Fix sphinx (#4644)
- Fix ArrayMap test (#4645)
- Remove copy of tifffile; install from pip (#4235)
- Refactor/move neighborhood utility functions in morphology (#4209)

Announcement: scikit-image 0.15.0
=================================

We're happy to announce the release of scikit-image v0.15.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org

0.15 is the first scikit-image release that is only compatible with Python 3.5
and above. Python 2.7 users should strongly consider upgrading to Python 3.5+,
or use the 0.14 long term support releases.


New Features
------------

- N-dimensional flood fill, with tolerance (#3245)
- Attribute operators (#2680)
- Extension of register_translation to enable subpixel precision in 3D and
  optionally disable error calculation (#2880)
- unsharp mask filtering (#2772)
- New options ``connectivity``, ``indices`` and ``allow_borders`` for
  ``skimage.morphology.local_maxima`` and ``local_minima``. (#3022)
- Image translation registration for masked data
  (``skimage.feature.masked_register_translation``) (#3334)
- Frangi (vesselness), Meijering (neuriteness), and Sato (tubeness) filters
  (#3515)
- Allow float->float conversion of any range (#3052)
- Let lower precision float arrays pass through ``img_as_float`` (#3110)
- Lazy apply_parallel (allows optimization of dask array operations) (#3121)
- Add range option for histogram. (#2479)
- Add histogram matching (#3568)


Improvements
------------

- Replace ``morphology.local_maxima`` with faster flood-fill based Cython
  version (#3022)
- ``skivi`` is now using ``qtpy`` for Qt4/Qt5/PySide/PySide2 compatibility (a
  new optional dependency).
- Performance is now monitored by
  `Airspeed Velocity <https://asv.readthedocs.io/en/stable/>`_. Benchmark
  results will appear at https://pandas.pydata.org/speed/ (#3137)
- Speed up inner loop of GLCM (#3378)
- Allow tuple to define kernel in threshold_niblack and threshold_sauvola (#3596)
- Add support for anisotropic blob detection in blob_log and blob_dog (#3690)


API Changes
-----------

- ``skimage.transform.seam_carve`` has been removed because the algorithm is
  patented. (#3751)
- Parameter ``dynamic_range`` in ``skimage.measure.compare_psnr`` has been
  removed. Use parameter ``data_range`` instead. (#3313)
- imageio is now the preferred plugin for reading and writing images. (#3126)
- imageio is now a dependency of scikit-image. (#3126)
- ``regular_grid`` now returns a tuple instead of a list for compatibility
  with numpy 1.15 (#3238)
- ``colorconv.separate_stains`` and ``colorconv.combine_stains`` now uses
  base10 instead of the natural logarithm as discussed in issue #2995. (#3146)
- Default value of ``clip_negative`` parameter in ``skimage.util.dtype_limits``
  has been set to ``False``.
- Default value of ``circle`` parameter in ``skimage.transform.radon``
  has been set to ``True``.
- Default value of ``circle`` parameter in ``skimage.transform.iradon``
  has been set to ``True``.
- Default value of ``mode`` parameter in ``skimage.transform.swirl``
  has been set to ``reflect``.
- Deprecated ``skimage.filters.threshold_adaptive`` has been removed.
  Use ``skimage.filters.threshold_local`` instead.
- Default value of ``multichannel`` parameter in
  ``skimage.restoration.denoise_bilateral`` has been set to ``False``.
- Default value of ``multichannel`` parameter in
  ``skimage.restoration.denoise_nl_means`` has been set to ``False``.
- Default value of ``mode`` parameter in ``skimage.transform.resize``
  and ``skimage.transform.rescale`` has been set to ``reflect``.
- Default value of ``anti_aliasing`` parameter in ``skimage.transform.resize``
  and ``skimage.transform.rescale`` has been set to ``True``.
- Removed the ``skimage.test`` function. This functionality can be achieved
  by calling ``pytest`` directly.
- ``morphology.local_maxima`` now returns a boolean array (#3749)


Bugfixes
--------

- Correct bright ridge detection for Frangi filter (#2700)
- ``skimage.morphology.local_maxima`` and ``skimage.morphology.local_minima``
  no longer raise an error if any dimension of the image is smaller 3 and
  the keyword ``allow_borders`` was false.
- ``skimage.morphology.local_maxima`` and ``skimage.morphology.local_minima``
  will return a boolean array instead of an array of 0s and 1s if the
  parameter ``indices`` was false.
- When ``compare_ssim`` is used with ``gaussian_weights`` set to True, the
  boundary crop used when computing the mean structural similarity will now
  exactly match the width of the Gaussian used. The Gaussian filter window is
  also now truncated at 3.5 rather than 4.0 standard deviations to exactly match
  the original publication on the SSIM. These changes should produce only a very
  small change in the computed SSIM value. There is no change to the existing
  behavior when ``gaussian_weights`` is False. (#3802)
- erroneous use of cython wrap around (#3481)
- Speed up block reduce by providing the appropriate parameters to numpy (#3522)
- Add urllib.request again (#3766)
- Repeat pixels in reflect mode when image has dimension 1 (#3174)
- Improve Li thresholding (#3402, 3622)


Deprecations
------------

- Python 2 support has been dropped. Users should have Python >= 3.5. (#3000)
- ``skimage.util.montage2d`` has been removed. Use ``skimage.util.montage`` instead.
- ``skimage.novice`` is deprecated and will be removed in 0.16.
- ``skimage.transform.resize`` and ``skimage.transform.rescale`` option
  ``anti_aliasing`` has been enabled by default.
- ``regionprops`` will use row-column coordinates in 0.16. You can start
  using them now with ``regionprops(..., coordinates='rc')``. You can silence
  warning messages, and retain the old behavior, with
  ``regionprops(..., coordinates='xy')``. However, that option will go away
  in 0.16 and result in an error. This change has a number of consequences.
  Specifically, the "orientation" region property will measure the
  anticlockwise angle from a *vertical* line, i.e. from the vector (1, 0) in
  row-column coordinates.
- ``skimage.morphology.remove_small_holes`` ``min_size`` argument is deprecated
  and will be removed in 0.16. Use ``area_threshold`` instead.
- ``skimage.filters.median`` will change behavior in the future to have an
  identical behavior as ``scipy.ndimage.median_filter``. This behavior can be
  set already using ``behavior='ndimage'``. In 0.16, it will be the default
  behavior and removed in 0.17 as well as the parameter of the previous
  behavior (i.e., ``mask``, ``shift_x``, ``shift_y``) will be removed.


Documentation improvements
--------------------------

- Correct rotate method's center parameter doc (#3341)
- Add Sphinx copybutton (#3530)
- Add glossary to the documentation (#3626)
- Add image of retina to our data (#3748)
- Add microaneurysms() to gallery (#3765)
- Better document remove_small_objects behaviour: int vs bool (#2830)
- Linking preserve_range parameter calls to docs (#3109)
- Update the documentation regarding datalocality (#3127)
- Specify conda-forge channel for scikit-image conda install (#3189)
- Turn DOIs into web links in docstrings (#3367)
- Update documentation for regionprops (#3602)
- DOC: Improve the RANSAC gallery example (#3554)
- DOC: "feature.peak_local_max" : explanation of multiple same-intensity peaks returned by the function; added details on ``exclude_border`` parameter  (#3600)


Improvements
------------

- MNT: handle a deprecation warning for np.linspace and floats for the num parameter (#3453)
- TST: numpy empty arrays are not inherently Falsy (#3455)
-  handle warning in scipy cdist for unused parameters (#3456)
- MNT: don't use filter_warnings in test suite. (#3459)
- Add doc notes on setting up the build environment (#3472)
- Release the GIL in numerous cython functions (#3490)
- Cython touchups to use float32 and float64 (#3493)
- rank_filters: Change how the bitdepth and max_bin are computed to ensure exact warnings. (#3501)
- Rank: Optimize OTSU filter (#3504)
- Rank - Fix rank entropy and OTSU tests (#3506)
- delay importing pyplot in manual segmentation (#3533)
- Get rid of the requirements-parser dependency (#3534)
- filter warning from ``correct_mesh_orientation`` in tests (#3549)
- cloudpickle is really a doc dependency, not a core one (#3634)
- optional dependencies on pip (#3645)
- Fewer test warnings in 3.7 (#3687)
- collections.abc nit (#3692)
- Streamlined issue template (#3697)
- Tighten the PR Template (#3701)
- Use language level to 3 in cython for future compatibility (#3707)
- Update ISSUE_TEMPLATE.md with info about numpy and skimage versions (#3730)
- Use relative imports for many cython modules (#3759)
- Pass tests that don't raise floating point exceptions on arm with soft-fp (#3337)


Other improvements
------------------

- BUG: Fix greycoprops correlation always returning 1 (#2532)
- Add section on API discovery via ``skimage.lookfor`` (#2539)
- Speedup 2D warping for affine transformations (#2902)
- Credit Reviewers in Release Notes (#2927)
- Added small galleries in the API (#2940)
- Use skimage gaussian filter to avoid integer rounding artifacts (#2983)
- Remove Python 2 compatibility (#3000)
- Add ``rectangle_perimeter`` feature to ``skimage.draw`` (#3069)
- Update installation instructions to reference existing requirements specification (#3113)
- Updated release notes with pre 0.13.1 phase (#3114)
- Release guidelines update (#3115)
- Ensure we are installing with / running on Python 3 (#3119)
- Hide warnings in test_unsharp_mask (#3130)
- Process 0.15 deprecations (#3132)
- Documentation: always use dev branch javascript (#3136)
- Add initial airspeed velocity (asv) framework (#3137)
- Supress warnings for flatten during io testing (#3143)
- Recover from exceptions in filters.try_all_threshold() (#3149)
- Fix skimage.test() to run the unittests (#3152)
- skivi: Use qtpy to handle different Qt versions (#3157)
- Refactor python version checking. (#3160)
- Move data_dir to within ``data/__init__.py`` (#3161)
- Move the definition of lookfor out of __init__.py (#3162)
- Normalize the package number to PEP440 (#3163)
- Remove skimage.test as it was never used. (#3164)
- Added a message about qtpy to the INSTALL.rst (#3168)
- Regression fix: Travis should fail if tests fail (#3170)
- Set minimum cython version to ``0.23.4`` (#3171)
- Add rgba2rgb to API docs (#3175)
- Minor doc formatting fixes in video.rst (#3176)
- Decrease the verbosity of the testing (#3182)
- Speedup rgb2gray using matrix multiply (#3187)
- Add instructions for meeseeksdev to PR template (#3194)
- Remove installation instructions for video packages (#3197)
- Big image labeling fix (#3202)
- Handle dask deprecation in cycle_spin (#3205)
- Fix Qt viewer painttool indexing (#3210)
- build_versions.py is no longer hard coded. (#3211)
- Remove dtype constructor call in exposure.rescale_intensity (#3213)
- Various updates to the ASV benchmarks (#3215)
- Add a link to stack overflow on github README (#3217)
- MAINT: remove encoding information in file headers (python 3) (#3219)
- Build tools: Dedicate a --pre build in appveyor and ensure other builds don't download --pre (#3222)
- Fix the human readable error message on a bad build. (#3223)
- Respect input array type in apply_parallel by default (#3225)
- Travis cleanup pip commands (#3227)
- Add benchmarks for morphology.watershed (#3234)
- Correcte docstring formatting so that code block is displayed as code (#3236)
- Defer skimage.io import of matplotlib.pyplot until needed (#3243)
- Add benchmark for Sobel filters (#3249)
- Remove cython md5 hashing since it breaks the build process (#3254)
- Fix typo in documentation. (#3262)
- Issue 3156: skimage/__init__.py Update docstring and fix import *  (#3265)
- Object detector module (#3267)
- Do not import submodules while building (#3270)
- Add benchmark suite for canny (#3271)
- improve segmentation.felzenszwalb document #3264 (#3272)
- Update _canny.py (#3276)
- Add benchmark suite for histogram equalization (#3285)
- fix link to equalist_hist blog reference (#3287)
- .gitignore: novice: Ignore save-demo.jpg (#3289)
- Guide the user of denoise_wavelet to choose an orthogonal wavelet. (#3290)
- Remove unused lib in skimage/__init__.py (#3291)
- BUILD: Add pyproject.toml to ensure cython is present (#3295)
- Handle intersphinx and mpl deprecation warnings in docs (#3300)
- Minor PEP8 fixes (#3305)
- cython: check for presence of cpp files during install from sdist (#3311)
- appveyor: don't upload any artifacts (#3315)
- Add benchmark suite for hough_line() (#3319)
- Novice skip url test (#3320)
- Remove benchmarks from wheel (#3321)
- Add license file to the wheel (binary) distribution (#3322)
- codecov: ignore build scripts in coverage and don't comment on PRs (#3326)
- Matplotlib 2.2.3 +  PyQt5.11 (#3345)
- Allow @hmaarrfk to mention MeeseeksDev to backport. (#3357)
- Add Python 3.7 to the test matrix (#3359)
- Fix deprecated keyword from dask (#3366)
- Incompatible modes with anti-aliasing in skimage.transform.resize (#3368)
- Missing cval parameter in threshold_local (#3370)
- Avoid Sphinx 1.7.8 (#3381)
- Show our data in the gallery (#3388)
- Minor updates to grammar in numpy images page (#3389)
- assert_all_close doesn't exist, make it ``assert_array_equal`` (#3391)
- Better behavior of Gaussian filter for arrays with a large number of dimensions (#3394)
- Allow import/execution with -OO (#3398)
- Mark tests known to fail on 32bit architectures with xfail (#3399)
- Hardcode the inputs to test_ssim_grad (#3403)
- TST: make test_wavelet_denoising_levels compatible with PyWavelets 1.0 (#3406)
- Allow tifffile.py to handle I/O. (#3409)
- Add explicit Trove classifier for Python 3 (#3415)
- Fix error in contribs.py (#3418)
- MAINT: remove pyside restriction since we don't support Python 3.4 anymore (#3421)
- Build tools: simplify how MPL_DIR is obtained. (#3422)
- Build tools: Don't run tests twice in travis. (#3423)
- Build tools: Add an OSX build with optional dependencies. (#3424)
- MAINT: Reverted the changes in #3300 that broke the MINIMIUM_REQUIREMENTS tests (#3427)
- MNT: Convert links using http to https (#3428)
- MAINT: Use upstream colormaps now that matplotlib has been upgraded (#3429)
- Build tools: Make pyamg an optional dependency and remove custom logic (#3431)
- Build tools: Fix PyQt installed in minimum requirements build (#3432)
- MNT: multiprocessing should always be available since we depend on python >=2.7 (#3434)
- MAINT Use np.full instead of cst*np.ones (#3440)
- DOC: Fix LaTeX build via ``make latexpdf``  (#3441)
- Update instructions et al for releases after 0.14.1 (#3442)
- Remove code specific to python 2 (#3443)
- Fix default value of ``methods`` in ``_try_all`` to avoid exception (#3444)
- Fix morphology.local_maxima for input with any dimension < 3 (#3447)
- Use raw strings to avoid unknown escape symbol warnings (#3450)
- Speed up xyz2rgb by clipping output in place (#3451)
- MNT; handle deprecation warnings in tifffile (#3452)
- Build tools: TST: filter away novice deprecation warnings during testing (#3454)
- Build tools: don't use the pytest.fixtures decorator anymore in class fixtures  (#3458)
- Preserving the fill_value of a masked array (#3461)
- Fix VisibleDeprecationWarning from np.histogram, normed=True (#3463)
- Build Tools: DOC: Document that now PYTHONOPTMIZE build is blocked by SciPy (#3470)
- DOC: Replace broken links by webarchive equivalent links (#3471)
- FIX: making the plot_marching_cubes example visible. (#3474)
- Avoid Travis failure regarding ``skimage.lookfor`` (#3477)
- Fix Python executable for sphinx-build in docs Makefile (#3478)
- Build Tools: Block specific Cython versions (#3479)
- Fix typos (#3480)
- Add "optional" indications to docstrings (#3495)
- Rename 'mnxc' (masked normalize cross-correlation) to something more descriptive (#3497)
- Random walker bug fix: no error should be raised when there is nothing to do (#3500)
- Various minor edits for active contour (#3508)
- Fix range for uint32 dtype in user guide (#3512)
- Raise meaningful exception in warping when image is empty (#3518)
- DOC: Development installation instructions for Ubuntu are missing tkinter (#3520)
- Better gallery examples and tests for masked translation registration (#3528)
- DOC: make more docstrings compliant with our standards (#3529)
- Build tools: Remove restriction on simpleitk for python 3.7 (#3535)
- Speedup and add benchmark for ``skeletonize_3d`` (#3536)
- Update requirements/README.md on justification of matplotlib 3.0.0 in favor of #3476 (#3542)
- Doc enhancements around denoising features. (#3553)
- Use 'getconf _NPROCESSORS_ONLN' as fallback for nproc in Makefile of docs (#3563)
- Fix matplotlib set_*lim API deprecations (#3564)
- Switched from np.power to np.cbrt (#3570)
- Filtered out DeprecationPendingWarning for matrix subclass (#3572)
- Add RGB to grayscale example to gallery (#3574)
- Build tools: Refactor check_sdist so that it takes a filename as a parameter (#3579)
- Turn dask to an optional requirement (#3582)
- _marching_cubes_lewiner_cy: mark char as signed (#3587)
- Hyperlink DOIs to preferred resolver (#3589)
- Missing parameter description in ``morphology.reconstruction`` docstring #3581 (#3591)
- Update chat location (#3598)
- Remove orphan code (skimage/filters/_ctmf.pyx). (#3601)
- More explicit example title, better list rendering in plot_cycle_spinning.py (#3606)
- Add rgb to hsv example in the gallery (#3607)
- Update documentation of ``perimeter`` and add input validation (#3608)
- Additionnal mask option to clear_border (#3610)
- Set up CI with Azure Pipelines (#3612)
- [MRG] EHN: median filters will accept floating image (#3616)
- Update Travis-CI to xcode 10.1 (#3617)
- Minor tweaks to _mean_std code (#3619)
- Add explicit ordering of gallery sections (#3627)
- Delete broken links (#3628)
- Build tools: Fix test_mpl_imshow for matplotlib 2.2.3 and numpy 1.16 (#3635)
- First draft of core dev guide (#3636)
- Add more details about the home page build process (#3639)
- Ensure images resources with long querystrings can be read (#3642)
- Delay matplotlib import in skimage/future/manual_segmentation.py (#3648)
- make the low contrast check optional when saving images (#3653)
- Correctly ignore release notes auto-generated for docs (#3656)
- Remove MANIFEST file when making the 'clean' target (#3657)
- Clarify return values in _overlap docstrings in feature/blob.py (#3660)
- Contribution script: allow specification of GitHub development branch (#3661)
- Update core dev guide: deprecation, contributor guide, required experience (#3662)
- Add release notes for 0.14.2 (#3664)
- FIX gallery: Add multichannel=True to match_histogram (#3672)
- MAINT Minor code style improvements (#3673)
- Pass parameters through tifffile plugin (#3675)
- DOC unusused im3d_t in example (#3677)
- Remove wrong cast of Py_ssize_t to int (#3682)
- Build tools: allow python 3.7 to fail, but travis to continue (#3683)
- Build tools: remove pyproject.toml (#3688)
- Fix ValueError: not enough values to unpack (#3703)
- Several fixes for heap.pyx (#3704)
- Enable the faulthandler module during testing (#3708)
- Build tools: Fix Python 3.7 builds on travis (#3709)
- Replace np.einsum with np.tensordot in _upsampled_dft (#3710)
- Fix potential use of NULL pointers (#3717)
- Fix potential memory leak (#3718)
- Fix potential use of NULL pointers (#3719)
- Fix and improve core_cy.pyx (#3720)
- Build tools: Downgrade Xcode to 9.4 on master (#3723)
- Improve visual_test.py (#3732)
- Updated painttool to work with color images and properly scale labels. (#3733)
- Add image.sc forum badge to README (#3738)
- Block PyQt 5.12.0 on Travis (#3743)
- Build tools: Fix matplotlib + qt 5.12 the same way upstream does it (#3744)
- gallery: remove xx or yy  sorted directory names (#3761)
- Allow for f-contiguous 2D arrays in convex_hull_image (#3762)
- Build tools: Set astropy minimum requirement to 1.2 to help the CIs. (#3767)
- Avoid NumPy warning while stacking arrays. (#3768)
- Set CC0 for microaneurysms (#3778)
- Unify LICENSE files for easier interpretation (#3791)
- Readme: Remove expectation for future fix from matplotlib (#3794)
- Improved documentation/test in ``flood()`` (#3796)
- Use ssize_t in denoise cython (#3800)
- Removed non-existent parameter in docstring (#3803)
- Remove redundant point in draw.polygon docstring example (#3806)
- Ensure watershed auto-markers respect mask (#3809)


75 authors added to this release [alphabetical by first name or login]
----------------------------------------------------------------------

- Abhishek Arya
- Adrian Roth
- alexis-cvetkov (Alexis Cvetkov-Iliev)
- Ambrose J Carr
- Arthur Imbert
- blochl (Leonid Bloch)
- Brian Smith
- Casper da Costa-Luis
- Christian Rauch
- Christoph Deil
- Christoph Gohlke
- Constantin Pape
- David Breuer
- Egor Panfilov
- Emmanuelle Gouillart
- fivemok
- François Boulogne
- François Cokelaer
- François-Michel De Rainville
- Genevieve Buckley
- Gregory R. Lee
- Gregory Starck
- Guillaume Lemaitre
- Hugo
- jakirkham (John Kirkham)
- Jan
- Jan Eglinger
- Jathrone
- Jeremy Metz
- Jesse Pangburn
- Johannes Schönberger
- Jonathan J. Helmus
- Josh Warner
- Jotham Apaloo
- Juan Nunez-Iglesias
- Justin
- Katrin Leinweber
- Kim Newell
- Kira Evans
- Kirill Klimov
- Lars Grueter
- Laurent P. René de Cotret
- Legodev
- mamrehn
- Marcel Beining
- Mark Harfouche
- Matt McCormick
- Matthias Bussonnier
- mrastgoo
- Nehal J Wani
- Nelle Varoquaux
- Onomatopeia
- Oscar Javier Hernandez
- Page-David
- PeterJackNaylor
- PinkFloyded
- R S Nikhil Krishna
- ratijas
- Rob
- robroooh
- Roman Yurchak
- Sarkis Dallakian
- Scott Staniewicz
- Sean Budd
- shcrela
- Stefan van der Walt
- Taylor D. Scott
- Thein Oo
- Thomas Walter
- Tom Augspurger
- Tommy Löfstedt
- Tony Tung
- Vilim Štih
- yangfl
- Zhanwen "Phil" Chen


46 reviewers added to this release [alphabetical by first name or login]
------------------------------------------------------------------------

- Abhishek Arya
- Adrian Roth
- Alexandre de Siqueira
- Ambrose J Carr
- Arthur Imbert
- Brian Smith
- Christian Rauch
- Christoph Gohlke
- David Breuer
- Egor Panfilov
- Emmanuelle Gouillart
- Evan Putra Limanto
- François Boulogne
- François Cokelaer
- Gregory R. Lee
- Grégory Starck
- Guillaume Lemaitre
- Ilya Flyamer
- jakirkham
- Jarrod Millman
- Johannes Schönberger
- Josh Warner
- Jotham Apaloo
- Juan Nunez-Iglesias
- Justin
- Lars Grueter
- Laurent P. René de Cotret
- Marcel Beining
- Mark Harfouche
- Matthew Brett
- Matthew Rocklin
- Matti Picus
- mrastgoo
- Onomatopeia
- PeterJackNaylor
- Rob
- Roman Yurchak
- Scott Staniewicz
- Stefan van der Walt
- Thein Oo
- Thomas A Caswell
- Thomas Walter
- Tom Augspurger
- Tomas Kazmar
- Tommy Löfstedt
- Vilim Štih
Announcement: scikit-image 0.11.0
=================================

We're happy to announce the release of scikit-image v0.11.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

http://scikit-image.org

Highlights
----------
For this release, we merged over 200 pull requests with bug fixes,
cleanups, improved documentation and new features.  Highlights
include:

- Region Adjacency Graphs
  - Color distance RAGs (#1031)
  - Threshold Cut on RAGs (#1031)
  - Similarity RAGs (#1080)
  - Normalized Cut on RAGs (#1080)
  - RAG drawing (#1087)
  - Hierarchical merging (#1100)
- Sub-pixel shift registration (#1066)
- Non-local means denoising (#874)
- Sliding window histogram (#1127)
- More illuminants in color conversion (#1130)
- Handling of CMYK images (#1360)
- `stop_probability` for RANSAC (#1176)
- Li thresholding (#1376)
- Signed edge operators (#1240)
- Full ndarray support for `peak_local_max` (#1355)
- Improve conditioning of geometric transformations (#1319)
- Standardize handling of multi-image files (#1200)
- Ellipse structuring element (#1298)
- Multi-line drawing tool (#1065), line handle style (#1179)
- Point in polygon testing (#1123)
- Rotation around a specified center (#1168)
- Add `shape` option to drawing functions (#1222)
- Faster regionprops (#1351)
- `skimage.future` package (#1365)
- More robust I/O module (#1189)

API Changes
-----------
- The ``skimage.filter`` subpackage has been renamed to ``skimage.filters``.
- Some edge detectors returned values greater than 1--their results are now
  appropriately scaled with a factor of ``sqrt(2)``.

Contributors to this release
----------------------------
(Listed alphabetically by last name)

- Fedor Baart
- Vighnesh Birodkar
- François Boulogne
- Nelson Brown
- Alexey Buzmakov
- Julien Coste
- Phil Elson
- Adam Feuer
- Jim Fienup
- Geoffrey French
- Emmanuelle Gouillart
- Charles Harris
- Jonathan Helmus
- Alexander Iacchetta
- Ivana Kajić
- Kevin Keraudren
- Almar Klein
- Gregory R. Lee
- Jeremy Metz
- Stuart Mumford
- Damian Nadales
- Pablo Márquez Neila
- Juan Nunez-Iglesias
- Rebecca Roisin
- Jasper St. Pierre
- Jacopo Sabbatini
- Michael Sarahan
- Salvatore Scaramuzzino
- Phil Schaf
- Johannes Schönberger
- Tim Seifert
- Arve Seljebu
- Steven Silvester
- Julian Taylor
- Matěj Týč
- Alexey Umnov
- Pratap Vardhan
- Stefan van der Walt
- Joshua Warner
- Tony S Yu
Announcement: scikit-image 0.16.2
=================================

We're happy to announce the release of scikit-image v0.16.2!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

This is a bug fix release that addresses several critical issues from 0.16.1.

Bug fixes
---------
- Migrate to networkx 2.x (#4236, #4237)
- Sync required numpy and dask to runtime versions (#4233, #4239)
- Fix wrong argument parsing in structural_similarity (#4246, #4247)
- Fix active contour gallery example after change to rc coordinates (#4257, #4262)

4 authors added to this release [alphabetical by first name or login]
---------------------------------------------------------------------
- François Boulogne
- Jarrod Millman
- Mark Harfouche
- Ondrej Pesek

6 reviewers added to this release [alphabetical by first name or login]
-----------------------------------------------------------------------
- Alexandre de Siqueira
- Egor Panfilov
- François Boulogne
- Juan Nunez-Iglesias
- Mark Harfouche
- Nelle Varoquaux


Announcement: scikit-image 0.16.1
=================================

We're happy to announce the release of scikit-image v0.16.1!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org

Starting from this release, scikit-image will follow the recently
introduced NumPy deprecation policy, `NEP 29
<https://github.com/numpy/numpy/blob/master/doc/neps/nep-0029-deprecation_policy.rst>__`.
Accordingly, scikit-image 0.16 drops support for Python 3.5.
This release of scikit-image officially supports Python 3.6 and 3.7.

Special thanks to Matthias Bussonnier for `Frappuccino
<https://github.com/Carreau/frappuccino>`__, which helped us catch all API
changes and nail down the APIs for new features.

New Features
------------
- New `skimage.evaluate` module containing simple metrics (mse,
  nrme, psd) and segmentation metrics (adapted rand error, variation of
  information) (#4025)
- n-dimensional TV-L1 optical flow algorithm for registration --
  `skimage.registration.optical_flow_tvl1` (#3983)
- Draw a line in an n-dimensional array -- `skimage.draw.line_nd`
  (#2043)
- 2D Farid & Simoncelli edge filters - `skimage.filters.farid`,
  `skimage.filters.farid_h`, and `skimage.filters.farid_v` (#3775)
- 2D majority voting filter assigning to each pixel the most commonly
  occurring value within its neighborhood -- `skimage.filters.majority`
  (#3836, #3839)
- Multi-level threshold "multi-Otsu" method, a thresholding algorithm
  used to separate the pixels of an input image into several classes by
  maximizing the variances between classes --
  `skimage.filters.threshold_multiotsu` (#3872, #4174)
- New example data -- `skimage.data.shepp_logan_phantom`, `skimage.data.colorwheel`,
  `skimage.data.brick`, `skimage.data.grass`, `skimage.data.roughwall`, `skimage.data.cell`
  (#3958, #3966)
- Compute and format image region properties as a table --
  `skimage.measure.regionprops_table` (#3959)
- Convert a polygon into a mask -- `skimage.draw.poly2mask`  (#3971, #3977)
- Visual image comparison helper `skimage.util.compare_images`,
  that returns an image showing the difference between two input images (#4089)
- `skimage.transform.warp_polar` to remap image into
  polar or log-polar coordinates. (#4097)

Improvements
------------

- RANSAC: new option to set initial samples selected for initialization (#2992)
- Better repr and str for `skimage.transform.ProjectiveTransform` (#3525,
  #3967)
- Better error messages and data type stability to
  `skimage.segmentation.relabel_sequential` (#3740)
- Improved compatibility with dask arrays in some image thresholding methods (#3823)
- `skimage.io.ImageCollection` can now receive lists of patterns (#3928)
- Speed up `skimage.feature.peak_local_max` (#3984)
- Better error message when incorrect value for keyword argument `kind` in
  `skimage.color.label2rgb` (#4055)
- All functions from `skimage.drawing` now supports multi-channel 2D images (#4134)

API Changes
-----------
- Deprecated subpackage ``skimage.novice`` has been removed.
- Default value of ``multichannel`` parameters has been set to False in
  `skimage.transform.rescale`, `skimage.transform.pyramid_reduce`,
  `skimage.transform.pyramid_laplacian`,
  `skimage.transform.pyramid_gaussian`, and
  `skimage.transform.pyramid_expand`. Guessing is no longer performed for 3D
  arrays.
- Deprecated argument ``visualise`` has been removed from
  `skimage.feature.hog`. Use ``visualize`` instead.¨
- `skimage.transform.seam_carve` has been completely removed from the
  library due to licensing restrictions.
- Parameter ``as_grey`` has been removed from `skimage.data.load` and
  `skimage.io.imread`. Use ``as_gray`` instead.
- Parameter ``min_size`` has been removed from
  `skimage.morphology.remove_small_holes`. Use ``area_threshold`` instead.
- Deprecated ``correct_mesh_orientation`` in `skimage.measure` has been
  removed.
- `skimage.measure._regionprops` has been completely switched to using
  row-column coordinates. Old x-y interface is not longer available.
- Default value of ``behavior`` parameter has been set to ``ndimage`` in
  `skimage.filters.median`.
- Parameter ``flatten`` in `skimage.io.imread` has been removed in
  favor of ``as_gray``.
- Parameters ``Hxx, Hxy, Hyy`` have been removed from
  `skimage.feature.corner.hessian_matrix_eigvals` in favor of ``H_elems``.
- Default value of ``order`` parameter has been set to ``rc`` in
  `skimage.feature.hessian_matrix`.
- ``skimage.util.img_as_*`` functions no longer raise precision and/or loss warnings.

Bugfixes
--------

- Corrected error with scales attribute in ORB.detect_and_extract (#2835)
  The scales attribute wasn't taking into account the mask, and thus was using
  an incorrect array size.
- Correct for bias in Inverse Randon Transform (`skimage.transform.irandon`) (#3067)
  Fixed by using the Ramp filter equation in the spatial domain as described
  in the reference
- Fix a rounding issue that caused  a rotated image to have a
  different size than the input (`skimage.transform.rotate`)  (#3173)
- RANSAC uses random subsets of the original data and not bootstraps. (#3901,
  #3915)
- Canny now produces the same output regardless of dtype (#3919)
- Geometry Transforms: avoid division by zero & some degenerate cases (#3926)
- Fixed float32 support in denoise_bilateral and denoise_tv_bregman (#3936)
- Fixed computation of Meijering filter and avoid ZeroDivisionError (#3957)
- Fixed `skimage.filters.threshold_li` to prevent being stuck on stationnary
  points, and thus at local minima or maxima (#3966)
- Edited `skimage.exposure.rescale_intensity` to return input image instead of
  nans when all 0 (#4015)
- Fixed `skimage.morphology.medial_axis`. A wrong indentation in Cython
  caused the function to not behave as intended. (#4060)
- Fixed `skimage.restoration.denoise_bilateral` by correcting the padding in
  the gaussian filter(#4080)
- Fixed `skimage.measure.find_contours` when input image contains NaN.
  Contours interesting NaN will be left open (#4150)
- Fixed `skimage.feature.blob_log` and `skimage.feature.blob_dog` for 3D
  images and anisotropic data (#4162)
- Fixed `skimage.exposure.adjust_gamma`, `skimage.exposure.adjust_log`,
  and `skimage.exposure.adjust_sigmoid` such that when provided with a 1 by
  1 ndarray, it returns 1 by 1 ndarrays and not single number floats (#4169)

Deprecations
------------
- Parameter ``neighbors`` in `skimage.measure.convex_hull_object` has been
  deprecated in favor of ``connectivity`` and will be removed in version 0.18.0.
- The following functions are deprecated in favor of the `skimage.metrics`
  module (#4025):

    - `skimage.measure.compare_mse`
    - `skimage.measure.compare_nrmse`
    - `skimage.measure.compare_psnr`
    - `skimage.measure.compare_ssim`

- The function `skimage.color.guess_spatial_dimensions` is deprecated and
  will be removed in 0.18 (#4031)
- The argument ``bc`` in `skimage.segmentation.active_contour` is
  deprecated.
- The function `skimage.data.load` is deprecated and will be removed in 0.18
  (#4061)
- The function `skimage.transform.match_histogram` is deprecated in favor of
  `skimage.exposure.match_histogram` (#4107)
- The parameter ``neighbors`` of `skimage.morphology.convex_hull_object` is
  deprecated. 
- The `skimage.transform.randon_tranform` function will convert input image
  of integer type to float by default in 0.18. To preserve current behaviour,
  set the new argument ``preserve_range`` to True. (#4131)


Documentation improvements
--------------------------

- DOC: Improve the documentation of transform.resize with respect to the anti_aliasing_sigma parameter (#3911)
- Fix URL for stain deconvolution reference (#3862)
- Fix doc for denoise guassian (#3869)
- DOC: various enhancements (cross links, gallery, ref...), mainly for corner detection (#3996)
- [DOC] clarify that the inertia_tensor may be nD in documentation (#4013)
- [DOC] How to test and write benchmarks (#4016)
- Spellcheck @CONTRIBUTING.txt (#4008)
- Spellcheck @doc/examples/segmentation/plot_watershed.py (#4009)
- Spellcheck @doc/examples/segmentation/plot_thresholding.py (#4010)
- Spellcheck @skimage/morphology/binary.py (#4011)
- Spellcheck @skimage/morphology/extrema.py (#4012)
- docs update for downscale_local_mean and N-dimensional images (#4079)
- Remove fancy language from 0.15 release notes (#3827)
- Documentation formatting / compilation fixes (#3838)
- Remove duplicated section in INSTALL.txt. (#3876)
- ENH: doc of ridge functions (#3933)
- Fix docstring for Threshold Niblack (#3917)
- adding docs to circle_perimeter_aa (#4155)
- Update link to NumPy docstring standard in Contribution Guide (replaces #4191) (#4192)
- DOC: Improve downscale_local_mean() docstring (#4180)
- DOC: enhance the result display in ransac gallery example (#4109)
- Gallery: use fstrings for better readability (#4110)
- MNT: Document stacklevel parameter in contribution guide (#4066)
- Fix minor typo (#3988)
- MIN: docstring improvements in canny functions (#3920)
- Minor docstring fixes for #4150 (#4184)
- Fix `full` parameter description in compare_ssim (#3860)
- State Bradley threshold equivalence in Niblack docstring (#3891)
- Add plt.show() to example-code for consistency. (#3908)
- CC0 is not equivalent to public domain. Fix the note of the horse image (#3931)
- Update the joblib link in tutorial_parallelization.rst (#3943)
- Fix plot_edge_filter.py references (#3946)
- Add missing argument to docstring of PaintTool (#3970)
- Improving documentation and tests for directional filters (#3956)
- Added new thorough examples on the inner working of
  ``skimage.filters.threshold_li`` (#3966)
- matplotlib: remove interpolation=nearest, none in our examples (#4002)
- fix URL encoding for wikipedia references in filters.rank.entropy and filters.rank.shannon_entropy docstring (#4007)
- Fixup integer division in examples (#4032)
- Update the links the installation guide (#4118)
- Gallery hough line transform (#4124)
- Cross-linking between function documentation should now be much improved! (#4188)
- Better documentation of the ``num_peaks`` of `skimage.feature.corner_peaks` (#4195)


Other Pull Requests
-------------------
- Add benchmark suite for exposure module (#3312)
- Remove precision and sign loss warnings from ``skimage.util.img_as_`` (#3575)
- Propose SKIPs and add mission/vision/values, governance (#3585)
- Use user-installed tifffile if available (#3650)
- Simplify benchmarks pinnings (#3711)
- Add project_urls to setup for PyPI and other services (#3834)
- Address deprecations for 0.16 release (#3841)
- Followup deprecations for 0.16 (#3851)
- Build and test the docs in Azure (#3873)
- Pin numpydoc to pre-0.8 to fix dev docs formatting (#3893)
- Change all HTTP links to HTTPS (#3896)
- Skip extra deps on OSX (#3898)
- Add location for Sphinx 2.0.1 search results; clean up templates (#3899)
- Fix CSS styling of Sphinx 2.0.1 + numpydoc 0.9 rendered docs (#3900)
- Travis CI: The sudo: tag is deprcated in Travis (#4164)
- MNT Preparing the 0.16 release (#4204)
- FIX generate_release_note when contributor_set contains None (#4205)
- Specify that travis should use Ubuntu xenial (14.04) not trusty (16.04) (#4082)
- MNT: set stack level accordingly in lab2xyz (#4067)
- MNT: fixup stack level for filters ridges (#4068)
- MNT: remove unused import `deprecated` from filters.thresholding (#4069)
- MNT: Set stacklevel correctly in io matplotlib plugin (#4070)
- MNT: set stacklevel accordingly in felzenszwalb_cython (#4071)
- MNT: Set stacklevel accordingly in img_as_* (convert) (#4072)
- MNT: set stacklevel accordingly in util.shape (#4073)
- MNT: remove extreneous matplotlib warning (#4074)
- Suppress warnings in tests for viewer (#4017)
- Suppress warnings in test suite regarding measure.label (#4018)
- Suppress warnings in test_rank due to type conversion (#4019)
- Add todo item for imread plugin testing (#3907)
- Remove matplotlib agg warning when using the sphinx gallery. (#3897)
- Forward-port release notes for 0.14.4 (#4137)
- Add tests for pathological arrays in threshold_li (#4143)
- setup.py: Fail gracefully when NumPy is not installed (#4181)
- Drop Python 3.5 support (#4102)
- Force imageio reader to return NumPy arrays (#3837)
- Fixing connecting to GitHub with SSH info. (#3875)
- Small fix to an error message of `skimage.measure.regionprops` (#3884)
- Unify skeletonize and skeletonize 3D APIs (#3904)
- Add location for Sphinx 2.0.1 search results; clean up templates (#3910)
- Pin numpy version forward (#3925)
- Replacing pyfits with Astropy to read FITS (#3930)
- Add warning for future dtype kwarg removal (#3932)
- MAINT: cleanup regionprop add PYTHONOPTIMIZE=2 to travis array (#3934)
- Adding complexity and new tests for filters.threshold_multiotsu (#3935)
- Fixup dtype kwarg warning in certain image plugins (#3948)
- don't cast integer to float before using it as integer in numpy logspace (#3949)
- avoid low contrast image save in a doctest. (#3953)
- MAINT: Remove unused _convert_input from filters._gaussian (#4001)
- Set minimum version for imread so that it compiles from source on linux in test builds (#3960)
- Cleanup plugin utilization in data.load and testsuite (#3961)
- Select minimum imageio such that it is compatible with pathlib (#3969)
- Remove pytest-faulthandler from test dependencies (#3987)
- Fix tifffile and __array_function__ failures in our CI (#3992)
- MAINT: Do not use assert in code, raise an exception instead. (#4006)
- Enable packagers to disable failures on warnings. (#4021)
- Fix numpy 117 rc and dask in thresholding filters (#4022)
- silence r,c  warnings when property does not depend on r,c (#4027)
- remove warning filter, fix doc wrt r,c (#4028)
- Import Iterable from collections.abc (#4033)
- Import Iterable from collections.abc in vendored tifffile code (#4034)
- Correction of typos after #4025 (#4036)
- Rename internal function called assert_* -> check_* (#4037)
- Improve import time (#4039)
- Remove .meeseeksdev.yml (#4045)
- Fix mpl deprecation on grid() (#4049)
- Fix gallery after deprecation from #4025 (#4050)
- fix mpl future deprecation normed -> density (#4053)
- Add shape= to circle perimeter in hough_circle example (#4047)
- Critical: address internal warnings in test suite related to metrics 4025 (#4063)
- Use functools instead of a real function for the internal warn function (#4062)
- Test rank capture warnings in threadsafe manner (#4064)
- Make use of FFTs more consistent across the library (#4084)
- Fixup region props test (#4099)
- Turn single backquotes to double backquotes in filters (#4127)
- Refactor radon transform module (#4136)
- Fix broken import of rgb2gray in benchmark suite (#4176)
- Fix doc building issues with SKIPs (#4182)
- Remove several __future__ imports (#4198)
- Restore deprecated coordinates arg to regionprops (#4144)
- Refactor/optimize threshold_multiotsu (#4167)
- Remove Python2-specific code (#4170)
- `view_as_windows` incorrectly assumes that a contiguous array is needed  (#4171)
- Handle case in which NamedTemporaryFile fails (#4172)
- Fix incorrect resolution date on SKIP1 (#4183)
- API updates before 0.16 (#4187)
- Fix conversion to float32 dtype (#4193)


Contributors to this release
----------------------------

- Abhishek Arya
- Alexandre de Siqueira
- Alexis Mignon
- Anthony Carapetis
- Bastian Eichenberger
- Bharat Raghunathan
- Christian Clauss
- Clement Ng
- David Breuer
- David Haberthür
- Dominik Kutra
- Dominik Straub
- Egor Panfilov
- Emmanuelle Gouillart
- Etienne Landuré
- François Boulogne
- Genevieve Buckley
- Gregory R. Lee
- Hadrien Mary
- Hamdi Sahloul
- Holly Gibbs
- Huang-Wei Chang
- i3v (i3v)
- Jarrod Millman
- Jirka Borovec
- Johan Jeppsson
- Johannes Schönberger
- Jon Crall
- Josh Warner
- Juan Nunez-Iglesias
- Kaligule (Kaligule)
- kczimm (kczimm)
- Lars Grueter
- Shachar Ben Harim
- Luis F. de Figueiredo
- Mark Harfouche
- Mars Huang
- Dave Mellert
- Nelle Varoquaux
- Ollin Boer Bohan
- Patrick J Zager
- Riadh Fezzani
- Ryan Avery
- Srinath Kailasa
- Stefan van der Walt
- Stuart Berg
- Uwe Schmidt


Reviewers for this release
--------------------------

- Alexandre de Siqueira
- Anthony Carapetis
- Bastian Eichenberger
- Clement Ng
- David Breuer
- Egor Panfilov
- Emmanuelle Gouillart
- Etienne Landuré
- François Boulogne
- Genevieve Buckley
- Gregory R. Lee
- Hadrien Mary
- Hamdi Sahloul
- Holly Gibbs
- Jarrod Millman
- Jirka Borovec
- Johan Jeppsson
- Johannes Schönberger
- Jon Crall
- Josh Warner
- jrmarsha
- Juan Nunez-Iglesias
- kczimm
- Lars Grueter
- leGIT-bot
- Mark Harfouche
- Mars Huang
- Dave Mellert
- Paul Müller
- Phil Starkey
- Ralf Gommers
- Riadh Fezzani
- Ryan Avery
- Sebastian Berg
- Stefan van der Walt
- Uwe Schmidt

scikit-image 0.13.1 is a bug-fix and compatibility update. See below for
the many new features in 0.13.0.

The main contribution in 0.13.1 is Jarrod Millman's valiant work to ensure
scikit-image works with both NetworkX 1.11 and 2.0 (#2766). Additional updates
include:

- Bug fix in similarity transform estimation, by GitHub user @zhongzyd (#2690)
- Bug fixes in ``skimage.util.plot_matches`` and ``denoise_wavelet``,
  by Gregory Lee (#2650, #2640)
- Documentation updates by Egor Panfilov (#2716) and Jirka Borovec (#2524)
- Documentation build fixes by Gregory Lee (#2666, #2731), Nelle
  Varoquaux (#2722), and Stéfan van der Walt (#2723, #2810)


Announcement: scikit-image 0.13.0
=================================

We're happy to (finally) announce the release of scikit-image v0.13.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

http://scikit-image.org

and our gallery of examples

http://scikit-image.org/docs/dev/auto_examples/

Highlights
----------

This release is the result of a year of work, with over 200 pull requests by
82 contributors. Highlights include:

- Improved n-dimensional image support. This release adds nD support to:

  * ``regionprops`` computation for centroids (#2083)
  * ``segmentation.clear_border`` (#2087)
  * Hessian matrix (#2194)

- In addition, the following new functions support nD images:

  * new wavelet denoising function, ``restoration.denoise_wavelet``
    (#1833, #2190, #2238, #2240, #2241, #2242, #2462)
  * new thresholding functions, ``filters.threshold_sauvola`` and
    ``filters.threshold_niblack`` (#2266, #2441)
  * new local maximum, local minimum, hmaxima, hminima functions (#2449)

- Grey level co-occurrence matrix (GLCM) now works with uint16 images
- ``filters.try_all_threshold`` to rapidly see output of various thresholding
  methods
- Frangi and Hessian filters (2D only) (#2153)
- New *compact watershed* algorithm in ``segmentation.watershed`` (#2211)
- New *shape index* algorithm in ``feature.shape_index`` (#2312)

New functions and features
--------------------------

- Add threshold minimum algorithm (#2104)
- Implement mean and triangle thresholding (#2126)
- Add Frangi and Hessian filters (#2153)
- add bbox_area to region properties (#2187)
- colorconv: Add rgba2rgb() (#2181)
- Lewiner marching cubes algorithm (#2052)
- image inversion (#2199)
- wavelet denoising (from #1833) (#2190)
- routine to estimate the noise standard deviation from an image (#1837)
- Add compact watershed and clean up existing watershed (#2211)
- Added the missing 'grey2rgb' function. (#2316)
- Shape index (#2312)
- Fundamental and essential matrix 8-point algorithm (#1357)
- Add YUV, YIQ, YPbPr, YCbCr colorspaces
- Detection of local extrema from morphology (#2449)
- shannon entropy (#2416)

Documentation improvements
--------------------------

- add details about github SSH keys in contributing page (#2073)
- Add example for felzenszwalb image segmentation (#2096)
- Sphinx gallery for example gallery (#2078)
- Improved region boundary RAG docs (#2106)
- Add gallery Lucy-Richardson deconvolution algorithm (#2376)
- Gallery: Use Horse to illustrate Convex Hull (#2431)
- Add working with OpenCV in user guide (#2519)

Code improvements
-----------------

- Remove lena image from test suite (#1985)
- Remove duplicate mean calculation in skimage.feature.match_template (#1980)
- Add nD support to clear_border (#2087)
- Add uint16 images support for co-occurrence matrix (#2095)
- Add default parameters for Gaussian and median filters (#2151)
- try_all to choose the best threshold algorithm (#2110)
- Add support for multichannel in Felzenszwalb segmentation (#2134)
- Improved SimilarityTransform, new EuclideanTransform class (#2044)
- ENH: Speed up Hessian matrix computation (#2194)
- add n-dimensional support to denoise_wavelet (#2242)
- Speedup ``inpaint_biharmonic`` (#2234)
- Update hessian matrix code to include order kwarg (#2327)
- Handle cases for label2rgb where input labels are negative and/or
  nonconsecutive (#2370)
- Added watershed_line parameter (#2393)

API Changes
-----------

- Remove deprecated ``filter`` module. Use ``filters`` instead. (#2023)
- Remove ``skimage.filters.canny`` links. Use ``feature.canny`` instead. (#2024)
- Removed Python 2.6 support and related checks (#2033)
- Remove deprecated {h/v}sobel, {h/v}prewitt, {h/v}scharr,
  roberts_{positive/negative} filters (#2159)
- Remove deprecated ``_mode_deprecations`` (#2156)
- Remove deprecated None defaults in ``rescale_intensity`` (#2161)
- Parameters ``ntiles_x`` and ``ntiles_y`` have been removed from
  ``exposure.equalize_adapthist``
- The minimum NumPy version is now 1.11, and the minimum SciPy version is now
  0.17

Deprecations
------------

- clip_negative will be set to false by default in version 0.15
  (func: dtype_limits) (#2228)
- Deprecate "dynamic_range" in favor of "data_range" (#2384)
- The default value of the ``circle`` argument to ``radon`` and ``iradon``
  transforms will be ``True`` in 0.15 (#2235)
- The default value of ``multichannel`` for ``denoise_bilateral`` and
  ``denoise_nl_means`` will be ``False`` in 0.15
- The default value of ``block_norm`` in ``feature.hog`` will be L2-Hysteresis in
  0.15.
- The ``threshold_adaptive`` function is deprecated. Use ``threshold_local``
  instead.
- The default value of ``mode`` in ``transform.swirl``, ``resize``, and ``rescale``
  will be "reflect" in 0.15.

Contributors to this release
----------------------------

- AbdealiJK
- Rodrigo Benenson
- Vighnesh Birodkar
- Jirka Borovec
- François Boulogne
- Matthew Brett
- Sarwat Fatima
- Rachel Finck
- Joe Futrelle
- Jeroen Van Goey
- Christoph Gohlke
- Roman Golovanov
- Emmanuelle Gouillart
- Anshita Gupta
- David Haberthür
- Jeff Hemmelgarn
- Hiyorimi
- Daniel Hyams
- Alex Izvorski
- Kyle Jackson
- Jirka
- JohnnyTeutonic
- Kevin Keraudren
- Almar Klein
- Yu Kobayashi
- Moriyoshi Koizumi
- Lachlan
- LachlanD
- George Laurent
- Gregory R. Lee
- Evan Limanto
- Ben Longo
- Victor MARTIN
- Oliver Mader
- Ken'ichi Matsui
- Jeremy Metz
- Jeyson Molina
- Michael Mueller
- Juan Nunez-Iglesias
- Egor Panfilov
- Paul
- PengchengAi
- Francisco de la Peña
- Pavlin Poličar
- Orion Poplawski
- Zoe Richards
- Todd V. Rovito
- Christian Sachs
- Sanya
- Johannes Schönberger
- Pavel Shevchuk
- Scott Sievert
- Steven Silvester
- Shaun Singh
- Sourav Singh
- Alexandre Fioravante de Siqueira
- Samuel St-Jean
- Noah Stier
- Ole Streicher
- Martin Thoma
- Matěj Týč
- Viraj
- Stefan van der Walt
- Josh Warner
- Olivia Wilson
- Robin Wilson
- Martin Zackrisson
- Yue Zheng
- Nick Zoghb
- alexandrejaguar
- almar
- cespenel
- danielballan
- dmesejo
- eli
- jwittenbach
- lgeorge
- mljli
- rjeli
- skrish13
- tseclaudia
- walter

Pull requests merged in this release
------------------------------------

- Warn if user tries to build with older Cython version (#1986)
- Remove lena image from test suite (#1985)
- Add inpaint to module init (#1987)
- Pre-calculate tempate mean (#1980)
- rgb2grey -> grey2rgb (#1989)
- Also expose rgb2gray as rgb2grey (#1990)
- Remove all .md5 files on clean (#1992)
- avoid deprecation warnings when calling compute_ssim with multichannel=True (#1994)
- DOC: Suggest multichannel=True in compute_ssim error (#1999)
- [DOC] add link to guide (#2001)
- Fix docs-->doc in CONTRIBUTING (#2009)
- Turn ``dask`` into an optional dependency (#2013)
- Correct regexp for catching mpl warnings (#2014)
- BUILD: Use --pre flag for Travis pip installs. (#1938)
- Github templates (#1954)
- added doc to PaintTool (#1934)
- skimage.segmentation.quickshift signature is missing from API docs (#2017)
- MAINT: Upgrade tifffile (#2016)
- Modified .gitignore to properly ignore auto_example files (#1966)
- MAINT: Switch from coveralls -> codecov in CI build (#2015)
- skimage.segmentation.quickshift signature is missing from API docs, third attempt (#2021)
- MAINT: Remove deprecated ``filter`` module (#2023)
- Remove ``skimage.filters.canny`` links (#2024)
- Document regionprops bbox property. (#2030)
- Fix URL to texturematch paper (#2031)
- Improved skimage.segmentation.active_contour input arguments' dtype support (#2032)
- Fix local test function (#2034)
- Removed Python 2.6 support and related checks (#2033)
- Test on OSX (#2038)
- Change coverage badge to codecov (#2055)
- TST: Speed up bilateral filter tests (#2061)
- Speed up colorconv._convert (#2064)
- FIX: Fix import of 'warn' in qt_plugin (#2070)
- Add YUV, YIQ, YPbPr, YCbCr colorspaces
- adding details about github SSH keys in contributing page (#2073)
- ENH: Pass np.random.RandomState to RANSAC (#2072)
- Handle IO objects with tifffile (#2046)
- Updated centroid to use coords - works in 3d (#2083)
- [WIP] Hierarchical Merging of Region Boundary RAGs (#2058)
- Add nD support to clear_border (#2087)
- DOC: update for new API (minor) (#2090)
- Add example for felzenszwalb image segmentation (#2096)
- DOC: add space before column on variable def (minor...) (#2102)
- DOC: Guide new contributors to HTTPS, not SSH (#2082)
- Add François Boulogne to the mailmap (#2117)
- Move skimage.filters.rank description and todos from README into docstring. (#2115)
- Fixing Error and documentation on Otsu Threshold (#2118)
- Add scuinto's second email address to mailmap (#2122)
- MAINT: around label and regionprops functions. (#2100)
- Add threshold minimum algorithm (#2104)
- Sphinx gallery for example gallery (#2078)
- DOC: make a title shorter in gallery (#2128)
- DOC: refactor axes with lists (#2129)
- DOC ENH + API fix on houghline transform (#2089)
- Fix indentation for example script (#2136)
- Implement mean and triangle thresholding (#2126)
- Move ``skimage.measure.label`` references to the docstring (#2143)
- Fix outdated GraphicsGems link (#2149)
- Docstring (#2145)
- Add uint16 images support for co-occurrence matrix (#2095)
- Remove deprecared {h/v}sobel, {h/v}prewitt, {h/v}scharr, roberts_{positive/negative} filters (#2159)
- Remove deprecated ``_mode_deprecations`` (#2156)
- Default parameters (#2151)
- ENH: try_all to choose the best threshold algorithm and DOC refactoring (#2110)
- BUGFIX: inverse_map should not be None (#2160)
- Switched felzenszwalb gray to multichannel version (#2134)
- Writing, style, and PEP8 fixes for greycomatrix (#2157)
- Add Frangi and Hessian filters (#2153)
- Improved SimilarityTransform, new EuclideanTransform class (#2044)
- color.colorconv: Fix documentation of rgb2gray() (#2169)
- fix region merging in ``segmentation.felzenszwalb`` (#2164)
- Remove deprecated None defaults in ``rescale_intensity`` (#2161)
- DOC: add a note to template_match (#2176)
- Added chapter title formatting for numpy_images.rst (#2177)
- Fix threshold_triangle to work with non-integer images. (#2171)
- Improved region boundary RAG docs (#2106)
- ENH add bbox_area to region properties (#2187)
- colorconv: Add rgba2rgb() (#2181)
- DOC: add DOI to references (#2188)
- remove local threshold in try_all_threshold (#2180)
- DOC: add a note on warning treatment (#2198)
- ENH: Speed up Hessian matrix computation (#2194)
- Add missing unittests for data and convert horse to binary (#2196)
- Fix ssim example (#2208)
- [MRG] MAINT: Replaced gaussian_filter with filters.gaussian (#2210)
- [MRG] DOC: corrected mssim docstring to return float (#2218)
- FEAT: Lewiner marching cubes algorithm (#2052)
- Fix bug in salt and pepper noise (#2223)
- TST: Updated AppVeyor to use Conda, added msvc_runtime (#2217)
- Improve docstrings for captions (#2185)
- Add task update version on wikipedia (#2230)
- NEW + DOC: image inversion (#2199)
- ENH: Implements wavelet denoising (from #1833) (#2190)
- TEST: define seed in setup() / Fix random test failure (#2227)
- add n-dimensional support to denoise_wavelet (#2242)
- API: clip_negative will be set to false by default in version 0.15 (func: dtype_limits) (#2228)
- Speedup ``inpaint_biharmonic`` (#2234)
- MAINT dtype.py (PEP8) (#2231)
- Removed unused extend_image (#2251)
- ENH:  routine to estimate the noise standard deviation from an image (#1837)
- Restrict sphinx builds to a single process.  Remove vendored numpydoc. (#2257)
- Added more specific check for image shape in threshold_otsu warning (#2259)
- Allow running ``setup.py egg_info`` without numpy installed. (#2260)
- Add compact watershed and clean up existing watershed (#2211)
- Use numpy.pad directly, removing most shipped code in util.pad (#2265)
- DOC: fix references (#2262)
- DOC: tiny fixes in gallery (#2226)
- DOC: fix typo (#2274)
- Update Manifest.in (#2255)
- Bugfix unbounded correlation -- Dhyams fix for match template (#2263)
- DOC: Refactor example skeletonize in the gallery (#2141)
- [MRG+1] Insert metadata in docstrings of images in skimage.data.* (#2236)
- MAINT: Radon (docstring, API, PEP8) (#2235)
- [MRG+2] MAINT: Fix numpy deprecation (#2283)
- Reduce whitespace around plots (#2144)
- [MRG+1] By default, clear_border is not inplace (#2285)
- Remove unused imports in ``transform.{pyx/pxd}`` (#2288)
- [MRG+1] Add community guidelines to doc navigation (#2287)
- Adding colors to the IHC (#2279)
- FIX: select num_peaks if labels is specified  (#2098)
- [MRG+1] Add felzenszwalb shape validation (#2286)
- [MRG+1] more closesly match the BayesShrink paper in _wavelet_threshold (#2241)
- Remove usages of ``subplots_adjust`` (#2289)
- [MRG+1] Change documentation page favicon (#2291)
- [MRG+1] TST: prefer ``assert_`` from numpy.testing over assert (#2298)
- TSTFIX: Bug fix for development version of scipy (#2302)
- Enhance ``compare_ssim`` docstring (#2314)
- Added the missing 'grey2rgb' function. (#2316)
- PEP8 (#2304)
- Made Python wrappers for public Cython functions (#2303)
- Update mailing list location (#2328)
- Shape Index (#2312)
- Add pywavelets to runtime requirements in DEPENDS.txt (#2238)
- Refactor variable names in ``skimage.draw`` (#2321)
- Fix display problem when printing error messages (#2326)
- Added catch for zero image in threshold_li (#2338)
- FIX: Modified peak_local_max to use relabel_sequential (#2341)
- Update favicon in _static (#2355)
- Remove incorrect input type assumption in doctrings for rgb2hsv and h… (#2354)
- Update the default boundary mode in transform.swirl (#2331)
- Update imread() document (#2358)
- Check for valid mode in random_walker(). (#2362)
- Fix 1 broken test in _shared not executed by nose/travis (#2229)
- Update hessian matrix code to include order kwarg (#2327)
- Clarify purpose of beta1 and beta2 parameters in documentations of sk… (#2382)
- Handle cases for label2rgb where input labels are negative and/or nonconsecutive (#2370)
- Update ``exposure.equalize_adapthist`` args and docstring (#2220)
- Fix (x, y) origin description in user guide (#2385)
- Update docstring for show_rag method (#2375)
- Fix display problem when printing error messages (#2372)
- Added a check for empty array in _shared.utils.py (#2364)
- Fix no peaks blob log (#2349)
- ENH: Extend draw.ellipse with orientation kwarg (#2366)
- Fundamental and essential matrix 8-point algorithm (#1357)
- Fix reference to travis notes (#2403)
- Fix deprecated option in sphinx that causes warning treated as error in travis (#2395)
- Update Travis Script (#2374)
- Remove the freeimage plugin (#1933)
- Fix shape type for histogram (#2417)
- Add illuminant and observer parameters to the rgb2lab and lab2rgb functions. (#2306)
- PEP8 (#2413)
- MAINT: merge lists of dtypes (#2420)
- Made (partially) ``pep8``-compliant (#2392)
- Added titles and text to make plot_brief.py example more clear (#2193)
- DOC: Add reference to standard illuminant (#2418)
- Added titles and text to the subplots to make it easier to new comers for plot_censure.py example (#2191)
- Deprecate "dynamic_range" in favor of "data_range" (#2384)
- Make PR 2266 n-D compatible (#4)
- Add new "thin" method based on Guo and Hall 1989 (#2294)
- local threshold niblack sauvola (from Jeysonmc PR) (#2266)
- stable ellipse fitting (#2394)
- Add gallery Lucy-Richardson deconvolution algorithm (#2376)
- Improve SIFT loader docstring according to comments and StackOverflow (#2404)
- Change to Javascript loading of search index (patch by Julian Taylor) (#2438)
- Fix segfault in connected components (patch by Yaroslav Halchenko) (#2437)
- Refactor ``util/dtype.py`` (#2425)
- ENH: Gallery, various little stylish corrections (DFT example). (#2430)
- Make peak_local_max return indices sorted, always (#2435)
- Correct comment of probabilistic_hough_line(). (#2448)
- Added watershed_line parameter (#2393)
- Solved Gaussian value range #2383 (#2388)
- Gallery: Use Horse to illustrate Convex Hull (#2431)
- MRG: update build matrix for Python 3.6 (#2451)
- Wavelet denoising in YCbCr color space (#2240)
- Gallery: Use gray cmap for coins (#2459)
- Bug fix for Sauvola and Niblack thresholding (#2441)
- MAINT: removes _wavelet_threshold docstring (#2460)
- BUG: fix denoise_wavelet for odd-length input (#2462)
- MAINT: warns for new multichannel default in denoise_{bilateral, nl_means} (#2467)
- Various enhancements in gallery for denoising (#2461)
- Tool for checking completeness of sdist (#2085)
- Add different ``skimage.hog`` blocks normalization methods (#2040)
- DOC: fix typos and add references (#2478)
- update sphinx gallery to 0.1.8 (#2474)
- DOC: Fix typo in gaussian filter docstring (#2487)
- Add threshold_local, deprecate old threshold_adaptive API (#2490)
- Default edge mode change for resize and rescale (#2484)
- Add ``dask[array]`` to optional requirements (#2494)
- DOC:  Adds an instruction to CONTRIBUTING.txt & Updates the git install link for Windows (#2495)
- ENH: generalize hough_peak functions (#2109)
- Fix gallery examples (#2504)
- Bump min scipy version (#2254)
- DOC: img_as_float add note about range if input dtype is float (#2499)
- Update tifffile for 2017.01.12 changes (#2497)
- Replace local_sum by block_reduce in docstrings. (#2498)
- MAINT: pass scipys truncate parameter to gaussian filter API (#2508)
- DOC: gallery: join segmentation: enhancement (#2507)
- Tidy up the deployment of dev docs (#2516)
-  Do not require cython for normal builds (#2509)
- Fix broken ``test_ncut_stable_subgraph`` for Python 3.6, enable Python 3.6 in Travis (#2511)
- Improved background labeling (#2381)
- For imread's load_func, make the img_num argument optional (#2054)
- Make compatible with current networkx master (#2455)
- Miscellaneous tidying in HOG code (#2526)
- BUG: Fix NumPy error when no descriptors are returned by ORB (#2537)
- BUG: ValueError in restoration.denoise_bilateral for zeros image (#2533)
- Fix link to Python XY (#2542)
- TST: fix ValueError with scipy-0.19.0rc2 (#2544)
- DOC: Update URL for data.coins() (#2548)
- Replace GRIN URL with Flickr URL (#2547)
- Have ``threshold_minimum`` return identical results on i686 and x86_64 (#2549)
- Minor Fix (Issue #2554) (#2556)
- Remove ``offset`` parameter from ``filters.threshold_sauvola`` docstring (#2566)
- Practical guide to reading video files (#1012)
- Remove dask from ``requirements.txt`` (#2572)
- Fix ``morphology.watershed`` error message (#2570)
- DOC: Added working with OpenCV in user guide (#2519)
- NEW: add shannon entropy (#2416)
- Fix typo in ylabel of GLCM demo (#2576)
- Detection of local extrema from morphology (#2449)
- Add extrema functions to ``__init__`` (#2588)

Announcement: scikit-image 0.9.0
================================

We're happy to announce the release of scikit-image v0.9.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

    http://scikit-image.org


New Features
------------

`scikit-image` now runs without translation under both Python 2 and 3.

In addition to several bug fixes, speed improvements and examples, the 204 pull
requests merged for this release include the following new features (PR number
in brackets):

Segmentation:

- 3D support in SLIC segmentation (#546)
- SLIC voxel spacing (#719)
- Generalized anisotropic spacing support for random_walker (#775)
- Yen threshold method (#686)

Transforms and filters:

- SART algorithm for tomography reconstruction (#584)
- Gabor filters (#371)
- Hough transform for ellipses (#597)
- Fast resampling of nD arrays (#511)
- Rotation axis center for Radon transforms with inverses. (#654)
- Reconstruction circle in inverse Radon transform (#567)
- Pixelwise image adjustment curves and methods (#505)

Feature detection:

- [experimental API] BRIEF feature descriptor (#591)
- [experimental API] Censure (STAR) Feature Detector (#668)
- Octagon structural element (#669)
- Add non rotation invariant uniform LBPs (#704)

Color and noise:

- Add deltaE color comparison and lab2lch conversion (#665)
- Isotropic denoising (#653)
- Generator to add various types of random noise to images (#625)
- Color deconvolution for immunohistochemical images (#441)
- Color label visualization (#485)

Drawing and visualization:

- Wu's anti-aliased circle, line, bezier curve (#709)
- Linked image viewers and docked plugins (#575)
- Rotated ellipse + bezier curve drawing (#510)
- PySide & PyQt4 compatibility in skimage-viewer (#551)

Other:

- Python 3 support without 2to3. (#620)
- 3D Marching Cubes (#469)
- Line, Circle, Ellipse total least squares fitting and RANSAC algorithm (#440)
- N-dimensional array padding (#577)
- Add a wrapper around `scipy.ndimage.gaussian_filter` with useful default behaviors. (#712)
- Predefined structuring elements for 3D morphology (#484)


API changes
-----------

The following backward-incompatible API changes were made between 0.8 and 0.9:

- No longer wrap ``imread`` output in an ``Image`` class
- Change default value of `sigma` parameter in ``skimage.segmentation.slic``
  to 0
- ``hough_circle`` now returns a stack of arrays that are the same size as the
  input image. Set the ``full_output`` flag to True for the old behavior.
- The following functions were deprecated over two releases:
  `skimage.filter.denoise_tv_chambolle`,
  `skimage.morphology.is_local_maximum`, `skimage.transform.hough`,
  `skimage.transform.probabilistic_hough`,`skimage.transform.hough_peaks`.
  Their functionality still exists, but under different names.


Contributors to this release
----------------------------

This release was made possible by the collaborative efforts of many
contributors, both new and old.  They are listed in alphabetical order by
surname:

- Ankit Agrawal
- K.-Michael Aye
- Chris Beaumont
- François Boulogne
- Luis Pedro Coelho
- Marianne Corvellec
- Olivier Debeir
- Ferdinand Deger
- Kemal Eren
- Jostein Bø Fløystad
- Christoph Gohlke
- Emmanuelle Gouillart
- Christian Horea
- Thouis (Ray) Jones
- Almar Klein
- Xavier Moles Lopez
- Alexis Mignon
- Juan Nunez-Iglesias
- Zachary Pincus
- Nicolas Pinto
- Davin Potts
- Malcolm Reynolds
- Umesh Sharma
- Johannes Schönberger
- Chintak Sheth
- Kirill Shklovsky
- Steven Silvester
- Matt Terry
- Riaan van den Dool
- Stéfan van der Walt
- Josh Warner
- Adam Wisniewski
- Yang Zetian
- Tony S Yu
Announcement: scikit-image 0.X.0
================================

We're happy to announce the release of scikit-image v0.X.0!

scikit-image is an image processing library for the scientific Python
ecosystem that includes algorithms for segmentation, geometric
transformations, feature detection, registration, color space
manipulation, analysis, filtering, morphology, and more.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org


New Features
------------



Improvements
------------



API Changes
-----------

- All references to EN-GB spelling for the word ``neighbour`` and others—e.g.,
  ``neigbourhood``, ``neighboring``, were changed to their EN-US spelling,
  ``neighbor``. With that, ``skimage.measure.perimeter` parameter ``neighbourhood``
  was deprecated in favor of ``neighborhood`` in 0.19.2.

Bugfixes
--------



Deprecations
------------



Contributors to this release
----------------------------
Announcement: scikits-image 0.7.0
=================================

We're happy to announce the 7th version of scikits-image!

Scikits-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website

  http://skimage.org


New Features
------------

It's been only 3 months since scikits-image 0.6 was released, but in that short
time, we've managed to add plenty of new features and enhancements, including

- Geometric image transforms
- 3 new image segmentation routines (Felsenzwalb, Quickshift, SLIC)
- Local binary patterns for texture characterization
- Morphological reconstruction
- Polygon approximation
- CIE Lab color space conversion
- Image pyramids
- Multispectral support in random walker segmentation
- Slicing, concatenation, and natural sorting of image collections
- Perimeter and coordinates measurements in regionprops
- An extensible image viewer based on Qt and Matplotlib, with plugins for edge
  detection, line-profiling, and viewing image collections

Plus, this release adds a number of bug fixes, new examples, and performance
enhancements.


Contributors to this release
----------------------------

This release was only possible due to the efforts of many contributors, both
new and old.

- Andreas Mueller
- Andreas Wuerl
- Andy Wilson
- Brian Holt
- Christoph Gohlke
- Dharhas Pothina
- Emmanuelle Gouillart
- Guillaume Gay
- Josh Warner
- James Bergstra
- Johannes Schonberger
- Jonathan J. Helmus
- Juan Nunez-Iglesias
- Leon Tietz
- Marianne Corvellec
- Matt McCormick
- Neil Yager
- Nicolas Pinto
- Nicolas Poilvert
- Pavel Campr
- Petter Strandmark
- Stefan van der Walt
- Tim Sheerman-Chase
- Tomas Kazmar
- Tony S Yu
- Wei Li
Announcement: scikit-image 0.12
===============================

We're happy to announce the release of scikit-image v0.12!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

http://scikit-image.org

and our gallery of examples

http://scikit-image.org/docs/dev/auto_examples/

Highlights and new features
---------------------------

For this release, we merged over 200 pull requests with bug fixes,
cleanups, improved documentation and new features.  Highlights
include:

- 3D skeletonization (#1923)
- Parallelization framework using ``dask``:``skimage.util.apply_parallel``
  (#1493)
- Laplacian operator ``filters.laplace`` (#1763)
- Synthetic 2-D and 3-D binary data with rounded blobs (#1485)
- Plugin for ``imageio`` library (#1575)
- Inpainting algorithm (#1804)
- New handling of background pixels for ``measure.label``: 0-valued
  pixels are considered as background by default, and the label of
  background pixels is 0.
- Partial support of 3-D images for ``skimage.measure.regionprops``
  (#1505)
- Multi-block local binary patterns (MB-LBP) for texture classification (#1536)
- Seam Carving (resizing without distortion) (#1459)
- Simple image comparison metrics (PSNR, NRMSE) (#1897)
- Region boundary based region adjacency graphs (RAG) (#1495)
- Construction of RAG from label images (#1826)
- ``morphology.remove_small_holes`` now complements
  ``morphology.remove_small_objects`` (#1689)
- Faster cython implementation of ``morphology.skeletonize``
- More appropriate default weights in
  ``restoration.denoise_tv_chambolle`` and ``feature.peak_local_max``
- Correction of bug in the computation of the Euler characteristic in
  ``measure.regionprops``.
- Replace jet by viridis as default colormap
- Alpha layer support for ``color.gray2rgb``
- The measure of structural similarity (``measure.compare_ssim``) is now
  n-dimensional and supports color channels as well.
- We fixed an issue related to incorrect propagation in plateaus in
  ``segmentation.watershed``

Documentation:

- New organization of gallery of examples in sections
- More frequent (automated) updates of online documentation

API Changes
-----------

- ``equalize_adapthist`` now takes a ``kernel_size`` keyword argument,
  replacing  the ``ntiles_*`` arguments.
- The functions ``blob_dog``, ``blob_log`` and ``blob_doh`` now return
  float  arrays instead of integer arrays.
- ``transform.integrate`` now takes lists of tuples instead of integers
  to define the window over which to integrate.
- `reverse_map` parameter in `skimage.transform.warp` has been removed.
- `enforce_connectivity` in `skimage.segmentation.slic` defaults to ``True``.
- `skimage.measure.fit.BaseModel._params`,
  `skimage.transform.ProjectiveTransform._matrix`,
  `skimage.transform.PolynomialTransform._params`,
  `skimage.transform.PiecewiseAffineTransform.affines_*` attributes
  have been removed.
- `skimage.filters.denoise_*` have moved to `skimage.restoration.denoise_*`.

Deprecations
------------

- ``filters.gaussian_filter`` has been renamed ``filters.gaussian``
- ``filters.gabor_filter`` has been renamed ``filters.gabor``
- ``restoration.nl_means_denoising`` has been renamed
  ``restoration.denoise_nl_means``
- ``measure.LineModel`` was deprecated in favor of ``measure.LineModelND``
- ``measure.structural_similarity`` has been renamed
  ``measure.compare_ssim``
- ``data.lena`` has been deprecated, and gallery examples use instead the
  ``data.astronaut()`` picture.

Contributors to this release
----------------------------
(Listed alphabetically by last name)

- K.-Michael Aye
- Connelly Barnes
- Sumit Binnani
- Vighnesh Birodkar
- François Boulogne
- Matthew Brett
- Steven Brown
- Arnaud De Bruecker
- Olivier Debeir
- Charles Deledalle
- endolith
- Eric Lubeck
- Victor Escorcia
- Ivo Flipse
- Joel Frederico
- Cédric Gilon
- Christoph Gohlke
- Korijn van Golen
- Emmanuelle Gouillart
- J. Goutin
- Blake Griffith
- M. Hawker
- Jonathan Helmus
- Dhruv Jawali
- Lee Kamentsky
- Kevin Keraudren
- Julius Bier Kirkegaard
- David Koeller
- Gustav Larsson
- Gregory R. Lee
- Guillaume Lemaitre
- Benny Lichtner
- Himanshu Mishra
- Juan Nunez-Iglesias
- Ömer Özak
- Leena P.
- Michael Pacer
- Daniil Pakhomov
- David Perez-Suarez
- Egor Panfilov
- David PS
- Sergio Pascual
- Ariel Rokem
- Nicolas Rougier
- Christian Sachs
- Kshitij Saraogi
- Martin Savc
- Johannes Schönberger
- Arve Seljebu
- Tim Sheerman-Chase
- Scott Sievert
- Steven Silvester
- Alexandre Fioravante de Siqueira
- Daichi Suzuo
- Noah Trebesch
- Pratap Vardhan
- Gael Varoquaux
- Stefan van der Walt
- Joshua Warner
- Josh Warner
- Warren Weckesser
- Daniel Wennberg
- John Wiggins
- Robin Wilson
- Olivia Wilson

Announcement: scikit-image 0.19.0
=================================

We're happy to announce the release of scikit-image v0.19.0!

scikit-image is an image processing toolbox for SciPy that includes algorithms
for segmentation, geometric transformations, color space manipulation,
analysis, filtering, morphology, feature detection, and more.

For more information, examples, and documentation, please visit our website:

https://scikit-image.org

A highlight of this release is the addition of the popular scale-invariant
feature transform (SIFT) feature detector and descriptor. This release also
introduces a perceptual blur metric, new pixel graph algorithms, and most
functions now operate in single-precision when single-precision inputs are
provided. Many other bug fixes, enhancements and performance improvements are
detailed below.

A significant change in this release is in the treatment of multichannel
images. The existing ``multichannel`` argument to functions has been deprecated
in favor of a new ``channel_axis`` argument. ``channel_axis`` can be used to
specify which axis of an array contains channel information (with
``channel_axis=None`` indicating a grayscale image).

scikit-image now uses "lazy loading", which enables users to access the
functions from all ``skimage`` submodules without the overhead of eagerly
importing all submodules. As a concrete example, after calling "import skimage"
a user can directly call a function such as ``skimage.transform.warp`` whereas
previously it would have been required to first "import skimage.transform".

An exciting change on the development side is the introduction of support for
Pythran as an alternative to Cython for generating compiled code. We plan to
keep Cython support as well going forward, so developers are free to use either
one as appropriate. For those curious about Pythran, a good overview was given
in the SciPy 2021 presentation, "Building SciPy Kernels with Pythran" (https://www.youtube.com/watch?v=6a9D9WL6ZjQ).


New Features
------------

- Added support for processing images with channels located along any array
  axis. This is in contrast to previous releases where channels were required
  to be the last axis of an image. See more info on the new ``channel_axis``
  argument under the API section of the release notes.
- A no-reference measure of perceptual blur was added
  (``skimage.measure.blur_effect``).
- Non-local means (``skimage.restoration.denoise_nl_means``) now supports
  3D multichannel, 4D and 4D multichannel data when ``fast_mode=True``.
- An n-dimensional Fourier-domain Butterworth filter
  (``skimage.filters.butterworth``) was added.
- Color conversion functions now have a new ``channel_axis`` keyword argument
  that allows specification of which axis of an array corresponds to channels.
  For backwards compatibility, this parameter defaults to ``channel_axis=-1``,
  indicating that channels are along the last axis.
- Added a new keyword only parameter ``random_state`` to
  ``morphology.medial_axis`` and ``restoration.unsupervised_wiener``.
- Seeding random number generators will not give the same results as the
  underlying generator was updated to use ``numpy.random.Generator``.
- Added ``saturation`` parameter to ``skimage.color.label2rgb``
- Added normalized mutual information metric
  ``skimage.metrics.normalized_mutual_information``
- threshold_local now supports n-dimensional inputs and anisotropic block_size
- New ``skimage.util.label_points`` function for assigning labels to points.
- Added nD support to several geometric transform classes
- Added ``skimage.metrics.hausdorff_pair`` to find points separated by the
  Hausdorff distance.
- Additional colorspace ``illuminants`` and ``observers`` parameter options
  were added to ``skimage.color.lab2rgb``, ``skimage.color.rgb2lab``,
  ``skimage.color.xyz2lab``, ``skimage.color.lab2xyz``,
  ``skimage.color.xyz2luv`` and ``skimage.color.luv2xyz``.
- ``skimage.filters.threshold_multiotsu`` has a new ``hist`` keyword argument
  to allow use with a user-supplied histogram. (gh-5543)
- ``skimage.restoration.denoise_bilateral`` added support for images containing
  negative values. (gh-5527)
- The ``skimage.feature`` functions ``blob_dog``, ``blob_doh`` and ``blob_log``
  now support a ``threshold_rel`` keyword argument that can be used to specify
  a relative threshold (in range [0, 1]) rather than an absolute one. (gh-5517)
- Implement lazy submodule importing (gh-5101)
- Implement weighted estimation of geometric transform matrices (gh-5601)
- Added new pixel graph algorithms in ``skimage.graph``:
  ``pixel_graph`` generates a graph (network) of pixels
  according to their adjacency, and ``central_pixel`` finds
  the geodesic center of the pixels. (gh-5602)
- scikit-image now supports use of Pythran in contributed code. (gh-3226)


Documentation
-------------

- A new doc tutorial presenting a 3D biomedical imaging example has been added
  to the gallery (gh-4946). The technical content benefited from conversations
  with Genevieve Buckley, Kevin Mader, and Volker Hilsenstein.
- New gallery example for 3D structure tensor.
- New gallery example displaying a 3D dataset.
- Extended rolling ball example with ECG data (1D).
- The stain unmixing gallery example was fixed and now displays proper
  separation of the stains.
- Documentation has been added to the contributing notes about how to submit a
  gallery example.
- Autoformat docstrings in morphology.
- Display plotly figures from gallery example even when running script at CLI.
- Single out docs-only PRs in review process.
- Use matplotlib's infinite axline to demonstrate hough transform.
- Clarify disk documentation inconsistency regarding 'shape'.
- docs: fix simple typo, convertions -> conversions.
- Fixes to linspace in example.
- Minor fixes to Hough line transform code and examples.
- Added 1/2 pixel bounds to extent of displayed images in several examples.
- Add release step on github to RELEASE.txt.
- Remove reference to opencv in threshold_local documentation.
- Update structure_tensor docstring to include per-axis sigma.
- Fix typo in _shared/utils.py docs.
- Proofread and crosslink examples with immunohistochemistry image.
- Spelling correction: witch -> which.
- Mention possible filters in radon_transform -> filtered-back-projection
- Fix dtype info in documentation for watershed.
- Proofread gallery example for Radon transform.
- Use internal function for noise + clarify code in Canny example.
- Make more comprehensive 'see also' sections in filters.
- Specify the release note version instead of the misleading `latest`.
- Remove misleading comment in ``plot_thresholding.py`` example.
- Fix sphinx layout to make the search engine work with recent sphinx versions.
- Draw node IDs in RAG example.
- Update sigma_color description in denoise_bilateral.
- Update intersphinx fallback inventories + add matplotlib fallback inventory.
- Fix numpy deprecation in ``plot_local_equalize.py``.
- Rename ``label`` variable in ``plot_regionprops.py`` to circumvent link issue
  in docs.
- Avoid duplicate API documentation for ImageViewer, CollectionViewer.
- Fix 'blog_dog' typo in ``gaussian`` docs.
- Update reference link documentation in the ``adjust_sigmoid`` function.
- Fix reference to multiscale_basic_features in TrainableSegmenter.
- Slight ``shape_index`` docstring modification to specify 2D array.
- Add stitching gallery example (gh-5365)
- Add draft SKIP3: transition to scikit-image 1.0 (gh-5475)
- Mention commit messages in the contribution guidelines. (gh-5504)
- Fix and standardize docstrings for blob detection functions. (gh-5547)
- Update the User Guide to reflect usage of ``channel_axis`` rather than
  ``multichannel``. (gh-5554)
- Update the user guide to use channel_axis rather than multichannel (gh-5556)
- Add hyperlinks to referenced documentation places. (gh-5560)
- Update branching instructions to change the location of the pooch repo.
  (gh-5565)
- Add Notes and References section to the Cascade class docstring. (gh-5568)
- Clarify 2D vs nD in skimage.feature.corner docstrings (gh-5569)
- Fix math formulas in plot_swirl.py example. (gh-5574)
- Update references in texture feature detectors docstrings (gh-5578)
- Update mailing list location to discuss.scientific-python.org forum (gh-5951)
- DOC: Fix docstring in rescale_intensity() (gh-5964)
- Fix slic documentation (gh-5975)
- Update docstring for dilation, which is now nD. (gh-5978)
- Change stitching gallery example thumbnail (gh-5985)
- Add circle and disk to glossary.md (gh-5590)
- Update pixel graphs example (gh-5991)
- Separate entries that have the same description in glossary.md (gh-5592)
- Do not use space before colon in directive name (gh-6002)


Improvements
------------

- Many more functions throughout the library now have single precision
  (float32) support.
- Biharmonic  inpainting (``skimage.restoration.inpaint_biharmonic``) was
  refactored and is orders of magnitude faster than before.
- Salt-and-pepper noise generation with ``skimage.util.random_noise`` is now
  faster.
- The performance of the SLIC superpixels algorithm
  (``skimage.segmentation.slice``) was improved for the case where a mask
  is supplied by the user (gh-4903). The specific superpixels produced by
  masked SLIC will not be identical to those produced by prior releases.
- ``exposure.adjust_gamma`` has been accelerated for ``uint8`` images thanks to
  a LUT (gh-4966).
- ``measure.label`` has been accelerated for boolean input images, by using
  ``scipy.ndimage``'s implementation for this case (gh-4945).
- ``util.apply_parallel`` now works with multichannel data (gh-4927).
- ``skimage.feature.peak_local_max`` supports now any Minkowski distance.
- Fast, non-Cython implementation for ``skimage.filters.correlate_sparse``.
- For efficiency, the histogram is now precomputed within
  ``skimage.filters.try_all_threshold``.
- Faster ``skimage.filters.find_local_max`` when given a finite ``num_peaks``.
- All filters in the ``skimage.filters.rank`` module now release the GIL,
  enabling multithreaded use.
- ``skimage.restoration.denoise_tv_bregman`` and
  ``skimage.restoration.denoise_bilateral`` now release the GIL, enabling
  multithreaded use.
- A ``skimage.color.label2rgb`` performance regression was addressed.
- Improve numerical precision in ``CircleModel.estimate``. (gh-5190)
- Add default keyword argument values to
  ``skimage.restoration.denoise_tv_bregman``, ``skimage.measure.block_reduce``,
  and ``skimage.filters.threshold_local``. (gh-5454)
- Make matplotlib an optional dependency (gh-5990)
- single precision support in skimage.filters (gh-5354)
- Support nD images and labels in label2rgb (gh-5550)
- Regionprops table performance refactor (gh-5576)
- add regionprops benchmark script (gh-5579)
- remove use of apply_along_axes from greycomatrix & greycoprops (gh-5580)
- refactor gabor_kernel for efficiency (gh-5582)
- remove need for channel_as_last_axis decorator in skimage.filters (gh-5584)
- replace use of scipy.ndimage.gaussian_filter with skimage.filters.gaussian
  (gh-5872)
- add channel_axis argument to quickshift (gh-5987)


API Changes
-----------

- The ``multichannel`` boolean argument has been deprecated. All functions with
  multichannel support now use an integer ``channel_axis`` to specify which
  axis corresponds to channels. Setting ``channel_axis`` to None is used to
  indicate that the image is grayscale. Specifically, existing code with
  ``multichannel=True`` should be updated to use ``channel_axis=-1`` and code
  with ``multichannel=False`` should now specify ``channel_axis=None``.
- Most functions now return float32 images when the input has float32 dtype.
- A default value has been added to ``measure.find_contours``, corresponding to
  the half distance between the min and max values of the image
  (gh-4862).
- ``data.cat`` has been introduced as an alias of ``data.chelsea`` for a more
  descriptive name.
- The ``level`` parameter of ``measure.find_contours`` is now a keyword
  argument, with a default value set to ``(max(image) - min(image)) / 2``.
- ``p_norm`` argument was added to ``skimage.feature.peak_local_max``
  to add support for Minkowski distances.
- ``skimage.transforms.integral_image`` now promotes floating point inputs to
  double precision by default (for accuracy). A new ``dtype`` keyword argument
  can be used to override this behavior when desired.
- Color conversion functions now have a new ``channel_axis`` keyword argument
  (see **New Features** section).
- SLIC superpixel segmentation outputs may differ from previous versions for
  data that was not already scaled to [0, 1] range. There is now an automatic
  internal rescaling of the input to [0, 1] so that the ``compactness``
  parameter has an effect that is independent of the input image's scaling.
- A bug fix to the phase normalization applied within
  ``skimage.register.phase_cross_correlation`` may result in a different result
  as compared to prior releases. The prior behavior of "unnormalized" cross
  correlation is still available by explicitly setting ``normalization=None``.
  There is no change to the masked cross-correlation case, which uses a
  different algorithm.


Bugfixes
--------

- Input ``labels`` argument renumbering in ``skimage.feature.peak_local_max``
  is avoided (gh-5047).
- fix clip bug in resize when anti_aliasing is applied (gh-5202)
- Nonzero values at the image edge are no longer incorrectly marked as a
  boundary when using ``find_bounaries`` with mode='subpixel' (gh-5447).
- Fix return dtype of ``_label2rgb_avg`` function.
- Ensure ``skimage.color.separate_stains`` does not return negative values.
- Prevent integer overflow in ``EllipseModel``.
- Fixed off-by one error in pixel bins in Hough line transform,
  ``skimage.transform.hough_line``.
- Handle 1D arrays properly in ``skimage.filters.gaussian``.
- Fix Laplacian matrix size bug in ``skimage.segmentation.random_walker``.
- Regionprops table (``skimage.measure.regionprops_table``) dtype bugfix.
- Fix ``skimage.transform.rescale`` when using a small scale factor.
- Fix ``skimage.measure.label`` segfault.
- Watershed (``skimage.segmentation.watershed``): consider connectivity when
  calculating markers.
- Fix ``skimage.transform.warp`` output dtype when order=0.
- Fix multichannel ``intensity_image`` extra_properties in regionprops.
- Fix error message for ``skimage.metric.structural_similarity`` when image is
  too small.
- Do not mark image edges in 'subpixel' mode of
  ``skimage.segmentation.find_boundaries``.
- Fix behavior of ``skimage.exposure.is_low_contrast`` for boolean inputs.
- Fix wrong syntax for the string argument of ValueError in
  ``skimage.metric.structural_similarity`` .
- Fixed NaN issue in ``skimage.filters.threshold_otsu``.
- Fix ``skimage.feature.blob_dog`` docstring example and normalization.
- Fix uint8 overflow in ``skimage.exposure.adjust_gamma``.
- Work with pooch 1.5.0 for fetching data (gh-5529).
- The ``offsets`` attribute of ``skimage.graph.MCP`` is now public. (gh-5547)
- Fix io.imread behavior with pathlib.Path inputs (gh-5543)
- Make scikit-image imports from Pooch, compatible with pooch >= 1.5.0.
  (gh-5529)
- Fix several broken doctests and restore doctesting on GitHub Actions.
  (gh-5505)
- Fix broken doctests in ``skimage.exposure.histogram`` and
  ``skimage.measure.regionprops_table``. (gh-5522)
- Rescale image consistently during SLIC superpixel segmentation. (gh-5518)
- Correct phase correlation in ``skimage.register.phase_cross_correlation``.
  (gh-5461)
- Fix hidden attribute 'offsets' in skimage.graph.MCP (gh-5551)
- fix phase_cross_correlation for 3D with reference masks (gh-5559)
- fix return shape of blob_log and blob_dog when no peaks are found (gh-5567)
- Fix find contours key error (gh-5577)
- Refactor measure.ransac and add warning when the estimated model is not valid
  (gh-5583)
- Restore integer image rescaling for edge filters (gh-5589)
- trainable_segmentation: re-raise in error case (gh-5600)
- allow regionprops_table to be called with deprecated property names (gh-5908)
- Fix weight calculation in fast mode of non-local means (gh-5923)
- fix for #5948: lower boundary 1 for kernel_size in equalize_adapthist
  (gh-5949)
- convert pathlib.Path to str in imsave (gh-5971)
- Fix slic spacing (gh-5974)
- Add small regularization to avoid zero-division in richardson_lucy (gh-5976)
- Fix benchmark suite (watershed function was moved) (gh-5982)
- catch QhullError and return empty array (``convex_hull``) (gh-6008)
- add property getters for all newly deprecated regionprops names (gh-6000)
- Fix the estimation of ellipsoid axis lengths in the 3D case (gh-6013)
- Fix peak local max segfault (gh-6035)
- Avoid circular import errors when EAGER_IMPORT=1 (gh-6042)
- remove all use of the deprecated distutils package (gh-6044)


Deprecations
------------

Completed deprecations from prior releases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- In ``measure.label``, the deprecated ``neighbors`` parameter has been
  removed (use ``connectivity`` instead).
- The deprecated ``skimage.color.rgb2grey`` and ``skimage.color.grey2rgb``
  functions have been removed (use ``skimage.color.rgb2gray`` and
  ``skimage.color.gray2rgb`` instead).
- ``skimage.color.rgb2gray`` no longer allows grayscale or RGBA inputs.
- The deprecated ``alpha`` parameter of ``skimage.color.gray2rgb`` has now been
  removed. Use ``skimage.color.gray2rgba`` for conversion to RGBA.
- Attempting to warp a boolean image with ``order > 0`` now raises a
  ValueError.
- When warping or rescaling boolean images, setting ``anti-aliasing=True`` will
  raise a ValueError.
- The ``bg_label`` parameter of ``skimage.color.label2rgb`` is now 0.
- The deprecated ``filter`` parameter of ``skimage.transform.iradon`` has now
  been removed (use ``filter_name`` instead).
- The deprecated ``skimage.draw.circle`` function has been removed (use
  ``skimage.draw.disk`` instead).
- The deprecated ``skimage.feature.register_translation`` function has
  been removed (use ``skimage.registration.phase_cross_correlation`` instead).
- The deprecated ``skimage.feature.masked_register_translation`` function has
  been removed (use ``skimage.registration.phase_cross_correlation`` instead).
- The deprecated ``skimage.measure.marching_cubes_classic`` function has
  been removed (use ``skimage.measure.marching_cubes`` instead).
- The deprecated ``skimage.measure.marching_cubes_lewiner`` function has
  been removed (use ``skimage.measure.marching_cubes`` instead).
- The deprecated ``skimage.segmentation.circle_level_set`` function has been
  removed (use ``skimage.segmentation.disk_level_set`` instead).
- The deprecated ``inplace`` parameter of ``skimage.morphology.flood_fill``
- The deprecated ``skimage.util.pad`` function has been removed (use
  ``numpy.pad`` instead).
  been removed (use ``in_place`` instead).
- The default ``mode`` in ``skimage.filters.hessian`` is now
  ``'reflect'``.
- The default boundary ``mode`` in ``skimage.filters.sato`` is now
  ``'reflect'``.
- The default boundary ``mode`` in ``skimage.measure.profile_line`` is now
  ``'reflect'``.
- The default value of ``preserve_range`` in
  ``skimage.restoration.denoise_nl_means`` is now False.
- The default value of ``start_label`` in ``skimage.segmentation.slic`` is now
  1.

Newly introduced deprecations:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The ``multichannel`` argument is now deprecated throughout the library and
  will be removed in 1.0. The new ``channel_axis`` argument should be used
  instead. Existing code with ``multichannel=True`` should be updated to use
  ``channel_axis=-1`` and code with ``multichannel=False`` should now specify
  ``channel_axis=None``.
- ``skimage.feature.greycomatrix`` and ``skimage.feature.greycoprops`` are
  deprecated in favor of ``skimage.feature.graycomatrix`` and
  ``skimage.feature.graycoprops``.
- The ``skimage.morphology.grey`` module has been renamed
  ``skimage.morphology.gray``. The old name is deprecated.
- The ``skimage.morphology.greyreconstruct`` module has been renamed
  ``skimage.morphology.grayreconstruct``. The old name is deprecated.
- see **API Changes** section regarding functions with deprecated argument
  names related to the number of iterations. ``num_iterations`` and
  ``max_num_iter`` are now used throughout the library.
- see **API Changes** section on deprecation of the ``selem`` argument in favor
  of ``footprint`` throughout the library
- Deprecate ``in_place`` in favor of the use of an explicit ``out`` argument
  in ``skimage.morphology.remove_small_objects``,
  ``skimage.morphology.remove_small_holes`` and
  ``skimage.segmentation.clear_border``
- The ``input`` argument of ``skimage.measure.label`` has been renamed
  ``label_image``. The old name is deprecated.
- standardize on ``num_iter`` for paramters describing the number of iterations
  and ``max_num_iter`` for parameters specifying an iteration limit. Functions
  where the old argument names have now been deprecated are::

    skimage.filters.threshold_minimum
    skimage.morphology.thin
    skimage.restoration.denoise_tv_bregman
    skimage.restoration.richardson_lucy
    skimage.segmentation.active_contour
    skimage.segmentation.chan_vese
    skimage.segmentation.morphological_chan_vese
    skimage.segmentation.morphological_geodesic_active_contour
    skimage.segmentation.slic

- The names of several parameters in ``skimage.measure.regionprops`` have been
  updated so that properties are better grouped by the first word(s) of the
  name. The old names will continue to work for backwards compatibility.
  The specific names that were updated are::

    ============================ ============================
    Old Name                     New Name
    ============================ ============================
    max_intensity                intensity_max
    mean_intensity               intensity_mean
    min_intensity                intensity_min

    bbox_area                    area_bbox
    convex_area                  area_convex
    filled_area                  area_filled

    convex_image                 image_convex
    filled_image                 image_filled
    intensity_image              image_intensity

    local_centroid               centroid_local
    weighted_centroid            centroid_weighted
    weighted_local_centroid      centroid_weighted_local

    major_axis_length            axis_major_length
    minor_axis_length            axis_minor_length

    weighted_moments             moments_weighted
    weighted_moments_central     moments_weighted_central
    weighted_moments_hu          moments_weighted_hu
    weighted_moments_normalized  moments_weighted_normalized

    equivalent_diameter          equivalent_diameter_area
    ============================ ============================

- The ``selem`` argument has been renamed to ``footprint`` throughout the
  library. The ``selem`` argument is now deprecated.


Development process
-------------------

- Test setup and teardown functions added to allow raising an error on any
  uncaught warnings via ``SKIMAGE_TEST_STRICT_WARNINGS_GLOBAL`` environment
  variable.
- Increase automation in release process.
- Release wheels before source
- update minimum supported Matplotlib, NumPy, SciPy and Pillow
- Pin pillow to !=8.3.0
- Rename `master` to `main` throughout
- Ensure that README.txt has write permissions for subsequent imports.
- Run face classification gallery example with a single thread
- Enable pip and skimage.data caching on Azure
- Fix CircleCI and Azure CI caching.
- Address Cython warnings.
- Disable calls to plotly.io.show when running on Azure.
- Remove legacy Travis-CI scripts and update contributor documentation
  accordingly.
- Increase cibuildwheel verbosity.
- Update pip during dev environment installation.
- Add benchmark checks to CI.
- Resolve stochastic rank filter test failures on CI.
- Ensure that README.txt has write permissions for subsequent imports.
- Decorators for helping with the transition between the keyword argument
  multichannel and channel_axis.
- Add missing import in lch2lab docstring example (gh-5998)
- Prefer importing build_py and sdist from setuptools (gh-6007)
- Reintroduce skimage.test utility (gh-5909)


Other Updates
-------------
- Refactor np.random.x to use np.random.Generator.
- Avoid warnings about use of deprecated ``scipy.linalg.pinv2``.
- Simplify resize implementation using new SciPy 1.6 zoom option.
- Fix duplicate test function names in ``test_unsharp_mask.py``.
- Benchmarks: ``fix ResizeLocalMeanSuite.time_resize_local_mean`` signature.
- Prefer use of new-style NumPy random API in tests (gh-5450)
- Add fixture enforcing SimpleITK I/O in test_simpleitk.py (gh-5526)
- MNT: Remove unused stat import from skimage data (gh-5566)
- MAINT: Remove unused imports (gh-5595)
- MAINT: Refactor duplicated tests, remove unnecessary assignments and
  variables (gh-5596)
- Remove obsolete lazy import (gh-5992)
- Lazily load data_dir into the top-level namespace (gh-5996)
- Update scipy requirement to 1.4.1 and use scipy.fft instead of scipy.fftpack
  (gh-5999)
- Remove lines generating Requires metadata (gh-6017)
- Update wheel builds to include Python 3.10 (gh-6021)
- Update pyproject.toml to handle Python 3.10 and Apple arm64 (gh-6022)
- Add python 3.10 test runs on GitHub Actions and Appveyor (gh-6027)
- Pin sphinx to <4.3 until new sphinx-gallery release is available (gh-6029)
- Relax a couple of equality tests causing i686 test failures on cibuildwheel
  (gh-6031)
- Avoid matplotlib import overhead during 'import skimage' (gh-6032)
- Update sphinx gallery pin (gh-6034)


Contributors to this release
----------------------------


80 authors added to this release [alphabetical by first name or login]
----------------------------------------------------------------------
- Abhinavmishra8960 (Abhinavmishra8960)
- abouysso
- Alessia Marcolini
- Alex Brooks
- Alexandre de Siqueira
- Andres Fernandez
- Andrew Hurlbatt
- andrewnags (andrewnags)
- Antoine Bierret
- BMaster123 (BMaster123)
- Boaz Mohar
- Bozhidar Karaargirov
- Carlos Andrés Álvarez Restrepo
- Christoph Gohlke
- Christoph Sommer
- Clement Ng
- cmarasinou
- Cris Luengo
- David Manthey
- Devanshu Shah
- Dhiraj Kumar Sah
- divyank agarwal
- Egor Panfilov
- Emmanuelle Gouillart
- Erik Reed
- erykoff (erykoff)
- Fabian Schneider
- Felipe Gutierrez-Barragan
- François Boulogne
- Fred Bunt
- Fukai Yohsuke
- Gregory R. Lee
- Hari Prasad
- Harish Venkataraman
- Harshit Dixit
- Ian Hunt-Isaak
- Jaime Rodríguez-Guerra
- Jan-Hendrik Müller
- Janakarajan Natarajan
- Jenny Vo
- john lee
- Jonathan Striebel
- Joseph Fox-Rabinovitz
- Juan Antonio Barragan Noguera
- Juan Nunez-Iglesias
- Julien Jerphanion
- Jurneo
- klaussfreire (klaussfreire)
- Larkinnjm1 (Larkinnjm1)
- Lars Grüter
- Mads Dyrmann
- Marianne Corvellec
- Marios Achilias
- Mark Boer
- Mark Harfouche
- Matthias Bussonnier
- Mauro Silberberg
- Max Frei
- michalkrawczyk (michalkrawczyk)
- Niels Cautaerts
- Pamphile ROY
- Pradyumna Rahul
- R
- Raphael
- Riadh Fezzani
- Robert Haase
- Sebastian Gonzalez Tirado
- Sebastián Vanrell
- serge-sans-paille (serge-sans-paille)
- Stefan van der Walt
- t.ae
- that1solodev (Xyno18)
- Thomas Walter
- Tim Gates
- Tom Flux
- Vinicius D. Cerutti
- Volker Hilsenstein
- WeiChungChang
- yacth
- Yash-10 (Yash-10)

63 reviewers added to this release [alphabetical by first name or login]
------------------------------------------------------------------------
- Abhinavmishra8960
- Alessia Marcolini
- Alex Brooks
- Alexandre de Siqueira
- Andres Fernandez
- Andrew Hurlbatt
- andrewnags
- BMaster123
- Boaz Mohar
- Carlos Andrés Álvarez Restrepo
- Clement Ng
- Cris Luengo
- Dan Schult
- David Manthey
- Egor Panfilov
- Emmanuelle Gouillart
- erykoff
- Fabian Schneider
- Felipe Gutierrez-Barragan
- François Boulogne
- Fukai Yohsuke
- Genevieve Buckley
- Gregory R. Lee
- Jan Eglinger
- Jan-Hendrik Müller
- Janakarajan Natarajan
- Jarrod Millman
- Jirka Borovec
- Joan Massich
- Johannes Schönberger
- john lee
- Jon Crall
- Joseph Fox-Rabinovitz
- Josh Warner
- Juan Nunez-Iglesias
- Julien Jerphanion
- Kenneth Hoste
- klaussfreire
- Larkinnjm1
- Lars Grüter
- Marianne Corvellec
- Mark Boer
- Mark Harfouche
- Matthias Bussonnier
- Max Frei
- michalkrawczyk
- Niels Cautaerts
- Pamphile ROY
- Pomax
- R
- Raphael
- Riadh Fezzani
- Robert Kern
- Ross Barnowski
- Sebastian Berg
- Sebastian Gonzalez Tirado
- Sebastian Wallkötter
- serge-sans-paille
- Stefan van der Walt
- t.ae
- Vinicius D. Cerutti
- Volker Hilsenstein
- Yash-10
