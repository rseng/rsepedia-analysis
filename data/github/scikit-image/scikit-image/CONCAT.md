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
