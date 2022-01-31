# kalepy: Kernel Density Estimation and Sampling

[![Build Status](https://travis-ci.org/lzkelley/kalepy.svg?branch=master)](https://travis-ci.org/lzkelley/kalepy)
[![codecov](https://codecov.io/gh/lzkelley/kalepy/branch/master/graph/badge.svg)](https://codecov.io/gh/lzkelley/kalepy)
[![Documentation Status](https://readthedocs.org/projects/kalepy/badge/?version=latest)](https://kalepy.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02784/status.svg)](https://doi.org/10.21105/joss.02784)
[![DOI](https://zenodo.org/badge/187267055.svg)](https://zenodo.org/badge/latestdoi/187267055)

![kalepy animated logo](https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/logo_anim_small.gif)

This package performs KDE operations on multidimensional data to: **1) calculate estimated PDFs** (probability distribution functions), and **2) resample new data** from those PDFs.

## Documentation

A number of examples (also used for continuous integration testing) are included in [the package notebooks](https://github.com/lzkelley/kalepy/tree/master/notebooks).  Some background information and references are included in [the JOSS paper](https://joss.theoj.org/papers/10.21105/joss.02784).

Full documentation is available on [kalepy.readthedocs.io](https://kalepy.readthedocs.io/en/latest/).

## README Contents

- [Installation](#Installation)
- Quickstart
    - [Basic Usage](#Basic-Usage)
    - [Fancy Usage](#Fancy-Usage)
- [Development & Contributions](#Development-&-Contributions)
- [Attribution (citation)](#Attribution)


## Installation

#### from pypi (i.e. via pip)

```bash
pip install kalepy
```

#### from source (e.g. for development)

```bash
git clone https://github.com/lzkelley/kalepy.git
pip install -e kalepy/
```

In this case the package can easily be updated by changing into the source directory, pulling, and rebuilding:

```bash
cd kalepy
git pull
pip install -e .
# Optional: run unit tests (using the `nosetests` package)
nosetests
```
# kalepy: Kernel Density Estimation and Sampling

[![Build Status](https://travis-ci.org/lzkelley/kalepy.svg?branch=master)](https://travis-ci.org/lzkelley/kalepy)
[![codecov](https://codecov.io/gh/lzkelley/kalepy/branch/master/graph/badge.svg)](https://codecov.io/gh/lzkelley/kalepy)
[![Documentation Status](https://readthedocs.org/projects/kalepy/badge/?version=latest)](https://kalepy.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02784/status.svg)](https://doi.org/10.21105/joss.02784)
[![DOI](https://zenodo.org/badge/187267055.svg)](https://zenodo.org/badge/latestdoi/187267055)

![kalepy animated logo](https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/logo_anim_small.gif)

This package performs KDE operations on multidimensional data to: **1) calculate estimated PDFs** (probability distribution functions), and **2) resample new data** from those PDFs.

## Documentation

A number of examples (also used for continuous integration testing) are included in [the package notebooks](https://github.com/lzkelley/kalepy/tree/master/notebooks).  Some background information and references are included in [the JOSS paper](https://joss.theoj.org/papers/10.21105/joss.02784).

Full documentation is available on [kalepy.readthedocs.io](https://kalepy.readthedocs.io/en/latest/).

## README Contents

- [Installation](#Installation)
- Quickstart
    - [Basic Usage](#Basic-Usage)
    - [Fancy Usage](#Fancy-Usage)
- [Development & Contributions](#Development-&-Contributions)
- [Attribution (citation)](#Attribution)


## Installation

#### from pypi (i.e. via pip)

```bash
pip install kalepy
```

#### from source (e.g. for development)

```bash
git clone https://github.com/lzkelley/kalepy.git
pip install -e kalepy/
```

In this case the package can easily be updated by changing into the source directory, pulling, and rebuilding:

```bash
cd kalepy
git pull
pip install -e .
# Optional: run unit tests (using the `nosetests` package)
nosetests
```


# Basic Usage


```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import kalepy as kale

from kalepy.plot import nbshow
```

Generate some random data, and its corresponding distribution function


```python
NUM = int(1e4)
np.random.seed(12345)
# Combine data from two different PDFs
_d1 = np.random.normal(4.0, 1.0, NUM)
_d2 = np.random.lognormal(0, 0.5, size=NUM)
data = np.concatenate([_d1, _d2])

# Calculate the "true" distribution
xx = np.linspace(0.0, 7.0, 100)[1:]
yy = 0.5*np.exp(-(xx - 4.0)**2/2) / np.sqrt(2*np.pi)
yy += 0.5 * np.exp(-np.log(xx)**2/(2*0.5**2)) / (0.5*xx*np.sqrt(2*np.pi))
```

### Plotting Smooth Distributions


```python
# Reconstruct the probability-density based on the given data points.
points, density = kale.density(data, probability=True)

# Plot the PDF
plt.plot(points, density, 'k-', lw=2.0, alpha=0.8, label='KDE')

# Plot the "true" PDF
plt.plot(xx, yy, 'r--', alpha=0.4, lw=3.0, label='truth')

# Plot the standard, histogram density estimate
plt.hist(data, density=True, histtype='step', lw=2.0, alpha=0.5, label='hist')

plt.legend()
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_8_0.png)


### resampling: constructing statistically similar values

Draw a new sample of data-points from the KDE PDF


```python
# Draw new samples from the KDE reconstructed PDF
samples = kale.resample(data)

# Plot new samples
plt.hist(samples, density=True, label='new samples', alpha=0.5, color='0.65', edgecolor='b')
# Plot the old samples
plt.hist(data, density=True, histtype='step', lw=2.0, alpha=0.5, color='r', label='input data')

# Plot the KDE reconstructed PDF
plt.plot(points, density, 'k-', lw=2.0, alpha=0.8, label='KDE')

plt.legend()
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_11_0.png)


### Multivariate Distributions


```python
reload(kale.plot)

# Load some random-ish three-dimensional data
np.random.seed(9485)
data = kale.utils._random_data_3d_02(num=3e3)

# Construct a KDE
kde = kale.KDE(data)

# Construct new data by resampling from the KDE
resamp = kde.resample(size=1e3)

# Plot the data and distributions using the builtin `kalepy.corner` plot
corner, h1 = kale.corner(kde, quantiles=[0.5, 0.9])
h2 = corner.clean(resamp, quantiles=[0.5, 0.9], dist2d=dict(median=False), ls='--')

corner.legend([h1, h2], ['input data', 'new samples'])

nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_13_0.png)



```python
# Resample the data (default output is the same size as the input data)
samples = kde.resample()


# ---- Plot the input data compared to the resampled data ----

fig, axes = plt.subplots(figsize=[16, 4], ncols=kde.ndim)

for ii, ax in enumerate(axes):
    # Calculate and plot PDF for `ii`th parameter (i.e. data dimension `ii`)
    xx, yy = kde.density(params=ii, probability=True)
    ax.plot(xx, yy, 'k--', label='KDE', lw=2.0, alpha=0.5)
    # Draw histograms of original and newly resampled datasets
    *_, h1 = ax.hist(data[ii], histtype='step', density=True, lw=2.0, label='input')
    *_, h2 = ax.hist(samples[ii], histtype='step', density=True, lw=2.0, label='resample')
    # Add 'kalepy.carpet' plots showing the data points themselves
    kale.carpet(data[ii], ax=ax, color=h1[0].get_facecolor())
    kale.carpet(samples[ii], ax=ax, color=h2[0].get_facecolor(), shift=ax.get_ylim()[0])

axes[0].legend()
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_14_0.png)


# Fancy Usage

### Reflecting Boundaries

What if the distributions you're trying to capture have edges in them, like in a uniform distribution between two bounds?  Here, the KDE chooses 'reflection' locations based on the extrema of the given data.


```python
# Uniform data (edges at -1 and +1)
NDATA = 1e3
np.random.seed(54321)
data = np.random.uniform(-1.0, 1.0, int(NDATA))

# Create a 'carpet' plot of the data
kale.carpet(data, label='data')
# Histogram the data
plt.hist(data, density=True, alpha=0.5, label='hist', color='0.65', edgecolor='k')

# ---- Standard KDE will undershoot just-inside the edges and overshoot outside edges
points, pdf_basic = kale.density(data, probability=True)
plt.plot(points, pdf_basic, 'r--', lw=3.0, alpha=0.5, label='KDE')

# ---- Reflecting KDE keeps probability within the given bounds
# setting `reflect=True` lets the KDE guess the edge locations based on the data extrema
points, pdf_reflect = kale.density(data, reflect=True, probability=True)
plt.plot(points, pdf_reflect, 'b-', lw=2.0, alpha=0.75, label='reflecting KDE')

plt.legend()
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_18_0.png)


Explicit reflection locations can also be provided (in any number of dimensions).


```python
# Construct random data, add an artificial 'edge'
np.random.seed(5142)
edge = 1.0
data = np.random.lognormal(sigma=0.5, size=int(3e3))
data = data[data >= edge]

# Histogram the data, use fixed bin-positions
edges = np.linspace(edge, 4, 20)
plt.hist(data, bins=edges, density=True, alpha=0.5, label='data', color='0.65', edgecolor='k')

# Standard KDE with over & under estimates
points, pdf_basic = kale.density(data, probability=True)
plt.plot(points, pdf_basic, 'r--', lw=4.0, alpha=0.5, label='Basic KDE')

# Reflecting KDE setting the lower-boundary to the known value
#    There is no upper-boundary when `None` is given.
points, pdf_basic = kale.density(data, reflect=[edge, None], probability=True)
plt.plot(points, pdf_basic, 'b-', lw=3.0, alpha=0.5, label='Reflecting KDE')

plt.gca().set_xlim(edge - 0.5, 3)
plt.legend()
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_20_0.png)


### Multivariate Reflection


```python
# Load a predefined dataset that has boundaries at:
#   x: 0.0 on the low-end
#   y: 1.0 on the high-end
data = kale.utils._random_data_2d_03()

# Construct a KDE with the given reflection boundaries given explicitly
kde = kale.KDE(data, reflect=[[0, None], [None, 1]])

# Plot using default settings
kale.corner(kde)

nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_22_0.png)


### Specifying Bandwidths and Kernel Functions


```python
# Load predefined 'random' data
data = kale.utils._random_data_1d_02(num=100)
# Choose a uniform x-spacing for drawing PDFs
xx = np.linspace(-2, 8, 1000)

# ------ Choose the kernel-functions and bandwidths to test -------  #
kernels = ['parabola', 'gaussian', 'box']                            #
bandwidths = [None, 0.9, 0.15]     # `None` means let kalepy choose  #
# -----------------------------------------------------------------  #

ylabels = ['Automatic', 'Course', 'Fine']
fig, axes = plt.subplots(figsize=[16, 10], ncols=len(kernels), nrows=len(bandwidths), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0.2, wspace=0.05)
for (ii, jj), ax in np.ndenumerate(axes):

    # ---- Construct KDE using particular kernel-function and bandwidth ---- #
    kern = kernels[jj]                                                       #
    bw = bandwidths[ii]                                                      #
    kde = kale.KDE(data, kernel=kern, bandwidth=bw)                          #
    # ---------------------------------------------------------------------- #

    # If bandwidth was set to `None`, then the KDE will choose the 'optimal' value
    if bw is None:
        bw = kde.bandwidth[0, 0]

    ax.set_title('{} (bw={:.3f})'.format(kern, bw))
    if jj == 0:
        ax.set_ylabel(ylabels[ii])

    # plot the KDE
    ax.plot(*kde.pdf(points=xx), color='r')
    # plot histogram of the data (same for all panels)
    ax.hist(data, bins='auto', color='b', alpha=0.2, density=True)
    # plot  carpet   of the data (same for all panels)
    kale.carpet(data, ax=ax, color='b')

ax.set(xlim=[-2, 5], ylim=[-0.2, 0.6])
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_24_0.png)


## Resampling

### Using different data `weights`


```python
# Load some random data (and the 'true' PDF, for comparison)
data, truth = kale.utils._random_data_1d_01()

# ---- Resample the same data, using different weightings ---- #
resamp_uni = kale.resample(data, size=1000)                       #
resamp_sqr  = kale.resample(data, weights=data**2, size=1000)      #
resamp_inv = kale.resample(data, weights=data**-1, size=1000)     #
# ------------------------------------------------------------ #


# ---- Plot different distributions ----

# Setup plotting parameters
kw = dict(density=True, histtype='step', lw=2.0, alpha=0.75, bins='auto')

xx, yy = truth
samples = [resamp_inv, resamp_uni, resamp_sqr]
yvals = [yy/xx, yy, yy*xx**2/10]
labels = [r'$\propto X^{-1}$', r'$\propto 1$', r'$\propto X^2$']

plt.figure(figsize=[10, 5])

for ii, (res, yy, lab) in enumerate(zip(samples, yvals, labels)):
    hh, = plt.plot(xx, yy, ls='--', alpha=0.5, lw=2.0)
    col = hh.get_color()
    kale.carpet(res, color=col, shift=-0.1*ii)
    plt.hist(res, color=col, label=lab, **kw)

plt.gca().set(xlim=[-0.5, 6.5])
# Add legend
plt.legend()
# display the figure if this is a notebook
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_27_0.png)


### Resampling while 'keeping' certain parameters/dimensions


```python
# Construct covariant 2D dataset where the 0th parameter takes on discrete values
xx = np.random.randint(2, 7, 1000)
yy = np.random.normal(4, 2, xx.size) + xx**(3/2)
data = [xx, yy]

# 2D plotting settings: disable the 2D histogram & disable masking of dense scatter-points
dist2d = dict(hist=False, mask_dense=False)

# Draw a corner plot
kale.corner(data, dist2d=dist2d)

nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_29_0.png)


A standard KDE resampling will smooth out the discrete variables, creating a smooth(er) distribution.  Using the `keep` parameter, we can choose to resample from the actual data values of that parameter instead of resampling with 'smoothing' based on the KDE.


```python
kde = kale.KDE(data)

# ---- Resample the data both normally, and 'keep'ing the 0th parameter values ---- #
resamp_stnd = kde.resample()                                                        #
resamp_keep = kde.resample(keep=0)                                                  #
# --------------------------------------------------------------------------------- #

corner = kale.Corner(2)
dist2d['median'] = False    # disable median 'cross-hairs'
h1 = corner.plot(resamp_stnd, dist2d=dist2d)
h2 = corner.plot(resamp_keep, dist2d=dist2d)

corner.legend([h1, h2], ['Standard', "'keep'"])
nbshow()
```


![png](https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_files/demo_31_0.png)


## Development & Contributions

Please visit the `github page <https://github.com/lzkelley/kalepy>`_ for issues or bug reports.  Contributions and feedback are very welcome.

Contributors:
* Luke Zoltan Kelley (@lzkelley)
* Zachary Hafen (@zhafen)

JOSS Paper:
* Kexin Rong (@kexinrong)
* Arfon Smith (@arfon)
* Will Handley (@williamjameshandley)


## Attribution

A JOSS paper has been published on the `kalepy` package.  If you have found this package useful in your research, please add a reference to the code paper:

.. code-block:: tex

    @article{Kelley2021,
      doi = {10.21105/joss.02784},
      url = {https://doi.org/10.21105/joss.02784},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {57},
      pages = {2784},
      author = {Luke Zoltan Kelley},
      title = {kalepy: a Python package for kernel density estimation, sampling and plotting},
      journal = {Journal of Open Source Software}
    }

## Development & Contributions

Please visit the `github page <https://github.com/lzkelley/kalepy>`_ for issues or bug reports.  Contributions and feedback are very welcome.

Contributors:
* Luke Zoltan Kelley (@lzkelley)
* Zachary Hafen (@zhafen)

JOSS Paper:
* Kexin Rong (@kexinrong)
* Arfon Smith (@arfon)
* Will Handley (@williamjameshandley)


## Attribution

A JOSS paper has been submitted.  If you have found this package useful in your research, please add a reference to the code paper:

.. code-block:: tex

    @article{kalepy,
      author = {Luke Zoltan Kelley},
      title = {kalepy: a python package for kernel density estimation and sampling},
      journal = {The Journal of Open Source Software},
      publisher = {The Open Journal},
    }
## To-Do / Known-Issues
- **Optimization needed**.  Things are done in (generally) the simplest ways, currently, need to be optimized for performance (both speed and memory [e.g. with reflections]).  Especially in the case of finite-support kernels, the calculations can be drastically sped up.  Can also use an approximation for infinite-support kernels, truncating at some threshold value of sigma (or percentile; etc).
- Try using `sp.stats.rv_continuous` as base-class for 'Distribution' to provide functionality like 'ppf' etc.
- Differences between covariance-matrix elements of numerous orders of magnitude can cause spurious results, in particular in the PDF marginalized over parameters.  See "KDE::Dynamic Range" docstrings.  Currently this is checked for in the `KDE._finalize()` method, at the end of initialization, and a warning is given if the dynamic range seems too large.
- Move all checking/sanitizing functionality to `KDE` and have `kernels` (etc) assume it's correct.
  - e.g. extrema, points, reflection, params, etc
- Add documenation/examples for base drawing functions in plotting submodule (e.g. `draw_contour2d`, `draw_hist1d`, etc).
- Add tests/documentation for `sample` submodule.

- `kalepy/`
    - Allow for calculating PDF and resampling in only particular dimensions/parameters.
        - FIX: Doesn't work quite right for non-fixed bandwidth, bandwidth needs to be re-calculated for different number of dimensions
    - `tests/`
        - No tests currently check that proper errors are raised.
        - Make sure tests check both cases of `num_points > num_data` and visa-versa (e.g. in PDF calculation).
    - `kernels.py`
        - See if `_resample_clear` and `_resample_reflect` can be combined.
    - `kde.py`
      - `KDE`
        - Explore more efficient ways of calculating the CDF using the underlying kernels instead of integrating the PDF.
        - Use different methods for `grid` edges in ND, instead of broadcasting and flattening (inefficient).
    - `plot.py`
      - Add some way of tracking how many carpet plots have been added to an axis to automatically offset them appropriately in different situations
      - Finish `scatter` method, including 1D carpet (optional)



## Current



## v1.2 - 2021/08/05

- Plotting Improvements
  - Allow an `origin` argument to be specified for corner plots, placing the location of the triangle at one of `bl` (default), `tl`, `tr`, `br`.

- New `kalepy.sample` submodule for sampling from continuously defined functions.
  - Functionality is defined in the `Sample_Grid` class, and the function `sample_grid()` provides a simple API to construct an instance and use it to construct sample points.  Accessible directly from top-level of package.
  - `sample_grid_proportional()` and `sample_outliers()`



## v1.1 - 2021/03/02

- Allow `covariance` to be manually specified in KDE constructor.
- New `KDE.from_hist()` method for constructing KDEs based on existing distributions (instead of finite points).
- Deprecated CDF functionality removed (for the time being).
- `Triweight` distribution works.
- Significant PDF evaluation speed improvements using numba
  - Sampling and evaluation code simplified.
- Use `abc.ABC` base classes
  - `Distribution(object)` ==> `_Distribution(abc.ABC)`
- Plotting improvements
  - BUG: fix incorrect label in rotate bottom-right panels of corner plots
  - Allow `target` lines to be drawn on corner plots using `Corner.target()`
  - Add arguments to limit the number of carpet and scatter points drawn


## v1.0.0 - 2021/01/21

- DOCS: significant expansion of documentation, both docstrings and sphinx (readthedocs.org).
- Allow `kale.plot.Corner` instances to accept externally-created axes.
- Simplify handling of `reflect` arguments.
- Improve bin-edge guessing.


## v0.5 - 2020/12/30

- Complete restructure of `kalepy.plot` submodule, particularly in the API.
- Extensive addition and improvements of documentation, both inline docstrings, and the addition of sphinx docs now available on [kalepy.readthedocs.io](kalepy.readthedocs.io).
  - This includes new demo/test notebooks which are automatically incorporated into the `README.md` and sphinx documentation.
  - Documentation, testing, and examples are now included for core plotting functionality.  More is needed for the base drawing functions (e.g. `draw_contour2d`, `draw_hist1d`, etc)
- `kalepy` paper
  - Fixed typos pointed out by JOSS referees.
  - Added citation and comparison to `GetDist` package.
- BUG: `weights` was not being passed correctly during resampling (without reflection).
- MAINT: fixed a huge number of deprecation warnings now raised by numpy about operations on jagged arrays.
  - Improved functionality of `kale.utils.jshape` and `kale.utils.really1d` functions to accommodate.
- General plotting improvements
  - The handling of colors and colormaps: plotting methods will automatically select the next colorcycle color, construct a matching colormap, and synchronize the color of all plot components.
  - The handling of quantiles for confidence and contour components: are now handles much more self-consistently and with a simpler API.
  - Drawing functions (e.g. `carpet`, `dist1d` and `dist2d`) will load the current, active axes by default.


## v0.4 - 2020/10/12

- Added paper submitted to JOSS


## v0.3.3 - 2020/07/27

- API:
  - Removed `KDE.pdf_grid` method, instead use `KDE.pdf(... grid=True)`.
  - `KDE.pdf(...)` just calls `KDE.density(..., probability=True)`
    - NOTE: this means that, like `density` the `pdf()` function now returns a (2,) tuple of the evaluation points in addition to the density values!
  - The arguments `reflect` and `params` can now be used in tandem.

- `kalepy/`
  - `kde.py`
    - `KDE`
      - `pdf(...)` is now identical to `density(..., probability=True)`
      - `pdf_grid()` [DELETED]
        - Call `pdf(..., grid=True)` instead.
  - `kernels.py`
    - `Kernel`
      - `density()` <== `pdf`, `_pdf_clear`, and `_pdf_reflect`
        - Combined latter functions into single new method.



## v0.3.2 - 2020/06/08

- `reflect` arguments: `True` can now be given (single value, or for a particular parameter/dimension), in which case the KDE will guess the reflection points based on the data extrema (in all dimensions, or only the target ones).  This happens in `kernels._check_reflect`.
- General bug fixes.
- Improve kwarg handling in plotting.

- API
  - `kalepy.density()`
    - BUG: fixed issue in 'grid' mode where output points didn't match values in shape.
    - Add `grid` kwarg.

- `kalepy/`
  - `kernels.py`
    - `_check_reflect()`
      - Added boolean functionality for `reflect` arguments, which are then replaced with data extrema as needed.
  - `plot.py`
    - General bug fixes, improvements in kwarg handling.
    - Return `handles` from plotting functions to allow for legends.
  - `utils.py`
    - New methods for checking / handling jagged arrays (`flatten()`, `flatlen()`, `isjagged()` and `jshape`)

- `notebooks/`
  - `api.ipynb`  [NEW-FILE]
    - New notebook for running API tests.

- `gen_readme.py`  [NEW-FILE]
  - Script to automatically assemble the `README.md` file based on an input template `_README.md` and the jupyter notebook `demo.ipynb`.  Automatically takes care of image files, and updating them with git.
- `README.md`
  - Updated (using `gen_readme.py`) to include new, cleaner examples (primarily using top-level API).


## v0.3.1 - 2020/04/30
- Improved how 'edges' (both for bins and PDF evaluation) are constructed, especially in multiple dimensions.  `KDE` constructs extrema from the given data and then calls `utils.parse_edges`.

- `kalepy/`
  - `__init__.py`
    - `corner()`  [NEW-METHOD]
      - New top-level API method for constructing corner plots using either a dataset or KDE instance.
    - `density()`  [NEW-METHOD]
      - Interface to `KDE.density()`
    - `resample()`  [NEW-METHOD]
      - Interface to `KDE.resample()`
  - `kde.py`  <==  `kde_base.py`  [RENAME]
    - Improved how 'edges' are constructed.  Constructs `extrema` based on input data, and uses `utils.parse_edges` to construct edges.
    - `_guess_edges()`  [REMOVED]
    - `KDE`
      - `density()`  [NEW-METHOD]
        - Calculate density using KDE, where 'density' can either be number density or probability density (i.e. 'pdf').
      - `pdf()`
        - Now calls `density()` using `probability=True`.
  - `kernels.py`
  - `plot.py`
    - Methods for constructing "corner" plots (based strongly on Dan Foreman-Mackey's `corner` package).
    - `Corner`
      - Class for managing corner plots and plotting scatter data or KDE PDFs.

    - `corner_data()`
      - Higher-level function for constructing a full corner plot given scatter-data.
    - `draw_carpet()`  <==  `draw_carpet_fuzz()`  [RENAME]
      - Add `rotate` argument to plot vertically instead of horizontally.
    - `hist()` [NEW-METHOD]
      - Calculate histogram using `utils.histogram()`, then draw it using `_draw_hist1d()`.
    - `utils()`
      - Add `positive` argument to filter by positive definite values.
    - `_get_smap()`  <==  `smap()`  [RENAME]
      - Add `log` argument to specify log-scaling.
  - `utils.py`
    - `histogram()`  [NEW-METHOD]
      - Calculate histograms with both `density` and `probability` parameters (instead of combined like in numpy).
    - `parse_edges()`
      - Allow `weights` to be passed for calculating effective number of data points and inter-quartile ranges
    - `quantiles()`  <==  `percentiles()`
    - `stats()`  [NEW-METHOD]
      - Combines `array_str()` and `stats_str()` output.
    - `_get_edges_1d()`
      - BUG: avoid negative bin-width for very small number of data points.

- `notebooks/`
  - `plotting.ipynb`  [NEW-FILE]
    - For testing and demonstration of plotting methods, especially corner plots.
  - `kde.ipynb`
    - Add corner plots using the `corner.py` submodule.

- `convert_notebook_tests.py`  <==  `build_notebook_tests.py`  [RENAME]



## v0.3.0 - 2020/04/07

- Started working on cleaning up the API (i.e. outward visible functions and structures).
  - New API Functions: `kalepy.pdf()`, `kalepy.cdf()`
- Cleanup variable naming conventions in KDE and Kernels.

- BUG: calculating PDF with `params` given would often result in an error from bad checking of edges/grid shapes.

- `kalepy/`
  - `__init__.py`
    - `pdf()`  [NEW-METHOD]
      - Convenience / API Method for constructing a quick PDF based on the given data.
    - `cdf()`  [NEW-METHOD]
      - Convenience / API Method for constructing a quick CDF based on the given data.
  - `kde_base.py`
    - `KDE`
      - BUG: when providing a scalar value for bandwidth, it was still being multiplied by the data covariance (as is needed for Scott and Silverman rules).  If scalar value(s) are provided do not rescale by covariance.
      - `cdf()`  [NEW-METHOD]
        - Calculate the CDF by integrating the KDE-derived CDF.  This could be done much better.
        - Seems to be working based on simple tests in 1D and 2D.
  - `plot.py` [NEW-FILE]
    - Plotting related functionality; not imported by default - primarily for internal usage.
    - `align_axes_loc()`  [NEW-METHOD]
      - Align a twin axes to a particular location of the base axes.
    - `draw_carpet_fuzz()`  [NEW-METHOD]
      - Draw a fuzz-style carpet plot
    - `nbshow()`  [moved from `utils.py`]
    - `save_fig()`  [moved from `utils.py`]
    - `smap()`  [NEW-METHOD]
      - Construct a ScalarMappable object (with colormap and normalization) for plotting.
    - `Plot_Control`  [moved from `utils.py`]
  - `utils.py`
    - Moved plotted related methods to `plot.py`
    - `assert_true()`  [NEW-METHOD]
      - Internal testing method.
    - `bins()`
      - Added some docstrings
    - `cumsum()` [NEW-METHOD]
      - Calculate cumulative sums along either a single axis, or all axes (unlike `numpy.cumsum`)
    - `cumtrapz()`
      - Added docstrings
    - `really1d()`  [NEW-METHOD]
      - Check if the given array is really one-dimensional (as opposed to a jagged array)
    - `run_if()`
      - Add `otherwise` argument for functions to run when negation
      - Applies to all `run_if_*` methods.
    - `spacing()`
      - BUG: convert `num` to integer before usage.
  - `tests/`
    - `test_utils.py`
      - Added tests for `cumsum()`
- `notebooks/`
  - Update and use the `init.ipy` for the initialization cell of each notebook.  Default save plots/files to notebooks/output


## v0.2.4 - 2020/03/25
- `Triweight` kernel temporarily disabled as it's having normalization problems in ND > 1.

- `kalepy/`
    - `kde_base.py`
        - `class KDE`
            - Addition (uncaught) keyword-arguments are passed from `KDE` initialization to `Kernel` initialization, so that additional arguments (e.g. `chunk`) can be passed along.
    - `kernels.py`
        - BUG: `Triweight` kernel is not working --> disabled kernel.
        - `class Kernel`
            - Implemented 'chunking' for resampling calculation.  Currently only reflection.
                - This produces an *extreme* memory and time performance increase.  For certain parameters, empirically a chunk size of ~ 1e5 seems to work best.
            - `resample()`
                - BUG: non-integer values of `size` would result in an error.
        - `class Distribution`
            - Significant improvements to the way CDFs are handled.
            - `ppf()`  [new-function]
                - "Percent point function" the inverse of the CDF (returns quantiles given cumulative-probabilities).
    - `utils.py`
        - `bound_indices()`
            - BUG: error in boolean logic.
        - `check_path()`  [new-function]
            - Create the given path if it does not already exist.
        - `cumtrapz()`  [new-function]
            - Cumulative summation using the trapezoid-rule.  Light wrapper around  the `trapz_dens_to_mass()` function.
        - `modify_exists()`  [new-function]
            - Modify the given filename if it already exists.
        - `run_if()`  [new-function]
            - New functions for running passed methods if the current environment is the target environment.
        - `save_fig()`  [new-function]
            - Save a `matplotlib` figure adding convenience features.
- `docs/`
    - `logo/`
        - Logo associated data files.
- `notebooks/`
    - `performance.ipynb`  [new-file]
        - New notebook for performance checks, comparisons and diagnostics.



## v0.2.3 - 2019/06/17
- Added code producing a `kalepy` logo, which is added to the attached media and README file.
- Updated notebooks to fix a few minor errors.



## v0.2.2 - 2019/06/11
- Significant improvement in memory and speed while resampling with reflecting boundaries by implementing chunking.



## v0.2.1 - 2019/06/09
- `kalepy/`
    - `__init__.py`
        - Import desired API methods into module namespace.  Use `__all__` in both `kernels.py` and `utils.py`.
    - `kde_base.py`
        - `class KDE`
            - Introduce `helper` argument upon initialization which determines if extra checks and verbose feedback are given.
            - Introcuce `bw_rescale` initialization argument to rescale the bw-matrix by some factor (matrix, or array).
            - `pdf_grid()`  [new-function]
                - Convenience / wrapper function to calculate the PDF given the edges of a grid.
    - `kernels.py`
        - Introduce `helper` parameter, see `class KDE`
        - Allow the `keep` parameter to be `True` in which case all parameters are kept, or `False` and none are kept (same as `None`).
        - `_check_reflect()`
            - Add additional checks for where the reflection boundaries are relative to the data-values and bandwidth.
        - `_resample_reflect()`
            - BUG: reflection was actually a periodic boundary (ish), instead of reflection.  Not sure why it was still behaving well in testing...
            - BUG: reflection was unnecessarily duplicating (already reflected) data, making fewer new points valid.
    - `utils.py`
        - `ave_std()`  [new-function]
            - Calculation of (optionally) *weighted* average and standard-deviation.
        - `bound_indices()`
            - Allow boundaries to be `None` (for no boundaries)
        - `percentiles()`  [new-function]
            - Copied from `zcode.math.statistic`, allows for weighted percentiles.
        - `stats_str()`
            - Copied function from `zcode.math.math_core` with more extended functionality.
        - `trapz_dens_to_mass()`
            - New argument `axis` to integrate only along target axes.
        - `trapz_nd()`
            - New argument `axis` to integrate only along target axes.
- `notebooks/`
    - `init.ipy`     [new-file]
        - Convenience script for setting up the imports in each notebook file
    - `utils.ipynb`  [new-file]
        - New notebook for testing/exploring the `utils.py` submodule.



## v0.2 – 2019/06/03
- Module renamed from `kdes` to `kalepy`.
- Notebooks are now included in travis unit testing.
- Added skeleton for sphinx documentation; not written yet.

- `README.md`
    - Added installation information and basic examples.
- `kalepy/`
    - `bandwidths.py`
    - `kde_base.py`  [new-file]
        - `class KDE`  [new-class]
            - Primary API for using the `kalepy` package.  Uses passed data and options to construct KDEs by interfacing with `Kernel` instances.
            - The `KDE` class calculates the bandwidth and constructs a `kernel` instance, and handles passing the data and covariance matrix to the kernel as needed.
            - `pdf()`
                - Interface to the kernel instance method: `kernel.pdf()`
            - `resample()`
                - Interface to the kernel instance method: `kernel.resample()`
    - `kernels.py`  [new-file]
        - Stores classes and methods for handling the kernels and their underlying distribution functions.
        - NOTE: some of the scaling and normalization does not work properly in multi-dimensions for all kernels.
        - `class Kernel`
            - Stores a covariance-matrix and uses it as needed with a `Distribution` class instance.
        - `class Distribution`
            - Subclassed to implement particular distribution functions to use in a kernel.
            - Agnostic of the data and covariance.  The `Kernel` class handles the covariance matrix and appropriately transforming the data.
        - `class Gaussian(Distribution)`
            - Gaussian/Normal distribution function with infinite support.
        - `class Box_Asym(Distribution)`
            - Boxcar/rectangle/uniform function with finite support.
        - `class Parabola(Distribution)`
            - Epanechnikov kernel-function with finite support.
        - `class Triweight`
            - Cubic kernel, similar to Parabola but with additional smooth derivatives.
            - WARNING: does not currently work in multiple-dimensions (normalization is off).
        - `get_all_distribution_classes()`
            - Method to retrieve a list of all `Distribution` sub-classes.  Mostly used for testing.
        - `get_distribution_class()`
            - Convert from the argument to a `Distribution` subclass as needed.  This argument can convert from a string specification of a distribution function to return the actual class.
    - `utils.py`
        - `class Test_Base`
            - Base-class to use in unittests.
        - `add_cov()`
            - Given a covariance matrix, use a Cholesky decomposition to transform the given data to have that covariance.
        - `allclose()`   [new-function]
            - Convenience function for unittests.
        - `alltrue()`    [new-function]
            - Convenience function for unittests.
        - `array_str()`  [new-function]
            - Format an array (or elements of) for printing.
        - `bins()`       [new-function]
            - Generate bin- edges, centers and widths all together.
        - `bound_indices()`
            - Find the indices of parameter space arrays within given bounds.
        - `cov_from_var_cor()`
            - Construct a covariance matrix given a set of variances of parameters, and the correlations between them.
        - `matrix_invert()`
            - Invert a matrix, following back to SVD if it initially fails.
        - `rem_cov()`
            - Given a covariance matrix, use a Cholesky decomposition to remove that covariance from the given data.
        - `stats_str()`  [new-function]
            - Method for calculating percentiles of given data and returning them as a str.
        - `trapz_dens_to_mass()`
            - Use the ndimensional trapezoid rule to convert from densities on a grid to masses (e.g. PDF to PMF).
    - `tests/`
        - `test_distributions.py`
            - Test the underlying distribution functions.
        - `test_kde.py`
            - Test the top-level KDE class and the accuracy of KDE calculation of PDFs and resampling.
        - `test_kernels.py` [new-file]
            - Tests of the kernels directly.
        - `test_utils.py`
            - Test the utility functions.

- `notebooks/`
    - `kernels.ipynb`  [new-file]
        - Examining / testing the behavior of different kernels specifically.
    - `demo.ipynb`     [new-file]
        - Currently includes the material used in the `README.rst`, should be expanded as a quick demonstration / tutorial of the package.



## v0.1 – 2019/05/19
- `kdes/`
    - `__init__.py`
        - `class KDE`
            - Base class for KDE calculations, modeled roughly on the `scipy.stats.gaussian_kde` class.
            - Allows for multidimensional PDF calculation and resampling of data, in multi-dimensional parameter spaces.
            - Reflecting boundary conditions are available in multiple dimensions, both for PDF calculation and resampling.
    - `utils.py`
        - General utility functions for the package.  Methods extracted from the `zcode` package.
        - `midpoints()`
            - Calculate the midpoints between values in an array, either in log or linear space.
        - `minmax()`
            - Calculate the extrema of a given dataset.  Allows for comparison with previous extrema, setting limits, or 'stretching' the return values by a given amount.
        - `spacing()`
            - Construct a linear or log spacing between the given extrema.
    - `tests/`
        - `test_kde.py`
            - Basic tests for the `KDE` base class and its operations.
        - `test_util.py`
            - Basic tests for the utility methods.

- `notebooks/`
    - `kde.ipynb`
        - Includes basic examples and tests with plots.  Mostly the same tests as in the `kdes/tests/` directory, but with plots.
---
title: 'kalepy: a Python package for kernel density estimation, sampling and plotting'
tags:
  - Python
  - astronomy
  - statistics
  - monte carlo methods
authors:
  - name: Luke Zoltan Kelley
    orcid: 0000-0002-6625-6450
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Center for Interdisciplinary Exploration and Research in Astrophysics (CIERA), Northwestern University, USA
   index: 1
 - name: Physics & Astronomy, Northwestern University, USA
   index: 2
date: 2020 October 11
bibliography: paper.bib

---

# Summary

'Kernel Density Estimation' or 'KDE' [@Rosenblatt-1956; @Parzen-1962] is a type of non-parametric density estimation [@Scott-2015] that improves upon the traditional 'histogram' approach by, for example, i) utilizing the exact location of each data point (instead of 'binning'), ii) being able to produce smooth distributions with continuous and meaningful derivatives, and iii) removing the arbitrary offset of an initial bin edge.  The `kalepy` package presents a Python KDE implementation designed for broad applicability by including numerous features absent in other packages.  `kalepy` provides optional weightings, reflecting boundary conditions, support for an arbitrary number of dimensions, numerous kernel (i.e., window) functions, built-in plotting, and built-in resampling.

# Statement of need

Numerous Python KDE implementations exist, for example in `scipy` (`scipy.stats.gaussian_kde`) [@SciPy-2020], `seaborn` (`seaborn.kdeplot`) [@seaborn-2020], `GetDist` [@GetDist-2019] and `KDEpy` [@kdepy-2018].  The `scipy` and `seaborn` tools are simple and accessible, but lack advanced functionality.  The `KDEpy` package provides excellent performance on large numbers of data points and dimensions, but does not include resampling, boundary conditions, or plotting tools.  The `GetDist` package offers extensive methods for plotting samples and utilizes numerous boundary treatments [@GetDist-2019], but lacks a standalone KDE interface or resampling functionality.  `kalepy` provides convenient access to both plotting and numerical results in the same package, including multiple kernel functions, built-in resampling, boundary conditions, and numerous plotting tools for 1D, 2D, and N-dimensional 'corner' plots.  `kalepy` is entirely class-based, and while focusing on ease of use, provides a highly extensible framework for modification and expansion in a range of possible applications.

While `kalepy` has no features specific to any particular field, it was designed for resampling from weighted astronomical datasets.  Consider a population of binaries derived from cosmological simulations.  If the initial population is costly to produce (e.g., requiring tens of millions of CPU hours), and as long as it accurately samples the parameter space of interest, it may be sufficiently accurate to produce larger populations by 'resampling with variation', e.g., using a KDE approach.  Depending on the details of the population, many of the parameters may be highly correlated and often abut a boundary: for example, the mass-ratio defined as the lower-mass component divided by the more massive component, is often highly correlated with the total mass of the binary, and is bounded to the unit interval i.e., $0 < q \equiv M_2 / M_1 \leq 1$.  Faithfully resampling from the population requires handling this discontinuity, while also preserving accurate covariances which may be distorted when transforming the variable, performing the KDE, and transforming back.

# Methods

Consider a $d$ dimensional parameter space, with $N$ data points given by $x_i = (x_{i1}, x_{i2}, ..., x_{id})$, with $i = \{1, ..., N\}$.  Each data point may have an associated 'weight' that is appropriately normalized, $\sum_i^N w_i = 1$.  The kernel density estimate at a general position $x = (x_1, x_2, ..., x_N)$ can be written as,
$$\hat{f}_H(x) = \sum_{i=1}^N w_i K_H(x - x_i),$$
where the kernel is typically expressed as,
$$K_H(x) = \|H\|^{-1/2} K\left(H^{-1/2} x \right).$$
Here $H$ is the 'bandwidth' (or covariance) matrix.  Choosing the kernel and bandwidth matrix produces most of the nuance and art of KDE.  The most common choice of kernel is likely the Gaussian, i.e.,
$$\hat{f}_H(x) = \sum_{i=1}^N \frac{w_i}{(2\pi)^{-d/2}\|H\|^{1/2}} \exp\{(x_j - x_{ij}) {H^j}_k (x^k - {x_i}^k)\}.$$
In the current implementation, the Gaussian, tri-weight, and box-car kernels are implemented, in addition to the Epanechnikov kernel [@Epanechnikov-1969] which in some cases has been shown to be statistically optimal but has discontinuous derivatives that can produce both numerical and aesthetic problems.
Often the bandwidth is chosen to be diagonal, and different rules-of-thumb are typically used to approximate a bandwidth that minimizes typical measures of error and/or bias.  For example, the so-called 'Silverman factor' [@Silverman-1978] bandwidth,
$$H_{ij} = \delta_{ij} \sigma_i \left[ \frac{4}{(d + 2)n}\right]^{1/(d+4)} \;\;\;\textrm{(summation not implied)},$$
where $\delta_{ij}$ is the Kronecker delta, and $\sigma_i$ is the standard deviation (or its estimate) for the $i$th parameter.  In the current implementation, both the Silverman and Scott factor [@Scott-1979] bandwidth estimators are included.

Reflecting boundary conditions can be used to improve reconstruction accuracy.  For example, with data drawn from a log-normal distribution, a standard KDE will produce 'leakage' outside of the domain.  To enforce the restriction that $f(x < 0) = 0$ (which must be known {\textit{a priori}}), the kernel is redefined such that $K_H(x < 0) = 0$, and re-normalized to preserve unitarity\footnote{Note that some implementations instead truncate and renormalize the resulting $\hat{f}_H$ which which incorrectly redistributes probability from near the boundaries to the whole domain.}.  This example is shown in \autoref{fig:one}, with histograms in the upper panel and KDEs on the bottom.

Resampling from the derived PDF can be done much more efficiently in the KDE framework than by the standard method of CDF inversion.  In particular, we can see that sampling from the PDF is identical to re-sampling with replacement from the weighted data points, while shifting each point based on the PDF of the Kernel at that location.

![Data drawn from a log-normal distribution is used to estimate the underlying PDF using histgrams (upper) and KDEs (lower).  The true distribution is shown in magenta.  In the upper panel, the default bins chosen by `matplotlib` are especially uninsightful (blue), while custom bins misrepresent the distributions position when the initial edge is poorly chosen (red).  The data is also included as a 'carpet' plot.  In the lower panel, a Gaussian KDE with no reflection (blue) is compared to one with a reflection at $x=0$, which better reproduces the true PDF.  Data resampled from the reflecting-KDE PDF is shown as the blue 'carpet' points which closely resemble the input data. \label{fig:one}](fig_one.png){width=75%}

`kalepy` has recently been used in astronomy and astrophysics, particularly in @Siwek+2020, @Kelley-2020, @Andrews-2020.

# Acknowledgements

We acknowledge very helpful consultations on statistical nuances from Christopher Berry and Diego Muñoz.  The support and 'beta-testing' performed by Magda Siwek and Jeff Andrews is also much appreciated.  `kalepy` utilizes tools and functionality from `numpy` [@numpy-2020], `matplotlib` [@matplotlib-2007], `scipy` [@SciPy-2020], `jupyter` notebooks [ipython-2007; @jupyter-2016], and `corner` [@Foreman-Mackey-2016].

# References
Development
===========

Please visit the `github page to make contributions to the package. <https://github.com/lzkelley/kalepy>`_  Particularly if you encounter any difficulties or bugs in the code, please `submit an issue <https://github.com/lzkelley/kalepy/issues>`_, which can also be used to ask questions about usage, or to submit general suggestions and feature requests.  Direct additions, fixes, or other contributions are very welcome which can be done by submitting `pull requests <https://github.com/lzkelley/kalepy/pulls>`_.  If you are considering making a contribution / pull-request, please open an issue first to make sure it won't clash with other changes in development or planned for the future.  Some known issues and indended future-updates are noted in the `change-log <https://github.com/lzkelley/kalepy/blob/master/CHANGES.md>`_ file.  If you are looking for ideas of where to contribute, this would be a good place to start.


Change-Log
----------

Updates and changes to the newest version of `kalepy` will not always be backwards compatible.  The package is consistently versioned, however, to ensure that functionality and compatibility can be maintained for dependencies.  Please consult the `change-log <https://github.com/lzkelley/kalepy/blob/master/CHANGES.md>`_ for summaries of recent changes.


Test Suite
----------

If you are making, or considering making, changes to the `kalepy` source code, the are a large number of built in continuous-integration tests, both in the `kalepy/tests <https://github.com/lzkelley/kalepy/tree/master/kalepy/tests>`_ directory, and in the `kalepy notebooks <https://github.com/lzkelley/kalepy/tree/master/notebooks>`_.  Many of the notebooks are automatically converted into test scripts, and run during continuous integration.  If you are working on a local copy of `kalepy`, you can run the tests using the `tester.sh script <https://github.com/lzkelley/kalepy/tree/master/tester.sh>`_ (i.e. running `$ bash tester.sh`), which will include the notebook tests.


Deploying to pypi (pip)
-----------------------

.. code-block:: bash

  $ python setup.py sdist bdist_wheel
  $ twine check dist/<PACKAGE> 
  $ twine upload dist/<PACKAGE>
kalepy.utils
============

.. automodule:: kalepy.utils
   :members:
   :undoc-members:
   :show-inheritance:

.. code:: ipython3

    import kalepy as kale
    import numpy as np
    import matplotlib.pyplot as plt

Top Level Functions
-------------------

kalepy.corner() and the kalepy.Corner class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the full documentation, see:

* `kalepy.plot.corner <kalepy_plot.html#kalepy.plot.corner>`_
* `kalepy.plot.Corner <kalepy_plot.html#kalepy.plot.Corner>`_
* `kalepy.plot.Corner.plot <kalepy_plot.html#kalepy.plot.Corner.plot>`_

Plot some three-dimensional data called ``data3`` with shape (3, N) with
``N`` data points.

.. code:: ipython3

    kale.corner(data3);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_9_0.png


Extensive modifications are possible with passed arguments, for example:

.. code:: ipython3

    # 1D plot settings: turn on histograms, and modify the confidence-interval quantiles
    dist1d = dict(hist=True, quantiles=[0.5, 0.9])
    # 2D plot settings: turn off the histograms, and turn on scatter
    dist2d = dict(hist=False, scatter=True)
    
    kale.corner(data3, labels=['a', 'b', 'c'], color='purple',
                dist1d=dist1d, dist2d=dist2d);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_11_0.png


The ``kalepy.corner`` method is a wrapper that builds a
``kalepy.Corner`` instance, and then plots the given data. For
additional flexibility, the ``kalepy.Corner`` class can be used
directly. This is particularly useful for plotting multiple
distributions, or using preconfigured plotting styles.

.. code:: ipython3

    # Construct a `Corner` instance for 3 dimensional data, modify the figure size
    corner = kale.Corner(3, figsize=[9, 9])
    
    # Plot two different datasets using the `clean` plotting style
    corner.clean(data3a)
    corner.clean(data3b);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_13_0.png


kalepy.dist1d and kalepy.dist2d
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Corner`` class ultimately calls the functions ``dist1d`` and
``dist2d`` to do the actual plotting of each figure panel. These
functions can also be used directly.

For the full documentation, see:

* `kalepy.plot.dist1d <kalepy_plot.html#kalepy.plot.dist1d>`_
* `kalepy.plot.dist2d <kalepy_plot.html#kalepy.plot.dist2d>`_


.. code:: ipython3

    # Plot a 1D dataset, shape: (N,) for `N` data points
    kale.dist1d(data1);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_17_0.png


.. code:: ipython3

    # Plot a 2D dataset, shape: (2, N) for `N` data points
    kale.dist2d(data2, hist=False);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_18_0.png


These functions can also be called on a ``kalepy.KDE`` instance, which
is particularly useful for utilizing the advanced KDE functionality like
reflection.

.. code:: ipython3

    # Construct a random dataset, and truncate it on the left at 1.0
    import numpy as np
    data = np.random.lognormal(sigma=0.5, size=int(3e3))
    data = data[data >= 1.0]
    
    # Construct a KDE, and include reflection (only on the lower/left side)
    kde_reflect = kale.KDE(data, reflect=[1.0, None])
    # plot, and include confidence intervals
    hr = kale.dist1d(kde_reflect, confidence=True);



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_20_0.png


.. code:: ipython3

    # Load a predefined 2D, 'random' dataset that includes boundaries on both dimensions
    data = kale.utils._random_data_2d_03(num=1e3)
    # Initialize figure
    fig, axes = plt.subplots(figsize=[10, 5], ncols=2)
    
    # Construct a KDE included reflection
    kde = kale.KDE(data, reflect=[[0, None], [None, 1]])
    
    # plot using KDE's included reflection parameters
    kale.dist2d(kde, ax=axes[0]);
    
    # plot data without reflection
    kale.dist2d(data, ax=axes[1], cmap='Reds')
    
    titles = ['reflection', 'no reflection']
    for ax, title in zip(axes, titles):
        ax.set(xlim=[-0.5, 2.5], ylim=[-0.2, 1.2], title=title)



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_plot_files/demo_plot_21_0.png

kalepy.kernels
==============

.. automodule:: kalepy.kernels
   :members:
   :undoc-members:
   :show-inheritance:
kalepy.plot
===========

Contents
--------

- `Corner <kalepy_plot.html#kalepy.plot.Corner>`_

   - `Corner.plot() <kalepy_plot.html#kalepy.plot.Corner.plot>`_
   - `Corner <kalepy_plot.html#kalepy.plot.Corner>`_


- `corner() <kalepy_plot.html#kalepy.corner>`_


`kalepy.plot` submodule
=======================

.. automodule:: kalepy.plot
   :members:
   :undoc-members:
   :show-inheritance:

Basic Usage
===========

.. code:: ipython3

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    import kalepy as kale
    
    from kalepy.plot import nbshow

Generate some random data, and its corresponding distribution function

.. code:: ipython3

    NUM = int(1e4)
    np.random.seed(12345)
    # Combine data from two different PDFs
    _d1 = np.random.normal(4.0, 1.0, NUM)
    _d2 = np.random.lognormal(0, 0.5, size=NUM)
    data = np.concatenate([_d1, _d2])
    
    # Calculate the "true" distribution
    xx = np.linspace(0.0, 7.0, 100)[1:]
    yy = 0.5*np.exp(-(xx - 4.0)**2/2) / np.sqrt(2*np.pi)
    yy += 0.5 * np.exp(-np.log(xx)**2/(2*0.5**2)) / (0.5*xx*np.sqrt(2*np.pi))

Plotting Smooth Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Reconstruct the probability-density based on the given data points.
    points, density = kale.density(data, probability=True)
    
    # Plot the PDF
    plt.plot(points, density, 'k-', lw=2.0, alpha=0.8, label='KDE')
    
    # Plot the "true" PDF
    plt.plot(xx, yy, 'r--', alpha=0.4, lw=3.0, label='truth')
    
    # Plot the standard, histogram density estimate
    plt.hist(data, density=True, histtype='step', lw=2.0, alpha=0.5, label='hist')
    
    plt.legend()
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_8_0.png


resampling: constructing statistically similar values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Draw a new sample of data-points from the KDE PDF

.. code:: ipython3

    # Draw new samples from the KDE reconstructed PDF
    samples = kale.resample(data)
    
    # Plot new samples
    plt.hist(samples, density=True, label='new samples', alpha=0.5, color='0.65', edgecolor='b')
    # Plot the old samples
    plt.hist(data, density=True, histtype='step', lw=2.0, alpha=0.5, color='r', label='input data')
    
    # Plot the KDE reconstructed PDF
    plt.plot(points, density, 'k-', lw=2.0, alpha=0.8, label='KDE')
    
    plt.legend()
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_11_0.png


Multivariate Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Load some random-ish three-dimensional data
    data = kale.utils._random_data_3d_02()
    
    # Construct a KDE
    kde = kale.KDE(data)
    
    # Plot the data and distributions using the builtin `kalepy.corner` plot
    kale.corner(kde)
    
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_13_0.png


.. code:: ipython3

    # Resample the data (default output is the same size as the input data)
    samples = kde.resample()
    
    
    # ---- Plot the input data compared to the resampled data ----
    
    fig, axes = plt.subplots(figsize=[16, 4], ncols=kde.ndim)
    
    for ii, ax in enumerate(axes):
        # Calculate and plot PDF for `ii`th parameter (i.e. data dimension `ii`)
        xx, yy = kde.density(params=ii, probability=True)
        ax.plot(xx, yy, 'k--', label='KDE', lw=2.0, alpha=0.5)
        # Draw histograms of original and newly resampled datasets
        *_, h1 = ax.hist(data[ii], histtype='step', density=True, lw=2.0, label='input')
        *_, h2 = ax.hist(samples[ii], histtype='step', density=True, lw=2.0, label='resample')
        # Add 'kalepy.carpet' plots showing the data points themselves
        kale.carpet(data[ii], ax=ax, color=h1[0].get_facecolor())
        kale.carpet(samples[ii], ax=ax, color=h2[0].get_facecolor(), shift=ax.get_ylim()[0])
    
    axes[0].legend()
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_14_0.png


Fancy Usage
===========

Reflecting Boundaries
~~~~~~~~~~~~~~~~~~~~~

What if the distributions you’re trying to capture have edges in them,
like in a uniform distribution between two bounds? Here, the KDE chooses
‘reflection’ locations based on the extrema of the given data.

.. code:: ipython3

    # Uniform data (edges at -1 and +1)
    NDATA = 1e3
    np.random.seed(54321)
    data = np.random.uniform(-1.0, 1.0, int(NDATA))
    
    # Create a 'carpet' plot of the data
    kale.carpet(data, label='data')
    # Histogram the data
    plt.hist(data, density=True, alpha=0.5, label='hist', color='0.65', edgecolor='k')
    
    # ---- Standard KDE will undershoot just-inside the edges and overshoot outside edges
    points, pdf_basic = kale.density(data, probability=True)
    plt.plot(points, pdf_basic, 'r--', lw=3.0, alpha=0.5, label='KDE')
    
    # ---- Reflecting KDE keeps probability within the given bounds
    # setting `reflect=True` lets the KDE guess the edge locations based on the data extrema
    points, pdf_reflect = kale.density(data, reflect=True, probability=True)
    plt.plot(points, pdf_reflect, 'b-', lw=2.0, alpha=0.75, label='reflecting KDE')
    
    plt.legend()
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_18_0.png


Explicit reflection locations can also be provided (in any number of
dimensions).

.. code:: ipython3

    # Construct random data, add an artificial 'edge'
    np.random.seed(5142)
    edge = 1.0
    data = np.random.lognormal(sigma=0.5, size=int(3e3))
    data = data[data >= edge]
    
    # Histogram the data, use fixed bin-positions
    edges = np.linspace(edge, 4, 20)
    plt.hist(data, bins=edges, density=True, alpha=0.5, label='data', color='0.65', edgecolor='k')
    
    # Standard KDE with over & under estimates
    points, pdf_basic = kale.density(data, probability=True)
    plt.plot(points, pdf_basic, 'r--', lw=4.0, alpha=0.5, label='Basic KDE')
    
    # Reflecting KDE setting the lower-boundary to the known value
    #    There is no upper-boundary when `None` is given.
    points, pdf_basic = kale.density(data, reflect=[edge, None], probability=True)
    plt.plot(points, pdf_basic, 'b-', lw=3.0, alpha=0.5, label='Reflecting KDE')
    
    plt.gca().set_xlim(edge - 0.5, 3)
    plt.legend()
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_20_0.png


Multivariate Reflection
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Load a predefined dataset that has boundaries at:
    #   x: 0.0 on the low-end
    #   y: 1.0 on the high-end
    data = kale.utils._random_data_2d_03()
    
    # Construct a KDE with the given reflection boundaries given explicitly
    kde = kale.KDE(data, reflect=[[0, None], [None, 1]])
    
    # Plot using default settings
    kale.corner(kde)
    
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_22_0.png


Specifying Bandwidths and Kernel Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Load predefined 'random' data
    data = kale.utils._random_data_1d_02(num=100)
    # Choose a uniform x-spacing for drawing PDFs
    xx = np.linspace(-2, 8, 1000)
    
    # ------ Choose the kernel-functions and bandwidths to test -------  #
    kernels = ['parabola', 'gaussian', 'box']                            #
    bandwidths = [None, 0.9, 0.15]     # `None` means let kalepy choose  #
    # -----------------------------------------------------------------  #
    
    ylabels = ['Automatic', 'Course', 'Fine']
    fig, axes = plt.subplots(figsize=[16, 10], ncols=len(kernels), nrows=len(bandwidths), sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0.2, wspace=0.05)
    for (ii, jj), ax in np.ndenumerate(axes):
        
        # ---- Construct KDE using particular kernel-function and bandwidth ---- #
        kern = kernels[jj]                                                       # 
        bw = bandwidths[ii]                                                      #
        kde = kale.KDE(data, kernel=kern, bandwidth=bw)                          #
        # ---------------------------------------------------------------------- #
        
        # If bandwidth was set to `None`, then the KDE will choose the 'optimal' value
        if bw is None:
            bw = kde.bandwidth[0, 0]
            
        ax.set_title('{} (bw={:.3f})'.format(kern, bw))
        if jj == 0:
            ax.set_ylabel(ylabels[ii])
    
        # plot the KDE
        ax.plot(*kde.pdf(points=xx), color='r')
        # plot histogram of the data (same for all panels)
        ax.hist(data, bins='auto', color='b', alpha=0.2, density=True)
        # plot  carpet   of the data (same for all panels)
        kale.carpet(data, ax=ax, color='b')
        
    ax.set(xlim=[-2, 5], ylim=[-0.2, 0.6])
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_24_0.png


Resampling
----------

Using different data ``weights``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Load some random data (and the 'true' PDF, for comparison)
    data, truth = kale.utils._random_data_1d_01()
    
    # ---- Resample the same data, using different weightings ---- #
    resamp_uni = kale.resample(data, size=1000)                       # 
    resamp_sqr  = kale.resample(data, weights=data**2, size=1000)      #
    resamp_inv = kale.resample(data, weights=data**-1, size=1000)     #
    # ------------------------------------------------------------ # 
    
    
    # ---- Plot different distributions ----
    
    # Setup plotting parameters
    kw = dict(density=True, histtype='step', lw=2.0, alpha=0.75, bins='auto')
    
    xx, yy = truth
    samples = [resamp_inv, resamp_uni, resamp_sqr]
    yvals = [yy/xx, yy, yy*xx**2/10]
    labels = [r'$\propto X^{-1}$', r'$\propto 1$', r'$\propto X^2$']
    
    plt.figure(figsize=[10, 5])
    
    for ii, (res, yy, lab) in enumerate(zip(samples, yvals, labels)):
        hh, = plt.plot(xx, yy, ls='--', alpha=0.5, lw=2.0)
        col = hh.get_color()
        kale.carpet(res, color=col, shift=-0.1*ii)
        plt.hist(res, color=col, label=lab, **kw)
    
    plt.gca().set(xlim=[-0.5, 6.5])
    # Add legend
    plt.legend()
    # display the figure if this is a notebook
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_27_0.png


Resampling while ‘keeping’ certain parameters/dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Construct covariant 2D dataset where the 0th parameter takes on discrete values
    xx = np.random.randint(2, 7, 1000)
    yy = np.random.normal(4, 2, xx.size) + xx**(3/2)
    data = [xx, yy]
    
    # 2D plotting settings: disable the 2D histogram & disable masking of dense scatter-points
    dist2d = dict(hist=False, mask_dense=False)
    
    # Draw a corner plot 
    kale.corner(data, dist2d=dist2d)
    
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_29_0.png


A standard KDE resampling will smooth out the discrete variables,
creating a smooth(er) distribution. Using the ``keep`` parameter, we can
choose to resample from the actual data values of that parameter instead
of resampling with ‘smoothing’ based on the KDE.

.. code:: ipython3

    kde = kale.KDE(data)
    
    # ---- Resample the data both normally, and 'keep'ing the 0th parameter values ---- #
    resamp_stnd = kde.resample()                                                        #
    resamp_keep = kde.resample(keep=0)                                                  #
    # --------------------------------------------------------------------------------- #
    
    corner = kale.Corner(2)
    dist2d['median'] = False    # disable median 'cross-hairs'
    h1 = corner.plot(resamp_stnd, dist2d=dist2d)
    h2 = corner.plot(resamp_keep, dist2d=dist2d)
    
    corner.legend([h1, h2], ['Standard', "'keep'"])
    nbshow()



.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/master/docs/media/demo_kde_files/demo_kde_31_0.png

kalepy
======

.. toctree::
   :maxdepth: 4

   kalepy
Introduction
============

*Multidimensional kernel density estimation for distribution functions, resampling, and plotting.*

`kalepy on github <https://github.com/lzkelley/kalepy>`_

|travis| |codecov| |rtd| |joss|

.. |travis| image:: https://travis-ci.org/lzkelley/kalepy.svg?branch=master
.. |codecov| image:: https://codecov.io/gh/lzkelley/kalepy/branch/master/graph/badge.svg
.. |rtd| image:: https://readthedocs.org/projects/kalepy/badge/?version=latest
.. |joss| image:: https://joss.theoj.org/papers/10.21105/joss.02784/status.svg
   :target: https://doi.org/10.21105/joss.02784
   
.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/logo_anim_small.gif

.. contents:: :local:

Installation
------------

.. code-block:: bash

    pip install kalepy


or from source, for development:

.. code-block:: bash

    git clone https://github.com/lzkelley/kalepy.git
    pip install -e kalepy



Quickstart
----------

| Basic examples are shown below.
| `The top-level API is documented here, with many KDE and plotting examples, <api.html>`_
| `The README file on github also includes installation and quickstart examples. <https://github.com/lzkelley/kalepy/blob/master/README.md>`_

One dimensional kernel density estimation:
******************************************

.. code-block:: python

   import kalepy as kale
   import matplotlib.pyplot as plt
   points, density = kale.density(data, points=None)
   plt.plot(points, density, 'k-', lw=2.0, alpha=0.8, label='KDE')

.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/kde1d.png
    :alt: my-picture1


One dimensional resampling:
***************************

.. code-block:: python

    # Draw new samples from the KDE reconstructed PDF
    samples = kale.resample(data)
    plt.hist(samples, density=True, alpha=0.5, label='new samples', color='0.65', edgecolor='b')

.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/resamp1d.png


Multi-dimensional kernel density estimation:
********************************************

.. code-block:: python

    # Construct a KDE instance from data, shaped (N, 3) for `N` data points, and 3 dimensions
    kde = kale.KDE(data)
    # Build a corner plot using the `kalepy` plotting submodule
    corner = kale.corner(kde)

.. image:: https://raw.githubusercontent.com/lzkelley/kalepy/dev/docs/media/kde3dresamp.png



Documentation
-------------

A number of examples are included in `the package notebooks <https://github.com/lzkelley/kalepy/tree/master/notebooks>`_, and the `readme file <https://github.com/lzkelley/kalepy/blob/master/README.md>`_.  Some background information and references are included in `the JOSS paper <https://joss.theoj.org/papers/10.21105/joss.02784>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   introduction <index>
   kalepy API <api>
   Full Package Documentation <kalepy>


Development & Contributions
---------------------------

Please visit the `github page to make contributions to the package. <https://github.com/lzkelley/kalepy>`_  Particularly if you encounter any difficulties or bugs in the code, please `submit an issue <https://github.com/lzkelley/kalepy/issues>`_, which can also be used to ask questions about usage, or to submit general suggestions and feature requests.  Direct additions, fixes, or other contributions are very welcome which can be done by submitting `pull requests <https://github.com/lzkelley/kalepy/pulls>`_.  If you are considering making a contribution / pull-request, please open an issue first to make sure it won't clash with other changes in development or planned for the future.  Some known issues and indended future-updates are noted in the `change-log <https://github.com/lzkelley/kalepy/blob/master/CHANGES.md>`_ file.  If you are looking for ideas of where to contribute, this would be a good place to start.

Updates and changes to the newest version of `kalepy` will not always be backwards compatible.  The package is consistently versioned, however, to ensure that functionality and compatibility can be maintained for dependencies.  Please consult the `change-log <https://github.com/lzkelley/kalepy/blob/master/CHANGES.md>`_ for summaries of recent changes.

Test Suite
^^^^^^^^^^

If you are making, or considering making, changes to the `kalepy` source code, the are a large number of built in continuous-integration tests, both in the `kalepy/tests <https://github.com/lzkelley/kalepy/tree/master/kalepy/tests>`_ directory, and in the `kalepy notebooks <https://github.com/lzkelley/kalepy/tree/master/notebooks>`_.  Many of the notebooks are automatically converted into test scripts, and run during continuous integration.  If you are working on a local copy of `kalepy`, you can run the tests using the `tester.sh script (i.e. '$ bash tester.sh') <https://github.com/lzkelley/kalepy/tree/master/tester.sh>`_, which will include the test notebooks.


Attribution
-----------

A JOSS paper has been published on the `kalepy` package.  If you have found this package useful in your research, please add a reference to the code paper:

.. code-block:: tex

    @article{Kelley2021,
      doi = {10.21105/joss.02784},
      url = {https://doi.org/10.21105/joss.02784},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {57},
      pages = {2784},
      author = {Luke Zoltan Kelley},
      title = {kalepy: a Python package for kernel density estimation, sampling and plotting},
      journal = {Journal of Open Source Software}
    }

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
kalepy full package documentation
=================================

.. contents:: :local:

kalepy.kde module
---------------------

.. automodule:: kalepy.kde
    :members:
    :undoc-members:
    :show-inheritance:

kalepy.kernels module
----------------------

.. automodule:: kalepy.kernels
    :members:
    :undoc-members:
    :show-inheritance:

kalepy.plot module
-------------------------

.. automodule:: kalepy.plot
    :members:
    :undoc-members:
    :show-inheritance:

kalepy.utils module
-----------------------

.. automodule:: kalepy.utils
    :members:
    :undoc-members:
    :show-inheritance:

.. Module contents
.. ---------------
.. 
.. .. automodule:: kalepy
..     :members:
..     :undoc-members:
..     :show-inheritance:
kalepy.kde
==========

.. contents:: :local:

Density estimation: `kalepy.density(...)`:
------------------------------------------

.. autofunction:: kalepy.density
    :noindex:

Data resampling: `kalepy.resample(...)`:
----------------------------------------

.. autofunction:: kalepy.resample
    :noindex:
    
The KDE class: `kalepy.KDE`
---------------------------

.. autoclass:: kalepy.KDE
    :members:
    :noindex:

Full submodule documentation
----------------------------

.. automodule:: kalepy.kde
   :members:
   :undoc-members:
   :show-inheritance:
==========
kalepy API
==========

.. contents:: :local:

Kernel Density Estimation
=========================

The primary API is two functions in the top level package: `kalepy.density` and `kalepy.resample`.  Additionally, `kalepy.pdf` is included which is a shorthand for `kalepy.density(..., probability=True)` --- i.e. a normalized density distribution.

Each of these functions constructs a `KDE` (kalepy.kde.KDE) instance, calls the corresponding member function, and returns the results.  If multiple operations will be done on the same data set, it will be more efficient to construct the `KDE` instance manually and call the methods on that.  i.e.

.. code-block:: python

    kde = kalepy.KDE(data)            # construct `KDE` instance
    points, density = kde.density()   # use `KDE` for density-estimation
    new_samples = kde.resample()      # use same `KDE` for resampling


.. include:: kde_api.rst



Plotting Distributions
======================

For more extended documentation, see the `kalepy.plot submodule documentation. <kalepy_plot.html>`_

.. include:: plot_api.rst


.. Module contents
.. ---------------
.. 
.. .. automodule:: kalepy
..    :members:
..    :undoc-members:
..    :show-inheritance:
