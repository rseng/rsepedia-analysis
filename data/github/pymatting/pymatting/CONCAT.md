### 1.1.5

- Add `relative_discard_threshold` for `ichol` preconditioner.

### 1.1.4

- Switch back to `njit` because ahead-of-time compilation caused too many issues with installation.

### 1.1.3

- Add optimization for `estimate_alpha_cf` which should reduce computation time if most pixels in the trimap are known.
- Allow sloppy trimaps.

### 1.1.2

- Recompile ahead-of-time-compiled modules if they are out of date.
- Add a gradient weighting term to `estimate_foreground_ml`.

### 1.1.1

- Compile on first import instead of during build to simplfy PyPI upload process.

### 1.1.0

- Replace just-in-time compilation with ahead-of-time compilation for faster import times.
* Thomas Germer
* Tobias Uelwer
* Joseph Adams
* Christian Clauss
# PyMatting: A Python Library for Alpha Matting
[![License: MIT](https://img.shields.io/github/license/pymatting/pymatting?color=brightgreen)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/workflow/status/pymatting/pymatting/tests?label=tests)](https://github.com/pymatting/pymatting/actions?query=workflow%3Atests)
[![PyPI](https://img.shields.io/pypi/v/pymatting)](https://pypi.org/project/PyMatting/)
[![JOSS](https://joss.theoj.org/papers/9766cab65bfbf07a70c8a835edd3875a/status.svg)](https://joss.theoj.org/papers/9766cab65bfbf07a70c8a835edd3875a)
[![Gitter](https://img.shields.io/gitter/room/pymatting/pymatting)](https://gitter.im/pymatting/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

We introduce the PyMatting package for Python which implements various methods to solve the alpha matting problem.

- **Website and Documentation:** [https://pymatting.github.io/](https://pymatting.github.io)
- **Benchmarks:**  [https://pymatting.github.io/benchmark.html](https://pymatting.github.io/benchmark.html)

![Lemur](https://github.com/pymatting/pymatting/raw/master/data/lemur/lemur_at_the_beach.png)

Given an input image and a hand-drawn trimap (top row), alpha matting estimates the alpha channel of a foreground object which can then be composed onto a different background (bottom row).

PyMatting provides:
- Alpha matting implementations for:
  - Closed Form Alpha Matting [[1]](#1)
  - Large Kernel Matting [[2]](#2)
  - KNN Matting [[3]](#3)
  - Learning Based Digital Matting [[4]](#4)
  - Random Walk Matting [[5]](#5)
- Foreground estimation implementations for:
  - Closed Form Foreground Estimation [[1]](#1)
  - Fast Multi-Level Foreground Estimation (CPU, CUDA and OpenCL) [[6]](#6)
- Fast multithreaded KNN search
- Preconditioners to accelerate the convergence rate of conjugate gradient descent:
  - The *incomplete thresholded Cholesky decomposition* (*Incomplete* is part of the name. The implementation is quite complete.)
  - The V-Cycle Geometric Multigrid preconditioner
- Readable code leveraging [NumPy](https://numpy.org/), [SciPy](https://www.scipy.org/scipylib/index.html) and [Numba](http://numba.pydata.org/)

## Getting Started

### Requirements

Minimal requiremens
* numpy>=1.16.0
* pillow>=5.2.0
* numba>=0.47.0
* scipy>=1.1.0

Additional requirements for GPU support
* cupy-cuda90>=6.5.0 or similar
* pyopencl>=2019.1.2

Requirements to run the tests
* pytest>=5.3.4

### Installation with PyPI

```bash
pip3 install pymatting
```

### Installation from Source

```bash
git clone https://github.com/pymatting/pymatting
cd pymatting
pip3 install .
```

## Example
```python
from pymatting import cutout

cutout(
    # input image path
    "data/lemur/lemur.png",
    # input trimap path
    "data/lemur/lemur_trimap.png",
    # output cutout path
    "lemur_cutout.png")
```

[More advanced examples](https://pymatting.github.io/examples.html)

## Trimap Construction

All implemented methods rely on trimaps which roughly classify the image into foreground, background and unknown reagions.
Trimaps are expected to be `numpy.ndarrays` of type `np.float64`  having the same shape as the input image with only one color-channel.
Trimap values of 0.0 denote pixels which are 100% background.
Similarly, trimap values of 1.0 denote pixels which are 100% foreground.
All other values indicate unknown pixels which will be estimated by the algorithm.


## Testing

Run the tests from the main directory:
```
 python3 tests/download_images.py
 pip3 install -r requirements_tests.txt
 pytest
```

Currently 89% of the code is covered by tests.

## Upgrade

```bash
pip3 install --upgrade pymatting
python3 -c "import pymatting"
```

## Bug Reports, Questions and Pull-Requests

Please, see [our community guidelines](https://github.com/pymatting/pymatting/blob/master/CONTRIBUTING.md).

## Authors

- **Thomas Germer**
- **Tobias Uelwer**
- **Stefan Conrad**
- **Stefan Harmeling**

See also the list of [contributors](https://github.com/pymatting/pymatting/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Citing

If you found PyMatting to be useful for your work, please consider citing our [paper](https://doi.org/10.21105/joss.02481):

```
@article{Germer2020,
  doi = {10.21105/joss.02481},
  url = {https://doi.org/10.21105/joss.02481},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2481},
  author = {Thomas Germer and Tobias Uelwer and Stefan Conrad and Stefan Harmeling},
  title = {PyMatting: A Python Library for Alpha Matting},
  journal = {Journal of Open Source Software}
}
```

## References

<a id="1">[1]</a> 
Anat Levin, Dani Lischinski, and Yair Weiss. A closed-form solution to natural image matting. IEEE transactions on pattern analysis and machine intelligence, 30(2):228–242, 2007.

<a id="2">[2]</a>
Kaiming He, Jian Sun, and Xiaoou Tang. Fast matting using large kernel matting laplacian matrices. In 2010 IEEE Computer Society Conference on Computer Vision and Pattern Recognition, 2165–2172. IEEE, 2010.

<a id="3">[3]</a>
Qifeng Chen, Dingzeyu Li, and Chi-Keung Tang. Knn matting. IEEE transactions on pattern analysis and machine intelligence, 35(9):2175–2188, 2013.

<a id="4">[4]</a>
Yuanjie Zheng and Chandra Kambhamettu. Learning based digital matting. In 2009 IEEE 12th international conference on computer vision, 889–896. IEEE, 2009.

<a id="5">[5]</a>
Leo Grady, Thomas Schiwietz, Shmuel Aharon, and Rüdiger Westermann. Random walks for interactive alpha-matting. In Proceedings of VIIP, volume 2005, 423–429. 2005.

<a id="6">[6]</a>
Germer, T., Uelwer, T., Conrad, S., & Harmeling, S. (2020). Fast Multi-Level Foreground Estimation. arXiv preprint arXiv:2006.14970.

Lemur image by Mathias Appel from https://www.flickr.com/photos/mathiasappel/25419442300/ licensed under [CC0 1.0 Universal (CC0 1.0) Public Domain License](https://creativecommons.org/publicdomain/zero/1.0/).
# Community Guidelines

## Bug Reports
For bug reports please open a new issue. Include a [minimal reproducible example](https://stackoverflow.com/help/minimal-reproducible-example) and state your operating system, as well as the versions of all relevant Python packages you are using.

**Important:** Before opening a new issue, make sure to check wether a similar issue already exists.

## Questions

If you have any other questions about the usage of the library feel free to open an issue or chat with us on [gitter](https://gitter.im/pymatting/community).

**Important:** Before opening a new issue, make sure to check wether a similar issue already exists.

## Contributing Source Code

Thank you for considering to contribute to the PyMatting project!

To make a contribution, follow those steps:

1. [Fork the repository.](https://guides.github.com/activities/forking/)
2. Make your changes.
3. Add tests if appropriate.
4. Add yourself to [CONTRIBUTORS.md](https://github.com/pymatting/pymatting/blob/master/CONTRIBUTORS.md) if you want to.
5. [Run the black code formatter.](https://pypi.org/project/black/)
6. [Run the tests.](https://github.com/pymatting/pymatting#testing)
7. Create a pull request.
8. State in your pull request
    * the changes made
    * if you agree to release the code under the MIT license

For larger changes, it is probably a good idea to open an issue and discuss them first.

### Licensing

* All contributed code must be released under the MIT license.
* You must be legally allowed to release the contributed code under the MIT license.
---
name: Bug report
about: Create a report to help us improve the PyMatting library.
title: "[BUG \U0001F41B]"
labels: ''
assignees: ''

---

**Bug description**

(Bug description here)

**To Reproduce**

(Include a [minimal reproducible example](https://stackoverflow.com/help/minimal-reproducible-example), input/output images and instructions on how to run it.)

**Expected behavior**

(What should happen instead?)

**Images**

(Add relevant images.)

**Library versions:**

(Run the following commands and paste the result here.)

```bash
python --version --version
python -c "import platform;print(platform.platform())"
python -c "import numpy; numpy.show_config()"
python -c "import scipy;scipy.show_config()"
python -c "import numba;print('Numba version:', numba.__version__)"
python -c "import PIL;print('PIL version:', PIL.__version__)"
python -c "from pymatting.__about__ import __version__;print('PyMatting version:', __version__)"

```
---
name: Question
about: Ask us a general question about the PyMatting library.
title: "[Question❓]"
labels: ''
assignees: ''

---


# Benchmarks

Running the benchmarks will take about a day, not including the time required to install all libraries.

```bash
# install required libraries for benchmark
sudo apt install build-essential unzip cmake libboost-all-dev libopenmpi-dev libmumps-dev petsc-dev libsuitesparse-dev swig
pip3 install psutil scikit-umfpack pyamg pymatting natsort

# download pymatting
git clone https://github.com/pymatting/pymatting
cd pymatting
# download test images
python3 tests/download_images.py

# build libraries
cd benchmarks
sh rebuild.sh

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

python3 calculate_laplacian_error.py
python3 plot_laplacian_error_per_image.py

python3 benchmark_image_sizes.py
python3 plot_results.py
```

If you should be unable to install a specific solver, you can disable it by removing the corresponding line from the list `SOLVER_NAMES` in `pymatting/benchmarks/config.py`.

For faster debugging, it might be helpful to uncomment `SCALES` in the config file.

Sometimes it helps to build packages from source instead of binaries:

```bash
pip3 uninstall <package>
pip3 install <package> --no-binary :all:
```
## Building the docs

1. Install the following packages:
```bash
pip3 install Sphinx sphinxcontrib-bibtex==1.0.0 nbsphinx sphinx_rtd_theme
```
2. Run `./build.sh`
3. The files will appear in `pymatting/build/html`.

## Doc string format

https://numpydoc.readthedocs.io/en/latest/format.html

## Doc string example:

```
    """This function splits the trimap into foreground pixels, background pixels and classifies each pixel as known or unknown. 

    Foreground pixels are pixels where the trimap has value 1.0. Background pixels are pixels where the trimap has value 0.

    Parameters
    ----------
    trimap: array_like
        Trimap with shape :math:`h\\times w`
    flatten: boolean
        If true np.flatten is called on the trimap

    Returns
    -------
    is_fg: np.ndarray
        Boolean array indicating which pixel belongs to the foreground
    is_bg: np.ndarray
        Boolean array indicating which pixel belongs to the background
    is_known: np.ndarray
        Boolean array indicating which pixel is known
    is_unknown: np.ndarray
        Boolean array indicating which pixel is unknown

    """
```
Examples
========

We provide different examples at different levels of abstraction.

.. _example-simple:

Simple Example
---------------

This simple example is intended for application-oriented users.
All parameters were set beforehand and should work well on most images.
The :code:`cutout()` method employs closed-form alpha matting :cite:`levin2007closed` and multi-level foreground extraction :cite:`germer2020multilevel`.

.. code-block:: python

    from pymatting import cutout

    cutout(
       # input image path
       "../data/lemur/lemur.png",
       # input trimap path
       "../data/lemur/lemur_trimap.png",
       # output cutout path
       "lemur_cutout.png")


Advanced Example
----------------

The following example demonstrates the use of the :code:`estimate_alpha_cf()` method as well as the :code:`estimate_foreground_ml()` method.
Both methods can be easily replaced by other methods from the :code:`pymatting.alpha` and from the :code:`pymatting.foreground` module, respectively.
Parameters can be tweaked by passing them to the corresponding function calls.

.. code-block:: python

    from pymatting import *
    import numpy as np

    scale = 1.0

    image = load_image("../data/lemur/lemur.png", "RGB", scale, "box")
    trimap = load_image("../data/lemur/lemur_trimap.png", "GRAY", scale, "nearest")

    # estimate alpha from image and trimap
    alpha = estimate_alpha_cf(image, trimap)

    # make gray background
    background = np.zeros(image.shape)
    background[:, :] = [0.5, 0.5, 0.5]

    # estimate foreground from image and alpha
    foreground = estimate_foreground_ml(image, alpha)

    # blend foreground with background and alpha, less color bleeding
    new_image = blend(foreground, background, alpha)

    # save results in a grid
    images = [image, trimap, alpha, new_image]
    grid = make_grid(images)
    save_image("lemur_grid.png", grid)

    # save cutout
    cutout = stack_images(foreground, alpha)
    save_image("lemur_cutout.png", cutout)

    # just blending the image with alpha results in color bleeding
    color_bleeding = blend(image, background, alpha)
    grid = make_grid([color_bleeding, new_image])
    save_image("lemur_color_bleeding.png", grid)


Expert Example
--------------

The third example provides an insight how PyMatting is working under-the-hood. The matting Laplacian matrix :code:`L` and the system of linear equations :code:`A x = b` are constructed manually. The solution vector :code:`x` is the flattened alpha matte.
The alpha matte :code:`alpha` is then calculated by solving the linear system using the :code:`cg()` method. The convergence of the :code:`cg()` method is accelerated with a preconditioner using the :code:`ichol()` method.
This example is intended for developers and (future) contributors to demonstrate the implementation of the different alpha matting methods.

.. code-block:: python

    from pymatting import *
    import numpy as np
    import scipy.sparse

    scale = 1.0

    image = load_image("../data/lemur/lemur.png", "RGB", scale, "box")
    trimap = load_image("../data/lemur/lemur_trimap.png", "GRAY", scale, "nearest")

    # height and width of trimap
    h, w = trimap.shape[:2]

    # calculate laplacian matrix
    L = cf_laplacian(image)

    # decompose trimap
    is_fg, is_bg, is_known, is_unknown = trimap_split(trimap)

    # constraint weight
    lambda_value = 100.0

    # build constraint pixel selection matrix
    c = lambda_value * is_known
    C = scipy.sparse.diags(c)

    # build constraint value vector
    b = lambda_value * is_fg

    # build linear system
    A = L + C

    # build ichol preconditioner for faster convergence
    A = A.tocsr()
    A.sum_duplicates()
    M = ichol(A)

    # solve linear system with conjugate gradient descent
    x = cg(A, b, M=M)

    # clip and reshape result vector
    alpha = np.clip(x, 0.0, 1.0).reshape(h, w)

    save_image("lemur_alpha.png", alpha)

***************
Getting Started
***************

Requirements
############

* numpy>=1.16.0
* pillow>=5.2.0
* numba>=0.44.0
* scipy>=1.1.0

Additional Requirements (for GPU support)
#########################################

* cupy-cuda90>=6.5.0 or similar
* pyopencl>=2019.1.2

Installation
############
To install PyMatting simply run:

.. code-block::
      
   git clone https://github.com/pymatting/pymatting
   cd pymatting
   pip3 install .

Testing
#######
Run the tests from the main directory:

.. code-block::

   python3 tests/download_images.py
   pip3 install -r requirements_tests.txt
   pytest

Pytest will throw a warning if PyOpenCL or CuPy are not available.
****************************
Benchmarks and Visualization
****************************


Quality
#######

To evaluate the performance of our implementation we calculate the mean squared error on the unknown pixels of the benchmark images of :cite:`rhemann2009perceptually`. 

.. _laplacians_quality_many_bars:
.. figure:: figures/laplacian_quality_many_bars.png
   :align: center
	    
   Figure 1: Mean squared error of the estimated alpha matte to the ground truth alpha matte.

.. _laplacians:
.. figure:: figures/laplacians.png
   :align: center

   Figure 2: Mean squared error across all images from the benchmark dataset.

Visualization
##############

The following videos show the iterates of the different methods. Note that the videos are timewarped.

   .. raw:: html
	    
      <table style="width:100%">
	 <tr align="center">
	 <td>
	 <embed>
	   <video width="320" height="180" loop autoplay muted playsinline>
	   <source src="https://github.com/pymatting/videos/blob/master/cf_web.mp4?raw=true" type="video/mp4">
	   </video>
	 </embed>
         </td>
	 <td>
	 <embed>
	   <video width="320" height="180" loop autoplay muted playsinline>
	   <source src="https://github.com/pymatting/videos/blob/master/knn_web.mp4?raw=true" type="video/mp4">
	 </video>
	 </embed>
	 </td> 
	 </tr>
	 <tr align="center">
	 <td>CF</td>
	 <td>KNN</td> 
	 </tr>
	 <tr align="center">
	 <td><embed>
	   <video width="320" height="180" loop autoplay muted playsinline>
	   <source src="https://github.com/pymatting/videos/blob/master/lkm_web.mp4?raw=true" type="video/mp4">
	   </video>
	 </embed>
	 </td>
	 <td>
	 <embed>
	   <video width="320" height="180" loop autoplay muted playsinline>
	   <source src="https://github.com/pymatting/videos/blob/master/rw_web.mp4?raw=true" type="video/mp4">
	   </video>
	 </embed>
	 </td> 
	 </tr>
	 <tr align="center">
	 <td>LKM</td>
	 <td>RW</td> 
	 </tr>
      </table>

Performance
###########

We compare the computational runtime of our solver with other solvers: pyAMG, UMFPAC, AMGCL, MUMPS, Eigen and SuperLU. Figure 3 shows that our implemented conjugate gradients method in combination with the incomplete Cholesky decomposition preconditioner outperforms the other methods by a large margin. For the iterative solver we used an absolute tolerance of :math:`10^{-7}`, which we scaled with the number of known pixels, i.e. pixels that are either marked as foreground or background in the trimap.


.. _time_image_size:
.. figure:: figures/time_image_size.png
   :align: center
	    
   Figure 3: Comparison of runtime for different image sizes.

.. _average_running_time:
.. figure:: figures/average_running_time.png
   :align: center
	    
   Figure 4: Peak memory for each solver usage in MB.

.. _average_preak_memory_usage:
.. figure:: figures/average_peak_memory_usage.png
   :align: center
	    
   Figure 5: Mean running time of each solver in seconds.
   
.. pymatting documentation master file, created by
   sphinx-quickstart on Tue Dec 17 14:50:36 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyMatting's documentation!
=====================================

The PyMatting package implements various methods for alpha matting and foreground estimation in Python.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
   start
   pymatting
   examples
   benchmark
   references
   

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
References
===========

.. bibliography:: pymatting.bib

		  
The lemur image was taken from https://www.flickr.com/photos/mathiasappel/25419442300/ `(CC0 1.0 Universal (CC0 1.0) Public Domain License) <https://creativecommons.org/publicdomain/zero/1.0/>`_ by Mathias Appel.
Introduction
============

Alpha Matting
-------------

For an image :math:`I` with foreground pixels :math:`F` and background :math:`B` the alpha matting problem aims to determine opacities :math:`\alpha`, such that the equality

.. math::
   I = \alpha F +(1-\alpha)B
   
holds. This problem is inherently ill-posed since for each pixel we have three equations with seven unknown variables. The alpha matte :math:`\alpha` determine how much a pixel contributes to the foreground and how much to the background of an image.

After estimating the alpha matte :math:`\alpha` the foreground pixels and background pixels can be estimated. We refer to this process as foreground estimation.

.. figure:: figures/lemur_at_the_beach.png
   :align: center

   Figure 1: Input image, input trimap, estimated alpha and extracted foreground.

To estimate the alpha matte Pymatting implements the following methods:

* Closed-form matting :cite:`levin2007closed`
* KNN matting :cite:`chen2013knn`
* Large kernel matting :cite:`he2010fast`
* Learning-based matting :cite:`zheng2009learning`
* Random-walk matting :cite:`grady2005random`


Foreground Extraction
---------------------

Simply multiplying the alpha matte with the input image results in halo artifacts. This motivates the developement of foreground extraction methods.

.. figure:: figures/lemur_color_bleeding.png
   :align: center

   Figure 2: Input image naively composed onto a grey background (left) and extracted foreground placed onto the same background (right).

The following foreground estimation methods are implemented in PyMatting:

* Closed-form foreground estimation :cite:`levin2007closed`
* Multilevel approach :cite:`germer2020multilevel`
