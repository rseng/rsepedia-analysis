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
