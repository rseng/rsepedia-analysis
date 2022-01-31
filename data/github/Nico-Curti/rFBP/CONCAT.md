# Authors #

----------
- Nico Curti - eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([nico.curti2@unibo.it](mailto:nico.curti2@unibo.it))
- Daniele Dall'Olio - Department of Physics and Astronomy, University of Bologna ([daniele.dallolio@studio.unibo.it](mailto:daniele.dallolio@studio.unibo.it))
- Enrico Giampieri - eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([enrico.giampieri@unibo.it](mailto:enrico.giampieri@unibo.it))
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.3] - 2020-10-19

This is the official release for the JOSS paper [![JOSS](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0/status.svg)](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0).

### Added

- Add the [pyproject.toml](https://github.com/Nico-Curti/rFBP/blob/master/pyproject.toml) for the build requirements according to PEP-518.
- Add the Pypi badge for the latest releases.
- Add the Doxygen build via CMake.
- Upload the package to Zenodo.

### Fixed

- Fix bibliography and paper document
- Fix minor issues in the documentation

## [1.0.2] - 2020-10-12

### Added

- Doxygen documentation of the C++ APIs
- Sphinx documentation of the Python APIs
- CI for the documentation build

### Fixed

- python setup using pip, fixing dependencies and installation paths

## [1.0.0] - 2020-07-03

### Added

- First C++ version of the BinaryCommitteeMachineFBP code written in Julia
- Add Cython wrap for the Python support
- Add scikit-learn compatibility
- Add scorer support for the score evaluation
- Add old std compatibility in C++

### Changed

- Change the file organization in C++ modules
- Change magnetization formulas for computational performance improvements

### Fixed

- Add a complete Python support for the magnetization variables
- New performances visualization plots
| **Authors**  | **Project** |  **Build Status** | **Code Quality** | **Coverage** |
|:------------:|:-----------:|:-----------------:|:----------------:|:------------:|
| [**N. Curti**](https://github.com/Nico-Curti) <br/> [**D. Dall'Olio**](https://github.com/DanieleDallOlio) <br/> [**E. Giampieri**](https://github.com/EnricoGiampieri)  |  **rFBP** <br/> [![JOSS](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0/status.svg)](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0) <br/> [![PyPI Latest Release](https://img.shields.io/pypi/v/ReplicatedFocusingBeliefPropagation.svg)](https://pypi.org/project/ReplicatedFocusingBeliefPropagation/)  | **Linux/MacOS** : [![Travis](https://travis-ci.com/Nico-Curti/rFBP.svg?token=7QqsqaQiuDHSyGDT3xek&branch=master)](https://travis-ci.com/Nico-Curti/rFBP) <br/> **Windows** : [![appveyor](https://ci.appveyor.com/api/projects/status/obuq56lhyd90pmup?svg=true)](https://ci.appveyor.com/project/Nico-Curti/rfbp) | **Codacy** : [![Codacy](https://api.codacy.com/project/badge/Grade/a6fdac990b6f4141a5bd9e8171ddaf53)](https://www.codacy.com/manual/Nico-Curti/rFBP?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Nico-Curti/rFBP&amp;utm_campaign=Badge_Grade) <br/> **Codebeat** : [![Codebeat](https://codebeat.co/badges/cc761a7c-79fa-4a66-984f-bef6fd145d34)](https://codebeat.co/projects/github-com-nico-curti-rfbp-master) | [![codecov](https://codecov.io/gh/Nico-Curti/rFBP/branch/master/graph/badge.svg)](https://codecov.io/gh/Nico-Curti/rFBP) |

[![rFBP C++ CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20C++%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+C%2B%2B+CI%22)
[![rFBP Python CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20Python%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+Python+CI%22)
[![rFBP Docs CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20Docs%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+Docs+CI%22)

[![docs](https://readthedocs.org/projects/rfbp/badge/?version=latest)](https://rfbp.readthedocs.io/en/latest/?badge=latest)
[![GitHub pull-requests](https://img.shields.io/github/issues-pr/Nico-Curti/rFBP.svg?style=plastic)](https://github.com/Nico-Curti/rFBP/pulls)
[![GitHub issues](https://img.shields.io/github/issues/Nico-Curti/rFBP.svg?style=plastic)](https://github.com/Nico-Curti/rFBP/issues)

[![GitHub stars](https://img.shields.io/github/stars/Nico-Curti/rFBP.svg?label=Stars&style=social)](https://github.com/Nico-Curti/rFBP/stargazers)
[![GitHub watchers](https://img.shields.io/github/watchers/Nico-Curti/rFBP.svg?label=Watch&style=social)](https://github.com/Nico-Curti/rFBP/watchers)

<a href="https://github.com/UniboDIFABiophysics">
  <div class="image">
    <img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90">
  </div>
</a>

# Replicated Focusing Belief Propagation algorithm

We propose a `C++` version of the [**Replicated Focusing Belief Propagation**](https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl) Julia package.
Our implementation optimizes and extends the original library including multi-threading support and an easy-to-use interface to the main algorithm.
To further improve the usage of our code, we propose also a `Python` wrap of the library with a full compatibility with the [`scikit-learn`](https://github.com/scikit-learn/scikit-learn) and [`scikit-optimize`](https://github.com/scikit-optimize/scikit-optimize) packages.

* [Overview](#overview)
* [Theory](#theory)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Efficiency](#efficiency)
* [Usage](#usage)
* [Testing](#testing)
* [Table of contents](#table-of-contents)
* [Contribution](#contribution)
* [References](#references)
* [Authors](#authors)
* [License](#license)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)

## Overview

The learning problem could be faced through statistical mechanic models joined with the so-called Large Deviation Theory.
In general, the learning problem can be split into two sub-parts: the classification problem and the generalization one.
The first aims to completely store a pattern sample, i.e a prior known ensemble of input-output associations (*perfect learning*).
The second one corresponds to compute a discriminant function based on a set of features of the input which guarantees a unique association of a pattern.

From a statistical point-of-view many Neural Network models have been proposed and the most promising seems to be spin-glass models based.
Starting from a balanced distribution of the system, generally based on Boltzmann distribution, and under proper conditions, we can prove that the classification problem becomes a NP-complete computational problem.
A wide range of heuristic solutions to that type of problems were proposed.

In this project we show one of these algorithms developed by Zecchina et al. [[BaldassiE7655](https://www.pnas.org/content/113/48/E7655)] and called *Replicated Focusing Belief Propagation* (`rFBP`).
The `rFBP` algorithm is a learning algorithm developed to justify the learning process of a binary neural network framework.
The model is based on a spin-glass distribution of neurons put on a fully connected neural network architecture.
In this way each neuron is identified by a spin and so only binary weights (-1 and 1) can be assumed by each entry.
The learning rule which controls the weight updates is given by the Belief Propagation method.

A first implementation of the algorithm was proposed in the original paper [[BaldassiE7655](https://www.pnas.org/content/113/48/E7655)] jointly with an open-source Github repository.
The original version of the code was written in `Julia` language and despite it is a quite efficient implementation the `Julia` programming language stays on difficult and far from many users.
To broaden the scope and use of the method, a `C++` implementation was developed with a joint `Cython` wrap for `Python` users.
The `C++` language guarantees better computational performances against the `Julia` implementation and the `Python` version enhances its usability.
This implementation is optimized for parallel computing and is endowed with a custom `C++` library called [`Scorer`](https://github.com/Nico-Curti/scorer)), which is able to compute a large number of statistical measurements based on a hierarchical graph scheme.
With this optimized implementation and its [`scikit-learn`](https://github.com/scikit-learn/scikit-learn) compatibility we try to encourage researchers to approach these alternative algorithms and to use them more frequently on real context.

As the `Julia` implementation also the `C++` one provides the entire `rFBP` framework in a single library callable via a command line interface.
The library widely uses template syntaxes to perform dynamic specialization of the methods between two magnetization versions of the algorithm.
The main object categories needed by the algorithm are wrapped in handy `C++` objects easy to use also from the `Python` interface.

## Theory

The `rFBP` algorithm derives from an out-of-equilibrium (non-Boltzmann) model of the learning process of binary neural networks [[DallAsta101103](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.77.031118)].
This model mimics a spin glass system whose realizations are equally likely to occur when sharing the same so-called entropy (not the same energy, i.e. out-of-equilibrium).
This entropy basically counts the number of solutions (zero-energy realizations) around a realization below a fixed-distance.

Within this out-of-equilibrium framework, the objective is to maximize the entropy instead of minimizing the energy.
From a machine learning standpoint, we aim at those weights sets that perfectly solve the learning process (zero-errors) and that are mathematically closed to each other.
To this end, the Belief Propagation method [[MézardMontanari](https://web.stanford.edu/~montanar/RESEARCH/book.html)] can be adopted as the underlying learning rule, although it must be properly adjusted to take into account the out-of-equilibrium nature of the model.

See [here](https://github.com/Nico-Curti/rFBP/blob/master/docs/model.md) for further details about the model.

## Prerequisites

C++ supported compilers:

![gcc version](https://img.shields.io/badge/gcc-4.8.5%20|%204.9.*%20|%205.*%20|%206.*%20|%207.*%20|%208.*%20|%209.*-yellow.svg)

![clang version](https://img.shields.io/badge/clang-3.6%20|3.9%20|5.*%20|%206.*%20|%207.*%20|-red.svg)

![msvc version](https://img.shields.io/badge/msvc-vs2017%20x86%20|%20vs2017%20x64|%20vs2019%20x86%20|%20vs2019%20x64-blue.svg)

The `rFBP` project is written in `C++` using a large amount of c++17 features.
To enlarge the usability of our package we provide also a retro-compatibility of all the c++17 modules reaching an usability (tested) of our code from gcc 4.8.5+.
The package installation can be performed via [`CMake`](https://github.com/Nico-Curti/rFBP/blob/master/CMakeLists.txt) or [`Makefile`](https://github.com/Nico-Curti/rFBP/blob/master/Makefile).

If you are using the `CMake` (recommended) installer the maximum version of C++ standard is automatic detected.
The `CMake` installer provides also the export of the library: after the installation you can use this library into other `CMake` projects using a simple `find_package` function.
The exported `CMake` library (`rFBP::rfbp`) is installed in the `share/rFBP` directory of the current project and the relative header files are available in the `rFBP_INCLUDE_DIR` variable.

The `CMake` installer provides also a `rFBP.pc`, useful if you want link to the `rFBP` using `pkg-config`.

You can also use the `rFBP` package in `Python` using the `Cython` wrap provided inside this project.
The only requirements are the following:

* numpy >= 1.15
* cython >= 0.29
* scipy >= 1.2.1
* scikit-learn >= 0.20.3
* requests >= 2.22.0

The `Cython` version can be built and installed via `CMake` enabling the `-DPYWRAP` variable.
The `Python` wrap guarantees also a good integration with the other common Machine Learning tools provided by `scikit-learn` `Python` package; in this way you can use the `rFBP` algorithm as an equivalent alternative also in other pipelines.
Like other Machine Learning algorithm also the `rFBP` one depends on many parameters, i.e its hyper-parameters, which has to be tuned according to the given problem.
The `Python` wrap of the library was written according to `scikit-optimize` `Python` package to allow an easy hyper-parameters optimization using the already implemented classical methods.

## Installation

Follow the instruction about your needs.

A complete list of instructions "for beginners" is also provided for both [cpp](https://github.com/Nico-Curti/rFBP/blob/master/docs/cpp_install.md) and [python](https://github.com/Nico-Curti/rFBP/blob/master/docs/python_install.md) versions.

### CMake C++ installation

We recommend to use `CMake` for the installation since it is the most automated way to reach your needs.
First of all make sure you have a sufficient version of `CMake` installed (3.9 minimum version required).
If you are working on a machine without root privileges and you need to upgrade your `CMake` version a valid solution to overcome your problems is provided [here](https://github.com/Nico-Curti/Shut).

With a valid `CMake` version installed first of all clone the project as:

```bash
git clone https://github.com/Nico-Curti/rFBP
cd rFBP
```

The you can build the `rFBP` package with

```bash
mkdir -p build
cd build && cmake .. && cmake --build . --target install
```

or more easily

```bash
./build.sh
```

if you are working on a Windows machine the right script to call is the [`build.ps1`](https://Nico-Curti/rFBP/blob/master/build.ps1).

**NOTE 1:** if you want enable the OpenMP support (*4.5 version is required*) compile the library with `-DOMP=ON`.

**NOTE 2:** if you want enable the Scorer support compile the library with `-DSCORER=ON`. If you want use a particular installation of the Scorer library or you have manually installed the library following the `README` instructions, we suggest to add the `-DScorer_DIR=/path/to/scorer/shared/scorer` in the command line.

**NOTE 3:** if you want enable the Cython support compile the library with `-DPYWRAP=ON`. The Cython packages will be compiled and correctly positioned in the `rFBP` Python package **BUT** you need to run also the setup before use it.

**NOTE 4:** if you use MagT configuration, please download the `atanherf coefficients` file before running any executable. You can find a downloader script inside the [scripts](https://github.com/Nico-Curti/rFBP/tree/master/scripts) folder. Enter in that folder and just run `python dowload_atanherf.py`.

### Make C++ installation

The `Make` installation requires more attention!
First of all the `Make` installation assumes that you compiler is able to support the c++17 standard: if it is not your case you have to change the `STD` variable into the `Makefile` script.

Then if you call just:

```bash
make
```

you can view the complete list of available examples.
With

```bash
make main
```

you can compile the main example and the `C++` library.

### Python installation

Python version supported : ![Python version](https://img.shields.io/badge/python-3.5|3.6|3.7|3.8-blue.svg)

The easiest way to install the package is to use `pip`

```bash
python -m pip install ReplicatedFocusingBeliefPropagation
```

> :warning: The setup file requires the `Cython` and `Numpy` packages, thus make sure to pre-install them!
> We are working on some workarounds to solve this issue.

The `Python` installation can be performed with or without the `C++` installation.
The `Python` installation is always executed using [`setup.py`](https://github.com/Nico-Curti/blob/master/setup.py) script.

If you have already built the `rFBP` `C++` library the installation is performed faster and the `Cython` wrap was already built using the `-DPYWRAP` definition.
Otherwise the full list of dependencies is build.

In both cases the installation steps are

```bash
python -m pip install -r ./requirements.txt
```

to install the prerequisites and then

```bash
python setup.py install
```

or for installing in development mode:

```bash
python setup.py develop --user
```

> :warning: The current installation via pip has no requirements about the version of `setuptools` package.
> If the already installed version of `setuptools` is `>= 50.*` you can find some troubles during the installation of our package (ref. [issue](https://github.com/Nico-Curti/rFBP/issues/5)).
> We suggest to temporary downgrade the `setuptools` version to `49.3.0` to workaround this `setuptools` issue.

## Efficiency

![Comparison of time performances between the original `Julia` implementation and our `Cython` one of the `rFBP` algorithm varying the input dimension sizes (number of samples, `M`, and number of features, `N`). For each input configuration 100 runs of both algorithm were performed and the results were normalized by the `Julia` implementation. In these cases we fixed the magnetization to **MagP64**.](./img/rfbp_magp_timing.svg)

![Comparison of time performances between the original `Julia` implementation and our `Cython` one of the `rFBP` algorithm varying the input dimension sizes (number of samples, `M`, and number of features, `N`). For each input configuration 100 runs of both algorithm were performed and the results were normalized by the `Julia` implementation. In these cases we fixed the magnetization to **MagT64**.](./img/rfbp_magt_timing.svg)

We test the computational efficiency of our implementation against the original `Julia` one ([update Jul 2020](https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl/tree/179443860083fc68c87e8e53a588bc22187c3ade)).
The tests were performed comparing our `Cython` version of the code (and thus with a slight overhead given by the `Python` interpreter) and the `Julia` implementation.
Varying the dimension sizes (number of samples, `M`, and number of features, `N`) we performed 100 runs of both the algorithms.
We divided our simulation according to the two possible types of magnetizations: `MagP64` and `MagT64`.
As described in the [original implementation](https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl), the `MagP64` type allows fast executions with inexact outcomes by neglecting all `tanh` operations.
In contrast, the `MagT64` exactly follows all theoretical equations with no further approximation, which necessarily causes slower executions.
The obtained results are showed in Fig. [[1](./img/rgbp_magp_timing.svg), [2](./img/rgbp_magt_timing.svg)], respectively.

As can be seen by the two simulations our implementation scales very well with the number of samples and it is quite stable in relation to the number of features.
However, we can not guarantee a perfect parallel execution of our version: also with multi-threading support the scalability of our implementation does not follow a linear trend with the number of available cores.
In our simulation, in fact, we used 32 cores against the single thread execution of the `Julia` implementation but we gained only a 4x and 2x of speedup for `MagT64` and `MagP64`, respectively.
The network training is a sequential process by definition and thus it is hard to obtain a relevant speedup using a parallel implementation.
In this case it is probably jointed to a not perfect parallelization strategy chosen which bring to a not efficient scalability of our algorithm version.
However, the improvements performed to the code allow us to use this algorithm with bigger dataset sizes.

## Usage

You can use the `rFBP` library into pure-Python modules or inside your C++ application.

### C++ Version

The easiest usage of `rFBP` library is given by the two examples provided in the [example](https://github.com/Nico-Curti/rFBP/blob/master/example) folder.
These two scripts include an easy-to-use command line support for both training and test procedure.

To train the model you can just use

```
./bin/train_main
Usage: ./train_main [-threads <std::remove_reference<int>> ] -f <std :: string> [-output <std :: string> ] [-bin <std::remove_reference<bool>> ] [-delimiter <std :: string> ] [-hidden <std::remove_reference<int>> ] [-iteration <std::remove_reference<int>> ] [-seed <std::remove_reference<int>> ] [-randfact <std::remove_reference<double>> ] [-damping <std::remove_reference<double>> ] [-accuracy <std :: string> ] [-protocol <std :: string> ] [-epsilon <std::remove_reference<double>> ] [-steps <std::remove_reference<int>> ] [-mag <std::remove_reference<int>> ] [-inmess <std :: string> ] [-outmess <std :: string> ] [-delmess <std :: string> ] [-binmess <std::remove_reference<bool>> ]

Training BeliefPropagation ${VERSION}

optional arguments:
        -t,   --threads                 Max number of threads exploitable
        -f,   --file                    Pattern Filename (with extension)
        -o,   --output                  Output Filename (with extension)
        -b,   --bin                     File format: (0) Textfile(default), (1) Binary
        -dl,  --delimiter               Delimiter for text files(default: "\t")
        -k,   --hidden                  Number of Hidden Layers(default:3)
        -i,   --iteration               Max Number of Iterations(default: 1000)
        -r,   --seed                    Seed random generator(default: 135)
        -g,   --randfact                Seed random generator of Cavity Messages(default: 0.1)
        -d,   --damping                 Damping parameter(default: 0.5)
        -a,   --accuracy                Accuracy of the messages computation at the hidden units level (choose between 'exact'(default), 'accurate', 'approx', 'none')
        -p,   --protocol                Specify protocol : scooping, pseudo_reinforcement (default), free_scoping, standard_reinforcement
        -e,   --epsilon                 Threshold for convergence(default: 0.1)
        -s,   --steps                   Max Number of Steps for chosen protocol(default: 101)
        -m,   --mag                     Specify Magnetization: (0) MagnetizationP (MagP64), (1) MagnetizationT (MagT64)
        -im,  --inmess                  Input Messages file
        -om,  --outmess                 Output Messages file
        -dm,  --delmess                 Delimiter for Messages files(default: "\t")
        -bm,  --binmess                 Messages files format: (0) Textfile(default), (1) Binary
```

and after training you can test your model using

```
./bin/test_main
Usage: ./test_main [-threads <std::remove_reference<int>> ] -f <std :: string> [-bin <std::remove_reference<bool>> ] -w <std :: string> [-delimiter <std :: string> ] [-output <std :: string> ]

Test BeliefPropagation ${VERSION}

optional arguments:
        -t,   --threads                 Max number of threads exploitable
        -f,   --file                    Pattern Filename (with extension)
        -b,   --bin                     File format: (0) Textfile(default), (1) Binary
        -w,   --weights                 Weights Matrix Filename (with extension)
        -dl,  --delimiter               Delimiter for text files(default: "\t")
        -o,   --output                  Output Filename (no extension)
```

If you are interested in using `rFBP` inside your code you can simply import the [`rfbp.hpp`](https://github.com/Nico-Curti/rFBP/blob/master/hpp/rfbp.hpp) and create a `ReplicatedFocusingBeliefPropagation` object.

Then all the work is performed by the `focusingBP` (template) function.
You can use it with `MagP64` type or `MagT64`.
We recommend the former when quick results are needed and the latter when the weights accuracy is top priority.

The input pattern must be wrapped into a `Pattern` object provided by the library.

```c++
#include <rfbp.hpp>

int main ()
{
  FocusingProtocol fp("pseudo_reinforcement", 101);
  Patterns patterns("patternsfile.csv", false, ",");

  long int ** bin_weights = focusingBP < MagP64 >(3,          // K,
                                                  patterns,   // patterns,
                                                  1000,       // max_iters,
                                                  101,        // max_steps,
                                                  42,         // seed,
                                                  0.5,        // damping,
                                                  "accurate", // accuracy1,
                                                  "exact",    // accuracy2,
                                                  0.1,        // randfact,
                                                  fp,         // fp,
                                                  0.1,        // epsil,
                                                  1,          // nth,
                                                  "",         // outfile,
                                                  "",         // outmess,
                                                  "",         // inmess,
                                                  false       // binmess
                                                  );

  return 0;
}
```

Then you can use the `nonbayes_test` function to predict your test set.

### Python Version

The `rfbp` object is totally equivalent to a `scikit-learn` classifier and thus it provides the member functions `fit` (to train your model) and `predict` (to test a trained model on new samples).

First of all you need to import the `rFBP` modules.

```python
from ReplicatedFocusingBeliefPropagation import MagT64
from ReplicatedFocusingBeliefPropagation import Pattern
from ReplicatedFocusingBeliefPropagation import ReplicatedFocusingBeliefPropagation as rFBP
```

If you want to run your script with multiple cores you can simply import also

```python
from ReplicatedFocusingBeliefPropagation import NTH
```

which is set to the maximum number of core in your computer.

You can start to try the package functionality using a random pattern

```python
N, M = (20, 101) # M must be odd
data = np.random.choice([-1, 1], p=[.5, .5], size=(N, M))
label = np.random.choice([-1, 1], p=[.5, .5], size=(N, ))
```

The input data must be composed by binary variables codified as `[-1, 1]`, since the model works only with spin-like variables.
The next step is the creation of the `Replicated Focusing Belief Propagation` model.

```python
rfbp = rFBP(mag=MagT64,
            hidden=3,
            max_iter=1000,
            seed=135,
            damping=0.5,
            accuracy=('accurate','exact'),
            randfact=0.1,
            epsil=0.5,
            protocol='pseudo_reinforcement',
            size=101,
            nth=NTH)
```

Now you can fit your model and predict:

```python
rfbp.fit(data, label)
predicted_labels = rfbp.predict(data)
```

which is clearly an overfitting! But it works as example :blush:

The internal implementation of the algorithm works with a custom data type called `Pattern` (ref. [here](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/rfbp/Patterns.py)).
You can explicitly use a `Pattern` object or convert your data to it

```python
n_sample, n_feature = (20, 101) # n_feature must be odd
data = np.random.choice(a=(-1, 1), p=(.5, .5), size=(n_sample, n_feature))
labels = np.random.choice(a=(-1, 1), p=(.5, .5), size=(n_sample, ))

pt = Pattern(X=data, y=labels)
# dimensions
assert pt.shape == (n_sample, n_feature)
# data
np.testing.assert_allclose(pt.data, data)
# labels
np.testing.assert_allclose(pt.labels, labels)
```

We suggest the usage of this data type if you have to load your data from file.
We check the consistency of the input variables into the `C++` code (**only** in DEBUG mode) and into the `Python` wrap.

In the [example](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/example/) folder you can find a training/test example using a pattern imported from file (a more realistic example).
Both the `fit` and `predict` functions work using either a `numpy` array and a `Pattern` object.

## Testing

`rFBP` uses CMake to build a full list of tests.
You can disable tests setting the `-DBUILD_TEST=OFF` during the building.
All the test are performed using the [`Catch2`](https://github.com/catchorg/Catch2/) (v2.11.0) library.

The test scripts can be found [here](https://github.com/Nico-Curti/rFBP/blob/master/test).

The Python version of the package is also tested using [`pytest`](https://docs.pytest.org/en/latest/).
To install the package in development mode you need to add also this requirement:

* pytest == 3.0.7

The full list of python test scripts can be found [here](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/rfbp/test).

## Table of contents

Description of the folders related to the `C++` version.

| **Directory**  |  **Description** |
|:--------------:|:-----------------|
| [example](https://github.com/Nico-Curti/rFBP/blob/master/example) | List of example usages for the C++ version of the code. In [train_main.cpp](https://github.com/Nico-Curti/rFBP/blob/master/example/train_main.cpp) we show how to build and train a C++ model and in [test_main.cpp](https://github.com/Nico-Curti/rFBP/blob/master/example/test_main.cpp) how to use this model to perform a prediction. |
| [hpp](https://github.com/Nico-Curti/rFBP/blob/master/hpp)         | Implementation of the C++ template functions and objects used in the `rFBP` library |
| [include](https://github.com/Nico-Curti/rFBP/blob/master/include) | Definition of the C++ function and objects used in the `rFBP` library |
| [src](https://github.com/Nico-Curti/rFBP/blob/master/src)         | Implementation of the C++ functions and objects used in the `rFBP` library |
| [test](https://github.com/Nico-Curti/rFBP/blob/master/test)       | Repository of tests for the C++ codes |

Description of the folders related to the `Python` version (base directory `ReplicatedFocusingBeliefPropagation`).

| **Directory**  |  **Description** |
|:--------------:|:-----------------|
| [example](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/example) | `Python` version of the `C++` examples. In [overall_example.py](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/example/overall_example.py) a full example (train + test) is showed using random patten. |
| [lib](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/lib)         | List of `Cython` definition files |
| [source](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/source)   | List of `Cython` implementation objects |
| [rfbp](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/rfbp)       | List of `Python` wraps |
| [rfbp/test](https://github.com/Nico-Curti/rFBP/blob/master/ReplicatedFocusingBeliefPropagation/rfbp/test) | List of test scripts for the `Python` wraps |

## Contribution

Any contribution is more than welcome :heart:. Just fill an [issue](https://github.com/Nico-Curti/rFBP/blob/master/.github/ISSUE_TEMPLATE/ISSUE_TEMPLATE.md) or a [pull request](https://github.com/Nico-Curti/rFBP/blob/master/.github/PULL_REQUEST_TEMPLATE/PULL_REQUEST_TEMPLATE.md) and we will check ASAP!

See [here](https://github.com/Nico-Curti/rFBP/blob/master/.github/CONTRIBUTING.md) for further informations about how to contribute with this project.

## References

<blockquote>1- D. Dall'Olio, N. Curti, G. Castellani, A. Bazzani, D. Remondini. "Classification of Genome Wide Association data by Belief Propagation Neural network", CCS Italy, 2019. </blockquote>

<blockquote>2- C. Baldassi, C. Borgs, J. T. Chayes, A. Ingrosso, C. Lucibello, L. Saglietti, and R. Zecchina. "Unreasonable effectiveness of learning neural networks: From accessible states and robust ensembles to basic algorithmic schemes", Proceedings of the National Academy of Sciences, 113(48):E7655-E7662, 2016. </blockquote>

<blockquote>3- C. Baldassi, A. Braunstein, N. Brunel, R. Zecchina. "Efficient supervised learning in networks with binary synapses", Proceedings of the National Academy of Sciences, 104(26):11079-11084, 2007. </blockquote>

<blockquote>4- A., Braunstein, R. Zecchina. "Learning by message passing in networks of discrete synapses". Physical Review Letters 96(3), 2006. </blockquote>

<blockquote>5- C. Baldassi, F. Gerace, C. Lucibello, L. Saglietti, R. Zecchina. "Learning may need only a few bits of synaptic precision", Physical Review E, 93, 2016 </blockquote>

<blockquote>6- A. Blum, R. L. Rivest. "Training a 3-node neural network is NP-complete", Neural Networks, 1992 </blockquote>

<blockquote>7- W. Krauth, M. Mezard. "Storage capacity of memory networks with binary coupling", Journal of Physics (France), 1989 </blockquote>

<blockquote>8- H. Huang, Y. Kabashima. "Origin of the computational hardness for learning with binary synapses", Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 2014 </blockquote>

<blockquote>9- C. Baldassi, A. Ingrosso, C. Lucibello, L. Saglietti, R. Zecchina. "Local entropy as a measure for sampling solutions in constraint satisfaction problems", Journal of Statistical Mechanics: Theory and Experiment, 2016 </blockquote>

<blockquote>10- R. Monasson, R. Zecchina. "Learning and Generalization Theories of Large Committee Machines", Modern Physics Letters B, 1995 </blockquote>

<blockquote>11- R. Monasson, R. Zecchina. "Weight space structure and internal representations: A direct approach to learning and generalization in multilayer neural networks", Physical Review Letters, 1995 </blockquote>

<blockquote>12- C. Baldassi, A. Braunstein. "A Max-Sum algorithm for training discrete neural networks", Journal of Statistical Mechanics: Theory and Experiment, 2015 </blockquote>

<blockquote>13- G. Parisi. "Mean field theory of spin glasses: statics and dynamics", arXiv, 2007 </blockquote>

<blockquote>14- L. Dall'Asta, A. Ramezanpour, R. Zecchina. "Entropy landscape and non-Gibbs solutions in constraint satisfaction problem", Physical Review E, 2008 </blockquote>

<blockquote>15- M. Mézard, A. Montanari. "Information, Physics and Computation", Oxford Graduate Texts, 2009 </blockquote>

<blockquote>16- C. Baldassi, A. Ingrosso, C. Lucibello, L. Saglietti, R. Zecchina. "Subdominant Dense Clusters Allow for Simple Learning and High Computational Performance in Neural Networks with Discrete Synapses", Physical Review Letters, 2015 </blockquote>

## Authors

* <img src="https://avatars0.githubusercontent.com/u/24650975?s=400&v=4" width="25px"> **Nico Curti** [git](https://github.com/Nico-Curti), [unibo](https://www.unibo.it/sitoweb/nico.curti2)
* <img src="https://avatars3.githubusercontent.com/u/23407684?s=400&v=4" width="25px"> **Daniele Dall'Olio** [git](https://github.com/DanieleDallOlio), [unibo](https://www.unibo.it/sitoweb/daniele.dallolio)
* <img src="https://avatars2.githubusercontent.com/u/25343321?s=400&v=4" width="25px"> **Daniel Remondini** [git](https://github.com/dremondini), [unibo](https://www.unibo.it/sitoweb/daniel.remondini)
* <img src="https://www.unibo.it/uniboweb/utils/UserImage.aspx?IdAnagrafica=236217&IdFoto=bf094429" width="25px"> **Gastone Castellani** [unibo](https://www.unibo.it/sitoweb/gastone.castellani)
* <img src="https://avatars2.githubusercontent.com/u/1419337?s=400&v=4" width="25px;"/> **Enrico Giampieri** [git](https://github.com/EnricoGiampieri), [unibo](https://www.unibo.it/sitoweb/enrico.giampieri)

See also the list of [contributors](https://github.com/Nico-Curti/rFBP/contributors) [![GitHub contributors](https://img.shields.io/github/contributors/Nico-Curti/rFBP.svg?style=plastic)](https://github.com/Nico-Curti/rFBP/graphs/contributors/) who participated in this project.

## License

The `rFBP` package is licensed under the MIT "Expat" License. [![License](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/Nico-Curti/rFBP/blob/master/LICENSE)

## Acknowledgments

Thanks goes to all contributors of this project.

We thank also the author(s) of [Catch2](https://github.com/catchorg/Catch2) library: we have used it in the testing procedure of our C++ version and it is amazing!

## Citation

If you have found `rFBP` helpful in your research, please consider citing the paper

```BibTeX
@misc{DallOlioCCS19,
  author = {Dall'Olio, Daniele and Curti, Nico and Castellani, Gastone and Bazzani, Armando and Remondini, Daniel},
  title = {Classification of Genome Wide Association data by Belief Propagation Neural network},
  year = {2019},
  conference = {Conference of Complex System}
}
```

or just this project repository

```BibTeX
@misc{ReplicatedFocusingBeliefPropagation,
  author = {Curti, Nico and Dall'Olio, Daniele and Giampieri, Enrico},
  title = {Replicated Focusing Belief Propagation},
  year = {2019},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/Nico-Curti/rFBP}},
}
```
---
title: 'rFBP: Replicated Focusing Belief Propagation algorithm'
tags:
- belief-propagation
- deep-neural-networks
- sping-glass
- statistical-mechanics
- learning-algorithm
- machine-learning-algorithms
- python3
- cpp17
authors:
- name: Nico Curti^[co-first author]
  orcid: 0000-0001-5802-1195
  affiliation: 1
- name: Daniele Dall'Olio^[co-first author]
  orcid: 0000-0003-0196-6870
  affiliation: 3
- name: Daniel Remondini
  orcid: 0000-0003-3185-7456
  affiliation: 3
- name: Gastone Castellani
  orcid: 0000-0003-4892-925X
  affiliation: 2
- name: Enrico Giampieri
  orcid: 0000-0003-2269-2338
  affiliation: 1
affiliations:
- name: eDIMESLab, Department of Experimental, Diagnostic and Specialty Medicine of Bologna University
  index: 1
- name: Department of Experimental, Diagnostic and Specialty Medicine of Bologna University
  index: 2
- name: Department of Physics and Astronomy of Bologna University
  index: 3
date: 06 October 2020
bibliography: paper.bib
---

# Summary

The `rFBP` project implements a `scikit-learn` compatible machine-learning binary classifier leveraging fully connected neural networks with a learning algorithm (*Replicated Focusing Belief Propagation*, rFBP) that is quickly converging and robust (less prone to brittle overfitting) for ill-posed datasets (very few samples compared to the number of features).
The current implementation works only with binary features such as one-hot encoding for categorical data.

This library has already been widely used to successfully predict *source attribution* starting from GWAS (*Genome Wide Association Studies*) data.
That study was trying to predict the animal origin for an infectious bacterial disease inside the H2020 European project COMPARE (Grant agreement ID: 643476).
A full description of the pipeline used in this study is available in the abstract and slides provided into the [publications](https://github.com/Nico-Curti/rFBP/blob/master/publications) folder of the project.

Algorithm application on real data:

- *Classification of Genome Wide Association data by Belief Propagation Neural network*, [CCS Italy 2019](https://github.com/Nico-Curti/rFBP/blob/master/publications/conference/ccs19.pdf), Conference paper

- *Classification of Genome Wide Association data by Belief Propagation Neural network*, [CCS Italy 2019](https://github.com/Nico-Curti/rFBP/blob/master/publications/presentation/ccs19.pdf), Conference slides

# Statement of need

The learning problem under ill-posed conditions can be tackled through statistical mechanic models joined with the so-called Large Deviation Theory [@parisi2007mean; @Baldassi_2015; @Monasson_1995; @Zecchina_1995; @Baldassi_2016_local].
In general, the learning problem can be split into two sub-parts: the classification problem and the generalization one.
The first aims to completely store a pattern sample, i.e., a prior known ensemble of input-output associations (*perfect learning*, Baldassi et al.) [@Baldassi_2016; @Krauth1989StorageCO].
The second one corresponds to compute a discriminant function based on a set of features of the input which guarantees a unique association of a pattern.

From a statistical point-of-view many Neural Network models have been proposed and spin-glass models have emerged as the most promising ones.
Starting from a balanced distribution of the system, generally based on Boltzmann distribution, and under proper conditions, we can prove that the classification problem becomes a NP-complete computational problem [@Blum_1992].
A wide range of heuristic solutions to that type of problems were proposed [@Huang_2014; @Braunstein_2006; @Baldassi_11079].

In this project we show one of these algorithms developed by Baldassi et al. [@Baldassi_2016] and called *Replicated Focusing Belief Propagation* (`rFBP`).
The `rFBP` algorithm is a learning algorithm developed to justify the learning process of a binary neural network framework.
The model is based on a spin-glass distribution of neurons put on a fully connected neural network architecture.
In this way each neuron is identified by a spin and so only binary weights (-1 and 1) can be assumed by each entry.
The learning rule which controls the weight updates is given by the Belief Propagation method.

A first implementation of the algorithm was proposed in the original paper [@Baldassi_2016] jointly with an open-source Github repository.
The original version of the algorithm was written in [`Julia` language](https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl).
`Julia` is certainly an efficient programming language but it is not part of most machine learning developers' tool of choice.
To broaden the scope and use of the method, a `C++` implementation was developed with a joint `Cython` wrap for `Python` users.
The `C++` language guarantees better computational performances against the `Julia` implementation and the `Python` version enhances its usability.
This implementation is optimized for parallel computing and is endowed with a custom `C++` library called [`Scorer`](https://github.com/Nico-Curti/scorer), which is able to compute a large number of statistical measurements based on a hierarchical graph scheme.
With this optimized implementation and its [`scikit-learn`](https://github.com/scikit-learn/scikit-learn) compatibility we try to encourage researchers to approach these alternative algorithms and to use them more frequently on real context.

As the `Julia` implementation also the `C++` one provides the entire `rFBP` framework in a single library callable via a command line interface.
The library widely uses template syntaxes to perform dynamic specialization of the methods between two magnetization versions of the algorithm.
The main object categories needed by the algorithm are wrapped in handy `C++` objects easy to use also from the `Python` interface.

# Acknowledgments

The authors acknowledge COMPARE n. 643476 EU Horizon 2020 (EU) Project.

# References
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at nico.curti2@unibo.it. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contribution

Any contribution is more than welcome :heart:. Just fill an [issue](https://github.com/Nico-Curti/rFBP/blob/master/.github/ISSUE_TEMPLATE/ISSUE_TEMPLATE.md) or a [pull request](https://github.com/Nico-Curti/rFBP/blob/master/.github/PULL_REQUEST_TEMPLATE/PULL_REQUEST_TEMPLATE.md) and we will check ASAP!

## Before start contributing

- Make sure you agree to contribute your code under `rFBP` [license](https://github.com/Nico-Curti/rFBP/blob/master/LICENSE)

- If you are going to fix a bug, check that it's still exists and that you are working with the latest version of `rFBP` library.

- Make sure that there is not someone working on the same issue. In this case you can provide support or suggestion in the issue.

- If you have any question about the library fill an issue with the same criteria described in the previous step.

## Pull Request and Contribution

Please consider the following step to create a pull request which can be processed as quickly as possible

1. Install Git.

2. Create your account on Github.

3. Fork [this](https://github.com/Nico-Curti/rFBP) repository.

4. Create your feature branch (under `dev` branch).

5. Choose a task for yourself. It could be a bugfix or some new code.

6. Add documentation (`docstring`) to your functions/methods if you are working on the Python version. Add comments about what you are doing in functions/methods if you are working on the C++ version.

7. Add tests for your contributions

9. Pass **ALL** CI tests.

10. Submit a Pull Request into `dev` (following the pull request [template](https://github.com/Nico-Curti/rFBP/blob/master/.github/PULL_REQUEST_TEMPLATE/PULL_REQUEST_TEMPLATE.md))

## Merging

1. As soon as possible the authors will review your pull request, answering to your message.

2. Make sure you got credits at the end of merge: we have a lot of project to follow and we may forget to add your name to the list of `AUTHORS.md`. Please, do not hesitate to remind us!
<!-- Please use this line to close one or multiple issues when this pullrequest gets merged
You can add another line right under the first one:
resolves #1234
resolves #1235
-->

#### This PR changes :grey_question:
<!-- Explain your changes -->

#### Any other comments?
<!--
This is a template helping you to create an issue which can be processed as quickly as possible. This is the bug reporting section for the rFBP library.
-->

#### Detailed Description :grey_question:
<!-- your description -->

#### Steps/Code to Reproduce :grey_question:
<!-- to add code example fence it with triple backticks and optional file extension
    ```.cpp
    // C++ code example
    ```
 or attach as .txt or .zip file
-->

#### Expected Behavior :grey_question:
<!-- Description of the expected result(s) -->

#### Actual Behavior :grey_question:
<!-- Description (possibly with some shell reports) of the actual result(s) -->

#### Operating System / Platform :grey_question:
<!-- Example
- OS: Windows 10 Pro
- System type: x64
- Processor: i7-6500U
- RAM: 8 GB
-->

#### Compiler (with version) :grey_question:
<!--
$ g++ --version
g++-9 (Ubuntu 9.2.1-17ubuntu1~16.04) 9.2.1 20191102
-->

#### Python Version :grey_question:
<!--
$ python --version
Python 3.7.3
-->

#### rFBP Version (`rFBP.__version__`) :grey_question:
<!--
$ python -c "import ReplicatedFocusingBeliefPropagation as rFBP; print(rFBP.__version__)"
'1.0.0'
-->
# Replicated Focusing Belief Propagation Documents

----------

In this folder we store all the papers and reports related to this library.

If we can not share source codes for copyright reasons we provide a link to them.

## Conference Paper

- *Classification of Genome Wide Association data by Belief Propagation Neural network*, [CCS Italy 2019](https://github.com/Nico-Curti/rFBP/blob/master/publications/conference/ccs19.pdf)

## Presentations

- *Classification of Genome Wide Association data by Belief Propagation Neural network*, [CCS Italy 2019](https://github.com/Nico-Curti/rFBP/blob/master/publications/presentation/ccs19.pdf)

## Thesis

- **ITA**: *Applicazione di un algoritmo d'apprendimento basato su sistemi fuori dall'equilibrio a dati di Genome Wide Association*, Master's degree thesis, Author: Dott. Daniele Dall'Olio, Supervisor: Prof. Gastone Castellani, Correlator: Dott. Nico Curti

| **Authors**  | **Project** |  **Build Status** | **Code Quality** | **Coverage** |
|:------------:|:-----------:|:-----------------:|:----------------:|:------------:|
| [**N. Curti**](https://github.com/Nico-Curti) <br/> [**D. Dall'Olio**](https://github.com/DanieleDallOlio) <br/> [**E. Giampieri**](https://github.com/EnricoGiampieri)  |  **rFBP** <br/> [![JOSS](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0/status.svg)](https://joss.theoj.org/papers/7643779111039dbc7776ff49d2a6b1b0) <br/> [![PyPI Latest Release](https://img.shields.io/pypi/v/ReplicatedFocusingBeliefPropagation.svg)](https://pypi.org/project/ReplicatedFocusingBeliefPropagation/)  | **Linux/MacOS** : [![Travis](https://travis-ci.com/Nico-Curti/rFBP.svg?token=7QqsqaQiuDHSyGDT3xek&branch=master)](https://travis-ci.com/Nico-Curti/rFBP) <br/> **Windows** : [![appveyor](https://ci.appveyor.com/api/projects/status/obuq56lhyd90pmup?svg=true)](https://ci.appveyor.com/project/Nico-Curti/rfbp) | **Codacy** : [![Codacy](https://api.codacy.com/project/badge/Grade/a6fdac990b6f4141a5bd9e8171ddaf53)](https://www.codacy.com/manual/Nico-Curti/rFBP?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Nico-Curti/rFBP&amp;utm_campaign=Badge_Grade) <br/> **Codebeat** : [![Codebeat](https://codebeat.co/badges/cc761a7c-79fa-4a66-984f-bef6fd145d34)](https://codebeat.co/projects/github-com-nico-curti-rfbp-master) | [![codecov](https://codecov.io/gh/Nico-Curti/rFBP/branch/master/graph/badge.svg)](https://codecov.io/gh/Nico-Curti/rFBP) |

[![rFBP C++ CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20C++%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+C%2B%2B+CI%22)
[![rFBP Python CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20Python%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+Python+CI%22)
[![rFBP Docs CI](https://github.com/Nico-Curti/rFBP/workflows/rFBP%20Docs%20CI/badge.svg)](https://github.com/Nico-Curti/rFBP/actions?query=workflow%3A%22rFBP+Docs+CI%22)

[![docs](https://readthedocs.org/projects/rfbp/badge/?version=latest)](https://rfbp.readthedocs.io/en/latest/?badge=latest)
[![GitHub pull-requests](https://img.shields.io/github/issues-pr/Nico-Curti/rFBP.svg?style=plastic)](https://github.com/Nico-Curti/rFBP/pulls)
[![GitHub issues](https://img.shields.io/github/issues/Nico-Curti/rFBP.svg?style=plastic)](https://github.com/Nico-Curti/rFBP/issues)

[![GitHub stars](https://img.shields.io/github/stars/Nico-Curti/rFBP.svg?label=Stars&style=social)](https://github.com/Nico-Curti/rFBP/stargazers)
[![GitHub watchers](https://img.shields.io/github/watchers/Nico-Curti/rFBP.svg?label=Watch&style=social)](https://github.com/Nico-Curti/rFBP/watchers)

<a href="https://github.com/UniboDIFABiophysics">
  <div class="image">
    <img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90">
  </div>
</a>

# Replicated Focusing Belief Propagation algorithm

We propose a `C++` version of the [**Replicated Focusing Belief Propagation**](https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl) Julia package.
Our implementation optimizes and extends the original library including multi-threading support and an easy-to-use interface to the main algorithm.
To further improve the usage of our code, we propose also a `Python` wrap of the library with a full compatibility with the [`scikit-learn`](https://github.com/scikit-learn/scikit-learn) and [`scikit-optimize`](https://github.com/scikit-optimize/scikit-optimize) packages.

- [C++ Install](./cpp_install.md)
- [Python Install](./python_install.md)
- [Model description](./model.md)
- [Authors](./authors.md)
# Authors

## Nico Curti

<img align="right" width="180" height="183" src="https://avatars0.githubusercontent.com/u/24650975?s=400&v=4"> About me:
- **Name**: Nico
- **Surname**: Curti
- **Profession**: PhD at eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University
- **University**: University of Bologna
- **Location**: Italy
- **Web page**: [git](https://github.com/Nico-Curti), [unibo](https://www.unibo.it/sitoweb/nico.curti2)
- **Contact me**: [email](mailto:nico.curti2@unibo.it)

## Daniele Dall'Olio

<img align="right" width="180" height="183" src="https://avatars3.githubusercontent.com/u/23407684?s=400&v=4"> About me:
- **Name**: Daniele
- **Surname**: Dall'Olio
- **Profession**: PhD student at Physics and Anstronomy Department in Bologna
- **University**: University of Bologna
- **Location**: Italy
- **Web page**: [git](https://github.com/DanieleDallOlio), [unibo](https://www.unibo.it/sitoweb/daniele.dallolio)
- **Contact me**: [email](mailto:daniele.dallolio@studio.unibo.it)

## Enrico Giampieri

<img align="right" width="180" height="183" src="https://avatars2.githubusercontent.com/u/1419337?s=400&v=4"> About me:
- **Name**: Enrico
- **Surname**: Giampieri
- **Profession**: PhD at eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine
- **University**: University of Bologna
- **Location**: Italy
- **Web page**: [git](https://github.com/EnricoGiampieri), [unibo](https://www.unibo.it/sitoweb/enrico.giampieri)
- **Contact me**: [email](mailto:enrico.giampieri@unibo.it)
# Replicated Focusing Belief Propagation

The *Replicated Focusing Belief Propagation* (`rFBP`) is an entropy-maximization based algorithm operating as the underlying learning rule of feed-forward binary neural networks.
Here, the entropy is defined as:

<img src="https://render.githubusercontent.com/render/math?math={S(\vec{w},\beta,\gamma) = \frac{1}{N} log \bigg ( \sum_{\vec{w}'}  e^{-\beta E(\vec{w}')} e^{\gamma \vec{w}' \cdot \vec{w}}	\bigg)}"> ,

where <img src="https://render.githubusercontent.com/render/math?math={\vec{w}}"> is the whole weights set, N is its size and <img src="https://render.githubusercontent.com/render/math?math={\beta}"> is the standard Boltzmann term.
Further, <img src="https://render.githubusercontent.com/render/math?math={E(\vec{w})}"> is the energy of such weights set, which is equal to the number of wrong predictions on the training set produced by <img src="https://render.githubusercontent.com/render/math?math={\vec{w}}">.
When <img src="https://render.githubusercontent.com/render/math?math={\beta \to inf}">, only those <img src="https://render.githubusercontent.com/render/math?math={\vec{w}'}"> with null energy, i.e. perfect solutions, sum up in the entropy.
At the time being, the rFBP only works with <img src="https://render.githubusercontent.com/render/math?math={\beta \to inf}">.

The realization of the rFBP is equivalent to evolve a spin-glass system composed of several interacting replicas with the Belief Propagation algorithm. The density distribution of the system is modelled by:

<img src="https://render.githubusercontent.com/render/math?math={P(\vec{w} | \beta,y,\gamma) = \frac{e^{yN S(\vec{w},\beta,\gamma)}}{Z(\beta,y,\gamma)}}"> ,

where Z is the partition function.

Such spin-glass model depends necessarily on two parameters: y and <img src="https://render.githubusercontent.com/render/math?math={\gamma}">.
The former is a temperature-alike related variable, similar to the one usually exploited by Gradient Descend approaches, but it can be also interpreted as the number of interacting replicas of the system.
The latter is the penalization term associated to the distance between two weights sets. Indeed, the term <img src="https://render.githubusercontent.com/render/math?math={e^{\gamma \vec{w}' \cdot \vec{w}}}"> in the entropy is larger, when <img src="https://render.githubusercontent.com/render/math?math={\vec{w}}"> and <img src="https://render.githubusercontent.com/render/math?math={\vec{w}'}"> are closer.

The Belief Propagation algorithm needs to be to adjusted by adding incoming extra messages for all weights, in order to involve the interacting replicas of the system.
This extra term is represented by:

<img src="https://render.githubusercontent.com/render/math?math={\hat{m}^{t_1}_{\star \to \w_i} = tanh \big[ (y-1) artanh ( m^{t_0}_{\w_i \to \star} tanh \gamma ) \big] tanh \gamma}">,

where <img src="https://render.githubusercontent.com/render/math?math={\w_i}"> and <img src="https://render.githubusercontent.com/render/math?math={\star}"> stand respectively for the i-th weight and a representation of all i-th replicas.

The `rFBP` is therefore an adjusted Belief Propagation algorithm, whose general procedure can be summarized as follows:
- set <img src="https://render.githubusercontent.com/render/math?math={\beta \to inf}">
- select protocols for y and <img src="https://render.githubusercontent.com/render/math?math={\gamma}">;
- set first values of y and <img src="https://render.githubusercontent.com/render/math?math={\gamma}"> and run the adjusted-BP method until convergence (<img src="https://render.githubusercontent.com/render/math?math={ < \epsilon}">) or up to a limited-number of iterations;
- step to the next pair values of y and <img src="https://render.githubusercontent.com/render/math?math={\gamma}"> with respect to the chosen protocols and re-run the adjusted-BP method;
- keep it going until a solution is reached or protocols end.

The `rFBP` algorithm focuses the replicated system to fall step by step into weights sets extremely closed to many perfect solutions (<img src="https://render.githubusercontent.com/render/math?math={\vec{w}}"> such that <img src="https://render.githubusercontent.com/render/math?math={E(\vec{w})=0}">), which ables them to well generalize out of the training set [[Baldassi101103](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.115.128101)].
## C++ Installation

1) Follow your system prerequisites (below)

2) Clone the `rFBP` package from this repository, or download a stable release

```bash
git clone https://github.com/Nico-Curti/rFBP.git
cd rFBP
```

3) `rFBP` could be built with CMake and Make or with the *build* scripts in the project.
Example:

**Unix OS:**
```bash
./build.sh
```

**Windows OS:**
```Powershell
PS \>                 ./build.ps1
```

## Installation

### Ubuntu

1) Define a work folder, which we will call WORKSPACE in this tutorial: this could be a "Code" folder in our home, a "c++" folder on our desktop, whatever you want. Create it if you don't already have, using your favourite method (mkdir in bash, or from the graphical interface of your distribution). We will now define an environment variable to tell the system where our folder is. Please note down the full path of this folder, which will look like `/home/$(whoami)/code/`

```bash
echo -e "\n export WORKSPACE=/full/path/to/my/folder \n" >> ~/.bashrc
source ~/.bashrc
```

2) Open a Bash terminal and type the following commands to install all the prerequisites.

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y gcc-8 g++-8

wget --no-check-certificate https://cmake.org/files/v3.13/cmake-3.13.1-Linux-x86_64.tar.gz
tar -xzf cmake-3.13.1-Linux-x86_64.tar.gz
export PATH=$PWD/cmake-3.13.1-Linux-x86_64/bin:$PATH

sudo apt-get install -y make git dos2unix ninja-build
git config --global core.autocrlf input
git clone https://github.com/physycom/sysconfig
```

3) Build the project with CMake (enable or disable OMP with the define **-DOMP**; enable or disable **Cython** building with the define **-DPYWRAP**; enable or disable the **Scorer** support using **-DSCORER**; enable or disable testing with the define **-DBUILD_TEST**):

```bash
cd $WORKSPACE
git clone https://github.com/Nico-Curti/rFBP
cd rFBP

mkdir -p build
cd build

cmake -DOMP=ON ..
make -j
cmake --build . --target install
cd ..
```

### macOS

1) If not already installed, install the XCode Command Line Tools, typing this command in a terminal:

```bash
xcode-select --install
```

2) If not already installed, install Homebrew following the [official guide](https://brew.sh/index_it.html).

3) Open the terminal and type these commands

```bash
brew update
brew upgrade
brew install gcc@8
brew install cmake make git ninja
```

4) Define a work folder, which we will call WORKSPACE in this tutorial: this could be a "Code" folder in our home, a "c++" folder on our desktop, whatever you want. Create it if you don't already have, using your favourite method (mkdir in bash, or from the graphical interface in Finder). We will now define an environment variable to tell the system where our folder is. Please note down the full path of this folder, which will look like /home/$(whoami)/code/

5) Open a Terminal and type the following command (replace /full/path/to/my/folder with the previous path noted down)

```bash
echo -e "\n export WORKSPACE=/full/path/to/my/folder \n" >> ~/.bash_profile
source ~/.bash_profile
```

6) Build the project with CMake (enable or disable OMP with the define **-DOMP**; enable or disable **Cython** building with the define **-DPYWRAP**; enable or disable the **Scorer** support using **-DSCORER**; enable or disable testing with the define **-DBUILD_TEST**):

```bash
cd $WORKSPACE
git clone https://github.com/Nico-Curti/rFBP
cd rFBP

mkdir -p build
cd build

cmake -DOMP=ON ..
make -j
cmake --build . --target install
cd ..
```

### Windows (7+)

1) Install Visual Studio 2017 from the [official website](https://www.visualstudio.com/)

2) Open your Powershell with Administrator privileges, type the following command and confirm it:

```PowerShell
PS \>                 Set-ExecutionPolicy unrestricted
```

3) If not already installed, please install chocolatey using the [official guide](http://chocolatey.org)

4) If you are not sure about having them updated, or even installed, please install `git`, `cmake` and an updated `Powershell`. To do so, open your Powershell with Administrator privileges and type

```PowerShell
PS \>                 cinst -y git cmake powershell
```

5) Restart the PC if required by chocolatey after the latest step

6) Install PGI 18.10 from the [official website](https://www.pgroup.com/products/community.htm) (the community edition is enough and is free; NOTE: install included MS-MPI, but avoid JRE and Cygwin)

7) Activate license for PGI 18.10 Community Edition (rename the file `%PROGRAMFILES%\PGI\license.dat-COMMUNITY-18.10` to `%PROGRAMFILES%\PGI\license.dat`) if necessary, otherwise enable a Professional License if available

8) Define a work folder, which we will call `WORKSPACE` in this tutorial: this could be a "Code" folder in our home, a "cpp" folder on our desktop, whatever you want. Create it if you don't already have, using your favourite method (mkdir in Powershell, or from the graphical interface in explorer). We will now define an environment variable to tell the system where our folder is. Please note down its full path. Open a Powershell (as a standard user) and type

```PowerShell
PS \>                 rundll32 sysdm.cpl,EditEnvironmentVariables
```

9) In the upper part of the window that pops-up, we have to create a new environment variable, with name `WORKSPACE` and value the full path noted down before.
If it not already in the `PATH` (this is possible only if you did it before), we also need to modify the "Path" variable adding the following string (on Windows 10 you need to add a new line to insert it, on Windows 7/8 it is necessary to append it using a `;` as a separator between other records):

```cmd
                      %PROGRAMFILES%\CMake\bin
```

10) If `vcpkg` is not installed, please follow the next procedure, otherwise please jump to #12

```PowerShell
PS \>                 cd $env:WORKSPACE
PS Code>              git clone https://github.com/Microsoft/vcpkg.git
PS Code>              cd vcpkg
PS Code\vcpkg>        .\bootstrap-vcpkg.bat
```

11) Open a Powershell with Administrator privileges and type

```PowerShell
PS \>                 cd $env:WORKSPACE
PS Code>              cd vcpkg
PS Code\vcpkg>        .\vcpkg integrate install
```

12) Open a Powershell and build `rFBP` using the `build.ps1` script

```PowerShell
PS \>                 cd $env:WORKSPACE
PS Code>              git clone https://github.com/Nico-Curti/rFBP
PS Code>              cd rFBP
PS Code\rFBP>         .\build.ps1
```
## Python Installation

First of all ensure that a right Python version is installed (Python >= 3.5 is required).
The [Anaconda/Miniconda](https://www.anaconda.com/) python version is recomended.

Download the project or the latest release:

```bash
git clone https://github.com/Nico-Curti/rFBP
cd rFBP
```

### Installing prerequisites

To install the prerequisites type:

```bash
pip install -r ./requirements.txt
```

### Installation from sources

In the `rFBP` directory execute:

```bash
python setup.py install
```

or for installing in development mode:

```bash
python setup.py develop --user
```

### Installation using pip

The latest release of the `rFBP` package can be installed using `pip`

```bash
pip install ReplicatedFocusingBeliefPropagation
```

The installation via `pip` requires to pre-install the `Cython` and `Numpy` packages, thus make sure to pre-install them!Theory
-------

The `rFBP` algorithm derives from an out-of-equilibrium (non-Boltzmann) model of the learning process of binary neural networks DallAsta101103_.
This model mimics a spin glass system whose realizations are equally likely to occur when sharing the same so-called entropy (not the same energy, i.e. out-of-equilibrium).
This entropy basically counts the number of solutions (zero-energy realizations) around a realization below a fixed-distance.

Within this out-of-equilibrium framework, the objective is to maximize the entropy instead of minimizing the energy.
From a machine learning standpoint, we aim at those weights sets that perfectly solve the learning process (zero-errors) and that are mathematically closed to each other.
To this end, the Belief Propagation method MézardMontanari_ can be adopted as the underlying learning rule, although it must be properly adjusted to take into account the out-of-equilibrium nature of the model.

The *Replicated Focusing Belief Propagation* (`rFBP`) is an entropy-maximization based algorithm operating as the underlying learning rule of feed-forward binary neural networks.
Here, the entropy is defined as:

:math:`S(\vec{w},\beta,\gamma) = \frac{1}{N} log \bigg ( \sum_{\vec{w}'}  e^{-\beta E(\vec{w}')} e^{\gamma \vec{w}' \cdot \vec{w}}	\bigg)`

where :math:`\vec{w}` is the whole weights set, N is its size and :math:`\beta` is the standard Boltzmann term.
Further, :math:`E(\vec{w})` is the energy of such weights set, which is equal to the number of wrong predictions on the training set produced by :math:`\vec{w}`.
When :math:`\beta \to inf`, only those :math:`\vec{w}` with null energy, i.e. perfect solutions, sum up in the entropy.
At the time being, the rFBP only works with :math:`\beta \to inf`.

The realization of the rFBP is equivalent to evolve a spin-glass system composed of several interacting replicas with the Belief Propagation algorithm. The density distribution of the system is modelled by:

:math:`P(\vec{w} | \beta,y,\gamma) = \frac{e^{yN S(\vec{w},\beta,\gamma)}}{Z(\beta,y,\gamma)}`

where Z is the partition function.

Such spin-glass model depends necessarily on two parameters: y and :math:`\gamma`.
The former is a temperature-alike related variable, similar to the one usually exploited by Gradient Descend approaches, but it can be also interpreted as the number of interacting replicas of the system.
The latter is the penalization term associated to the distance between two weights sets. Indeed, the term :math:`e^{\gamma \vec{w}' \cdot \vec{w}}` in the entropy is larger, when :math:`\vec{w}` and :math:`\vec{w}` are closer.

The Belief Propagation algorithm needs to be to adjusted by adding incoming extra messages for all weights, in order to involve the interacting replicas of the system.
This extra term is represented by:

:math:`\hat{m}^{t_1}_{\star \to w_i} = tanh \big[ (y-1) artanh ( m^{t_0}_{w_i \to \star} tanh \gamma ) \big] tanh \gamma`

where :math:`w_i` and :math:`\star` stand respectively for the i-th weight and a representation of all i-th replicas.

The `rFBP` is therefore an adjusted Belief Propagation algorithm, whose general procedure can be summarized as follows:

- set :math:`\beta \to inf`
- select protocols for y and :math:`\gamma`
- set first values of y and :math:`\gamma` and run the adjusted-BP method until convergence (:math:`\epsilon` or up to a limited-number of iterations;
- step to the next pair values of y and :math:`\gamma` with respect to the chosen protocols and re-run the adjusted-BP method;
- keep it going until a solution is reached or protocols end.

The `rFBP` algorithm focuses the replicated system to fall step by step into weights sets extremely closed to many perfect solutions (:math:`\vec{w}` such that :math:`E(\vec{w})=0`), which ables them to well generalize out of the training set `Zecchina et al`_.

.. _`Zecchina et al`: https://www.pnas.org/content/113/48/E7655
.. _MézardMontanari: https://web.stanford.edu/~montanar/RESEARCH/book.html
.. _DallAsta101103 : https://journals.aps.org/pre/abstract/10.1103/PhysRevE.77.031118
CMake C++ Installation
======================

We recommend to use `CMake` for the installation since it is the most automated way to reach your needs.
First of all make sure you have a sufficient version of `CMake` installed (3.9 minimum version required).
If you are working on a machine without root privileges and you need to upgrade your `CMake` version a valid solution to overcome your problems is provided shut_.

With a valid `CMake` version installed first of all clone the project as:

.. code-block:: bash

	git clone https://github.com/Nico-Curti/rFBP
	cd rFBP


The you can build the `rFBP` package with

.. code-block:: bash

	mkdir -p build
	cd build && cmake .. && cmake --build . --target install

or more easily

.. code-block:: bash

	./build.sh

if you are working on a Windows machine the right script to call is the `build.ps1`_.

.. note::
  If you want enable the OpenMP support (*4.5 version is required*) compile the library with `-DOMP=ON`.

.. note::
	If you want enable the Scorer support compile the library with `-DSCORER=ON`. If you want use a particular installation of the Scorer library or you have manually installed the library following the `README` instructions, we suggest to add the `-DScorer_DIR=/path/to/scorer/shared/scorer` in the command line.

.. note::
	If you want enable the Cython support compile the library with `-DPYWRAP=ON`. The Cython packages will be compiled and correctly positioned in the `rFBP` Python package **BUT** you need to run also the setup before use it.

.. note::
	If you use MagT configuration, please download the `atanherf coefficients` file before running any executable. You can find a downloader script inside the scripts_ folder. Enter in that folder and just run `python dowload_atanherf.py`.

.. _shut: https://github.com/Nico-Curti/Shut
.. _`build.ps1`: https://Nico-Curti/rFBP/blob/master/build.ps1
.. _scripts: https://github.com/Nico-Curti/rFBP/tree/master/scriptsInstallation guide
==================

C++ supported compilers:

|gcc version|

|clang version|

|msvc version|

The `rFBP` project is written in `C++` using a large amount of c++17 features.
To enlarge the usability of our package we provide also a retro-compatibility of all the c++17 modules reaching an usability (tested) of our code from gcc 4.8.5+.
The package installation can be performed via `CMake` or `Makefile`.

If you are using the `CMake` (recommended) installer the maximum version of C++ standard is automatic detected.
The `CMake` installer provides also the export of the library: after the installation you can use this library into other `CMake` projects using a simple `find_package` function.
The exported `CMake` library (`rFBP::rfbp`) is installed in the `share/rFBP` directory of the current project and the relative header files are available in the `rFBP_INCLUDE_DIR` variable.

The `CMake` installer provides also a `rFBP.pc`, useful if you want link to the `rFBP` using `pkg-config`.

You can also use the `rFBP` package in `Python` using the `Cython` wrap provided inside this project.
The only requirements are the following:

* numpy >= 1.15
* cython >= 0.29
* scipy >= 1.2.1
* scikit-learn >= 0.20.3
* requests >= 2.22.0

The `Cython` version can be built and installed via `CMake` enabling the `-DPYWRAP` variable.
The `Python` wrap guarantees also a good integration with the other common Machine Learning tools provided by `scikit-learn` `Python` package; in this way you can use the `rFBP` algorithm as an equivalent alternative also in other pipelines.
Like other Machine Learning algorithm also the `rFBP` one depends on many parameters, i.e its hyper-parameters, which has to be tuned according to the given problem.
The `Python` wrap of the library was written according to `scikit-optimize` `Python` package to allow an easy hyper-parameters optimization using the already implemented classical methods.


.. |gcc version| image:: https://img.shields.io/badge/gcc-4.8.5%20|%204.9.*%20|%205.*%20|%206.*%20|%207.*%20|%208.*%20|%209.*-yellow.svg
.. |clang version| image:: https://img.shields.io/badge/clang-3.6%20|3.9%20|5.*%20|%206.*%20|%207.*%20|-red.svg
.. |msvc version| image:: https://img.shields.io/badge/msvc-vs2017%20x86%20|%20vs2017%20x64|%20vs2019%20x86%20|%20vs2019%20x64-blue.svg
.. _CMake: https://github.com/Nico-Curti/rFBP/blob/master/CMakeLists.txt
.. _Makefile: https://github.com/Nico-Curti/rFBP/blob/master/Makefile

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   CMake
   Python
Python Installation
===================

Python version supported : |Python version|

The easiest way to install the package is using `pip`

.. code-block:: bash

	python -m pip install ReplicatedFocusingBeliefPropagation

.. warning::

	The setup file requires the `Cython` and `Numpy` packages, thus make sure to pre-install them!
	We are working on some workarounds to solve this issue.

The `Python` installation can be performed with or without the `C++` installation.
The `Python` installation is always executed using `setup.py`_ script.

If you have already built the `rFBP` `C++` library the installation is performed faster and the `Cython` wrap was already built using the `-DPYWRAP` definition.
Otherwise the full list of dependencies is build.

In both cases the installation steps are

.. code-block:: bash

	python -m pip install -r ./requirements.txt

to install the prerequisites and then

.. code-block:: bash

	python setup.py install

or for installing in development mode:

.. code-block:: bash

	python setup.py develop --user

.. warning::

	The current installation via pip has no requirements about the version of `setuptools` package.
	If the already installed version of `setuptools` is `>= 50.*` you can find some troubles during the installation of our package (ref. issue_).
	We suggest to temporary downgrade the `setuptools` version to `49.3.0` to workaround this `setuptools` issue.


.. |Python version| image:: https://img.shields.io/badge/python-3.5|3.6|3.7|3.8-blue.svg
.. _`setup.py`: https://github.com/Nico-Curti/blob/master/setup.py
.. _issue: https://github.com/Nico-Curti/rFBP/issues/5.. Replicated Focusing Belief Propagation algorithm documentation master file, created by
   sphinx-quickstart on Fri Oct  2 12:42:24 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Replicated Focusing Belief Propagation's documentation!
==================================================================

The **Replicated Focusing Belief Propagation** package is inspired by the original BinaryCommitteeMachineFBP_ verion written in Julia.
In our implementation we optimize and extend the original library inclu multi-threading support and an easy-to-use interface to the main algorithm.
To further improve the usage of our code, we propose also a `Python` wrap of the library with a full compatibility with the `scikit-learn`_ and `scikit-optimize`_ packages.

Overview
========

The learning problem could be faced through statistical mechanic models joined with the so-called Large Deviation Theory.
In general, the learning problem can be split into two sub-parts: the classification problem and the generalization one.
The first aims to completely store a pattern sample, i.e a prior known ensemble of input-output associations (*perfect learning*).
The second one corresponds to compute a discriminant function based on a set of features of the input which guarantees a unique association of a pattern.

From a statistical point-of-view many Neural Network models have been proposed and the most promising seems to be spin-glass models based.
Starting from a balanced distribution of the system, generally based on Boltzmann distribution, and under proper conditions, we can prove that the classification problem becomes a NP-complete computational problem.
A wide range of heuristic solutions to that type of problems were proposed.

In this project we show one of these algorithms developed by `Zecchina et al`_ and called *Replicated Focusing Belief Propagation* (`rFBP`).
The `rFBP` algorithm is a learning algorithm developed to justify the learning process of a binary neural network framework.
The model is based on a spin-glass distribution of neurons put on a fully connected neural network architecture.
In this way each neuron is identified by a spin and so only binary weights (-1 and 1) can be assumed by each entry.
The learning rule which controls the weight updates is given by the Belief Propagation method.

A first implementation of the algorithm was proposed in the original paper (`Zecchina et al`_) jointly with an open-source Github repository.
The original version of the code was written in `Julia` language and despite it is a quite efficient implementation the `Julia` programming language stays on difficult and far from many users.
To broaden the scope and use of the method, a `C++` implementation was developed with a joint `Cython` wrap for `Python` users.
The `C++` language guarantees better computational performances against the `Julia` implementation and the `Python` version enhances its usability.
This implementation is optimized for parallel computing and is endowed with a custom `C++` library called scorer_, which is able to compute a large number of statistical measurements based on a hierarchical graph scheme.
With this optimized implementation and its `scikit-learn`_ compatibility we try to encourage researchers to approach these alternative algorithms and to use them more frequently on real context.

As the `Julia` implementation also the `C++` one provides the entire `rFBP` framework in a single library callable via a command line interface.
The library widely uses template syntaxes to perform dynamic specialization of the methods between two magnetization versions of the algorithm.
The main object categories needed by the algorithm are wrapped in handy `C++` objects easy to use also from the `Python` interface.


Usage example
-------------

The `rfbp` object is totally equivalent to a `scikit-learn` classifier and thus it provides the member functions `fit` (to train your model) and `predict` (to test a trained model on new samples).

.. code-block:: python

   import numpy as np
   from sklearn.model_selection import train_test_split
   from ReplicatedFocusingBeliefPropagation import MagT64
   from ReplicatedFocusingBeliefPropagation import Pattern
   from ReplicatedFocusingBeliefPropagation import ReplicatedFocusingBeliefPropagation as rFBP

   N, M = (20, 101) # M must be odd
   X = np.random.choice([-1, 1], p=[.5, .5], size=(N, M))
   y = np.random.choice([-1, 1], p=[.5, .5], size=(N, ))

   X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

   rfbp = rFBP(mag=MagT64,
               hidden=3,
               max_iter=1000,
               seed=135,
               damping=0.5,
               accuracy=('accurate','exact'),
               randfact=0.1,
               epsil=0.5,
               protocol='pseudo_reinforcement',
               size=101,
               nth=1)

   rfbp.fit(X_train, y=y_train)
   y_pred = rfbp.predict(X_test)

The same code could be easily translated also in a pure C++ application as

.. code-block:: c++

   #include <rfbp.hpp>

   int main ()
   {
     const int N = 20;
     const int M = 101; // M must be odd
     const int K = 3;

     FocusingProtocol fp("pseudo_reinforcement", M);
     Patterns patterns(N, M);

     long int ** bin_weights = focusingBP < MagP64 >(K,          // hidden,
                                                     patterns,   // patterns,
                                                     1000,       // max_iters,
                                                     101,        // max_steps,
                                                     42,         // seed,
                                                     0.5,        // damping,
                                                     "accurate", // accuracy1,
                                                     "exact",    // accuracy2,
                                                     0.1,        // randfact,
                                                     fp,         // fp,
                                                     0.1,        // epsil,
                                                     1,          // nth,
                                                     "",         // outfile,
                                                     "",         // outmess,
                                                     "",         // inmess,
                                                     false       // binmess
                                                     );

     // It is clearly an overfitting! But it works as example
     long int ** y_pred = nonbayes_test(bin_weights, patterns, K);

     return 0;
   }


.. _BinaryCommitteeMachineFBP: https://github.com/carlobaldassi/BinaryCommitteeMachineFBP.jl
.. _`scikit-learn`: https://github.com/scikit-learn/scikit-learn
.. _`scikit-optimize`: https://github.com/scikit-optimize/scikit-optimize
.. _`Zecchina et al`: https://www.pnas.org/content/113/48/E7655
.. _scorer: https://github.com/Nico-Curti/scorer


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   theory
   installation
   cppAPI/modules
   pyAPI/modules
   references

References
----------

- D. Dall'Olio, N. Curti, G. Castellani, A. Bazzani, D. Remondini. "Classification of Genome Wide Association data by Belief Propagation Neural network", CCS Italy, 2019.
- C. Baldassi, C. Borgs, J. T. Chayes, A. Ingrosso, C. Lucibello, L. Saglietti, and R. Zecchina. "Unreasonable effectiveness of learning neural networks: From accessible states and robust ensembles to basic algorithmic schemes", Proceedings of the National Academy of Sciences, 113(48):E7655-E7662, 2016
- C. Baldassi, A. Braunstein, N. Brunel, R. Zecchina. "Efficient supervised learning in networks with binary synapses", Proceedings of the National Academy of Sciences, 104(26):11079-11084, 2007.
- A., Braunstein, R. Zecchina. "Learning by message passing in networks of discrete synapses". Physical Review Letters 96(3), 2006.
- C. Baldassi, F. Gerace, C. Lucibello, L. Saglietti, R. Zecchina. "Learning may need only a few bits of synaptic precision", Physical Review E, 93, 2016
- A. Blum, R. L. Rivest. "Training a 3-node neural network is NP-complete", Neural Networks, 1992
- W. Krauth, M. Mezard. "Storage capacity of memory networks with binary coupling", Journal of Physics (France), 1989
- H. Huang, Y. Kabashima. "Origin of the computational hardness for learning with binary synapses", Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 2014
- C. Baldassi, A. Ingrosso, C. Lucibello, L. Saglietti, R. Zecchina. "Local entropy as a measure for sampling solutions in constraint satisfaction problems", Journal of Statistical Mechanics: Theory and Experiment, 2016
- R. Monasson, R. Zecchina. "Learning and Generalization Theories of Large Committee Machines", Modern Physics Letters B, 1995
- R. Monasson, R. Zecchina. "Weight space structure and internal representations: A direct approach to learning and generalization in multilayer neural networks", Physical Review Letters, 1995
- C. Baldassi, A. Braunstein. "A Max-Sum algorithm for training discrete neural networks", Journal of Statistical Mechanics: Theory and Experiment, 2015
- G. Parisi. "Mean field theory of spin glasses: statics and dynamics", arXiv, 2007
- L. Dall'Asta, A. Ramezanpour, R. Zecchina. "Entropy landscape and non-Gibbs solutions in constraint satisfaction problem", Physical Review E, 2008
- M. Mézard, A. Montanari. "Information, Physics and Computation", Oxford Graduate Texts, 2009
- C. Baldassi, A. Ingrosso, C. Lucibello, L. Saglietti, R. Zecchina. "Subdominant Dense Clusters Allow for Simple Learning and High Computational Performance in Neural Networks with Discrete Synapses", Physical Review Letters, 2015.
MagT64
------------------

.. automodule:: rfbp.MagT64
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __mod__, __xor__
Magnetization functions
-----------------------

.. automodule:: rfbp.magnetization
   :members:
   :undoc-members:
   :show-inheritance:
MagP64
------------------

.. automodule:: rfbp.MagP64
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __mod__, __xor__
Atanherf
--------------------

.. automodule:: rfbp.atanherf
   :members:
   :undoc-members:
   :show-inheritance:
ReplicatedFocusingBeliefPropagation
-----------------------------------------------

.. automodule:: rfbp.ReplicatedFocusingBeliefPropagation
   :members:
   :undoc-members:
   :show-inheritance:

Patterns
--------------------

.. automodule:: rfbp.Patterns
   :members:
   :undoc-members:
   :show-inheritance:
Mag type
==========

.. toctree::
   :maxdepth: 4

   MagP64
   MagT64
   magnetization

Mag
---------------

.. automodule:: rfbp.Mag
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __add__, __truediv__, __mul__, __sub__, __neg__, __eq__, __ne__
Python API
==========

.. toctree::
   :maxdepth: 4

   FocusingProtocol
   Mag
   ReplicatedFocusingBeliefPropagation
   Patterns
   atanherf
FocusingProtocol
----------------------------

.. automodule:: rfbp.FocusingProtocol
   :members:
   :undoc-members:
   :show-inheritance:
Cavity Message
--------------

.. doxygenclass:: Cavity_Message
   :project: Cavity_Message
   :members:

Params
------

.. doxygenclass:: Params
   :project: Params
   :members:

MagT64
------

.. doxygenclass:: MagT64
   :project: MagT64
   :members:

Magnetization functions
-----------------------

.. doxygennamespace:: mag
   :project: magnetization

MagP64
------

.. doxygenclass:: MagP64
   :project: MagP64
   :members:

Atanherf
--------------------

.. doxygennamespace:: AtanhErf
   :project: AtanhErf

Replicated Focusing Belief Propagation
--------------------------------------

.. doxygenfile:: rfbp.h
   :project: ReplicatedFocusingBeliefPropagation

Patterns
---------

.. doxygenclass:: Patterns
   :project: Patterns
   :members:Mag type
==========

.. toctree::
   :maxdepth: 4

   MagP64
   MagT64
   magnetization
C++ API
-------

.. toctree::
   :maxdepth: 4

   FocusingProtocol
   Mag
   ReplicatedFocusingBeliefPropagation
   Patterns
   params
   cavity_message
   atanherf
FocusingProtocol
----------------

.. doxygenclass:: FocusingProtocol
   :project: FocusingProtocol
   :members: