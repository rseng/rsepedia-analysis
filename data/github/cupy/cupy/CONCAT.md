<div align="center"><img src="https://raw.githubusercontent.com/cupy/cupy/master/docs/image/cupy_logo_1000px.png" width="400"/></div>

# CuPy : NumPy & SciPy for GPU

[![pypi](https://img.shields.io/pypi/v/cupy.svg)](https://pypi.python.org/pypi/cupy)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/cupy.svg)](https://anaconda.org/conda-forge/cupy)
[![GitHub license](https://img.shields.io/github/license/cupy/cupy.svg)](https://github.com/cupy/cupy)
[![coveralls](https://img.shields.io/coveralls/cupy/cupy.svg)](https://coveralls.io/github/cupy/cupy)
[![Gitter](https://badges.gitter.im/cupy/community.svg)](https://gitter.im/cupy/community)
[![Twitter](https://img.shields.io/twitter/follow/CuPy_Team?label=%40CuPy_Team)](https://twitter.com/CuPy_Team)

[**Website**](https://cupy.dev/)
| [**Install**](https://docs.cupy.dev/en/stable/install.html)
| [**Tutorial**](https://docs.cupy.dev/en/stable/user_guide/basic.html)
| [**Examples**](https://github.com/cupy/cupy/tree/master/examples)
| [**Documentation**](https://docs.cupy.dev/en/stable/)
| [**API Reference**](https://docs.cupy.dev/en/stable/reference/)
| [**Forum**](https://groups.google.com/forum/#!forum/cupy)

CuPy is a NumPy/SciPy-compatible array library for GPU-accelerated computing with Python.
CuPy acts as a [drop-in replacement](https://docs.cupy.dev/en/stable/reference/comparison.html) to run existing NumPy/SciPy code on NVIDIA CUDA or AMD ROCm platforms.

```py
>>> import cupy as cp
>>> x = cp.arange(6).reshape(2, 3).astype('f')
>>> x
array([[ 0.,  1.,  2.],
       [ 3.,  4.,  5.]], dtype=float32)
>>> x.sum(axis=1)
array([  3.,  12.], dtype=float32)
```

CuPy also provides access to low-level CUDA features.
You can pass `ndarray` to existing CUDA C/C++ programs via [RawKernels](https://docs.cupy.dev/en/stable/user_guide/kernel.html#raw-kernels), use [Streams](https://docs.cupy.dev/en/stable/reference/cuda.html) for performance, or even call [CUDA Runtime APIs](https://docs.cupy.dev/en/stable/reference/cuda.html#runtime-api) directly.

## Installation

Wheels (precompiled binary packages) are available for Linux (x86_64) and Windows (amd64).
Choose the right package for your platform.

| Platform      | Command                       |
| ------------- | ----------------------------- |
| CUDA 10.2     | `pip install cupy-cuda102`    |
| CUDA 11.0     | `pip install cupy-cuda110`    |
| CUDA 11.1     | `pip install cupy-cuda111`    |
| CUDA 11.2     | `pip install cupy-cuda112`    |
| CUDA 11.3     | `pip install cupy-cuda113`    |
| CUDA 11.4     | `pip install cupy-cuda114`    |
| CUDA 11.5     | `pip install cupy-cuda115`    |
| ROCm 4.0 (*)  | `pip install cupy-rocm-4-0`   |
| ROCm 4.2 (*)  | `pip install cupy-rocm-4-2`   |
| ROCm 4.3 (*)  | `pip install cupy-rocm-4-3`   |

(\*) ROCm support is an experimental feature. Refer to the [docs](https://docs.cupy.dev/en/latest/install.html#using-cupy-on-amd-gpu-experimental) for details.

Use `-f https://pip.cupy.dev/pre` option to install pre-releases (e.g., `pip install cupy-cuda114 -f https://pip.cupy.dev/pre`).
See the [Installation Guide](https://docs.cupy.dev/en/stable/install.html) if you are using Conda/Anaconda or building from source.

## Run on Docker

Use [NVIDIA Container Toolkit](https://github.com/NVIDIA/nvidia-docker) to run CuPy image with GPU.

```
$ docker run --gpus all -it cupy/cupy
```

## More information

- [Release Notes](https://github.com/cupy/cupy/releases)
- [Projects using CuPy](https://github.com/cupy/cupy/wiki/Projects-using-CuPy)
- [Contribution Guide](https://docs.cupy.dev/en/stable/contribution.html)

## License

MIT License (see `LICENSE` file).

CuPy is designed based on NumPy's API and SciPy's API (see `docs/LICENSE_THIRD_PARTY` file).

CuPy is being maintained and developed by [Preferred Networks Inc.](https://preferred.jp/en/) and [community contributors](https://github.com/cupy/cupy/graphs/contributors).

## Reference

Ryosuke Okuta, Yuya Unno, Daisuke Nishino, Shohei Hido and Crissman Loomis.
**CuPy: A NumPy-Compatible Library for NVIDIA GPU Calculations.**
*Proceedings of Workshop on Machine Learning Systems (LearningSys) in The Thirty-first Annual Conference on Neural Information Processing Systems (NIPS)*, (2017).
[[PDF](http://learningsys.org/nips17/assets/papers/paper_16.pdf)]

```bibtex
@inproceedings{cupy_learningsys2017,
  author       = "Okuta, Ryosuke and Unno, Yuya and Nishino, Daisuke and Hido, Shohei and Loomis, Crissman",
  title        = "CuPy: A NumPy-Compatible Library for NVIDIA GPU Calculations",
  booktitle    = "Proceedings of Workshop on Machine Learning Systems (LearningSys) in The Thirty-first Annual Conference on Neural Information Processing Systems (NIPS)",
  year         = "2017",
  url          = "http://learningsys.org/nips17/assets/papers/paper_16.pdf"
}
```
# CuPy Code of Conduct

CuPy follows the [NumFOCUS Code of Conduct][homepage] available at https://numfocus.org/code-of-conduct.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at `dlfw@preferred.jp`.

[homepage]: https://numfocus.org/
# "include" directory

All files and directories in this directory will be copied to the distribution (sdist and wheel).
Note that items starting with `.` (e.g., `cub/.git`) are excluded.
See `setup.py` for details.

## CUB

The `cub` folder is a git submodule for the CUB project.
Including the CUB headers as a submodule enables not only building the `cupy.cuda.cub` module,
but also easier maintenance.
For further information on CUB, see the [CUB Project Website](http://nvlabs.github.com/cub).

## Jitify
The `Jitify` folder is a git submodule for the Jitify project.
Including the Jitify header as a submodule for building the `cupy.cuda.jitify` module.
For further information on Jitify, see the [Jitify repo](https://github.com/NVIDIA/jitify).

## DLPack
The `dlpack` folder stores the DLPack header for building the `cupy._core.dlpack` module,
see `README.md` therein.
For further information on DLPack, see the [DLPack repo](https://github.com/dmlc/dlpack).
## DLPack header

The header `dlpack.h` is downloaded from https://github.com/dmlc/dlpack/blob/main/include/dlpack/dlpack.h.
The commit is [`2775088`](https://github.com/dmlc/dlpack/commit/277508879878e0a5b5b43599b1bea11f66eb3c6c).
These files are copied from thrust project and are modified.

    http://thrust.github.io/These files are copied from CUDA Toolkit Distribution and redistributed under the following license:

https://docs.nvidia.com/cuda/archive/9.2/eula/#nvidia-cuda-toolkit-license-agreement

For CUDA 11.2+, we enable CUDA Enhanced Compatibility by hosting the fp16 headers from the latest
CUDA Toolkit release and place them in the ``cuda-11`` folder.

Currenly, the headers in ``cuda-11`` are from CTK 11.5.
Please refer to [our Contribution Guide](https://docs.cupy.dev/en/stable/contribution.html).

# CuPy ROCm Docker Image

```
docker run -it --device=/dev/kfd --device=/dev/dri --group-add video cupy/cupy-rocm
```
# SGEMM example

This example contains implementation of single-precision general matrix-multiplication  (SGEMM).
The implementation is based on the one in [MAGMA](http://icl.cs.utk.edu/magma/).


### How to demo
The demo contains a script that calculates matrix multiplication of A (m x k) and B (k x n).
The demo can be run by the following command.

```
python sgemm.py [--gpu GPU_ID] [--m m] [--n n] [--k k]
```


### What this demo contains

In this example, we work on a SGEMM kernel that requires a complete interface to `cuLaunchKernel` (e.g. grid size and size of shared memory), which is not provided by `cupy.ElementwiseKernel`.
CuPy arrays work regardless of the underlying memory layouts thanks to `ndarray` abstraction.
As is the case for this example, `ndarray` abstraction does not need to be used if the underlying memory layouts of arrays match the ones expected by a kernel.
The SGEMM kernel expects input and output arrays to be in Fortran contiguous memory layout, and this layout is enforced by `cupy.asfortranarray`.

#### How to dynamically compile and launch a kernel function written in CUDA C

For compilation, `cupy.RawKernel` class is used to compile a CUDA code written in `sgemm.cu`.
The class takes a text of code and name of the kernel as an constructor argument.
The instance is a callable; the CUDA code will be compiled and then invoked when it is called.
The compiled code is cached, and it avoids the compilation process after the first time.
Also, the CUDA code can be modified at Python level because it is simply a text.
In this example, C macros that determine a distribution of data to threads are specified at runtime.
Note that `"extern C"` needs to be put on top of the kernel that is called.

#### How to supply grid size, block size and shared memory size on launching a kernel function

`cupy.RawKernel` object allows you to call the kernel with CUDA's `cuLaunchKernel` interface.
In other words, you have control over grid size, block size, shared memory size and stream.
At this level of interface, it becomes straightforward to replace host `.cu` that calls CUDA kernels with Python code.
# kmeans example

This example contains implementation of K-means clustering.


### How to demo
The demo contains a script that partitions data into groups using K-means clustering.
The demo can be run by the following command.

```
python kmeans.py [--gpu-id GPU_ID] [--n-clusters N_CLUSTERS] [--num NUM]
                 [--max-iter MAX_ITER] [--use-custom-kernel]
                 [--output-image OUTPUT_IMAGE]
```

If you run this script on environment without matplotlib renderers (e.g., non-GUI environment), setting the environmental variable `MPLBACKEND` to `Agg` may be required to use `matplotlib`. For example,

```
MPLBACKEND=Agg python kmeans.py ...
```
# GMM example

This example contains implementation of Gaussian Mixture Model (GMM).


### How to demo
The demo contains a script that partitions data into groups using Gaussian Mixture Model.
The demo can be run by the following command.

```
python gmm.py [--gpu-id GPU_ID] [--num NUM] [--dim DIM]
              [--max-iter MAX_ITER] [--tol TOL] [--output-image OUTPUT]
```

If you run this script on environment without matplotlib renderers (e.g., non-GUI environment), setting the environmental variable `MPLBACKEND` to `Agg` may be required to use `matplotlib`. For example,

```
MPLBACKEND=Agg python gmm.py ...
```
# Custom user structure examples

This folder contains examples of custom user structures in `cupy.RawKernel` (see [https://docs.cupy.dev/en/stable/tutorial/kernel.html](https://docs.cupy.dev/en/stable/tutorial/kernel.html) for corresponding documentation).

This folder provides three scripts ranked by increasing complexity:

1. `builtins_vectors.py` shows how to use CUDA builtin vectors such as `float4` both as scalar parameter (pass by value from host) and array parameter in RawKernels.
2. `packed_matrix.py` demonstrates how to create and use templated packed structures in RawModules.
3. `complex_struct.py` illustrates the possibility to recursively build complex NumPy dtypes matching device structure memory layout.

All examples can be run as simple python scripts: `python3.x example_name.py`.
<!-- AUTO GENERATED: DO NOT EDIT! -->

# CuPy CI Test Coverage

| Param                 |              | Test                        |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | #   |
| --------------------- | ------------ | --------------------------- | --------------------------------- | --------------------------- | --------------------------------- | --------------------------- | --------------------------------- | --------------------------- | --------------------------------- | --------------------------- | --------------------------------- | ------------------------------ | ------------------------------------ | ------------------------------ | ------------------------------------ | ------------------------------- | ------------------------------- | ------------------------------- | -------------------------------- | ----------------------------------- | -------------------------------- | ------------------------------------------ | --- |
|                       | System       | linux                       | linux                             | linux                       | linux                             | linux                       | linux                             | linux                       | linux                             | linux                       | linux                             | linux                          | linux                                | linux                          | linux                                | linux                           | linux                           | linux                           | linux                            | linux                               | linux                            | linux                                      |     |
|                       | Target       | [cuda102][t0][🐳][d0][📜][s0] | [cuda102.multi][t1][🐳][d1][📜][s1] | [cuda110][t2][🐳][d2][📜][s2] | [cuda110.multi][t3][🐳][d3][📜][s3] | [cuda111][t4][🐳][d4][📜][s4] | [cuda111.multi][t5][🐳][d5][📜][s5] | [cuda112][t6][🐳][d6][📜][s6] | [cuda112.multi][t7][🐳][d7][📜][s7] | [cuda113][t8][🐳][d8][📜][s8] | [cuda113.multi][t9][🐳][d9][📜][s9] | [cuda114][t10][🐳][d10][📜][s10] | [cuda114.multi][t11][🐳][d11][📜][s11] | [cuda115][t12][🐳][d12][📜][s12] | [cuda115.multi][t13][🐳][d13][📜][s13] | [rocm-4-0][t14][🐳][d14][📜][s14] | [rocm-4-2][t15][🐳][d15][📜][s15] | [rocm-4-3][t16][🐳][d16][📜][s16] | [cuda-slow][t17][🐳][d17][📜][s17] | [cuda-example][t18][🐳][d18][📜][s18] | [cuda-head][t19][🐳][d19][📜][s19] | [cuda11x-cuda-python][t20][🐳][d20][📜][s20] |     |
|                       |              |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            |     |
| system                | linux        | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    | ✅                               | ✅                               | ✅                               | ✅                                | ✅                                   | ✅                                | ✅                                          | 21  |
|                       | windows      |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
| os                    | ubuntu:18.04 | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 4   |
|                       | ubuntu:20.04 |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    | ✅                               | ✅                               | ✅                               | ✅                                | ✅                                   | ✅                                | ✅                                          | 13  |
|                       | centos:7     |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | centos:8     |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | ws:2016      |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
| cuda                  | null         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     |                                  |                                            | 3   |
|                       | 10.2         | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 11.0         |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 11.1         |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 11.2         |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 11.3         |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 11.4         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    |                                |                                      |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                |                                            | 5   |
|                       | 11.5         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      | ✅                              | ✅                                    |                                 |                                 |                                 |                                  |                                     |                                  | ✅                                          | 3   |
| rocm                  | null         | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                | ✅                                          | 18  |
|                       | 4.0          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               |                                 |                                 |                                  |                                     |                                  |                                            | 1   |
|                       | 4.2          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 | ✅                               |                                 |                                  |                                     |                                  |                                            | 1   |
|                       | 4.3          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 | ✅                               |                                  |                                     |                                  |                                            | 1   |
| nccl                  | null         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     |                                  |                                            | 3   |
|                       | 2.8          | ✅                           | ✅                                 |                             |                                   | ✅                           | ✅                                 | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 6   |
|                       | 2.9          |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 4   |
|                       | 2.10         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                |                                            | 3   |
|                       | 2.11         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 |                                  |                                     |                                  | ✅                                          | 5   |
| cutensor              | null         | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     |                                  |                                            | 5   |
|                       | 1.4          |                             |                                   | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                | ✅                                          | 16  |
| cusparselt            | null         | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     |                                  |                                            | 11  |
|                       | 0.1.0        |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                | ✅                                          | 10  |
| cudnn                 | null         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     |                                  |                                            | 3   |
|                       | 7.6          | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 8.0          |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 8.1          |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 8.2          |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 4   |
|                       | 8.3          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 | ✅                                | ✅                                   | ✅                                | ✅                                          | 8   |
| python                | 3.7          | ✅                           | ✅                                 |                             |                                   | ✅                           | ✅                                 | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      | ✅                               |                                 |                                 |                                  |                                     |                                  |                                            | 7   |
|                       | 3.8          |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 | ✅                               |                                 |                                  | ✅                                   |                                  |                                            | 4   |
|                       | 3.9          |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 | ✅                               | ✅                                |                                     | ✅                                |                                            | 5   |
|                       | 3.10         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 |                                  |                                     |                                  | ✅                                          | 5   |
|                       | pre          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
| numpy                 | 1.18         | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                           | ✅                                 | ✅                           | ✅                                 |                                |                                      |                                |                                      | ✅                               |                                 |                                 |                                  |                                     |                                  |                                            | 7   |
|                       | 1.19         |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 1.20         |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 | ✅                               |                                 |                                  | ✅                                   |                                  |                                            | 4   |
|                       | 1.21         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    |                                |                                      |                                 |                                 | ✅                               | ✅                                |                                     |                                  | ✅                                          | 5   |
|                       | 1.22         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      | ✅                              | ✅                                    |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | pre          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     | ✅                                |                                            | 1   |
| scipy                 | null         |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 1.4          | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      | ✅                               |                                 |                                 |                                  |                                     |                                  |                                            | 3   |
|                       | 1.5          |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 2   |
|                       | 1.6          |                             |                                   |                             |                                   |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 | ✅                               |                                 |                                  |                                     |                                  |                                            | 3   |
|                       | 1.7          |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 | ✅                               | ✅                                | ✅                                   |                                  | ✅                                          | 10  |
|                       | pre          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     | ✅                                |                                            | 1   |
| optuna                | null         |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
|                       | 2            | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    | ✅                               | ✅                               | ✅                               | ✅                                | ✅                                   |                                  | ✅                                          | 20  |
|                       | pre          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     | ✅                                |                                            | 1   |
| cython                | 0.29         | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    | ✅                               | ✅                               | ✅                               | ✅                                | ✅                                   | ✅                                | ✅                                          | 21  |
|                       | pre          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
| cuda-python           | null         | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                           | ✅                                 | ✅                              | ✅                                    | ✅                              | ✅                                    | ✅                               | ✅                               | ✅                               | ✅                                | ✅                                   | ✅                                |                                            | 20  |
|                       | 11           |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  | ✅                                          | 1   |
| env:CUPY_ACCELERATORS | null         | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                                |                                      |                                |                                      | ✅                               | ✅                               | ✅                               |                                  | ✅                                   |                                  |                                            | 8   |
|                       | cub          |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
|                       | cutensor     |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 0 🚨 |
|                       | cub,cutensor |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                           | ✅                                 |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  |                                     | ✅                                |                                            | 5   |
|                       | cutensor,cub |                             |                                   |                             |                                   | ✅                           | ✅                                 |                             |                                   |                             |                                   | ✅                              | ✅                                    | ✅                              | ✅                                    |                                 |                                 |                                 | ✅                                |                                     |                                  | ✅                                          | 8   |
| test                  | unit         | ✅                           |                                   | ✅                           |                                   | ✅                           |                                   | ✅                           |                                   | ✅                           |                                   | ✅                              |                                      | ✅                              |                                      | ✅                               | ✅                               | ✅                               |                                  |                                     | ✅                                | ✅                                          | 12  |
|                       | unit-multi   |                             | ✅                                 |                             | ✅                                 |                             | ✅                                 |                             | ✅                                 |                             | ✅                                 |                                | ✅                                    |                                | ✅                                    |                                 |                                 |                                 |                                  |                                     |                                  |                                            | 7   |
|                       | unit-slow    |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 | ✅                                |                                     |                                  |                                            | 1   |
|                       | example      |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                             |                                   |                                |                                      |                                |                                      |                                 |                                 |                                 |                                  | ✅                                   |                                  |                                            | 1   |

[t0]:https://ci.preferred.jp/cupy.linux.cuda102/
[d0]:linux/tests/cuda102.Dockerfile
[s0]:linux/tests/cuda102.sh
[t1]:https://ci.preferred.jp/cupy.linux.cuda102.multi/
[d1]:linux/tests/cuda102.multi.Dockerfile
[s1]:linux/tests/cuda102.multi.sh
[t2]:https://ci.preferred.jp/cupy.linux.cuda110/
[d2]:linux/tests/cuda110.Dockerfile
[s2]:linux/tests/cuda110.sh
[t3]:https://ci.preferred.jp/cupy.linux.cuda110.multi/
[d3]:linux/tests/cuda110.multi.Dockerfile
[s3]:linux/tests/cuda110.multi.sh
[t4]:https://ci.preferred.jp/cupy.linux.cuda111/
[d4]:linux/tests/cuda111.Dockerfile
[s4]:linux/tests/cuda111.sh
[t5]:https://ci.preferred.jp/cupy.linux.cuda111.multi/
[d5]:linux/tests/cuda111.multi.Dockerfile
[s5]:linux/tests/cuda111.multi.sh
[t6]:https://ci.preferred.jp/cupy.linux.cuda112/
[d6]:linux/tests/cuda112.Dockerfile
[s6]:linux/tests/cuda112.sh
[t7]:https://ci.preferred.jp/cupy.linux.cuda112.multi/
[d7]:linux/tests/cuda112.multi.Dockerfile
[s7]:linux/tests/cuda112.multi.sh
[t8]:https://ci.preferred.jp/cupy.linux.cuda113/
[d8]:linux/tests/cuda113.Dockerfile
[s8]:linux/tests/cuda113.sh
[t9]:https://ci.preferred.jp/cupy.linux.cuda113.multi/
[d9]:linux/tests/cuda113.multi.Dockerfile
[s9]:linux/tests/cuda113.multi.sh
[t10]:https://ci.preferred.jp/cupy.linux.cuda114/
[d10]:linux/tests/cuda114.Dockerfile
[s10]:linux/tests/cuda114.sh
[t11]:https://ci.preferred.jp/cupy.linux.cuda114.multi/
[d11]:linux/tests/cuda114.multi.Dockerfile
[s11]:linux/tests/cuda114.multi.sh
[t12]:https://ci.preferred.jp/cupy.linux.cuda115/
[d12]:linux/tests/cuda115.Dockerfile
[s12]:linux/tests/cuda115.sh
[t13]:https://ci.preferred.jp/cupy.linux.cuda115.multi/
[d13]:linux/tests/cuda115.multi.Dockerfile
[s13]:linux/tests/cuda115.multi.sh
[t14]:https://jenkins.preferred.jp/job/chainer/job/cupy_master/TEST=rocm-4-0,label=mnj-mi50/
[d14]:linux/tests/rocm-4-0.Dockerfile
[s14]:linux/tests/rocm-4-0.sh
[t15]:https://jenkins.preferred.jp/job/chainer/job/cupy_master/TEST=rocm-4-2,label=mnj-mi50/
[d15]:linux/tests/rocm-4-2.Dockerfile
[s15]:linux/tests/rocm-4-2.sh
[t16]:https://jenkins.preferred.jp/job/chainer/job/cupy_master/TEST=rocm-4-3,label=mnj-mi50/
[d16]:linux/tests/rocm-4-3.Dockerfile
[s16]:linux/tests/rocm-4-3.sh
[t17]:https://ci.preferred.jp/cupy.linux.cuda-slow/
[d17]:linux/tests/cuda-slow.Dockerfile
[s17]:linux/tests/cuda-slow.sh
[t18]:https://ci.preferred.jp/cupy.linux.cuda-example/
[d18]:linux/tests/cuda-example.Dockerfile
[s18]:linux/tests/cuda-example.sh
[t19]:https://ci.preferred.jp/cupy.linux.cuda-head/
[d19]:linux/tests/cuda-head.Dockerfile
[s19]:linux/tests/cuda-head.sh
[t20]:https://ci.preferred.jp/cupy.linux.cuda11x-cuda-python/
[d20]:linux/tests/cuda11x-cuda-python.Dockerfile
[s20]:linux/tests/cuda11x-cuda-python.sh
# CuPy CI

CuPy uses two infrastructures for GPU tests.

* FlexCI (`pfn-public-ci`): GCP, Linux/Windows, CUDA only
* Jenkins: on-premise, Linux only, CUDA/ROCm

Currently most of test configurations are managed by [chainer-test](http://github.com/chainer/chainer-test), which contains a set of scripts for Jenkins, but we are gradually migrating to new tooling in this directory.
We are also gradually migrating from Jenkins to FlexCI for better performance; eventually Jenkins will only be used for ROCm tests.

This directory contains the test matrix definition, and a tool to generate test environment from the matrix.

* `schema.yaml` defines all the possible values for each test axis, and constraints between them.
* `matrix.yaml` defines the configuration of each matrix.
* `generate.py` generates the test environment (Dockerfile/shell script for Linux, PowerShell script for Windows) for each matrix from the schema and the matrix.
  This program also generates `coverage.md` to see the configuration coverage.

## Usage

To generate `linux/tests/*.Dockerfile`, `linux/tests/*.sh` and `coverage.md`:

```
pip install PyYAML
./generate.py -s schema.yaml -m matrix.yaml
```

## Future work

* Support generating Windows test environment.
* Test notifications to Gitter.
* Generate shuffle tests from `schema.yaml`.
* Support using OS-provided Python binary packages instead of pyenv.
* Support coverage reporting.
* Support installation tests.
# Linux CI Scripts

This directory contains assets used for CI.

All tests here are isolated by Docker so that developers and contributors can reproduce the same environment as the CI.
You can use the `run.sh` tool to build the Docker image and run unit tests in the image.
The current (local) codebase is read-only mounted to the container and used for testing.

`./run.sh` takes a TARGET and one or more STAGEs as arguments.
Run `./run.sh` without arguments for the detailed usage and list of all STAGEs.

Here are some examples:

```
# Target: cuda114
# Stages: Build the docker image for testing, then run the unit test.
./run.sh cuda114 build test

# Target: rocm-4-0
# Stages: Only build the docker image.
./run.sh rocm-4-0 build
```

## Test Targets

`tests/` directory contains Dockerfiles and bootstrap shell scripts with TARGET name prefixed.
For example, when the target is `rocm-4-2`, `rocm-4-2.Dockerfile` and `rocm-4-2.sh` are used for testing.

These files are generated by the `generate.py` tool, except for the following targets:

* `cuda-rapids`: Run unit tests for components requiring RAPIDS (e.g., `cupyx.scipy.sparse.csgraph`).
* `array-api`: Run [array-api-tests](https://github.com/data-apis/array-api-tests) for `cupy.array_api` module.
=============
Upgrade Guide
=============

This page covers changes introduced in each major version that users should know when migrating from older releases.
Please see also the :ref:`compatibility_matrix` for supported environments of each major version.

CuPy v10
========

Dropping CUDA 9.2 / 10.0 / 10.1 Support
---------------------------------------

CUDA 10.1 or earlier is no longer supported.
Use CUDA 10.2 or later.

Dropping NCCL v2.4 / v2.6 / v2.7 Support
----------------------------------------

NCCL v2.4, v2.6, and v2.7 are no longer supported.

Dropping Python 3.6 Support
---------------------------

Python 3.6 is no longer supported.

Dropping NumPy 1.17 Support
---------------------------

NumPy 1.17 is no longer supported.

Change in :class:`cupy.cuda.Device` Behavior
--------------------------------------------

Current device set via ``use()`` will not be restored when exiting ``with`` block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current device set via :func:`cupy.cuda.Device.use()` will not be reactivated when exiting a device context manager. An existing code mixing ``with device:`` block and ``device.use()`` may get different results between CuPy v10 and v9.

.. code-block:: py

   with cupy.cuda.Device(1) as d1:
       d2 = cupy.cuda.Device(0).use()
       with d1:
           pass
       cupy.cuda.Device()  # -> CuPy v10 returns device 1 instead of device 0

Changes in :class:`cupy.cuda.Stream` Behavior
---------------------------------------------

Stream is now managed per-device
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Previoulys, it was users' responsibility to keep the current stream consistent with the current CUDA device.
For example, the following code raises an error in CuPy v9 or earlier:

.. code-block:: py

   import cupy

   with cupy.cuda.Device(0):
       # Create a stream on device 0.
       s0 = cupy.cuda.Stream()

   with cupy.cuda.Device(1):
       with s0:
           # Try to use the stream on device 1
           cupy.arange(10)  # -> CUDA_ERROR_INVALID_HANDLE: invalid resource handle

CuPy v10 manages the current stream per-device, thus eliminating the need of switching the stream every time the active device is changed.
When using CuPy v10, the above example behaves differently because whenever a stream is created, it is automatically associated with the current device and will be ignored when switching devices. 
In early versions, trying to use `s0` in device 1 raises an error because `s0` is associated with device 0. However, in v10, this `s0` is ignored and the default stream for device 1 will be used instead.

Current stream set via ``use()`` will not be restored when exiting ``with`` block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Samely as the change of :class:`cupy.cuda.Device` above, the current stream set via :func:`cupy.cuda.Stream.use` will not be reactivated when exiting a stream context manager.
An existing code mixing ``with stream:`` block and ``stream.use()`` may get different results between CuPy v10 and v9.

.. code-block:: py

   s1 = cupy.cuda.Stream()
   s2 = cupy.cuda.Stream()
   s3 = cupy.cuda.Stream()
   with s1:
       s2.use()
       with s3:
           pass
       cupy.cuda.get_current_stream()  # -> CuPy v10 returns `s1` instead of `s2`.

Streams can now be shared between threads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same :class:`cupy.cuda.Stream` instance can now safely be shared between multiple threads.

To achieve this, CuPy v10 will not destroy the stream (``cudaStreamDestroy``) if the stream is the current stream of any thread.

Big-Endian Arrays Automatically Converted to Little-Endian
----------------------------------------------------------

:func:`cupy.array`, :func:`cupy.asarray` and its variants now always transfer the data to GPU in little-endian byte order.

Previously CuPy was copying the given :class:`numpy.ndarray` to GPU as-is, regardless of the endianness.
In CuPy v10, big-endian arrays are converted to little-endian before the transfer, which is the native byte order on GPUs.
This change eliminates the need to manually change the array endianness before creating the CuPy array.

Baseline API Update
-------------------

Baseline API has been bumped from NumPy 1.20 and SciPy 1.6 to NumPy 1.21 and SciPy 1.7.
CuPy v10 will follow the upstream products' specifications of these baseline versions.

API Changes
-----------

* Device synchronize detection APIs (:func:`cupyx.allow_synchronize` and :class:`cupyx.DeviceSynchronized`), introduced as an experimental feature in CuPy v8, have been marked as deprecated because it is impossible to detect synchronizations reliably.

* An *internal* API :func:`cupy.cuda.compile_with_cache` has been marked as deprecated as there are better alternatives (see :class:`~cupy.RawModule` added since CuPy v7 and :class:`~cupy.RawKernel` since v5). While it has a longstanding history, this API has never been meant to be public. We encourage downstream libraries and users to migrate to the aforementioned public APIs. See :doc:`./user_guide/kernel` for their tutorials.

* The DLPack routine :func:`cupy.fromDlpack` is deprecated in favor of :func:`cupy.from_dlpack`, which addresses potential data race issues.

* A new module :mod:`cupyx.profiler` is added to host all profiling related APIs in CuPy. Accordingly, the following APIs are relocated to this module as follows. The old routines are deprecated.

    * :func:`cupy.prof.TimeRangeDecorator` -> :func:`cupyx.profiler.time_range`
    * :func:`cupy.prof.time_range` -> :func:`cupyx.profiler.time_range`
    * :func:`cupy.cuda.profile` -> :func:`cupyx.profiler.profile`
    * :func:`cupyx.time.repeat` -> :func:`cupyx.profiler.benchmark`

* :func:`cupy.ndarray.__pos__` now returns a copy (samely as :func:`cupy.positive`) instead of returning ``self``.

Note that deprecated APIs may be removed in the future CuPy releases.

Update of Docker Images
-----------------------

CuPy official Docker images (see :doc:`install` for details) are now updated to use CUDA 11.4 and ROCm 4.3.

CuPy v9
=======

Dropping Support of CUDA 9.0
----------------------------

CUDA 9.0 is no longer supported.
Use CUDA 9.2 or later.

Dropping Support of cuDNN v7.5 and NCCL v2.3
--------------------------------------------

cuDNN v7.5 (or earlier) and NCCL v2.3 (or earlier) are no longer supported.

Dropping Support of NumPy 1.16 and SciPy 1.3
--------------------------------------------

NumPy 1.16 and SciPy 1.3 are no longer supported.

Dropping Support of Python 3.5
------------------------------

Python 3.5 is no longer supported in CuPy v9.

NCCL and cuDNN No Longer Included in Wheels
-------------------------------------------

NCCL and cuDNN shared libraires are no longer included in wheels (see `#4850 <https://github.com/cupy/cupy/issues/4850>`_ for discussions). 
You can manually install them after installing wheel if you don't have a previous installation; see :doc:`install` for details.

cuTENSOR Enabled in Wheels
--------------------------

cuTENSOR can now be used when installing CuPy via wheels.

``cupy.cuda.{nccl,cudnn}`` Modules Needs Explicit Import
--------------------------------------------------------

Previously ``cupy.cuda.nccl`` and ``cupy.cuda.cudnn`` modules were automatically imported.
Since CuPy v9, these modules need to be explicitly imported (i.e., ``import cupy.cuda.nccl`` / ``import cupy.cuda.cudnn``.)

Baseline API Update
-------------------

Baseline API has been bumped from NumPy 1.19 and SciPy 1.5 to NumPy 1.20 and SciPy 1.6.
CuPy v9 will follow the upstream products' specifications of these baseline versions.

Following NumPy 1.20, aliases for the Python scalar types (``cupy.bool``, ``cupy.int``, ``cupy.float``, and ``cupy.complex``) are now deprecated.
``cupy.bool_``, ``cupy.int_``, ``cupy.float_`` and ``cupy.complex_`` should be used instead when required.

Update of Docker Images
-----------------------

CuPy official Docker images (see :doc:`install` for details) are now updated to use CUDA 11.2 and Python 3.8.


CuPy v8
=======

Dropping Support of CUDA 8.0 and 9.1
------------------------------------

CUDA 8.0 and 9.1 are no longer supported.
Use CUDA 9.0, 9.2, 10.0, or later.

Dropping Support of NumPy 1.15 and SciPy 1.2
--------------------------------------------

NumPy 1.15 (or earlier) and SciPy 1.2 (or earlier) are no longer supported.

Update of Docker Images
-----------------------

* CuPy official Docker images (see :doc:`install` for details) are now updated to use CUDA 10.2 and Python 3.6.
* SciPy and Optuna are now pre-installed.

CUB Support and Compiler Requirement
------------------------------------

CUB module is now built by default.
You can enable the use of CUB by setting ``CUPY_ACCELERATORS="cub"`` (see :envvar:`CUPY_ACCELERATORS` for details).

Due to this change, g++-6 or later is required when building CuPy from the source.
See :doc:`install` for details.

The following environment variables are no longer effective:

* ``CUB_DISABLED``: Use :envvar:`CUPY_ACCELERATORS` as aforementioned.
* ``CUB_PATH``: No longer required as CuPy uses either the CUB source bundled with CUDA (only when using CUDA 11.0 or later) or the one in the CuPy distribution.

API Changes
-----------

* ``cupy.scatter_add``, which was deprecated in CuPy v4, has been removed. Use :func:`cupyx.scatter_add` instead.
* ``cupy.sparse`` module has been deprecated and will be removed in future releases. Use :mod:`cupyx.scipy.sparse` instead.
* ``dtype`` argument of :func:`cupy.ndarray.min` and :func:`cupy.ndarray.max` has been removed to align with the NumPy specification.
* :func:`cupy.allclose` now returns the result as 0-dim GPU array instead of Python bool to avoid device synchronization.
* :class:`cupy.RawModule` now delays the compilation to the time of the first call to align the behavior with :class:`cupy.RawKernel`.
* ``cupy.cuda.*_enabled`` flags (``nccl_enabled``, ``nvtx_enabled``, etc.) has been deprecated. Use ``cupy.cuda.*.available`` flag (``cupy.cuda.nccl.available``, ``cupy.cuda.nvtx.available``, etc.) instead.
* ``CHAINER_SEED`` environment variable is no longer effective. Use ``CUPY_SEED`` instead.


CuPy v7
=======

Dropping Support of Python 2.7 and 3.4
--------------------------------------

Starting from CuPy v7, Python 2.7 and 3.4 are no longer supported as it reaches its end-of-life (EOL) in January 2020 (2.7) and March 2019 (3.4).
Python 3.5.1 is the minimum Python version supported by CuPy v7.
Please upgrade the Python version if you are using affected versions of Python to any later versions listed under :doc:`install`.


CuPy v6
=======

Binary Packages Ignore ``LD_LIBRARY_PATH``
------------------------------------------

Prior to CuPy v6, ``LD_LIBRARY_PATH`` environment variable can be used to override cuDNN / NCCL libraries bundled in the binary distribution (also known as wheels).
In CuPy v6, ``LD_LIBRARY_PATH`` will be ignored during discovery of cuDNN / NCCL; CuPy binary distributions always use libraries that comes with the package to avoid errors caused by unexpected override.


CuPy v5
=======

``cupyx.scipy`` Namespace
-------------------------

:mod:`cupyx.scipy` namespace has been introduced to provide CUDA-enabled SciPy functions.
:mod:`cupy.sparse` module has been renamed to :mod:`cupyx.scipy.sparse`; :mod:`cupy.sparse` will be kept as an alias for backward compatibility.

Dropped Support for CUDA 7.0 / 7.5
----------------------------------

CuPy v5 no longer supports CUDA 7.0 / 7.5.

Update of Docker Images
-----------------------

CuPy official Docker images (see :doc:`install` for details) are now updated to use CUDA 9.2 and cuDNN 7.

To use these images, you may need to upgrade the NVIDIA driver on your host.
See `Requirements of nvidia-docker <https://github.com/NVIDIA/nvidia-docker/wiki/CUDA#requirements>`_ for details.


CuPy v4
=======

.. note::

   The version number has been bumped from v2 to v4 to align with the versioning of Chainer.
   Therefore, CuPy v3 does not exist.

Default Memory Pool
-------------------

Prior to CuPy v4, memory pool was only enabled by default when CuPy is used with Chainer.
In CuPy v4, memory pool is now enabled by default, even when you use CuPy without Chainer.
The memory pool significantly improves the performance by mitigating the overhead of memory allocation and CPU/GPU synchronization.

.. attention::

   When you monitor GPU memory usage (e.g., using ``nvidia-smi``), you may notice that GPU memory not being freed even after the array instance become out of scope.
   This is expected behavior, as the default memory pool "caches" the allocated memory blocks.

To access the default memory pool instance, use :func:`get_default_memory_pool` and :func:`get_default_pinned_memory_pool`.
You can access the statistics and free all unused memory blocks "cached" in the memory pool.

.. code-block:: py

   import cupy
   a = cupy.ndarray(100, dtype=cupy.float32)
   mempool = cupy.get_default_memory_pool()

   # For performance, the size of actual allocation may become larger than the requested array size.
   print(mempool.used_bytes())   # 512
   print(mempool.total_bytes())  # 512

   # Even if the array goes out of scope, its memory block is kept in the pool.
   a = None
   print(mempool.used_bytes())   # 0
   print(mempool.total_bytes())  # 512

   # You can clear the memory block by calling `free_all_blocks`.
   mempool.free_all_blocks()
   print(mempool.used_bytes())   # 0
   print(mempool.total_bytes())  # 0

You can even disable the default memory pool by the code below.
Be sure to do this before any other CuPy operations.

.. code-block:: py

   import cupy
   cupy.cuda.set_allocator(None)
   cupy.cuda.set_pinned_memory_allocator(None)

Compute Capability
------------------

CuPy v4 now requires NVIDIA GPU with Compute Capability 3.0 or larger.
See the `List of CUDA GPUs <https://developer.nvidia.com/cuda-gpus>`_ to check if your GPU supports Compute Capability 3.0.


CUDA Stream
-----------

As CUDA Stream is fully supported in CuPy v4, ``cupy.cuda.RandomState.set_stream``, the function to change the stream used by the random number generator, has been removed.
Please use :func:`cupy.cuda.Stream.use` instead.

See the discussion in `#306 <https://github.com/cupy/cupy/pull/306>`_ for more details.

``cupyx`` Namespace
-------------------

``cupyx`` namespace has been introduced to provide features specific to CuPy (i.e., features not provided in NumPy) while avoiding collision in future.
See :doc:`reference/ext` for the list of such functions.

For this rule, :func:`cupy.scatter_add` has been moved to :func:`cupyx.scatter_add`.
:func:`cupy.scatter_add` is still available as an alias, but it is encouraged to use :func:`cupyx.scatter_add` instead.

Update of Docker Images
-----------------------

CuPy official Docker images (see :doc:`install` for details) are now updated to use CUDA 8.0 and cuDNN 6.0.
This change was introduced because CUDA 7.5 does not support NVIDIA Pascal GPUs.

To use these images, you may need to upgrade the NVIDIA driver on your host.
See `Requirements of nvidia-docker <https://github.com/NVIDIA/nvidia-docker/wiki/CUDA#requirements>`_ for details.

CuPy v2
=======

Changed Behavior of count_nonzero Function
------------------------------------------

For performance reasons, :func:`cupy.count_nonzero` has been changed to return zero-dimensional :class:`ndarray` instead of `int` when `axis=None`.
See the discussion in `#154 <https://github.com/cupy/cupy/pull/154>`_ for more details.


.. _compatibility_matrix:

Compatibility Matrix
====================

.. list-table::
   :header-rows: 1

   * - CuPy
     - CC [1]_
     - CUDA
     - ROCm
     - cuTENSOR
     - NCCL
     - cuDNN
     - Python
     - NumPy
     - SciPy
     - Baseline API Spec.
     - Docs
   * - v11
     -
     -
     -
     -
     -
     -
     -
     -
     -
     -
     - `latest <https://docs.cupy.dev/en/stable/install.html>`__
   * - v10
     - 3.0~
     - 10.2~
     - 4.0~
     - 1.3~
     - 2.8~
     - 7.6~
     - 3.7~
     - 1.18~
     - 1.4~
     - NumPy 1.21 & SciPy 1.7
     - `stable <https://docs.cupy.dev/en/stable/install.html>`__
   * - v9
     - 3.0~8.x
     - 9.2~11.5
     - 3.5~4.3
     - 1.2~1.3
     - 2.4 & 2.6~2.11
     - 7.6~8.2
     - 3.6~3.9
     - 1.17~1.21
     - 1.4~1.7
     - NumPy 1.20 & SciPy 1.6
     - `v9.6.0 <https://docs.cupy.dev/en/v9.6.0/install.html>`__
   * - v8
     - 3.0~8.x
     - 9.0 & 9.2~11.2
     - 3.x [2]_
     - 1.2
     - 2.0~2.8
     - 7.0~8.1
     - 3.5~3.9
     - 1.16~1.20
     - 1.3~1.6
     - NumPy 1.19 & SciPy 1.5
     - `v8.6.0 <https://docs.cupy.dev/en/v8.6.0/install.html>`__
   * - v7
     - 3.0~8.x
     - 8.0~11.0
     - 2.x [2]_
     - 1.0
     - 1.3~2.7
     - 5.0~8.0
     - 3.5~3.8
     - 1.9~1.19
     - (not specified)
     - (not specified)
     - `v7.8.0 <https://docs.cupy.dev/en/v7.8.0/install.html>`__
   * - v6
     - 3.0~7.x
     - 8.0~10.1
     - n/a
     - n/a
     - 1.3~2.4
     - 5.0~7.5
     - 2.7 & 3.4~3.8
     - 1.9~1.17
     - (not specified)
     - (not specified)
     - `v6.7.0 <https://docs.cupy.dev/en/v6.7.0/install.html>`__
   * - v5
     - 3.0~7.x
     - 8.0~10.1
     - n/a
     - n/a
     - 1.3~2.4
     - 5.0~7.5
     - 2.7 & 3.4~3.7
     - 1.9~1.16
     - (not specified)
     - (not specified)
     - `v5.4.0 <https://docs.cupy.dev/en/v5.4.0/install.html>`__
   * - v4
     - 3.0~7.x
     - 7.0~9.2
     - n/a
     - n/a
     - 1.3~2.2
     - 4.0~7.1
     - 2.7 & 3.4~3.6
     - 1.9~1.14
     - (not specified)
     - (not specified)
     - `v4.5.0 <https://docs.cupy.dev/en/v4.5.0/install.html>`__

.. [1] CUDA Compute Capability
.. [2] Highly experimental support with limited features.
.. _overview:

Overview
========

`CuPy <https://github.com/cupy/cupy>`__ is a NumPy/SciPy-compatible array library for GPU-accelerated computing with Python.
CuPy acts as a drop-in replacement to run existing NumPy/SciPy code on `NVIDIA CUDA <https://developer.nvidia.com/cuda-toolkit>`__ or `AMD ROCm <https://www.amd.com/en/graphics/servers-solutions-rocm>`__ platforms.

CuPy provides a ``ndarray``, sparse matrices, and the associated routines for GPU devices, all having the same API as NumPy and SciPy:

* **N-dimensional array** (``ndarray``): :doc:`cupy.ndarray <reference/ndarray>`

  * Data types (dtypes): boolean (``bool_``), integer (``int8``, ``int16``, ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64``), float (``float16``, ``float32``, ``float64``), and complex (``complex64``, ``complex128``)
  * Supports the semantics identical to :class:`numpy.ndarray`, including basic / advanced indexing and broadcasting

* **Sparse matrices**: :doc:`cupyx.scipy.sparse <reference/scipy_sparse>`

  * 2-D sparse matrix: ``csr_matrix``, ``coo_matrix``, ``csc_matrix``, and ``dia_matrix``

* **NumPy Routines**

  * :doc:`Module-level Functions <reference/routines>` (``cupy.*``)
  * :doc:`Linear Algebra Functions <reference/linalg>` (``cupy.linalg.*``)
  * :doc:`Fast Fourier Transform <reference/fft>` (``cupy.fft.*``)
  * :doc:`Random Number Generator <reference/random>` (``cupy.random.*``)

* **SciPy Routines**

  * :doc:`Discrete Fourier Transforms <reference/scipy_fft>` (``cupyx.scipy.fft.*`` and ``cupyx.scipy.fftpack.*``)
  * :doc:`Advanced Linear Algebra <reference/scipy_linalg>` (``cupyx.scipy.linalg.*``)
  * :doc:`Multidimensional Image Processing <reference/scipy_ndimage>` (``cupyx.scipy.ndimage.*``)
  * :doc:`Sparse Matrices <reference/scipy_sparse>` (``cupyx.scipy.sparse.*``)
  * :doc:`Sparse Linear Algebra <reference/scipy_sparse_linalg>` (``cupyx.scipy.sparse.linalg.*``)
  * :doc:`Special Functions <reference/scipy_special>` (``cupyx.scipy.special.*``)
  * :doc:`Signal Processing <reference/scipy_signal>` (``cupyx.scipy.signal.*``)
  * :doc:`Statistical Functions <reference/scipy_stats>` (``cupyx.scipy.stats.*``)

Routines are backed by CUDA libraries (cuBLAS, cuFFT, cuSPARSE, cuSOLVER, cuRAND), Thrust, CUB, and cuTENSOR to provide the best performance.

It is also possible to easily implement :doc:`custom CUDA kernels <user_guide/kernel>` that work with ``ndarray`` using:

* **Kernel Templates**: Quickly define element-wise and reduction operation as a single CUDA kernel
* **Raw Kernel**: Import existing CUDA C/C++ code
* **Just-in-time Transpiler (JIT)**: Generate CUDA kernel from Python source code
* **Kernel Fusion**: Fuse multiple CuPy operations into a single CUDA kernel

CuPy can run in multi-GPU or cluster environments. The distributed communication package (:mod:`cupyx.distributed`) provides collective and peer-to-peer primitives for ``ndarray``, backed by NCCL.

For users who need more fine-grain control for performance, accessing :doc:`low-level CUDA features <user_guide/cuda_api>` are available:

* **Stream and Event**: CUDA stream and per-thread default stream are supported by all APIs
* **Memory Pool**: Customizable memory allocator with a built-in memory pool
* **Profiler**: Supports profiling code using CUDA Profiler and NVTX
* **Host API Binding**: Directly call CUDA libraries, such as NCCL, cuDNN, cuTENSOR, and cuSPARSELt APIs from Python

CuPy implements standard APIs for data exchange and interoperability, such as `DLPack <https://github.com/dmlc/dlpack>`__, `CUDA Array Interface <https://numba.readthedocs.io/en/stable/cuda/cuda_array_interface.html>`__, ``__array_ufunc__`` (`NEP 13 <https://numpy.org/neps/nep-0013-ufunc-overrides.html>`__), ``__array_function__`` (`NEP 18 <https://numpy.org/neps/nep-0018-array-function-protocol.html>`__), and `Array API Standard <https://data-apis.org/array-api/latest/>`__.
Thanks to these protocols, CuPy easily :doc:`integrates <user_guide/interoperability>` with NumPy, PyTorch, TensorFlow, MPI4Py, and any other libraries supporting the standard.

Under AMD ROCm environment, CuPy automatically translates all CUDA API calls to ROCm HIP (hipBLAS, hipFFT, hipSPARSE, hipRAND, hipCUB, hipThrust, RCCL, etc.), allowing code written using CuPy to run on both NVIDIA and AMD GPU without any modification.

Project Goal
------------

The goal of the CuPy project is to provide Python users GPU acceleration capabilities, without the in-depth knowledge of underlying GPU technologies.
The CuPy team focuses on providing:

* A complete NumPy and SciPy API coverage to become a full drop-in replacement, as well as advanced CUDA features to maximize the performance.
* Mature and quality library as a fundamental package for all projects needing acceleration, from a lab environment to a large-scale cluster.
============================================================
CuPy -- NumPy & SciPy for GPU
============================================================

.. module:: cupy

.. toctree::
   :maxdepth: 2

   overview
   install
   user_guide/index
   reference/index

.. toctree::
   :maxdepth: 2
   :caption: Development

   contribution

.. toctree::
   :maxdepth: 2
   :caption: Misc Notes

   upgrade
   license
Installation
============

Requirements
------------

* `NVIDIA CUDA GPU <https://developer.nvidia.com/cuda-gpus>`_ with the Compute Capability 3.0 or larger.

* `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`_: v10.2 / v11.0 / v11.1 / v11.2 / v11.3 / v11.4 / v11.5

    * If you have multiple versions of CUDA Toolkit installed, CuPy will automatically choose one of the CUDA installations.
      See :ref:`install_cuda` for details.

    * This requirement is optional if you install CuPy from ``conda-forge``. However, you still need to have a compatible
      driver installed for your GPU. See :ref:`install_cupy_from_conda_forge` for details.

* `Python <https://python.org/>`_: v3.7.0+ / v3.8.0+ / v3.9.0+ / v3.10.0+

.. note::

   Currently, CuPy is tested against  `Ubuntu <https://www.ubuntu.com/>`_ 18.04 LTS / 20.04 LTS (x86_64), `CentOS <https://www.centos.org/>`_ 7 / 8 (x86_64) and Windows Server 2016 (x86_64).

Python Dependencies
~~~~~~~~~~~~~~~~~~~

NumPy/SciPy-compatible API in CuPy v10 is based on NumPy 1.22 and SciPy 1.7, and has been tested against the following versions:

* `NumPy <https://numpy.org/>`_: v1.18 / v1.19 / v1.20 / v1.21 / v1.22

* `SciPy <https://scipy.org/>`_ (*optional*): v1.4 / v1.5 / v1.6 / v1.7

    * Required only when using :doc:`../reference/scipy` (``cupyx.scipy``).

* `Optuna <https://optuna.org/>`_ (*optional*): v2.x

    * Required only when using :ref:`kernel_param_opt`.

.. note::

   SciPy and Optuna are optional dependencies and will not be installed automatically.

.. note::

   Before installing CuPy, we recommend you to upgrade ``setuptools`` and ``pip``::

    $ python -m pip install -U setuptools pip

Additional CUDA Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~

Part of the CUDA features in CuPy will be activated only when the corresponding libraries are installed.

* `cuTENSOR <https://developer.nvidia.com/cutensor>`_: v1.4

    * The library to accelerate tensor operations. See :doc:`../reference/environment` for the details.

* `NCCL <https://developer.nvidia.com/nccl>`_: v2.8 / v2.9 / v2.10 / v2.11

    * The library to perform collective multi-GPU / multi-node computations.

* `cuDNN <https://developer.nvidia.com/cudnn>`_: v7.6 / v8.0 / v8.1 / v8.2 / v8.3

    * The library to accelerate deep neural network computations.

* `cuSPARSELt <https://docs.nvidia.com/cuda/cusparselt/>`_: v0.1.0

    * The library to accelerate sparse matrix-matrix multiplication.


Installing CuPy
---------------

Installing CuPy from PyPI
~~~~~~~~~~~~~~~~~~~~~~~~~

Wheels (precompiled binary packages) are available for Linux (x86_64) and Windows (amd64).
Package names are different depending on your CUDA Toolkit version.

.. list-table::
   :header-rows: 1

   * - CUDA
     - Command
   * - v10.2
     - ``$ pip install cupy-cuda102``
   * - v11.0
     - ``$ pip install cupy-cuda110``
   * - v11.1
     - ``$ pip install cupy-cuda111``
   * - v11.2
     - ``$ pip install cupy-cuda112``
   * - v11.3
     - ``$ pip install cupy-cuda113``
   * - v11.4
     - ``$ pip install cupy-cuda114``
   * - v11.5
     - ``$ pip install cupy-cuda115``

.. note::

   To enable features provided by additional CUDA libraries (cuTENSOR / NCCL / cuDNN), you need to install them manually.
   If you installed CuPy via wheels, you can use the installer command below to setup these libraries in case you don't have a previous installation::

    $ python -m cupyx.tools.install_library --cuda 11.2 --library cutensor

.. note::

   Use ``pip install cupy-cudaXXX -f https://pip.cupy.dev/pre`` to install pre-release (development) versions.


When using wheels, please be careful not to install multiple CuPy packages at the same time.
Any of these packages and ``cupy`` package (source installation) conflict with each other.
Please make sure that only one CuPy package (``cupy`` or ``cupy-cudaXX`` where XX is a CUDA version) is installed::

  $ pip freeze | grep cupy


.. _install_cupy_from_conda_forge:

Installing CuPy from Conda-Forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Conda/Anaconda is a cross-platform package management solution widely used in scientific computing and other fields.
The above ``pip install`` instruction is compatible with ``conda`` environments. Alternatively, for both Linux (x86_64,
ppc64le, aarch64-sbsa) and
Windows once the CUDA driver is correctly set up, you can also install CuPy from the ``conda-forge`` channel::

    $ conda install -c conda-forge cupy

and ``conda`` will install a pre-built CuPy binary package for you, along with the CUDA runtime libraries
(``cudatoolkit``). It is not necessary to install CUDA Toolkit in advance.

Conda has a built-in mechanism to determine and install the latest version of ``cudatoolkit`` supported by your driver.
However, if for any reason you need to force-install a particular CUDA version (say 11.0), you can do::

    $ conda install -c conda-forge cupy cudatoolkit=11.0

.. note::

    cuDNN, cuTENSOR, and NCCL are available on ``conda-forge`` as optional dependencies. The following command can install them all at once::

        $ conda install -c conda-forge cupy cudnn cutensor nccl

    Each of them can also be installed separately as needed.

.. note::

    If you encounter any problem with CuPy installed from ``conda-forge``, please feel free to report to `cupy-feedstock
    <https://github.com/conda-forge/cupy-feedstock/issues>`_, and we will help investigate if it is just a packaging
    issue in ``conda-forge``'s recipe or a real issue in CuPy.

.. note::

    If you did not install CUDA Toolkit by yourself, the ``nvcc`` compiler might not be available, as
    the ``cudatoolkit`` package from ``conda-forge`` does not include the ``nvcc`` compiler toolchain. If you would like to use
    it from a local CUDA installation, you need to make sure the version of CUDA Toolkit matches that of ``cudatoolkit`` to
    avoid surprises.


.. _install_cupy_from_source:

Installing CuPy from Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use of wheel packages is recommended whenever possible.
However, if wheels cannot meet your requirements (e.g., you are running non-Linux environment or want to use a version of CUDA / cuDNN / NCCL not supported by wheels), you can also build CuPy from source.

.. note::

   CuPy source build requires ``g++-6`` or later.
   For Ubuntu 18.04, run ``apt-get install g++``.
   For Ubuntu 16.04, CentOS 6 or 7, follow the instructions :ref:`here <install_gcc6>`.

.. note::

   When installing CuPy from source, features provided by additional CUDA libraries will be disabled if these libraries are not available at the build time.
   See :ref:`install_cudnn` for the instructions.

.. note::

   If you upgrade or downgrade the version of CUDA Toolkit, cuDNN, NCCL or cuTENSOR, you may need to reinstall CuPy.
   See :ref:`install_reinstall` for details.

You can install the latest stable release version of the `CuPy source package <https://pypi.python.org/pypi/cupy>`_ via ``pip``.

::

  $ pip install cupy

If you want to install the latest development version of CuPy from a cloned Git repository::

  $ git clone --recursive https://github.com/cupy/cupy.git
  $ cd cupy
  $ pip install .

.. note::

   Cython 0.29.22 or later is required to build CuPy from source.
   It will be automatically installed during the build process if not available.


Uninstalling CuPy
-----------------

Use ``pip`` to uninstall CuPy::

  $ pip uninstall cupy

.. note::

   If you are using a wheel, ``cupy`` shall be replaced with ``cupy-cudaXX`` (where XX is a CUDA version number).

.. note::

   If CuPy is installed via ``conda``, please do ``conda uninstall cupy`` instead.


Upgrading CuPy
---------------

Just use ``pip install`` with ``-U`` option::

  $ pip install -U cupy

.. note::

   If you are using a wheel, ``cupy`` shall be replaced with ``cupy-cudaXX`` (where XX is a CUDA version number).


.. _install_reinstall:


Reinstalling CuPy
-----------------

To reinstall CuPy, please uninstall CuPy and then install it.
When reinstalling CuPy, we recommend using ``--no-cache-dir`` option as ``pip`` caches the previously built binaries::

  $ pip uninstall cupy
  $ pip install cupy --no-cache-dir

.. note::

   If you are using a wheel, ``cupy`` shall be replaced with ``cupy-cudaXX`` (where XX is a CUDA version number).


Using CuPy inside Docker
------------------------

We are providing the `official Docker images <https://hub.docker.com/r/cupy/cupy/>`_.
Use `NVIDIA Container Toolkit <https://github.com/NVIDIA/nvidia-docker>`_ to run CuPy image with GPU.
You can login to the environment with bash, and run the Python interpreter::

  $ docker run --gpus all -it cupy/cupy /bin/bash

Or run the interpreter directly::

  $ docker run --gpus all -it cupy/cupy /usr/bin/python3


FAQ
---

.. _install_error:

``pip`` fails to install CuPy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please make sure that you are using the latest ``setuptools`` and ``pip``::

  $ pip install -U setuptools pip

Use ``-vvvv`` option with ``pip`` command.
This will display all logs of installation::

  $ pip install cupy -vvvv

If you are using ``sudo`` to install CuPy, note that ``sudo`` command does not propagate environment variables.
If you need to pass environment variable (e.g., ``CUDA_PATH``), you need to specify them inside ``sudo`` like this::

  $ sudo CUDA_PATH=/opt/nvidia/cuda pip install cupy

If you are using certain versions of conda, it may fail to build CuPy with error ``g++: error: unrecognized command line option ‘-R’``.
This is due to a bug in conda (see `conda/conda#6030 <https://github.com/conda/conda/issues/6030>`_ for details).
If you encounter this problem, please upgrade your conda.

.. _install_cudnn:

Installing cuDNN and NCCL
~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend installing cuDNN and NCCL using binary packages (i.e., using ``apt`` or ``yum``) provided by NVIDIA.

If you want to install tar-gz version of cuDNN and NCCL, we recommend installing it under the ``CUDA_PATH`` directory.
For example, if you are using Ubuntu, copy ``*.h`` files to ``include`` directory and ``*.so*`` files to ``lib64`` directory::

  $ cp /path/to/cudnn.h $CUDA_PATH/include
  $ cp /path/to/libcudnn.so* $CUDA_PATH/lib64

The destination directories depend on your environment.

If you want to use cuDNN or NCCL installed in another directory, please use ``CFLAGS``, ``LDFLAGS`` and ``LD_LIBRARY_PATH`` environment variables before installing CuPy::

  $ export CFLAGS=-I/path/to/cudnn/include
  $ export LDFLAGS=-L/path/to/cudnn/lib
  $ export LD_LIBRARY_PATH=/path/to/cudnn/lib:$LD_LIBRARY_PATH

.. _install_cuda:

Working with Custom CUDA Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have installed CUDA on the non-default directory or multiple CUDA versions on the same host, you may need to manually specify the CUDA installation directory to be used by CuPy.

CuPy uses the first CUDA installation directory found by the following order.

#. ``CUDA_PATH`` environment variable.
#. The parent directory of ``nvcc`` command. CuPy looks for ``nvcc`` command from ``PATH`` environment variable.
#. ``/usr/local/cuda``

For example, you can build CuPy using non-default CUDA directory by ``CUDA_PATH`` environment variable::

  $ CUDA_PATH=/opt/nvidia/cuda pip install cupy

.. note::

   CUDA installation discovery is also performed at runtime using the rule above.
   Depending on your system configuration, you may also need to set ``LD_LIBRARY_PATH`` environment variable to ``$CUDA_PATH/lib64`` at runtime.

CuPy always raises ``cupy.cuda.compiler.CompileException``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If CuPy raises a ``CompileException`` for almost everything, it is possible that CuPy cannot detect CUDA installed on your system correctly.
The followings are error messages commonly observed in such cases.

* ``nvrtc: error: failed to load builtins``
* ``catastrophic error: cannot open source file "cuda_fp16.h"``
* ``error: cannot overload functions distinguished by return type alone``
* ``error: identifier "__half_raw" is undefined``

Please try setting ``LD_LIBRARY_PATH`` and ``CUDA_PATH`` environment variable.
For example, if you have CUDA installed at ``/usr/local/cuda-9.2``::

  $ export CUDA_PATH=/usr/local/cuda-9.2
  $ export LD_LIBRARY_PATH=$CUDA_PATH/lib64:$LD_LIBRARY_PATH

Also see :ref:`install_cuda`.

.. _install_gcc6:

Build fails on Ubuntu 16.04, CentOS 6 or 7
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to build CuPy from source on systems with legacy GCC (g++-5 or earlier), you need to manually set up g++-6 or later and configure ``NVCC`` environment variable.

On Ubuntu 16.04::

  $ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  $ sudo apt update
  $ sudo apt install g++-6
  $ export NVCC="nvcc --compiler-bindir gcc-6"

On CentOS 6 / 7::

  $ sudo yum install centos-release-scl
  $ sudo yum install devtoolset-7-gcc-c++
  $ source /opt/rh/devtoolset-7/enable
  $ export NVCC="nvcc --compiler-bindir gcc"


Using CuPy on AMD GPU (experimental)
====================================

CuPy has an experimental support for AMD GPU (ROCm).

Requirements
------------

* `AMD GPU supported by ROCm <https://github.com/RadeonOpenCompute/ROCm#Hardware-and-Software-Support>`_

* `ROCm <https://rocmdocs.amd.com/en/latest/index.html>`_: v4.0 / v4.2 / v4.3
    * See the `ROCm Installation Guide <https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html>`_ for details.

The following ROCm libraries are required:

::

  $ sudo apt install hipblas hipsparse rocsparse rocrand rocthrust rocsolver rocfft hipcub rocprim rccl

Environment Variables
---------------------

When building or running CuPy for ROCm, the following environment variables are effective.

* ``ROCM_HOME``: directory containing the ROCm software (e.g., ``/opt/rocm``).

Docker
------

You can try running CuPy for ROCm using Docker.

::

  $ docker run -it --device=/dev/kfd --device=/dev/dri --group-add video cupy/cupy-rocm

.. _install_hip:

Installing Binary Packages
--------------------------

Wheels (precompiled binary packages) are available for Linux (x86_64).
Package names are different depending on your ROCm version.

.. list-table::
   :header-rows: 1

   * - ROCm
     - Command
   * - v4.0
     - ``$ pip install cupy-rocm-4-0``
   * - v4.2
     - ``$ pip install cupy-rocm-4-2``
   * - v4.3
     - ``$ pip install cupy-rocm-4-3``

Building CuPy for ROCm From Source
----------------------------------

To build CuPy from source, set the ``CUPY_INSTALL_USE_HIP``, ``ROCM_HOME``, and ``HCC_AMDGPU_TARGET`` environment variables.
(``HCC_AMDGPU_TARGET`` is the ISA name supported by your GPU.
Run ``rocminfo`` and use the value displayed in ``Name:`` line (e.g., ``gfx900``).
You can specify a comma-separated list of ISAs if you have multiple GPUs of different architectures.)

::

  $ export CUPY_INSTALL_USE_HIP=1
  $ export ROCM_HOME=/opt/rocm
  $ export HCC_AMDGPU_TARGET=gfx906
  $ pip install cupy

.. note::

  If you don't specify the ``HCC_AMDGPU_TARGET`` environment variable, CuPy will be built for the GPU architectures available on the build host.
  This behavior is specific to ROCm builds; when building CuPy for NVIDIA CUDA, the build result is not affected by the host configuration.

Limitations
-----------

The following features are not available due to the limitation of ROCm or because that they are specific to CUDA:

* CUDA Array Interface
* cuTENSOR
* Handling extremely large arrays whose size is around 32-bit boundary (HIP is known to fail with sizes `2**32-1024`)
* Atomic addition in FP16 (``cupy.ndarray.scatter_add`` and ``cupyx.scatter_add``)
* Multi-GPU FFT and FFT callback
* Some random number generation algorithms
* Several options in RawKernel/RawModule APIs: Jitify, dynamic parallelism
* Per-thread default stream
* Random generation API (``cupy.random.Generator``) for ROCm versions older than 4.3

The following features are not yet supported:

* Sparse matrices (``cupyx.scipy.sparse``)
* cuDNN (hipDNN)
* Hermitian/symmetric eigenvalue solver (``cupy.linalg.eigh``)
* Polynomial roots (uses Hermitian/symmetric eigenvalue solver)

The following features may not work in edge cases (e.g., some combinations of dtype):

.. note::
   We are investigating the root causes of the issues. They are not necessarily
   CuPy's issues, but ROCm may have some potential bugs.

* ``cupy.ndarray.__getitem__`` (`#4653 <https://github.com/cupy/cupy/pull/4653>`_)
* ``cupy.ix_`` (`#4654 <https://github.com/cupy/cupy/pull/4654>`_)
* Some polynomial routines (`#4758 <https://github.com/cupy/cupy/pull/4758>`_, `#4759 <https://github.com/cupy/cupy/pull/4759>`_)
* ``cupy.broadcast`` (`#4662 <https://github.com/cupy/cupy/pull/4662>`_)
* ``cupy.convolve`` (`#4668 <https://github.com/cupy/cupy/pull/4668>`_)
* ``cupy.correlate`` (`#4781 <https://github.com/cupy/cupy/pull/4781>`_)
* Some random sampling routines (``cupy.random``, `#4770 <https://github.com/cupy/cupy/pull/4770>`_)
* ``cupy.linalg.einsum``
* ``cupyx.scipy.ndimage`` and ``cupyx.scipy.signal`` (`#4878 <https://github.com/cupy/cupy/pull/4878>`_, `#4879 <https://github.com/cupy/cupy/pull/4879>`_, `#4880 <https://github.com/cupy/cupy/pull/4880>`_)
.. _contrib:

Contribution Guide
==================

This is a guide for all contributions to CuPy.
The development of CuPy is running on `the official repository at GitHub <https://github.com/cupy/cupy>`_.
Anyone that wants to register an issue or to send a pull request should read through this document.


Classification of Contributions
-------------------------------

There are several ways to contribute to CuPy community:

1. Registering an issue
2. Sending a pull request (PR)
3. Sending a question to `CuPy's Gitter channel <https://gitter.im/cupy/community>`_, `CuPy User Group <https://groups.google.com/forum/#!forum/cupy>`_, or `StackOverflow <https://stackoverflow.com/questions/tagged/cupy>`_
4. Open-sourcing an external example
5. Writing a post about CuPy

This document mainly focuses on 1 and 2, though other contributions are also appreciated.


Development Cycle
-----------------

This section explains the development process of CuPy.
Before contributing to CuPy, it is strongly recommended to understand the development cycle.

Versioning
~~~~~~~~~~

The versioning of CuPy follows `PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_ and a part of `Semantic versioning <https://semver.org/>`_.
The version number consists of three or four parts: ``X.Y.Zw`` where ``X`` denotes the **major version**, ``Y`` denotes the **minor version**, ``Z`` denotes the **revision number**, and the optional ``w`` denotes the prelease suffix.
While the major, minor, and revision numbers follow the rule of semantic versioning, the pre-release suffix follows PEP 440 so that the version string is much friendly with Python eco-system.

**Note that a major update basically does not contain compatibility-breaking changes from the last release candidate (RC).**
This is not a strict rule, though; if there is a critical API bug that we have to fix for the major version, we may add breaking changes to the major version up.

As for the backward compatibility, see :doc:`user_guide/compatibility`.


.. _contrib-release-cycle:

Release Cycle
~~~~~~~~~~~~~

The first one is the track of **stable versions**, which is a series of revision updates for the latest major version.
The second one is the track of **development versions**, which is a series of pre-releases for the upcoming major version.

Consider that ``X.0.0`` is the latest major version and ``Y.0.0``, ``Z.0.0`` are the succeeding major versions.
Then, the timeline of the updates is depicted by the following table.

========== =========== =========== ============
   Date       ver X       ver Y       ver Z
========== =========== =========== ============
  0 weeks    X.0.0rc1    --         --
  4 weeks    X.0.0       Y.0.0a1    --
  8 weeks    X.1.0*      Y.0.0b1    --
 12 weeks    X.2.0*      Y.0.0rc1   --
 16 weeks    --          Y.0.0      Z.0.0a1
========== =========== =========== ============

(* These might be revision releases)

The dates shown in the left-most column are relative to the release of ``X.0.0rc1``.
In particular, each revision/minor release is made four weeks after the previous one of the same major version, and the pre-release of the upcoming major version is made at the same time.
Whether these releases are revision or minor is determined based on the contents of each update.

Note that there are only three stable releases for the versions ``X.x.x``.
During the parallel development of ``Y.0.0`` and ``Z.0.0a1``, the version ``Y`` is treated as an **almost-stable version** and ``Z`` is treated as a development version.

If there is a critical bug found in ``X.x.x`` after stopping the development of version ``X``, we may release a hot-fix for this version at any time.

We create a milestone for each upcoming release at GitHub.
The GitHub milestone is basically used for collecting the issues and PRs resolved in the release.

.. _contrib-git-branches:

Git Branches
~~~~~~~~~~~~

The ``master`` branch is used to develop pre-release versions.
It means that **alpha, beta, and RC updates are developed at the** ``master`` **branch**.
This branch contains the most up-to-date source tree that includes features newly added after the latest major version.

The stable version is developed at the individual branch named as ``vN`` where "N" reflects the version number (we call it a *versioned branch*).
For example, v1.0.0, v1.0.1, and v1.0.2 will be developed at the ``v1`` branch.

**Notes for contributors:**
When you send a pull request, you basically have to send it to the ``master`` branch.
If the change can also be applied to the stable version, a core team member will apply the same change to the stable version so that the change is also included in the next revision update.

If the change is only applicable to the stable version and not to the ``master`` branch, please send it to the versioned branch.
We basically only accept changes to the latest versioned branch (where the stable version is developed) unless the fix is critical.

If you want to make a new feature of the ``master`` branch available in the current stable version, please send a *backport PR* to the stable version (the latest ``vN`` branch).
See the next section for details.

*Note: a change that can be applied to both branches should be sent to the* ``master`` *branch.*
*Each release of the stable version is also merged to the development version so that the change is also reflected to the next major version.*

Feature Backport PRs
~~~~~~~~~~~~~~~~~~~~

We basically do not backport any new features of the development version to the stable versions.
If you desire to include the feature to the current stable version and you can work on the backport work, we welcome such a contribution.
In such a case, you have to send a backport PR to the latest ``vN`` branch.
**Note that we do not accept any feature backport PRs to older versions because we are not running quality assurance workflows (e.g. CI) for older versions so that we cannot ensure that the PR is correctly ported.**

There are some rules on sending a backport PR.

- Start the PR title from the prefix **[backport]**.
- Clarify the original PR number in the PR description (something like "This is a backport of #XXXX").
- (optional) Write to the PR description the motivation of backporting the feature to the stable version.

Please follow these rules when you create a feature backport PR.

Note: PRs that do not include any changes/additions to APIs (e.g. bug fixes, documentation improvements) are usually backported by core dev members.
It is also appreciated to make such a backport PR by any contributors, though, so that the overall development proceeds more smoothly!

Issues and Pull Requests
------------------------

In this section, we explain how to file issues and send pull requests (PRs).

Issue/PR Labels
~~~~~~~~~~~~~~~

Issues and PRs are labeled by the following tags:

* **Bug**: bug reports (issues) and bug fixes (PRs)
* **Enhancement**: implementation improvements without breaking the interface
* **Feature**: feature requests (issues) and their implementations (PRs)
* **NoCompat**: disrupts backward compatibility
* **Test**: test fixes and updates
* **Document**: document fixes and improvements
* **Example**: fixes and improvements on the examples
* **Install**: fixes installation script
* **Contribution-Welcome**: issues that we request for contribution (only issues are categorized to this)
* **Other**: other issues and PRs

Multiple tags might be labeled to one issue/PR.
**Note that revision releases cannot include PRs in Feature and NoCompat categories.**

How to File an Issue
~~~~~~~~~~~~~~~~~~~~

On registering an issue, write precise explanations on how you want CuPy to be.
Bug reports must include necessary and sufficient conditions to reproduce the bugs.
Feature requests must include **what** you want to do (and **why** you want to do, if needed) with CuPy.
You can contain your thoughts on **how** to realize it into the feature requests, though **what** part is most important for discussions.

.. warning::

   If you have a question on usages of CuPy, it is highly recommended to send a post to `CuPy's Gitter channel <https://gitter.im/cupy/community>`_, `CuPy User Group <https://groups.google.com/forum/#!forum/cupy>`_ or `StackOverflow <https://stackoverflow.com/questions/tagged/cupy>`_ instead of the issue tracker.
   The issue tracker is not a place to share knowledge on practices.
   We may suggest these places and immediately close how-to question issues.

How to Send a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you can write code to fix an issue, we encourage to send a PR.

First of all, before starting to write any code, do not forget to confirm the following points.

- Read through the :ref:`coding-guide` and :ref:`testing-guide`.
- Check the appropriate branch that you should send the PR following :ref:`contrib-git-branches`.
  If you do not have any idea about selecting a branch, please choose the ``master`` branch.

In particular, **check the branch before writing any code.**
The current source tree of the chosen branch is the starting point of your change.

After writing your code **(including unit tests and hopefully documentations!)**, send a PR on GitHub.
You have to write a precise explanation of **what** and **how** you fix;
it is the first documentation of your code that developers read, which is a very important part of your PR.

Once you send a PR, it is automatically tested on ``GitHub Actions``.
After the automatic test passes, core developers will start reviewing your code.
Note that this automatic PR test only includes CPU tests.

.. note::

   We are also running continuous integration with GPU tests for the ``master`` branch and the versioned branch of the latest major version.
   Since this service is currently running on our internal server, we do not use it for automatic PR tests to keep the server secure.

If you are planning to add a new feature or modify existing APIs, **it is recommended to open an issue and discuss the design first.**
The design discussion needs lower cost for the core developers than code review.
Following the consequences of the discussions, you can send a PR that is smoothly reviewed in a shorter time.

Even if your code is not complete, you can send a pull request as a *work-in-progress PR* by putting the ``[WIP]`` prefix to the PR title.
If you write a precise explanation about the PR, core developers and other contributors can join the discussion about how to proceed the PR.
WIP PR is also useful to have discussions based on a concrete code.


.. _coding-guide:

Coding Guidelines
-----------------

.. note::

   Coding guidelines are updated at v5.0.
   Those who have contributed to older versions should read the guidelines again.

We use `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ and a part of `OpenStack Style Guidelines <https://docs.openstack.org/developer/hacking/>`_ related to general coding style as our basic style guidelines.

You can use ``autopep8`` and ``flake8`` commands to check your code.

In order to avoid confusion from using different tool versions, we pin the versions of those tools.
Install them with the following command (from within the top directory of CuPy repository)::

  $ pip install -e '.[stylecheck]'

And check your code with::

  $ autopep8 path/to/your/code.py
  $ flake8 path/to/your/code.py

To check Cython code, use ``.flake8.cython`` configuration file::

  $ flake8 --config=.flake8.cython path/to/your/cython/code.pyx

The ``autopep8`` supports automatically correct Python code to conform to the PEP 8 style guide::

  $ autopep8 --in-place path/to/your/code.py

The ``flake8`` command lets you know the part of your code not obeying our style guidelines.
Before sending a pull request, be sure to check that your code passes the ``flake8`` checking.

Note that ``flake8`` command is not perfect.
It does not check some of the style guidelines.
Here is a (not-complete) list of the rules that ``flake8`` cannot check.

* Relative imports are prohibited. [H304]
* Importing non-module symbols is prohibited.
* Import statements must be organized into three parts: standard libraries, third-party libraries, and internal imports. [H306]

In addition, we restrict the usage of *shortcut symbols* in our code base.
They are symbols imported by packages and sub-packages of ``cupy``.
For example, ``cupy.cuda.Device`` is a shortcut of ``cupy.cuda.device.Device``.
**It is not allowed to use such shortcuts in the ``cupy`` library implementation**.
Note that you can still use them in :tree:`tests` and :tree:`examples` directories.

Once you send a pull request, your coding style is automatically checked by `GitHub Actions`.
The reviewing process starts after the check passes.

The CuPy is designed based on NumPy's API design. CuPy's source code and documents contain the original NumPy ones.
Please note the followings when writing the document.

* In order to identify overlapping parts, it is preferable to add some remarks
  that this document is just copied or altered from the original one. It is
  also preferable to briefly explain the specification of the function in a
  short paragraph, and refer to the corresponding function in NumPy so that
  users can read the detailed document. However, it is possible to include a
  complete copy of the document with such a remark if users cannot summarize
  in such a way.
* If a function in CuPy only implements a limited amount of features in the
  original one, users should explicitly describe only what is implemented in
  the document.

For changes that modify or add new Cython files, please make sure the pointer types follow these guidelines (`#1913 <https://github.com/cupy/cupy/issues/1913>`_).

* Pointers should be ``void*`` if only used within Cython, or ``intptr_t`` if exposed to the Python space.
* Memory sizes should be ``size_t``.
* Memory offsets should be ``ptrdiff_t``.

.. note::

     We are incrementally enforcing the above rules, so some existing code may not follow the above guidelines, but please ensure all new contributions do.

.. _testing-guide:

Unit Testing
------------

Testing is one of the most important part of your code.
You must write test cases and verify your implementation by following our testing guide.

Note that we are using pytest and mock package for testing, so install them before writing your code::

  $ pip install pytest mock

How to Run Tests
~~~~~~~~~~~~~~~~

In order to run unit tests at the repository root, you first have to build Cython files in place by running the following command::

  $ pip install -e .

.. note::

  When you modify ``*.pxd`` files, before running ``pip install -e .``, you must clean ``*.cpp`` and ``*.so`` files once with the following command, because Cython does not automatically rebuild those files nicely::

    $ git clean -fdx

Once Cython modules are built, you can run unit tests by running the following command at the repository root::

  $ python -m pytest

CUDA must be installed to run unit tests.

Some GPU tests require cuDNN to run.
In order to skip unit tests that require cuDNN, specify ``-m='not cudnn'`` option::

  $ python -m pytest path/to/your/test.py -m='not cudnn'

Some GPU tests involve multiple GPUs.
If you want to run GPU tests with insufficient number of GPUs, specify the number of available GPUs to ``CUPY_TEST_GPU_LIMIT``.
For example, if you have only one GPU, launch ``pytest`` by the following command to skip multi-GPU tests::

  $ export CUPY_TEST_GPU_LIMIT=1
  $ python -m pytest path/to/gpu/test.py

Following this naming convention, you can run all the tests by running the following command at the repository root::

  $ python -m pytest

Or you can also specify a root directory to search test scripts from::

  $ python -m pytest tests/cupy_tests     # to just run tests of CuPy
  $ python -m pytest tests/install_tests  # to just run tests of installation modules

If you modify the code related to existing unit tests, you must run appropriate commands.

Test File and Directory Naming Conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tests are put into the :tree:`tests/cupy_tests` directory.
In order to enable test runner to find test scripts correctly, we are using special naming convention for the test subdirectories and the test scripts.

* The name of each subdirectory of ``tests`` must end with the ``_tests`` suffix.
* The name of each test script must start with the ``test_`` prefix.

When we write a test for a module, we use the appropriate path and file name for the test script whose correspondence to the tested module is clear.
For example, if you want to write a test for a module ``cupy.x.y.z``, the test script must be located at ``tests/cupy_tests/x_tests/y_tests/test_z.py``.

How to Write Tests
~~~~~~~~~~~~~~~~~~

There are many examples of unit tests under the :tree:`tests` directory, so reading some of them is a good and recommended way to learn how to write tests for CuPy.
They simply use the :mod:`unittest` package of the standard library, while some tests are using utilities from :mod:`cupy.testing`.

In addition to the :ref:`coding-guide` mentioned above, the following rules are applied to the test code:

* All test classes must inherit from :class:`unittest.TestCase`.
* Use :mod:`unittest` features to write tests, except for the following cases:

    * Use ``assert`` statement instead of ``self.assert*`` methods (e.g., write ``assert x == 1`` instead of ``self.assertEqual(x, 1)``).
    * Use ``with pytest.raises(...):`` instead of ``with self.assertRaises(...):``.

.. note::

   We are incrementally applying the above style.
   Some existing tests may be using the old style (``self.assertRaises``, etc.), but all newly written tests should follow the above style.

Even if your patch includes GPU-related code, your tests should not fail without GPU capability.
Test functions that require CUDA must be tagged by the ``cupy.testing.attr.gpu``::

  import unittest
  from cupy.testing import attr

  class TestMyFunc(unittest.TestCase):
      ...

      @attr.gpu
      def test_my_gpu_func(self):
          ...

The functions tagged by the ``gpu`` decorator are skipped if ``CUPY_TEST_GPU_LIMIT=0`` environment variable is set.
We also have the ``cupy.testing.attr.cudnn`` decorator to let ``pytest`` know that the test depends on cuDNN.
The test functions decorated by ``cudnn`` are skipped if ``-m='not cudnn'`` is given.

The test functions decorated by ``gpu`` must not depend on multiple GPUs.
In order to write tests for multiple GPUs, use ``cupy.testing.attr.multi_gpu()`` decorators instead::

  import unittest
  from cupy.testing import attr

  class TestMyFunc(unittest.TestCase):
      ...

      @attr.multi_gpu(2)  # specify the number of required GPUs here
      def test_my_two_gpu_func(self):
          ...

If your test requires too much time, add ``cupy.testing.attr.slow`` decorator.
The test functions decorated by ``slow`` are skipped if ``-m='not slow'`` is given::

  import unittest
  from cupy.testing import attr

  class TestMyFunc(unittest.TestCase):
      ...

      @attr.slow
      def test_my_slow_func(self):
          ...

.. note::
   If you want to specify more than two attributes, use ``and`` operator like ``-m='not cudnn and not slow'``.
   See detail in `the document of pytest <https://docs.pytest.org/en/latest/example/markers.html>`_.

Once you send a pull request, `Travis-CI <https://travis-ci.org/cupy/cupy/>`_ automatically checks if your code meets our coding guidelines described above.
Since Travis-CI does not support CUDA, we cannot run unit tests automatically.
The reviewing process starts after the automatic check passes.
Note that reviewers will test your code without the option to check CUDA-related code.

.. note::
   Some of numerically unstable tests might cause errors irrelevant to your changes.
   In such a case, we ignore the failures and go on to the review process, so do not worry about it!


Documentation
-------------

When adding a new feature to the framework, you also need to document it in the reference.

.. note::

   If you are unsure about how to fix the documentation, you can submit a pull request without doing so.
   Reviewers will help you fix the documentation appropriately.

The documentation source is stored under `docs directory <https://github.com/cupy/cupy/tree/master/docs>`_ and written in `reStructuredText <http://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_ format.

To build the documentation, you need to install `Sphinx <http://www.sphinx-doc.org/>`_::

  $ pip install -r docs/requirements.txt

Then you can build the documentation in HTML format locally::

  $ cd docs
  $ make html

HTML files are generated under ``build/html`` directory.
Open ``index.html`` with the browser and see if it is rendered as expected.

.. note::

   Docstrings (documentation comments in the source code) are collected from the installed CuPy module.
   If you modified docstrings, make sure to install the module (e.g., using `pip install -e .`) before building the documentation.


Tips for Developers
-------------------

Here are some tips for developers hacking CuPy source code.

Install as Editable
~~~~~~~~~~~~~~~~~~~

During the development we recommend using ``pip`` with ``-e`` option to install as editable mode::

  $ pip install -e .

Please note that even with ``-e``, you will have to rerun ``pip install -e .`` to regenerate C++ sources using Cython if you modified Cython source files (e.g., ``*.pyx`` files).

Use ccache
~~~~~~~~~~

``NVCC`` environment variable can be specified at the build time to use the custom command instead of ``nvcc`` .
You can speed up the rebuild using `ccache <https://ccache.dev/>`_ (v3.4 or later) by::

  $ export NVCC='ccache nvcc'

Limit Architecture
~~~~~~~~~~~~~~~~~~

Use ``CUPY_NVCC_GENERATE_CODE`` environment variable to reduce the build time by limiting the target CUDA architectures.
For example, if you only run your CuPy build with NVIDIA P100 and V100, you can use::

  $ export CUPY_NVCC_GENERATE_CODE=arch=compute_60,code=sm_60;arch=compute_70,code=sm_70

See :doc:`reference/environment` for the description.
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   ..
      Methods

{% block methods %}

   .. rubric:: Methods

   ..
      Special methods

{% for item in ('__call__', '__enter__', '__exit__', '__getitem__', '__setitem__', '__len__', '__next__', '__iter__', '__copy__') %}
{% if item in all_methods or item in all_attributes %}
   .. automethod:: {{ item }}
{% endif %}
{%- endfor %}

   ..
      Ordinary methods

{% for item in methods %}
{% if item not in ('__init__',) %}
   .. automethod:: {{ item }}
{% endif %}
{%- endfor %}

   ..
      Special methods

{% for item in ('__eq__', '__ne__', '__lt__', '__le__', '__gt__', '__ge__', '__nonzero__', '__bool__') %}
{% if item in all_methods %}
   .. automethod:: {{ item }}
{% endif %}
{%- endfor %}
{% endblock %}

   ..
      Atributes

{% block attributes %} {% if attributes %}

   .. rubric:: Attributes

{% for item in attributes %}
   .. autoattribute:: {{ item }}
{%- endfor %}
{% endif %} {% endblock %}
Memory Management
=================

CuPy uses *memory pool* for memory allocations by default.
The memory pool significantly improves the performance by mitigating the overhead of memory allocation and CPU/GPU synchronization.

There are two different memory pools in CuPy:

* Device memory pool (GPU device memory), which is used for GPU memory allocations.
* Pinned memory pool (non-swappable CPU memory), which is used during CPU-to-GPU data transfer.

.. attention::

   When you monitor the memory usage (e.g., using ``nvidia-smi`` for GPU memory or ``ps`` for CPU memory), you may notice that memory not being freed even after the array instance become out of scope.
   This is an expected behavior, as the default memory pool "caches" the allocated memory blocks.

See :doc:`../reference/cuda` for the details of memory management APIs.

For using pinned memory more conveniently, we also provide a few high-level APIs in the ``cupyx`` namespace,
including :func:`cupyx.empty_pinned`, :func:`cupyx.empty_like_pinned`, :func:`cupyx.zeros_pinned`, and
:func:`cupyx.zeros_like_pinned`. They return NumPy arrays backed by pinned memory. If CuPy's pinned memory pool
is in use, the pinned memory is allocated from the pool.

.. note::

    CuPy v8 and above provides a :ref:`FFT plan cache <fft_plan_cache>` that could use a portion of device memory if FFT and related functions are used.
    The memory taken can be released by shrinking or disabling the cache.


Memory Pool Operations
----------------------

The memory pool instance provides statistics about memory allocation.
To access the default memory pool instance, use :func:`cupy.get_default_memory_pool` and :func:`cupy.get_default_pinned_memory_pool`.
You can also free all unused memory blocks hold in the memory pool.
See the example code below for details:

.. code-block:: py

   import cupy
   import numpy

   mempool = cupy.get_default_memory_pool()
   pinned_mempool = cupy.get_default_pinned_memory_pool()

   # Create an array on CPU.
   # NumPy allocates 400 bytes in CPU (not managed by CuPy memory pool).
   a_cpu = numpy.ndarray(100, dtype=numpy.float32)
   print(a_cpu.nbytes)                      # 400

   # You can access statistics of these memory pools.
   print(mempool.used_bytes())              # 0
   print(mempool.total_bytes())             # 0
   print(pinned_mempool.n_free_blocks())    # 0

   # Transfer the array from CPU to GPU.
   # This allocates 400 bytes from the device memory pool, and another 400
   # bytes from the pinned memory pool.  The allocated pinned memory will be
   # released just after the transfer is complete.  Note that the actual
   # allocation size may be rounded to larger value than the requested size
   # for performance.
   a = cupy.array(a_cpu)
   print(a.nbytes)                          # 400
   print(mempool.used_bytes())              # 512
   print(mempool.total_bytes())             # 512
   print(pinned_mempool.n_free_blocks())    # 1

   # When the array goes out of scope, the allocated device memory is released
   # and kept in the pool for future reuse.
   a = None  # (or `del a`)
   print(mempool.used_bytes())              # 0
   print(mempool.total_bytes())             # 512
   print(pinned_mempool.n_free_blocks())    # 1

   # You can clear the memory pool by calling `free_all_blocks`.
   mempool.free_all_blocks()
   pinned_mempool.free_all_blocks()
   print(mempool.used_bytes())              # 0
   print(mempool.total_bytes())             # 0
   print(pinned_mempool.n_free_blocks())    # 0

See :class:`cupy.cuda.MemoryPool` and :class:`cupy.cuda.PinnedMemoryPool` for details.

Limiting GPU Memory Usage
-------------------------

You can hard-limit the amount of GPU memory that can be allocated by using ``CUPY_GPU_MEMORY_LIMIT`` environment variable (see :doc:`../reference/environment` for details).

.. code-block:: py

   # Set the hard-limit to 1 GiB:
   #   $ export CUPY_GPU_MEMORY_LIMIT="1073741824"

   # You can also specify the limit in fraction of the total amount of memory
   # on the GPU. If you have a GPU with 2 GiB memory, the following is
   # equivalent to the above configuration.
   #   $ export CUPY_GPU_MEMORY_LIMIT="50%"

   import cupy
   print(cupy.get_default_memory_pool().get_limit())  # 1073741824

You can also set the limit (or override the value specified via the environment variable) using :meth:`cupy.cuda.MemoryPool.set_limit`.
In this way, you can use a different limit for each GPU device.

.. code-block:: py

   import cupy

   mempool = cupy.get_default_memory_pool()

   with cupy.cuda.Device(0):
       mempool.set_limit(size=1024**3)  # 1 GiB

   with cupy.cuda.Device(1):
       mempool.set_limit(size=2*1024**3)  # 2 GiB

.. note::

   CUDA allocates some GPU memory outside of the memory pool (such as CUDA context, library handles, etc.).
   Depending on the usage, such memory may take one to few hundred MiB.
   That will not be counted in the limit.

Changing Memory Pool
--------------------

You can use your own memory allocator instead of the default memory pool by passing the memory allocation function to :func:`cupy.cuda.set_allocator` / :func:`cupy.cuda.set_pinned_memory_allocator`.
The memory allocator function should take 1 argument (the requested size in bytes) and return :class:`cupy.cuda.MemoryPointer` / :class:`cupy.cuda.PinnedMemoryPointer`.

CuPy provides two such allocators for using managed memory and stream ordered memory on GPU,
see :func:`cupy.cuda.malloc_managed` and :func:`cupy.cuda.malloc_async`, respectively, for details.
To enable a memory pool backed by managed memory, you can construct a new :class:`~cupy.cuda.MemoryPool` instance with its allocator
set to :func:`~cupy.cuda.malloc_managed` as follows

.. code-block:: py

    import cupy

    # Use managed memory
    cupy.cuda.set_allocator(cupy.cuda.MemoryPool(cupy.cuda.malloc_managed).malloc)

Note that if you pass :func:`~cupy.cuda.malloc_managed` directly to :func:`~cupy.cuda.set_allocator` without constructing
a :class:`~cupy.cuda.MemoryPool` instance, when the memory is freed it will be released back to the system immediately,
which may or may not be desired.

Stream Ordered Memory Allocator is a new feature added since CUDA 11.2. CuPy provides an *experimental* interface to it.
Similar to CuPy's memory pool, Stream Ordered Memory Allocator also allocates/deallocates memory *asynchronously* from/to
a memory pool in a stream-ordered fashion. The key difference is that it is a built-in feature implemented in the CUDA
driver by NVIDIA, so other CUDA applications in the same processs can easily allocate memory from the same pool.

To enable a memory pool that manages stream ordered memory, you can construct a new :class:`~cupy.cuda.MemoryAsyncPool`
instance:

.. code-block:: py

    import cupy

    # Use asynchronous stream ordered memory
    cupy.cuda.set_allocator(cupy.cuda.MemoryAsyncPool().malloc)

    # Create a custom stream
    s = cupy.cuda.Stream()

    # This would allocate memory asynchronously on stream s
    with s:
        a = cupy.empty((100,), dtype=cupy.float64)

Note that in this case we do not use the :class:`~cupy.cuda.MemoryPool` class. The :class:`~cupy.cuda.MemoryAsyncPool` takes
a different input argument from that of :class:`~cupy.cuda.MemoryPool` to indicate which pool to use.
Please refer to :class:`~cupy.cuda.MemoryAsyncPool`'s documentation for further detail.

Note that if you pass :func:`~cupy.cuda.malloc_async` directly to :func:`~cupy.cuda.set_allocator` without constructing
a :class:`~cupy.cuda.MemoryAsyncPool` instance, the device's *current* memory pool will be used.

When using stream ordered memory, it is important that you maintain a correct stream semantics yourselves using, for example,
the :class:`~cupy.cuda.Stream` and :class:`~cupy.cuda.Event` APIs (see :ref:`cuda_stream_event` for details); CuPy does not
attempt to act smartly for you. Upon deallocation, the memory is freed asynchronously either on the stream it was
allocated (first attempt), or on any current CuPy stream (second attempt). It is permitted that the stream on which the
memory was allocated gets destroyed before all memory allocated on it is freed.

In addition, applications/libraries internally use ``cudaMalloc`` (CUDA's default, synchronous allocator) could have unexpected
interplay with Stream Ordered Memory Allocator. Specifically, memory freed to the memory pool might not be immediately visible
to ``cudaMalloc``, leading to potential out-of-memory errors. In this case, you can either call :meth:`~cupy.cuda.MemoryAsyncPool.free_all_blocks()`
or just manually perform a (event/stream/device) synchronization, and retry.

Currently the :class:`~cupy.cuda.MemoryAsyncPool` interface is *experimental*. In particular, while its API is largely identical
to that of :class:`~cupy.cuda.MemoryPool`, several of the pool's methods require a sufficiently new driver (and of course, a
supported hardware, CUDA version, and platform) due to CUDA's limitation.

You can even disable the default memory pool by the code below.
Be sure to do this before any other CuPy operations.

.. code-block:: py

   import cupy

   # Disable memory pool for device memory (GPU)
   cupy.cuda.set_allocator(None)

   # Disable memory pool for pinned memory (CPU).
   cupy.cuda.set_pinned_memory_allocator(None)
Interoperability
================

CuPy can be used in conjunction with other libraries.


CUDA functionalities
--------------------

Under construction. For using CUDA streams created in foreign libraries in CuPy, see :ref:`cuda_stream_event`.


NumPy
-----

:class:`cupy.ndarray` implements ``__array_ufunc__`` interface (see `NEP 13 — A Mechanism for Overriding Ufuncs <http://www.numpy.org/neps/nep-0013-ufunc-overrides.html>`_ for details).
This enables NumPy ufuncs to be directly operated on CuPy arrays.
``__array_ufunc__`` feature requires NumPy 1.13 or later.

.. code:: python

    import cupy
    import numpy

    arr = cupy.random.randn(1, 2, 3, 4).astype(cupy.float32)
    result = numpy.sum(arr)
    print(type(result))  # => <class 'cupy._core.core.ndarray'>

:class:`cupy.ndarray` also implements ``__array_function__`` interface (see `NEP 18 — A dispatch mechanism for NumPy’s high level array functions <http://www.numpy.org/neps/nep-0018-array-function-protocol.html>`_ for details).
This enables code using NumPy to be directly operated on CuPy arrays.
``__array_function__`` feature requires NumPy 1.16 or later; As of NumPy 1.17, ``__array_function__`` is enabled by default.


Numba
-----

`Numba <https://numba.pydata.org/>`_ is a Python JIT compiler with NumPy support.

:class:`cupy.ndarray` implements ``__cuda_array_interface__``, which is the CUDA array interchange interface compatible with Numba v0.39.0 or later (see `CUDA Array Interface <https://numba.readthedocs.io/en/stable/cuda/cuda_array_interface.html>`_ for details).
It means you can pass CuPy arrays to kernels JITed with Numba.
The following is a simple example code borrowed from `numba/numba#2860 <https://github.com/numba/numba/pull/2860>`_:

.. code:: python

	import cupy
	from numba import cuda

	@cuda.jit
	def add(x, y, out):
		start = cuda.grid(1)
		stride = cuda.gridsize(1)
		for i in range(start, x.shape[0], stride):
			out[i] = x[i] + y[i]

	a = cupy.arange(10)
	b = a * 2
	out = cupy.zeros_like(a)

	print(out)  # => [0 0 0 0 0 0 0 0 0 0]

	add[1, 32](a, b, out)

	print(out)  # => [ 0  3  6  9 12 15 18 21 24 27]

In addition, :func:`cupy.asarray` supports zero-copy conversion from Numba CUDA array to CuPy array.

.. code:: python

    import numpy
    import numba
    import cupy

    x = numpy.arange(10)  # type: numpy.ndarray
    x_numba = numba.cuda.to_device(x)  # type: numba.cuda.cudadrv.devicearray.DeviceNDArray
    x_cupy = cupy.asarray(x_numba)  # type: cupy.ndarray

.. warning::

    ``__cuda_array_interface__`` specifies that the object lifetime must be managed by the user, so it is an undefined behavior if the
    exported object is destroyed while still in use by the consumer library.

.. note::

    CuPy uses two environment variables controlling the exchange behavior: :envvar:`CUPY_CUDA_ARRAY_INTERFACE_SYNC` and :envvar:`CUPY_CUDA_ARRAY_INTERFACE_EXPORT_VERSION`.


mpi4py
------

`MPI for Python (mpi4py) <https://mpi4py.readthedocs.io/en/latest/>`_ is a Python wrapper for the Message Passing Interface (MPI) libraries.

MPI is the most widely used standard for high-performance inter-process communications. Recently several MPI vendors, including MPICH, Open MPI and MVAPICH, have extended their support beyond the MPI-3.1 standard to enable "CUDA-awareness"; that is, passing CUDA device pointers directly to MPI calls to avoid explicit data movement between the host and the device.

With the ``__cuda_array_interface__`` (as mentioned above) and ``DLPack`` data exchange protocols (see :ref:`dlpack` below) implemented in CuPy, mpi4py now provides (experimental) support for passing CuPy arrays to MPI calls, provided that mpi4py is built against a CUDA-aware MPI implementation. The following is a simple example code borrowed from `mpi4py Tutorial <https://mpi4py.readthedocs.io/en/latest/tutorial.html>`_:

.. code:: python

    # To run this script with N MPI processes, do
    # mpiexec -n N python this_script.py

    import cupy
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    size = comm.Get_size()

    # Allreduce
    sendbuf = cupy.arange(10, dtype='i')
    recvbuf = cupy.empty_like(sendbuf)
    comm.Allreduce(sendbuf, recvbuf)
    assert cupy.allclose(recvbuf, sendbuf*size)

This new feature is added since mpi4py 3.1.0. See the `mpi4py website <https://mpi4py.readthedocs.io/en/latest/>`_ for more information.


PyTorch
-------

`PyTorch <https://pytorch.org/>`_ is a machine learning framefork that provides high-performance, differentiable tensor operations.

PyTorch also supports ``__cuda_array_interface__``, so zero-copy data exchange between CuPy and PyTorch can be achieved at no cost.
The only caveat is PyTorch by default creates CPU tensors, which do not have the ``__cuda_array_interface__`` property defined, and
users need to ensure the tensor is already on GPU before exchanging.

.. code:: python

    >>> import cupy as cp
    >>> import torch
    >>>
    >>> # convert a torch tensor to a cupy array
    >>> a = torch.rand((4, 4), device='cuda')
    >>> b = cp.asarray(a)
    >>> b *= b
    >>> b
    array([[0.8215962 , 0.82399917, 0.65607935, 0.30354425],
           [0.422695  , 0.8367199 , 0.00208597, 0.18545236],
           [0.00226746, 0.46201342, 0.6833052 , 0.47549972],
           [0.5208748 , 0.6059282 , 0.1909013 , 0.5148635 ]], dtype=float32)
    >>> a
    tensor([[0.8216, 0.8240, 0.6561, 0.3035],
            [0.4227, 0.8367, 0.0021, 0.1855],
            [0.0023, 0.4620, 0.6833, 0.4755],
            [0.5209, 0.6059, 0.1909, 0.5149]], device='cuda:0')
    >>> # check the underlying memory pointer is the same
    >>> assert a.__cuda_array_interface__['data'][0] == b.__cuda_array_interface__['data'][0]
    >>>
    >>> # convert a cupy array to a torch tensor
    >>> a = cp.arange(10)
    >>> b = torch.as_tensor(a, device='cuda')
    >>> b += 3
    >>> b
    tensor([ 3,  4,  5,  6,  7,  8,  9, 10, 11, 12], device='cuda:0')
    >>> a
    array([ 3,  4,  5,  6,  7,  8,  9, 10, 11, 12])
    >>> assert a.__cuda_array_interface__['data'][0] == b.__cuda_array_interface__['data'][0]

PyTorch also supports zero-copy data exchange through ``DLPack`` (see :ref:`dlpack` below):

.. code:: python

	import cupy
	import torch

	from torch.utils.dlpack import to_dlpack
	from torch.utils.dlpack import from_dlpack

	# Create a PyTorch tensor.
	tx1 = torch.randn(1, 2, 3, 4).cuda()

	# Convert it into a DLPack tensor.
	dx = to_dlpack(tx1)

	# Convert it into a CuPy array.
	cx = cupy.from_dlpack(dx)

	# Convert it back to a PyTorch tensor.
	tx2 = from_dlpack(cx.toDlpack())

`pytorch-pfn-extras <https://github.com/pfnet/pytorch-pfn-extras/>`_ library provides additional integration features with PyTorch, including memory pool sharing and stream sharing:

.. code:: python

   >>> import cupy
   >>> import torch
   >>> import pytorch_pfn_extras as ppe
   >>>
   >>> # Perform CuPy memory allocation using the PyTorch memory pool.
   >>> ppe.cuda.use_torch_mempool_in_cupy()
   >>> torch.cuda.memory_allocated()
   0
   >>> arr = cupy.arange(10)
   >>> torch.cuda.memory_allocated()
   512
   >>>
   >>> # Change the default stream in PyTorch and CuPy:
   >>> stream = torch.cuda.Stream()
   >>> with ppe.cuda.stream(stream):
   ...     ...


Using custom kernels in PyTorch
*******************************

With the DLPack protocol, it becomes very simple to implement functions in PyTorch using CuPy user-defined kernels. Below is the example of a PyTorch autograd function
that computes the forward and backward pass of the logarithm using :class:`cupy.RawKernel` s.

.. code:: python

    import cupy
    import torch
    
    
    cupy_custom_kernel_fwd = cupy.RawKernel(
        r"""
    extern "C" __global__
    void cupy_custom_kernel_fwd(const float* x, float* y, int size) {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        if (tid < size)
            y[tid] = log(x[tid]);
    }
    """,
        "cupy_custom_kernel_fwd",
    )
    
    
    cupy_custom_kernel_bwd = cupy.RawKernel(
        r"""
    extern "C" __global__
    void cupy_custom_kernel_bwd(const float* x, float* gy, float* gx, int size) {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        if (tid < size)
            gx[tid] = gy[tid] / x[tid];
    }
    """,
        "cupy_custom_kernel_bwd",
    )
    
    
    class CuPyLog(torch.autograd.Function):
        @staticmethod
        def forward(ctx, x):
            ctx.input = x
            # Enforce contiguous arrays to simplify RawKernel indexing.
            cupy_x = cupy.ascontiguousarray(cupy.from_dlpack(x.detach()))
            cupy_y = cupy.empty(cupy_x.shape, dtype=cupy_x.dtype)
            x_size = cupy_x.size
            bs = 128
            cupy_custom_kernel_fwd(
                (bs,), ((x_size + bs - 1) // bs,), (cupy_x, cupy_y, x_size)
            )
            # the ownership of the device memory backing cupy_y is implicitly
            # transferred to torch_y, so this operation is safe even after
            # going out of scope of this function.
            torch_y = torch.from_dlpack(cupy_y)
            return torch_y
    
        @staticmethod
        def backward(ctx, grad_y):
            # Enforce contiguous arrays to simplify RawKernel indexing.
            cupy_input = cupy.from_dlpack(ctx.input.detach()).ravel()
            cupy_grad_y = cupy.from_dlpack(grad_y.detach()).ravel()
            cupy_grad_x = cupy.zeros(cupy_grad_y.shape, dtype=cupy_grad_y.dtype)
            gy_size = cupy_grad_y.size
            bs = 128
            cupy_custom_kernel_bwd(
                (bs,),
                ((gy_size + bs - 1) // bs,),
                (cupy_input, cupy_grad_y, cupy_grad_x, gy_size),
            )
            # the ownership of the device memory backing cupy_grad_x is implicitly
            # transferred to torch_y, so this operation is safe even after
            # going out of scope of this function.
            torch_grad_x = torch.from_dlpack(cupy_grad_x)
            return torch_grad_x

.. note::

   Directly feeding a ``torch.Tensor`` to :func:`cupy.from_dlpack` is only supported in the (new) DLPack data exchange protocol added in CuPy v10+ and PyTorch 1.10+.
   For earlier versions, you will need to wrap the ``Tensor`` with ``torch.utils.dlpack.to_dlpack()`` as shown in the above examples.

RMM
---

`RMM (RAPIDS Memory Manager) <https://docs.rapids.ai/api/rmm/stable/index.html>`_ provides highly configurable memory allocators.

RMM provides an interface to allow CuPy to allocate memory from the RMM memory pool instead of from CuPy's own pool. It can be set up
as simple as:

.. code:: python

    import cupy
    import rmm
    cupy.cuda.set_allocator(rmm.rmm_cupy_allocator)

Sometimes, a more performant allocator may be desirable. RMM provides an option to switch the allocator:

.. code:: python

    import cupy
    import rmm
    rmm.reinitialize(pool_allocator=True)  # can also set init pool size etc here
    cupy.cuda.set_allocator(rmm.rmm_cupy_allocator)

For more information on CuPy's memory management, see :doc:`./memory`.


.. _dlpack:

DLPack
------

`DLPack <https://github.com/dmlc/dlpack>`__ is a specification of tensor structure to share tensors among frameworks.

CuPy supports importing from and exporting to DLPack data structure (:func:`cupy.from_dlpack` and :func:`cupy.ndarray.toDlpack`).

Here is a simple example:

.. code:: python

	import cupy

	# Create a CuPy array.
	cx1 = cupy.random.randn(1, 2, 3, 4).astype(cupy.float32)

	# Convert it into a DLPack tensor.
	dx = cx1.toDlpack()

	# Convert it back to a CuPy array.
	cx2 = cupy.from_dlpack(dx)

`TensorFlow <https://www.tensorflow.org>`_ also supports DLpack, so zero-copy data exchange between CuPy and TensorFlow through
DLPack is possible:

.. code:: python

    >>> import tensorflow as tf
    >>> import cupy as cp
    >>>
    >>> # convert a TF tensor to a cupy array
    >>> with tf.device('/GPU:0'):
    ...     a = tf.random.uniform((10,))
    ...
    >>> a
    <tf.Tensor: shape=(10,), dtype=float32, numpy=
    array([0.9672388 , 0.57568085, 0.53163004, 0.6536236 , 0.20479882,
           0.84908986, 0.5852566 , 0.30355775, 0.1733712 , 0.9177849 ],
          dtype=float32)>
    >>> a.device
    '/job:localhost/replica:0/task:0/device:GPU:0'
    >>> cap = tf.experimental.dlpack.to_dlpack(a)
    >>> b = cp.from_dlpack(cap)
    >>> b *= 3
    >>> b
    array([1.4949363 , 0.60699713, 1.3276931 , 1.5781245 , 1.1914308 ,
           2.3180873 , 1.9560868 , 1.3932796 , 1.9299742 , 2.5352407 ],
          dtype=float32)
    >>> a
    <tf.Tensor: shape=(10,), dtype=float32, numpy=
    array([1.4949363 , 0.60699713, 1.3276931 , 1.5781245 , 1.1914308 ,
           2.3180873 , 1.9560868 , 1.3932796 , 1.9299742 , 2.5352407 ],
          dtype=float32)>
    >>>
    >>> # convert a cupy array to a TF tensor
    >>> a = cp.arange(10)
    >>> cap = a.toDlpack()
    >>> b = tf.experimental.dlpack.from_dlpack(cap)
    >>> b.device
    '/job:localhost/replica:0/task:0/device:GPU:0'
    >>> b
    <tf.Tensor: shape=(10,), dtype=int64, numpy=array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])>
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

Be aware that in TensorFlow all tensors are immutable, so in the latter case any changes in ``b`` cannot be reflected in the CuPy array ``a``.

Note that as of DLPack v0.5 for correctness the above approach (implicitly) requires users to ensure that such conversion (both importing and exporting a CuPy array) must happen on the same CUDA/HIP stream. If in doubt, the current CuPy stream in use can be fetched by, for example, calling :func:`cupy.cuda.get_current_stream`. Please consult the other framework's documentation for how to access and control the streams.

DLPack data exchange protocol
*****************************

To obviate user-managed streams and DLPack tensor objects, the `DLPack data exchange protocol <https://data-apis.org/array-api/latest/design_topics/data_interchange.html>`_ provides a mechanism to shift the responsibility from users to libraries. Any compliant objects (such as :class:`cupy.ndarray`) must implement a pair of methods ``__dlpack__`` and ``__dlpack_device__``. The function :func:`cupy.from_dlpack` accepts such object and returns a :class:`cupy.ndarray` that is safely accessible on CuPy's current stream. Likewise, :class:`cupy.ndarray` can be exported via any compliant library's ``from_dlpack()`` function.

.. note::

    CuPy uses :envvar:`CUPY_DLPACK_EXPORT_VERSION` to control how to handle tensors backed by CUDA managed memory.
Performance Best Practices
==========================

Here we gather a few tricks and advices for improving CuPy's performance.

Benchmarking
------------

It is utterly important to first identify the performance bottleneck before making any attempt to optimize
your code. To help set up a baseline benchmark, CuPy provides a useful utility :func:`cupyx.profiler.benchmark`
for timing the elapsed time of a Python function on both CPU and GPU:

.. doctest::

    >>> from cupyx.profiler import benchmark
    >>> 
    >>> def my_func(a):
    ...     return cp.sqrt(cp.sum(a**2, axis=-1))
    ... 
    >>> a = cp.random.random((256, 1024))
    >>> print(benchmark(my_func, (a,), n_repeat=20))  # doctest: +SKIP
    my_func             :    CPU:   44.407 us   +/- 2.428 (min:   42.516 / max:   53.098) us     GPU-0:  181.565 us   +/- 1.853 (min:  180.288 / max:  188.608) us

Because GPU executions run asynchronously with respect to CPU executions, a common pitfall in GPU programming is to mistakenly
measure the elapsed time using CPU timing utilities (such as :py:func:`time.perf_counter` from the Python Standard Library
or the ``%timeit`` magic from IPython), which have no knowledge in the GPU runtime. :func:`cupyx.profiler.benchmark` addresses
this by setting up CUDA events on the :ref:`current_stream` right before and after the function to be measured and
synchronizing over the end event (see :ref:`cuda_stream_event` for detail). Below we sketch what is done internally in :func:`cupyx.profiler.benchmark`:

.. doctest::

    >>> import time
    >>> start_gpu = cp.cuda.Event()
    >>> end_gpu = cp.cuda.Event()
    >>>
    >>> start_gpu.record()
    >>> start_cpu = time.perf_counter()
    >>> out = my_func(a)
    >>> end_cpu = time.perf_counter()
    >>> end_gpu.record()
    >>> end_gpu.synchronize()
    >>> t_gpu = cp.cuda.get_elapsed_time(start_gpu, end_gpu)
    >>> t_cpu = end_cpu - start_cpu

Additionally, :func:`cupyx.profiler.benchmark` runs a few warm-up runs to reduce timing fluctuation and exclude the overhead in first invocations.


One-Time Overheads
~~~~~~~~~~~~~~~~~~

Be aware of these overheads when benchmarking CuPy code.

Context Initialization
......................

It may take several seconds when calling a CuPy function for the first time in a process.
This is because CUDA driver creates a CUDA context during the first CUDA API call in CUDA applications.

Kernel Compilation
..................

CuPy uses on-the-fly kernel synthesis. When a kernel call is required, it compiles a kernel code optimized for the dimensions and dtypes of the given arguments, sends them to the GPU device, and executes the kernel.

CuPy caches the kernel code sent to GPU device within the process, which reduces the kernel compilation time on further calls.

The compiled code is also cached in the directory ``${HOME}/.cupy/kernel_cache`` (the path can be overwritten by setting the :envvar:`CUPY_CACHE_DIR` environment variable).
This allows reusing the compiled kernel binary across the process.


In-depth profiling
------------------

Under construction. To mark with NVTX/rocTX ranges, you can use the :func:`cupyx.profiler.time_range` API. To start/stop the profiler, you can use the :func:`cupyx.profiler.profile` API.


Use CUB/cuTENSOR backends for reduction and other routines
----------------------------------------------------------

For reduction operations (such as :func:`~cupy.sum`, :func:`~cupy.prod`, :func:`~cupy.amin`, :func:`~cupy.amax`, :func:`~cupy.argmin`, :func:`~cupy.argmax`) and many more routines built upon them, CuPy ships with our own implementations so that things just work out of the box. However, there are dedicated efforts to further accelerate these routines, such as `CUB <https://github.com/NVIDIA/cub>`_ and `cuTENSOR <https://developer.nvidia.com/cutensor>`_.

In order to support more performant backends wherever applicable, starting v8 CuPy introduces an environment variable :envvar:`CUPY_ACCELERATORS` to allow users to specify the desired backends (and in what order they are tried). For example, consider summing over a 256-cubic array:

.. doctest::

    >>> from cupyx.profiler import benchmark
    >>> a = cp.random.random((256, 256, 256), dtype=cp.float32)
    >>> print(benchmark(a.sum, (), n_repeat=100))  # doctest: +SKIP
    sum                 :    CPU:   12.101 us   +/- 0.694 (min:   11.081 / max:   17.649) us     GPU-0:10174.898 us   +/-180.551 (min:10084.576 / max:10595.936) us

We can see that it takes about 10 ms to run (on this GPU). However, if we launch the Python session using ``CUPY_ACCELERATORS=cub python``, we get a ~100x speedup for free (only ~0.1 ms):

.. doctest::

    >>> print(benchmark(a.sum, (), n_repeat=100))  # doctest: +SKIP
    sum                 :    CPU:   20.569 us   +/- 5.418 (min:   13.400 / max:   28.439) us     GPU-0:  114.740 us   +/- 4.130 (min:  108.832 / max:  122.752) us

CUB is a backend shipped together with CuPy.
It also accelerates other routines, such as inclusive scans (ex: :func:`~cupy.cumsum`), histograms,
sparse matrix-vector multiplications (not applicable in CUDA 11), and :class:`~cupy.ReductionKernel`.
cuTENSOR offers optimized performance for binary elementwise ufuncs, reduction and tensor contraction.
If cuTENSOR is installed, setting ``CUPY_ACCELERATORS=cub,cutensor``, for example, would try CUB first and fall back to cuTENSOR if CUB does not provide the needed support. In the case that both backends are not applicable, it falls back to CuPy's default implementation.

Note that while in general the accelerated reductions are faster, there could be exceptions
depending on the data layout. In particular, the CUB reduction only supports reduction over
contiguous axes.
In any case, we recommend to perform some benchmarks to determine whether CUB/cuTENSOR offers
better performance or not.


Overlapping work using streams
------------------------------

Under construction.


Use JIT compiler
----------------

Under construction. For now please refer to :ref:`jit_kernel_definition` for a quick introduction.


Prefer float32 over float64
---------------------------

Under construction.
Difference between CuPy and NumPy
=================================

The interface of CuPy is designed to obey that of NumPy.
However, there are some differences.


Cast behavior from float to integer
-----------------------------------

Some casting behaviors from float to integer are not defined in C++ specification.
The casting from a negative float to unsigned integer and infinity to integer is one of such examples.
The behavior of NumPy depends on your CPU architecture.
This is the result on an Intel CPU:

  >>> np.array([-1], dtype=np.float32).astype(np.uint32)
  array([4294967295], dtype=uint32)
  >>> cupy.array([-1], dtype=np.float32).astype(np.uint32)
  array([0], dtype=uint32)

  >>> np.array([float('inf')], dtype=np.float32).astype(np.int32)
  array([-2147483648], dtype=int32)
  >>> cupy.array([float('inf')], dtype=np.float32).astype(np.int32)
  array([2147483647], dtype=int32)


Random methods support dtype argument
-------------------------------------

NumPy's random value generator does not support a `dtype` argument and instead always returns a ``float64`` value.
We support the option in CuPy because cuRAND, which is used in CuPy, supports both ``float32`` and ``float64``.


  >>> np.random.randn(dtype=np.float32)
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  TypeError: randn() got an unexpected keyword argument 'dtype'
  >>> cupy.random.randn(dtype=np.float32)    # doctest: +SKIP
  array(0.10689262300729752, dtype=float32)


Out-of-bounds indices
---------------------
CuPy handles out-of-bounds indices differently by default from NumPy when
using integer array indexing.
NumPy handles them by raising an error, but CuPy wraps around them.

  >>> x = np.array([0, 1, 2])
  >>> x[[1, 3]] = 10
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  IndexError: index 3 is out of bounds for axis 1 with size 3
  >>> x = cupy.array([0, 1, 2])
  >>> x[[1, 3]] = 10
  >>> x
  array([10, 10,  2])


Duplicate values in indices
---------------------------
CuPy's ``__setitem__`` behaves differently from NumPy when integer arrays
reference the same location multiple times.
In that case, the value that is actually stored is undefined.
Here is an example of CuPy.

  >>> a = cupy.zeros((2,))
  >>> i = cupy.arange(10000) % 2
  >>> v = cupy.arange(10000).astype(np.float32)
  >>> a[i] = v
  >>> a  # doctest: +SKIP
  array([ 9150.,  9151.])

NumPy stores the value corresponding to the
last element among elements referencing duplicate locations.

  >>> a_cpu = np.zeros((2,))
  >>> i_cpu = np.arange(10000) % 2
  >>> v_cpu = np.arange(10000).astype(np.float32)
  >>> a_cpu[i_cpu] = v_cpu
  >>> a_cpu
  array([9998., 9999.])


Zero-dimensional array
-----------------------------------------------

Reduction methods
~~~~~~~~~~~~~~~~~

NumPy's reduction functions (e.g. :func:`numpy.sum`) return scalar values (e.g. :class:`numpy.float32`).
However CuPy counterparts return zero-dimensional :class:`cupy.ndarray` s.
That is because CuPy scalar values (e.g. :class:`cupy.float32`) are aliases of NumPy scalar values and are allocated in CPU memory.
If these types were returned, it would be required to synchronize between GPU and CPU.
If you want to use scalar values, cast the returned arrays explicitly.

  >>> type(np.sum(np.arange(3))) == np.int64
  True
  >>> type(cupy.sum(cupy.arange(3))) == cupy._core.core.ndarray
  True


Type promotion
~~~~~~~~~~~~~~

CuPy automatically promotes dtypes of :class:`cupy.ndarray` s in a function with two or more operands, the result dtype is determined by the dtypes of the inputs.
This is different from NumPy's rule on type promotion, when operands contain zero-dimensional arrays.
Zero-dimensional :class:`numpy.ndarray` s are treated as if they were scalar values if they appear in operands of NumPy's function,
This may affect the dtype of its output, depending on the values of the "scalar" inputs.

  >>> (np.array(3, dtype=np.int32) * np.array([1., 2.], dtype=np.float32)).dtype
  dtype('float32')
  >>> (np.array(300000, dtype=np.int32) * np.array([1., 2.], dtype=np.float32)).dtype
  dtype('float64')
  >>> (cupy.array(3, dtype=np.int32) * cupy.array([1., 2.], dtype=np.float32)).dtype
  dtype('float64')


Matrix type (:class:`numpy.matrix`)
-----------------------------------

SciPy returns :class:`numpy.matrix` (a subclass of :class:`numpy.ndarray`) when dense matrices are computed from sparse matrices (e.g., ``coo_matrix + ndarray``). However, CuPy returns :class:`cupy.ndarray` for such operations.

There is no plan to provide :class:`numpy.matrix` equivalent in CuPy.
This is because the use of :class:`numpy.matrix` is no longer recommended since NumPy 1.15.


Data types
----------

Data type of CuPy arrays cannot be non-numeric like strings or objects.
See :ref:`overview` for details.


Universal Functions only work with CuPy array or scalar
-------------------------------------------------------

Unlike NumPy, Universal Functions in CuPy only work with CuPy array or scalar.
They do not accept other objects (e.g., lists or :class:`numpy.ndarray`).

  >>> np.power([np.arange(5)], 2)
  array([[ 0,  1,  4,  9, 16]])

  >>> cupy.power([cupy.arange(5)], 2)
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  TypeError: Unsupported type <class 'list'>


Random seed arrays are hashed to scalars
----------------------------------------

Like Numpy, CuPy's RandomState objects accept seeds either as numbers or as
full numpy arrays.

  >>> seed = np.array([1, 2, 3, 4, 5])
  >>> rs = cupy.random.RandomState(seed=seed)

However, unlike Numpy, array seeds will be hashed down to a single number and
so may not communicate as much entropy to the underlying random number
generator.


NaN (not-a-number) handling
---------------------------

By default CuPy's reduction functions (e.g., :func:`cupy.sum`) handle NaNs in complex numbers differently from NumPy's
counterparts:

  >>> a = [0.5 + 3.7j, complex(0.7, np.nan), complex(np.nan, -3.9), complex(np.nan, np.nan)]
  >>>
  >>> a_np = np.asarray(a)
  >>> print(a_np.max(), a_np.min())
  (0.7+nanj) (0.7+nanj)
  >>>
  >>> a_cp = cp.asarray(a_np)
  >>> print(a_cp.max(), a_cp.min())
  (nan-3.9j) (nan-3.9j)

The reason is that internally the reduction is performed in a strided fashion, thus it does not ensure a proper
comparison order and cannot follow NumPy's rule to always propagate the first-encountered NaN.
.. _udkernel:

User-Defined Kernels
====================

CuPy provides easy ways to define three types of CUDA kernels: elementwise kernels, reduction kernels and raw kernels.
In this documentation, we describe how to define and call each kernels.


Basics of elementwise kernels
-----------------------------

An elementwise kernel can be defined by the :class:`~cupy.ElementwiseKernel` class.
The instance of this class defines a CUDA kernel which can be invoked by the ``__call__`` method of this instance.

A definition of an elementwise kernel consists of four parts: an input argument list, an output argument list, a loop body code, and the kernel name.
For example, a kernel that computes a squared difference :math:`f(x, y) = (x - y)^2` is defined as follows:

.. doctest::

   >>> squared_diff = cp.ElementwiseKernel(
   ...    'float32 x, float32 y',
   ...    'float32 z',
   ...    'z = (x - y) * (x - y)',
   ...    'squared_diff')

The argument lists consist of comma-separated argument definitions.
Each argument definition consists of a *type specifier* and an *argument name*.
Names of NumPy data types can be used as type specifiers.

.. note::
   ``n``, ``i``, and names starting with an underscore ``_`` are reserved for the internal use.

The above kernel can be called on either scalars or arrays with broadcasting:

.. doctest::

   >>> x = cp.arange(10, dtype=np.float32).reshape(2, 5)
   >>> y = cp.arange(5, dtype=np.float32)
   >>> squared_diff(x, y)
   array([[ 0.,  0.,  0.,  0.,  0.],
          [25., 25., 25., 25., 25.]], dtype=float32)
   >>> squared_diff(x, 5)
   array([[25., 16.,  9.,  4.,  1.],
          [ 0.,  1.,  4.,  9., 16.]], dtype=float32)

Output arguments can be explicitly specified (next to the input arguments):

.. doctest::

   >>> z = cp.empty((2, 5), dtype=np.float32)
   >>> squared_diff(x, y, z)
   array([[ 0.,  0.,  0.,  0.,  0.],
          [25., 25., 25., 25., 25.]], dtype=float32)


Type-generic kernels
--------------------

If a type specifier is one character, then it is treated as a **type placeholder**.
It can be used to define a type-generic kernels.
For example, the above ``squared_diff`` kernel can be made type-generic as follows:

.. doctest::

   >>> squared_diff_generic = cp.ElementwiseKernel(
   ...     'T x, T y',
   ...     'T z',
   ...     'z = (x - y) * (x - y)',
   ...     'squared_diff_generic')

Type placeholders of a same character in the kernel definition indicate the same type.
The actual type of these placeholders is determined by the actual argument type.
The ElementwiseKernel class first checks the output arguments and then the input arguments to determine the actual type.
If no output arguments are given on the kernel invocation, then only the input arguments are used to determine the type.

The type placeholder can be used in the loop body code:

.. doctest::

   >>> squared_diff_generic = cp.ElementwiseKernel(
   ...     'T x, T y',
   ...     'T z',
   ...     '''
   ...         T diff = x - y;
   ...         z = diff * diff;
   ...     ''',
   ...     'squared_diff_generic')

More than one type placeholder can be used in a kernel definition.
For example, the above kernel can be further made generic over multiple arguments:

.. doctest::

   >>> squared_diff_super_generic = cp.ElementwiseKernel(
   ...     'X x, Y y',
   ...     'Z z',
   ...     'z = (x - y) * (x - y)',
   ...     'squared_diff_super_generic')

Note that this kernel requires the output argument explicitly specified, because the type ``Z`` cannot be automatically determined from the input arguments.


Raw argument specifiers
-----------------------

The ElementwiseKernel class does the indexing with broadcasting automatically, which is useful to define most elementwise computations.
On the other hand, we sometimes want to write a kernel with manual indexing for some arguments.
We can tell the ElementwiseKernel class to use manual indexing by adding the ``raw`` keyword preceding the type specifier.

We can use the special variable ``i`` and method ``_ind.size()`` for the manual indexing.
``i`` indicates the index within the loop.
``_ind.size()`` indicates total number of elements to apply the elementwise operation.
Note that it represents the size **after** broadcast operation.

For example, a kernel that adds two vectors with reversing one of them can be written as follows:

.. doctest::

   >>> add_reverse = cp.ElementwiseKernel(
   ...     'T x, raw T y', 'T z',
   ...     'z = x + y[_ind.size() - i - 1]',
   ...     'add_reverse')

(Note that this is an artificial example and you can write such operation just by ``z = x + y[::-1]`` without defining a new kernel).
A raw argument can be used like an array.
The indexing operator ``y[_ind.size() - i - 1]`` involves an indexing computation on ``y``, so ``y`` can be arbitrarily shaped and strode.

Note that raw arguments are not involved in the broadcasting.
If you want to mark all arguments as ``raw``, you must specify the ``size`` argument on invocation, which defines the value of ``_ind.size()``.


Texture memory
--------------
Texture objects (:class:`~cupy.cuda.texture.TextureObject`) can be passed to :class:`~cupy.ElementwiseKernel` with their type marked by a unique type placeholder distinct from any other types used in the same kernel, as its actual datatype is determined when populating the texture memory. The texture coordinates can be computed in the kernel by the per-thread loop index ``i``.


Reduction kernels
-----------------

Reduction kernels can be defined by the :class:`~cupy.ReductionKernel` class.
We can use it by defining four parts of the kernel code:

1. Identity value: This value is used for the initial value of reduction.
2. Mapping expression: It is used for the pre-processing of each element to be reduced.
3. Reduction expression: It is an operator to reduce the multiple mapped values.
   The special variables ``a`` and ``b`` are used for its operands.
4. Post mapping expression: It is used to transform the resulting reduced values.
   The special variable ``a`` is used as its input.
   Output should be written to the output parameter.

ReductionKernel class automatically inserts other code fragments that are required for an efficient and flexible reduction implementation.

For example, L2 norm along specified axes can be written as follows:

.. doctest::

   >>> l2norm_kernel = cp.ReductionKernel(
   ...     'T x',  # input params
   ...     'T y',  # output params
   ...     'x * x',  # map
   ...     'a + b',  # reduce
   ...     'y = sqrt(a)',  # post-reduction map
   ...     '0',  # identity value
   ...     'l2norm'  # kernel name
   ... )
   >>> x = cp.arange(10, dtype=np.float32).reshape(2, 5)
   >>> l2norm_kernel(x, axis=1)
   array([ 5.477226 , 15.9687195], dtype=float32)

.. note::
   ``raw`` specifier is restricted for usages that the axes to be reduced are put at the head of the shape.
   It means, if you want to use ``raw`` specifier for at least one argument, the ``axis`` argument must be ``0`` or a contiguous increasing sequence of integers starting from ``0``, like ``(0, 1)``, ``(0, 1, 2)``, etc.

.. note::
   Texture memory is not yet supported in :class:`~cupy.ReductionKernel`.


Raw kernels
-----------

Raw kernels can be defined by the :class:`~cupy.RawKernel` class.
By using raw kernels, you can define kernels from raw CUDA source.

:class:`~cupy.RawKernel` object allows you to call the kernel with CUDA's ``cuLaunchKernel`` interface.
In other words, you have control over grid size, block size, shared memory size and stream.

.. doctest::

   >>> add_kernel = cp.RawKernel(r'''
   ... extern "C" __global__
   ... void my_add(const float* x1, const float* x2, float* y) {
   ...     int tid = blockDim.x * blockIdx.x + threadIdx.x;
   ...     y[tid] = x1[tid] + x2[tid];
   ... }
   ... ''', 'my_add')
   >>> x1 = cp.arange(25, dtype=cp.float32).reshape(5, 5)
   >>> x2 = cp.arange(25, dtype=cp.float32).reshape(5, 5)
   >>> y = cp.zeros((5, 5), dtype=cp.float32)
   >>> add_kernel((5,), (5,), (x1, x2, y))  # grid, block and arguments
   >>> y
   array([[ 0.,  2.,  4.,  6.,  8.],
          [10., 12., 14., 16., 18.],
          [20., 22., 24., 26., 28.],
          [30., 32., 34., 36., 38.],
          [40., 42., 44., 46., 48.]], dtype=float32)

Raw kernels operating on complex-valued arrays can be created as well:

.. doctest::

   >>> complex_kernel = cp.RawKernel(r'''
   ... #include <cupy/complex.cuh>
   ... extern "C" __global__
   ... void my_func(const complex<float>* x1, const complex<float>* x2,
   ...              complex<float>* y, float a) {
   ...     int tid = blockDim.x * blockIdx.x + threadIdx.x;
   ...     y[tid] = x1[tid] + a * x2[tid];
   ... }
   ... ''', 'my_func')
   >>> x1 = cupy.arange(25, dtype=cupy.complex64).reshape(5, 5)
   >>> x2 = 1j*cupy.arange(25, dtype=cupy.complex64).reshape(5, 5)
   >>> y = cupy.zeros((5, 5), dtype=cupy.complex64)
   >>> complex_kernel((5,), (5,), (x1, x2, y, cupy.float32(2.0)))  # grid, block and arguments
   >>> y
   array([[ 0. +0.j,  1. +2.j,  2. +4.j,  3. +6.j,  4. +8.j],
          [ 5.+10.j,  6.+12.j,  7.+14.j,  8.+16.j,  9.+18.j],
          [10.+20.j, 11.+22.j, 12.+24.j, 13.+26.j, 14.+28.j],
          [15.+30.j, 16.+32.j, 17.+34.j, 18.+36.j, 19.+38.j],
          [20.+40.j, 21.+42.j, 22.+44.j, 23.+46.j, 24.+48.j]],
         dtype=complex64)

Note that while we encourage the usage of ``complex<T>`` types for complex numbers (available by including ``<cupy/complex.cuh>`` as shown above), for CUDA codes already written using functions from ``cuComplex.h`` there is no need to make the conversion yourself: just set the option ``translate_cucomplex=True`` when creating a :class:`~cupy.RawKernel` instance.

The CUDA kernel attributes can be retrieved by either accessing the :attr:`~cupy.RawKernel.attributes` dictionary,
or by accessing the :class:`~cupy.RawKernel` object's attributes directly; the latter can also be used to set certain
attributes:

.. doctest::

   >>> add_kernel = cp.RawKernel(r'''
   ... extern "C" __global__
   ... void my_add(const float* x1, const float* x2, float* y) {
   ...     int tid = blockDim.x * blockIdx.x + threadIdx.x;
   ...     y[tid] = x1[tid] + x2[tid];
   ... }
   ... ''', 'my_add')
   >>> add_kernel.attributes  # doctest: +SKIP
   {'max_threads_per_block': 1024, 'shared_size_bytes': 0, 'const_size_bytes': 0, 'local_size_bytes': 0, 'num_regs': 10, 'ptx_version': 70, 'binary_version': 70, 'cache_mode_ca': 0, 'max_dynamic_shared_size_bytes': 49152, 'preferred_shared_memory_carveout': -1}
   >>> add_kernel.max_dynamic_shared_size_bytes  # doctest: +SKIP
   49152
   >>> add_kernel.max_dynamic_shared_size_bytes = 50000  # set a new value for the attribute  # doctest: +SKIP
   >>> add_kernel.max_dynamic_shared_size_bytes  # doctest: +SKIP
   50000

Dynamical parallelism is supported by :class:`~cupy.RawKernel`. You just need to provide the linking flag (such as ``-dc``) to :class:`~cupy.RawKernel`'s ``options`` argument. The static CUDA device runtime library (``cudadevrt``) is automatically discovered by CuPy. For further detail, see `CUDA Toolkit's documentation`_.

.. _CUDA Toolkit's documentation: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compiling-and-linking

Accessing texture (surface) memory in :class:`~cupy.RawKernel` is supported via CUDA Runtime's Texture (Surface) Object API, see the documentation for :class:`~cupy.cuda.texture.TextureObject` (:class:`~cupy.cuda.texture.SurfaceObject`) as well as CUDA C Programming Guide. For using the Texture Reference API, which is marked as deprecated as of CUDA Toolkit 10.1, see the introduction to :class:`~cupy.RawModule` below.

If your kernel relies on the C++ std library headers such as ``<type_traits>``, it is likely you will encounter compilation errors. In this case, try enabling CuPy's `Jitify <https://github.com/NVIDIA/jitify>`_ support by setting ``jitify=True`` when creating the :class:`~cupy.RawKernel` instance. It provides basic C++ std support to remedy common errors.

.. note::
    The kernel does not have return values.
    You need to pass both input arrays and output arrays as arguments.

.. note::
    When using ``printf()`` in your CUDA kernel, you may need to synchronize the stream to see the output.
    You can use ``cupy.cuda.Stream.null.synchronize()`` if you are using the default stream.

.. note::
    In all of the examples above, we declare the kernels in an ``extern "C"`` block,
    indicating that the C linkage is used. This is to ensure the kernel names are not
    mangled so that they can be retrived by name.

Kernel arguments
----------------
Python primitive types and NumPy scalars are passed to the kernel by value.
Array arguments (pointer arguments) have to be passed as CuPy ndarrays.
No validation is performed by CuPy for arguments passed to the kernel, including types and number of arguments.

Especially note that when passing a CuPy :class:`~cupy.ndarray`, its ``dtype`` should match with the type of the argument declared in the function signature of the CUDA source code (unless you are casting arrays intentionally). 

As an example, ``cupy.float32`` and ``cupy.uint64`` arrays must be passed to the argument typed as ``float*`` and ``unsigned long long*``, respectively. CuPy does not directly support arrays of non-primitive types such as ``float3``, but nothing prevents you from casting a ``float*`` or ``void*`` to a ``float3*`` in a kernel.

Python primitive types, ``int``, ``float``, ``complex`` and ``bool`` map to ``long long``, ``double``, ``cuDoubleComplex`` and ``bool``, respectively.

NumPy scalars (``numpy.generic``) and NumPy arrays (``numpy.ndarray``) **of size one** 
are passed to the kernel by value.
This means that you can pass by value any base NumPy types such as ``numpy.int8`` or ``numpy.float64``, provided the kernel arguments match in size. You can refer to this table to match CuPy/NumPy dtype and CUDA types:

+-----------------+-----------------------------------------------+------------------+
| CuPy/NumPy type | Corresponding kernel types                    | itemsize (bytes) |
+=================+===============================================+==================+
| bool            | bool                                          | 1                |
+-----------------+-----------------------------------------------+------------------+
| int8            | char, signed char                             | 1                |
+-----------------+-----------------------------------------------+------------------+
| int16           | short, signed short                           | 2                |
+-----------------+-----------------------------------------------+------------------+
| int32           | int, signed int                               | 4                |
+-----------------+-----------------------------------------------+------------------+
| int64           | long long, signed long long                   | 8                |
+-----------------+-----------------------------------------------+------------------+
| uint8           | unsigned char                                 | 1                |
+-----------------+-----------------------------------------------+------------------+
| uint16          | unsigned short                                | 2                |
+-----------------+-----------------------------------------------+------------------+
| uint32          | unsigned int                                  | 4                |
+-----------------+-----------------------------------------------+------------------+
| uint64          | unsigned long long                            | 8                |
+-----------------+-----------------------------------------------+------------------+
| float16         | half                                          | 2                |
+-----------------+-----------------------------------------------+------------------+
| float32         | float                                         | 4                |
+-----------------+-----------------------------------------------+------------------+
| float64         | double                                        | 8                |
+-----------------+-----------------------------------------------+------------------+
| complex64       | float2, cuFloatComplex, complex<float>        | 8                |
+-----------------+-----------------------------------------------+------------------+
| complex128      | double2, cuDoubleComplex, complex<double>     | 16               |
+-----------------+-----------------------------------------------+------------------+

The CUDA standard guarantees that the size of fundamental types on the host and device always match.
The itemsize of ``size_t``, ``ptrdiff_t``, ``intptr_t``, ``uintptr_t``, 
``long``, ``signed long`` and ``unsigned long`` are however platform dependent. 
To pass any CUDA vector builtins such as ``float3`` or any other user defined structure 
as kernel arguments (provided it matches the device-side kernel parameter type), see :ref:`custom_user_structs` below.

.. _custom_user_structs:

Custom user types
-----------------

It is possible to use custom types (composite types such as structures and structures of structures)
as kernel arguments by defining a custom NumPy dtype.
When doing this, it is your responsibility to match host and device structure memory layout.
The CUDA standard guarantees that the size of fundamental types on the host and device always match.
It may however impose device alignment requirements on composite types.
This means that for composite types the struct member offsets may be different from what you might expect.

When a kernel argument is passed by value, the CUDA driver will copy exactly ``sizeof(param_type)`` bytes starting from the beginning of the NumPy object data pointer, where ``param_type`` is the parameter type in your kernel. 
You have to match ``param_type``'s memory layout (ex: size, alignment and struct padding/packing) 
by defining a corresponding `NumPy dtype <https://numpy.org/doc/stable/reference/arrays.dtypes.html>`_.

For builtin CUDA vector types such as ``int2`` and ``double4`` and other packed structures with 
named members you can directly define such NumPy dtypes as the following:

.. doctest::

    >>> import numpy as np
    >>> names = ['x', 'y', 'z']
    >>> types = [np.float32]*3
    >>> float3 = np.dtype({'names': names, 'formats': types})
    >>> arg = np.random.rand(3).astype(np.float32).view(float3)
    >>> print(arg)  # doctest: +SKIP
    [(0.9940819, 0.62873816, 0.8953669)]
    >>> arg['x'] = 42.0
    >>> print(arg)  # doctest: +SKIP
    [(42., 0.62873816, 0.8953669)]

Here ``arg`` can be used directly as a kernel argument.
When there is no need to name fields you may prefer this syntax to define packed structures such as 
vectors or matrices:

.. doctest::

    >>> import numpy as np
    >>> float5x5 = np.dtype({'names': ['dummy'], 'formats': [(np.float32,(5,5))]}) 
    >>> arg = np.random.rand(25).astype(np.float32).view(float5x5)
    >>> print(arg.itemsize)
    100

Here ``arg`` represents a 100-byte scalar (i.e. a NumPy array of size 1)
that can be passed by value to any kernel.
Kernel parameters are passed by value in a dedicated 4kB memory bank which has its own cache with broadcast.
Upper bound for total kernel parameters size is thus 4kB
(see `this link <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#function-parameters>`_).
It may be important to note that this dedicated memory bank is not shared with the device ``__constant__`` memory space.

For now, CuPy offers no helper routines to create user defined composite types. 
Such composite types can however be built recursively using NumPy dtype `offsets` and `itemsize` capabilities,
see `cupy/examples/custum_struct <https://github.com/cupy/cupy/tree/master/examples/custom_struct>`_ for examples of advanced usage.

.. warning::
    You cannot directly pass static arrays as kernel arguments with the ``type arg[N]`` syntax where N is a compile time constant. The signature of ``__global__ void kernel(float arg[5])`` is seen as ``__global__ void kernel(float* arg)`` by the compiler. If you want to pass five floats to the kernel by value you need to define a custom structure ``struct float5 { float val[5]; };`` and modify the kernel signature to ``__global__ void kernel(float5 arg)``.


Raw modules
-----------

For dealing a large raw CUDA source or loading an existing CUDA binary, the :class:`~cupy.RawModule` class can be more handy. It can be initialized either by a CUDA source code, or by a path to the CUDA binary. It accepts most of the arguments as in :class:`~cupy.RawKernel`. The needed kernels can then be retrieved by calling the :meth:`~cupy.RawModule.get_function` method, which returns a :class:`~cupy.RawKernel` instance that can be invoked as discussed above.

.. doctest::

    >>> loaded_from_source = r'''
    ... extern "C"{
    ...
    ... __global__ void test_sum(const float* x1, const float* x2, float* y, \
    ...                          unsigned int N)
    ... {
    ...     unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
    ...     if (tid < N)
    ...     {
    ...         y[tid] = x1[tid] + x2[tid];
    ...     }
    ... }
    ...
    ... __global__ void test_multiply(const float* x1, const float* x2, float* y, \
    ...                               unsigned int N)
    ... {
    ...     unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
    ...     if (tid < N)
    ...     {
    ...         y[tid] = x1[tid] * x2[tid];
    ...     }
    ... }
    ...
    ... }'''
    >>> module = cp.RawModule(code=loaded_from_source)
    >>> ker_sum = module.get_function('test_sum')
    >>> ker_times = module.get_function('test_multiply')
    >>> N = 10
    >>> x1 = cp.arange(N**2, dtype=cp.float32).reshape(N, N)
    >>> x2 = cp.ones((N, N), dtype=cp.float32)
    >>> y = cp.zeros((N, N), dtype=cp.float32)
    >>> ker_sum((N,), (N,), (x1, x2, y, N**2))   # y = x1 + x2
    >>> assert cp.allclose(y, x1 + x2)
    >>> ker_times((N,), (N,), (x1, x2, y, N**2)) # y = x1 * x2
    >>> assert cp.allclose(y, x1 * x2)

The instruction above for using complex numbers in :class:`~cupy.RawKernel` also applies to :class:`~cupy.RawModule`.

For CUDA kernels that need to access global symbols, such as constant memory, the :meth:`~cupy.RawModule.get_global` method can be used, see its documentation for further detail.

CuPy also supports the Texture Reference API. A handle to the texture reference in a module can be retrieved by name via :meth:`~cupy.RawModule.get_texref`. Then, you need to pass it to :class:`~cupy.cuda.texture.TextureReference`, along with a resource descriptor and texture descriptor, for binding the reference to the array. (The interface of :class:`~cupy.cuda.texture.TextureReference` is meant to mimic that of :class:`~cupy.cuda.texture.TextureObject` to help users make transition to the latter, since as of CUDA Toolkit 10.1 the former is marked as deprecated.)

To support C++ template kernels, :class:`~cupy.RawModule` additionally provide a ``name_expressions`` argument. A list of template specializations should be provided, so that the corresponding kernels can be generated and retrieved by type:

.. doctest::

    >>> code = r'''
    ... template<typename T>
    ... __global__ void fx3(T* arr, int N) {
    ...     unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    ...     if (tid < N) {
    ...         arr[tid] = arr[tid] * 3;
    ...     }
    ... }
    ... '''
    >>>
    >>> name_exp = ['fx3<float>', 'fx3<double>']
    >>> mod = cp.RawModule(code=code, options=('-std=c++11',),
    ...     name_expressions=name_exp)
    >>> ker_float = mod.get_function(name_exp[0])  # compilation happens here
    >>> N=10
    >>> a = cp.arange(N, dtype=cp.float32)
    >>> ker_float((1,), (N,), (a, N))
    >>> a
    array([ 0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.], dtype=float32)
    >>> ker_double = mod.get_function(name_exp[1])
    >>> a = cp.arange(N, dtype=cp.float64)
    >>> ker_double((1,), (N,), (a, N))
    >>> a
    array([ 0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.])

.. note::

    The name expressions used to both initialize a :class:`~cupy.RawModule` instance and retrieve the kernels are
    the original (*un-mangled*) kernel names with all template parameters unambiguously specified. The name mangling
    and demangling are handled under the hood so that users do not need to worry about it.

.. _kernel_fusion:

Kernel fusion
--------------------

:func:`cupy.fuse` is a decorator that fuses functions.  This decorator can be used to define an elementwise or reduction kernel more easily than :class:`~cupy.ElementwiseKernel` or :class:`~cupy.ReductionKernel`.

By using this decorator, we can define the ``squared_diff`` kernel as follows:

.. doctest::

   >>> @cp.fuse()
   ... def squared_diff(x, y):
   ...     return (x - y) * (x - y)

The above kernel can be called on either scalars, NumPy arrays or CuPy arrays likes the original function.

.. doctest::

   >>> x_cp = cp.arange(10)
   >>> y_cp = cp.arange(10)[::-1]
   >>> squared_diff(x_cp, y_cp)
   array([81, 49, 25,  9,  1,  1,  9, 25, 49, 81])
   >>> x_np = np.arange(10)
   >>> y_np = np.arange(10)[::-1]
   >>> squared_diff(x_np, y_np)
   array([81, 49, 25,  9,  1,  1,  9, 25, 49, 81])

At the first function call, the fused function analyzes the original function based on the abstracted information of arguments (e.g. their dtypes and ndims) and creates and caches an actual CUDA kernel.  From the second function call with the same input types, the fused function calls the previously cached kernel, so it is highly recommended to reuse the same decorated functions instead of decorating local functions that are defined multiple times.

:func:`cupy.fuse` also supports simple reduction kernel.

.. doctest::

   >>> @cp.fuse()
   ... def sum_of_products(x, y):
   ...     return cp.sum(x * y, axis = -1)

You can specify the kernel name by using the ``kernel_name`` keyword argument as follows:

.. doctest::

   >>> @cp.fuse(kernel_name='squared_diff')
   ... def squared_diff(x, y):
   ...     return (x - y) * (x - y)

.. note::
   Currently, :func:`cupy.fuse` can fuse only simple elementwise and reduction operations.  Most other routines (e.g. :func:`cupy.matmul`, :func:`cupy.reshape`) are not supported.

.. _jit_kernel_definition:

JIT kernel definition
---------------------

The :class:`cupyx.jit.rawkernel` decorator can create raw CUDA kernels from Python functions.

In this section, a Python function wrapped with the decorator is called a *target function*.

A target function consists of elementary scalar operations, and users have to manage how to parallelize them. CuPy's array operations which automatically parallelize operations (e.g., :func:`~cupy.add`, :func:`~cupy.sum`) are not supported. If a custom kernel based on such array functions is desired, please refer to the :ref:`kernel_fusion` section.

Basic Usage
^^^^^^^^^^^

Here is a short example for how to write a :class:`cupyx.jit.rawkernel` to copy the values from ``x`` to ``y`` using a grid-stride loop:

.. doctest::

   >>> from cupyx import jit
   >>>
   >>> @jit.rawkernel()
   ... def elementwise_copy(x, y, size):
   ...     tid = jit.blockIdx.x * jit.blockDim.x + jit.threadIdx.x
   ...     ntid = jit.gridDim.x * jit.blockDim.x
   ...     for i in range(tid, size, ntid):
   ...         y[i] = x[i]

   >>> size = cupy.uint32(2 ** 22)
   >>> x = cupy.random.normal(size=(size,), dtype=cupy.float32)
   >>> y = cupy.empty((size,), dtype=cupy.float32)

   >>> elementwise_copy((128,), (1024,), (x, y, size))  # RawKernel style
   >>> assert (x == y).all()

   >>> elementwise_copy[128, 1024](x, y, size)  #  Numba style
   >>> assert (x == y).all()

The above two kinds of styles to launch the kernel are supported, see the documentation of :class:`cupyx.jit._interface._JitRawKernel` for details.

The compilation will be deferred until the first function call. CuPy's JIT compiler infers the types of arguments at the call time, and will cache the compiled kernels for speeding up any subsequent calls.

See :doc:`../reference/kernel` for a full list of API.

Basic Design
^^^^^^^^^^^^

CuPy's JIT compiler generates CUDA code via Python AST. We decided not to use Python bytecode to analyze the target function to avoid perforamance degradation. The CUDA source code generated from the Python bytecode will not effectively optimized by CUDA compiler, because for-loops and other control statements of the target function are fully transformed to jump instruction when converting the target function to bytecode.

Typing rule
^^^^^^^^^^^

The types of local variables are inferred at the first assignment in the function. The first assignment must be done at the top-level of the function; in other words, it must *not* be in ``if``/``else`` bodies or ``for``-loops.

Limitations
^^^^^^^^^^^

JIT does not work inside Python's interactive interpreter (REPL) as the compiler needs to get the source code of the target function.
Fast Fourier Transform with CuPy
================================

CuPy covers the full Fast Fourier Transform (FFT) functionalities provided in NumPy (:mod:`cupy.fft`) and a
subset in SciPy (:mod:`cupyx.scipy.fft`). In addition to those high-level APIs that can be used
as is, CuPy provides additional features to

1. access advanced routines that `cuFFT`_ offers for NVIDIA GPUs,
2. control better the performance and behavior of the FFT routines.

Some of these features are *experimental* (subject to change, deprecation, or removal, see :doc:`./compatibility`)
or may be absent in `hipFFT`_/`rocFFT`_ targeting AMD GPUs.

.. _cuFFT: https://docs.nvidia.com/cuda/cufft/index.html
.. _hipFFT: https://hipfft.readthedocs.io/en/latest/
.. _rocFFT: https://rocfft.readthedocs.io/en/latest/


.. _scipy_fft_backend:

SciPy FFT backend
-----------------

Since SciPy v1.4 a backend mechanism is provided so that users can register different FFT backends and use SciPy's API to perform the actual transform
with the target backend, such as CuPy's :mod:`cupyx.scipy.fft` module. For a one-time only usage, a context manager :func:`scipy.fft.set_backend` can be used:

.. code-block:: python

    import cupy as cp
    import cupyx.scipy.fft as cufft
    import scipy.fft

    a = cp.random.random(100).astype(cp.complex64)
    with scipy.fft.set_backend(cufft):
        b = scipy.fft.fft(a)  # equivalent to cufft.fft(a)

However, such usage can be tedious. Alternatively, users can register a backend through :func:`scipy.fft.register_backend` or :func:`scipy.fft.set_global_backend`
to avoid using context managers:

.. code-block:: python

    import cupy as cp
    import cupyx.scipy.fft as cufft
    import scipy.fft
    scipy.fft.set_global_backend(cufft)

    a = cp.random.random(100).astype(cp.complex64)
    b = scipy.fft.fft(a)  # equivalent to cufft.fft(a)

.. note::

    Please refer to `SciPy FFT documentation`_ for further information.

.. note::
    To use the backend together with an explicit ``plan`` argument requires SciPy version 1.5.0 or higher.
    See below for how to create FFT plans.

.. _SciPy FFT documentation: https://docs.scipy.org/doc/scipy/reference/fft.html#backend-control


User-managed FFT plans
----------------------

For performance reasons, users may wish to create, reuse, and manage the FFT plans themselves. CuPy provides a high-level *experimental* API :func:`~cupyx.scipy.fftpack.get_fft_plan` for this need. Users specify the transform to be performed as they would with most of the high-level FFT APIs, and a plan will be generated based on the input.

.. code-block:: python

    import cupy as cp
    from cupyx.scipy.fft import get_fft_plan

    a = cp.random.random((4, 64, 64)).astype(cp.complex64)
    plan = get_fft_plan(a, axes=(1, 2), value_type='C2C')  # for batched, C2C, 2D transform

The returned plan can be used either explicitly as an argument with the :mod:`cupyx.scipy.fft` APIs:

.. code-block:: python

    import cupyx.scipy.fft

    # the rest of the arguments must match those used when generating the plan
    out = cupyx.scipy.fft.fft2(a, axes=(1, 2), plan=plan)

or as a context manager for the :mod:`cupy.fft` APIs:

.. code-block:: python

    with plan:
        # the arguments must match those used when generating the plan
        out = cp.fft.fft2(a, axes=(1, 2))


.. _fft_plan_cache:

FFT plan cache
--------------

However, there are occasions when users may *not* want to manage the FFT plans by themselves. Moreover, plans could also be reused internally in CuPy's routines, to which user-managed plans would not be applicable. Therefore, starting CuPy v8 we provide a built-in plan cache, enabled by default. The plan cache is done on a *per device, per thread* basis, and can be retrieved by the :func:`~cupy.fft.config.get_plan_cache` API.

.. code-block:: python

    >>> import cupy as cp
    >>>
    >>> cache = cp.fft.config.get_plan_cache()
    >>> cache.show_info()
    ------------------- cuFFT plan cache (device 0) -------------------
    cache enabled? True
    current / max size   : 0 / 16 (counts)
    current / max memsize: 0 / (unlimited) (bytes)
    hits / misses: 0 / 0 (counts)
    
    cached plans (most recently used first):
    
    >>> # perform a transform, which would generate a plan and cache it
    >>> a = cp.random.random((4, 64, 64))
    >>> out = cp.fft.fftn(a, axes=(1, 2))
    >>> cache.show_info()  # hit = 0
    ------------------- cuFFT plan cache (device 0) -------------------
    cache enabled? True
    current / max size   : 1 / 16 (counts)
    current / max memsize: 262144 / (unlimited) (bytes)
    hits / misses: 0 / 1 (counts)
    
    cached plans (most recently used first):
    key: ((64, 64), (64, 64), 1, 4096, (64, 64), 1, 4096, 105, 4, 'C', 2, None), plan type: PlanNd, memory usage: 262144
    
    >>> # perform the same transform again, the plan is looked up from cache and reused
    >>> out = cp.fft.fftn(a, axes=(1, 2))
    >>> cache.show_info()  # hit = 1
    ------------------- cuFFT plan cache (device 0) -------------------
    cache enabled? True
    current / max size   : 1 / 16 (counts)
    current / max memsize: 262144 / (unlimited) (bytes)
    hits / misses: 1 / 1 (counts)
    
    cached plans (most recently used first):
    key: ((64, 64), (64, 64), 1, 4096, (64, 64), 1, 4096, 105, 4, 'C', 2, None), plan type: PlanNd, memory usage: 262144
    
    >>> # clear the cache
    >>> cache.clear()
    >>> cp.fft.config.show_plan_cache_info()  # = cache.show_info(), for all devices
    =============== cuFFT plan cache info (all devices) ===============
    ------------------- cuFFT plan cache (device 0) -------------------
    cache enabled? True
    current / max size   : 0 / 16 (counts)
    current / max memsize: 0 / (unlimited) (bytes)
    hits / misses: 0 / 0 (counts)
    
    cached plans (most recently used first):
    

The returned :class:`~cupy.fft._cache.PlanCache` object has other methods for finer control, such as setting the cache size (either by counts or by memory usage). If the size is set to 0, the cache is disabled. Please refer to its documentation for more detail.

.. note::

    As shown above each FFT plan has an associated working area allocated. If an out-of-memory error happens, one may want to inspect, clear, or limit the plan cache.

.. note::

    The plans returned by :func:`~cupyx.scipy.fftpack.get_fft_plan` are not cached.


FFT callbacks
-------------

`cuFFT`_ provides FFT callbacks for merging pre- and/or post- processing kernels with the FFT routines so as to reduce the access to global memory.
This capability is supported *experimentally* by CuPy. Users need to supply custom load and/or store kernels as strings, and set up a context manager
via :func:`~cupy.fft.config.set_cufft_callbacks`. Note that the load (store) kernel pointer has to be named as ``d_loadCallbackPtr`` (``d_storeCallbackPtr``).

.. code-block:: python

    import cupy as cp

    # a load callback that overwrites the input array to 1
    code = r'''
    __device__ cufftComplex CB_ConvertInputC(
        void *dataIn,
        size_t offset,
        void *callerInfo,
        void *sharedPtr)
    {
        cufftComplex x;
        x.x = 1.;
        x.y = 0.;
        return x;
    }
    __device__ cufftCallbackLoadC d_loadCallbackPtr = CB_ConvertInputC;
    '''

    a = cp.random.random((64, 128, 128)).astype(cp.complex64)

    # this fftn call uses callback
    with cp.fft.config.set_cufft_callbacks(cb_load=code):
        b = cp.fft.fftn(a, axes=(1,2))

    # this does not use
    c = cp.fft.fftn(cp.ones(shape=a.shape, dtype=cp.complex64), axes=(1,2))

    # result agrees
    assert cp.allclose(b, c)

    # "static" plans are also cached, but are distinct from their no-callback counterparts
    cp.fft.config.get_plan_cache().show_info()


.. note::

    Internally, this feature requires recompiling a Python module *for each distinct pair* of load and store kernels. Therefore, the first invocation will be very slow, and this cost is amortized if the callbacks can be reused in the subsequent calculations. The compiled modules are cached on disk, with a default position ``$HOME/.cupy/callback_cache`` that can be changed by the environment variable ``CUPY_CACHE_DIR``.


Multi-GPU FFT
-------------

CuPy currently provides two kinds of *experimental* support for multi-GPU FFT.

.. warning::

    Using multiple GPUs to perform FFT is not guaranteed to be more performant. The rule of thumb is if the transform fits in 1 GPU, you should avoid using multiple.

The first kind of support is with the high-level :func:`~cupy.fft.fft` and :func:`~cupy.fft.ifft` APIs, which requires the input array to reside on one of the participating GPUs. The multi-GPU calculation is done under the hood, and by the end of the calculation the result again resides on the device where it started. Currently only 1D complex-to-complex (C2C) transform is supported; complex-to-real (C2R) or real-to-complex (R2C) transforms (such as :func:`~cupy.fft.rfft` and friends) are not. The transform can be either batched (batch size > 1) or not (batch size = 1).

.. code-block:: python

    import cupy as cp

    cp.fft.config.use_multi_gpus = True
    cp.fft.config.set_cufft_gpus([0, 1])  # use GPU 0 & 1

    shape = (64, 64)  # batch size = 64
    dtype = cp.complex64
    a = cp.random.random(shape).astype(dtype)  # reside on GPU 0

    b = cp.fft.fft(a)  # computed on GPU 0 & 1, reside on GPU 0

If you need to perform 2D/3D transforms (ex: :func:`~cupy.fft.fftn`) instead of 1D (ex: :func:`~cupy.fft.fft`), it would likely still work, but in this particular use case it loops over the transformed axes under the hood (which is exactly what is done in NumPy too), which could lead to suboptimal performance.

The second kind of usage is to use the low-level, *private* CuPy APIs. You need to construct a :class:`~cupy.cuda.cufft.Plan1d` object and use it as if you are programming in C/C++ with `cuFFT`_. Using this approach, your input array can reside on the host as a :class:`numpy.ndarray` so that its size can be much larger than what a single GPU can accommodate, which is one of the main reasons to run multi-GPU FFT.

.. code-block:: python

    import numpy as np
    import cupy as cp

    # no need to touch cp.fft.config, as we are using low-level API

    shape = (64, 64)
    dtype = np.complex64
    a = np.random.random(shape).astype(dtype)  # reside on CPU

    if len(shape) == 1:
        batch = 1
        nx = shape[0]
    elif len(shape) == 2:
        batch = shape[0]
        nx = shape[1]

    # compute via cuFFT
    cufft_type = cp.cuda.cufft.CUFFT_C2C  # single-precision c2c
    plan = cp.cuda.cufft.Plan1d(nx, cufft_type, batch, devices=[0,1])
    out_cp = np.empty_like(a)  # output on CPU
    plan.fft(a, out_cp, cufft.CUFFT_FORWARD)

    out_np = numpy.fft.fft(a)  # use NumPy's fft
    # np.fft.fft alway returns np.complex128
    if dtype is numpy.complex64:
        out_np = out_np.astype(dtype)

    # check result
    assert np.allclose(out_cp, out_np, rtol=1e-4, atol=1e-7)

For this use case, please consult the `cuFFT`_ documentation on multi-GPU transform for further detail.

.. note::

    The multi-GPU plans are cached if auto-generated via the high-level APIs, but not if manually generated via the low-level APIs.


Half-precision FFT
------------------

`cuFFT`_ provides ``cufftXtMakePlanMany`` and ``cufftXtExec`` routines to support a wide range of FFT needs, including 64-bit indexing and half-precision FFT. CuPy provides an *experimental* support for this capability via the new (though *private*) :class:`~cupy.cuda.cufft.XtPlanNd` API. For half-precision FFT, on supported hardware it can be twice as fast than its single-precision counterpart. NumPy does not yet provide the necessary infrastructure for half-precision complex numbers (i.e., ``numpy.complex32``), though, so the steps for this feature is currently a bit more involved than common cases.

.. code-block:: python

    import cupy as cp
    import numpy as np


    shape = (1024, 256, 256)  # input array shape
    idtype = odtype = edtype = 'E'  # = numpy.complex32 in the future

    # store the input/output arrays as fp16 arrays twice as long, as complex32 is not yet available
    a = cp.random.random((shape[0], shape[1], 2*shape[2])).astype(cp.float16)
    out = cp.empty_like(a)

    # FFT with cuFFT
    plan = cp.cuda.cufft.XtPlanNd(shape[1:],
                                  shape[1:], 1, shape[1]*shape[2], idtype,
                                  shape[1:], 1, shape[1]*shape[2], odtype,
                                  shape[0], edtype,
                                  order='C', last_axis=-1, last_size=None)

    plan.fft(a, out, cp.cuda.cufft.CUFFT_FORWARD)

    # FFT with NumPy
    a_np = cp.asnumpy(a).astype(np.float32)  # upcast
    a_np = a_np.view(np.complex64)
    out_np = np.fft.fftn(a_np, axes=(-2,-1))
    out_np = np.ascontiguousarray(out_np).astype(np.complex64)  # downcast
    out_np = out_np.view(np.float32)
    out_np = out_np.astype(np.float16)

    # don't worry about accruacy for now, as we probably lost a lot during casting
    print('ok' if cp.mean(cp.abs(out - cp.asarray(out_np))) < 0.1 else 'not ok')

The 64-bit indexing support for all high-level FFT APIs is planned for a future CuPy release.
User Guide
==========

This user guide provides an overview of CuPy and explains its important features; details are found in :ref:`CuPy API Reference <cupy_reference>`.

.. toctree::
   :maxdepth: 1

   basic
   kernel
   cuda_api
   fft
   memory
   performance
   interoperability
   difference
   compatibility
Accessing CUDA Functionalities
==============================

.. _cuda_stream_event:

Streams and Events
------------------

In this section we discuss basic usages for CUDA streams and events. For the API reference please see
:ref:`stream_event_api`. For their roles in the CUDA programming model, please refer to `CUDA Programming Guide`_.

CuPy provides high-level Python APIs :class:`~cupy.cuda.Stream` and :class:`~cupy.cuda.Event` for creating
streams and events, respectively. Data copies and kernel launches are enqueued onto the :ref:`current_stream`,
which can be queried via :func:`~cupy.cuda.get_current_stream` and changed either by setting up a context
manager:

.. doctest::

    >>> import numpy as np
    >>>
    >>> a_np = np.arange(10)
    >>> s = cp.cuda.Stream()
    >>> with s:
    ...     a_cp = cp.asarray(a_np)  # H2D transfer on stream s
    ...     b_cp = cp.sum(a_cp)      # kernel launched on stream s
    ...     assert s == cp.cuda.get_current_stream()
    ...
    >>> # fall back to the previous stream in use (here the default stream)
    >>> # when going out of the scope of s

or by using the :meth:`~cupy.cuda.Stream.use` method:

.. doctest::

    >>> s = cp.cuda.Stream()
    >>> s.use()  # any subsequent operations are done on steam s  # doctest: +ELLIPSIS
    <Stream ... (device ...)>
    >>> b_np = cp.asnumpy(b_cp)
    >>> assert s == cp.cuda.get_current_stream()
    >>> cp.cuda.Stream.null.use()  # fall back to the default (null) stream
    <Stream 0 (device -1)>
    >>> assert cp.cuda.Stream.null == cp.cuda.get_current_stream()

Events can be created either manually or through the :meth:`~cupy.cuda.Stream.record` method.
:class:`~cupy.cuda.Event` objects can be used for timing GPU activities (via :func:`~cupy.cuda.get_elapsed_time`)
or setting up inter-stream dependencies:

.. doctest::

    >>> e1 = cp.cuda.Event()
    >>> e1.record()
    >>> a_cp = b_cp * a_cp + 8
    >>> e2 = cp.cuda.get_current_stream().record()
    >>>
    >>> # set up a stream order
    >>> s2 = cp.cuda.Stream()
    >>> s2.wait_event(e2)
    >>> with s2:
    ...     # the a_cp is guaranteed updated when this copy (on s2) starts
    ...     a_np = cp.asnumpy(a_cp)
    >>>
    >>> # timing
    >>> e2.synchronize()
    >>> t = cp.cuda.get_elapsed_time(e1, e2)  # only include the compute time, not the copy time

Just like the :class:`~cupy.cuda.Device` objects, :class:`~cupy.cuda.Stream` and :class:`~cupy.cuda.Event`
objects can also be used for synchronization.

.. note::

    In CuPy, the :class:`~cupy.cuda.Stream` objects are managed on the per thread, per device basis.

.. note::

    On NVIDIA GPUs, there are two stream singleton objects :obj:`~cupy.cuda.Stream.null` and
    :obj:`~cupy.cuda.Stream.ptds`, referred to as the *legacy* default stream and the *per-thread* default
    stream, respectively. CuPy uses the former as default when no user-defined stream is in use. To
    change this behavior, set the environment variable ``CUPY_CUDA_PER_THREAD_DEFAULT_STREAM`` to 1,
    see :ref:`environment`. This is not applicable to AMD GPUs.

.. _CUDA Programming Guide: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html

To interoperate with streams created in other Python libraries, CuPy provides the :class:`~cupy.cuda.ExternalStream`
API to wrap an existing stream pointer (given as a Python `int`). In this case, the stream lifetime is not managed
by CuPy. In addition, you need to make sure the :class:`~cupy.cuda.ExternalStream` object is used on the device
where the stream was created, either manually or by explicitly setting the optional `device_id` argument. But the
created :class:`~cupy.cuda.ExternalStream` object can otherwise be used like a :class:`~cupy.cuda.Stream` object.

CUDA Driver and Runtime API
---------------------------

Under construction. Please see :ref:`runtime_api` for the API reference.
API Compatibility Policy
========================

This document expresses the design policy on compatibilities of CuPy APIs.
Development team should obey this policy on deciding to add, extend, and change APIs and their behaviors.

This document is written for both users and developers.
Users can decide the level of dependencies on CuPy’s implementations in their codes based on this document.
Developers should read through this document before creating pull requests that contain changes on the interface.
Note that this document may contain ambiguities on the level of supported compatibilities.


Versioning and Backward Compatibilities
---------------------------------------

The updates of CuPy are classified into three levels: major, minor, and revision.
These types have distinct levels of backward compatibilities.

- **Major update** contains disruptive changes that break the backward compatibility.
- **Minor update** contains additions and extensions to the APIs that keep the backward compatibility supported.
- **Revision update** contains improvements on the API implementations without changing any API specifications.

Note that we do not support full backward compatibility, which is almost infeasible for Python-based APIs, since there is no way to completely hide the implementation details.


Processes to Break Backward Compatibilities
-------------------------------------------

Deprecation, Dropping, and Its Preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any APIs may be *deprecated* at some minor updates.
In such a case, the deprecation note is added to the API documentation, and the API implementation is changed to fire a deprecation warning (if possible).
There should be another way to reimplement the same functionality previously written using the deprecated APIs.

Any APIs may be marked as *to be dropped in the future*.
In such a case, the dropping is stated in the documentation with the major version number on which the API is planned to be dropped, and the API implementation is changed to fire a future warning (if possible).

The actual dropping should be done through the following steps:

- Make the API deprecated.
  At this point, users should not use the deprecated API in their new application codes.
- After that, mark the API as *to be dropped in the future*.
  It must be done in the minor update different from that of the deprecation.
- At the major version announced in the above update, drop the API.

Consequently, it takes at least two minor versions to drop any APIs after the first deprecation.

API Changes and Its Preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any APIs may be marked as *to be changed in the future* for changes without backward compatibility.
In such a case, the change is stated in the documentation with the version number on which the API is planned to be changed, and the API implementation is changed to fire the future warning on the certain usages.

The actual change should be done in the following steps:

- Announce that the API will be changed in the future.
  At this point, the actual version of change need not be accurate.
- After the announcement, mark the API as *to be changed in the future* with version number of planned changes.
  At this point, users should not use the marked API in their new application codes.
- At the major update announced in the above update, change the API.


Supported Backward Compatibility
--------------------------------

This section defines backward compatibilities that minor updates must maintain.

Documented Interface
~~~~~~~~~~~~~~~~~~~~

CuPy has an official API documentation.
Many applications can be written based on the documented features.
We support backward compatibilities of documented features.
In other words, codes only based on the documented features run correctly with minor-/revision- updated versions.

Developers are encouraged to use apparent names for objects of implementation details.
For example, attributes outside of the documented APIs should have one or more underscores at the prefix of their names.

.. _undocumented_behavior:

Undocumented behaviors
~~~~~~~~~~~~~~~~~~~~~~

Behaviors of CuPy implementation not stated in the documentation are undefined.
Undocumented behaviors are not guaranteed to be stable between different minor/revision versions.

Minor update may contain changes to undocumented behaviors.
For example, suppose an API X is added at the minor update.
In the previous version, attempts to use X cause AttributeError.
This behavior is not stated in the documentation, so this is undefined.
Thus, adding the API X in minor version is permissible.

Revision update may also contain changes to undefined behaviors.
Typical example is a bug fix.
Another example is an improvement on implementation, which may change the internal object structures not shown in the documentation.
As a consequence, **even revision updates do not support compatibility of pickling, unless the full layout of pickled objects is clearly documented.**

Documentation Error
~~~~~~~~~~~~~~~~~~~

Compatibility is basically determined based on the documentation, though it sometimes contains errors.
It may make the APIs confusing to assume the documentation always stronger than the implementations.
We therefore may fix the documentation errors in any updates that may break the compatibility in regard to the documentation.

.. note::
   Developers MUST NOT fix the documentation and implementation of the same functionality at the same time in revision updates as "bug fix".
   Such a change completely breaks the backward compatibility.
   If you want to fix the bugs in both sides, first fix the documentation to fit it into the implementation, and start the API changing procedure described above.

Object Attributes and Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Object attributes and properties are sometimes replaced by each other at minor updates.
It does not break the user codes, except for the codes depending on how the attributes and properties are implemented.

Functions and Methods
~~~~~~~~~~~~~~~~~~~~~

Methods may be replaced by callable attributes keeping the compatibility of parameters and return values in minor updates.
It does not break the user codes, except for the codes depending on how the methods and callable attributes are implemented.

Exceptions and Warnings
~~~~~~~~~~~~~~~~~~~~~~~

The specifications of raising exceptions are considered as a part of standard backward compatibilities.
No exception is raised in the future versions with correct usages that the documentation allows, unless the API changing process is completed.

On the other hand, warnings may be added at any minor updates for any APIs.
It means minor updates do not keep backward compatibility of warnings.


Installation Compatibility
--------------------------

The installation process is another concern of compatibilities.
We support environmental compatibilities in the following ways.

- Any changes of dependent libraries that force modifications on the existing environments must be done in major updates.
  Such changes include following cases:

  - dropping supported versions of dependent libraries (e.g. dropping cuDNN v2)
  - adding new mandatory dependencies (e.g. adding h5py to setup_requires)

- Supporting optional packages/libraries may be done in minor updates (e.g. supporting h5py in optional features).

.. note::
   The installation compatibility does not guarantee that all the features of CuPy correctly run on supported environments.
   It may contain bugs that only occurs in certain environments.
   Such bugs should be fixed in some updates.
Basics of CuPy
==============

.. currentmodule:: cupy

In this section, you will learn about the following things:

* Basics of :class:`cupy.ndarray`
* The concept of *current device*
* host-device and device-device array transfer


Basics of cupy.ndarray
----------------------

CuPy is a GPU array backend that implements a subset of NumPy interface.
In the following code, ``cp`` is an abbreviation of ``cupy``, following the standard convention of abbreviating ``numpy`` as ``np``:

.. doctest::

   >>> import numpy as np
   >>> import cupy as cp

The :class:`cupy.ndarray` class is at the core of ``CuPy`` and is a replacement class for ``NumPy``'s :class:`numpy.ndarray`.

.. doctest::

   >>> x_gpu = cp.array([1, 2, 3])

``x_gpu`` above is an instance of :class:`cupy.ndarray`.
As one can see, CuPy's syntax here is identical to that of NumPy.
The main difference between :class:`cupy.ndarray` and :class:`numpy.ndarray` is that
the CuPy arrays are allocated on the *current device*, which we will talk about later.

Most of the array manipulations are also done in the way similar to NumPy.
Take the Euclidean norm (a.k.a L2 norm), for example.
NumPy has :func:`numpy.linalg.norm` function that calculates it on CPU.

.. doctest::

   >>> x_cpu = np.array([1, 2, 3])
   >>> l2_cpu = np.linalg.norm(x_cpu)

Using CuPy, we can perform the same calculations on GPU in a similar way:

.. doctest::

   >>> x_gpu = cp.array([1, 2, 3])
   >>> l2_gpu = cp.linalg.norm(x_gpu)

CuPy implements many functions on :class:`cupy.ndarray` objects.
See the :ref:`reference <cupy_reference>` for the supported subset of NumPy API.
Knowledge of NumPy will help you utilize most of the CuPy features.
We, therefore, recommend you familiarize yourself with the `NumPy documentation <https://numpy.org/doc/stable/index.html>`_.


Current Device
--------------

CuPy has a concept of a *current device*, which is the default GPU device on which
the allocation, manipulation, calculation, etc., of arrays take place.
Suppose ID of the current device is 0.
In such a case, the following code would create an array ``x_on_gpu0`` on GPU 0.

.. doctest::

   >>> x_on_gpu0 = cp.array([1, 2, 3, 4, 5])

To switch to another GPU device, use the :class:`~cupy.cuda.Device` context manager:

.. doctest::

   >>> with cp.cuda.Device(1):
   ...    x_on_gpu1 = cp.array([1, 2, 3, 4, 5])
   >>> x_on_gpu0 = cp.array([1, 2, 3, 4, 5])

All CuPy operations (except for multi-GPU features and device-to-device copy) are performed on the currently active device.

In general, CuPy functions expect that the array is on the same device as the current one.
Passing an array stored on a non-current device may work depending on the hardware configuration but is generally discouraged as it may not be performant.

.. note::
  If the array's device and the current device mismatch, CuPy functions try to establish `peer-to-peer memory access <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#peer-to-peer-memory-access>`_ (P2P) between them so that the current device can directly read the array from another device.
  Note that P2P is available only when the topology permits it.
  If P2P is unavailable, such an attempt will fail with ``ValueError``.

``cupy.ndarray.device`` attribute indicates the device on which the array is allocated.

.. doctest::

   >>> with cp.cuda.Device(1):
   ...    x = cp.array([1, 2, 3, 4, 5])
   >>> x.device
   <CUDA Device 1>

.. note::

   When only one device is available, explicit device switching is not needed.


.. _current_stream:

Current Stream
--------------

Associated with the concept of current devices are *current streams*, which help avoid explicitly passing streams
in every single operation so as to keep the APIs pythonic and user-friendly. In CuPy, all CUDA operations
such as data transfer (see the :ref:`data-transfer-basics` section) and kernel launches are enqueued onto the current stream,
and the queued tasks on the same stream will be executed in serial (but *asynchronously* with respect to the host).

The default current stream in CuPy is CUDA's null stream (i.e., stream 0). It is also known as the *legacy*
default stream, which is unique per device. However, it is possible to change the current stream using the
:class:`cupy.cuda.Stream` API, please see :doc:`cuda_api` for example. The current stream in CuPy can be
retrieved using :func:`cupy.cuda.get_current_stream`.

It is worth noting that CuPy's current stream is managed on a *per thread, per device* basis, meaning that on different
Python threads or different devices the current stream (if not the null stream) can be different.

.. _data-transfer-basics:

Data Transfer
-------------

Move arrays to a device
~~~~~~~~~~~~~~~~~~~~~~~

:func:`cupy.asarray` can be used to move a :class:`numpy.ndarray`, a list, or any object
that can be passed to :func:`numpy.array` to the current device:

.. doctest::

   >>> x_cpu = np.array([1, 2, 3])
   >>> x_gpu = cp.asarray(x_cpu)  # move the data to the current device.

:func:`cupy.asarray` can accept :class:`cupy.ndarray`, which means we can
transfer the array between devices with this function.

.. doctest::

   >>> with cp.cuda.Device(0):
   ...     x_gpu_0 = cp.ndarray([1, 2, 3])  # create an array in GPU 0
   >>> with cp.cuda.Device(1):
   ...     x_gpu_1 = cp.asarray(x_gpu_0)  # move the array to GPU 1

.. note::

   :func:`cupy.asarray` does not copy the input array if possible.
   So, if you put an array of the current device, it returns the input object itself.

   If we do copy the array in this situation, you can use :func:`cupy.array` with `copy=True`.
   Actually :func:`cupy.asarray` is equivalent to `cupy.array(arr, dtype, copy=False)`.

Move array from a device to the host
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Moving a device array to the host can be done by :func:`cupy.asnumpy` as follows:

.. doctest::

   >>> x_gpu = cp.array([1, 2, 3])  # create an array in the current device
   >>> x_cpu = cp.asnumpy(x_gpu)  # move the array to the host.

We can also use :meth:`cupy.ndarray.get()`:

.. doctest::

   >>> x_cpu = x_gpu.get()


Memory management
-----------------

Check :doc:`./memory` for a detailed description of how memory is managed in CuPy
using memory pools.


How to write CPU/GPU agnostic code
----------------------------------

CuPy's compatibility with NumPy makes it possible to write CPU/GPU agnostic code.
For this purpose, CuPy implements the :func:`cupy.get_array_module` function that
returns a reference to :mod:`cupy` if any of its arguments resides on a GPU
and :mod:`numpy` otherwise.
Here is an example of a CPU/GPU agnostic function that computes ``log1p``:

.. doctest::

   >>> # Stable implementation of log(1 + exp(x))
   >>> def softplus(x):
   ...     xp = cp.get_array_module(x)  # 'xp' is a standard usage in the community
   ...     print("Using:", xp.__name__)
   ...     return xp.maximum(0, x) + xp.log1p(xp.exp(-abs(x)))

When you need to manipulate CPU and GPU arrays, an explicit data
transfer may be required to move them to the same location -- either CPU or GPU.
For this purpose, CuPy implements two sister methods called :func:`cupy.asnumpy`  and
:func:`cupy.asarray`. Here is an example that demonstrates the use of both methods:

.. doctest::

   >>> x_cpu = np.array([1, 2, 3])
   >>> y_cpu = np.array([4, 5, 6])
   >>> x_cpu + y_cpu
   array([5, 7, 9])
   >>> x_gpu = cp.asarray(x_cpu)
   >>> x_gpu + y_cpu
   Traceback (most recent call last):
   ...
   TypeError: Unsupported type <class 'numpy.ndarray'>
   >>> cp.asnumpy(x_gpu) + y_cpu
   array([5, 7, 9])
   >>> cp.asnumpy(x_gpu) + cp.asnumpy(y_cpu)
   array([5, 7, 9])
   >>> x_gpu + cp.asarray(y_cpu)
   array([5, 7, 9])
   >>> cp.asarray(x_gpu) + cp.asarray(y_cpu)
   array([5, 7, 9])

The :func:`cupy.asnumpy` method returns a NumPy array (array on the host),
whereas :func:`cupy.asarray` method returns a CuPy array (array on the current device).
Both methods can accept arbitrary input, meaning that they can be applied to any data that
is located on either the host or device and can be converted to an array.
Miscellaneous routines
======================

.. Hint:: `NumPy API Reference: Miscellaneous routines <https://numpy.org/doc/stable/reference/routines.other.html>`_

.. currentmodule:: cupy

Memory ranges
-------------

.. autosummary::
   :toctree: generated/

   shares_memory
   may_share_memory

Utility
-------------

.. autosummary::
   :toctree: generated/

   show_config

Matlab-like Functions
---------------------

.. autosummary::
   :toctree: generated/

   who
Binary operations
=================

.. Hint:: `NumPy API Reference: Binary operations <https://numpy.org/doc/stable/reference/routines.bitwise.html>`_

.. currentmodule:: cupy

Elementwise bit operations
--------------------------

.. autosummary::
   :toctree: generated/

   bitwise_and
   bitwise_or
   bitwise_xor
   invert
   left_shift
   right_shift


Bit packing
-----------

.. autosummary::
   :toctree: generated/

   packbits
   unpackbits


Output formatting
-----------------

.. autosummary::
   :toctree: generated/

   binary_repr
Indexing routines
=================

.. Hint:: `NumPy API Reference: Indexing routines <https://numpy.org/doc/stable/reference/routines.indexing.html>`_

.. currentmodule:: cupy

Generating index arrays
-----------------------

.. autosummary::
   :toctree: generated/

   c_
   r_
   nonzero
   where
   indices
   mask_indices
   tril_indices
   tril_indices_from
   triu_indices
   triu_indices_from
   ix_
   ravel_multi_index
   unravel_index
   diag_indices
   diag_indices_from


Indexing-like operations
------------------------

.. autosummary::
   :toctree: generated/

   take
   take_along_axis
   choose
   compress
   diag
   diagonal
   select
   lib.stride_tricks.as_strided


Inserting data into arrays
--------------------------

.. autosummary::
   :toctree: generated/

   place
   put
   putmask
   fill_diagonal


Iterating over arrays
---------------------

.. autosummary::
   :toctree: generated/

   flatiter
:orphan:


DLPack helper
-------------

.. autosummary::
   :toctree: generated/

   cupy.fromDlpack

Time range
----------

.. autosummary::
   :toctree: generated/

   cupy.prof.TimeRangeDecorator
   cupy.prof.time_range

Timing helper
-------------

.. autosummary::
   :toctree: generated/

   cupyx.time.repeat

Device synchronization detection
--------------------------------

.. warning::

   These APIs are deprecated in CuPy v10 and will be removed in future releases.

.. autosummary::
   :toctree: generated/

   cupyx.allow_synchronize
   cupyx.DeviceSynchronized
Window functions
================

.. Hint:: `NumPy API Reference: Window functions <https://numpy.org/doc/stable/reference/routines.window.html>`_

.. currentmodule:: cupy

Various windows
---------------

.. autosummary::
   :toctree: generated/

   bartlett
   blackman
   hamming
   hanning
   kaiser
.. module:: cupyx.scipy.ndimage

Multidimensional image processing (:mod:`cupyx.scipy.ndimage`)
==============================================================

.. Hint:: `SciPy API Reference: Multidimensional image processing (scipy.ndimage) <https://docs.scipy.org/doc/scipy/reference/ndimage.html>`_


Filters
-------

.. autosummary::
   :toctree: generated/

   convolve
   convolve1d
   correlate
   correlate1d
   gaussian_filter
   gaussian_filter1d
   gaussian_gradient_magnitude
   gaussian_laplace
   generic_filter
   generic_filter1d
   generic_gradient_magnitude
   generic_laplace
   laplace
   maximum_filter
   maximum_filter1d
   median_filter
   minimum_filter
   minimum_filter1d
   percentile_filter
   prewitt
   rank_filter
   sobel
   uniform_filter
   uniform_filter1d


Fourier filters
---------------

.. autosummary::
   :toctree: generated/

   fourier_ellipsoid
   fourier_gaussian
   fourier_shift
   fourier_uniform


Interpolation
-------------

.. autosummary::
   :toctree: generated/

   affine_transform
   map_coordinates
   rotate
   shift
   spline_filter
   spline_filter1d
   zoom


Measurements
------------

.. autosummary::
   :toctree: generated/

   center_of_mass
   extrema
   histogram
   label
   labeled_comprehension
   maximum
   maximum_position
   mean
   median
   minimum
   minimum_position
   standard_deviation
   sum_labels
   variance


Morphology
----------

.. autosummary::
   :toctree: generated/

   binary_closing
   binary_dilation
   binary_erosion
   binary_fill_holes
   binary_hit_or_miss
   binary_opening
   binary_propagation
   black_tophat
   generate_binary_structure
   grey_closing
   grey_dilation
   grey_erosion
   grey_opening
   iterate_structure
   morphological_gradient
   morphological_laplace
   white_tophat


OpenCV mode
-----------
:mod:`cupyx.scipy.ndimage` supports additional mode, ``opencv``.
If it is given, the function performs like `cv2.warpAffine <https://docs.opencv.org/master/da/d54/group__imgproc__transform.html#ga0203d9ee5fcd28d40dbc4a1ea4451983>`_ or `cv2.resize <https://docs.opencv.org/master/da/d54/group__imgproc__transform.html#ga47a974309e9102f5f08231edc7e7529d>`_. Example:


.. code:: python

   import cupyx.scipy.ndimage
   import cupy as cp
   import cv2

   im = cv2.imread('TODO') # pls fill in your image path

   trans_mat = cp.eye(4)
   trans_mat[0][0] = trans_mat[1][1] = 0.5

   smaller_shape = (im.shape[0] // 2, im.shape[1] // 2, 3)
   smaller = cp.zeros(smaller_shape) # preallocate memory for resized image

   cupyx.scipy.ndimage.affine_transform(im, trans_mat, output_shape=smaller_shape,
                                        output=smaller, mode='opencv')

   cv2.imwrite('smaller.jpg', cp.asnumpy(smaller)) # smaller image saved locally

Logic functions
===============

.. Hint:: `NumPy API Reference: Logic functions <https://numpy.org/doc/stable/reference/routines.logic.html>`_

.. currentmodule:: cupy

Truth value testing
-------------------

.. autosummary::
   :toctree: generated/

   all
   any
   union1d


Array contents
--------------

.. autosummary::
   :toctree: generated/

   isfinite
   isinf
   isnan


Array type testing
------------------

.. autosummary::
   :toctree: generated/

   iscomplex
   iscomplexobj
   isfortran
   isreal
   isrealobj
   isscalar


Logic operations
----------------

.. autosummary::
   :toctree: generated/

   logical_and
   logical_or
   logical_not
   logical_xor


Comparison
----------

.. autosummary::
   :toctree: generated/

   allclose
   isclose
   array_equal
   array_equiv
   greater
   greater_equal
   less
   less_equal
   equal
   not_equal
Functional programming
======================

.. Hint:: `NumPy API Reference: Functional programming <https://numpy.org/doc/stable/reference/routines.functional.html>`_

.. currentmodule:: cupy

.. note::

   :class:`cupy.vectorize` applies JIT compiler to the given Python function.
   See :ref:`jit_kernel_definition` for details.

.. autosummary::
   :toctree: generated/

   apply_along_axis
   vectorize
   piecewise
Custom kernels
==============

.. autosummary::
   :toctree: generated/

   cupy.ElementwiseKernel
   cupy.ReductionKernel
   cupy.RawKernel
   cupy.RawModule
   cupy.fuse


JIT kernel definition
---------------------

.. autosummary::
   :toctree: generated/

   cupyx.jit.rawkernel
   cupyx.jit.threadIdx
   cupyx.jit.blockDim
   cupyx.jit.blockIdx
   cupyx.jit.gridDim
   cupyx.jit.grid
   cupyx.jit.gridsize
   cupyx.jit.laneid
   cupyx.jit.warpsize
   cupyx.jit.syncthreads
   cupyx.jit.syncwarp
   cupyx.jit.shfl_sync
   cupyx.jit.shfl_up_sync
   cupyx.jit.shfl_down_sync
   cupyx.jit.shfl_xor_sync
   cupyx.jit.shared_memory
   cupyx.jit.atomic_add
   cupyx.jit.atomic_sub
   cupyx.jit.atomic_exch
   cupyx.jit.atomic_min
   cupyx.jit.atomic_max
   cupyx.jit.atomic_inc
   cupyx.jit.atomic_dec
   cupyx.jit.atomic_cas
   cupyx.jit.atomic_and
   cupyx.jit.atomic_or
   cupyx.jit.atomic_xor
   cupyx.jit._interface._JitRawKernel


Kernel binary memoization
-------------------------

.. autosummary::
   :toctree: generated/

   cupy.memoize
   cupy.clear_memo
Array API Functions
===================

This section is a full list of implemented APIs. For the detailed documentation, see the
`array API specification <https://data-apis.org/array-api/latest/API_specification/index.html>`_.

.. automodule:: cupy.array_api
   :members:
:orphan:

This document has been moved to :doc:`scipy_sparse`.
Comparison Table
================

Here is a list of NumPy / SciPy APIs and its corresponding CuPy implementations.

``-`` in CuPy column denotes that CuPy implementation is not provided yet.
We welcome contributions for these functions.

.. include:: comparison_table.rst.inc
.. module:: cupy.random

Random sampling (:mod:`cupy.random`)
====================================

Differences between :mod:`cupy.random` and :mod:`numpy.random`:

* Most functions under :mod:`cupy.random` support the ``dtype`` option, which do not exist in the corresponding NumPy APIs.
  This option enables generation of float32 values directly without any space overhead.
* :func:`cupy.random.default_rng` uses XORWOW bit generator by default.
* Random states cannot be serialized. See the description below for details.
* CuPy does not guarantee that the same number generator is used across major versions.
  This means that numbers generated by :mod:`cupy.random` by new major version may not be the same as the previous one, even if the same seed and distribution are used.

.. currentmodule:: cupy.random

New Random Generator API
------------------------

.. Hint:: `NumPy API Reference: Random sampling (numpy.random) <https://numpy.org/doc/stable/reference/random/>`_

Random Generator
~~~~~~~~~~~~~~~~

.. Hint:: `NumPy API Reference: Random Generator <https://numpy.org/doc/stable/reference/random/generator.html>`_

.. autosummary::
   :toctree: generated/

   default_rng
   Generator

Bit Generators
~~~~~~~~~~~~~~

.. Hint:: `NumPy API Reference: Bit Generators <https://numpy.org/doc/stable/reference/random/bit_generators/index.html>`_

.. autosummary::
   :toctree: generated/

   BitGenerator

CuPy provides the following bit generator implementations:

.. autosummary::
   :toctree: generated/

   XORWOW
   MRG32k3a
   Philox4x3210

Legacy Random Generation
------------------------

.. Hint::
   * `NumPy API Reference: Legacy Random Generation <https://numpy.org/doc/stable/reference/random/legacy.html>`_
   * `NumPy 1.16 Reference <https://numpy.org/doc/1.16/reference/routines.random.html>`_

.. autosummary::
   :toctree: generated/

   RandomState

Functions in :mod:`cupy.random`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   beta
   binomial
   bytes
   chisquare
   choice
   dirichlet
   exponential
   f
   gamma
   geometric
   gumbel
   hypergeometric
   laplace
   logistic
   lognormal
   logseries
   multinomial
   multivariate_normal
   negative_binomial
   noncentral_chisquare
   noncentral_f
   normal
   pareto
   permutation
   poisson
   power
   rand
   randint
   randn
   random
   random_integers
   random_sample
   ranf
   rayleigh
   sample
   seed
   shuffle
   standard_cauchy
   standard_exponential
   standard_gamma
   standard_normal
   standard_t
   triangular
   uniform
   vonmises
   wald
   weibull
   zipf

CuPy does not provide ``cupy.random.get_state`` nor ``cupy.random.set_state`` at this time.
Use the following CuPy-specific APIs instead.
Note that these functions use :class:`cupy.random.RandomState` instance to represent the internal state, which cannot be serialized.

.. autosummary::
   :toctree: generated/

   get_random_state
   set_random_state
.. module:: cupyx.scipy.stats

Statistical functions (:mod:`cupyx.scipy.stats`)
================================================

.. Hint:: `SciPy API Reference: Statistical functions (scipy.stats) <https://docs.scipy.org/doc/scipy/reference/stats.html>`_


Summary statistics
------------------

.. autosummary::
   :toctree: generated/

   trim_mean
   entropy
----------------
Routines (NumPy)
----------------

The following pages describe NumPy-compatible routines.
These functions cover a subset of
`NumPy routines <https://numpy.org/doc/stable/reference/routines.html>`_.

.. toctree::
   :maxdepth: 2

   creation
   manipulation
   binary
   dtype
   fft
   functional
   indexing
   io
   linalg
   logic
   math
   misc
   pad
   polynomials
   random
   set
   sorting
   statistics
   testing
   window
.. module:: cupyx.scipy.fft

Discrete Fourier transforms (:mod:`cupyx.scipy.fft`)
====================================================

.. Hint:: `SciPy API Reference: Discrete Fourier transforms (scipy.fft) <https://docs.scipy.org/doc/scipy/reference/fft.html>`_

.. seealso:: :doc:`../user_guide/fft`

Fast Fourier Transforms (FFTs)
------------------------------

.. autosummary::
   :toctree: generated/

   fft
   ifft
   fft2
   ifft2
   fftn
   ifftn
   rfft
   irfft
   rfft2
   irfft2
   rfftn
   irfftn
   hfft
   ihfft
   hfft2
   ihfft2
   hfftn
   ihfftn

Discrete Cosine and Sine Transforms (DST and DCT)
-------------------------------------------------

.. autosummary::
   :toctree: generated/

   dct
   idct
   dctn
   idctn
   dst
   idst
   dstn
   idstn

Helper functions
----------------

.. autosummary::
   :toctree: generated/

   fftshift
   ifftshift
   fftfreq
   rfftfreq
   next_fast_len


Code compatibility features
---------------------------
1. As with other FFT modules in CuPy, FFT functions in this module can take advantage of an existing cuFFT plan (returned by :func:`~cupyx.scipy.fftpack.get_fft_plan`) to accelerate the computation. The plan can be either passed in explicitly via the keyword-only ``plan`` argument or used as a context manager. One exception to this are the DCT and DST transforms, which do not currently support a plan argument.

2. The boolean switch ``cupy.fft.config.enable_nd_planning`` also affects the FFT functions in this module, see :doc:`./fft`. This switch is neglected when planning manually using :func:`~cupyx.scipy.fftpack.get_fft_plan`.

3. Like in ``scipy.fft``, all FFT functions in this module have an optional argument ``overwrite_x`` (default is ``False``), which has the same semantics as in ``scipy.fft``: when it is set to ``True``, the input array ``x`` *can* (not *will*) be overwritten arbitrarily. For this reason, when an in-place FFT is desired, the user should always reassign the input in the following manner: ``x = cupyx.scipy.fftpack.fft(x, ..., overwrite_x=True, ...)``.

4. The ``cupyx.scipy.fft`` module can also be used as a backend for ``scipy.fft`` e.g. by installing with ``scipy.fft.set_backend(cupyx.scipy.fft)``. This can allow ``scipy.fft`` to work with both ``numpy`` and ``cupy`` arrays. For more information, see :ref:`scipy_fft_backend`.

5. The boolean switch :data:`cupy.fft.config.use_multi_gpus` also affects the FFT functions in this module, see :doc:`./fft`. Moreover, this switch is *honored* when planning manually using :func:`~cupyx.scipy.fftpack.get_fft_plan`.

6. Both type II and III DCT and DST transforms are implemented. Type I and IV transforms are currently unavailable.
CuPy-specific functions
=======================

CuPy-specific functions are placed under ``cupyx`` namespace.

.. TODO(kmaehashi): use module:: cupyx
.. autosummary::
   :toctree: generated/

   cupyx.rsqrt
   cupyx.scatter_add
   cupyx.scatter_max
   cupyx.scatter_min
   cupyx.empty_pinned
   cupyx.empty_like_pinned
   cupyx.zeros_pinned
   cupyx.zeros_like_pinned

Profiling utilities
-------------------

.. autosummary::
   :toctree: generated/

   cupyx.profiler.benchmark
   cupyx.profiler.time_range
   cupyx.profiler.profile

DLPack utilities
----------------

Below are helper functions for creating a :class:`cupy.ndarray` from either a DLPack tensor
or any object supporting the DLPack data exchange protocol.
For further detail see :ref:`dlpack`.

.. autosummary::
   :toctree: generated/

   cupy.from_dlpack


.. _kernel_param_opt:

Automatic Kernel Parameters Optimizations (:mod:`cupyx.optimizing`)
-------------------------------------------------------------------

.. module:: cupyx.optimizing
.. autosummary::
   :toctree: generated/

   cupyx.optimizing.optimize
Low-level CUDA support
======================

.. _device_management:

Device management
-----------------

.. autosummary::
   :toctree: generated/

   cupy.cuda.Device


Memory management
-----------------

.. autosummary::
   :toctree: generated/

   cupy.get_default_memory_pool
   cupy.get_default_pinned_memory_pool
   cupy.cuda.Memory
   cupy.cuda.MemoryAsync
   cupy.cuda.ManagedMemory
   cupy.cuda.UnownedMemory
   cupy.cuda.PinnedMemory
   cupy.cuda.MemoryPointer
   cupy.cuda.PinnedMemoryPointer
   cupy.cuda.malloc_managed
   cupy.cuda.malloc_async
   cupy.cuda.alloc
   cupy.cuda.alloc_pinned_memory
   cupy.cuda.get_allocator
   cupy.cuda.set_allocator
   cupy.cuda.using_allocator
   cupy.cuda.set_pinned_memory_allocator
   cupy.cuda.MemoryPool
   cupy.cuda.MemoryAsyncPool
   cupy.cuda.PinnedMemoryPool
   cupy.cuda.PythonFunctionAllocator
   cupy.cuda.CFunctionAllocator


Memory hook
-----------

.. autosummary::
   :toctree: generated/

   cupy.cuda.MemoryHook
   cupy.cuda.memory_hooks.DebugPrintHook
   cupy.cuda.memory_hooks.LineProfileHook


.. _stream_event_api:

Streams and events
------------------

.. autosummary::
   :toctree: generated/

   cupy.cuda.Stream
   cupy.cuda.ExternalStream
   cupy.cuda.get_current_stream
   cupy.cuda.Event
   cupy.cuda.get_elapsed_time


.. _graph_api:

Graphs
------

.. autosummary::
   :toctree: generated/

   cupy.cuda.Graph


Texture and surface memory
--------------------------

.. autosummary::
   :toctree: generated/

   cupy.cuda.texture.ChannelFormatDescriptor
   cupy.cuda.texture.CUDAarray
   cupy.cuda.texture.ResourceDescriptor
   cupy.cuda.texture.TextureDescriptor
   cupy.cuda.texture.TextureObject
   cupy.cuda.texture.SurfaceObject
   cupy.cuda.texture.TextureReference


Profiler
--------

.. autosummary::
   :toctree: generated/

   cupy.cuda.profile
   cupy.cuda.profiler.initialize
   cupy.cuda.profiler.start
   cupy.cuda.profiler.stop
   cupy.cuda.nvtx.Mark
   cupy.cuda.nvtx.MarkC
   cupy.cuda.nvtx.RangePush
   cupy.cuda.nvtx.RangePushC
   cupy.cuda.nvtx.RangePop


NCCL
----

.. autosummary::
   :toctree: generated/

   cupy.cuda.nccl.NcclCommunicator
   cupy.cuda.nccl.get_build_version
   cupy.cuda.nccl.get_version
   cupy.cuda.nccl.get_unique_id
   cupy.cuda.nccl.groupStart
   cupy.cuda.nccl.groupEnd


.. _runtime_api:

Runtime API
-----------

CuPy wraps CUDA Runtime APIs to provide the native CUDA operations.
Please check the `CUDA Runtime API documentation <https://docs.nvidia.com/cuda/cuda-runtime-api/index.html>`_
to use these functions.

.. autosummary::
   :toctree: generated/

   cupy.cuda.runtime.driverGetVersion
   cupy.cuda.runtime.runtimeGetVersion
   cupy.cuda.runtime.getDevice
   cupy.cuda.runtime.getDeviceProperties
   cupy.cuda.runtime.deviceGetAttribute
   cupy.cuda.runtime.deviceGetByPCIBusId
   cupy.cuda.runtime.deviceGetPCIBusId
   cupy.cuda.runtime.deviceGetDefaultMemPool
   cupy.cuda.runtime.deviceGetMemPool
   cupy.cuda.runtime.deviceSetMemPool
   cupy.cuda.runtime.memPoolTrimTo
   cupy.cuda.runtime.getDeviceCount
   cupy.cuda.runtime.setDevice
   cupy.cuda.runtime.deviceSynchronize
   cupy.cuda.runtime.deviceCanAccessPeer
   cupy.cuda.runtime.deviceEnablePeerAccess
   cupy.cuda.runtime.deviceGetLimit
   cupy.cuda.runtime.deviceSetLimit
   cupy.cuda.runtime.malloc
   cupy.cuda.runtime.mallocManaged
   cupy.cuda.runtime.malloc3DArray
   cupy.cuda.runtime.mallocArray
   cupy.cuda.runtime.mallocAsync
   cupy.cuda.runtime.hostAlloc
   cupy.cuda.runtime.hostRegister
   cupy.cuda.runtime.hostUnregister
   cupy.cuda.runtime.free
   cupy.cuda.runtime.freeHost
   cupy.cuda.runtime.freeArray
   cupy.cuda.runtime.freeAsync
   cupy.cuda.runtime.memGetInfo
   cupy.cuda.runtime.memcpy
   cupy.cuda.runtime.memcpyAsync
   cupy.cuda.runtime.memcpyPeer
   cupy.cuda.runtime.memcpyPeerAsync
   cupy.cuda.runtime.memcpy2D
   cupy.cuda.runtime.memcpy2DAsync
   cupy.cuda.runtime.memcpy2DFromArray
   cupy.cuda.runtime.memcpy2DFromArrayAsync
   cupy.cuda.runtime.memcpy2DToArray
   cupy.cuda.runtime.memcpy2DToArrayAsync
   cupy.cuda.runtime.memcpy3D
   cupy.cuda.runtime.memcpy3DAsync
   cupy.cuda.runtime.memset
   cupy.cuda.runtime.memsetAsync
   cupy.cuda.runtime.memPrefetchAsync
   cupy.cuda.runtime.memAdvise
   cupy.cuda.runtime.pointerGetAttributes
   cupy.cuda.runtime.streamCreate
   cupy.cuda.runtime.streamCreateWithFlags
   cupy.cuda.runtime.streamDestroy
   cupy.cuda.runtime.streamSynchronize
   cupy.cuda.runtime.streamAddCallback
   cupy.cuda.runtime.streamQuery
   cupy.cuda.runtime.streamWaitEvent
   cupy.cuda.runtime.launchHostFunc
   cupy.cuda.runtime.eventCreate
   cupy.cuda.runtime.eventCreateWithFlags
   cupy.cuda.runtime.eventDestroy
   cupy.cuda.runtime.eventElapsedTime
   cupy.cuda.runtime.eventQuery
   cupy.cuda.runtime.eventRecord
   cupy.cuda.runtime.eventSynchronize
   cupy.cuda.runtime.ipcGetMemHandle
   cupy.cuda.runtime.ipcOpenMemHandle
   cupy.cuda.runtime.ipcCloseMemHandle
   cupy.cuda.runtime.ipcGetEventHandle
   cupy.cuda.runtime.ipcOpenEventHandle
.. module:: cupyx.scipy.signal

Signal processing (:mod:`cupyx.scipy.signal`)
=============================================

.. Hint:: `SciPy API Reference: Signal processing (scipy.signal) <https://docs.scipy.org/doc/scipy/reference/signal.html>`_

Convolution
-----------

.. autosummary::
   :toctree: generated/

   convolve
   correlate
   fftconvolve
   oaconvolve
   convolve2d
   correlate2d
   sepfir2d
   choose_conv_method


Filtering
---------

.. autosummary::
   :toctree: generated/

   order_filter
   medfilt
   medfilt2d
   wiener
.. _environment:

Environment variables
=====================

For runtime
-----------

Here are the environment variables that CuPy uses at runtime.

.. envvar:: CUDA_PATH

  Path to the directory containing CUDA.
  The parent of the directory containing ``nvcc`` is used as default.
  When ``nvcc`` is not found, ``/usr/local/cuda`` is used.
  See :ref:`install_cuda` for details.

.. envvar:: CUPY_CACHE_DIR

  Default: ``${HOME}/.cupy/kernel_cache``

  Path to the directory to store kernel cache.
  See :doc:`../user_guide/performance` for details.

.. envvar:: CUPY_CACHE_SAVE_CUDA_SOURCE

  Default: ``0``

  If set to ``1``, CUDA source file will be saved along with compiled binary in the cache directory for debug purpose.
  Note: the source file will not be saved if the compiled binary is already stored in the cache.

.. envvar:: CUPY_CACHE_IN_MEMORY

  Default: ``0``

  If set to ``1``, :envvar:`CUPY_CACHE_DIR` and :envvar:`CUPY_CACHE_SAVE_CUDA_SOURCE` will be ignored, and the cache is in memory.
  This environment variable allows reducing disk I/O, but is ignoed when ``nvcc`` is set to be the compiler backend.

.. envvar:: CUPY_DUMP_CUDA_SOURCE_ON_ERROR

  Default: ``0``

  If set to ``1``, when CUDA kernel compilation fails,
  CuPy dumps CUDA kernel code to standard error.

.. envvar:: CUPY_CUDA_COMPILE_WITH_DEBUG

  Default: ``0``

  If set to ``1``, CUDA kernel will be compiled with debug information (``--device-debug`` and ``--generate-line-info``).

.. envvar:: CUPY_GPU_MEMORY_LIMIT

  Default: ``0`` (unlimited)

  The amount of memory that can be allocated for each device.
  The value can be specified in absolute bytes or fraction (e.g., ``"90%"``) of the total memory of each GPU.
  See :doc:`../user_guide/memory` for details.

.. envvar:: CUPY_SEED

  Set the seed for random number generators.

.. envvar:: CUPY_EXPERIMENTAL_SLICE_COPY

  Default: ``0``
  
  If set to ``1``, the following syntax is enabled::

    cupy_ndarray[:] = numpy_ndarray

.. envvar:: CUPY_ACCELERATORS

  Default: ``""`` (no accelerators)

  A comma-separated string of backend names (``cub`` or ``cutensor``) which indicates the acceleration backends used in CuPy operations and its priority.
  All accelerators are disabled by default.

.. envvar:: CUPY_TF32

  Default: ``0``

  If set to ``1``, it allows CUDA libraries to use Tensor Cores TF32 compute for 32-bit floating point compute.

.. envvar:: CUPY_CUDA_ARRAY_INTERFACE_SYNC

  Default: ``1``

  This controls CuPy's behavior as a Consumer.
  If set to ``0``, a stream synchronization will *not* be performed when a device array provided by an external library that implements the CUDA Array Interface is being consumed by CuPy.
  For more detail, see the `Synchronization`_ requirement in the CUDA Array Interface v3 documentation.

.. envvar:: CUPY_CUDA_ARRAY_INTERFACE_EXPORT_VERSION

  Default: ``3``

  This controls CuPy's behavior as a Producer.
  If set to ``2``, the CuPy stream on which the data is being operated will not be exported and thus the Consumer (another library) will not perform any stream synchronization.
  For more detail, see the `Synchronization`_ requirement in the CUDA Array Interface v3 documentation.

.. envvar:: CUPY_DLPACK_EXPORT_VERSION

  Default: ``0.6``

  This controls CuPy's DLPack support. Currently, setting a value smaller than 0.6 would disguise managed memory as normal device memory, which enables data exchanges with libraries that have not updated their DLPack support, whereas starting 0.6 CUDA managed memory can be correctly recognized as a valid device type.

.. envvar:: NVCC

  Default: ``nvcc``

  Define the compiler to use when compiling CUDA source.
  Note that most CuPy kernels are built with NVRTC; this environment variable is only effective for :class:`~cupy.RawKernel`/:class:`~cupy.RawModule` with the ``nvcc`` backend or when using ``cub`` as the accelerator.

.. envvar:: CUPY_CUDA_PER_THREAD_DEFAULT_STREAM

  Default: ``0``

  If set to ``1``, CuPy will use the CUDA per-thread default stream, effectively causing each host thread to automatically execute in its own stream, unless the CUDA default (``null``) stream or a user-created stream is specified.
  If set to ``0`` (default), the CUDA default (``null``) stream is used, unless the per-thread default stream (``ptds``) or a user-created stream is specified.

.. envvar:: CUPY_COMPILE_WITH_PTX

  Default: ``0``

  By default, CuPy directly compiles kernels into SASS (CUBIN) to support `CUDA Enhanced Compatibility <https://docs.nvidia.com/deploy/cuda-compatibility/>`_
  If set to ``1``, CuPy instead compiles kernels into PTX and lets CUDA Driver assemble SASS from PTX.
  This option is only effective for CUDA 11.1 or later; CuPy always compiles into PTX on earlier CUDA versions. Also, this option only applies when NVRTC is selected as the compilation backend. NVCC backend always compiles into SASS (CUBIN).

CUDA Toolkit Environment Variables
  In addition to the environment variables listed above, as in any CUDA programs, all of the CUDA environment variables listed in the `CUDA Toolkit Documentation`_ will also be honored.

.. note::

  When :envvar:`CUPY_ACCELERATORS` or :envvar:`NVCC` environment variables are set, g++-6 or later is required as the runtime host compiler.
  Please refer to :ref:`install_cupy_from_source` for the details on how to install g++.

.. _CUDA Toolkit Documentation: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#env-vars

.. _Synchronization: https://numba.readthedocs.io/en/latest/cuda/cuda_array_interface.html#synchronization


For installation
----------------

These environment variables are used during installation (building CuPy from source).

.. envvar:: CUTENSOR_PATH

  Path to the cuTENSOR root directory that contains ``lib`` and ``include`` directories. (experimental)

.. envvar:: CUPY_INSTALL_USE_HIP

  Default: ``0``

  If set to ``1``, CuPy is built for AMD ROCm Platform (experimental).
  For building the ROCm support, see :ref:`install_hip` for further detail.

.. envvar:: CUPY_USE_CUDA_PYTHON

  Default: ``0``

  If set to ``1``, CuPy is built using `CUDA Python <https://github.com/NVIDIA/cuda-python>`_.

.. envvar:: CUPY_NVCC_GENERATE_CODE

  Build CuPy for a particular CUDA architecture. For example::

    CUPY_NVCC_GENERATE_CODE="arch=compute_60,code=sm_60"

  For specifying multiple archs, concatenate the ``arch=...`` strings with semicolons (``;``).
  If ``current`` is specified, then it will automatically detect the currently installed GPU architectures in build time.
  When this is not set, the default is to support all architectures.

.. envvar:: CUPY_NUM_BUILD_JOBS

  Default: ``4``

  To enable or disable parallel build, sets the number of processes used to build the extensions in parallel.


.. envvar:: CUPY_NUM_NVCC_THREADS

  Default: ``2``

  To enable or disable nvcc parallel compilation, sets the number of threads used to compile files using nvcc.

Additionally, the environment variables :envvar:`CUDA_PATH` and :envvar:`NVCC` are also respected at build time.
The N-dimensional array (:class:`ndarray <cupy.ndarray>`)
=========================================================

:class:`cupy.ndarray` is the CuPy counterpart of NumPy :class:`numpy.ndarray`.
It provides an intuitive interface for a fixed-size multidimensional array which resides
in a CUDA device.

For the basic concept of ``ndarray``\s, please refer to the `NumPy documentation <https://numpy.org/doc/stable/reference/arrays.ndarray.html>`_.


.. TODO(kmaehashi): use currentmodule:: cupy
.. autosummary::
   :toctree: generated/

   cupy.ndarray


Conversion to/from NumPy arrays
-------------------------------

:class:`cupy.ndarray` and :class:`numpy.ndarray` are not implicitly convertible to each other.
That means, NumPy functions cannot take :class:`cupy.ndarray`\s as inputs, and vice versa.

- To convert :class:`numpy.ndarray` to :class:`cupy.ndarray`, use :func:`cupy.array` or :func:`cupy.asarray`.
- To convert :class:`cupy.ndarray` to :class:`numpy.ndarray`, use :func:`cupy.asnumpy` or :meth:`cupy.ndarray.get`.

Note that converting between :class:`cupy.ndarray` and :class:`numpy.ndarray` incurs data transfer between
the host (CPU) device and the GPU device, which is costly in terms of performance.


.. TODO(kmaehashi): use currentmodule:: cupy
.. autosummary::
   :toctree: generated/

   cupy.array
   cupy.asarray
   cupy.asnumpy


Code compatibility features
---------------------------

:class:`cupy.ndarray` is designed to be interchangeable with :class:`numpy.ndarray` in terms of code compatibility as much as possible.
But occasionally, you will need to know whether the arrays you're handling are :class:`cupy.ndarray` or :class:`numpy.ndarray`.
One example is when invoking module-level functions such as :func:`cupy.sum` or :func:`numpy.sum`.
In such situations, :func:`cupy.get_array_module` can be used.

.. autosummary::
   :toctree: generated/

   cupy.get_array_module

.. autosummary::
   :toctree: generated/

   cupyx.scipy.get_array_module
.. module:: cupyx.scipy.special

Special functions (:mod:`cupyx.scipy.special`)
===============================================

.. Hint:: `SciPy API Reference: Special functions (scipy.special) <https://docs.scipy.org/doc/scipy/reference/special.html>`_

Bessel functions
----------------

.. autosummary::
   :toctree: generated/

   j0
   j1
   y0
   y1
   yn
   i0
   i1


Raw statistical functions
-------------------------

.. seealso:: :mod:`cupyx.scipy.stats`

.. autosummary::
   :toctree: generated/

   ndtr
   logit
   expit
   log_expit


Information Theory functions
----------------------------

.. autosummary::
   :toctree: generated/

   entr
   rel_entr
   kl_div
   huber
   pseudo_huber


Gamma and related functions
---------------------------

.. autosummary::
   :toctree: generated/

   gamma
   gammaln
   gammainc
   gammaincinv
   gammaincc
   gammainccinv
   psi
   polygamma
   digamma
   poch


Error function and Fresnel integrals
------------------------------------

.. autosummary::
   :toctree: generated/

   erf
   erfc
   erfcx
   erfinv
   erfcinv


Legendre functions
---------------------------

.. autosummary::
   :toctree: generated/

   lpmv
   sph_harm


Other special functions
-----------------------

.. autosummary::
   :toctree: generated/

   zeta


Convenience functions
-----------------------

.. autosummary::
   :toctree: generated/

   cbrt
   exp10
   exp2
   radian
   cosdg
   sindg
   tandg
   cotdg
   log1p
   expm1
   round
   xlogy
   xlog1py
   sinc
Polynomials
===========

.. Hint:: `NumPy API Reference: Polynomials <https://numpy.org/doc/stable/reference/routines.polynomials.html>`_

Power Series (:mod:`cupy.polynomial.polynomial`)
------------------------------------------------

.. Hint:: `NumPy API Reference: Power Series (numpy.polynomial.polynomial) <https://numpy.org/doc/stable/reference/routines.polynomials.polynomial.html>`_

Misc Functions
~~~~~~~~~~~~~~

.. module:: cupy.polynomial.polynomial
.. autosummary::
   :toctree: generated/

   polyvander
   polycompanion


Polyutils
---------

.. Hint:: `NumPy API Reference: Polyutils <https://numpy.org/doc/stable/reference/routines.polynomials.polyutils.html>`_

Functions
~~~~~~~~~

.. module:: cupy.polynomial.polyutils
.. autosummary::
   :toctree: generated/

   as_series
   trimseq
   trimcoef


Poly1d
------

.. Hint:: `NumPy API Reference: Poly1d <https://numpy.org/doc/stable/reference/routines.polynomials.poly1d.html>`_

.. currentmodule:: cupy

Basics
~~~~~~

.. autosummary::
   :toctree: generated/

   poly1d
   polyval
   roots


Fitting
~~~~~~~

.. autosummary::
   :toctree: generated/

   polyfit


Arithmetic
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   polyadd
   polysub
   polymul
.. module:: cupy.linalg

Linear algebra (:mod:`cupy.linalg`)
===================================

.. Hint:: `NumPy API Reference: Linear algebra (numpy.linalg) <https://numpy.org/doc/stable/reference/routines.linalg.html>`_

.. seealso:: :doc:`scipy_linalg`

.. currentmodule:: cupy

Matrix and vector products
--------------------------

.. autosummary::
   :toctree: generated/

   dot
   vdot
   inner
   outer
   matmul
   tensordot
   einsum
   linalg.matrix_power
   kron

Decompositions
--------------

.. autosummary::
   :toctree: generated/

   linalg.cholesky
   linalg.qr
   linalg.svd

Matrix eigenvalues
------------------

.. autosummary::
   :toctree: generated/

   linalg.eigh
   linalg.eigvalsh

Norms and other numbers
-----------------------

.. autosummary::
   :toctree: generated/

   linalg.norm
   linalg.det
   linalg.matrix_rank
   linalg.slogdet
   trace


Solving equations and inverting matrices
----------------------------------------

.. autosummary::
   :toctree: generated/

   linalg.solve
   linalg.tensorsolve
   linalg.lstsq
   linalg.inv
   linalg.pinv
   linalg.tensorinv
:orphan:


cuFFT Plan Cache
----------------

.. autoclass:: cupy.fft._cache.PlanCache
   :members:
   :undoc-members:
.. module:: cupyx.scipy.linalg

Linear algebra (:mod:`cupyx.scipy.linalg`)
==========================================

.. Hint:: `SciPy API Reference: Linear algebra (scipy.linalg) <https://docs.scipy.org/doc/scipy/reference/linalg.html>`_

Basics
------

.. autosummary::
   :toctree: generated/

   solve_triangular
   tril
   triu


Decompositions
--------------

.. autosummary::
   :toctree: generated/

   lu
   lu_factor
   lu_solve


Special Matrices
----------------

.. autosummary::
   :toctree: generated/

   block_diag
   circulant
   companion
   convolution_matrix
   dft
   fiedler
   fiedler_companion
   hadamard
   hankel
   helmert
   hilbert
   kron
   leslie
   toeplitz
   tri
----------------
Distributed
----------------

The following pages describe the APIs used to easily perform communication
between different processes in CuPy.


.. module:: cupyx.distributed

.. autosummary::
   :toctree: generated/

   init_process_group
   NCCLBackend
.. module:: cupyx.scipy.sparse.csgraph

Compressed sparse graph routines (:mod:`cupyx.scipy.sparse.csgraph`)
====================================================================

.. note::

   The ``csgraph`` module uses ``pylibcugraph`` as a backend.
   You need to install `pylibcugraph package <https://anaconda.org/rapidsai/pylibcugraph>` from ``rapidsai`` Conda channel to use features listed on this page.

.. note::
   Currently, the ``csgraph`` module is not supported on AMD ROCm platforms.

.. Hint:: `SciPy API Reference: Compressed sparse graph routines (scipy.sparse.csgraph) <https://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html>`_

Contents
---------------------

.. autosummary::
   :toctree: generated/

   connected_components
Statistics
==========

.. Hint:: `NumPy API Reference: Statistics <https://numpy.org/doc/stable/reference/routines.statistics.html>`_

.. currentmodule:: cupy

Order statistics
----------------

.. autosummary::
   :toctree: generated/

   amin
   amax
   nanmin
   nanmax
   ptp
   percentile
   quantile


Averages and variances
----------------------

.. autosummary::
   :toctree: generated/

   median
   average
   mean
   std
   var
   nanmedian
   nanmean
   nanstd
   nanvar


Correlations
------------

.. autosummary::
   :toctree: generated/

   corrcoef
   correlate
   cov


Histograms
----------

.. autosummary::
   :toctree: generated/

   histogram
   histogram2d
   histogramdd
   bincount
   digitize
Set routines
============

.. Hint:: `NumPy API Reference: Set routines <https://numpy.org/doc/stable/reference/routines.set.html>`_

.. currentmodule:: cupy

Making proper sets
------------------

.. autosummary::
   :toctree: generated/

   unique


Boolean operations
------------------

.. autosummary::
   :toctree: generated/

   in1d
   isin
Padding arrays
==============

.. Hint:: `NumPy API Reference: Padding arrays <https://numpy.org/doc/stable/reference/routines.padding.html>`_

.. currentmodule:: cupy
.. autosummary::
   :toctree: generated/

   pad
Universal functions (:class:`cupy.ufunc`)
=========================================

.. Hint:: `NumPy API Reference: Universal functions (numpy.ufunc) <https://numpy.org/doc/stable/reference/ufuncs.html>`_

.. currentmodule:: cupy

CuPy provides universal functions (a.k.a. ufuncs) to support various elementwise operations.
CuPy's ufunc supports following features of NumPy's one:

- Broadcasting
- Output type determination
- Casting rules

CuPy's ufunc currently does not provide methods such as ``reduce``, ``accumulate``, ``reduceat``, ``outer``, and ``at``.


ufunc
-----

.. autosummary::
   :toctree: generated/

   ufunc


Available ufuncs
----------------

Math operations
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   add
   subtract
   multiply
   matmul
   divide
   logaddexp
   logaddexp2
   true_divide
   floor_divide
   negative
   positive
   power
   remainder
   mod
   fmod
   absolute
   rint
   sign
   conj
   conjugate
   exp
   exp2
   log
   log2
   log10
   expm1
   log1p
   sqrt
   square
   cbrt
   reciprocal
   gcd
   lcm


Trigonometric functions
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   sin
   cos
   tan
   arcsin
   arccos
   arctan
   arctan2
   hypot
   sinh
   cosh
   tanh
   arcsinh
   arccosh
   arctanh
   degrees
   radians
   deg2rad
   rad2deg


Bit-twiddling functions
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   bitwise_and
   bitwise_or
   bitwise_xor
   invert
   left_shift
   right_shift


Comparison functions
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   greater
   greater_equal
   less
   less_equal
   not_equal
   equal
   logical_and
   logical_or
   logical_xor
   logical_not
   maximum
   minimum
   fmax
   fmin


Floating functions
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   isfinite
   isinf
   isnan
   signbit
   copysign
   nextafter
   modf
   ldexp
   frexp
   fmod
   floor
   ceil
   trunc


ufunc.at
--------

Currently, CuPy does not support ``at`` for ufuncs in general.
However, :func:`cupyx.scatter_add` can substitute ``add.at`` as both behave identically.


Generalized Universal Functions
-------------------------------

.. currentmodule:: cupyx

In addition to regular ufuncs, CuPy also provides a wrapper class to convert
regular cupy functions into Generalized Universal Functions as in NumPy `<https://numpy.org/doc/stable/reference/c-api/generalized-ufuncs.html>`_.
This allows to automatically use keyword arguments such as ``axes``, ``order``, ``dtype``
without needing to explicitly implement them in the wrapped function.


.. autosummary::
   :toctree: generated/

   GeneralizedUFunc
.. module:: cupyx.scipy.sparse.linalg

Sparse linear algebra (:mod:`cupyx.scipy.sparse.linalg`)
========================================================

.. Hint:: `SciPy API Reference: Sparse linear algebra (scipy.sparse.linalg) <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html>`_

Abstract linear operators
-------------------------

.. autosummary::
   :toctree: generated/

   LinearOperator
   aslinearoperator


Matrix norms
------------

.. autosummary::
   :toctree: generated/

   norm


Solving linear problems
-----------------------

Direct methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   spsolve
   spsolve_triangular
   factorized

Iterative methods for linear equation systems:

.. autosummary::
   :toctree: generated/

   cg
   gmres
   cgs
   minres

Iterative methods for least-squares problems:

.. autosummary::
   :toctree: generated/

   lsqr
   lsmr


Matrix factorizations
---------------------

Eigenvalue problems:

.. autosummary::
   :toctree: generated/

   eigsh
   lobpcg

Singular values problems:

.. autosummary::
   :toctree: generated/

   svds

Complete or incomplete LU factorizations:

.. autosummary::
   :toctree: generated/

   splu
   spilu
   SuperLU
Array API Complaint Object
==========================

:class:`~cupy.array_api._array_object.Array` is a wrapper class built upon :class:`cupy.ndarray`
to enforce strict complaince with the array API standard. See the
`documentation <https://data-apis.org/array-api/latest/API_specification/array_object.html>`_
for detail.

This object should not be constructed directly. Rather, use one of the
`creation functions <https://data-apis.org/array-api/latest/API_specification/creation_functions.html>`_,
such as :func:`cupy.array_api.asarray`.

.. currentmodule:: cupy.array_api._array_object

.. autosummary::
   :toctree: generated/

   Array
.. module:: cupy.fft

Discrete Fourier Transform (:mod:`cupy.fft`)
============================================

.. Hint:: `NumPy API Reference: Discrete Fourier Transform (numpy.fft) <https://numpy.org/doc/stable/reference/routines.fft.html>`_

.. seealso:: :doc:`scipy_fft`, :doc:`../user_guide/fft`

Standard FFTs
-------------

.. autosummary::
   :toctree: generated/

   fft
   ifft
   fft2
   ifft2
   fftn
   ifftn


Real FFTs
---------

.. autosummary::
   :toctree: generated/

   rfft
   irfft
   rfft2
   irfft2
   rfftn
   irfftn


Hermitian FFTs
--------------

.. autosummary::
   :toctree: generated/

   hfft
   ihfft


Helper routines
---------------

.. autosummary::
   :toctree: generated/

   fftfreq
   rfftfreq
   fftshift
   ifftshift

CuPy-specific APIs
------------------

See the description below for details.

.. autosummary::
   :toctree: generated/

   config.set_cufft_callbacks
   config.set_cufft_gpus
   config.get_plan_cache
   config.show_plan_cache_info


Normalization
-------------
The default normalization (``norm`` is ``"backward"`` or ``None``) has the direct transforms unscaled and the inverse transforms scaled by :math:`1/n`.
If the keyword argument ``norm`` is ``"forward"``, it is the exact opposite of ``"backward"``:
the direct transforms are scaled by :math:`1/n` and the inverse transforms are unscaled.
Finally, if the keyword argument ``norm`` is ``"ortho"``, both transforms are scaled by :math:`1/\sqrt{n}`.

Code compatibility features
---------------------------
FFT functions of NumPy always return numpy.ndarray which type is ``numpy.complex128`` or ``numpy.float64``.
CuPy functions do not follow the behavior, they will return ``numpy.complex64`` or ``numpy.float32`` if the type of the input is ``numpy.float16``, ``numpy.float32``, or ``numpy.complex64``.

Internally, ``cupy.fft`` always generates a *cuFFT plan* (see the `cuFFT documentation`_ for detail) corresponding to the desired transform. When possible, an n-dimensional plan will be used, as opposed to applying separate 1D plans for each axis to be transformed. Using n-dimensional planning can provide better performance for multidimensional transforms, but requires more GPU memory than separable 1D planning. The user can disable n-dimensional planning by setting ``cupy.fft.config.enable_nd_planning = False``. This ability to adjust the planning type is a deviation from the NumPy API, which does not use precomputed FFT plans.

Moreover, the automatic plan generation can be suppressed by using an existing plan returned by :func:`cupyx.scipy.fftpack.get_fft_plan` as a context manager. This is again a deviation from NumPy.

Finally, when using the high-level NumPy-like FFT APIs as listed above, internally the cuFFT plans are cached for possible reuse. The plan cache can be retrieved by :func:`~cupy.fft.config.get_plan_cache`, and its current status can be queried by :func:`~cupy.fft.config.show_plan_cache_info`. For finer control of the plan cache, see :doc:`plan_cache`.


Multi-GPU FFT
-------------
:mod:`cupy.fft` can use multiple GPUs. To enable (disable) this feature, set :data:`cupy.fft.config.use_multi_gpus` to ``True`` (``False``). Next, to set the number of GPUs or the participating GPU IDs, use the function :func:`cupy.fft.config.set_cufft_gpus`. All of the limitations listed in the `cuFFT documentation`_ apply here. In particular, using more than one GPU does not guarantee better performance.


.. _cuFFT documentation: https://docs.nvidia.com/cuda/cufft/index.html
Array creation routines
=======================

.. Hint:: `NumPy API Reference: Array creation routines <https://numpy.org/doc/stable/reference/routines.array-creation.html>`_

.. currentmodule:: cupy

Ones and zeros
--------------

.. autosummary::
   :toctree: generated/

   empty
   empty_like
   eye
   identity
   ones
   ones_like
   zeros
   zeros_like
   full
   full_like


From existing data
------------------

.. autosummary::
   :toctree: generated/

   array
   asarray
   asanyarray
   ascontiguousarray
   copy
   frombuffer
   fromfile
   fromfunction
   fromiter
   fromstring
   loadtxt


Numerical ranges
----------------

.. autosummary::
   :toctree: generated/

   arange
   linspace
   logspace
   meshgrid
   mgrid
   ogrid


Building matrices
-----------------

.. autosummary::
   :toctree: generated/

   diag
   diagflat
   tri
   tril
   triu
   vander
Data type routines
==================

.. Hint:: `NumPy API Reference: Data type routines <https://numpy.org/doc/stable/reference/routines.dtype.html>`_

.. currentmodule:: cupy

.. autosummary::
   :toctree: generated/

   can_cast
   result_type
   common_type

.. csv-table::
   :align: left

   ``promote_types`` (alias of :func:`numpy.promote_types`)
   ``min_scalar_type`` (alias of :func:`numpy.min_scalar_type`)
   ``obj2sctype`` (alias of :func:`numpy.obj2sctype`)

Creating data types
-------------------

.. csv-table::
   :align: left

   ``dtype`` (alias of :class:`numpy.dtype`)
   ``format_parser`` (alias of :class:`numpy.format_parser`)

Data type information
---------------------

.. csv-table::
   :align: left

   ``finfo`` (alias of :class:`numpy.finfo`)
   ``iinfo`` (alias of :class:`numpy.iinfo`)
   ``MachAr`` (alias of :class:`numpy.MachAr`)

Data type testing
-----------------

.. csv-table::
   :align: left

   ``issctype`` (alias of :func:`numpy.issctype`)
   ``issubdtype`` (alias of :func:`numpy.issubdtype`)
   ``issubsctype`` (alias of :func:`numpy.issubsctype`)
   ``issubclass_`` (alias of :func:`numpy.issubclass_`)
   ``find_common_type`` (alias of :func:`numpy.find_common_type`)

Miscellaneous
-------------

.. csv-table::
   :align: left

   ``typename`` (alias of :func:`numpy.typename`)
   ``sctype2char`` (alias of :func:`numpy.sctype2char`)
   ``mintypecode`` (alias of :func:`numpy.mintypecode`)
.. _cupy_reference:

*************
API Reference
*************

* :ref:`genindex`
* :ref:`modindex`

----

.. currentmodule:: cupy

..
  For NumPy/SciPy-compatible APIs (ndarray, ufunc, routines, scipy),
  omit the module name prefix in the API list, following the convension in
  NumPy/SciPy documentation.
  For CuPy-specific APIs, use fully-qualified names.

.. toctree::
   :maxdepth: 2

   ndarray
   ufunc
   routines
   scipy
   ext
   cuda
   kernel
   distributed
   environment
   comparison
   array_api
Input and output
================

.. Hint:: `NumPy API Reference: Input and output <https://numpy.org/doc/stable/reference/routines.io.html>`_

.. currentmodule:: cupy


NumPy binary files (NPY, NPZ)
-----------------------------

.. autosummary::
   :toctree: generated/

   load
   save
   savez
   savez_compressed

Text files
-----------------------------

.. autosummary::
   :toctree: generated/

   loadtxt
   savetxt
   genfromtxt
   fromstring

String formatting
-----------------

.. autosummary::
   :toctree: generated/

   array2string
   array_repr
   array_str
   format_float_positional


Base-n representations
----------------------

.. autosummary::
   :toctree: generated/

   binary_repr
   base_repr
Python Array API Support
========================

The `Python array API standard <https://data-apis.org/array-api/latest/>`_ aims to provide a coherent set of
APIs for array and tensor libraries developed by the community to build upon. This solves the API fragmentation
issue across the community by offering concrete function signatures, semantics and scopes of coverage, enabling
writing backend-agnostic codes for better portability.

CuPy provides **experimental** support based on NumPy's `NEP-47 <https://numpy.org/neps/nep-0047-array-api-standard.html>`_,
which is in turn based on the draft standard to be finalized in 2021. All of the functionalities can be accessed
through the :mod:`cupy.array_api` namespace.

The key difference between NumPy and CuPy is that we are a GPU-only library, therefore CuPy users should be aware
of potential `device management <https://data-apis.org/array-api/latest/design_topics/device_support.html>`_ issues.
Same as in regular CuPy codes, the GPU-to-use can be specified via the :class:`~cupy.cuda.Device` objects, see
:ref:`device_management`.

.. toctree::
   :maxdepth: 2

   array_api_functions
   array_api_array
:orphan:

This document has been moved to :doc:`scipy_stats`.
:orphan:

This document has been moved to :doc:`scipy_signal`.
:orphan:

This document has been moved to :doc:`scipy_fftpack`.
:orphan:


This document has been moved to :doc:`scipy_special`.
----------------
Routines (SciPy)
----------------

The following pages describe SciPy-compatible routines.
These functions cover a subset of
`SciPy routines <https://docs.scipy.org/doc/scipy/reference/#api-reference>`_.


.. module:: cupyx.scipy

.. toctree::
   :maxdepth: 2

   scipy_fft
   scipy_fftpack
   scipy_linalg
   scipy_ndimage
   scipy_signal
   scipy_sparse
   scipy_sparse_linalg
   scipy_sparse_csgraph
   scipy_special
   scipy_stats
.. module:: cupy.testing

Test support (:mod:`cupy.testing`)
==================================

.. Hint:: `NumPy API Reference: Test support (numpy.testing) <https://numpy.org/doc/stable/reference/routines.testing.html>`_

Asserts
-------

.. Hint:: These APIs can accept both :class:`numpy.ndarray` and :class:`cupy.ndarray`.

.. autosummary::
   :toctree: generated/

   assert_array_almost_equal
   assert_allclose
   assert_array_almost_equal_nulp
   assert_array_max_ulp
   assert_array_equal
   assert_array_less

CuPy-specific APIs
------------------

Asserts
~~~~~~~

.. autosummary::
   :toctree: generated/

   assert_array_list_equal


NumPy-CuPy Consistency Check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following decorators are for testing consistency
between CuPy's functions and corresponding NumPy's ones.

.. autosummary::
   :toctree: generated/

   numpy_cupy_allclose
   numpy_cupy_array_almost_equal
   numpy_cupy_array_almost_equal_nulp
   numpy_cupy_array_max_ulp
   numpy_cupy_array_equal
   numpy_cupy_array_list_equal
   numpy_cupy_array_less


Parameterized dtype Test
~~~~~~~~~~~~~~~~~~~~~~~~

The following decorators offer the standard way for
parameterized test with respect to single or the
combination of dtype(s).

.. autosummary::
   :toctree: generated/

   for_dtypes
   for_all_dtypes
   for_float_dtypes
   for_signed_dtypes
   for_unsigned_dtypes
   for_int_dtypes
   for_complex_dtypes
   for_dtypes_combination
   for_all_dtypes_combination
   for_signed_dtypes_combination
   for_unsigned_dtypes_combination
   for_int_dtypes_combination


Parameterized order Test
~~~~~~~~~~~~~~~~~~~~~~~~

The following decorators offer the standard way to parameterize tests with
orders.

.. autosummary::
   :toctree: generated/

   for_orders
   for_CF_orders
Mathematical functions
======================

.. Hint:: `NumPy API Reference: Mathematical functions <https://numpy.org/doc/stable/reference/routines.math.html>`_

.. currentmodule:: cupy

Trigonometric functions
-----------------------

.. autosummary::
   :toctree: generated/

   sin
   cos
   tan
   arcsin
   arccos
   arctan
   hypot
   arctan2
   degrees
   radians
   unwrap
   deg2rad
   rad2deg


Hyperbolic functions
--------------------

.. autosummary::
   :toctree: generated/

   sinh
   cosh
   tanh
   arcsinh
   arccosh
   arctanh


Rounding
--------

.. autosummary::
   :toctree: generated/

   around
   round_
   rint
   fix
   floor
   ceil
   trunc


Sums, products, differences
---------------------------

.. autosummary::
   :toctree: generated/

   prod
   sum
   nanprod
   nansum
   cumprod
   cumsum
   nancumprod
   nancumsum
   diff
   gradient
   ediff1d
   cross
   trapz


Exponents and logarithms
------------------------

.. autosummary::
   :toctree: generated/

   exp
   expm1
   exp2
   log
   log10
   log2
   log1p
   logaddexp
   logaddexp2


Other special functions
-----------------------

.. autosummary::
   :toctree: generated/

   i0
   sinc


Floating point routines
-----------------------

.. autosummary::
   :toctree: generated/

   signbit
   copysign
   frexp
   ldexp
   nextafter


Rational routines
-----------------

.. autosummary::
   :toctree: generated/

   lcm
   gcd


Arithmetic operations
---------------------

.. autosummary::
   :toctree: generated/

   add
   reciprocal
   positive
   negative
   multiply
   divide
   power
   subtract
   true_divide
   floor_divide
   fmod
   mod
   modf
   remainder
   divmod


Handling complex numbers
------------------------

.. autosummary::
   :toctree: generated/

   angle
   real
   imag
   conj
   conjugate


Miscellaneous
-------------

.. autosummary::
   :toctree: generated/

   convolve
   clip
   sqrt
   cbrt
   square
   absolute
   fabs
   sign
   maximum
   minimum
   fmax
   fmin
   nan_to_num
   interp
.. module:: cupyx.scipy.fftpack

Legacy discrete fourier transforms (:mod:`cupyx.scipy.fftpack`)
===============================================================

.. note::

   As of SciPy version 1.4.0, :mod:`scipy.fft` is recommended over
   :mod:`scipy.fftpack`. Consider using :mod:`cupyx.scipy.fft` instead.

.. Hint:: `SciPy API Reference: Legacy discrete Fourier transforms (scipy.fftpack) <https://docs.scipy.org/doc/scipy/reference/fftpack.html>`_

Fast Fourier Transforms (FFTs)
------------------------------

.. autosummary::
   :toctree: generated/

   fft
   ifft
   fft2
   ifft2
   fftn
   ifftn
   rfft
   irfft
   get_fft_plan


Code compatibility features
---------------------------
1. As with other FFT modules in CuPy, FFT functions in this module can take advantage of an existing cuFFT plan (returned by :func:`~cupyx.scipy.fftpack.get_fft_plan`) to accelarate the computation. The plan can be either passed in explicitly via the ``plan`` argument or used as a context manager. The argument ``plan`` is currently experimental and the interface may be changed in the future version. The :func:`~cupyx.scipy.fftpack.get_fft_plan` function has no counterpart in ``scipy.fftpack``.

2. The boolean switch :data:`cupy.fft.config.enable_nd_planning` also affects the FFT functions in this module, see :doc:`./fft`. This switch is neglected when planning manually using :func:`~cupyx.scipy.fftpack.get_fft_plan`.

3. Like in ``scipy.fftpack``, all FFT functions in this module have an optional argument ``overwrite_x`` (default is ``False``), which has the same semantics as in ``scipy.fftpack``: when it is set to ``True``, the input array ``x`` *can* (not *will*) be overwritten arbitrarily. For this reason, when an in-place FFT is desired, the user should always reassign the input in the following manner: ``x = cupyx.scipy.fftpack.fft(x, ..., overwrite_x=True, ...)``.

4. The boolean switch :data:`cupy.fft.config.use_multi_gpus` also affects the FFT functions in this module, see :doc:`./fft`. Moreover, this switch is *honored* when planning manually using :func:`~cupyx.scipy.fftpack.get_fft_plan`.
Sorting, searching, and counting
================================

.. Hint:: `NumPy API Reference: Sorting, searching, and counting <https://numpy.org/doc/stable/reference/routines.sort.html>`_

.. currentmodule:: cupy

Sorting
-------

.. autosummary::
   :toctree: generated/

   sort
   lexsort
   argsort
   msort
   sort_complex
   partition
   argpartition

.. seealso:: :func:`cupy.ndarray.sort`

Searching
---------

.. autosummary::
   :toctree: generated/

   argmax
   nanargmax
   argmin
   nanargmin
   argwhere
   nonzero
   flatnonzero
   where
   searchsorted
   extract

Counting
--------

.. autosummary::
   :toctree: generated/

   count_nonzero
:orphan:

This document has been moved to :doc:`scipy_ndimage`.
.. module:: cupyx.scipy.sparse

Sparse matrices (:mod:`cupyx.scipy.sparse`)
===========================================

.. Hint:: `SciPy API Reference: Sparse matrices (scipy.sparse) <https://docs.scipy.org/doc/scipy/reference/sparse.html>`_

CuPy supports sparse matrices using `cuSPARSE <https://developer.nvidia.com/cusparse>`_.
These matrices have the same interfaces of `SciPy's sparse matrices <https://docs.scipy.org/doc/scipy/reference/sparse.html>`_.

Conversion to/from SciPy sparse matrices
----------------------------------------

``cupyx.scipy.sparse.*_matrix`` and ``scipy.sparse.*_matrix`` are not implicitly convertible to each other.
That means, SciPy functions cannot take ``cupyx.scipy.sparse.*_matrix`` objects as inputs, and vice versa.

- To convert SciPy sparse matrices to CuPy, pass it to the constructor of each CuPy sparse matrix class.
- To convert CuPy sparse matrices to SciPy, use :func:`get <cupyx.scipy.sparse.spmatrix.get>` method of each CuPy sparse matrix class.

Note that converting between CuPy and SciPy incurs data transfer between
the host (CPU) device and the GPU device, which is costly in terms of performance.

Conversion to/from CuPy ndarrays
--------------------------------

- To convert CuPy ndarray to CuPy sparse matrices, pass it to the constructor of each CuPy sparse matrix class.
- To convert CuPy sparse matrices to CuPy ndarray, use ``toarray`` of each CuPy sparse matrix instance (e.g., :func:`cupyx.scipy.sparse.csr_matrix.toarray`).

Converting between CuPy ndarray and CuPy sparse matrices does not incur data transfer; it is copied inside the GPU device.

Contents
--------

Sparse matrix classes
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   coo_matrix
   csc_matrix
   csr_matrix
   dia_matrix
   spmatrix


Functions
~~~~~~~~~

Building sparse matrices:

.. autosummary::
   :toctree: generated/

   eye
   identity
   kron
   kronsum
   diags
   spdiags
   tril
   triu
   bmat
   hstack
   vstack
   rand
   random


Sparse matrix tools:

.. autosummary::
   :toctree: generated/

   find

Identifying sparse matrices:

.. autosummary::
   :toctree: generated/

   issparse
   isspmatrix
   isspmatrix_csc
   isspmatrix_csr
   isspmatrix_coo
   isspmatrix_dia


Submodules
~~~~~~~~~~

.. autosummary::

   csgraph - Compressed sparse graph routines
   linalg - Sparse linear algebra routines

Exceptions
~~~~~~~~~~

* :class:`scipy.sparse.SparseEfficiencyWarning`
* :class:`scipy.sparse.SparseWarning`
Array manipulation routines
===========================

.. Hint:: `NumPy API Reference: Array manipulation routines <https://numpy.org/doc/stable/reference/routines.array-manipulation.html>`_

.. currentmodule:: cupy

Basic operations
----------------

.. autosummary::
   :toctree: generated/

   copyto
   shape


Changing array shape
--------------------

.. autosummary::
   :toctree: generated/

   reshape
   ravel

.. seealso:: :attr:`cupy.ndarray.flat` and :func:`cupy.ndarray.flatten`

Transpose-like operations
-------------------------

.. autosummary::
   :toctree: generated/

   moveaxis
   rollaxis
   swapaxes
   transpose

.. seealso:: :attr:`cupy.ndarray.T`

Changing number of dimensions
-----------------------------

.. autosummary::
   :toctree: generated/

   atleast_1d
   atleast_2d
   atleast_3d
   broadcast
   broadcast_to
   broadcast_arrays
   expand_dims
   squeeze


Changing kind of array
----------------------

.. autosummary::
   :toctree: generated/

   asarray
   asanyarray
   asfarray
   asfortranarray
   ascontiguousarray
   asarray_chkfinite
   require


Joining arrays
--------------

.. autosummary::
   :toctree: generated/

   concatenate
   stack
   vstack
   hstack
   dstack
   column_stack
   row_stack


Splitting arrays
----------------

.. autosummary::
   :toctree: generated/

   split
   array_split
   dsplit
   hsplit
   vsplit


Tiling arrays
-------------

.. autosummary::
   :toctree: generated/

   tile
   repeat


Adding and removing elements
----------------------------

.. autosummary::
   :toctree: generated/

   append
   resize
   unique
   trim_zeros


Rearranging elements
--------------------

.. autosummary::
   :toctree: generated/

   flip
   fliplr
   flipud
   reshape
   roll
   rot90
