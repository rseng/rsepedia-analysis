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
