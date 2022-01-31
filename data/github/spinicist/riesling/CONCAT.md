![Logo](riesling-logo.png)

[![Build](https://github.com/spinicist/riesling/workflows/Build/badge.svg)](https://github.com/spinicist/riesling/actions)
[![DOI](https://zenodo.org/badge/317237623.svg)](https://zenodo.org/badge/latestdoi/317237623)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03500/status.svg)](https://doi.org/10.21105/joss.03500)

## Radial Interstices Enable Speedy Low-Volume imagING

This is a reconstruction toolbox optimised for 3D ZTE MR images. There are many high quality MR recon toolboxes available, e.g. [BART](http://mrirecon.github.io/bart/), but these are mostly optimised for 2D sequences. 3D non-cartesian sequences present unique challenges for efficient reconstruction, so we wrote our own.

This toolbox was presented at ISMRM 2020.

## Authors

Tobias C Wood, Emil Ljungberg, Florian Wiesinger.

## Installation

Pre-compiled executables are provided for Linux and Mac OS X in a .tar.gz 
archive from http://github.com/spinicist/riesling/releases. Download the 
archive and extract it with `tar -xzf riesling-platform.tar.gz`. Then, move the 
resulting `riesling` executable to somewhere on your `$PATH`, for instance 
`/usr/local/bin`. That's it.

- MacOS Catalina or higher users should use `curl` to download the binary, i.e. 
  `curl -L https://github.com/spinicist/riesling/releases/download/v1.0/riesling-macos.tar.gz`. 
  This is because Safari now sets the quarantine attribute of all downloads, 
  which prevents them being run as the binary is unsigned. It is possible to 
  remove the quarantine flag with `xattr`, but downloading with `curl` is more 
  straightforward.
- The Linux executable is compiled on Ubuntu 16.04 with GLIBC version 2.3 and a 
  statically linked libc++. This means it will hopefully run on most modern 
  Linux distributions. Let us know if it doesn't.

## Compilation

If you wish to compile RIESLING yourself, compilation should hopefully be 
straightforward as long as you have access to a C++17 compiler (GCC 8 or higher,
Clang 7 or higher). RIESLING relies on `vcpkg` for dependency management. To 
download and compile RIESLING, follow these steps:

### 0. MacOS Dependencies
Install the [MacOS vcpkg dependencies](https://github.com/microsoft/vcpkg#installing-macos-developer-tools). This includes:
1. XCode from the AppStore
2. Run `$ xcode-select --install` in the terminal

You may also need to install `pkg-config` depending on your macOS version. This is easily installed with [Homebrew](https://brew.sh/) using
```
$ brew install pkg-config
```

**_(17/01/2022): The Apple Silicon Architecture (M1) is currently not supported. We are working on it though, stay tuned!)_**

### 0. Linux Dependencies
Install the [Linux vcpkg dependencies](https://github.com/microsoft/vcpkg#installing-linux-developer-tools). These includes:

These are `cmake`, `tar`, `curl`, `zip`, `unzip`, `pkg-config` \& `build-essentials`. You may be surprised by which distributions do not include these by default.

### 1. Clone repository
```
$ git clone https://github.com/spinicist/riesling
```

### 2. Compile
In the `riesling` folder execute
```
$ ./bootstraph.sh
```

## Usage

RIESLING comes as a single executable file with multiple commands, similar to 
`git` or `bart`. Type `riesling` to see a list of all the available commands. If you run a RIESLING command without any additional parameter RIESLING will output all available options for the given command.

RIESLING uses HDF5 (.h5) files but can also output NIFTI (.nii). To create an 
example digital phantom, use `riesling phantom`. RIESLING will append suffixes 
to input filenames when writing outputs to indicate which command was executed.

There are several reconstruction algorithms currently provided in RIESLING. 
Simple non-iterative recon is available with (`riesling recon`).

A separate examples repository https://github.com/spinicist/riesling-examples
contains Jupyter notebooks demonstrating most functionality. These can also be
run on the mybinder website.

## Documentation & Help

Documentation is available at https://riesling.readthedocs.io.

If you can't find an answer there or in the help strings, 
you can open an [issue](https://github.com/spinicist/riesling/issues), or find
the developers on Twitter ([@spinicist](https://twitter.com/spinicist)) or
e-mail tobias.wood@kcl.ac.uk.
# Contributing to RIESLING

Welcome to RIESLING. Any help is very much appreciated.

These guidelines exist to make it easy to get involved and ensure that RIESLING is welcoming to everyone.

## Asking Questions / Reporting Bugs

If you have general questions about how to use RIESLING programs, you can either open an [issue](https://github.com/spinicist/riesling/issues).

If you think you have found a bug, then opening a Github issue is the preferred avenue. Please check whether an identical or similar issue already exists first. When opening a new issue, please give as much information as you can, including
- The version of riesling you are using. Either give the version number if you downloaded the binaries, or a branch/git commit id if you compiled from source.
- The operating system you are running on.
- A description of the bug, including the full shell command you ran and any input files.

## Contributing Changes

If you want to edit the RIESLING code yourself to fix a problem or implement a new feature, you are very welcome to! Please follow this model for submitting a Pull Request:
- (Fork)[https://help.github.com/articles/fork-a-repo/] the riesling repo to your Github profile.
- Clone this repo, then create a new branch.
- Make your changes on this branch.
- Run clang-format. Riesling has a .clang-format file to avoid wasting time arguing over style.
- Squash/rebase your changes into a minimal number of commits.
- Make sure your copy of the `main` branch is up to date, and if necessary `rebase` your branch to the latest commit.
- Push your branch to Github.
- Open a Pull Request.

## Code of Conduct

- Anyone who participates in the development of RIESLING is expected to show respect and courtesy to all other community members at all times.
- Harrasment in any form towards any members will not be tolerated.
- All communication should be appropraite for a professional audience including people of different backgrounds.
- Be kind to others. Do not insult or put down other contributers.

## Thanks

Above all, thank you for using or contributing to RIESLING. It really is appreciated.
---
title: "Radial Interstices Enable Speedy Low-volume Imaging"
tags:
  - mri
  - reconstruction
  - ZTE
  - cpp
authors:
  - name: Tobias C Wood
    orcid: 0000-0001-7640-5520
    affiliation: 1
  - name: Emil Ljungberg
    orcid: 0000-0003-1456-7967
    affiliation: 1
  - name: Florian Wiesinger
    orcid: 0000-0002-5597-6057
    affiliation: "1,2"
affiliations:
  - name: Department of Neuroimaging, King's College London
    index: 1
  - name: GE Healthcare
    index: 2
date: 2021-06-25
bibliography: paper.bib
---

# Summary

- Standard MRI methods acquire Fourier-encoded data on a regularly spaced Cartesian grid. Noncartesian MRI methods abandon this grid and instead acquire data along arbitrary trajectories in k-space, which can lead to advantages such as motion robustness and more flexible pulse sequence design.
- Zero Echo Time (ZTE) imaging is a specialised form of MRI with an inherently noncartesian three-dimensional radial trajectory [@Ljungberg]. ZTE imaging has many interesting benefits, such as near silent operation, but also technical challenges such as dead-time gap artefacts [@Froidevaux].
- We have developed a toolbox, named Radial Interstices Enable Speedy Low-volume imagING (RIESLING), tuned for high performance reconstruction of MRI data acquired with 3D noncartesian k-space trajectories. While our group has a focus on ZTE imaging, RIESLING is suitable for all 3D noncartesian trajectories.

# Statement of Need

3D noncartesian trajectories can be challenging to reconstruct compared to other MRI trajectories. The major issues are:

- **Inseparability of the trajectory dimensions.** In Cartesian imaging it is often possible to consider one of the image dimensions separately to the others, leading to algorithmic complexity reduction and memory savings. 3D noncartesian must consider all dimensions simultaneously.
- **Oversampling requirements.** In order to obtain good image quality, it is necessary to oversample the reconstruction grid. Coupled with the inseparability of the 3D problem, this leads to large memory requirements.
- **Sample density compensation.** Noncartesian trajectories can lead to uneven k-space sampling, the correction of which is equivalent to preconditioning a linear system [@cgSENSE]. For 2D noncartesian trajectories, it is generally not necessary to correct for this in the context of an iterative reconstruction, as the problem will still converge. For 3D noncartesian trajectories, such preconditioning is essential to ensure reasonable convergence properties.

Existing MRI reconstruction toolboxes such as BART [@BART] or SigPy [@SigPy] provide high-quality implementations of many reconstruction algorithms but have prioritised applications other than 3D noncartesian reconstruction. We hence created a dedicated toolbox for this, including specific features such as:

- A configurable grid oversampling factor
- A thread-safe gridding implementation suitable for multi-core CPUs
- Integrated sample density compensation functions [@Zwart]

RIESLING (http://github.com/spinicist/riesling) is written with a modern C++ toolchain and utilizes Eigen [@Eigen] for all core operations, providing high performance. Data is stored in the HDF5 format allowing easy integration with other tools. We provide multiple reconstruction strategies, including classic conjugate-gradient SENSE [@cgSENSE] and Total-Generalized Variation [@TGV]. We are actively using RIESLING for our own studies and hope that it will be of interest to other groups using such sequences.

# References
# MATLAB Integration
Includes two main functions for reading and writing riesling k-space data

```
varargout = riesling_read(fname);
```

```
riesling_write(fname, data, traj, info, meta);
```

## Info struct

You can create an empty header info struct with

info = riesling_info();

## Examples
One examples is included:

1. `riesling_matlab_demo.m`: Demonstrating basic data input, output, manipulation and plotting.
.. RIESLING documentation master file, created by
   sphinx-quickstart on Sun May  9 11:26:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RIESLING!
====================================

.. image:: riesling-logo.png
   :alt: RIESLING Logo

RIESLING is a reconstruction toolbox tuned for 3D radial MRI - the name stands for Radial Interstices Enable Speedy Low-volume imagING - but it can process other trajectories as well.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   docs/contents

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Notes for Developers
====================

These notes are intended to highlight to give a brief overview of the RIESLING code. They are not intended as a comprehensive API reference, but may evolve into that over time.

Unit Tests
----------

RIESLING includes some unit tests which cover the core functionality. They do not cover everything, and should likely be expanded. The tests use the `Catch2 <https://github.com/catchorg/Catch2>`_ framework.

The tests are built with the ``riesling-tests`` target, which is included in the ``all`` target. ``riesling-tests`` does not have an install location set, so it will be built but remain in the build directory if the ``install`` target is used.

To run the tests, simply run ``./riesling-tests`` within the build directory. By default, all cases are tested. A list of individual test cases can be obtained with ``riesling-tests --list-tests``.

The test case code is located in ``/test``. Generally the files are named after the corresponding file in ``/src``, unless it makes sense to sub-divide the tests further.

The ultimate set of tests are the various notebooks contained in the `riesling-examples` repository.

Code Structure
--------------

RIESLING is built as a monolithic binary. There is a simple ``main.cpp`` file which specifies the available commands, each of which is contained in a ``main_command.cpp`` file. The commands can roughly be split into simple utilities (such as ``riesling hdr``) and the more complicated commands that perform a particular reconstruction method (such as ``riesling cg``).

The most important command to explain is ``riesling cg`` / ``main_cg.cpp``. The top of the ``main_cg`` function is fairly straight-forward - a set of flags is declared, including those common across recon methods which are defined by a macro defined in ``parse_args.h``. Once the command-line has been parsed, the trajectory and ``Info`` header struct are read from the input file.

The next few lines initialise the gridding kernel and the ``GridOp`` object. The ``Trajectory`` and ``GridOp`` objects are the work-horses of RIESLING. ``Trajectory`` can calculate an efficient ``Mapping`` between non-Cartesian and Cartesian co-ordinates. This consists of lists of matching integer non-Cartesian and Cartesian co-ordinates, the floating-point offset from the Cartesian grid-point, and a list of indices sorted by Cartesian grid location. This ``Mapping`` depends on the chosen over-sampling factor. This enables fast, thread-safe, non-Cartesian to Cartesian gridding as each thread can work on a section of the Cartesian grid without conflicting writes. The ``Mapping`` is used to construct the ``GridOp`` object, which contains the interpolation code.

After the ``GridOp`` is constructed, the Sample Density Compensation is either calculated or loaded from a file on disk. Then the necessary cropping between the oversampled reconstruction grid and the output image is calculated, followed by the apodization required to correct for any apodization introduced by the gridding kernel. Next the required FFT for the reconstruction grid is planned using the FFTW library.

The next key step is calculating the sensitivity maps, or loading them from disk. Calculating them currently requires constructing a second ``GridOp`` object within the ``DirectSENSE`` function. This uses a very low resolution, so constructing and sorting the grid co-ordinates is usually fast.

If a basis has been specified for a low-rank/subspace reconstruction, then the next step is to replace the GridOp with one including the basis.

At this point we can create the ReconOp object which combines gridding, FFT, and SENSE channel combination. This object is then passed to the iterative optimization routines. If the user wants, Töplitz embedding can be enabled.

We then loop through all the volumes in the input file, load the data, apply the adjoint operator to get a starting image, and then apply conjugate gradients to optimize the image. We crop to the desired output FoV and apply a Tukey filter if requested by the user. Finally the images are written to disk.
Reconstruction
==============

There are currently three different reconstruction commands provided in RIESLING - ``recon``, ``cg`` and ``tgv``. More may be added in future.

Non-iterative
-------------

The ``recon`` command provides basic non-iterative reconstructions. This is useful when you want to run a quick reconstruction to ensure that the data file is in the correct format, but is unlikely to yield optimal image quality. However it will be useful to describe some of the command options here, because they are shared with ``cg`` and ``tgv``.

By default ``recon`` will output complex-valued images in ``.h5`` format that have been combined using channel sensitivities extracted from the scan. If you are only interested in the magnitude images, add ``--mag`` to the command-line. If you only want a root-sum-squares reconstruction, add ``--rss``. If you would like NIFTI images as output, add ``--oft=nii``.

Non-cartesian MRI data requires a Non-Uniform FFT (NUFFT) instead of a simple FFT for conversion between k-space and image space. The NUFFT consists of a gridding step and then a normal FFT [1, 2]. The gridding step in RIESLING is controlled by the oversampling and kernel options.

``--os=X`` controls the grid oversampling ratio. The default value is 2, but this is very memory intensive (as it is a 3D grid, 2x oversampling requires 8 times the memory of the native matrix size). Reducing the over-sampling to a value of 1.3 leads to a reduced memory foot-print for little impact in image quality [3]. However, oversampling factors below 2 do not work with Töplitz embedding.

``--kernel=KB3`` selects a width 3 Kaiser-Bessel (KB) kernel instead of the default nearest-neighbour (NN)kernel [3, 4]. Kaiser-Bessel is the default in most other toolboxes and may become the default in RIESLING in the future. KB5 is also available.

The final gridding option is ``--fast-grid``. RIESLING uses a parallelized gridding operation. In order to avoid multiple threads writing to the same Cartesian k-space point, the trajectory is sorted by Cartesian k-space location and each thread processes its own chunk of the sorted co-ordinates. However, there may still be race conditions at the edge of each chunk, particularly for small images with highly oversampled k-space centers. To prevent this, each thread uses its own workspace, and then these are combined into the final grid at the end. This process is thread-safe but doubles the memory requirements for the gridding operations. The ``--fast-grid`` option makes the threads write directly into the final grid, reducing memory consumption but at the risk of conflicting writes into the grid. When gridding high resolution images on a small number of threads, e.g. 10 or fewer, the probability of a race condition is vanishingly small. Use at your own risk.

The gridding step can also compensate for the increased density of samples in the oversampled central k-space region with most non-Cartesian trajectories. This step is often omitted in 2D non-Cartesian iterative reconstructions, but is essential for reasonable convergence in 3D non-Cartesian reconstruction. The ``--sdc`` option controls the sample density compensation method - valid values are "none" to turn it off, "pipe" for the iterative method of Pipe et al [5], and "pipenn" for an approximate but fast version of the Pipe method. In the ``riesling sdc`` command analytic radial weights are also available. The default is "pipenn". It is also possible to pre-calculate the densities of a given trajectory using the ``riesling sdc`` command and then pass in the resulting file, i.e. ``--sdc=file-sdc.h5``. The weighting of the compensation can be reduced using the ``--sdcPow`` option - the power is applied to all k-space densities equally and values between 0 and 1 make sense.

Due to the oversampled central region in most non-Cartesian trajectories reasonable channel sensitivities can be extracted directly from the data. This is the default option. Tikhonov regularization can be applied to the sensitivities using the ``--lambda`` option (the value should be approximately the same as the background noise in the ``--rss`` reconstruction). In a multi-volume reconstruction, the sensitivities are taken from the last volume by default but can be specified using ``--senseVolume``, or can be taken from an external file using ``--sense``.

You can apply a basic Tukey filter to the final image k-space using ``--tukey_start``, ``--tukey_end`` and ``--tukey_height``. The start and end are defined as the fractional radius in k-space, i.e. 0.5 and 1. The height option is specified at the end radius and should be between 0 (heavy filtering) and 1 (no filtering). Finally, if you want to expand (or contract) the field-of-view of an image, for instance with a read-oversampled acquisition, then use the ``--fov`` option.

Iterative
---------

The workhorse reconstruction tool in RIESLING is ``cg``, which runs an un-regularized cgSENSE reconstruction. For speed RIESLING uses a Töplitz embedding strategy [6]. This uses the gridding method to calculate the k-space transfer function on the Cartesian grid. After the initial gridding from non-Cartesian to Cartesian grid, each iteration only requires SENSE combination/expansion, the forwards/reverse FFT, and a multiplication of Cartesian k-space by the transfer function.

The additional options added for ``cg`` control the iterations strategy. ``--max_its`` specifies the maximum number of iterations. cgSENSE image quality often benefits from early stopping of the iterations, which is an implicit form of regularization as it prevents the algorithm from over-fitting noise. The default value is 8, with correct density compensation reasonable images can often be obtained in only 4. You can also specify a threshold to terminate the iterations using ``--thresh``. The default value is 1e-10 which is very strict and rarely reached. Values of ``1e-3`` or so would lead to early stopping.

Finally, ``cg`` adds an additional ``--iter_fov`` option which controls the field-of-view cropping used during the iterations. This needs to be larger than the final FOV to avoid aliasing and edge effects. The default value is 256 mm which is sufficient for most brain reconstructions. Note that if you pre-compute sensitivities, their FOV must match this value.

The ``riesling admm`` command uses the Alternating Directions Method-of-Multipliers, also known as an Augmented Lagrangian method, to add regularizers to the cgSENSE method. Currently the only regularizer available is Locally Low-Rank. Additional options are available to: control the number of outer iterations ``--admm_its``, the regularization strength ``--reg``, the patch size for LLR ``--patch``, and how tightly coupled the regularizer and data-fidelity terms are ``--rho`` [7].

The ``tgv`` command uses Total Generalized Variation regularization to improve image quality [8]. It uses a different optimization algorithm to ``cg`` which is noticeable slower, but still reasonable. It adds three more options. ``--alpha`` controls the initial regularization level. The default is 1e-5, better results can often be obtained with 2e-5. ``--reduce`` will reduce the regularization over the course of the iterations, which can prevent over-smoothing. ``--step`` controls the gradient-descent step size and is specified as an inverse, i.e. a value of 8 results in a step-size of 1/8th the gradient. Smaller values (larger step sizes) give faster convergence but can lead to artefacts.

References
----------

1. Pruessmann, K. P., Weiger, M., Börnert, P. & Boesiger, P. Advances in sensitivity encoding with arbitrary k-space trajectories. Magn. Reson. Med. 46, 638–651 (2001).
2. JI Jackson, C. H. Meyer, D. G. Nishimura, and A. Macovski, ‘Selection of a convolution function for Fourier inversion using gridding (computerised tomography application)’, IEEE Transactions on Medical Imaging, vol. 10, no. 3, pp. 473–478, Sep. 1991, doi: 10.1109/42.97598.
3. Beatty, P. J., Nishimura, D. G. & Pauly, J. M. Rapid gridding reconstruction with a minimal oversampling ratio. IEEE Transactions on Medical Imaging 24, 799–808 (2005).
4. Oesterle, C., Markl, M., Strecker, R., Kraemer, F. M. & Hennig, J. Spiral reconstruction by regridding to a large rectilinear matrix: A practical solution for routine systems. Journal of Magnetic Resonance Imaging 10, 84–92 (1999).
5. Zwart, N. R., Johnson, K. O. & Pipe, J. G. Efficient sample density estimation by combining gridding and an optimized kernel: Efficient Sample Density Estimation. Magn. Reson. Med. 67, 701–710 (2012).
6. CA Baron, N. Dwork, J. M. Pauly, and D. G. Nishimura, ‘Rapid compressed sensing reconstruction of 3D non-Cartesian MRI’, Magnetic Resonance in Medicine, vol. 79, no. 5, pp. 2685–2692, May 2018, doi: 10.1002/mrm.26928.
7. J. I. Tamir et al., ‘T2 shuffling: Sharp, multicontrast, volumetric fast spin‐echo imaging’, vol. 77, pp. 180–195, 2017.
8. Knoll, F., Bredies, K., Pock, T. & Stollberger, R. Second order total generalized variation (TGV) for MRI. Magnetic Resonance in Medicine 65, 480–491 (2011).
Data Format
===========

RIESLING uses `HDF5 <https://www.hdfgroup.org/solutions/hdf5>`_ to store input and intermediate data. Output data is by default written to HDF5 (``.h5``), but optionally can also be output to NiFTI (``.nii``).

RIESLING mandates that the non-cartesian data is stored in "spokes", which could equally be called frames or traces. Here we will treat the data as being stored in `S` spokes, with `N` data-points per spoke, and `C` k-space channels.

Header
------

To be considered valid RIESLING input, the HDF5 file must contain the header information datastructure, stored as a compound data-type in a dataset with the name ``info``. We reserve the right to change these fields of the header structure between versions of RIESLING. For the canonical definition of the header, see ``src/info.h``. A pseudo-code version of the header is given here for clarity:

.. code-block:: c

  struct Info {
    long type;
    long matrix[3];

    long channels;
    long read_points;
    long spokes;

    long volumes;
    long echoes;

    float tr;
    float voxel_size[3];
    float origin[3];
    float direction[3][3];
  };

* ``type`` defines the kind of acquisition. Currently two values are supported - 1 means the acquisition is fully 3D, while 2 means the acquisition is a 3D stack-of-stars or stack-of-spirals type acquisition, with cartesian phase-encoding blips for the third axis.
* ``matrix`` defines the nominal matrix size for the scan - i.e. it determines the matrix size of the final reconstructed image (unless the `--fov` option is used).
* ``channels`` defines the number of k-space channels / coil-elements.
* ``read_points`` sets the number of data-points in the readout direction, i.e. how many readout points per spoke.
* ``spokes`` sets the number of spokes in the non-cartesian k-space acquisition.
* ``volumes`` indicates how many volumes or time-points were acquired in the acquisition.
* ``echoes`` specifies how many separate echoes were acquired per time-point.

The final four fields specify the TR and image orientation as required to build a valid NIfTI output file.

* ``tr`` The repetition time. Should be specified in milliseconds as per NIfTI convention.
* ``voxel_size`` The nominal voxel-size. Should be specified in millimeters as per NIfTI/ITK convention.
* ``origin`` The physical-space location of the center of the voxel at index 0,0,0, as per ITK convention.
* ``direction`` The physical-space axes directions, as per ITK convention.

Trajectory
----------

The trajectory should be stored as a float array in a dataset with the name ``trajectory`` with dimensions ``SxNx3``. HDF5 uses a row-major convention, if your software is column major (RIESLING is internally) then this will be ``3xNxS``. The 3 co-ordinates correspond to the x, y & z locations within the k-space volume. For a full 3D acquisition these should be scaled such that the nominal edge of k-space in each direction is 0.5. Hence, for radial spokes the k-space locations go between 0 and 0.5, and for diameter spokes between -0.5 and 0.5. For a 3D stack trajectory, the z co-ordinate should be the slice/stack position.

Echoes
------

If the dataset contains multiple echoes, or other temporal points (e.g. cardiac or respiratory phases) which should be reconstructed together, then an additional dataset should be added to the input H5 file called ``echoes``. This should be an integer valued, one-dimensional array with the number of entries equal to the number of spokes specified in the ``info`` structure.

Non-cartesian Data
------------------

The non-cartesian data itself must be stored in a complex-valued float-precision dataset named ``noncartesian`` with dimensions ``VxSxNxC`` where V is the number of volumes. HDF5 does not have a native complex-valued datatype, hence a compound datatype with a ``r`` and ``i`` members corresponding to the real and imaginary parts must be used. In contrast to other toolkits RIESLING stores the channels as the fastest-varying index, i.e. the data for each k-space point across all channels is stored contiguously.

Cartesian Data
--------------

The ``riesling grid`` command will produce a complex-valued dataset named ``cartesian`` containing the gridded cartesian data for all channels. The dimensions will depend on the reconstruction settings (notably the oversampling factor).

Image Data
----------

The output of a reconstruction command will write a complex-valued dataset named ``image``, unless the ``--mag`` command is specified in which case the dataset will be real-valued. The dimensions will be ``VxZxYxXxE`` where V is the number of volumes, X, Y & Z are the matrix size, and E is either the number of echoes or the number of basis-vectors if a low-rank reconstruction has been used.

Density Compensation
--------------------

``riesling sdc`` pre-calculates Sample Density Correction factors. It produces a real-valued dataset ``sdc`` of dimension ``SxN``.

Meta-Information
----------------

RIESLING is capable of storing additional meta-information and passing it through the processing chain. This should be stored in an HDF5 group named ``meta``, and consist of key-value pairs where the key is the dataset name and the value is a single floating-point number.
Utilities
=========

RIESLING includes a number of utilities that either implement functionality that is not strictly reconstruction or assist with basic file management.

Sensitivities
-------------

Estimating channel sensitivities is a key step in modern reconstruction methods. The ``sense`` command allows you to run the direct sensitivity extraction step that is incorporated into the main reconstruction commands and save the results, either for inspection or for use across multiple images. As described in :doc:`recon`, Tikhonov regularization can be applied with ``--lambda``. You must set the ``--fov`` option to match the FOV that the channels will be combined at - for the iterative recon methods this is equal to the ``--iter_fov`` value, which defaults to 256 mm.

``riesling espirit`` implements the popular ESPIRiT method for estimating channel sensitivities that exploits the correlations between channels in k-space. Important options are ``--kRad``, which sets the initial k-space kernel radius and ``--calRad`` which determines the size of the calibration region. RIESLING defines the calibration region as the cubic region with "radius" of twice the kernel radius plus the calibration radius, i.e. the calibration radius is expanded by the kernel width. The default value is one, or one plus the dead-time gap if the data has a gap. ``--thresh`` defines the threshold for retaining kernels after the first SVD step. Any kernel/singular vector with a singular value (as a fraction of the first singular value) above the threshold will be retained.

Pre-calculation
---------------

In addition to pre-calculating sensitivities, it is also possible to pre-calculate sample densities using the ``--sdc`` command. All the reconstruction options (``--kb``, ``--kw``, ``--os``) must match to the settings you will use for the actual reconstruction or you will get artefacts. This command is useful if you will be running many reconstructions with the same trajectory.

The first time you run RIESLING with a new trajectory it will be fairly slow while it "plans" the Fourier transforms. This is a feature of the FFTW library that RIESLING uses internally. In short, FFTW attempts different strategies for any given FFT, measures which is fastest, and then saves this "wisdom" for future use. RIESLING stores these in a file called ``.riesling-wisdom`` in the user's home directory. Similarly to sample densities, the FFT settings can be re-used for any future FFT of the same size. The size of the FFT in RIESLING is controlled by the trajectory and the ``--os`` option. You can run ``riesling plan`` with a particular ``.h5`` file and ``--os`` value to force the planning before you run any real reconstructions.

Compression
-----------

Modern MR scanners are often equipped with multi-channel receive coils with a high number of elements. This increases memory requirements, both on disk and in RAM, and can make steps in the reconstruction ill-conditioned. It is hence advisable to compress the the raw data to a smaller number of virtual channels before running the reconstruction. This step should be carried out first before any subsequent operations, e.g. sensitivity estimation. The ``compress`` command in RIESLING implements basic PCA coil-compression, with the number of output virtual channels specified by the ``--cc`` option. We have not implemented other methods, such as Geometric Coil Compression, because they rely on properties of Cartesian sequences that we cannot rely on for non-Cartesian.

Gridding
--------

``riesling grid`` will carry out only the first step of the NUFFT, i.e. it will grid non-Cartesian k-space to Cartesian (or vice versa) and save the result. This can be useful to check that a dataset has been acquired correctly.

To diagnose trajectory and sample density issues, you can instead use ``riesling traj``. This will apply the gridding and sample density compensation to a set of ones, allowing you to see if the Cartesian k-space has an even weighting.

Simulations
-----------

``riesling phantom`` will produce a simulated image, useful for experimenting with reconstruction settings. Currently spherical and Shepp-Logan phantoms with simple multi-channel coil sensitivities are implemented. It is advisable to use a high grid oversampling rate to minimise rasterization errors. The trajectory can be read from a ``.h5`` file with the ``--traj`` command, otherwise an Archimedean 3D spiral trajectory will be used.

Data Tools
----------

``riesling split`` will split out a single volume from a multi-volume ``.h5`` file, and will separate the low- and high-resolution k-space trajectories if they are present.

Dead-time Gap Filling
---------------------

``riesling zinfandel`` implements an experimental ZTE dead-time gap filling method based on 1D GRAPPA. This will be the subject of a future publication.

References
----------

1. Yeh, E. N. et al. Inherently self-calibrating non-cartesian parallel imaging. Magnetic Resonance in Medicine 54, 1–8 (2005).
2. Uecker, M. et al. ESPIRiT-an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA. Magnetic Resonance in Medicine 71, 990–1001 (2014).
3. Zwart, N. R., Johnson, K. O. & Pipe, J. G. Efficient sample density estimation by combining gridding and an optimized kernel: Efficient Sample Density Estimation. Magn. Reson. Med. 67, 701–710 (2012).
4. Wong, S. T. S. & Roos, M. S. A strategy for sampling on a sphere applied to 3D selective RF pulse design. Magnetic Resonance in Medicine 32, 778–784 (1994).
Contents
========

.. toctree::
   :maxdepth: 2
   
   intro
   recon
   util
   data
   devel
Introduction
============

RIESLING is a tool for reconstructing non-cartesian MRI scans. It has been tuned specifically for reconstructing 3D radial center-out trajectories associated with Zero Echo-Time (ZTE) or Ultrashort Echo-Time (UTE) sequences. These trajectories provide unique challenges in efficient reconstruction compared to Cartesian trajectories (both 2D and 3D).

RIESLING is provided as a single executable file, similar to ``bart``. The ``riesling`` executable provides multiple individual commands, which vary from basic utilities to interrogate files to complete reconstruction pipelines. To see a full list of commands currently available, run ``riesling``. The full list is not repeated here as they are subject to change. However, the core commands are:

- ``riesling hdr`` Prints the header information from compatible ``.h5`` files
- ``riesling recon`` Performs a non-iterative reconstruction with either root-sum-squares (``--rss``) channel combination or a sensitivity-based combination (the default)
- ``riesling sense`` Extract channel sensitivities from a dataset for future use
- ``riesling sdc`` Pre-calculate sample densities
- ``riesling plan`` Pre-plan FFTs
- ``riesling cg`` Iterative conjugate-gradients SENSE reconstruction
- ``riesling admm`` Regularized cgSENSE reconstruction (only Locally-Low-Rank available currently)
- ``riesling tgv`` TGV-regularized iterative reconstruction

A `tutorial notebook <https://github.com/spinicist/riesling-examples/tutorial.ipynb>`_ can be run interactively at on `MyBinder <https://mybinder.org/v2/gh/spinicist/riesling-examples/HEAD?filepath=tutorial.ipynb>`_. This explains the various steps required to generate a simulated phantom dataset and then reconstruct it. You will need to reduce the matrix size to 64 to run with MyBinder's RAM limit.

An important step with using RIESLING is providing data in the correct ``.h5`` format. Details of this format can be found in :doc:`data`. Users of the ZTE sequence on GE platforms should contact the authors to discuss conversion strategies.

Further details about the reconstruction tools can be found in :doc:`recon`.
