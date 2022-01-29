

| **master**  | **dev** |
| ------------- | ------------- |
| [![Build Status](https://travis-ci.org/nlesc-dirac/sagecal.svg?branch=master)](https://travis-ci.org/nlesc-dirac/sagecal)  | [![Build Status](https://travis-ci.org/nlesc-dirac/sagecal.svg?branch=dev)](https://travis-ci.org/nlesc-dirac/sagecal)  |

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1289316.svg)](https://doi.org/10.5281/zenodo.1289316)
[![Documentation](https://codedocs.xyz/nlesc-dirac/sagecal.svg)](https://codedocs.xyz/nlesc-dirac/sagecal/)

# SAGECAL


## Features

- Levenberg-Marquardt, (stochastic) LBFGS, Riemannian Trust Region, Nesterov's accelerated gradient descent algorithms
- GPU acceleration using CUDA
- Fast and accurate interferometric calibration
- Gaussian and Student's t noise models
- Shapelet source models
- CASA MS data format supported
- Distributed calibration using MPI - consensus optimization with data multiplexing
- Tools to build sky models and restore sky models to images
- Adaptive update of ADMM penalty (Barzilai-Borwein a.k.a. Spectral method)
- Bandpass calibration and unprecedented RFI mitigation with stochastic LBFGS
- Stochastic calibration for handling data at highest resolution (with federated averaging and consensus optimization)


Please read INSTALL.md for installation instructions, but 'cmake' should work in most cases. We give a brief guide to use SAGECal here but there is extensive documentation here.

## Contributing
Read the [contributing guide](https://github.com/nlesc-dirac/sagecal/blob/master/CONTRIBUTING.md)


## Code documentation
Code documentation can be found [here](https://codedocs.xyz/nlesc-dirac/sagecal).


## Step by Step Introduction:

### 1) Input Data
Input to sagecal must be in CASA MS format, make sure to create a column in the MS to write output data as well. The data can be in raw or averaged form, also initial calibration using other software can be also applied.


### 2) Sky Model:
#### 2a) Make an image of your MS (using ExCon/casapy). 
Use Duchamp to create a mask for the image. Use buildsky to create a sky model. (see the README file on top level directory). Also create a proper cluster file.
Special options to buildsky: -o 1 (NOTE: not -o 2)

Alternatively, create these files by hand according to the following formats.

#### 2b) Cluster file format:
cluster_id chunk_size source1 source2 ...
e.g.
```

0 1 P0C1 P0C2
2 3 P11C2 P11C1 P13C1

```

Note: putting -ve values for cluster_id will not subtract them from data.
chunk_size: find hybrid solutions during one solve run. Eg. if -t 120 is used 
to select 120 timeslots, cluster 0 will find a solution using the full 120 timeslots while cluster 2 will solve for every 120/3=40 timeslots.

#### 2c) Sky model format:
```
#name h m s d m s I Q U V spectral_index RM extent_X(rad) extent_Y(rad) pos_angle(rad) freq0
```

or

```
#name h m s d m s I Q U V spectral_index1 spectral_index2 spectral_index3 RM extent_X(rad) extent_Y(rad) pos_angle(rad) freq0
```

e.g.:

```
P1C1 0 12 42.996 85 43 21.514 0.030498 0 0 0 -5.713060 0 0 0 0 115039062.0
P5C1 1 18 5.864 85 58 39.755 0.041839 0 0 0 -6.672879 0 0 0 0 115039062.0
#A Gaussian mjor,minor 0.1375,0.0917 deg diameter -> radius(rad), PA 43.4772 deg (-> rad)
#Position Angle: "West from North (counter-clockwise)" (0 deg = North, 90 deg = West). 
#Note: PyBDSM and BBS use "North from East (counter-clockwise)" (0 deg = East, 90 deg = North). 
G0  5 34 31.75 22 00 52.86 100 0 0 0 0.00 0 0.0012  0.0008 -2.329615801 130.0e6
#A Disk radius=0.041 deg
D01 23 23 25.67 58 48 58 80 0 0 0 0 0 0.000715 0.000715 0 130e6
#A Ring radius=0.031 deg
R01 23 23 25.416 58 48 57 70 0 0 0 0 0 0.00052 0.00052 0 130e6
#A shapelet ('S3C61MD.fits.modes' file must be in the current directory)
S3C61MD 2 22 49.796414 86 18 55.913266 0.135 0 0 0 -6.6 0 1 1 0.0 115000000.0
```


Note: Comments starting with a '#' are allowed for both sky model and cluster files.
Note: 3rd order spectral indices are also supported, use -F 1 option in sagecal.
Note: Spectral indices use natural logarithm, ```exp(ln(I0) + p1 * ln(f/f0) + p2 * ln(f/f0)^2 + ..)``` so if you have a model with common logarithms like ```10^(log(J0) + q1*log(f/f0) + q2*log(f/f0)^2 + ..)``` then, conversion is

```
ln(I0)+p1*ln(f/f0)+p2*ln(f/f0)^2+... = ln(10)*(log(J0)+q1*log(f/f0)+q2*log(f/f0))^2)+...)
=ln(10)*(ln(J0)/ln(10)+q1*ln(f/f0)/ln(10)+q2*ln(f/f0)^2/ln(10)^2+...)
```
so
```
I0=J0
p1=q1
p2=q2/ln(10)
p3=q3/ln(10)^2
...
```

### 3) Run sagecal
Optionally: Make sure your machine has (1/2 working NVIDIA GPU cards or Intel Xeon Phi MICs) to use sagecal.
Recommended usage: (with GPUs)

```
sagecal -d my_data.MS -s my_skymodel -c my_clustering -n no.of.threads -t 60 -p my_solutions -e 3 -g 2 -l 10 -m 7 -w 1 -b 1
```

Use your solution interval (-t 60) so that its big enough to get a decent solution and not too big to make the parameters vary too much. (about 20 minutes per solution is reasonable).

Note: It is also possible to calibrate more than one MS together. See section 4 below.
Note: To fully use GPU acceleration use ```-E 1``` option.

Simulations:
With ```-a 1```, only a simulation of the sky model is done.
With ```-a 1``` and ```-p``` 'solutions_file', simulation is done with the sky model corrupted with solutions in 'solutions_file'.
With ```-a 1``` and ```-p``` 'solutions_file' and ```-z``` 'ignore_file', simulation is done with the solutions in the 'solutions_file', but ignoring the cluster ids in the 'ignore_file'.
E.g., If you need to ignore cluster ids '-1', '10', '999', create a text file :

```
-1
10
999
```

and use it as the 'ignore_file'.

Bandpass correction using **stochastic** calibration with consensus:
Use ```-N 1``` combined with options for ```-M```,```-w``` (see also section 4 below).


### 4) Distributed calibration

Use mpirun to run sagecal-mpi, example:
```
 mpirun  -np 11 -hostfile ./machines --map-by node --cpus-per-proc 8 
 --mca yield_when_idle 1 -mca orte_tmpdir_base /scratch/users/sarod 
 /full/path/to/sagecal-mpi -f 'MS*pattern' -A 30 -P 2 -r 5 
 -s sky.txt -c cluster.txt -n 16 -t 1 -e 3 -g 2 -l 10 -m 7 -x 10 -F 1 -j 5
```

Specific options : 
```-np 11``` : 11 processes : starts 10 workers + 1 master

```./machines``` : will list the host names of the 11 (or fewer) nodes used ( 1st name is the master ) : normally the node where you invoke mpirun

```-f 'MS*pattern'``` : Search MS names that match this pattern and calibrate all of them together. The total number of MS being calibrated can be higher than the actual number of slaves (multiplexing).

```-A 30``` : 30 ADMM iterations.

```-P 2``` : polynomial in frequency has 2 terms.

```-Q``` : can change the type of polynomial used (```-Q 2``` gives Bernstein polynomials).

```-r 5``` : regularization factor is 5.0.

```-G textfile```: each cluster can have a different regularization factor, instead of using ```-r``` option when the regularization is the same for all clusters.

```-N 1```: enable **stochastic** calibration (minibatches of data), combined with options ```-M```, ```-w``` and ```-u```.

```-U 1```: use global solution instead of local solution for residual calculation.

MPI specific options:

```/scratch/users/sarod``` : this is where MPI stores temp files (default is probably ```/tmp```).

```--mca*```: various options to tune the networking and scheduling.

Note: the number of workers (-np option) can be lower than the number of MS calibrated. The program will divide the workload among the number of available workers.


The rest of the options are similar to sagecal.

### 5) Spatial regularization
Spatial regularization (with distributed multi-directional calibration) enables the use of solutions along all directions in the sky as spatial regularization for calibration along each single direction. In other words, we construct a model of the systematic errors covering the full field of view (not explicitly but implicitly) and use this as an additional constraint. To turn this on, the following options are needed:

```-X L2 penalty,L1 penalty,model order,FISTA iterations,update cadence``` : For example ```-X 0.01,1e-4,3,40,10``` will use 0.01 as the L2 penalty and 1e-4 as the L1 penalty for solving an elastic net regression to build a spatial model. The model will have 3x3=9 spatial basis functions. The elastic net model will be found using 40 FISTA iterations. The spatial model will be updated at every 10 ADMM iterations.

```-u alpha``` : The regularization factor for the spatial constraint while solving the consensus problem. This factor is scaled according to the value of rho given to each cluster (see section 5 above). The cluster with the highest rho value will have the value of alpha given by ```-u```.


### 6) Solution format
All SAGECal solutions are stored as text files. Lines starting with '#' are comments.
The first non-comment line includes some general information, i.e.
freq(MHz) bandwidth(MHz) time_interval(min) stations clusters effective_clusters

The remaining lines contain solutions for each cluster as a single column, the first column is just a counter. 
Let's say there are K effective clusters and N directions. Then there will be K+1 columns, the first column will start from 0 and increase to 8N-1, 
which can be used to count the row number. It will keep repeating this, for each time interval.
The rows 0 to 7 belong to the solutions for the 1st station. The rows 8 to 15 for the 2nd station and so on. 
Each 8 rows of any given column represent the 8 values of a 2x2 Jones matrix. Lets say these are ```S0,S1,S2,S3,S4,S5,S6``` and ```S7```. Then the Jones matrix is ```[S0+j*S1, S4+j*S5; S2+j*S3, S6+j*S7]``` (the ';' denotes the 1st row of the 2x2 matrix).

When a cluster has a chunk size > 1, there will be more than 1 solution per given time interval. 
So for this cluster, there will be more than 1 column in the solution file, the exact number of columns being equal to the chunk size.



### Additional Info
See a [Tutorial](http://sagecal.sourceforge.net/tutorial/html/index.html)
and the [LOFAR Cookbook Chapter](https://support.astron.nl/LOFARImagingCookbook/sagecal.html).
# Code of Conduct

## 1. Purpose

A primary goal of SAGECal is to be inclusive to the largest number of contributors, with the most varied and diverse backgrounds possible. As such, we are committed to providing a friendly, safe and welcoming environment for all, regardless of gender, sexual orientation, ability, ethnicity, socioeconomic status, and religion (or lack thereof).

This code of conduct outlines our expectations for all those who participate in our community, as well as the consequences for unacceptable behavior.

We invite all those who participate in SAGECal to help us create safe and positive experiences for everyone.

## 2. Open Source Citizenship

A supplemental goal of this Code of Conduct is to increase open source citizenship by encouraging participants to recognize and strengthen the relationships between our actions and their effects on our community.

Communities mirror the societies in which they exist and positive action is essential to counteract the many forms of inequality and abuses of power that exist in society.

If you see someone who is making an extra effort to ensure our community is welcoming, friendly, and encourages all participants to contribute to the fullest extent, we want to know.

## 3. Expected Behavior

The following behaviors are expected and requested of all community members:

*   Participate in an authentic and active way. In doing so, you contribute to the health and longevity of this community.
*   Exercise consideration and respect in your speech and actions.
*   Attempt collaboration before conflict.
*   Refrain from demeaning, discriminatory, or harassing behavior and speech.
*   Be mindful of your surroundings and of your fellow participants. Alert community leaders if you notice a dangerous situation, someone in distress, or violations of this Code of Conduct, even if they seem inconsequential.
*   Remember that community event venues may be shared with members of the public; please be respectful to all patrons of these locations.

## 4. Unacceptable Behavior

The following behaviors are considered harassment and are unacceptable within our community:

*   Violence, threats of violence or violent language directed against another person.
*   Sexist, racist, homophobic, transphobic, ableist or otherwise discriminatory jokes and language.
*   Posting or displaying sexually explicit or violent material.
*   Posting or threatening to post other people’s personally identifying information ("doxing").
*   Personal insults, particularly those related to gender, sexual orientation, race, religion, or disability.
*   Inappropriate photography or recording.
*   Inappropriate physical contact. You should have someone’s consent before touching them.
*   Unwelcome sexual attention. This includes, sexualized comments or jokes; inappropriate touching, groping, and unwelcomed sexual advances.
*   Deliberate intimidation, stalking or following (online or in person).
*   Advocating for, or encouraging, any of the above behavior.
*   Sustained disruption of community events, including talks and presentations.

## 5. Consequences of Unacceptable Behavior

Unacceptable behavior from any community member, including sponsors and those with decision-making authority, will not be tolerated.

Anyone asked to stop unacceptable behavior is expected to comply immediately.

If a community member engages in unacceptable behavior, the community organizers may take any action they deem appropriate, up to and including a temporary ban or permanent expulsion from the community without warning (and without refund in the case of a paid event).

## 6. Reporting Guidelines

If you are subject to or witness unacceptable behavior, or have any other concerns, please notify a community organizer as soon as possible. yatawatta@astron.nl, f.diblen@esciencecenter.nl.



Additionally, community organizers are available to help community members engage with local law enforcement or to otherwise help those experiencing unacceptable behavior feel safe. In the context of in-person events, organizers will also provide escorts as desired by the person experiencing distress.

## 7. Addressing Grievances

If you feel you have been falsely or unfairly accused of violating this Code of Conduct, you should notify NLeSC Dirac team with a concise description of your grievance. Your grievance will be handled in accordance with our existing governing policies.



## 8. Scope

We expect all community participants (contributors, paid or otherwise; sponsors; and other guests) to abide by this Code of Conduct in all community venues–online and in-person–as well as in all one-on-one communications pertaining to community business.

This code of conduct and its related procedures also applies to unacceptable behavior occurring outside the scope of community activities when such behavior has the potential to adversely affect the safety and well-being of community members.

## 9. Contact info

yatawatta@astron.nl
f.diblen@esciencecenter.nl
h.spreeuw@esciencecenter.nl

## 10. License and attribution

This Code of Conduct is distributed under a [Creative Commons Attribution-ShareAlike license](http://creativecommons.org/licenses/by-sa/3.0/).

Portions of text derived from the [Django Code of Conduct](https://www.djangoproject.com/conduct/) and the [Geek Feminism Anti-Harassment Policy](http://geekfeminism.wikia.com/wiki/Conference_anti-harassment/Policy).

Retrieved on November 22, 2016 from [http://citizencodeofconduct.org/](http://citizencodeofconduct.org/)
# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a 
   build.
2. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters.
3. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
4. You may merge the Pull Request in once you have the sign-off of two other developers, or if you 
   do not have permission to do that, you may request the second reviewer to merge it for you.

## Code of Conduct

see [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

## Requirements for your pull request
 - [ ] Read the [contributing guide](https://github.com/nlesc-dirac/sagecal/blob/master/CONTRIBUTING.md) before creating a Pull request
 - [ ] Pull request should be motivated, i.e. what does it fix, why, and if relevant how
 - [ ] If possible / relevant include an example in the description, that could help all readers
       including project members to get a better picture of the change
 - [ ] Avoid other runtime dependencies
 - [ ] Fix all the merge conflicts
 - [ ] Meaningful commit history ; intention is important please rebase your commit history so that each
       commit is meaningful and help the people that will explore a change in later
 - [ ] The pull request follows coding style
 - [ ] Mention `Fixes #<issue number>` in the description _if relevant_
 - [ ] At least one commit should mention `Fixes #<issue number>` _if relevant_
 
Each pull request needs an approval of two developers and it should pass build test.
wo 10 nov 2021  8:52:47 CET
# SAGECal Installation

## Cmake Build
#### Ubuntu 20.04 (quick install)
```
 sudo apt-get install -y git cmake g++ pkg-config libcfitsio-bin libcfitsio-dev libopenblas-base libopenblas-dev wcslib-dev wcslib-tools libglib2.0-dev libcasa-casa4 casacore-dev casacore-data casacore-tools gfortran libopenmpi-dev libfftw3-dev

```
Run cmake (with GPU support) for example like
```
 mkdir build && cd build
 cmake .. -DHAVE_CUDA=ON -DCMAKE_CXX_FLAGS='-DMAX_GPU_ID=0' -DCMAKE_CXX_COMPILER=g++-8  -DCMAKE_C_FLAGS='-DMAX_GPU_ID=0' -DCMAKE_C_COMPILER=gcc-8 -DCUDA_NVCC_FLAGS='-gencode arch=compute_75,code=sm_75'
```
where *MAX_GPU_ID=0* is when there is only one GPU (ordinal 0). If you have more GPUs, increase this number to 1,2, and so on. This will produce *sagecal_gpu* and *sagecal-mpi_gpu* binary files (after running *make* of course).

CPU only version can be build as
```
 cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_C_COMPILER=gcc-8
```
which will produce *sagecal* and *sagecal-mpi*.

If you get **-lgfortran is not found** error, run the following in the build directory
```
 cd dist/lib
 ln -s /usr/lib/x86_64-linux-gnu/libgfortran.so.5 libgfortran.so
```
to make a symbolic link to libgfortran.so.5 or whatever version that is installed.

To only build *libdirac* library, use *-DLIB_ONLY=1* option. This library can be used with pkg-config using *lib/pkgconfig/libdirac.pc*.

### Requirements for older installations
#### das5

Load the modules below before compiling SageCal.
```
module load cmake/3.8.2
module load mpich/ge/gcc/64/3.2
module load gcc/4.9.3
module load casacore/2.3.0-gcc-4.9.3
module load wcslib/5.13-gcc-4.9.3
module load wcslib/5.16-gcc-4.9.3
module load cfitsio/3.410-gcc-4.9.3
```

checkout the source code and compile it with the instructions below(in source folder):
```
git clone https://github.com/nlesc-dirac/sagecal.git

cd sagecal && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
make
make install
```
$INSTALL_PATH is where you want to install SageCal.

#### Other systems

- Install equivalent packages for your distribution
    - g++
    - cmake
    - git
    - pkg-config
    - openblas
    - libglib2.0-dev
    - follow the instructions at
[https://github.com/casacore/casacore](https://github.com/casacore/casacore) to install casacore.
    - Additional packages (not essential, but recommended): MPI (openmpi), FFTW



### Building
- Clone the repository
```
    git clone -b master https://git@github.com/nlesc-dirac/sagecal.git

```

- Build SAGECal
```
    mkdir build && cd build
    cmake ..
```

**OPTIONAL:** You can also define a custom casacore path:

```
    cmake .. -DCASACORE_ROOT_DIR=/opt/soft/casacore
```
**OPTIONAL:** You can also define a custom paths to everything:

```
    cmake -DCFITSIO_ROOT_DIR=/cm/shared/package/cfitsio/3380-gcc-4.9.3 -DCASACORE_ROOT_DIR=/cm/shared/package/casacore/v2.3.0-gcc-4.9.3 -DWCSLIB_INCLUDE_DIR=/cm/shared/package/wcslib/5.13-gcc-4.9.3/include -DWCSLIB_LIBRARY=/cm/shared/package/wcslib/5.13-gcc-4.9.3/lib/libwcs.so -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_LINKER=/cm/shared/package/gcc/4.9.3/bin/gcc -DCMAKE_CXX_FLAGS=-L/cm/shared/package/cfitsio/3380-gcc-4.9.3/lib -DCMAKE_C_FLAGS=-L/cm/shared/package/cfitsio/3380-gcc-4.9.3/lib ..
```

    Compile with:
```
    make
```
    Install at your favorite place
```
    make DEST=/path/to/sagecal/dir install
```

- The sagecal executable can be found in **/path/to/sagecal/dir/usr/local/bin**, also **sagecal-mpi**,**buildsky** and **restore** might be installed depending on the availability of MPI and WCSLIB/FFTW.

### MPI support
MPI support is automatically detected, otherwise, it can be forced with:
```
cmake -DENABLE_MPI=ON
```

## GPU Support

### Loading modules on DAS5
See scripts folder for the modules.
```
source ./scripts/load_das5_modules_gcc6.sh
```

### Compiling with GPU support
```
mkdir -p build && cd build
cmake -DCUDA_DEBUG=ON -DDEBUG=ON -DHAVE_CUDA=ON ..
make VERBOSE=1
```



## Installation via Anaconda (WIP)
```
    conda install -c sagecal=0.6.0
```



## Manual installation
For expert users, and for custom architectures (GPU), the manual install is recommended.
### 1 Prerequisites:
 - CASACORE http://casacore.googlecode.com/
 - glib http://developer.gnome.org/glib
 - BLAS/LAPACK
   Highly recommended is OpenBLAS http://www.openblas.net/
   Also, to avoid any linking issues (and to get best performance), build OpenBLAS from source and link SAGECal with the static library (libopenblas***.a) and NOT libopenblas***.so
 - Compilers gcc/g++ or Intel icc/icpc
 - If you have NVIDIA GPUs,
  -- CUDA/CUBLAS/CUSOLVER and nvcc
  -- NVML Nvidia management library
 - If you are using Intel Xeon Phi MICs.
  -- Intel MKL and other libraries
 - Get the source for SAGECal
```
    git clone -b master https://git@github.com/nlesc-dirac/sagecal.git
```

### 2 The basic way to build is
  1.a) go to ./src/lib/Dirac and ./src/lib/Radio  and run make (which will create libdirac.a and libradio.a)
  1.b) go to ./src/MS and run make (which will create the executable)


### 3 Build settings
In ./src/lib/Dirac and ./src/lib/Radio and ./src/MS you MUST edit the Makefiles to suit your system. Some common items to edit are:
 - LAPACK: directory where LAPACK/OpenBLAS is installed
 - GLIBI/GLIBL: include/lib files for glib
 - CASA_LIBDIR/CASA_INCDIR/CASA_LIBS : casacore include/library location and files:
  Note with new CASACORE might need two include paths, e.g.
    -I/opt/casacore/include/ -I/opt/casacore/include/casacore
 - CUDAINC/CUDALIB : where CUDA/CUBLAS/CUSOLVER is installed
 - NVML_INC/NVML_LIB : NVML include/lib path
 - NVCFLAGS : flags to pass to nvcc, especially -arch option to match your GPU
 - MKLROOT : for Intel MKL

 Example makefiles:
   Makefile : plain build
   Makefile.gpu: with GPU support
   Note: Edit ./lib/Radio/Radio.h MAX_GPU_ID to match the number of available GPUs, e.g., for 2 GPUs, MAX_GPU_ID=1



## SAGECAL-MPI Manual Installation
This is for manually installing the distributed version of sagecal (sagecal-mpi), the cmake build will will work for most cases.
## 1 Prerequsites:
 - Same as for SAGECal.
 - MPI (e.g. OpenMPI)

## 2 Build ./src/lib/Dirac ./src/lib/Radio as above (using mpicc -DMPI_BUILD)

## 3 Build ./src/MPI using mpicc++



## BUILDSKY Installation

  - See INSTALL in ./src/buildsky


## RESTORE Installation

  - See INSTALL in ./src/restore



---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Using libdirac in general optimization
This directory includes some simple demo programs illustrating how to use optimization routines in libdirac for solving any general optimization problem. We higlight the use of (stochastic) limited-meory Broyden Fletcher Goldfarb Shanno (LBFGS) algorithm.

  * `demo.c`: full batch LBFGS
  * `demo_stochastic.c`: minibatch (stochastic) LBFGS
  * `demo_stochastic_cuda.c`: minibatch LBFGS with full GPU acceleration

Use `Makefile` to build `demo` and `demo_stochastic`. Use `Makefile.cuda` to build `demo_stochastic_cuda`.

Note that for the CUDA demo, libdirac should be built with CUDA support, by using `-DHAVE_CUDA=ON` cmake option.# Calibration examples
This directory includes a few examples of calibration using SAGECal. Direction independent calibration, direction dependent calibration as well as bandpass calibration is possible. Two shell scripts are included:

  * `dosage.sh`: calibration in a single computer
  * `dosage-mpi.sh`: calibration using a cluster of computers, using MPI

A small test data file `sm.ms.tar` is included, and you need to untar it before running the examples. The input model is given by:
  * `3c196.sky.txt`: sky model 
  * `3c196.sky.txt.cluster`: cluster file that specifies the directions being calibrated

## Distributed calibration
Here is a step-by-step guide to get going with distributed calibration. You only need one computer to test this, but using MPI configuration options, the same steps can be carried out over a cluster.

1 Untar `sm.ms.tar` like
```console
tar xvf sm.ms.tar
```
 then you will have `sm.ms` dataset (it is actually a directory).

2 Copy this to multiple files such as
```console
cp -r sm.ms sm1.ms
cp -r sm.ms sm2.ms
...
```

3 After this, use [this script](https://github.com/nlesc-dirac/sagecal/blob/master/test/Calibration/Change_freq.py) to change the frequency of each data file, like
```console
python2 ./Change_freq.py sm1.ms 110e6
python2 ./Change_freq.py sm2.ms 120e6
python2 ./Change_freq.py sm3.ms 130e6
...
```
 Each data file should have a distinct frequency, and the above method changes the frequency of each dataset to a unique value.

4 After this step, you can run `dosage-mpi.sh` to calibrate all datasets matching `*.ms` in the current directory.
