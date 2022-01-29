# HH-suite3 for sensitive sequence searching

(C) Johannes Soeding, Markus Meier, Martin Steinegger, Milot Mirdita, Michael Remmert, Andreas Hauser, Andreas Biegert

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/hhsuite.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/hhsuite)
[![Biocontainer Pulls](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dhhsuite)](https://biocontainers.pro/#/tools/hhsuite)
[![Github All Releases](https://img.shields.io/github/downloads/soedinglab/hh-suite/total.svg)](https://github.com/soedinglab/hh-suite/releases/latest)
[![Docker Pulls](https://img.shields.io/docker/pulls/soedinglab/hh-suite.svg)](https://hub.docker.com/r/soedinglab/hh-suite)
[![Build Status](https://dev.azure.com/themartinsteinegger/hhsuite/_apis/build/status/soedinglab.hh-suite?branchName=master)](https://dev.azure.com/themartinsteinegger/hhsuite/_build/latest?definitionId=4&branchName=master)

The HH-suite is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

## Documentation

We provide an extensive [user guide](https://github.com/soedinglab/hh-suite/wiki) with many usage examples, frequently asked questions and guides to build your own databases. 

### Installation

HH-suite3 can also be installed by downloading a statically compiled version, [conda](https://github.com/conda/conda) or [Docker](https://github.com/moby/moby). HH-suite3 requires a 64-bit system (check with `uname -a | grep x86_64`). On AMD/Intel CPUs it requires at least support for the SSE2 instruction set (check by executing `cat /proc/cpuinfo | grep sse2` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE2` on macOS). `AVX2` is roughly 2x faster compared to SSE2. HH-suite3 also works on Linux systems with ARM64 and PPC64LE CPUs. Precompiled binaries for all supported systems can be found at [mmseqs.com/hhsuite](https://mmseqs.com/hhsuite).

```
# install via conda
conda install -c conda-forge -c bioconda hhsuite 
# install docker
docker pull soedinglab/hh-suite
# static SSE2 build
wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-SSE2-Linux.tar.gz; tar xvfz hhsuite-3.3.0-SSE2-Linux.tar.gz; export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
# static AVX2 build
wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-AVX2-Linux.tar.gz; tar xvfz hhsuite-3.3.0-AVX2-Linux.tar.gz; export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```
:exclamation: Only the self-compiled HH-suite3 version includes MPI support, since MPI configuration is specific to the local environment.

### Available Databases
List of available database for HH-suite3: 
  - [Uniclust30](https://uniclust.mmseqs.com) [[pub]](https://doi.org/10.1093/nar/gkw1081)
  - [BFD](https://bfd.mmseqs.com) (consists of 2.5 billion, mostly enviromental, protein sequences) [[pub]](https://doi.org/10.1038/s41592-019-0437-4)
  - [Pfam/SCOP/PDB70/dbCAN](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)


### Compilation
To compile from source, you will need a recent C/C++ compiler (at least GCC 4.8 or Clang 3.6) and [CMake](http://cmake.org/) 2.8.12 or later.

To download the source code and compile the HH-suite execute the following commands:
```
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

:exclamation: To compile HH-suite3 on macOS, first install the `gcc` compiler from [Homebrew](https://brew.sh). The default macOS `clang` compiler does not support OpenMP and HH-suite3 will only be able to use a single thread. Then replace the `cmake` call above with the following one:

```
CC="$(brew --prefix)/bin/gcc-10" CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_INSTALL_PREFIX=. ..
```    


## Usage
For performing a single search iteration of HHblits, run HHblits with the following command:
```
hhblits -i <input-file> -o <result-file> -n 1 -d <database-basename>
```

For generating an alignment of homologous sequences:
```
hhblits -i <input-file> -o <result-file> -oa3m <result-alignment> -d <database-basename>
```

A detailed list of options for HHblits is available by running HHblits with the `-h` parameter.

## Reference

Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J (2019)
HH-suite3 for fast remote homology detection and deep protein annotation, *BMC Bioinformatics*, 473. [doi: 10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7)

## Links

* [Soeding lab](http://www.mpibpc.mpg.de/soeding)
* [Databases for the HH-suite](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/)
# FFindex - A database wrapped around mmap

FFindex is a very simple index/database for huge amounts of small files. The
files are stored concatenated in one big data file, seperated by '\0'. A second
file contains a plain text index, giving name, offset and length of of the
small files. The lookup is currently done with a binary search on an array made
from the index file. The attatched binaries (see Usage below) and their source
code shall give an impression of how to use the functions supported by the library in C/C++ code.
 
## Copyright

FFindex was written by Andreas Hauser <hauser@genzentrum.lmu.de>.

Please add your name here if you distribute modified versions.
* Martin Steinegger <martin.steinegger@mpibpc.mpg.de>
* Markus Meier <markus.meier@mpibpc.mpg.de>

FFindex is provided under the Create Commons license "Attribution-ShareAlike 4.0",
which basically captures the spirit of the Gnu Public License (GPL).

See:
http://creativecommons.org/licenses/by-sa/4.0/

## Thanks

Thanks to Laszlo Kajan for creating and maintaining Debian packages
and many suggestions to improve the build and user experience.



## Installation

### Compilation
With the sourcecode ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE_DIR} ..
	make
	make install


**Please use a sensible value for ${INSTALL_BASE_DIR}, e.g. /usr/local or /opt/ffindex or $HOME/ffindex**


### Setting environment variables

Please note that before querying or unlinking entries a ffindex must be
sorted, although you can add to it without. So either specify -s with
ffindex_build or sorted later with ffindex_modify -s.

Setup environment:

	export PATH="${INSTALL_BASE_DIR}/bin:${PATH}"
	export LD_LIBRARY_PATH="${INSTALL_BASE_DIR}/lib:${LD_LIBRARY_PATH}"
On OS X set DYLD_LIBRARY_PATH instead of LD_LIBRARY_PATH.

## Usage

Build index from files in test/data and test/data2.

	ffindex_build -s /tmp/test.data /tmp/test.ffindex test/data test/data2

Retrieve three entries:

	ffindex_get  /tmp/test.data /tmp/test.ffindex a b foo

Unlink (Remove reference from index) an entry:

	ffindex_modify -u /tmp/test.ffindex b

Retrieve three entries, "b" should now be missing:

	ffindex_get /tmp/test.data /tmp/test.ffindex a b foo

Convert a Fasta file to ffindex, entry names are incerental IDs starting from 1:

	ffindex_from_fasta -s fasta.ffdata fasta.ffindex NC_007779.ffn

Get first entry by name:

	ffindex_get fasta.ffdata fasta.ffindex 1

Get first and third entry by entry index, this a little faster:

	ffindex_get fasta.ffdata fasta.ffindex -n 1 3

Count the characters including header in each entry:

	ffindex_apply fasta.ffdata fasta.ffindex wc -c

Count the number of characters in each sequence, without the header:

	ffindex_apply fasta.ffdata fasta.ffindex perl -ne '$x += length unless(/^>/); END{print "$x\n"}'

Parallel version for counting the characters including header in each entry:

	mpirun -np 4 ffindex_apply_mpi fasta.ffdata fasta.ffindex -- wc -c

Parallel version for counting the characters including header in each entry and
saving the output to a new ffindex:

	mpirun -np 4 ffindex_apply_mpi fasta.ffdata fasta.ffindex -i out-wc.ffindex -o out-wc.ffdata -- wc -c
# SIMDe Without Test Cases

This repository contains only the core of
[SIMDe](https://github.com/simd-everywhere/simde).
It is generated automatically for every commit to master, and is
intended to be used as a submodule in projects which don't want to
include the (rather large) test cases.

All development work happens in the main repository, please do not
file issues or create pull requests against this repository.
:exclamation: Make to check out our [User Guide](https://github.com/soedinglab/hh-suite/wiki).

## Expected Behavior

## Current Behavior

## Steps to Reproduce (for bugs)
Please make sure to execute the reproduction steps.

## HH-suite Output (for bugs)
Please make sure to post the complete output of the tool you called. Please use [gist.github.com](https://gist.github.com).

## Context
Providing context helps us come up with a solution and improve our documentation for the future.

## Your Environment
Include as many relevant details about the environment you experienced the issue in.
* Version/Git commit used:
* Server specifications (especially CPU support for AVX2/SSE and amount of system memory):
* Operating system and version:
