# Release 2.0.5  2022-01-25

This release fixes fixes a major (though rare) bug where the reference phase could take forever on some inputs.  It includes a newer version of lastz which seems to fix some crashes as well.

- Debug symbols no longer stripped from `cactus_consolidated` binary in Release.
- Update to Toil 3.5.6
- Update examples to use Python 3.8 specifically (was previously just python3, but this is often python3.6, support for which was dropped in Toil 3.5.6)
- cactus-prepare WDL output can now batch up `hal_append_subtree` jobs
- Fix bug where "reference" phase within cactus_consolidated could take ages on some input
- Fix bug where --realTimeLogging flag would cause infinite loop after cactus_consolidated within some Docker invocations.
- Upgrade to more recent version of lastz

# Release 2.0.4   2021-11-12

This release fixes a bug introduced in 2.0.0 where ancestral sequences could not be specified in the input, which prevented the recommended producedure for updating existing alignments from working. 

- Fix `Assertion `cap_getSequence(cap) == sequence' failed` error when ancestral fasta provided in input seqfile.
- Several minor pangenome updates, mostly in `cactus-graphmap-join`

# Release 2.0.3   2021-07-22

This release fixes some issues in pangenome normalization and CAF running time.

- Fix new regression that caused CAF's secondary filter to sometimes take forever.  This code has been causing occaisional slowdowns for some time, but should finally be fixed once and for all.
- Fix cactus-preprocess to work on zipped fasta inputs even when not running dna-brnn.
- Fix normalization in cactus-graphmap-join
- Update to abPOA v1.2.5


# Release 2.0.2   2021-07-07

This release primarily addresses stability issues during pangenome construction.

Changelog:
- Use latest abpoa, which fixes bug where aligning >1024 sequences would lead to a segfault
- Update to Toil 5.4.0
- More consistently apply filters to minimap2 output in the fallback stage of graphmap-split
- Build abpoa with AVX2 SIMD extensions instead of SSE4.1 in order to work around instability when building pangenomes.  This ups the hardware requirements for releases, unfortunately, as AVX2 is slightly newer.
- Clean up CAF config parameters
- Fix CAF secondary filter worst-case runtime issue.  It was very rare but could add days to runtime.
- Slightly tune minimap2 thresholds used for chromosome splitting
- Normalization option added to cactus-graphmap-join (should be used to work around soon-to-be addressed underalignment bug)

# Release 2.0.1   2021-06-19

This a patch release that fixes an issue where the new `--consCores` option could not be used with `--maxCores` (Thanks @RenzoTale88).  It also reverts some last minute CAF parameter changes to something more tested (though known to be slow in some cases with large numbers of secondary alignments)

Changelog:
- Fix bug where `cactus` doesn't work when both `--maxCores` and `--consCores` are specified.
- Static binaries script more portable.
- Revert CAF trimming parameters to their previous defaults. 

# Release 2.0.0   2021-06-18

This release includes a major update to the Cactus workflow which should dramatically improve both speed and robustness. Previously, Cactus used a multiprocess architecture for all cactus graph operations (everything after the "blast" phase).  Each process was run in its own Toil job, and they would communicate via the CactusDisk database that ran as its own separate service process (ktserver by default). Writing to and from the database was often a bottleneck, and it would fail sporadically on larger inputs with frustrating "network errors". This has all now been changed to run as a single multithreaded executable, `cactus_consolidated`.  Apart from saving on database I/O, `cactus_consolidated` now uses the much-faster, SIMD-accelerated abPOA by default instead of cPecan for performing multiple sequence alignments within the BAR phase.

Cactus was originally designed for a heterogeneous compute environment where a handful of large memory machines ran a small number of jobs, and much of the compute could be farmed off to a large number of smaller machines.  While lastz jobs from the "preprocess" and "blast" phases (or `cactus-preprocess` and `cactus-blast`) can still be farmed out to smaller machines, the rest of cactus (`cactus-align`) can now only be run on more powerful systems.  The exact requirements depend as usual on genome size and divergence, but roughly 64 cores / 512G RAM are required for distant mammals. 

This release also contains several fixes and usability improvements for the pangenome pipeline, and finally includes `halPhyloP`.

Changlelog:
- Fold all post-blast processing into single binary executable,`cactus_consolidated`
- New option, `--consCores`, to control the number of threads for each `cactus_consolidated` process.
- Cactus database (ktserver) no longer used.
- abPOA now default base aligner, replacing cPecan
- cPecan updated to include multithreading support via MUM anchors (as opposed to spawning lastz processes), and can be toggled on in the config
- Fix bug in how `cactus-prepare` transmits Toil size parameters
- `cactus-prepare-join` tool added to combine and index chromosome output from `cactus-align-batch`
- `cactus-graphmap-split` fixes
- Update to latest Segalign
- Update to Toil 5.3
- Update HAL
- Add `halPhyloP` to binary release and docker images

# Release 1.3.0   2021-02-11

This release introduces the [Cactus Pangenome Pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md), which can be used to align together samples from the same species in order to create a pangenome graph:

- `cactus_bar` now has a POA-mode via the abpoa aligner, which scales better than Pecan for large numbers of sequences and is nearly as accurate if the sequences are highly similar
- `cactus-refmap` tool added to produce cactus alignment anchors with all-to-reference minimap2 alignments instead of all-to-all lastz
- `cactus-graphmap` tool added to produce cactus alignment anchors with all-to-reference-graph minigraph alignments instead of all-to-all lastsz
- `--maskAlpha` option added to `cactus-preprocess` to softmask (or clip out) satellite sequence using `dna-brnn.
- `cactus_bar` now has an option to ignore masked sequence with a given length threshold.
- `cactus-graphmap-split` tool added to split input fasta sequences into chromosomes using minigraph alignments in order to create alignment subproblems and improve scaling.
- `cactus-align-batch` tool added to align several chromsomes at once using one job per chromosome. (`--batch` option added to `cactus-align` to achieve the same using many jobs per chromosome)
- `--outVG` and `outGFA` options added to `cactus-align` to output pangenome graphs in addtion to hal.

Other changes:
- `cactus-prepare` scheduling bug fix
- `--database redis` option added to use Redis instead of Kyoto Tycoon
- `cactus-blast --restart` bug fix
- "Legacy" binary release provided for those whose hardware is too old to run the normal release. 

# Release 1.2.3   2020-10-05

- Fix bug where `cactus_fasta_softmask_intervals.py` was expecting 1-based intervals from GPU lastz repeatmasker.

hal2vg version included: [v1.0.1](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.1/hal2vg)
GPU Lastz version used in GPU-enabled Docker image: [8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7)

# Release 1.2.2   2020-10-02

- hal2vg updated to version 1.0.1
- GPU lastz updated for more disk-efficient repeat masking and better error handling
- Fixed memory in `cactus_convertAlignmentsToInternalNames`
- CAF fixes targeted towards pangenome construction

hal2vg version included: [v1.0.1](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.1/hal2vg)
GPU Lastz version used in GPU-enabled Docker image: [8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7)

# Release 1.2.1   2020-08-31

This release fixes bugs related to GPU lastz

- Cactus fixed to correctly handle GPU lastz repeatmasking output, as well to keep temporary files inside Toil's working directory.
- `cactus-prepare --wdl` updated to support `--preprocessorBatchSize > 1`
- `cactus_covered_intervals` bug fix and speedup
- GPU lastz updated to fix crash

GPU Lastz version used in GPU-enabled Docker image: [12af3c295da7e1ca87e01186ddf5b0088cb29685](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/12af3c295da7e1ca87e01186ddf5b0088cb29685)
hal2vg version included: [v1.0.0](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.0/hal2vg)


# Release 1.2.0   2020-08-21

This release adds GPU lastz repeatmasking support in the preprocessor, and includes hal2vg

Notable Changes:
 - GPU lastz repeat masking is about an order of magnitude faster than CPU masking, and should provide better results.  It's toggled on in the config file or by using the GPU-enabled Docker image.
 - hal2vg (pangenome graph export) included in the binary release as well as docker images.
 - update hal to [f8f3fa2dada4751b642f0089b2bf30769967e68a](https://github.com/ComparativeGenomicsToolkit/hal/commit/f8f3fa2dada4751b642f0089b2bf30769967e68a)

GPU Lastz version used in GPU-enabled Docker image: [f84a94663bbd6c42543f63b50c5843b0b5025dda](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/f84a94663bbd6c42543f63b50c5843b0b5025dda)
hal2vg version included: [v1.0.0](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.0/hal2vg)

# Release 1.1.1   2020-07-31

This release fixes how Kent tools required for `hal2assemblyHub.py` were packaged in 1.1.0 (thanks @nathanweeks).  

Notable Changes:
 - The required shared libaries to run the Kent tools are added to the Docker Image
 - The same Kent tools are removed from the binary release.  They were included under the assumption that they were statically built and fully standalone, but they are not.  Instead, instrucitons are provided to guide interested users to installing them on their own. 

GPU Lastz version used in GPU-enabled Docker image: [3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa)

# Release 1.1.0   2020-07-30

This release contains some important bug fixes, as well as major improvements to `cactus-prepare` functionality.

Notable Changes:
 - `cactus-prepare` improvements including:
    - WDL / Terra support
    - GPU lastz support
    - bug fixes
- Upgrade to Toil 4.1.0
- Fix bug causing `cactus-reference` to run forever in presence of 0-length branches
- Major speed improvement for `cactus-caf` by fixing secondary alignment filter.
- Include HAL python tools in Docker images and binary release
- Fix static binary for `cPecanRealign`
- Provide GPU-enabled Docker image for release
- Included HAL tools contains several crash fixes

GPU Lastz version used in GPU-enabled Docker image: [3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa)

# Release 1.0.0   2020-04-19

This is the first official release of Cactus.  The goal is to provide a
stable, track-able version of Cactus to users.  The releases is provided in
source, static-compiled binaries, and Docker formats.

Notable Changes:
 - Kyoto Cabinet and Typhoon are now included as a sub-module.  This is due to
   the lack of consistent, stable releases and the difficulty in compiling it.
 - Cactus is now available as a tar file of static-linked binary executables,
   along with a wheel of the Cactus Python packages.  This avoids compilation and dependency problems. These should work on most Linux platforms when Docker is not used.
 - Added support to run Cactus in individual steps. This works around problems with using Toil in some distributed environments by dividing up alignment into tasks that can be run manually on separate machines.
   See *Running step by step (experimental)* in `README.md`.
 - Conversion to Python 3, allowing Toil to drop Python 2 support.


# Notes on developing and debugging cactus

## Overriding make settings
A file include.local.mk can be created in the root directory
to override make variables, including setting environment variables.
This should not be committed.

## Environment variables controlling how cactus is run
- CACTUS_BINARIES_MODE - how are cactus programs found?
  - docker <default>
  - singularity
  - local
- CACTUS_DOCKER_MODE - is Docker being used?
  - 1 <default>
  - 0
- CACTUS_USE_LOCAL_IMAGE - is Docker image on local server?
  - 0 <default>
  - 1

## Environment variables controlling tests
- SON_TRACE_DATASETS location of test data set, currently available with
    git clone https://github.com/ComparativeGenomicsToolkit/cactusTestData

- SONLIB_TEST_LENGTH  filters tests by maximum run time length category (case-insensitive)
  - SHORT - tests taking less than ~10 seconds, with some exceptions <default>
  - MEDIUM - tests taking less than ~100 seconds
  - LONG - test taking less than ~1000 seconds
  - VERG_LONG - test taking even longer

- CACTUS_TEST_LOG_LEVEL - Set log-level used for the test, may not set it for all test, but very useful for Toil.
  


## Running tests with docker in single machine mode
    make docker
    export CACTUS_USE_LOCAL_IMAGE=1
    make test

## Debugging hints
   - The main Cactus Python process will print out a stack trace of all of the Python
     threads if sent a SIGUSR1 signal.  They will then continue execution.  This
     maybe useful in determining the state of cactus.
   
   

# Installation of the Cactus binary distribution 

This describes the steps require to install the Cactus
pre-compile binary, static linked distribution.

## Extracting
If you have not already extract the distribution:
```
tar -xzf "cactus-bin-${REL_TAG}.tar.gz"
cd "cactus-bin-${REL_TAG}"
```

## Setup

To build a python virtualenv and activate, do the following steps:
```
virtualenv -p python3.6 venv
source venv/bin/activate
pip install -U setuptools pip
pip install -U -r ./toil-requirement.txt
pip install -U .
export PATH=$(pwd)/bin:$PATH
export PYTHONPATH=$(pwd)/lib:$PYTHONPATH
```

Some tools required for `hal2assemblyHub.py` are not included and must be downloaded separately.
They are `wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent`.  More information
can be found [here](https://hgdownload.cse.ucsc.edu/admin/exe/).  Note that some may require
a license for commercial use.  Static binaries are not available, but the following command
should set them up successfully on many 64 bit Linux systems:
```
cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done
```

## Testing

To test Cactus, the following will run a moderately sized alignment.  It may
take several hours, depending on your system.
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
``
# Cactus
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/cactus.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/cactus)

Cactus is a reference-free whole-genome multiple alignment program. Please cite the [Progressive Cactus paper](https://doi.org/10.1038/s41586-020-2871-y) when using Cactus.  Additional descriptions of the core algorithms can be found [here](https://doi.org/10.1101/gr.123356.111) and [here](https://doi.org/10.1089/cmb.2010.0252).

Please subscribe to the [cactus-announce](https://groups.google.com/d/forum/cactus-announce) low-volume mailing list so we may reach about releases and other announcements.

## Acknowledgements

Cactus uses many different algorithms and individual code contributions, principally from Joel Armstrong, Glenn Hickey, Mark Diekhans and Benedict Paten. We are particularly grateful to:

- Yung H. Tsin and Nima Norouzi for contributing their 3-edge connected components program code, which is crucial in constructing the cactus graph structure, see: Tsin,Y.H., "A simple 3-edge-connected component algorithm," Theory of Computing Systems, vol.40, No.2, 2007, pp.125-142.
- Bob Harris for providing endless support for his [LastZ](https://github.com/lastz/lastz) pairwise, blast-like genome alignment tool.
- Sneha Goenka and Yatish Turakhia for [SegAlign](https://github.com/gsneha26/SegAlign), the GPU-accelerated version of LastZ.
- Yan Gao et al. for [abPOA](https://github.com/yangao07/abPOA)
- Heng Li for [minigraph](https://github.com/lh3/minigraph), [minimap2](https://github.com/lh3/minimap2), [gfatools](https://github.com/lh3/gfatools) and [dna-brnn](https://github.com/lh3/dna-rnn)
- Dany Doerr for [GFAffix](https://github.com/marschall-lab/GFAffix), used to optionally clean pangenome graphs.
- The vg team for [vg](https://github.com/vgteam/vg), used to process pangenome graphs.

## Setup

### System requirements
We regularly test on Ubuntu 18.04 (Bionic) and to a more limited degree on Mac OS X (using Docker).

Cactus requires Python >= 3.7.

Cactus uses substantial resources. For primate-sized genomes (3 gigabases each), you should expect Cactus to use approximately 120 CPU-days of compute per genome, with about 120 GB of RAM used at peak. The requirements scale roughly quadratically, so aligning two 1-megabase bacterial genomes takes only 1.5 CPU-hours and 14 GB RAM.

Note that to run even the very small evolverMammals example, you may need up to 12 GB RAM. The actual resource requirements are much less, but the individual jobs have resource estimates based on much larger alignments, so the jobs will refuse to run unless there are enough resources to meet their estimates.

IMPORTANT:  It is highly recommend that one **not** run Cactus using the Toil Grid Engine-like batch systems (GridEngine, HTCondor, LSF, SLURM, or Torque).  Cactus creates a very large number of small jobs, which can overwhelm these systems.  There is a work-around described [here](#running-the-step-by-step-workflow-direclty-in-toil) for clusters with large compute nodes available.  **Update**:  Cactus version >= 2.0 with GPU enabled will spawn far fewer jobs which, in theory, should make this less of an issue.

NEW:  Cactus can now align individuals from the same species without a tree using the [Cactus Pangenome Pipeline](doc/pangenome.md).

### Installation Overview

There are many different ways to install and run Cactus:
* [Docker Image](#docker-image)
* [Precompiled Binaries](#precompiled-binaries)
* [Build From Source](#build-from-source)
* [Python Install with Docker Binaries](#python-install-with-docker-binaries)

#### Docker Image

Cactus docker images are hosted on [quay](https://quay.io/repository/comparative-genomics-toolkit/cactus).  The image for the latest release is listed on the [Releases Page](https://github.com/ComparativeGenomicsToolkit/cactus/releases).  Here is an command line to run the included evolver mammals example with release 2.0.5
```
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt
docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.0.5 cactus /data/jobStore /data/evolverMammals.txt /data/evolverMammals.hal --root mr --binariesMode local

```

#### Precompiled Binaries

Precompiled binaries can be found on the [Releases Page](https://github.com/ComparativeGenomicsToolkit/cactus/releases).  Download by clicking the `cactus-bin-vx.x.x.tar.gz` and install following the instructions in the `BIN-INSTALL.md` link.

#### Build From Source

##### Install the Cactus Python package and its dependencies

To avoid problems with conflicting versions of dependencies on your system, we strongly recommend installing Cactus inside a Python 3 [virtual environment](https://virtualenv.pypa.io/en/stable/). To install the `virtualenv` command, if you don't have it already, run:
```
python3 -m pip install virtualenv
```

To set up a virtual environment in the directory `cactus_env`, run:
```
python3 -m virtualenv -p python3.8 cactus_env
```

Then, to enter the virtualenv, run:
```
source cactus_env/bin/activate
```

You can always exit out of the virtualenv by running `deactivate`.


To install Cactus in Python, clone it and **its submodules with --recursive** from github and install it with pip:
```
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
python3.8 -m pip install --upgrade setuptools
python3.8 -m pip install --upgrade -r toil-requirement.txt
python3.8 -m pip install --upgrade .
```

##### Build the Cactus Binaries

Several binaries are required to run Cactus.  They can be built as follows.

Compile time settings can be overridden by creating a make include file in the top level cactus directory.  
```
cactus/include.local.mk
```

Cactus has several dependencies that need to be installed on the system, including HDF5. HDF5 is available through most package managers (`apt-get install libhdf5-dev`) or can be manual installed from source files at [The HDF Group](https://www.hdfgroup.org/).   HDF5 should be configured with the `--enable-cxx` option. If you've installed it in a non-standard location, have the `h5c++` command in your `PATH` or add this to `include.local.mk`:
`export PATH := <hdf5 bin dir>:${PATH}`

You can use the the [Dockerfile](Dockerfile) as a guide to see how all dependencies are installed with `apt` on Ubuntu.

In the top level cactus directory.  The binaries can then be built with
```
make -j $(nproc)
```
and added to the PATH with
```
export PATH=$(pwd)/bin:$PATH
```

In order to run the [Cactus Pangenome Pipeline](doc/pangenome.md), additional tools must be installed with:
```
build-tools/downloadPangenomeTools
```

To use HAL python scripts such as `hal2mafMP.py`, add the submodules directory to the PYTHONPATH with
```
export PYTHONPATH=$(pwd)/submodules:$PYTHONPATH
```

#### Python Install With Docker Binaries

Cactus can be setup and used in a virtual environment as in the [previous section](#build-from-source), without compiling the binaries.  When used like this (which will happen automatically when running `cactus` without the appropriate binaries in the `PATH` environment variable), a Docker image will be automatically pulled to run commands as needed.  The main use case for this is running with Toils AWS provisioner as [described here](doc/running-in-aws.md).

Singularity binaries can be used in place of docker binaries with the `--binariesMode singularity` flag.  Note, you must use Singularity 2.3 - 2.6 or Singularity 3.1.0+. Singularity 3 versions below 3.1.0 are incompatible with cactus (see [issue #55](https://github.com/ComparativeGenomicsToolkit/cactus/issues/55) and [issue #60](https://github.com/ComparativeGenomicsToolkit/cactus/issues/60)).

The `--binariesMode local` flag can be used to force `cactus` to run local binaries -- this is the default behavior if they are found.

## Running
To run Cactus, the basic format is:
```
cactus <jobStorePath> <seqFile> <outputHal>
```

Note: alternative ways of running include the [step-by-step interface](#running-step-by-step) and the [Cactus Pangenome Pipeline](doc/pangenome.md).

The `jobStorePath` is where intermediate files, as well as job metadata, will be stored. It must be accessible to all worker systems.

When first testing out Cactus on a new system or cluster, before running anything too large, try running the small (5 600kb genomes) simulated example in `examples/evolverMammals.txt`. It should take less than an hour to run on a modern 4-core system. That example, even though it's small, should be enough to expose any major problems Cactus may have with your setup. The command you should run is:
```
cactus jobStore examples/evolverMammals.txt examples/evolverMammals.hal --root mr
```

Within an hour at most (on modern computers), you should have a [HAL](https://github.com/ComparativeGenomicsToolkit/hal) file which relates simulated mouse and rat genomes.


### seqFile: the input file
The input file, called a "seqFile", is just a text file containing the locations of the input sequences as well as their phylogenetic tree. The tree will be used to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parents in a bottom-up fashion. Polytomies in the tree are allowed, though the amount of computation required for a sub-alignment rises quadratically with the degree of the polytomy.  Cactus uses the predicted branch lengths from the tree to determine appropriate pairwise alignment parameters, allowing closely related species to be aligned more quickly with no loss in accuracy. The file is formatted as follows:

    NEWICK tree
    name1 path1
    name2 path2
    ...
    nameN pathN

An optional * can be placed at the beginning of a name to specify that its assembly is of reference quality. This implies that it can be used as an outgroup for sub-alignments. If no genomes are marked in this way, all genomes are assumed to be of reference quality. The star should only be placed on the name-path lines and not inside the tree.

* The tree must be on a single line. All leaves must be labeled and these labels must be unique. Ancestors may be named, or left blank (in which case the ancestors in the final output will automatically be labeled Anc0, Anc1, etc.) Labels must not contain any spaces.
* Branch lengths that are not specified are assumed to be 1.
* Lines beginning with # are ignored. 
* Sequence paths must point to either a FASTA file or a directory containing 1 or more FASTA files.
* Sequence paths must not contain spaces.
* Each name / path pair must be on its own line
* `http://`, `s3://`, etc. URLs may be used.


Please ensure your genomes are *soft*-masked with RepeatMasker. We do some basic masking as a preprocessing step to ensure highly repetitive elements are masked when repeat libraries are incomplete, but genomes that aren't properly masked can still take tens of times longer to align that those that are masked. Hard-masking (totally replacing repeats with stretches of Ns) isn't necessary, and is strongly discouraged (you will miss a *lot* of alignments!).

Example:

	  # Sequence data for progressive alignment of 4 genomes
	  # human, chimp and gorilla are flagged as good assemblies.
	  # since orang isn't, it will not be used as an outgroup species.
     (((human:0.006,chimp:0.006667):0.0022,gorilla:0.008825):0.0096,orang:0.01831);
     *human /data/genomes/human/human.fa
     *chimp /data/genomes/chimp/
     *gorilla /data/genomes/gorilla/gorilla.fa
     orang /cluster/home/data/orang/

### Running on a cluster
Cactus (through Toil) supports many batch systems in theory, including LSF, SLURM, GridEngine, Parasol, and Torque. To run on a cluster, add `--batchSystem <batchSystem>`, e.g. `--batchSystem gridEngine`. If your batch system needs additional configuration, Toil exposes some [environment variables](http://toil.readthedocs.io/en/3.10.1/developingWorkflows/batchSystem.html#batch-system-enivronmental-variables) that can help.

IMPORTANT:  It is highly recommend that one **not** run Cactus in its default mode using the Toil Grid Engine-like batch systems (GridEngine, HTCondor, LSF, SLURM, or Torque).  Cactus creates a very large number of small jobs, which can overwhelm these systems.  The work-around described [here](#running-the-step-by-step-workflow-direclty-in-toil) for clusters with large compute nodes available must be used instead.  **Update**:  Cactus version >= 2.0 with GPU enabled will spawn far fewer jobs which, in theory, should make this less of an issue.

### Running on the cloud
Cactus supports running on AWS, Azure, and Google Cloud Platform using [Toil's autoscaling features](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html). For more details on running in AWS, check out [these instructions](doc/running-in-aws.md) (other clouds are similar).

### Running step by step

#### Printing a list of commands to run locally

Breaking Cactus up into smaller jobs can be practical, both for development and debugging, and managing larger workflows.  Here is an example of how to break the Evolver Mammals example up into three steps: 1) Preprocessing 2) Blast 3) Multiple Aligment:
```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore
```

It will print the sequence of commands to run the alignment step-by-step.  Blocks of commands within each alignment run can be run in parallel

`cactus-prepare` can also be used to simplify preprocessing sequences without decomposing the remaining workflow:

```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore --preprocessOnly
```

#### Creading a WDL script to run on Cromwell or Terra

The `--wdl` option in `cactus-prepare` can be used to generate a bespoke [WDL](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md) script for running the alignment from the input seqFile.  Here is an example on how to run locally in [Cromwell](https://github.com/broadinstitute/cromwell)
```
cactus-prepare examples/evolverMammals.txt --wdl > evolver.wdl
wget https://github.com/broadinstitute/cromwell/releases/download/49/cromwell-49.jar
javac -jar ./cromwell-49.jar run evolver.wdl

```

To run on [Terra](https://terra.bio/), use the `--noLocalInputs` option to make sure no local files are embedded in the script.  Also, care must be taken to specify some minimum resource requirements.

```
cactus-prepare examples/evolverMammals.txt --wdl --noLocalInputs --alignCores 2 --defaultMemory 16G > evolver_terra.wdl

```

Then in Terra's [workspace menu](https://app.terra.bio/#workspaces):
* Create a new workspace if necessary with the "+" button
* Click on the workspace
* Click on the "DATA" tab in the workspace menu and use the "Files" link to upload `examples/evolverMammals.txt` to Goggle Cloud
* Click on the "WORKFLOWS" tab
* Click the "+" button to add a workflow
* Click the link in the bottom right to the "Broad Methods Repository"
* Click the "Create New Method... +" button
* Choose and namespace and name, then either upload or paste `evolver_terra.wdl` as created above and click "Upload"
* If this WDL is valid, you can use the "Export To Workspace" button to link it to the Terra Workspace (using a blank configuration)
* You can select the option to go back to the Terra Workspace, otherwise the workflow should now appear as a card in the Terra "workflows" tab the next time you navigate there or refresh
* To run it, click the workflow then click the "INPUTS" tab, and select the `evolverMammals.txt` file in the Attribute field for Task=`cactus_prepare` Variable=`prep_seq_file`
* Tick "Run workflow with inputs defined by file paths"
* Save and click "RUN ANALYSIS"

In the evolver example, all input sequences are specified in public URLs.  If sequences are not specified as URLs in the seqfile, then they must be uploaded in similar fashion to how the evolverMammals.txt was uploaded and selected in the example above.

Here is an example of some settings that have worked on a mammalian-sized genome alignment on Terra:

```
cactus-prepare --wdl mammals.txt --noLocalInputs --preprocessBatchSize 5 --alignDisk 3000G --halAppendDisk 3000G --preprocessDisk 3000G --defaultDisk 1000G --defaultCores 64 --gpu --gpuCount 8 --defaultMemory 385G > mammals.wdl

```

If the workflow fails for whatever reason, it can be edited (to, say, increase job requirements) then resumed as follows:
* In the Workflows tab, click the scripts link beside "Source:" to go back to the Firecloud page to edit the WDL script
* Edit it and "Save a New Snapshot"
* Back in the Terra Workflows tab for the workflow, refresh the page, and select the new snapshot from the "Snapshots" menu.
* Click the "Save" button, ensure that "Use call caching" is ticked, then "Run Analysis" again to resume the workflow.

#### Running the step-by-step workflow direclty in Toil

`cactus-prepare-toil` shares the interface of `cactus-prepare` except instead of printing the command lines or WDL script, it runs them directly from Toil.  An example use case of this is within UCSC's kubernetes cluster.  Like many computing environments, the number of jobs that can be scheduled is limited, so running `cactus` directly using Toil's `kubernetes` batch system will swamp the cluster.  But if the computation can be broken up into a handful of steps, and a job is only created for each step (as in the Cromwell/WDL method), then it can run through.  So `cactus-prepare-toil` will run as a high-level Toil workflow on the specified batch system, and it will launch jobs for each command (`cactus-preprocess, cactus-blast, cactus-align`), and each one of these jobs will get scheduled on a node and run its command with the `singleMachine` batch system.  Here is an example invocation for kubernetes:

```
cactus-prepare-toil aws:us-west-2:<JOBSTORE-NAME> examples/evolverMammals.txt --binariesMode singularity --batchSystem kubernetes --realTimeLogging --outHal s3://<BUCKET-NAME>/out.hal --defaultDisk 20G --defaultMemory 12G --defaultCores 4
```

### Pangenome Pipeline

[The Cactus Pangenome Pipeline is described here](doc/pangenome.md)

## GPU Acceleration

[SegAlign](https://github.com/ComparativeGenomicsToolkit/SegAlign), a GPU-accelerated version of lastz, can be used in the "preprocess" and "blast" phases to speed up the runtime considerably, provided the right hardware is available. Unlike lastz, the input sequences do not need to be chunked before running SegAlign, so it also reduces the number of Toil jobs substantially.  The [GPU-enabled Docker releases](https://github.com/ComparativeGenomicsToolkit/cactus/releases) have SegAlign turned on by default and require no extra options from the user.  Otherwise, it is possible to [manually install it](https://github.com/ComparativeGenomicsToolkit/SegAlign#-dependencies) and then enable it in `cactus` using the `--gpu` command line option. One effective way of ensuring that only GPU-enabled parts of the workflow are run on GPU nodes is on Terra with `cactus-prepare --gpu --wdl` (see above example).

For mammal-sized genomes, we've tested SegAlign with 8 V100 GPUs / 64 CPU cores / 400G RAM. 

Please [cite SegAlign](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043).
 

## Using the output
Cactus outputs its alignments in the [HAL](https://github.com/ComparativeGenomicsToolkit/hal) format. This format represents the alignment in a reference-free, indexed way, but isn't readable by many tools. To export a MAF (which by its nature is usually reference-based), you can use the `hal2maf` tool to export the alignment from any particular genome: `hal2maf <hal> --refGenome <reference> <maf>`.

You can use the alignment to generate gene annotatations for your assemblies, using the [Comparative Annotation Toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit).

You can also [convert the HAL alignment into a Pangenome Graph](https://github.com/ComparativeGenomicsToolkit/hal#pangenome-graph-export-gfa-and-vg).  `hal2vg` is now included in the Cactus Docker images and binary release.

Please [cite HAL](https://doi.org/10.1093/bioinformatics/btt128).

## Updating existing alignments
Cactus supports incrementally updating existing alignments to add, remove, or update genomes. The process involves minor surgery on the output HAL files. See [this document](doc/updating-alignments.md) for details.
# Frequently Asked Questions
Q: I'm running under macOS using the Docker functionality and get an error from Docker: `docker: Error response from daemon: Mounts denied: [...]`

A: Go to your Docker preferences. In the "File Sharing" tab, double-click the last entry ("/path/to/exported/directory") and type in `/var/folders`. (Don't use the `+` button, it won't work because it resolves symlinks before adding).

The reason you have to do this is that the Docker VM requires explicitly listing the directories that can be bind-mounted. The default temp directory on macOS (`/var/folders/...`) is *symlinked* to a directory that is already listed as bind-mountable, but Docker checks the listing before resolving the symlink, returning an error.
The Cactus Pangenome Pipeline
===

## Introduction

[Cactus](../README.md) uses a phylogenetic tree as a guide in order to progressively create multiple alignments.  This heuristic allows Cactus scale linearly with the number of input genomes, by decomposiing the work into one alignment per internal (ancestral) node of the tree.  If the guide tree is fully resolved (binary), only two genomes (plus up to three outgroups) are aligned in each subproblem.

Progressively aligning up a guide tree makes sense when the evolution of the input genomes can be explained by a tree.  It is robust to small errors in the tree, as well as small numbers of non-treelike events (ex incomplete lineage sorting or horizontal genen transfer), making it [suitable for alignments of different vertebrate species](https://doi.org/10.1038/s41586-020-2871-y).

But the tree-like assumption breaks down when considering an alignment of individuals from the *same species*.  Such within-population genome alignments are increasingly in demand as high-quality assemblies become more available (ex: [HPP](https://humanpangenome.org/)), given their potential to better identify and represent structural variation than more traditional reference-based re-sequencing approaches.

The Cactus Pangenome Pipeline adapts [Cactus](../README.md) to no longer rely on a guide tree, by taking advantage of the relative similarity of the input sequence to use minimizer sketches to determine initial anchors, then partial order alignments to refine them.  It also provides the options to generate output in standard pangenome graph formats such as [vg](https://github.com/vgteam/vg) and [GFA](https://github.com/GFA-spec/GFA-spec), in addition to the usual HAL. 

*This is a work in progress and is not yet published.* 

## Overview

The interface is very similar to the [default utilization of Cactus](../README.md), and is dependent on a [seqfile mapping genome names to fasta locations](seqFile-the-input-file).  The main difference here is that a tree need not be provided at the top of the seqfile.  If a tree is present, it must be a star tree (all leaves connected to one root).

Also, a [minigraph](https://github.com/lh3/minigraph) GFA is required.  It can be constructed by running `minigraph -xggs`.  It is strongly suggested but not required to use the fasta files from the seqfile as input here.  Note that the first sequence passed to minigraph will be considered the "reference", and its paths will be acyclic in the graph output.  

The following two Cactus commands are run to produce an alignment and pangenome graph from the minigraph GFA and fasta input:

1. `cactus-graphmap`: Align each input fasta sequence to the minigraph
2. `cactus-align --pangenome --pafInput --outVG`: Run cactus in pangenome mode to produce a HAL alignment and vg graph from the minigraph alignments.

The pipeline can be run without the minigraph GFA by using `cactus-refmap` (which will require identifying one of the inputs as the reference) instead of `cactus-graphmap`, but this will at some cost to sensitivity.

While lastz preprocessing with `cactus-preprocess` is strongly recommended for the Progressive Cactus pipeline, it is not necessary here.  But masking centromeres with the new `cactus-preprocess --alphaMask` can be useful in some cases.  

`cactus-align` does not yet scale to whole human genomes, but this can be worked around by decomposing into chromosomes as [described in the example below](hprc-graph-construction).

**Important**: If you built from source, you must run `build-tools/downloadPangenomeTools` in order to install the extra binary dependencies required for the pangenome pipeline.

## Evolver Simulation Example

This is a very small example along the same lines as Cactus's "evolverMammals" test.  Begin by constructing the minigraph:
```
# download the fasta sequences
for i in `awk '{print $2}' examples/evolverPrimates.txt; do wget $i; done
# make the minigraph with human as the reference
minigraph -xggs simHuman.chr6 simChimp.chr6 simGorilla.chr6 simOrang.chr6 > primates.gfa
```

Align the sequences back to the graph.  Note that the `--outFasta` option is required.  It will be used to update the seqfile with an entry for the minigraph node sequences.  
```
cactus-graphmap ./jobstore ./examples/evolverPrimates.txt primates.gfa primates.paf --outputFasta primates.gfa.fa --realTimeLogging
```

Create the cactus alignment from the seqfile and PAF, and export the output to vg.
```
cactus-align ./jobstore ./examples/evolverPrimates.txt primates.paf primates.hal --pangenome --pafInput --realTimeLogging --outVG --reference simChimp
```

## HPRC Chromosome-by-Chromosome Graph Construction

This [script](https://github.com/glennhickey/pg-stuff/blob/main/cactus-pangenome.sh) was used to make the "year 1" minigraph-cactus pangenome graphs on AWS. 
```
./cactus-pangenome.sh -j aws:us-west-2:glennhickey-jobstore -s ./hprc-year1-f1g.fix.HG02080.1.brnn.leaveout.seqfile -m ftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPRC-f1g/GRCh38-f1g-90.gfa.gz  -o s
3://vg-k8s/vgamb/wg/cactus/GRCh38-f1g-90/aug11 -n GRCh38-f1g-90-mc-aug11 -r GRCh38 -d s3://vg-k8s/vgamb/wg/fasta/hs38d1.decoys.only.vg  -g  -F  -C -M 100000 -K 10000  2>> stderr.aug11.log > /dev/null
```

The steps are described in a bit more detail below:

### Input

*Todo: replace with public-facing links once available*

* `seqfile` containing links to gzipped fasta files for 90 haplotypes from `2020AUG26` bucket + hg38 (no alts) + [CHM13](https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz)
* minigraph GFA `GRCh38-freeze1.gfa.gz`

### Startup

Make a node on AWS with `toil launch-cluster --leaderNodeType t2.medium --leaderStorage 32 ...` and set up cactus as [described here](./running-in-aws.md).

### Preprocessing (about 7 hours)

In order to limit spurious alignments and graph complexity in difficult regions, it is recommended to mask these regions out from the get-go.  In an abundance of caution, we use two methods to detect these regions.  First, align the sequences to create a PAF.

```
cactus-graphmap <jobstore> seqfile.txt <minigraph GFA> s3://<bucket>/GRCh38-freeze1-orig.paf --realTimeLogging --outputFasta s3://<bucket>/GRCh38-freeze1.gfa.fa --logFile graphmap-1.log --batchSystem mesos --provisioner aws --nodeTypes r3.8xlarge:0.7 --maxNodes 20 --defaultPreemptable --betaInertia 0 --targetTime 1 
```

Next, we combine the coverage gaps in the PAF with dna-brnn to produce the final masking, using a 100kb length threshold (default in cactus config).  Note that we need to create a `seqfile.masked.txt` seqfile specifying the locations of the masked fasta sequences.  

```
cactus-preprocess <jobstore> seqfile.txt seqfile.masked.txt --realTimeLogging --logFile preprocess.log  --batchSystem mesos --provisioner aws --nodeTypes r3.8xlarge:0.7 --maxNodes 25 --defaultPreemptable --betaInertia 0 --targetTime 1  --maskFile s3://<bucket>/GRCh38-freeze1-orig.paf  --maskAlpha --brnnCores 8
```

Note that instead of softmasking sequences and leaving them unaligned, they can be clipped out entirely using `--clipAlpha` instead of `--maskAlpha` above.

### Map contigs to minigraph (about 2 hours)

Use `cactus-graphmap` to do the mapping.  It will produce a PAF file of pairwise alignments for cactus to build on.  It will be much faster this time than above, as it will ignore the masked sequences.  A future version of Cactus will clean up the interface in order to make this second call unnecessary! 

```
cactus-graphmap <jobstore> `seqfile.masked.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --maskFilter 100000 --realTimeLogging --logFile graphmap2.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeTypes r3.8xlarge:0.7 --maxNodes 20 --outputFasta s3://<bucket>/GRCh38-freeze1.gfa.fa
```

### Split the sequences by reference chromosome (about 3 hours)

Cactus hasn't yet been tested on the entire input set at once.  In the meantime, we split into chromosomes and make a graph for each.  All output will be written in bucket specified by `--outDir`.  A separate cactus project will be written for each chromosome (`--refContigs`).  The decomposition is determined from mapping: an input contig maps to the reference chromosome to which it has the most alignment in the PAF. Input contigs that cannot be assigned confidently to a reference chromosome will be flagged as `_AMBIGUOUS_`. It's currently important to use `--reference hg38` in order for it not to be filtered out with ambiguity filter.

```
cactus-graphmap-split <jobstore> `seqfile.masked.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --refContigs $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr$i; done) --reference hg38 --maskFilter 100000  --outDir s3://<bucket>/chroms --realTimeLogging --logFile graphmap-split.log --batchSystem mesos --provisioner aws --defaultPreemptable --betaInertia 0 --targetTime 1 --nodeTypes r3.8xlarge:0.7 --maxNodes 25
```

### Make a pangenome graph for each chromosome (about 10 hours)

Begin by downloading the chromfile and seqfiles that were created above (they must be input from local folders)
```
pip install awscli
aws s3 sync s3://<bucket>/chroms/seqfiles ./seqfiles
aws s3 cp s3://<bucket>/chroms/chromfile.txt ./
```

Use `cactus-align-batch` to align each chromosome on its own aws instance.  The `--alignCoresOverrides` option is used to ensure that the larger chromosomes get run on bigger instances.

```
cactus-align-batch  <jobstore> ./chromfile.txt s3://<bucket>/align-batch --alignCores 32  --alignOptions "--pafInput --pangenome --outVG --barMaskFilter 100000 --realTimeLogging --reference GRCh38"  --batchSystem mesos --provisioner aws --defaultPreemptable  --nodeTypes r3.8xlarge:0.7 --nodeStorage 1000 --maxNodes 25 --betaInertia 0 --targetTime 1 --logFile align-batch.log --realTimeLogging
```

### Combining the output into a single graph (about 5 hours)

The chromosome graphs can then be merged with `cactus-graphmap-join`:
```
cactus-graphmap-join <jobstore> --outDir s3://<bucket/join --outName GRCh38.mc --reference GRCh38 --realTimeLogging --rename "CHM13>CHM13.0" --clipLength 100000 --wlineSep . --vg {list of s3 vg paths output above} --hal {list of s3 hal paths output above} --nodeTypes r4.16xlarge--nodeStorage 1000 --maxNodes 1 --logFile align-batch.log --realTimeLogging
```


Running Cactus on AWS
===
Cactus supports running on AWS with an auto-scaling cluster using [Toil](https://toil.readthedocs.io/en/latest/). Check out the [Toil docs on running in the cloud](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html) for the full story, but here's a short walkthrough of running on AWS.

## AWS setup
If you have a fresh AWS account or haven't used EC2 before, you'll need to go through some initial setup.
### Keypair
Make sure you have an AWS keypair ready. [This document](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html) will tell you how to create an AWS keypair that will allow you to log into the instances you create.
### Access keys
You'll also need to have your AWS access credentials set up in `~/.aws/credentials` or the typical AWS environment variables (`AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`). [Here is AWS's documentation on setting up your access key](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html), and [this guide](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html) will help you set up `~/.aws/credentials`.
### Instance limits
By default, AWS will restrict you to running only a few small instances at a time. If you have a new AWS account, or you're not sure what your limits are, you will probably need to increase them. See [this guide](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-resource-limits.html) for information on how to check your existing limits and how to increase them if necessary. AWS support takes only one or two business days to process your request if you have the default "basic" support package.

[Look below](#estimating-the-maximum-number-of-worker-instances-youll-need) for tips on how many instances you may need for your specific alignment. Keep in mind that spot-market instance limits are separate from "on-demand" (non-preemptable) instance limits. (You'll probably want to request a slightly higher limit than you think you need, just in case you want to tweak the number of instances later on.)
## Installing Toil on your local machine
Follow the steps in the README, making sure to install toil with its extra AWS support:
```
git clone https://github.com/comparativegenomicstoolkit/cactus.git
cd cactus
virtualenv -p python 3.8 venv
source venv/bin/activate
pip install -r toil-requirement.txt
```
## Estimating the maximum number of worker instances you'll need
The cluster will automatically scale up and down, but you'll want to set a maximum number of nodes so the scaler doesn't get overly aggressive and waste money, or go over your AWS limits. We typically use `c4.8xlarge` on the spot market for most jobs, and `r4.8xlarge` on-demand for database jobs. Here are some very rough estimates of what we typically use for the maximum of each type (round up):

- `N` mammal-size genomes (~2-4Gb): `(N / 2) * 20` c4.8xlarge on the spot market, `(N / 2)` r3.8xlarge on-demand
- `N` bird-size genomes (~1-2Gb): `(N / 2) * 10` c4.8xlarge on the spot market, `(N / 4)` r3.8xlarge on-demand
- `N` nematode-size genomes (~100-300Mb): `(N / 2)` c4.8xlarge on the spot market, `(N / 10)` r3.8xlarge on-demand
- For anything less than 100Mb, the computational requirements are so small that you may be better off running it on a single machine than using an autoscaling cluster.
## Launching the "leader" instance
Make sure you have your AWS keypair active in your `ssh-agent`, unless you used your existing SSH key. 
If `ssh-agent` isn't started by your operating system, you may have to run `eval $(ssh-agent)` to start it. You can activate the AWS keypair for a session by running:
```
ssh-add path/to/your_aws_ssh_keypair
```


Then launch the leader like so:

```
toil launch-cluster -z us-west-2b <clusterName> --keyPairName <yourKeyPairName> --leaderNodeType t2.medium
```
## Transfer over local input data (or use URLs in your seqfile)
You need to get your actual data to the leader somehow. You can use `http://` and `s3://` to specify FASTA locations in your input file, or you can rsync your data over like so:
```
toil rsync-cluster -z us-west-2b my-cactus-cluster -avP seqFile.txt input1.fa input2.fa :/
```
## Log into the leader and set up the cluster's Cactus installation
Log into the leader by running:
```
toil ssh-cluster -z us-west-2b <clusterName>
```

You should now get a prompt like:
```
root@ip-172-31-34-148:/#
```
indicating that you're on the leader.

Install Cactus in a virtual environment on the leader:
```
apt update
apt install -y git tmux
virtualenv --system-site-packages -p python3.8 venv
source venv/bin/activate
git clone https://github.com/comparativegenomicstoolkit/cactus.git --recursive
cd cactus
pip install --upgrade .
cd /
```
## Run the alignment (on the leader)
The key parameters you'll care about changing (besides the usual Cactus parameters) are the autoscaling parameters `--nodeTypes`, `--minNodes`, `--maxNodes` and the jobStore location, `aws:<region>:<jobStoreName>`.

You must use the AWS jobstore (or other cloud jobstore, though others may incur data egress charges), not a directory jobstore, because there is no shared filesystem in the cluster. Set the region to whatever region you're running the leader in. The jobStoreName must be globally unique.

The `--nodeTypes` option lets you specify the list of instances you want as well as the price you're willing to pay for spot instances. For example, `c4.8xlarge:0.6` says that we want a c4.8xlarge instance on the spot market, and we're willing to pay up to 60 cents an hour for it. A value like `r3.8xlarge` indicates that we also want on-demand r3.8xlarge instances (for which we pay exactly the on-demand price, which is usually substantially higher). The `--maxNodes` option will let you specify a list containing the maximum number of nodes of each type to use at any given time (in the same order as `--nodeTypes`. The `--minNodes` option should probably be left at 0 unless you really know what you're doing and you're having substantial difficulties with the autoscaler.

```
cactus --nodeTypes c4.8xlarge:0.6,r3.8xlarge --minNodes 0,0 --maxNodes 20,2 --provisioner aws --batchSystem mesos --metrics aws:us-west-2:<jobstoreName> seqFile.txt output.hal
```

This will take a while. You'll want to run this inside something that will preserve your session, like `tmux` or `screen`, so the command doesn't terminate when you disconnect.

## Restarting after failure
If you cancel your run, or it fails for some reason, you can start where you left off by running:
```
cactus [all your usual options...] --restart
```
## Shut down your leader
When the alignment is done, your leader will still be active, costing you a little bit of money every day. Make sure you get rid of it when you're done:
```
toil destroy-cluster -z us-west-2b <yourClusterName>
```
# Updating Cactus alignments
Because of the hierarchical structure of the alignments Cactus produces, you can add, remove, or update genomes to an existing alignment without recomputing the entire alignment from scratch. This document briefly outlines the steps required. Ensure you have the [HAL toolkit](https://github.com/ComparativeGenomicsToolkit/cactus) installed; it comes with alignment-modification tools which we'll be using.

**Important: back up your alignment file before attempting any of these operations. If something goes wrong, the file may be completely unrecoverable.**
## Deleting a genome
Deleting a genome is easy. Just run:
```sh
halRemoveGenome <hal file> <genome to delete>
```

The parent ancestral genome of the genome you deleted will still hang around, since it's needed to establish the alignment relationships between the other genomes in the parent's subtree and its supertree. Only leaf genomes can be deleted, though reconstructed genomes may be deleted if all their children have already been deleted (since in that case they would now be a leaf genome).

The HAL file will remain the same size, since the genome data offsets cannot easily be shifted. If you need a smaller HAL file, check out the documentation for the `h5repack` command, or use the `halExtract` tool from the HAL toolkit (which is a bit slower than `h5repack`.
## Adding a new genome
Adding a new genome can be fairly complicated. The first thing to do is to determine *how* the genome should be added into the species tree: by splitting an existing branch (**add-to-branch**) or by adding the new genome as a child of an existing ancestral genome (**add-to-node**). These two ways allow you to maintain whatever species tree shape you want. Adding a genome to a node is slightly easier (it involves creating only one alignment rather than two), but may be more expensive since the alignment involves more genomes than in the add-to-branch case. Check out this handy diagram, showing the different ways that a genome "6" can be added:

![Different ways of adding a genome](add-genome-fig-github.png)

### Adding to a branch
Splitting a branch involves creating two alignments. This involves creating two separate Cactus alignments, one after the other. You'll have to start with the bottom half, since you're inferring a new ancestor. Here's what that might look like, where `genome6`, `genome7` etc. reflect 6, 7, etc. from the diagram (your names and FASTA locations should reflect the genomes in your tree):

"bottom_half.txt"
```
((genome6, genome4)genome7, genome3)genome5;

genome6 genome6.fa
genome4 genome4.fa
[...]
```

Then, run that alignment as normal, i.e. `cactus jobStore bottom_half.txt bottom_half.hal --root genome7 [--options]`. The `--root genome7` option is included here because we only want genome3 as outgroup information to help us infer genome7 more accurately.

After that's complete, you can work on the top half. First, you'll need to extract the FASTA file for the genome you just inferred:

```
hal2fasta bottom_half.hal genome7 > genome7.fa
```

Then create a file that looks something like (again, filling in whatever makes sense in your tree):

"top_half.txt"
```
(genome7, genome3)genome5;

genome7 genome7.fa
genome3 genome3.fa
genome5 genome5.fa
```

Use that to run Cactus again, creating `top_half.hal`.

Now you can modify your HAL file (ensuring you've backed it up):
```
halAddToBranch <your original alignment hal file> bottom_half.hal top_half.hal genome5 genome7 genome4 genome6 <length of genome5-genome7 branch> <length of genome7-genome6 branch>
```
which isn't exactly obvious. The meaning behind specifying these genomes is to orient the process, showing which genomes actually need to be updated or inserted in which places.

The branch lengths don't impact the alignment at this stage, they only impact the Newick tree that HAL displays.
### Adding to a node
Adding a genome as a new child of an existing ancestor requires creating only a single alignment (though that alignment will be more computationally expensive than usual, since it will involve more genomes than in the add-to-branch case). Again following our example, the sequence file should look something like this, with names and paths modified to suit your specific situation:

"add_to_node.txt"
```
(genome4, genome6, genome3)genome5;

genome6 genome6.fa
genome5 genome5.fa
[...]
```

Run this alignment as you would normally to create "add_to_node.hal". Then you can run:
```
halReplaceGenome alignment_to_add_to.hal --bottomAlignmentFile add_to_node.hal --topAlignmentFile alignment_to_add_to.hal genome5
```

After that, `alignment_to_add_to.hal` should contain the new genome.
## Replacing a genome
A genome can be replaced (for example to update an assembly version) by removing it and then following the "add-to-node" procedure to add the new version back as a child of its parent.
## Validating the file
It's a good idea to run a validation against the resulting file to make sure the operation went OK. Run:
```sh
halValidate --genome <genome you modified> <hal file>
```
If it doesn't raise an error, you should be good to go! If it does raise an error, please raise an issue on the [HAL github page](https://github.com/ComparativeGenomicsToolkit/hal).
