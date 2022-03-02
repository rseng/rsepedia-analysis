Our hope is that octopus will become a community driven project, so please feel tree to contribute! 

## General guidance

Before contributing please read the documentation, check existing issues and branches. Please ask questions if you're unsure of anything.

We follow [this](http://nvie.com/posts/a-successful-git-branching-model/) git branch model. In particular:

* The `master` branch is for stable versions only. Usually only tagged releases will be commited to `master`.
* The `develop` branch is the main working branch. Only small changes should be commited directly to this branch.
* New major features should branch from `develop` and be named like `feature/name`.
* Minor bug fixes should branch from the relevant branch and be named like `fix/name`.
* Experimental work should branch from the relevant branch and be named like `exp/name`.

## Style

Please try to adhere to the existing code style, in particular:

#### Naming

* Variable names are all lower case and seperated with underscores (e.g. `foo_bar`).
* Type names (including template types and aliases) are UpperCamelCase (e.g. `FooBar`).
* Function names are all lower case and seperated_with_underscores (e.g. `foor_bar()`).
* `enum` and `enum class` members are all lower case and seperated with underscores (e.g. `foo_bar`).
* `constexpr` variables are lowerCamelCase (e.g. `fooBar`).
* Global variables are lowerCamelCase (e.g. `fooBar`).

#### Brackets

* Classes and functions start with an open bracket on the next line, e.g:

```cpp
class Foo
{
    // code
};

int bar(int x)
{
    return 0;
}
```

* Loops and if statements open the bracket on the same line, e.g:

```cpp
for (;;) {
    // stuff
}

if (something) {
    // stuff
} else {
    // other stuff
}
```

* Small `if` statements can appear on the same line if it improves readability:

```cpp
if (something) return true;
```

#### Namespaces

* All code goes in the `octopus namespace`.
* The open bracket appears on the same line as the `namespace` name (e.g. `namespace octopus {`).

#### Misc

* Use 4 spaces for tab indentation.
* There is no strict limit for line lengths; use best judgement. 
![Octopus Logo](logo.png)

[![Website](https://img.shields.io/website?url=https%3A%2F%2Fluntergroup.github.io%2Foctopus%2F)](https://luntergroup.github.io/octopus/)
[![Build Status](https://travis-ci.org/luntergroup/octopus.svg?branch=master)](https://travis-ci.org/luntergroup/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![Gitter](https://badges.gitter.im/octopus-caller/Lobby.svg)](https://gitter.im/octopus-caller/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/luntergroup/octopus)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/octopus/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/dancooke/octopus?label=docker)](https://hub.docker.com/r/dancooke/octopus)

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly consider all possible haplotypes at a given loci in reasonable time.

There are currently six calling models available:

- [individual](https://luntergroup.github.io/octopus/docs/guides/models/individual): call germline variants in a single healthy individual.
- [population](https://luntergroup.github.io/octopus/docs/guides/models/population): jointly call germline variants in small cohorts.
- [trio](https://luntergroup.github.io/octopus/docs/guides/models/trio): call germline and _de novo_ mutations in a parent-offspring trio.
- [cancer](https://luntergroup.github.io/octopus/docs/guides/models/cancer): call germline and somatic mutations tumour samples.
- [polyclone](https://luntergroup.github.io/octopus/docs/guides/models/polyclone): call variants in samples with an unknown mixture of haploid clones, such a bacteria or viral samples.
- [cell](https://luntergroup.github.io/octopus/docs/guides/models/cell): call variants in a set of single cell samples from the same individual.

Octopus calls SNVs, small-medium sized indels, and small complex rearrangements in [VCF 4.3](https://luntergroup.github.io/octopus/docs/guides/advanced/vcf).

## Quick start

Install Octopus:

```shell
$ git clone https://github.com/luntergroup/octopus.git
$ octopus/scripts/install.py --dependencies --forests
$ echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile
$ source ~/.bash_profile
```

Call some variants:

```shell
$ FOREST="$(pwd)/octopus/resources/forests/germline.v0.7.4.forest"
$ octopus -R hs37d5.fa -I NA12878.bam -T 1 to MT -o NA12878.octopus.vcf.gz --forest $FOREST --threads 8
```

## Documentation

Documentation is hosted on [GitHub pages](https://luntergroup.github.io/octopus/).

## Support

Please report any bugs or feature requests to the [octopus issue tracker](https://github.com/luntergroup/octopus/issues). General chat is hosted on [Gitter](https://gitter.im/octopus-caller/Lobby).

## Contributing

Contributions are very welcome! Please consult the [contribution guidelines](CONTRIBUTING.md).

## Authors

Daniel Cooke and Gerton Lunter

## Citing

See [publications](https://luntergroup.github.io/octopus/docs/publications) associated with Octopus.

## License

Octopus is distributed under the [MIT LICENSE](LICENSE).
---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**Version**

```shell
$ octopus --version
// PASTE OUTPUT HERE
```

**Command**
Command line to install octopus:
```shell
$ 
```

Command line to run octopus:
```shell
$ octopus
```

**Additional context**
Add any other context about the problem here, e.g.

 - Reference [e.g. hg19]. Please provide link to fasta if atypical.
 - BAM files (if possible).
This folder contains all of the tests for octopus. These can be divded into three categories:

NOTE: Many of the tests use real data. In order to run the tests the files specified in 'test_common.h' must be present in your system.

1. Component unit tests: these tests cover functionality requirments of the major components of octopus. They are designed to ensure expected functionality, especially at edge cases, and avoid common bugs (e.g. off-by-one errors). Note many of the tests here are run on real data.
2. Benchmarks: these tests contain benchmarks for various key components. Generally these are tests that have directed design decisions (e.g. using virtual methods).
3. Data: these are tests on real data, usually 1000G. They are designed to measure and improve calling performance.
# Website

This website is built using [Docusaurus 2](https://docusaurus.io/), a modern static website generator.

## Installation

```console
yarn install
```

## Local Development

```console
yarn start
```

This command starts a local development server and opens up a browser window. Most changes are reflected live without having to restart the server.

## Build

```console
yarn build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

## Deployment

```console
GIT_USER=<Your GitHub username> USE_SSH=true yarn deploy
```

If you are using GitHub pages for hosting, this command is a convenient way to build the website and push to the `gh-pages` branch.
---
title: Markdown page example
---

# Markdown page example

You don't need React to write simple standalone pages.
---
id: publications
title: Publications
---

Here are the current papers associated with Octopus.

## [A unified haplotype-based method for accurate and comprehensive variant calling](https://www.nature.com/articles/s41587-021-00861-3)

This is the main Octopus paper describing the method with benchmarks for germline and somatic variant calling.

#### Abstract

Almost all haplotype-based variant callers were designed specifically for detecting common germline variation in diploid populations, and give suboptimal results in other scenarios. Here we present Octopus, a variant caller that uses a polymorphic Bayesian genotyping model capable of modeling sequencing data from a range of experimental designs within a unified haplotype-aware framework. Octopus combines sequencing reads and prior information to phase-called genotypes of arbitrary ploidy, including those with somatic mutations. We show that Octopus accurately calls germline variants in individuals, including single nucleotide variants, indels and small complex replacements such as microinversions. Using a synthetic tumor data set derived from clean sequencing data from a sample with known germline haplotypes and observed mutations in a large cohort of tumor samples, we show that Octopus is more sensitive to low-frequency somatic variation, yet calls considerably fewer false positives than other methods. Octopus also outputs realigned evidence BAM files to aid validation and interpretation.

#### Bibtex

```tex
@article{octopus,
   author = {Cooke, Daniel P. and Wedge, David C. and Lunter, Gerton},
   title = {A unified haplotype-based method for accurate and comprehensive variant calling},
   journal = {Nature Biotechnology},
   ISSN = {1546-1696},
   DOI = {10.1038/s41587-021-00861-3},
   url = {https://doi.org/10.1038/s41587-021-00861-3
https://www.nature.com/articles/s41587-021-00861-3.pdf},
   year = {2021},
   type = {Journal Article}
}
```

## [Benchmarking small-variant genotyping in polyploids](https://www.biorxiv.org/content/10.1101/2021.03.29.436766v1)

In this paper, we benchmark Octopus on polyploid samples.

#### Abstract

Genotyping from sequencing is the basis of emerging strategies in the molecular breeding of polyploid plants. However, compared with the situation for diploids, where genotyping accuracies are confidently determined with comprehensive benchmarks, polyploids have been neglected; there are no benchmarks measuring genotyping error rates for small variants using real sequencing reads. We previously introduced a variant calling method – Octopus – that accurately calls germline variants in diploids and somatic mutations in tumors. Here, we evaluate Octopus and other popular tools on whole-genome tetraploid and hexaploid datasets created using in silico mixtures of diploid Genome In a Bottle samples. We find that genotyping errors are abundant for typical sequencing depths, but that Octopus makes 25% fewer errors than other methods on average. We supplement our benchmarks with concordance analysis in real autotriploid banana datasets.

#### Bibtex

```tex
@article{octopus_polyploid,
   author = {Cooke, Daniel P. and Wedge, David C. and Lunter, Gerton},
   title = {Benchmarking small-variant genotyping in polyploids},
   journal = {bioRxiv},
   pages = {2021.03.29.436766},
   DOI = {10.1101/2021.03.29.436766},
   url = {https://www.biorxiv.org/content/biorxiv/early/2021/03/29/2021.03.29.436766.full.pdf},
   year = {2021},
   type = {Journal Article}
}
```

## [Accurate genotyping of single cells with Octopus](https://www.researchsquare.com/article/rs-583831/v1)

#### Absract

We describe an extension to our variant calling tool, Octopus (https://github.com/luntergroup/octopus), for single-cell DNA sequencing data. Octopus jointly genotypes cells from a lineage, accounting for amplification stochasticity and sequencing error with a haplotype-based Bayesian model. Octopus is considerably more accurate at genotyping single cells than existing methods.

#### Bibtex

```tex
@article{RN862,
   author = {Daniel, Cooke P. and Gerton, Lunter and David C., Wedge},
   title = {Accurate genotyping of single cells with Octopus},
   journal = {Research Square},
   ISSN = {2693-5015},
   DOI = {10.21203/rs.3.rs-583831/v1},
   url = {https://doi.org/10.21203/rs.3.rs-583831/v1},
   year = {2021},
   type = {Journal Article}
}

```---
id: installation
title: Installation
---

Octopus can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested.

The recommend way to install Octopus for most users is:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git
$ octopus/scripts/install.py --dependencies --forests
```

You can then optionally add `octopus` to your `PATH`:

```shell
$ echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile
$ source ~/.bash_profile
```

Then check the installation was successful:

```shell
$ octopus --version
```

## Requirements

* A [C++14](https://isocpp.org/wiki/faq/cpp14) compiler and compatibility standard library. Either [GCC](https://gcc.gnu.org) (version >= 9.3) or [Clang](https://clang.llvm.org) (version >= 11.0) are recommended.
* [Git](https://git-scm.com) version >= 2.5
* [Boost](https://www.boost.org) version >= 1.65
* [htslib](https://github.com/samtools/htslib) version >= 1.4; version != 1.12
* [GMP](https://gmplib.org) version >= 5.1.0
* [CMake](https://cmake.org) version >= 3.9
* Optional:
    * [Python](https://www.python.org) version >= 3 plus the [distro](https://pypi.org/project/distro/) package

:::important

Octopus uses [SIMD](https://en.wikipedia.org/wiki/SIMD) instructions for performance reasons. The instruction set used (minimum [SSE2](https://en.wikipedia.org/wiki/SSE2)) is built statically, so if you compile with [AVX2](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2), you won't be able to use the resulting binary on machines that doesn't support AVX2.

:::

## Python

First clone the git repository in your preferred directory:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git && cd octopus
```

The easiest way to install Octopus from source is with the Python3 installer script. To see the options available to this script run `scripts/install.py --help`.

If all the requirements are accessible on your `PATH` then simply run

```shell
$ scripts/install.py
```

otherwise you can specify paths to each dependency, for example, to set the compiler you'd use 

```shell
$ scripts/install.py --cxx_compiler /path/to/cpp/compiler
```

By default, this installs to `/bin` relative to where octopus is installed. To install to a different location (e.g. `/usr/local/bin`) use:

```shell
$ scripts/install.py --prefix /user/local/bin
```

You can also request all dependencies to be installed locally:

```shell
$ scripts/install.py --dependencies
```

:::tip

If a build isn't working after an update then try adding `--clean` to the install command.

:::

### Setting the build architecture

By default, the binary is optimised for the build machine architecture. If you need to run Octopus on another machine with a different architecture then use the `--architecture` option:

```shell
$ scripts/install.py --architecture haswell
```

This is passed to the [-march](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html) compiler option. 

## CMake

If Python3 isn't available, Octopus can be installed directly with [CMake](https://cmake.org):

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git
$ cd octopus/build
$ cmake .. && make install
```

CMake will try to find a suitable compiler on your system, if you'd like you use a specific compiler use the `-D` option, for example:

```shell
$ cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ ..
```

## Docker

Pre-built Docker images are available on [DockerHub](https://hub.docker.com/r/dancooke/octopus):

```shell
$ docker pull dancooke/octopus
$ docker run dancooke/octopus -h
```

:::important

The Octopus images on DockerHub are built to a [Haswell](https://en.wikipedia.org/wiki/Haswell_(microarchitecture)) architecture. This means that they will only work on  Haswell (with AVX2) or newer machines.

:::

You can also build a new image from the Dockerfile:

```shell
$ git clone https://github.com/luntergroup/octopus.git && cd octopus
$ docker build -t octopus .
$ docker run octopus -h
```

This is especially useful if you need to build to a specific architecture:

```shell
$ docker build -t octopus --build-args architecture=sandybridge .
```

## Singularity

To build a [Singularity](https://singularity.hpcng.org) container directly from the DockerHub images use

```shell
$ singularity build octopus.sif docker://dancooke/octopus
```

## Conda

Octopus is available [pre-built for Linux](https://anaconda.org/bioconda/octopus) as part of [Bioconda](https://bioconda.github.io/):

```shell
$ conda install -c bioconda octopus
```

:::important

The Octopus package on Bioconda is built to a [Haswell](https://en.wikipedia.org/wiki/Haswell_(microarchitecture)) architecture. This means that it will only work on  Haswell (with AVX2) or newer machines. If you need another architecture then consider using [conda-build](https://docs.conda.io/projects/conda-build/en/latest/).

:::---
id: cli
title: Command Line Reference
---

## General options

### `--help`

Command `--help` (short `-h`) prints the list of available command line options.

```shell
$ octopus --help
Octopus command line options:

General:
  -h [ --help ]                         Report detailed option information
  --version                             Report detailed version information
  --config arg                          Config file to populate command line 
                                        options
  --options                             Log all command line option values at 
                                        startup
  --debug [=arg(="octopus_debug.log")]  Create log file for debugging
  --trace [=arg(="octopus_trace.log")]  Create very verbose log file for 
                                        debugging
  -w [ --working-directory ] arg        Sets the working directory
  --threads [=arg(=0)]                  Maximum number of threads to be used. 
                                        If no argument is provided unlimited 
                                        threads are assumed
<OUTPUT STRIPPED>
```

### `--version`

Command `--version` prints the version of the current binary, including the git branch and commit. The command also prints some system information about the binary.

```shell
$ octopus --version
octopus version 0.6.3-beta (develop e7205439)
Target: x86_64 Darwin 18.7.0
SIMD extension: AVX2
Compiler: AppleClang 10.0.1.10010046
Boost: 1_70
```

### `--config`

Option `--config` is used to specify [configuration files](https://github.com/luntergroup/octopus/wiki/How-to:-Use-configs).

```shell
$ octopus ---config myconfig.txt
```


### `--trace`

Option `--trace` results in a log file containing detailed debug information being created.

```shell
$ octopus -R ref.fa -I reads.bam --trace # default name is "octopus_trace.log"
$ octopus -R ref.fa -I reads.bam --trace =/foo/bar/trace.log
```

**Warning** Trace files can get very large, so only use on small inputs..

### `--working-directory`

Option `--working-directory` (short `-w`) is used to set the working directory of the run. All output and temporary files will be relative to the working directory, unless absolute file paths are provided.

```shell
$ octopus -w ~/vcf -o octopus.vcf
```

The output will be in `~/vcf/octopus.vcf`.

**Notes**

* If no working directory is specified then the current directory is used.

### `--resolve-symlinks`

Command `--resolve-symlinks` forces all symbolic link input paths to be replaced with their resolved targets during program initialisation.

```shell
$ octopus -R ref.fa -I reads.bam --resolve-symlinks
```

### `--threads`

Option `--threads` is used to enable [multithreaded execution](https://github.com/luntergroup/octopus/wiki/How-to:-Use-multiple-threads). The option accepts a positive integer argument that specifies the maximum number of threads to use:

```shell
$ octopus -R ref.fa -I reads.bam --threads 4 # use maximum of 4 threads
```

If no argument is provided then automatic thread handling is used. This may be greater than the number of system threads.

```shell
$ octopus -R ref.fa -I reads.bam --threads #automatic thread handling
```

### `--max-reference-cache-memory`

Option `--max-reference-cache-memory` (short `-X`) controls the size of the buffer used for reference caching, and is therefore one way to [control memory use](https://github.com/luntergroup/octopus/wiki/How-to:-Adjust-memory-consumption). The option accepts a non-negative integer argument in bytes, and an optional unit specifier.

```shell
$ octopus -R ref.fa -I reads.bam -X 0 # disable reference caching
$ octopus -R ref.fa -I reads.bam -X 1g # reference cache is 1 gigabyte
```

**Notes**

* Capitisation of the units is ignored.
* An argument of `0` disables reference caching.
* The buffer size never exceeds the size of the reference.

### `--target-read-buffer-memory`

Option `--target-read-buffer-memory` (short `-B`) controls the size of the memory buffer used for input sequencing reads, and is therefore one way to [control memory use](https://github.com/luntergroup/octopus/wiki/How-to:-Adjust-memory-consumption). The option accepts a positive integer argument in bytes, and an optional unit specifier.

```shell
$ octopus -R ref.fa -I reads.bam -B 20G # reference cache is 20 gigabyte
```

**Notes**

* Capitisation of the units is ignored.
* The minimum read buffer size is `50Mb`; arguments less than this are ignored.

### `--target-working-memory`

Option `--target-working-memory` sets the target amount of working memory for computation, and is therefore one way to [control memory use](https://github.com/luntergroup/octopus/wiki/How-to:-Adjust-memory-consumption). The option accepts a positive integer argument in bytes, and an optional unit specifier. The option is not strictly enforced, but is sometimes used to decide whether to switch to lower-memory versions of some methods (possibly at the cost of additional runtime).

```shell
$ octopus -R ref.fa -I reads.bam --target-working-memory 40Gb
```

**Notes**

* Capitisation of the units is ignored.

### `--temp-directory-prefix`

Option `--temp-directory-prefix` sets the Octopus working temporary directory prefix.

```shell
$ octopus -R ref.fa -I reads.bam --temp-directory-prefix ~/octopus_tmp
```

**Notes**

* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.
* If a directory with the prefix already exists, then `-N` is appended to the prefix where `N` is the smallest integer such that the resulting path does not exist.

### `--reference`

Option `--reference` (short `-R`) sets the reference FASTA file used for variant calling. This is a required option.

```shell
$ octopus -R ref.fa -I reads.bam
```

**Notes**

* In addition to the FASTA file, a reference FASTA index file (extension `.fai`) is also required. This must be present in the same directory as the FASTA file.
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.

### `--reads`

Option `--reads` (short `-I`) specifies the input read files (extension `.bam` or `.cram`) used for variant calling. The argument is a list of file paths. This is a required option (unless `--reads-file` is specified).

```shell
$ octopus -R ref.fa -I reads1.bam reads2.bam
```

**Notes**

* Each read file must have a paired index file in the same directory.
* The option can be specified multiple times on the command line (arguments are concatenated).
* Can be used in conjunction with `--reads-file` (arguments are concatenated).
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.


### `--reads-file`

Option `--reads-file` (short `-i`) specifies a file containing a list of input read files (one per line).

```shell
$ octopus -R ref.fa -i bams.txt
```

### `--regions`

Option `--regions` (short `-T`) specifies a list of genomic [regions to call](https://github.com/luntergroup/octopus/wiki/How-to:-Call-targeted-regions).

```shell
$ octopus -R ref.fa -I reads.bam -T chr1 # all of chr1
$ octopus -R ref.fa -I reads.bam -T chr1 chr2 # all of chr1 and chr2
$ octopus -R ref.fa -I reads.bam -T chr1:0-1,000 # first 1,000bp of chr1
$ octopus -R ref.fa -I reads.bam -T chr1:100- # everything after position 100 of chr1
$ octopus -R ref.fa -I reads.bam -T chr1 to chrM # as specified in reference index
```

**Notes**

* If this option is not specified (and neither is `--regions-file`), then all contigs present in the reference genome index are used.
* Commas in the position tokens are ignored.
* The option can be specified multiple times on the command line (arguments are concatenated).

### `--regions-file`

Option `--regions-file` (short `-t`) specifies a file containing a list of genomic regions to call (one per line). The file format can either be plain text (using the same input format as the `--regions` option), or BED format (i.e. tab-separated). GZipped files are also accepted.

```shell
$ octopus -R ref.fa -I reads.bam -t regions.bed
```

### `--skip-regions`

Option `--skip-regions` (short `-K`) specifies a list of genomic regions to ignore during calling.

```shell
$ octopus -R ref.fa -I reads.bam -K X Y MT
```

**Notes**

* The option can be specified multiple times on the command line (arguments are concatenated).

### `--skip-regions-file`

Option `--skip-regions-file` (short `-k`) specifies a file containing a list of genomic regions to skip (one per line). The file format can either be plain text (using the same input format as the `--regions` option), or BED format (i.e. tab-separated). GZipped files are also accepted.

### `--one-based-indexing`

Command `--one-based-indexing` directs the program to read all input regions (i.e. those specified in options `--regions`, `--regions-file`, `--skip-regions`, and `--skip-regions-file`) using 1-based indexing rather than 0-based indexing.

```shell
$ octopus -R ref.fa -I reads.bam -T chr1:100 --one-based-indexing
$ octopus -R ref.fa -I reads.bam -T chr1:99 # equivalent to above
```

### `--ignore-unmapped-contigs`

Command `--ignore-unmapped-contigs` can be used to force execution if there is a mismatch between the input reference genome and the one used to map input reads (as specified in the SAM header). In particular, if there are contigs present in the reference genome but not in the reference used for read mapping.

```shell
$ octopus -R ref.fa -I reads.bam
<INFO> ------------------------------------------------------------------------
<INFO> octopus v0.6.1-beta (develop 9c97253c)
<INFO> Copyright (c) 2015-2019 University of Oxford
<INFO> ------------------------------------------------------------------------
<EROR> A user error has occurred:
<EROR> 
<EROR>     Some or all of the contigs in the reference genome (ref) are not
<EROR>     present in the read files.
<EROR> 
<EROR> To help resolve this error ensure the reference genome used for mapping
<EROR> is the same as the one used for calling (ref) and all input contigs are
<EROR> present in the read headers.
<INFO> ------------------------------------------------------------------------

$ octopus -R ref.fa -I reads.bam --ignore-unmapped-contigs # ignore error above

```

### `--samples`

Option `--samples` (short `-S`) specifies a list of samples to use for variant calling. Each sample must be present in the input read files (using the `SM` SAM tag).

```shell
$ octopus -R ref.fa -I reads.bam -S sample1 sample2
```

**Notes**

* If no samples are specified, all samples present in the input read files are used.
* The option can be specified multiple times on the command line (arguments are concatenated).

### `--samples-file`

Option `--samples-file` (short `-s`) specifies a file containing a list (one per line) of samples to use for variant calling.

```shell
$ octopus -R ref.fa -I reads.bam -s samples.txt
```

### `--output`

Option `--output` (short `-O`) sets the output destination. The file extension is used to determine the output file type (valid types are `.vcf`, `.vcf.gz`, and `.bcf`). If the file type is not recognised then uncompressed VCF is used.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.bcf
```

**Notes**

* If no output is specified, `stdout` is used.
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.

### `--contig-output-order`

Option `--contig-output-order` specifies the order that records will be processed and written to the output. Possible options are: `lexicographicalAscending`, `lexicographicalDescending`, `contigSizeAscending`, `contigSizeDescending`, `asInReferenceIndex`, `asInReferenceIndexReversed`, `unspecified`

```shell
$ octopus -R ref.fa -I reads.bam --contig-output-order asInReferenceIndexReversed
```

### `--sites-only`

Command `--sites-only` removes genotype information from the final output (i.e. drops the VCF `FORMAT` and sample columns).

```shell
$ octopus -R ref.fa -I reads.bam --sites-only
```


### `--bamout`

Option `--bamout` is used to produce [realigned evidence BAMs](https://github.com/luntergroup/octopus/wiki/How-to:-Make-evidence-BAMs). The option input is a file prefix where the BAMs should be written.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --bamout reads.realigned.bam # if 'reads.bam' contains one sample
$ octopus -R ref.fa -I readsA.bam readsB.bam -o calls.vcf --bamout minibams
```

**Notes**

* The `--bamout` command is only allowed if the final output is written to file (i.e. with `--output`).
* The default behaviour of `--bamout` is to realign only those reads supporting variant haplotypes (see also `--full-bamout`).


### `--bamout-type`

Option `--bamout-type` is used to select the type of realigned evidence BAM to produce when using the `--bamout` option. If `FULL` is chosen then all reads in the input BAM are written to the evidence BAM. If `MINI` is chosen (default) then only reads supporting variant sites are written.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --bamout reads.realigned.bam --bamout-type FULL
```

### `--pedigree`

Option `--pedigree` is used to input a [pedigree file](http://www.helsinki.fi/~tsjuntun/autogscan/pedigreefile.html) (`.ped`). Some calling models may use pedigree information to improve calling accuracy.

```shell
$ octopus -R ref.fa -I mother.bam father.bam child.bam --pedigree family.ped 
```

where `family.ped` might contain, for example

```text
my_trio	FATHER	0	0	1	0
my_trio	MOTHER	0	0	2	0
my_trio	CHILD	FATHER	MOTHER	2	0
``` 

**Notes**

* Specifying this option with a `.ped` file defining a trio relationship between three input samples will automatically activate the `trio` calling model.

### `--fast`

Command `--fast` sets up Octopus for fast variant calling.

```shell
$ octopus -R ref.fa -I reads.bam --fast
```

**Warning** This command may reduce calling accuracy.

### `--very-fast`

Command `--very-fast` sets up Octopus for very fast variant calling.

```shell
$ octopus -R ref.fa -I reads.bam --very-fast

```

**Warning** This command may reduce calling accuracy

### `--data-profile`

Option `--data-profile` is used to generate a profile of the input read data, which may be used to make new [sequence error models](https://github.com/luntergroup/octopus/wiki/How-to:-Use-error-models).

```shell
$ octopus -R ref.fa -I reads.bam -o calls.bcf --data-profile reads.profile.csv

```

### `--bad-region-tolerance`

Option `--bad-region-tolerance` specifies the user tolerance for regions that may be 'uncallable' (e.g. due to mapping errors) and slow down calling. The possible arguments are:

* `LOW` Low tolerance to bad regions.
* `NORMAL` Default tolerance to bad regions.
* `HIGH` High tolerance to bad regions.
* `UNLIMITED` Turn off bad region detection.

```shell
$ octopus -R ref.fa -I reads.bam --bad-region-tolerance UNLIMITED

```

## Read pre-processing options

### `--disable-read-preprocessing`

Command `--disable-read-preprocessing` can be used to disable all optional read preprocessing - all viable raw input alignments will be used.

```shell
$ octopus -R ref.fa -I reads.bam --disable-read-preprocessing
```

### `--max-base-quality`

Option `--max-base-quality` caps the base quality of all input reads.

```shell
$ octopus -R ref.fa -I reads.bam --max-base-quality 40
```

### `--mask-low-quality-tails`

Option `--mask-low-quality-tails` is used to mask (assign base quality zero) read tail bases that have low base quality scores. The value provided to the option is the threshold used to define a 'low' quality score.

```shell
$ octopus -R ref.fa -I reads.bam --mask-low-quality-tails 5 // use 5 as minimum quality score
$ octopus -R ref.fa -I reads.bam --mask-low-quality-tails // use implicit 'low' score (see -h)
```

### `--mask-tails`

Option `--mask-tails` is used to unconditionally mask (assign base quality zero) read tail bases. The value provided to the option is the number of bases to mask.

```shell
$ octopus -R ref.fa -I reads.bam --mask-tails 3 // mask the last 3 bases of all reads
```

### `--mask-soft-clipped-bases`

Command `--mask-soft-clipped-bases` is used to mask (assign base quality zero) all read bases that are soft clipped in the input alignments.

```shell
$ octopus -R ref.fa -I reads.bam --mask-soft-clipped-bases
```

### `--soft-clip-mask-threshold`

Option `--soft-clip-mask-threshold` makes the option `--soft-clip-masking` conditional on the base quality score; only bases below the value given to the option will be masked.

```shell
$ octopus -R ref.fa -I reads.bam --soft-clip-masking yes --soft-clip-mask-threshold 10
```

### `--mask-soft-clipped-boundary-bases`

Option `--mask-soft-clipped-boundary-bases` makes the option `--soft-clip-masking` mask addition bases adjacent to soft clipped bases. The value provided to the option is the number of additional bases to mask.

```shell
$ octopus -R ref.fa -I reads.bam --soft-clip-masking yes --mask-soft-clipped-boundary-bases 5
```

### `--mask-inverted-soft-clipping`

Command `--mask-inverted-soft-clipping` is used to mask (assign base quality zero) all read bases that are soft clipped **and** are inverted copies of nearby sequence.

```shell
$ octopus -R ref.fa -I reads.bam --mask-inverted-soft-clipping
```

### `--mask-3prime-shifted-soft-clipped-heads`

Command `--mask-3prime-shifted-soft-clipped-heads` is used to mask (assign base quality zero) all read bases that are soft clipped **and** are shifted copies of nearby 3' sequence.

```shell
$ octopus -R ref.fa -I reads.bam --mask-3prime-shifted-soft-clipped-heads
```

### `--disable-adapter-masking`

Command `--disable-adapter-masking` is used to mask (assign base quality zero) read bases that are considered to be adapter contamination.

```shell
$ octopus -R ref.fa -I reads.bam --disable-adapter-masking
```

**Notes**

* The algorithm used to detect adapter contamination depends only on the input alignment mapping information; no library of adapter sequences are used.

### `--disable-overlap-masking`

Command `--disable-overlap-masking` is used to mask (assign base quality zero) read bases of read templates that contain overlapping segments, in order to remove non-independent base observations. 

```shell
$ octopus -R ref.fa -I reads.bam --disable-overlap-masking
```

**Notes**

* All but of the overlapping read bases are masked, leaving one base untouched. If two segments are overlapping, then half of the 5' bases of each segment are masked.

### `--split-long-reads`

Command `--split-long-reads` stipulates that reads longer than [`--max-read-length`](#option---max-read-length) should be split into smaller linked reads.

```shell
$ octopus -R ref.fa -I reads.bam --max-read-length 200 --split-long-reads
```

### `--consider-unmapped-reads`

Command `--consider-unmapped-reads` turns off the read filter that removes reads marked as unmapped in the input alignments.

```shell
$ octopus -R ref.fa -I reads.bam --consider-unmapped-reads
```

### `--min-mapping-quality`

Option `--min-mapping-quality` specifies the minimum mapping quality that reads must have to be considered; reads with mapping quality below this will be filtered and not used for analysis.

```shell
$ octopus -R ref.fa -I reads.bam --min-mapping-quality 10
```

### `--good-base-quality`

Option `--good-base-quality` defines the minimum quality of a 'good' base for the options `--min-good-base-fraction` and `--min-good-bases`.

```shell
$ octopus -R ref.fa -I reads.bam --good-base-quality 10
```

### `--min-good-base-fraction`

Option `--min-good-base-fraction` specifies the fraction of 'good' (see `--good-base-quality`) base qualities a read must have in order to be considered; reads with a fraction of 'good' base qualities less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-good-base-fraction 0.3
```

### `--min-good-bases`

Option `--min-good-bases` specifies the number of 'good' (see `--good-base-quality`) base qualities a read must have in order to be considered; reads with 'good' base qualities less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-good-bases 50
```

### `--allow-qc-fails`

Command `--allow-qc-fails` permits reads marked as QC fail in the input alignments to be considered, otherwise they will be filtered.

```shell
$ octopus -R ref.fa -I reads.bam --allow-qc-fails
```

### `--min-read-length`

Option `--min-read-length` specifies the minimum length (number of sequence bases) of reads to be considered; reads with length less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-read-length 50
```

### `--max-read-length`

Option `--max-read-length` specifies the maximum length (number of sequence bases) of reads to be considered; reads with length greater than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --max-read-length 500
```

### `--allow-marked-duplicates`

Command `--allow-marked-duplicates` disables the filter that removes reads marked as duplicate in the input alignments. The default behaviour is to remove reads marked duplicate.

```shell
$ octopus -R ref.fa -I reads.bam --allow-marked-duplicates
```

### `--allow-octopus-duplicates`

Command `--allow-octopus-duplicates` disables the filter that removes reads that Octopus considers to be duplicates. The default behaviour is to remove duplicate reads.

```shell
$ octopus -R ref.fa -I reads.bam --allow-octopus-duplicates
```

### `--duplicate-read-detection-policy`

Option `--duplicate-read-detection-policy` specifies approach to use for detecting [duplicate reads](https://github.com/luntergroup/octopus/wiki/Frequently-Asked-Questions#should-i-deduplicate-my-input-bams). Possible arguments are:

* `RELAXED` Require 5' mapping co-ordinate matches **and** identical cigar strings.
* `AGGRESSIVE` Require 5' mapping co-ordinate matches only.

```shell
$ octopus -R ref.fa -I reads.bam --duplicate-read-detection-policy AGGRESSIVE
```

### `--allow-secondary-alignments`

Command `--allow-secondary-alignments` disables the filter that removes reads that are marked as secondary alignments in the input alignments. The default behaviour is to remove reads marked secondary.

```shell
$ octopus -R ref.fa -I reads.bam --allow-secondary-alignments
```

### `--allow-supplementary-alignments`

Command `--allow-supplementary-alignments` disables the filter that removes reads that are marked as supplementary alignments in the input alignments. The default behaviour is to remove reads marked supplementary.

```shell
$ octopus -R ref.fa -I reads.bam --allow-supplementary-alignments
```

### `--max-decoy-supplementary-alignment-mapping-quality`

Option `--max-decoy-supplementary-alignment-mapping-quality` removes any reads with supplementary alignments (i.e. `SA` SAM tag) to decoy contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-decoy-supplementary-alignment-mapping-quality 20
```

### `--max-unplaced-supplementary-alignment-mapping-quality`

Command `--max-unplaced-supplementary-alignment-mapping-quality removes reads with supplementary alignments (i.e. `SA` SAM tag) to unplaced contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-unplaced-supplementary-alignment-mapping-quality 20
```

### `--max-unlocalized-supplementary-alignment-mapping-quality`

Command `--max-unlocalized-supplementary-alignment-mapping-quality` removes reads with supplementary alignments (i.e. `SA` SAM tag) to unlocalized contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-unlocalized-supplementary-alignment-mapping-quality 20
```

### `--no-reads-with-unmapped-segments`

Command `--no-reads-with-unmapped-segments` removes reads that are marked as having unmapped segments in the input alignments. The default behaviour is to allow reads with unmapped segments.

```shell
$ octopus -R ref.fa -I reads.bam --no-reads-with-unmapped-segments
```

### `--no-reads-with-distant-segments`

Command `--no-reads-with-distant-segments` removes reads that have segments mapped to different contigs. The default behaviour is to allow reads with distant segments.

```shell
$ octopus -R ref.fa -I reads.bam --no-reads-with-distant-segments
```

### `--no-adapter-contaminated-reads`

Command `--no-adapter-contaminated-reads` removes reads that are considered to have adapter contamination (i.e. sequence bases from the adapter). The default behaviour is to allow reads with adapter contamination.

```shell
$ octopus -R ref.fa -I reads.bam --no-adapter-contaminated-reads
```

### `--disable-downsampling`

Command `--disable-downsampling` turns off downsampling.

```shell
$ octopus -R ref.fa -I reads.bam --disable-downsampling
```

### `--downsample-above`

Option `--downsample-above` specifies the read depth required to mark a position as a candidate for the downsampler.

```shell
$ octopus -R ref.fa -I reads.bam --downsample-above 5000 # depths up to 5000x are permitted
```

### `--downsample-target`

Option `--downsample-target` specifies the target read depth for the downsampler for all candidate sites. Reads will be removed from the input alignments until all positions have read depth not greater than this.

```shell
$ octopus -R ref.fa -I reads.bam --downsample-target 100 # downsample to 100x
```

### `--use-same-read-profile-for-all-samples`

Command `--use-same-read-profile-for-all-samples` specifies that the same input read profile should be used for all samples, rather than generating one for each sample. This essentially means that the same read distribution is assumed for all samples. 

```shell
$ octopus -R ref.fa -I reads.bam --use-same-read-profile-for-all-samples
```

## Variant discovery options

### `--variant-discovery-mode`

Option `--variant-discovery-mode` specifies the thresholds used for candidate variant discovery, which affects the overall sensitivity of the generators. Possible values are:

* `ILLUMINA` The default mode, for Illumina quality reads.
* `PACBIO` For PacBio quality reads - requires more observations to propose a candidate, particular indel candidates.

```shell
$ octopus -R ref.fa -I reads.bam --variant-discovery-mode PACBIO
```

### `--disable-denovo-variant-discovery`

Command `--disable-denovo-variant-discovery` disables the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-denovo-variant-discovery
```

### `--disable-repeat-candidate-generator`

Command `--disable-repeat-candidate-generator` disables the tandem repeat candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-repeat-candidate-generator
```

### `--disable-assembly-candidate-generator`

Command `--disable-assembly-candidate-generator` disables the local *de novo* assembly candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-assembly-candidate-generator
```

### `--source-candidates`

Option `--source-candidates` (short `-c`) accepts a list of VCF files, the contents of which will be added to the candidate variant set.

```shell
$ octopus -R ref.fa -I reads.bam --source-candidates calls.vcf.gz
```

**Notes**

* The option accepts files in the `.vcf`, `.vcf.gz`, and `.bcf` formats. Files in the `.vcf.gz` and `.bcf` format must be indexed. For larger

### `--source-candidates-file`

Option `--source-candidates-file` can be used to provide one or more files containing lists (one per line) of files in the `.vcf`, `.vcf.gz`, and `.bcf` formats (see notes on `--source-candidates`).

```shell
$ octopus -R ref.fa -I reads.bam --source-candidates-file vcfs.txt
```

where `vcfs.txt` may contain, for example

```text
/path/to/vcfs/callsA.vcf.gz
/path/to/vcfs/callsB.bcf
```

### `--min-source-candidate-quality`

Option `--min-source-candidate-quality` specifies the minimum `QUAL` score for a user provided source variant (see `--source-candidates` and `--source-candidates-file`) to be added to the final candidate variant list; any records with `QUAL` less than this will not be considered.

```shell
$ octopus -R ref.fa -I reads.bam -c calls.vcf.gz --min-source-candidate-quality 10
```

### `--use-filtered-source-candidates`

Command `--use-filtered-source-candidates` specifies allows variants in the user-provided candidate variants (see `--source-candidates` and `--source-candidates-file`) that are marked as filtered (according to the `FILTER` column) should be added to the final candidate variant list. The default behaviour is to remove filtered variants.

```shell
$ octopus -R ref.fa -I reads.bam -c calls.vcf.gz --use-filtered-source-candidates
```

### `--min-pileup-base-quality`

Option `--min-pileup-base-quality` specifies the minimum base quality that a SNV in the input alignments must have in ordered to be considered by the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --min-pileup-base-quality 10
```

### `--min-supporting-reads`

Option `--min-supporting-reads` specifies the minimum number of supporting reads a variant must have in the input alignments to be included in the candidate variant list from the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --min-supporting-reads 3
```

### `--allow-pileup-candidates-from-likely-misaligned-reads`

Command `--allow-pileup-candidates-from-likely-misaligned-reads` stops the pileup candidate generator filtering out candidates that are considered likely to originate from read misalignment, which can otherwise result in many false candidates.

```shell
$ octopus -R ref.fa -I reads.bam --allow-pileup-candidates-from-likely-misaligned-reads
```

### `--max-variant-size`

Option `--max-variant-size` specifies the maximum size (w.r.t reference region) a candidate variant can have to be considered by the calling algorithm; candidate variants with size greater than this will be removed and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --max-variant-size 500
```

### `--kmer-sizes`

Option `--kmer-sizes` specifies the default kmer sizes to try for local *de novo* assembly. Assembly graphs will be constructed for each kmer size in all active assembler regions, and the union of candidate variants from each assembler graph used for the final candidate set.

```shell
$ octopus -R ref.fa -I reads.bam --kmer-sizes 5 10 15 20 25 30 35 40 45 50
```

### `--num-fallback-kmers`

Option `--num-fallback-kmers` specifies the number of additional kmer sizes to try in the case that none of the default (i.e. kmer sizes specified in `--kmer-sizes`) is able to construct a valid assembly graph (often the case for smaller kmer sizes).

```shell
$ octopus -R ref.fa -I reads.bam --num-fallback-kmers 20
```

### `--fallback-kmer-gap`

Option `--fallback-kmer-gap` speifies the increment between fallback kmer sizes (see `--num-fallback-kmers`)

```shell
$ octopus -R ref.fa -I reads.bam --fallback-kmer-gap 5
$ octopus -R ref.fa -I reads.bam --kmer-sizes 5 10 --num-fallback-kmers 2 --fallback-kmer-gap 5 # always try kmer sizes 5 and 10, and then try 15 and 20 if 5 and 10 don't work
```

### `--max-region-to-assemble`

Option `--max-region-to-assemble` specifies the maximum reference region size that can be used to construct a local *de novo* assembly graph. Larger values enable detection of larger variants (e.g. large deletions), but increase the likelihood that small kmer sizes will result in an invalid assembly graph, and therefore decrease sensitivity for smaller variation.

```shell
$ octopus -R ref.fa -I reads.bam --max-region-to-assemble 1000
```

### `--max-assemble-region-overlap`

Option `--max-assemble-region-overlap` specifies the maximum overlap between reference regions used to build local *de novo* assembly graphs. Larger overlaps result in more assembly graphs being constructed and may increase sensitivity for variation, at the expense of compute time. 

```shell
$ octopus -R ref.fa -I reads.bam --max-assemble-region-overlap 100
```

### `--assemble-all`

Command `--assemble-all` forces all reference regions to be used for local *de novo* assembly, rather than only regions considered likely to contain variation.

```shell
$ octopus -R ref.fa -I reads.bam --assemble-all
```

**Warning** This command may result in substantially longer runtimes. 

### `--assembler-mask-base-quality`

Option `--assembler-mask-base-quality` specifies the minimum base quality an aligned base should have to avoid being 'masked' (i.e. converted to reference) before being inserted into the local *de novo* assembly graph. Higher values reduce sensitivity to noise and increase the likelihood of finding bubbles in graphs with larger kmer sizes, as the cost of decreased sensitivity at lower kmer sizes. 

```shell
$ octopus -R ref.fa -I reads.bam --assembler-mask-base-quality 20
```

### `--allow-cycles`

Command `--allow-cycles` forces the assembler to consider assembly graphs containing non-reference cycles, which usually result in false candidates.

```shell
$ octopus -R ref.fa -I reads.bam --allow-cycles
```

### `--min-kmer-prune`

Option `--min-kmer-prune` specifies the minimum number of kmer observations that must be present in the local *de novo* assembly graph for the kmer to stay in the graph before variant extraction. Lower values increase sensitivity but lower specificity.

```shell
$ octopus -R ref.fa -I reads.bam --min-kmer-prune 3
```

### `--max-bubbles`

Option `--max-bubbles` specifies the maximum number of bubbles in the final local *de novo* assembly graph that may be explored for candidate variant generation. Higher values increase sensitivity but reduce specificity.

```shell
$ octopus -R ref.fa -I reads.bam --max-bubbles 100
```

### `--min-bubble-score`

Option `--min-bubble-score` specifies the minimum score a bubble explored in a local *de novo* assembly graph must have for the bubble to be extracted as a candidate variant. The 'bubble score' is the average of all kmer observations along the bubble edge, scaled by the probability the kmer observations along the bubble edge originate from a single read strand. The bubble score can be viewed as a proxy for the number of 'good' read observations for the bubble (i.e. variant), and therefore higher values increase specificity but reduce sensitivity.

```shell
$ octopus -R ref.fa -I reads.bam --min-bubble-score 5
```

### `--min-candidate-credible-vaf-probability`

Option `--min-candidate-credible-vaf-probability` sets the minimum probability mass above `--min-credible-somatic-frequency` required to 'discover' a variant when using the `cancer` calling model. Smaller values increase the number of candidate variants generated, potentially improving sensitivity but also increasing computational complexity.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-candidate-credible-vaf-probability 0.5
```

## Haplotype generation options

### `--max-haplotypes`

Option `--max-haplotypes` specifies the maximum number of haplotypes that are considered by the calling model. It also specifies the target number of haplotypes that the haplotype generator should produce on each iteration of the algorithm. If the haplotype generator is unable to satisfy the request (i.e. produces a greater number of haplotypes), then the set of haplotypes is reduced to this size by removing haplotypes considered unlikely using several likelihood based statistics.

Increasing `--max-haplotypes` reduces the chance that a true haplotype is incorrectly filtered before evaluation by the calling model. It also increases the number of candidate variants that may be considered on each iteration of the algorithm, potentially improving calling accuracy. However, increasing this value will usually increase runtimes - sometimes substantially. 

```shell
$ octopus -R ref.fa -I reads.bam --max-haplotypes 500
```

### `--haplotype-holdout-threshold`

Option `--haplotype-holdout-threshold` specifies the number of haplotypes that the haplotype generator can produce before some active candidate variants are added to the holdout stack. The value must not be less than `--max-haplotypes`.

```shell
$ octopus -R ref.fa -I reads.bam --haplotype-holdout-threshold 1000
```

### `--haplotype-overflow`

Option `--haplotype-overflow` specifies the number of haplotypes that the haplotype generator can produce before the current active region must be skipped. The value must not be less than `--haplotype-holdout-threshold`.

```shell
$ octopus -R ref.fa -I reads.bam --haplotype-overflow 5000
```

### `--max-holdout-depth`

Option `--max-holdout-depth` specifies the maximum size of the holdout stack.

```shell
$ octopus -R ref.fa -I reads.bam --max-holdout-depth 0 # no holdouts
```

### `--extension-level`

Option `--extension-level` specifies the condition for extending the active haplotype tree with novel alleles. The possible values are:

* `MINIMAL` Only include novel alleles overlapping with overlapping reads. 
* `NORMAL` Include novel alleles with reads overlapping the rightmost included allele.
* `AGGRESSIVE` No conditions on extension other than the number of alleles.

More aggressive extension levels result in larger haplotype blocks, which may improve phase lengths and accuracy. However, aggressive extension increases the possibility of including novel alleles that cannot be phased with active alleles which will increase compute time.

```shell
$ octopus -R ref.fa -I reads.bam --extension-level AGGRESSIVE
```

### `--lagging-level`

Option `--lagging-level` specifies the extent to which active alleles remain active in the next algorithm iteration, increasing the length of haplotype blocks.

* `NONE` Disable lagging; each iteration of the algorithm will evaluate a novel set of alleles.
* `NORMAL` Consider previous active alleles if there are reads overlapping them with the next active allele set.
* `AGGRESSIVE` Consider previous active alleles if there are overlapping reads that span them, and the next active alleles.

```shell
$ octopus -R ref.fa -I reads.bam --lagging-level AGGRESSIVE
```

### `--backtrack-level`

Option `--backtrack-level` specifies the extent to which 

* `NONE` Disable all backtracking.
* `NORMAL` Enables backtracking.
* `AGGRESSIVE` Currently the same as `NORMAL`.

```shell
$ octopus -R ref.fa -I reads.bam --backtrack-level NORMAL
```

### `--min-protected-haplotype-posterior`

Option `--min-protected-haplotype-posterior` specifies the minimum posterior probability that a haplotype is present in the samples (according to the calling model) for the haplotype to avoid being removed from consideration; haplotypes with posterior probability less than this may be filtered. Increasing the value of this option results in a greater number of haplotypes being filtered, allowing the haplotype tree to grow to include more candidate alleles, and reducing computational complexity. However, larger values also increase the chance that a true haplotype is incorrectly discarded. 

```shell
$ octopus -R ref.fa -I reads.bam --min-protected-haplotype-posterior 1e-5
```

### `--dont-protect-reference-haplotype`

Command `--dont-protect-reference-haplotype` disables protection of the reference haplotype during haplotype filtering, and therefore ensures that the reference haplotype is always considered by the calling model.

```shell
$ octopus -R ref.fa -I reads.bam --dont-protect-reference-haplotype
```

### `--bad-region-tolerance`

Option `--bad-region-tolerance` specifies the 'tolerance' for spending time calling variants in regions that are unlikely to be callable (e.g. due to low complexity sequence). Such regions tend to be computationally difficult and skipping them can save a lot of wasted computation time. However, identification of such regions is not perfect and could result in skipping regions with real variation. Possible argument values are:

* `LOW` Skip region that show any signs of being uncallable.
* `NORMAL` Skip region that show reasonable signs of being uncallable.
* `HIGH` Skip region that show strong signs of being uncallable.

```shell
$ octopus -R ref.fa -I reads.bam --bad-region-tolerance LOW
```

## Common variant calling options

### `--caller`

Option `--caller` (short '-C') specifies the [calling model](https://github.com/luntergroup/octopus/wiki/Calling-models) to be used. The option must only be set if the calling model is not automatically determined from other options.

```shell
$ octopus -R ref.fa -I reads.bam --caller cancer # e.g. for tumour-only
```

### `--organism-ploidy`

Option `--organism-ploidy` (short '-P') specifies the default ploidy of all input samples. All contigs will be assumed to have this ploidy unless specified otherwise in `--contig-ploidies`.

```shell
$ octopus -R ref.fa -I reads.bam --organism-ploidy 3 # triploid calling
```

### `--contig-ploidies`

Option `--contig-ploidies` (short `-p`) can be used to specify the ploidy of contigs (format `contig=ploidy`) that are not the same as the `--organism-ploidy`, or the ploidy of individual sample contigs (format `sample:contig=ploidy`).

```shell
$ octopus -R ref.fa -I reads.bam --contig-ploidies X=1 # Contig "X" has ploidy 1 for all samples
$ octopus -R ref.fa -I reads.bam --contig-ploidies HG002:X=1 # Contig "X" has ploidy 1 for sample "HG002"
```

### `--contig-ploidies-file`

Option `--contig-ploidies-file` can be used to provide a file specifying contig or sample contig ploidies (see `--contig-ploidies`), one per line.  

```shell
$ octopus -R ref.fa -I reads.bam --contig-ploidies-file ploidies.txt
```

where `ploidies.txt` could contain

```text
X=1
HG002:X=1
```

### `--min-variant-posterior`

Option `--min-variant-posterior` specifies the minimum posterior probability (Phred scale) required to call a variant. Candidate variants with posterior probability less than this will not be reported in the final call set.

```shell
$ octopus -R ref.fa -I reads.bam --min-variant-posterior 0.5
```

**Notes**

* For calling models with more than one class of variation, this option refers to the calling models default variant class. For example, in the `cancer` calling model there are both germline and somatic variants, and this option refers to germline variants (see the option `--min-somatic-posterior` for somatic variants).

### `--refcall`

Option `--refcall` can be used to enable reference calling, meaning Octopus will generate a 'gVCF' format output. There are two possible arguments for this option:

* `BLOCKED` Merge adjacent called reference positions with similar quality (see `--refcall-block-merge-quality`) a single gVCF record.
* `POSITIONAL` Emit a gVCF record for all positions.

If the option is specified but no argument is provide, `BLOCKED` is assumed.

```shell
$ octopus -R ref.fa -I reads.bam --refcall BLOCKED
$ octopus -R ref.fa -I reads.bam --refcall # same as above
$ octopus -R ref.fa -I reads.bam --refcall POSITIONAL
```

### `--refcall-block-merge-quality`

Option `--refcall-block-merge-quality` specifies the quality (Phred scale) threshold for merging adjacent called reference positions when `BLOCKED` refcalls are requested; adjacent reference positions with an absolute quality difference less or equal than this will be merged into a block.

```shell
$ octopus -R ref.fa -I reads.bam --refcall --refcall-block-merge-quality 20
```

### `--min-refcall-posterior`

Option `--min-refcall-posterior` specifies the minimum posterior probability required to call a position as homozygous reference; positions with posterior probability less than this will not be reported in the output gVCF.

```shell
$ octopus -R ref.fa -I reads.bam --min-refcall-posterior 10
```

### `--max-refcall-posterior`

Option `--max-refcall-posterior` caps the `QUAL` of all reference calls, which may result in larger reference blocks and smaller gVCF file sizes.

```shell
$ octopus -R ref.fa -I reads.bam --max-refcall-posterior 100
```

### `--snp-heterozygosity`

Option `--snp-heterozygosity` specifies the SNV heterozygosity parameter for the Coalesent mutation model, used to assign prior probabilities.

```shell
$ octopus -R ref.fa -I reads.bam --snp-heterozygosity 0.1
```

### `--snp-heterozygosity-stdev`

Option `--snp-heterozygosity-stdev` specifies the standard deviation of the SNV heterozygosity parameter (see `--snp-heterozygosity`).

```shell
$ octopus -R ref.fa -I reads.bam --snp-heterozygosity-stdev 0.1
```

### `--indel-heterozygosity`

Option `--indel-heterozygosity` specifies the INDEL heterozygosity parameter for the Coalesent mutation model, used to assign prior probabilities.

```shell
$ octopus -R ref.fa -I reads.bam --indel-heterozygosity 0.001
```

### `--max-genotypes`

Option `--max-genotypes` specifies the maximum number of candidate genotypes that must be evaluated by the calling model. If there are more possible candidate genotypes than this value then the algorithm may decide to remove some candidate genotypes using heuristics. The number of candidate genotypes to evaluate can have a substantial impact on runtime for some calling models (e.g. `cancer`).

```shell
$ octopus -R ref.fa -I reads.bam --max-genotypes 1000
```

### `--max-genotype-combinations`

Option `--max-genotype-combinations` specifies the maximum number of candidate joint genotype vectors that must be evaluated by the calling model. If there are more possible candidate joint genotype vectors than this value the algorithm may decide to remove some candidate genotype vectors using heuristics. The number of candidate joint genotype vectors to evaluate can have a substantial impact on runtime. 

```shell
$ octopus -R ref.fa -I reads.bam --max-genotype-combinations 100000
```

**Notes**

* This option is only used by calling models that consider joint genotypes (i.e. `population` and `trio`)

### `--use-uniform-genotype-priors`

Command `--use-uniform-genotype-priors` indicates that the uniform genotype prior model should be used for calling genotypes and variants.

```shell
$ octopus -R ref.fa -I reads.bam --use-uniform-genotype-priors
```

### `--use-independent-genotype-priors`

Command `--use-independent-genotype-priors` indicates that an independent genotype prior model should be used for evaluating joint genotype vectors.

```shell
$ octopus -R ref.fa -I reads.bam --use-independent-genotype-priors
```

### `--model-posterior`

Option `--model-posterior` enables or disables model posterior evaluation at each called variant site.

```shell
$ octopus -R ref.fa -I reads.bam --model-posterior yes
```

### `--disable-inactive-flank-scoring`

Command `--disable-inactive-flank-scoring` disables an additional step during haplotype likelihood calculation that attempts to correct low likelihoods caused by inactive variation in the flanking regions of the haplotype under evaluation.

```shell
$ octopus -R ref.fa -I reads.bam --disable-inactive-flank-scoring
```

### `--dont-model-mapping-quality`

Command `--dont-model-mapping-quality` disables consideration of mapping quality in the haplotype likelihood calculation. This can improve calling accuracy if read mapping qualities are well calibrated.

```shell
$ octopus -R ref.fa -I reads.bam --dont-model-mapping-quality
```

### `--sequence-error-model`

Option `--sequence-error-model` specifies the [sequence error model](https://github.com/luntergroup/octopus/wiki/How-to:-Use-error-models) to use.

```shell
$ octopus -R ref.fa -I reads.bam --sequence-error-model PCR.NOVASEQ # built-in error model
$ octopus -R ref.fa -I reads.bam --sequence-error-model /path/to/my/error.model # custom error model
```

**Notes**

* The same error model is used for all input reads. 

### `--max-indel-errors`

Option `--max-indel-errors` specifies the maximum number of indel errors in an individual read fragment that can be accuretely modelled by the haplotype likelihood model. Larger values usually require greater computational resources (determined by your systems available SIMD instructions).

```shell
$ octopus -R ref.fa -I reads.bam --max-indel-errors 32
```

### `--use-wide-hmm-scores`

Command `--use-wide-hmm-scores` sets the score variable computed by the pair HMM for haplotype likelihoods to 32 bits (the default is 16 bits). This can avoid score overflow in long noisy reads, but will slow does the computation.

```shell
$ octopus -R ref.fa -I reads.bam --use-wide-hmm-scores
```

### `--max-vb-seeds`

Option `--max-vb-seeds` specifies the maximum number of seeds that Variational Bayes models can use for posterior evaluation. Increasing the number of seeds increases the likelihood that a posterior mode will be identified, but results in more computation time.

```shell
$ octopus -R ref.fa -I reads.bam --max-vb-seeds 50
```

### `--read-linkage`

Option `--read-linkage` specifies how reads are linked in the input alignments. Read linkage information is used by the haplotype likelihood calculation and can improve calling accuracy and increase phase lengths.

* `NONE` Reads are not linked in any way.
* `PAIRED` Reads may be paired, with pairs having identical read names.
* `LINKED` Reads may be linked or paired, with linked reads having identical `BX` tags.

```shell
$ octopus -R ref.fa -I reads.bam --read-linkage LINKED
```

### `--min-phase-score`

Option `--min-phase-score` specifies the minimum phase score (`PQ` in VCF) required to emit adjacent variant calls in the same [phase set](https://github.com/luntergroup/octopus/wiki/VCF-format#phasing). Increasing this value results in less sites being phased, but reduces the phase false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --min-phase-score 5
```

### `--disable-early-phase-detection`

Command `--disable-early-phase-detection` prevents the phasing algorithm being applied to partially resolved haplotype blocks, which can lead to removal of complete phased segments from the head of the current haplotype block. This heuristic can prevent discontiguous phase blocks being resolved, which are more likely in some data (e.g. linked reads).

```shell
$ octopus -R ref.fa -I reads.bam --disable-early-phase-detection
```

## Cancer variant calling options

### `--normal-samples`

Option `--normal-samples` specifies which of the input samples are normal samples for tumour-normal paired analysis.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL
```

**Notes**

* Specifying this option will automatically activate the `cancer` calling model.

### `--max-somatic-haplotypes`

Option `--max-somatic-haplotypes` specifies the maximum number of unique haplotypes containing somatic variation that can be modelled by the somatic genotype model. If there are more true somatic haplotypes present in the input data than this value then the model will not accuretely fit the data and true somatic variants may not be called, however, larger values substantially increase the computational complexity of the model and potentially increase the false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --max-somatic-haplotypes 3
```

### `--somatic-snv-prior`

Option `--somatic-snv-prior` specifies the somatic SNV mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --somatic-snv-prior 1e-7
```

### `--somatic-indel-prior`

Option `--somatic-indel-prior` specifies the somatic INDEL mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --somatic-indel-prior 1e-7
```

### `--min-expected-somatic-frequency`

Option `--min-expected-somatic-frequency` specifies the minimum expected Variant Allele Frequency (VAF) for somatic mutations in the samples. This value is used by the `cancer` calling model as the lower bound for the VAF posterior marginalisation used to compute the posterior probability that a candidate variant is a somatic mutation. Decreasing this value increases sensitivity for somatic mutations with smaller VAFs, but also increases sensitivity to noise.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-expected-somatic-frequency 0.05
```

### `--min-credible-somatic-frequency`

Option `--min-credible-somatic-frequency` specifies the minimum credible Variant Allele Frequency (VAF) for somatic mutations in the samples. This value is used for candidate variant discovery, and also by the `cancer` calling model as the lower-bound on candidate somatic mutation VAF credible regions (i.e. if the credible region computed for a mutation contains VAFs less than this value then the mutation will not be called as somatic). Decreasing this option increases sensitivity for somatic mutations with lower VAFs, but also increases sensitivity to noise and increases computational complexity as more candidate variants will be generated.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-credible-somatic-frequency 0.001
```

**Notes**

* If the value provided to this option is greater than the value specified by `--min-expected-somatic-frequency`, then this value is used for both options.

### `--tumour-germline-concentration`

Option `--tumour-germline-concentration` sets the Dirichlet concentration parameter for the germline haplotypes of tumour samples. Larger values concentrate more prior probability mass on equal frequencies of germline haplotypes and also increase prior mass on low VAFs for somatic haplotypes. This can help the model correctly classify somatic variation when the normal sample is not informative (or not present), but also reduces sensitivity to somatic variation with larger VAFs.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --tumour-germline-concentration 5
```

### `--somatic-credible-mass`

Option `--somatic-credible-mass` specifies the probability mass to used to compute the credible interval for determining whether to call a variant as somatic (see also `--min-credible-somatic-frequency`). Larger values result in wider credible regions, increasing sensitivity for lower VAF somatic mutations, but also increase sensitivity to noise. 

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatic-credible-mass 0.99
```

### `--min-somatic-posterior`

Option `--min-somatic-posterior` specifies the minimum posterior probability (Phred scale) required to call a candidate variant as a somatic mutation.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-somatic-posterior 1
```

### `--normal-contamination-risk`

Option `--normal-contamination-risk` indicates the risk that the normal sample contains contamination from the tumour samples. There are two possible values:

* `LOW` The algorithm will not consider normal contamination when generating candidate genotypes.
* `HIGH` The algorithm will consider normal contamination when generating candidate genotypes.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --normal-contamination-risk HIGH
```

### `--somatics-only`

Command `--somatics-only` indicates that only variant sites called as somatic mutations (i.e. tagged `SOMATIC`) should appear in the final output.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatics-only
```

**Warning** Using this command will produce a VCF file that cannot be re-filtered using Octopus.

## Trio variant calling options

### `--maternal-sample`

Option `--maternal-sample` (short '-M') indicates which of the input samples is the mother of the proband. If this option is specified then `--paternal-sample` must also be specified.

```shell
$ octopus -R ref.fa -I reads.bam -M mother -F father
```

**Notes**

* Specifying this option will automatically activate the `trio` calling model.

### `--paternal-sample`

Option `--paternal-sample` (short '-F') indicates which of the input samples is the father of the proband. If this option is specified then `--maternal-sample` must also be specified.

```shell
$ octopus -R ref.fa -I reads.bam -M mother -F father
```

**Notes**

* Specifying this option will automatically activate the `trio` calling model.

### `--denovo-snv-prior`

Option `--denovo-snv-prior` specifies the *de novo* SNV mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-snv-prior 1e-7
```

### `--denovo-indel-prior`

Option `--denovo-indel-prior` specifies the *de novo* INDEL mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-indel-prior 1e-8
```

### `--min-denovo-posterior`

Option `--min-denovo-posterior` specifies the minimum posterior probability (Phred scale) required to call a candidate variant as a *de novo* mutation.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --min-denovo-posterior 2
```

### `--denovos-only`

Command `--denovos-only` indicates that only variant sites called as *de novo* mutations (i.e. tagged `DENOVO`) should appear in the final output.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovos-only
```

**Warning** Using this command will produce a VCF file that cannot be re-filtered using Octopus.

## Polyclone variant calling options

### `--max-clones`

Option `--max-clones` specifies the maximum number of clones that can be modelled by the `polyclone` calling model. If there are more clones than this present in the input data then the model will not fit the data well and some true variation may not be called. However, larger values increase the computational complexity of the model and also increase sensitivity to noise.

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --max-clones 5
```

### `--min-clone-frequency`

Option `--min-clone-frequency` specifies the lower-bound Variant Allele Frequency (VAF) to use when computing the posterior probability for a variant. Smaller values increase sensitivity and the false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --min-clone-frequency 0.05
```

### `--clone-prior`

Option `--clone-prior` sets the prior probability for each new clone proposed by the model

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --clone-prior 0.1
```

### `--clone-concentration`

Option `--clone-concentration ` sets the concentration parameter for the symmetric Dirichlet distribution used to model clone frequencies

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --clone-concentration 10
```

## Cell variant calling options

### `--max-copy-loss`

Option `--max-copy-loss` specifies the maximum number of germline haplotypes losses that can be considered by the model.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --max-copy-loss 1
```

### `--max-copy-gains`

Option `--max-copy-gains ` specifies the maximum number of haplotypes gains that can be considered by the model.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --max-copy-gains 1
```

### `--somatic-cnv-prior`

Option `--somatic-cnv-prior ` sets the prior probability of a loci having a copy-number change.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --somatic-cnv-prior 1e-10
```

### `--dropout-concentration`

Option `--dropout-concentration` sets the default Dirichlet concentration prior on haplotype frequencies. The higher the concentration parameter, the more probability mass around the centre of the distribution (0.5 for diploid). Lowering the concentration parameter means there’s more mass on extreme frequencies, which in turn means the model is less sensitive to dropout, but also less sensitive to real somatic variation. Setting to unity would put a uniform prior on frequencies.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --dropout-concentration 10
```

### `--sample-dropout-concentration`

Option `--sample-dropout-concentration` sets the Dirichlet concentration prior on haplotype frequencies for a specific sample, otherwise `--dropout-concentration` is used.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --sample-dropout-concentration NORMAL=100
```

### `--phylogeny-concentration`

Option `--phylogeny-concentration` sets the symmetric Dirichlet concentration prior on group mixture proportions in each phylogeny. A larger concentration parameter implies a more even distribution of samples across the tree.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --phylogeny-concentration 10
```

## Variant filtering options

### `--disable-call-filtering`

Command `--disable-call-filtering` disables variant call filtering.

```shell
$ octopus -R ref.fa -I reads.bam --disable-call-filtering
```

### `--filter-expression`

Option `--filter-expression` sets the threshold filter expression to use for filtering variants not tagged with `SOMATIC` or `DENOVO`.

```shell
$ octopus -R ref.fa -I reads.bam --filter-expression "QUAL < 10" # PASS calls with QUAL >= 10
```

### `--somatic-filter-expression`

Option `--somatic-filter-expression` sets the threshold filter expression to use for filtering `SOMATIC` variants.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatic-filter-expression "QUAL < 20 | PP < 20"
```

### `--denovo-filter-expression`

Option `--denovo-filter-expression` sets the threshold filter expression to use for filtering `DENOVO` variants.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-filter-expression "QUAL < 20 | PP < 20"
```

### `--refcall-filter-expression`

Option `--refcall-filter-expression` sets the threshold filter expression to use for filtering sites called homozygous reference.

```shell
$ octopus -R ref.fa -I reads.bam --refcall --refcall-filter-expression "QUAL < 10"
```

### `--use-preprocessed-reads-for-filtering`

Command `--use-preprocessed-reads-for-filtering` forces use of the same read pre-processing steps used for calling variants for filtering variants; otherwise all well-formed reads are used for filtering.

```shell
$ octopus -R ref.fa -I reads.bam --use-preprocessed-reads-for-filtering
```

### `--keep-unfiltered-calls`

Command `--keep-unfiltered-calls` requests that Octopus keep a copy of the VCF file produced before filtering is applied. The copy has `unfiltered` appended to the final output name.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --keep-unfiltered-calls # writes calls.unfiltered.vcf
```

### `--annotations`

Option `--annotations` requests that the values of a sub-set of the measures used for filtering are reported in the final VCF output.

```shell
$ octopus -R ref.fa -I reads.bam --annotations SB # adds SB FORMAT field to each record
```

### `--filter-vcf`

Option `--filter-vcf` specifies an Octopus VCF file to filter. No calling is performed.

```shell
$ octopus -R ref.fa -I reads.bam --filter-vcf calls.bcf
```

### `--forest-model`

Option `--forest-model` enables [random forest variant filtering](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-Random-Forest). The argument to the option is a ranger forest file.

```shell
$ octopus -R ref.fa -I reads.bam \
	--forest-model octopus/resources/forests/germline.v0.6.3-beta.forest
```

### `--somatic-forest-model`

Option `--somatic-forest-model` enables [random forest variant filtering](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-Random-Forest) for somatic variants. The argument to the option is a ranger forest file.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL \
	--forest-model octopus/resources/forests/germline.v0.6.3-beta.forest
	--somatic-forest-model octopus/resources/forests/somatic.v0.6.3-beta.forest
```

### `--min-forest-quality`

Option `--min-forest-quality` specifies the random forest minimum quality score (phred scale) required to PASS a variant call (`RFGQ_ALL`), and each samples genotype calls (`RFGQ`).

```shell
$ octopus -R ref.fa -I reads.bam --min-forest-quality 7
```

### `--use-germline-forest-for-somatic-normals`

Command `--use-germline-forest-for-somatic-normals` specifies that the forest model given to `--forest-model` should be used to score normal sample genotypes for in somatic records, rather than the forest model given to `--somatic-forest-model`.

```shell
$ octopus -R ref.fa -I normal.bam tumour.bam -N NORMAL --forest germline.forest --somatic-forest somatic.forest -o calls.vcf --use-germline-forest-for-somatic-normals
```---
id: introduction
title: Introduction
description: Octopus is a haplotype-based variant caller with multiple calling modes.
sidebar_position: 1
---

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Each calling model is designed for a particular kind of experimental design. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly consider all possible haplotypes at a given loci in reasonable time.---
id: germline
title: Germline WGS
---

This case study considers whole-genome germline variant calling in an individual. We will use the syntheic-diploid sample [CHM1-CHM13](https://www.nature.com/articles/s41592-018-0054-7).

This tutorial is [included](../../static/snakemake/germline.smk) as a [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflow. Run it with

```shell
$ snakemake --snakefile germline.smk -j 16 --use-singularity
```

## Prerequisites

- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [BWA](https://github.com/lh3/bwa)
- [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools)
- [Octopus](https://github.com/luntergroup/octopus) (with [random forests](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-Random-Forest) installed).

This tutorial assumes a directory structure like

```
.
└───data
│   └───references
│   └───reads
│   │   └───raw
│   │   └───mapped
│   └───truth
└───results
    └───calls
    └───eval
```

We'll go ahead and create this upfront:

```shell
$ mkdir -p data/{references,truth} data/reads/{raw,mapped}
$ mkdir -p results/{calls,eval}
```

## Download data

First, download WGS CHM1-CHM13 PCR-free Illumina HiSeq-X10 reads ([PRJEB13208](https://www.ebi.ac.uk/ena/browser/view/PRJEB13208)) from EBI:

```shell
$ curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_1.fastq.gz | gzip > data/reads/raw/CHM1-CHM13_1.fastq.gz
$ curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_2.fastq.gz | gzip > data/reads/raw/CHM1-CHM13_2.fastq.gz
```

Next, we need a reference genome. In this example we'll use GRCh38 with ALT and decoys contigs (hs38DH). The simplest way to get this is with `bwakit`:

```shell
$ run-gen-ref hs38DH
$ mv hs38DH.fa* data/references
```

We'll also need the truth variants to evaluate our calls, which we can get from the [CHM-eval GitHub](https://github.com/lh3/CHM-eval):

```shell
$ curl -L https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar | tar xf - > data/truth/syndip
```

## Map reads

First, we need to index the reference genome for alignment:

```shell
$ samtools faidx data/references/hs38DH.fa
$ bwa index data/references/hs38DH.fa
```

Next, we map the raw reads to the reference genome with `bwa mem`:

```shell
$ bwa mem -t 16 \
     -R "@RG\tID:0\tSM:NA128CHM1-CHM1378\tLB:Broad\tPU:Illumina" \
     data/reference/hs38DH.fa \
     data/reads/raw/CHM1-CHM13_1.fastq.gz data/reads/raw/CHM1-CHM13_2.fastq.gz | \
     samtools view -bh | \
     samtools sort -@ 4 -o data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam -
$ samtools index data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam
```

This should take a few hours.

## Call variants

Now we can call variants with `octopus`. Since this is single-sample germline calling, we'll be using the [individual](../guides/models/individual.md) calling model. We'll also set the [sequence error model](../guides/errorModels.md) to reflect the PCR-free library design of this sample and Illlumina X10 sequencer, and use [random forest filtering](../guides/filtering/forest.md). Finally, we [restrict calling](../guides/advanced/targeted.md) to the primary chromosomes and use the built-in [multithreading](../guides/advanced/threading.md):

```shell
$ octopus \
     -R data/reference/hs38DH.fa \
     -I data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam \
     -T chr1 to chrM \
     --sequence-error-model PCRF.X10 \
     --forest resources/forests/germline.v0.7.4.forest \
     -o results/calls/CHM1-CHM13.hs38DH.bwa-mem.octopus.vcf.gz \
     --threads 16
```

This should complete in a few hours.

## Evaluate variants

Finally, we will evaluate our calls with RTG Tools `vcfeval`. This tool requires the reference sequence to be preprocessed:

```shell
$ rtg format -o ~/reference/hs37d5_sdf ~/reference/hs37d5.fa
```

Then evaluate the calls with `vcfeval`:

```shell
$ rtg vcfeval \
     -t data/reference/hs38DH.sdf \
     -b data/truth/syndip/full.38.vcf.gz \
     --evaluation-regions data/truth/syndip/full.38.bed.gz \
     -c results/calls/CHM1-CHM13.hs38DH.bwa-mem.octopus.vcf.gz \
     -o results/eval/CHM1-CHM13.hs38DH.bwa-mem.octopus.pass.vcfeval \
     -f RFGQ \
     -m annotate \
     --ref-overlap
```

We should see the following results:

```shell
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
```---
id: bamout
title: Realigned BAMs
---

Octopus can generate realigned BAMs that provide visual evidence for why a call has been made. Realigned BAMs are particularly helpful for confirming complex variation where the mapper alignments are incorrect, as can be seen in the IGV pileups below 

<iframe width="750" height="600" src="https://igv.org/web/release/2.8.5/embed.html?sessionURL=blob:rVTbitswEP2VRU8teGU7VytQCm3ptlD6sNvSh1LK2B472tiSI8l2LuTfd.x19sKmEMqGxMTM6Mw5ZzSzZw0aK7ViCzbiEZ8yj9mlbm.grAr8DiVatsigsOgxgxkaVAmyxZ7JlE4s83FEBxSl0duXugR18ebq.uNyHPld7C0FM7AOfl5_69Kdq.zC9.2YQwk7raC1PNGlL_OGx0ZDKpV10tUOuTa5n6PSRMC3uO7h.gfPgFClSnHzuqj0lYScbJ2OQaVngBMaP6L1SKCUduDITut3OB8I52uKmruN4_mO0KGQYP8Hunv87U9zBzE7eKzQSU3NYcnShIvIm80ibz4NLrt_wguigKo5A8mKcn7vmdtWXY9IdN230GPapGjY4lIEwTwUYjSdzCeBEOHB27PaFE8Ytm3LU6OrWG96gtYXo1rd5rjTEwHBcuXfbNUnWfG4hcsSSx5D.T4t3nUUhrvxMuF0E1.WMtpUKs_WmRJ6OytPlKKfPJbLtCnBEeB9hUE1.ZarEpXrO1BggZ8Nrn8sDdJdL.gmB3z0YMjkHAPW8rax2Xwz1jJcResjK504XdX23wY8STjXAJxGtRgF1Bqduiw5UeqVDZieY4COZ5GATbJp1Dhah0dWTZLRRT.t_T52rmyx263Kqhq7KhcwC54X4C5.qZhCj4obMBJ6vYOqWadqIHSNGQ3CxRUqWm9PEGjD0dB1w_Fc_nkjqpJYEnKH.jjwvVYkh4cl2kgrY1lIt_1FId3SAIbdai11A3FB5Ia8gXUY9J_HRj5sGHb4c7gD">
</iframe><br /><br />

Evidence BAMs are requested using the `--bamout` option. The argument to `--bamout` changes slightly depending on whether you're calling one or more samples: If you're only calling a single sample then the argument to `--bamout` is a file path to write the BAM to, e.g.:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout octopus.bam
```

For multiple samples the argument to `--bamout` is a directory path, e.g.:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout minibams
```

Realigned BAMs with the same names as the input BAMs will be written to this directory, so this cannot be a directory where any of the input BAMs are located.

:::important

Realigned BAMs are only available for single-sample BAMs and when `--output` is specified (i.e. no stdout output).

:::

Octopus adds several useful annotations to realigned reads:

| Name        | Description           |
| ------------- |:-------------|
| `HP`      | A list (`,` separated) of haplotype IDs that the read is inferred to originate from. A haplotype ID, which is zero-indexed, corresponds to column in the `GT` field of the affiliated phased VCF. A haplotype ID indicates that the read was unambiguously assigned to the haplotype, while multiple values indicate that the read could equally well be assigned to any of the listed haplotype. |
| `MD`      | Reference free alignments. As defined in the [SAM specficiation](https://samtools.github.io/hts-specs/SAMtags.pdf) |
| `md`      | Like MD but alignment is relative to the inferred haplotype rather than the reference (i.e. mismatches are inferred sequencing errors). | 
| `hc`      | The CIGAR alignment to the inferred haplotype. |
| `PS`      | The phase set the read was assigned to. |

:::tip

The `HP` tag is useful for colouring and grouping alignments in IGV.

:::

By default, only reads supporting regions containing called variation are realigned. However, Octopus can also copy reads overlapping regions where no variation was called using the `--bamout-type FULL` command. Only primary reads are used for BAM realignment.

:::caution

Reads are assigned and realigned to haplotypes called in the `--output` VCF. This means that read-pairs in different phase sets can appear discordant, and reads that are not completely spanned by a phase set (or overlap multiple phase sets) may have poor alignments. Consider trying to [increase haplotype lengths](phasing.md) if this occurs.

:::---
id: haplotypes
title: Haplotype Proposal
---

## Growth

## Lagging

## Backtracking---
id: preprocessing
title: Read Preprocessing
---

The basic idea of read pre-processing is to remove or modify reads containing artifact sequences that are likely to mislead variant calling or increase runtime by introducing spurious candidates.

## Deduplication

Read duplicates can arise during sequencing in several ways, and are library-prep and technology dependent:

![Docusaurus](/img/guides/duplicates.png)

See [here](https://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/) and [here](http://core-genomics.blogspot.com/2016/05/increased-read-duplication-on-patterned.html) for more detailed discussions on how duplicates arise. 

Duplicates can be problematic for variant calling as they can introduce systematic error (e.g. copying errors during PCR). Removing them is usually recommended for WGS libraries, but this remains somewhat [controversial](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3).




---
id: phasing
title: Phasing
------
id: errorModels
title: Error Models
---

Octopus accounts for SNV and indel sequencing errors with a context aware error model. The parameterisation of this model is conditional on the library preparations and sequencing technology used, and can have consequences on calling accuracy, particular for indel errors in tandem repeat regions. Octopus comes packaged with parameter sets for several common library preparation and sequencing combinations, and also allows custom sequence error models to be used.

Built-in error models are selected using the `--sequence-error-model` option, which accepts inputs of the form `[library preparation]<.sequencer>`. library preparation is selected from: `PCR`, `PCR-FREE`, or `10X`. sequencer is selected from: `HISEQ-2000`, `HISEQ-2500`, `HISEQ-4000`, `X10`, `NOVASEQ`, `BGISEQ-5000`. For example, `PCR.NOVASEQ` would select the sequence error model parametrised for a `PCR` library preparation and a `NOVASEQ` sequencer. If no sequencer is provided then the default is used (see `octopus --help`).

Custom error models can be used by providing a path to a valid Octopus error model file. These can be produced using the [`profiler.py`](https://github.com/luntergroup/octopus/blob/develop/scripts/profiler.py) Python script in the scripts top level directory. The script creates error model files given the output of the `--data-profile` command line option.
---
id: discovery
title: Variant Discovery
---

## Pileups

## Local reassembly

## Input VCF---
id: cell
title: Cell
---

The `cell` calling model is used to call germline and somatic variants in single cell and minibatch cell sequencing data. The model attempts to infer local phylogenies for the cells and accounts for allelic biases and dropout often observed in single cell sequencing data.

## Usage: basic

If all of the samples are single cells and none are control cells:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    -o cells.vcf
```

## Usage: with controls

If the experiment includes control cells (e.g. for tumour-normals) then provide the control cell sample names (`--normal-samples`; `-N`):

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --normal-samples CONTROL1 CONTROL2 ... CONTROLM \
    -o cells.vcf
```

All normal cells are assumed to originate from the root (i.e. founder) node of the phylogeny relating cells, and are therefore assumed to all have the same genotype.

## Usage: with minibatchs

If any of the samples are derived from minibatches of cells then specify high dropout concentrations (`--sample-dropout-concentration`) for these samples:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --sample-dropout-concentration MINIBATCH1=100 MINIBATCH2=100 .. MINIBATCHM=100 \
    -o cells.vcf
```

The argument for each minibatch sample may reflect the number of cells contained in the minibatch; the larger the number of cells, the larger the argument value.

The usual use for minibatch samples is for better controls, in which case the minibatches will be normal samples:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --normal-samples MINIBATCH1 MINIBATCH2 ... MINIBATCHM \
    --sample-dropout-concentration MINIBATCH1=100 MINIBATCH2=100 .. MINIBATCHM=100 \
    -o cells.vcf
```

#### VCF output

There are several annotations included in the VCF output:

| Name        | INFO/FORMAT           | Description  |
| ------------- |:-------------:| :-----|
| SOMATIC      | INFO | Indicates that a somatic mutation was inferred (i.e. the phylogeny contains more than one node). |
| PY      | INFO | The MAP phylogeny inferred for the variant loci. This annotation is only added for SOMATIC calls. |
| PPP      | INFO | Posterior probability (Phred) for the MAP phylogeny. |
| PSPP      | INFO | Posterior probabilities (Phred) that the local phylogeny contains `0`,`1`,... nodes |
| PNAP      | FORMAT | Posterior probabilities (Phred) that this sample is assigned to node ID `0`,`1`,.. in the MAP phylogeny (`PY`). |

##### `PY` notation

The phylogeny is serialised using the following algorithm:

```
def serialise(result, node):
    if (node != NULL):
        result += '(' + str(node.id)
        for child in node:
            serialise(result, child)
        result += ')'
```

The algorithm is called with the root node of the phylogeny `serialise("", ROOT)`. Examples:

```
     0
   /   \
  1     2

(0(1)(2))
```

```
     0
       \
        1
         \
          2

(0(1(2)))
```

#### CNV calling

The model can try to identify local copy changes (i.e. deletions or gains of haplotypes). This will result in some samples having called genotypes with different ploidies to the default ploidy. The maximum number of gains and losses is specified with the `--max-copy-gain` and `--max-copy-loss` options, respectively. For example, to identify up to one copy gain or loss:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --max-copy-gain 1 --max-copy-loss 1 \
    -o cells.vcf
```
**Warning** calling copy gains is currently computationally very expensive.

#### Performance considerations

A critical parameter for this calling model is the maximum size of the phylogeny (`--max-clones`). Copy loss and gain calling are also computationally expensive.

It is recommended to allow automatic thread usage with this calling model (use `--threads` option without an argument).---
id: cancer
title: Cancer
---

The `cancer` calling model is for calling germline variation and somatic mutations in tumours. The model can jointly genotype multiple tumours from the same individual, and make use of a normal sample for improved classification power.

## Usage: tumour-normal pairs

To call germline and somatic mutations in a paired tumour-normal sample, just specify which sample is the normal (`--normal-sample`; `-N`):

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam --normal-sample NORMAL
```

It is also possible to genotype multiple tumours from the same individual jointly:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumourA.bam tumourB.bam -N NORMAL
```

## Usage: tumour only

If a normal sample is not present the cancer calling model must be invoked explicitly:

```shell
$ octopus -R hs37d5.fa -I tumourA.bam tumourB.bam -C cancer
```

Be aware that without a normal sample, somatic mutation classification power is significantly reduced.

## VCF output

By default both germline and somatic variants are called, somatic mutations are tagged with the `SOMATIC` INFO field. The `GT` fields for `SOMATIC` variants and any other variants in the same phase set (`PS`) are augmented with the number of unique somatic haplotypes inferred (somatic haplotypes are identified with the `HSS` tag). For example, in the following VCF:

```
##fileformat=VCFv4.3
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Record includes a somatic mutation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set">
##FORMAT=<ID=HSS,Number=.,Type=Integer,Description="Somatic status for each haplotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOUR
1 50  . G  A   . . SOMATIC GT:PS:HSS   0|0:50:0,0    1|0|0:50:1,0,0
1 100 . A  C   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|0|1:100:0,0,1,1
1 150 . G  T   . . .       GT:PS:HSS   1|0:100:0,0   1|0|1|0:100:0,0,1,1
1 200 . A  AC  . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,0,1,1
1 250 . TA T   . . .       GT:PS:HSS   1|1:100:0,0   1|1|1|1:100:0,0,1,1
2 300 . C  A,T . . SOMATIC GT:PS:HSS   0|2:300:0,0   1|0|2:300:1,0,0
2 350 . G  T   . . .       GT:PS:HSS   0|1:300:0,0   0|0|1:300:1,0,0
2 400 . T  C   . . SOMATIC GT:PS:HSS   0|0:300:0,0   1|0|0:300:1,0,0
3 100 . A  C   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|1:100:0,0,1,1
3 150 . C  CC  . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,0,1,1
3 200 . G  T   . . .       GT:PS:HSS   0|1:100:0,0   0|1|1|1:100:0,0,1,1
4 100 . C  T   . . .       GT:PS:HSS   0|1:100:0,0   0|1|1|1:100:0,1,1,0
4 150 . A  G   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,1,1,0
4 200 . A  G   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|1|0|0:100:0,1,1,0
```

The first phase set `1:50` includes a simple somatic mutation; it is not phased with any other variants. Downstream of this is another phase block starting at `1:100` that includes 2 germline and 2 somatic variants. The first somatic mutation in this phase set at `1:100` is phased onto the germline haplotype including the reference allele at the germline variant at `1:150`. The second somatic mutation in this phase set is phased with the alternate allele of this germline variant. In the third phase set stating at `2:300` there is a somatic mutation that segregates with a germline variant at the same position, and another somatic mutation which is phased onto the same germline haplotype - the reference. The somatic allele in the multiallelic record is determined by looking at the `HSS` flags (1 indicates somatic), so "C>A" is the somatic mutation in this case. The third phase set starting at `3:100` includes 2 somatic mutations phased onto the same germline haplotype, but the first somatic mutation was inferred to segregate with both the germline allele and somatic allele at `3:150`, suggests linear progression of the C>CC somatic mutation. The last phase set starting at `4:100` includes 2 somatic mutations phased onto the same germline haplotype, but not with each other.

## `QUAL` vs `PP`

For both paired and tumour-only calling, octopus reports two quality scores for each call (both germline and somatic):

  * `QUAL` is the posterior probability the variant is *segregating* in the sample regardless of somatic classification.
  * `PP` (an `INFO` field) is the posterior probability the variant is segregating **and** classified correctly.

The difference between `QUAL` and `PP` indicates the uncertainty in the calls classification; a call may have high `QUAL` but low `PP` if the classification is uncertain (common in tumour-only calling or if the normal coverage is low). `PP` should always be less than `QUAL` in theory.

## `HPC`, `MAP_HF` and `HF_CR`

Octopus infers a probability distribution over haplotype frequencies, including any somatic haplotypes. For each variant in a phase-block with a `SOMATIC` variant, three statistics are reported that relate to haplotype frequency (per haplotype).

  * `HPC` is the Dirichlet posterior pseudo-count. Note that this count includes the prior count so should not be taken to mean the empirical count.
  * `MAP_HF` (`FORMAT`) is the [Maximum a Posteriori](https://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation) haplotype frequency point estimate.
  * `HF_CR` (`FORMAT`) is a [credible interval](https://en.wikipedia.org/wiki/Credible_interval) of the haplotypes frequency. The mass of the credible interval is specified by `--credible-mass`. The credible interval gives you an indication of how certain the MAP_HF estimate is; A very narrow interval means the MAP estimate is very certain, a wide interval means it is uncertain.

## Performance considerations

* The number of genotypes considered by the model `--max-genotype` has a significant impact on overall runtime.
* The parameter `--max-somatic-haplotypes` controls the maximum number of unique segregating somatic haplotypes to be modelled. There must be at-least one somatic haplotype, but adding more can resolve somatic mutations falling on different germline haplotypes or multiple distinct haplotypes due to sub clonal evolution. 

## SOMATICs only

To report only `SOMATIC` calls just add the `--somatics-only` command. This will only work if calls are beings filtered already (i.e. it will have no effect if `-f off`). It is generally not recommended to use this option until you are 100% satisfied with your calls, as call filtering will not work correctly if germline variants have been filtered from the VCF; you will not be able to re-filter your SOMATIC-only VCF.---
id: trio
title: Trio
---

The `trio` calling model is for calling germline variation and _de novo_ mutations in parent-offspring trios.

## Usage

To call germline and de novo mutations in a trio, either specify both maternal (`--maternal-sample`; `-M`) and paternal (`--paternal-sample`; `-F`) samples:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891
```

or provide a PED file which defines the trio:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam --pedigree ceu_trio.ped
```

## Setting ploidies

The trio calling model is currently only fully supported for diploid samples. You can set sample contig ploidies with the `--contig-ploidies` (`-p`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891 -p NA12891:X=1
``` 

## *De novo* mutations only

To call only `DENOVO` mutations, just add the `--denovos-only` command:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam --ped ceu_trio.ped --denovos-only
```

Note this only works for filtered calls, and should only be used if you do not plan on re-filtering the calls, since this will not be possible once the germline mutations are removed.

## Performance consideration

Like for the population model, `--max-joint-genotypes` has a large baring on overall accuracy and runtime. Increasing this parameter can improve sensitivity, but will also result in longer runtimes. The default value is set with accuracy in mind.---
id: polyclone
title: Polyclone
---

The `polyclone` calling model is for calling variation in a pooled sample of haploid clones where the number and mixture composition of clones is unknown, as is sometimes the case in bacteria or virus sequencing studies. 

## Usage

If your sample contains an unknown mix of haploid clones (e.g. some bacteria or viral samples), use the `polyclone` calling model:

```shell
$ octopus -R H37Rv.fa -I mycobacterium_tuberculosis.bam -C polyclone
```

This model will automatically detect the number of subclones in your sample (up to the maximum given by `--max-clones`).

## Performance considerations

The most important parameter in this model is `--max-clones` which determines the maximum number of haplotypes that octopus will try to use to fit the data. This is a bit like setting the 'ploidy' for the sample. Higher values of `--max-clones` may lead to longer runtimes and more memory usage - we do not recommend values above 10.---
id: individual
title: Individual
---

The `individual` calling model is used to call germline variants in a single sample with known ploidy. It is the simplest model Octopus offers.

## Usage

If the file `NA12878.bam` contains a single sample, to call germline variants in all regions use:

```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam
```

or less verbosely:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam
```

By default, octopus automatically detects and calls all samples contained in the input read files. If your BAM file contains multiple samples, and you want to call just one of these, use the `--samples` (`-S`) option:

```shell
$ octopus -R hs37d5.fa -I multi-sample.bam -S NA12878
```

## Setting the ploidy

Octopus assume diploid samples by default. If your sample is not diploid you can set the ploidy with the `--organism-ploidy` (`-P`) option:

```shell
$ octopus -R hs37d5.fa -I haploid.bam -P 1
```

You can also set contig specific policies with the `--contig-ploidies` (`-p`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -p Y=1
```

Octopus automatically sets contigs `Y` and `chrY` to ploidy 1.

_Performance note_: Larger ploidies require greater computational resources. In general, ploidies above 4 are currently intractable. If you wish to call tetraploid sample, you may find that you need to tweak other performance related parameters to get reasonable run times.---
id: population
title: Population
---

The `population` calling model is for jointly calling germline variation in small cohorts of samples with known ploidy but unknown pedigree structure. 

## Usage

Multiple samples from the same population, without pedigree information, can be called jointly:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam
```

Joint calling samples may increase calling power, especially for low coverage sequencing.

## Setting ploidies

Octopus can currently joint call samples that have the same ploidies for all contigs. You can set ploidies with the `--organism-ploidy` (`-P`):

```shell
$ octopus -R hs37d5.fa -I haploid1.bam haploid2.bam haploid3.bam -P 1
```

and `--contig-ploidies` (`-p`) options

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -p X=1
```

Octopus sets the organism ploidy to 2 by default and also sets chromosome Y to ploidy 1. X is left as two, so if your samples are all male, you could set `-p X=1`.

## Performance considerations

For joint calling, the single most important parameter which will determine performance is [`--max-genotype-combinations`](cli.md#--max-genotype-combinations). This parameter determines the maximum number of joint genotype combinations that octopus can consider. The number of possibl genotype combinations at a given site is exponential in the number of samples, so this value will usually be significantly lower than the number of possible joint genotypes. Octopus will select the 'best' `--max-genotype-combinations` genotype combinations using approximation. Increasing this value may improve accuracy, but result in longer runtimes and memory use.

## Calling large cohorts, and the n + 1 problem

Joint calling many samples (> 10), especially high coverage ones, is computationally expensive. This is because true joint calling necessarily requires considering all possible genotypes for each sample, which grow approximately linearly with the number of samples (due to noise). The number of possible genotype combinations is also polynomial in the number of possible genotypes.

There is also the so-called 'n + 1 problem' of adding an extra sample to a cohort already jointly called. Ideally this would be done cheaply rather than having to recall all samples again. Other methods have addressed this problem by individually calling all samples, then jointly calling all samples by evaluating previously computed genotype likelihoods under a joint model. However, this approach will never correctly call cases where a true allele was not proposed for a specific individual (since no genotype likelihood will be computed). This method also suffers the same 'representation problems' as typical merge procedures as it does not consider haplotypes.

Although octopus does not specifically address the n + 1 problem, we suggest a procedure for jointly calling many samples:

1. Call each sample individually.
2. Group the samples into small subsets (< 20), joint call these using only the variants called previously.
3. Repeat step 2 until there is only one group left (containing all samples).

Adding another sample is achieved in a similar way:

1. Call the new sample individually.
2. Joint call all samples using the variants called in previous joint calling and the new sample.

The reason these procedures is cheaper than joint calling de-novo is the number of candidate alleles is likely significantly reduced as most noise should be removed during individual calling.

For example, three samples could be jointly called with this procedure as follows:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.NA12878.bcf
$ octopus -R hs37d5.fa -I NA12891.bam -o octopus.NA12891.bcf
$ octopus -R hs37d5.fa -I NA12892.bam -o octopus.NA12892.bcf
$ octopus -R hs37d5.fa -I NA12892.bam NA12891.bam NA12892.bam \
 --disable-denovo-variant-discovery \
 -c octopus.NA12878.bcf octopus.NA12891.bcf octopus.NA12892.bcf \
 -o octopus.joint.bcf
```

Note the last step disables candidate generation from raw alignments (`-g`) and local reassembly (`-a`), and only uses provided source variants (`-c`).---
id: forest
title: Random Forest
---

Octopus provides a powerful way to classify variant calls with random forests using the [Ranger](https://github.com/imbs-hl/ranger) random forest library. Pre-trained random forests are available on [Google Cloud](https://console.cloud.google.com/storage/browser/luntergroup/octopus/forests/?project=parabolic-eon-208710). Currently there are forests for germline and somatic variants. You can easily obtain the forests using the Python install script:

```shell
$ ./scripts/install.py --forests
```

will download the forests into `octopus/resources/forests/`. The provided forests are suitable for calling typical diploid germline and cancer data. They *may* work well for other types of samples but this is untested. 

## Using random forest filtering

### Germline variant random forest filtering

To filter germline variants using the germline random forest, just specify the path to the forest in the `--forest-file` option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    --forest resources/forests/germline.forest
```

All calls will be annotated with the `FORMAT` field `RFGQ`, which is the phred-scaled probability of the `GT` call being correct according to the forest model. Each record will also be annotated with an `INFO` field, `RFGQ_ALL`, which is the phred-scalled probability that all `GT` calls in the record are correct, this probability is used to filter calls (see [`--min-forest-quality`](https://github.com/luntergroup/octopus/wiki/Command-line-reference#option---min-forest-quality)) with the `RF` filter.

### Somatic variant random forest filtering

When calling germline and somatic variants, you need to provide both random forests:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \
    --forest resources/forests/germline.forest \
    --somatic-forest resources/forests/somatic.forest
```

However, if you are calling somatic variants only (i.e. using the `--somatics-only` flag), then you just need to provide the somatic forest:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \
    --somatics-only \
    --somatic-forest resources/forests/somatic.forest
```

### Re-filtering an Octopus VCF with random forests

Random forest filtering can be applied to a VCF file produced by Octopus without recalling with the `--filter-vcf` option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    --filter-vcf unfiltered.vcf.gz \
    --forest resources/forests/germline.forest \
    -o filtered.vcf.gz
```

Note this will overwrite any existing filtering information. Note, the VCF file cannot be a `--legacy` VCF file, and should contain all variant calls that Octopus would make by default (i.e. do not use a VCF produced with the `--somatics-only` option).

## Training random forests

The supplied random forest should perform well for the majority of users, however, it is possible that performance can be improved by training a new forest model on your own data.

In addition to Octopus requirements, training new forest models requires the following:

1. A truth set of variants and high-confidence regions (e.g., GIAB or SynDip).
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/).
3. [RTG Tools](https://www.realtimegenomics.com/products/rtg-tools).

Forest models can - and ideally should - be trained on multiple examples runs (i.e., a run of Octopus). Each example can itself be calls for a single sample, or joint calls for multiple samples. For joint calls, each sample is used to generate training data independently, and a subset of samples can be selected.

Training new forests is done using the bundled Snakemake script [forest.smk](https://github.com/luntergroup/octopus/blob/develop/scripts/forest.smk). Basic usage is like:

```shell
$ snakemake \
    --snakefile forest.smk \
    --configfile config.yaml \
    --cores 20
```

You can of course use any [Snakemake option](https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html). Note that the RTG Tools binary `rtg` must be in your path.

The configuration file (YAML or JSON format) provided to Snakemake specifies the training data and setup. The format of the config file is:

* `truths` is a dictionary of truth sets, the key being some unique label (e.g. `SAMPLE1.truth`). The label does **not** need to correspond to the sample name. Each truth set contains a sub dictionary with `vcf` and `bed` value pairs corresponding to the truth variants and high confidence regions, respectively.
* examples is a list of examples to use. Each example is dictionary with the following fields:
    - `name`: the filename prefix to use for this example. If this is not specified then the name is automatically set by concatenating input BAM names. However, this can result in long filenames that may cause OS errors if many inputs are used.
    - `reference`: the reference fasta file. [**required**]
    - 'reads': A list of BAM files to use. The list may be omitted if a single BAM is used. [default none]
    - `truth`: Either the label of the truth set for single sample calling, or a dictionary of sample-truth pairs, where the key is the sample name (in the BAM file), and the value is one of the keys in `truths`. [**required**]
    - `regions`: A bed file to use for calling (i.e. sets option `--regions-file`). [default none]
    - `options`: Additional options to pass to the octopus command.
    - `tp_fraction`: The fraction of true positive calls to use for training. [default 1]
    - `fp_fraction`: The fraction of false positive calls to use for training. [default 1]
    - `threads`: How many threads to use for analysis in this run.
* `training` is a dictionary of forest training setting:
    - `training_fraction`: the fraction of calls to use for training.
    - `hyperparameters`:the random forest hyperparameters to use for training:
        * `trees`: the number of trees to use (ranger option `--ntrees`).
        * `min_node_size`: minimum size of a node in the tree (ranger option `--targetpartitionsize`).
        * `maxdepth`: maximum depth of each tree (ranger option `--maxdepth`).

A minimal example is:

```yaml
examples:
    -
        reference: /path/to/reference.fa
        reads: /path/to/reads.bam
        truth: GIAB//GRCh38//HG001
```

Here, `truth` is set to a special label `GIAB//GRCh38//HG001`. The script accepts special labels like this for any of the GIAB truth sets - it will download all the required truth data. A more detailed example is: 

```yaml
truths:
    SAMPLE1.truth:
        vcf: /path/to/truth/variants/sample1.vcf.gz
        bed: /path/to/truth/variants/sample1.bed
    SAMPLE2.truth:
        vcf: /path/to/truth/variants/sample2.vcf.gz
        bed: /path/to/truth/variants/sample2.bed
examples:
    -
        reference: /path/to/reference.fa
        reads: /path/to/reads1.bam
        truth: SAMPLE1.truth
        tp_fraction: 0.5
        fp_fraction: 1
        threads: 10
    -
        name: joint
        reference: /path/to/reference.fa
        reads:
            - /path/to/reads1.bam
            - /path/to/reads2.bam
        regions: /path/to/calling/regions.bed
        truth:
            SAMPLE1: SAMPLE1.truth
            SAMPLE2: SAMPLE2.truth
        options: --config /path/to/octopus/octopus.config
        threads: 20
training:
    training_fraction: 0.25
    hyperparameters:
        -
            trees: 200
            min_node_size: 20
```

The default behaviour of the script is to use all calls for training, however if the option `--kind=somatic` is used then only variants called `SOMATIC` are used for training. This can be used to generate somatic forests.---
id: thresholds
title: Hard Thresholds
---

Octopus currently provides simple threshold based filtering. A number of measures that can be used to define a Boolean filter expression. Currently, Octopus only supports expressions with OR operations and the less than (<) and greater than (>) comparators. Furthermore, measure name must appear on the left hand side of each condition in the expression.---
id: annotations
title: Introduction
---

Octopus calls variants and genotypes using a Bayesian model. Like any model, there are some aspects of real error that are not fully captured in the model which can lead to false calls. We therefore recommended filtering calls to improve precision.

There are currently two approaches available in Octopus for filtering variants:

1. [Hard coded thresholds](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-thresholds)
2. [Random forests](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-Random-Forest)

Both methods use [annotations](#annotations) computed by Octopus.

The random forest approach is preferred when sufficient training data is available (e.g. typical germline and somatic calling). Hard coded thresholds are appropriate for other cases (e.g. UMI low frequency calling).

## Annotations

Octopus provides various annotations for filtering variants. To list available annotations, use the command 

```shell
$ octopus --help --annotations
```

For example, the current annotations are:

```shell
Name	Kind	Number	Type	Description
AC	FORMAT	A	Integer	"Number of non-reference alleles called for each sample"	
AD	FORMAT	R	Integer	"Empirical allele depth"	
ADP	FORMAT	1	Integer	"Number of reads overlapping the position that could be assigned to an allele"	
AF	FORMAT	R	Float	"Empirical allele frequency (AD / ADP)"	
AFB	FORMAT	R	Float	"Absolute difference between empirical allele frequency (AF) and expected allele frequency given genotype"	
AMQ	FORMAT	R	Integer	"Median mapping quality of reads assigned to each allele"	
ARF	FORMAT	1	Float	"Fraction of reads overlapping the call that cannot be assigned to a unique haplotype"	
BMC	FORMAT	1	Integer	"Number of base mismatches at variant position in reads supporting variant haplotype"	
BMF	FORMAT	1	Float	"Fraction of base mismatches at variant position"	
BMQ	FORMAT	1	Integer	"Median quality of base mismatches in reads assigned to a unique allele"	
BQ	FORMAT	R	Integer	"Median base quality of reads supporting each allele"	
CC	INFO	1	Float	"PP divided by QUAL"	
CRF	INFO	1	Float	"Fraction of clipped reads covering the call"	
DAD	FORMAT	R	Integer	"Number of realigned reads supporting ALT alleles identified as duplicates"	
DAF	FORMAT	R	Float	"Fraction of realigned reads supporting ALT alleles identified as duplicates"	
DC	INFO	1	Float	"Number of reads supporting a de novo haplotype in the normal"	
DENOVO	FORMAT	1	Integer	"DENOVO status of each sample"	
DP	FORMAT	1	Integer	"Number of read overlapping the call"	
DPC	FORMAT	1	Float	"Concordance of allele support from duplicate reads"	
ER	FORMAT	R	Float	"Error rate in supporting reads"	
ERS	FORMAT	R	Float	"Error rate standard deviation in supporting reads"	
FRF	FORMAT	1	Float	"Fraction of reads filtered for calling"	
GC	INFO	1	Float	"GC bias of the reference in a window centred on the call"	
GQ	INFO	A	Integer	"Number of non-reference alleles called for each sample"	
GQD	FORMAT	1	Float	"GQ divided by DP"	
ITV	INFO	A	Flag	"Is the variant a transversion"	
MC	FORMAT	1	Integer	"Number of mismatches at variant position in reads supporting variant haplotype"	
MF	FORMAT	1	Float	"Fraction of reads with mismatches at variant position"	
MHL	FORMAT	.	Integer	"Mean likelihood (Phreds) of reads overlapping the site assigned to each haplotype"	
MP	INFO	1	Float	"Model posterior for this haplotype block"	
MQ	INFO	1	Float	"Mean mapping quality of reads overlapping the call"	
MQ0	INFO	1	Integer	"Number of reads overlapping the call with mapping quality zero"	
MQD	FORMAT	1	Integer	"Maximum pairwise difference in median mapping qualities of reads supporting each haplotype"	
MRC	FORMAT	1	Integer	"Number of reads supporting the call that appear misaligned"	
MRL	FORMAT	1	Integer	"Maximum read length overlapping the site"	
NC	INFO	1	Float	"Fraction of overlapping reads supporting a somatic haplotype in the normal"	
PLN	FORMAT	1	Integer	"Length of the phase block for the call"	
PP	INFO	1	Float	"Call posterior probability"	
PPD	INFO	1	Float	"PP divided by DP"	
QD	INFO	1	Float	"QUAL divided by DP"	
QUAL	INFO	A	Integer	"Number of non-reference alleles called for each sample"	
REB	FORMAT	1	Float	"Probability allele occurs at the end (head or tail) of supporting reads"	
REFCALL	FORMAT	1	Integer	"REFCALL status of each sample"	
RSB	FORMAT	1	Float	"Bias of variant side (head or tail half) in supporting reads"	
RTB	FORMAT	1	Float	"Probability allele occurs in the tail of supporting reads"	
SB	FORMAT	1	Float	"Strand bias of reads based on haplotype support"	
SD	FORMAT	1	Float	"Strand bias of reads overlapping the site; probability mass in tails of Beta distribution"	
SF	FORMAT	1	Float	"Max fraction of reads supporting ALT alleles that are supplementary"	
SHC	FORMAT	1	Integer	"Number of called somatic haplotypes"	
SMQ	FORMAT	1	Integer	"Median mapping quality of reads assigned to called somatic haplotypes"	
SOMATIC	FORMAT	1	Integer	"SOMATIC status of each sample"	
STRL	INFO	1	Integer	"Length of overlapping STR"	
STRP	INFO	1	Integer	"Period of overlapping STR"	
VL	FORMAT	1	Integer	"Maximum length of called alleles"
```---
id: threading
title: Multithreading
---

Octopus has built in multithreading capabilities, activated with the `--threads` option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads
```

This will let octopus automatically decide how many threads to use (usually the number of available cores). However, the number of threads to use can also be specified explicitly:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads 4
```

## Discussion

Using more than one threads is the simplest way to speed up variant calling. However, there is not always a linear payoff in runtime and the number of threads. Optimising thread throughput is challenging as it is highly data dependent. The best case scenario is to have each thread assigned single tasks of equal complexity so that they all finish simultaneously. Internally, octopus divides the calling region into chunks (genomic regions) with approximately the same number of reads. In particular each chunk has [`--target-read-buffer-memory`](https://github.com/luntergroup/octopus/wiki/Command-line-reference#option---target-read-buffer-memory) / `--threads` worth of reads. Allowing octopus more read memory therefore increases the size of each chunk (assuming constant number of threads). Here are some consequences of larger chunks to consider:

* Less thread management overhead (good).
* Less IO accesses (good).
* Higher chance of unbalanced workload (bad).
* More thread blocking due to IO constraint (bad).

Some experimentation may be required to find the best thread/memory combination for your particular data.---
id: targeted
title: Targeted Calling
---

By default, octopus will call all contigs specified in the reference index. However, calling can be restricted to a subset of regions by providing a list of regions to call or skip. Note all input regions are assumed to be zero-indexed. If you're using one-indexed regions then add the `--one-based-indexing` option (applied to **all** input regions).

#### `--regions` (`-T`) and `--skip-regions` (`-K`)

Provide a list of regions directly to the command line to call or not call. The format is `<chr>[:start][-[end]]`. So the following are valid:

  * `chr1`: all of `chr1`.
  * `chr2:10,000,000`: the single position `10000000` in `chr2`.
  * `chr3:5,000,000-`: everything from `chr3:5,000,000` onwards.
  * `chr4:100,000,000-200,000,000`: everything between `chr4:100,000,000` and `chr4:200,000,000`. The interval is half open so position `chr4:200,000,000` is **not** included.

You can provide multiple regions with this option, for example:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    -T 1 2:30,000,000- 3:10,000,000-20,000,000
```

Conversely the `--skip-regions` is for providing a list of regions **not** to call. The format is exactly the same as for `--regions`, so:

```console
$ octopus -R hs37d5.fa -I NA12878.bam -K Y
```

Will call all contigs in the reference index other than `Y`. You can provide both `--regions` and `--skip-regions` together, in which case the complement of `--region` and `--skip-regions` will be used.

:::tip

The `--region` option accepts a special ranged argument in the form `<lhs> to <rhs>` where `lhs` and `rhs` are contig names or positions. The reference index is used to expand the range. For example `chr1 to chr3` is expanded to `chr1` `chr2` `chr3`. You can also specify positional start or end points, e.g., `1:10,000,000 to X`. Only one region range can be specified. If it so happens that one of your reference contigs is called `to` then you cannot use this feature!

:::

#### `--regions-file` (`-t`) and `--skip-regions-file` (`-k`)

These commands accept a file containing a line separated list of regions to call or skip. The regions can either be in the same format as `--regions`, or in BED format.

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -t autosomes.txt
$ octopus -R hs37d5.fa -I NA12878.bam -k avoid.bed
```
---
id: configs
title: Config files
---

The `--config` command line option is a handy way configure Octopus. The argument is a text file containing settings for any subset of options (other than `config` itself). Each line of the file contains a parameter setting in the format

```shell
parameter = argument
```

Comment lines are allowed in the config file are proceeded with `#`.

It is perfectly fine to specify a config file and explicit command line options, e.g.:

```shell
$ octopus --config my-config.config -R reference.fa -I reads.bam -o octopus.vcf
```

Explicit command line options that are in the config file are ignored.

The [configs](https://github.com/luntergroup/octopus/tree/develop/configs) directory in the main source tree will be used to store useful config files.---
id: memory
title: Memory Use
---

#### File based factors

Octopus has been designed to minimise disk accesses by buffering frequently accessed resources in main memory. In other words, rather than making lots of short disk reads and using little main memory, octopus prefers to making few large disk reads, and use more main memory. This is beneficial for a number of reasons, especially in a cluster setting. There are essentially two resources which octopus actively buffers: read data and reference sequence. Both these can be controlled by the user:

* `--max-reference-cache-footprint` (`-X`) controls the reference buffer size.
* `--target-read-buffer-footprint` (`-B`) controls the read buffer size. This is not a hard limit, but for most normal samples octopus will respect this.

Both options are specified in all the standard memory units (e.g. `500Mb`, `2G`, `4GB`, etc). Note both buffer sizes are shared between threads, so `--target-read-buffer-footprint=2Gb` using two threads would mean 1GB per thread. This is important is it determines the size of thread 'jobs', and can therefore have a significant impact on throughput and overall runtime.

#### Other factors

The other factors to consider when optimising memory usage are:

* Multithreading: More threads means more memory overhead. A good starting point is to 'budget' an extra 100MB per additional thread, however, this will also depend on your data and the factors below, so some trial and error may be required.
* Calling model: Simpler models use less memory.
* Model setup: The parameterisation of the calling model can have a large influence on short term memory usage. For example setting `--max-joint-genotypes` (in the trio or population model) to a very large memory may mean requests for large memory blocks.---
id: vcf
title: VCF Format
---

Octopus reports variants in [VCF 4.3 format](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).

## Haplotypes

Octopus always reports phased genotypes (`GT` separated with `|` rather than `/`). The extent of phasing is provided in the `PS` `FORMAT` field, which refers to the `POS` of a previous record, so:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	.	A	C	.	.	.	GT:PS	1|0:100
1	200	.	G	T	.	.	.	GT:PS	0|1:100
```

indicates the variants at positions `100` and `200` are phased. While

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	.	A	C	.	.	.	GT:PS	1|0:100
1	200	.	G	T	.	.	.	GT:PS	0|1:200
```

indicates the variants are unphased. Note that phase sets may not be contiguous.

## Spanning alleles

To represent complex multi-allelic loci, Octopus prefers a decompose alleles into multiple records and use the `*` allele to resolve conflicts. In particular, Octopus always splits variants that have unique `REF` alleles into multiple VCF records. For example, two overlapping deletions are represented like:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	102738191	.	ATTATTTAT	A,*	.	.	.	GT	1|2
1	102738191	.	ATTATTTATTTAT	A,*	.	.	.	GT	2|1
```

In contrast, some other tools may choose to represent these variant with

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTAT	A	.	.	.	GT	1/0
1	102738191	.	ATTATTTATTTAT	A	.	.	.	GT	1/0
```

which is inconsistent as the reference is deduced in each record, or

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTATTTAT	ATTAT,A	.	.	.	GT	1/2
```

which is consistent, but rapidly becomes unmanageable as the length and number of overlapping variants increases.

### Multiple `*` records

In some cases, a site may overlap with multiple distinct alleles. Octopus represents these sites using multiple `*` alleles in the `ALT` field.

Consider the [realignments](../bamout.md) from the [Ashkenazim trio](https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials):

<iframe width="750" height="700" src="https://igv.org/web/release/2.8.5/embed.html?sessionURL=blob:rZRva9swEMa_StGrDVzZcZLZDozBBmsHYy.yjb0YY8jy2VYiS64k2_lDvnvPatJ0bIUwGhKTSKdfnufudHvSg7FCK7IgMU3pnATE1nr4yppWwhfWgCWLkkkLATFQggHFgSz2RBR4oq6mKR5QGIa_bruGqatXN8sP9TQNx73XuFky69j35ecx3LnWLsLQTilr2E4rNljKdROKqqe50awQyjrhOgdUmyqsQGkUEFq48zj_oCVDqlAFbF6Wim.BZL51OmequACONHqieRJTSjvmMJ02HDnvkfOpAE3dxtFqh3QmBbP_gx4fv_1p6lhODgGRmndYHMJrM1tMsiCexkESpdenr2mC_.cM42uM.rknbtuOVULbnS9iQLQpwJDFdRZFySTL4vksmUVZNjkEe9IZ.UTjMAy0MLrN9cZLtOFGDNtVtlt1nZg3fR_e3kRRTHPWvCvk2.hJU5zW_120v8Hrsl8XQ5mUd0NVRuYMxo84wUttGuaQ8wA.WsP0VKoB5XyiJUj4aODuW20AW1piw0Y0fnQ9u8RllXUy1jYxepOUzYPL6TMu_fqlLlet4Mk0Wsluu2tbdga_sMv5JS53XaPrbcpz3oGsuRcze8alX7_Upd7VrQFd5yu.XsflGfzCLt.MLo8il1Bii1_dgMLR9YSO0wsv1Nj2f6bjsuuneC6QPFLPl9nnAFDLcUD2wopcSOG2P3BLD3i1JuPYbHTPconijnFH1ZPIv86WH6cHOfw63AM-">
</iframe><br /><br />

The variants called by Octopus in this region are:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002	HG003	HG004
chr4	19232687	.	A	G	50	PASS	AC=2;AN=6	GT:PS	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232728	.	ATCTG	A	50	PASS	AC=2;AN=6	GT:PS:PQ	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232732	.	GTCTGTCTATCTA	G,*,*	50	PASS	AC=2,2,2;AN=6	GT:PS	2|3:19232687	1|2:19232687	1|3:19232687
chr4	19232736	.	G	GTCTATCTA,*	50	PASS	AC=2,2;AN=6	GT:PS	1|0:19232687	2|19232687	2|0:19232687
chr4	19232762	.	C	CTA	550	PASS	AC=2;AN=6	GT:PS	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232762	.	CTATATA	C	50	PASS	AC=2;AN=6	GT:PS	0|0:19232687	1|0:19232687	1|0:19232687
chr4	19232858	.	A	ACT	50	PASS	AC=4;AN=6	GT:PS	0|1:19232687	1|0:19232687	1|1:19232687
```

Notably, Octopus reports two `*` alleles at the `chr4:19232732` record. This indicates that there are two distinct alleles overlapping this site (at `chr4:19232728` and `chr4:19232736`). 

:::caution

In some rare cases, Octopus may report double `*` records in a single sample, resulting in calls like

```
1	100	.	ACGT	*,*	50	PASS	AC=1,1;AN=2	GT	2|3
```

This generally indicates a homozygous deletion directly upstream of a heterozygous deletion, like:

```
1	99	.	TA	T	50	PASS	AC=2;AN=2	GT:PS	1|1:99
1	100	.	ACGT	*,*	50	PASS	AC=1,1;AN=2	GT:PS	2|3:99
```

:::

---
slug: welcome
title: Welcome
author: Daniel Cooke
author_title: Octopus author and Bioinformatics Engineer @ Invitae
author_url: https://github.com/dancooke
author_image_url: https://github.com/dancooke.png
tags: [welcome, octopus]
---

Welcome to the new Octopus website!
