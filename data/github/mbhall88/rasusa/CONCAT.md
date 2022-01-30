# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.1]

### Added

- Warning if the actual coverage of the file(s) is less than the requested coverage
  [[#36][36]]
- JOSS manuscript

### Changed

- Use `rasusa` as the entry command for docker container [[#35][35]]

## [0.6.0]

### Addedd

- `--bases` option to allow for manually setting the target number of bases to keep
  [[#30][30]]
- `--genome-size` can now take a FASTA/Q index file and the sum of all reference
  sequences will be used as the genome size [[#31][31]]

## [0.5.0]

### Added

- Support for LZMA, Bzip, and Gzip output compression (thanks to
  [`niffler`](https://github.com/luizirber/niffler/)). This is either inferred from the
  file extension or manually via the `-O` option.
- Option to specify the compression level for the output via `-l`

### Changed

- Use a `Vec<bool>` instead of `HashSet` to store the indices of reads to keep. This
  gives a nice little speedup (see [#28][28]), A big thank you to
  [@natir](https://github.com/natir) for this.

### Fixed

- Restore compression of output files [[#27][27]]

## [0.4.2]

### Fixed

- I had stupidly forgetten to merge the fix for [#22][22] onto master 🤦

## [0.4.1]

### Fixes

- Releasing cross-compiled binaries didn't work for version 0.4.0
- Docker image is now correctly built

## [0.4.0]

### Changed

- Switch from using `snafu` and `failure` for error handling to `anyhow` and
  `thiserror`. Based on the procedure outlined in [this excellent blog
  post][error-blog].
- Switched fasta/q parsing to use [needletail](https://github.com/onecodex/needletail)
  instead of rust-bio. See [benchmark] for improvement in runtimes.
- Changed the way Illumina paired reads are subsampled. Previously, there was an
  assumption made that the reads of a pair were both the same length as the R1 read. We
  are now more careful and look at each read's length individually [[#22][22]]
- Moved container hosting to quay.io

## [0.3.0]

Version 0.3.0 may give different results to previous versions. If so, the differences
will likely be a handful of extra reads (possibly none). The reason for this is
`--coverage` is now treated as a float. Previously we immediately round coverage down to
the nearest integer. As the number of reads to keep is based on the target total number
of bases, which is coverage * genome size. So if coverage is 10.7 and genome size is
100, previously our target number of bases would have been 1000, whereas now, it would
be 1070.

### Changed

- `--coverage` is now treated as a `f32` instead of being converted immediately to an
  integer [#19][19].
- Updated `rust-bio` to version 0.31.0. This means `rasusa` now handles wrapped fastq
  files.
- Preallocate fastx records instead of using iterator. Gives marginal speedup.
- Added `bash` to the docker image b47a8b75943098bdd845b7758cf2eab01ef5a3d8

## [0.2.0]

### Added

- Support paired Illumina [#15](https://github.com/mbhall88/rasusa/issues/15)

[0.2.0]: https://github.com/mbhall88/rasusa/releases/tag/0.2.0
[0.3.0]: https://github.com/mbhall88/rasusa/releases/tag/0.3.0
[0.4.0]: https://github.com/mbhall88/rasusa/releases/tag/0.4.0
[0.4.1]: https://github.com/mbhall88/rasusa/releases/tag/0.4.1
[0.4.2]: https://github.com/mbhall88/rasusa/releases/tag/0.4.2
[0.5.0]: https://github.com/mbhall88/rasusa/releases/tag/0.5.0
[0.6.0]: https://github.com/mbhall88/rasusa/releases/tag/0.6.0
[0.6.1]: https://github.com/mbhall88/rasusa/releases/tag/0.6.1
[19]: https://github.com/mbhall88/rasusa/issues/19
[22]: https://github.com/mbhall88/rasusa/issues/22
[27]: https://github.com/mbhall88/rasusa/issues/27
[28]: https://github.com/mbhall88/rasusa/pull/28
[30]: https://github.com/mbhall88/rasusa/issues/30
[31]: https://github.com/mbhall88/rasusa/issues/31
[35]: https://github.com/mbhall88/rasusa/issues/35
[36]: https://github.com/mbhall88/rasusa/issues/36
[benchmark]: https://github.com/mbhall88/rasusa#benchmark
[error-blog]: https://nick.groenen.me/posts/rust-error-handling/
[unreleased]: https://github.com/mbhall88/rasusa/compare/0.6.1...HEAD

![rasusa](img/logo.png)

[![Build Status](https://github.com/mbhall88/rasusa/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/mbhall88/rasusa/actions/workflows/rust-ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/rasusa/branch/master/graph/badge.svg)](https://codecov.io/gh/mbhall88/rasusa)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![github release version](https://img.shields.io/github/v/release/mbhall88/rasusa)](https://github.com/mbhall88/rasusa/releases)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03941/status.svg)](https://doi.org/10.21105/joss.03941)

**Ra**ndomly **su**b**sa**mple sequencing reads to a specified coverage.

[TOC]: #
## Table of Contents
- [Motivation](#motivation)
- [Install](#install)
  - [`cargo`](#cargo)
  - [`conda`](#conda)
  - [Container](#container)
  - [`homebrew`](#homebrew)
  - [Release binaries](#release-binaries)
  - [Build locally](#build-locally)
- [Usage](#usage)
  - [Basic usage](#basic-usage)
  - [Required parameters](#required-parameters)
  - [Optional parameters](#optional-parameters)
  - [Full usage](#full-usage)
  - [Snakemake](#snakemake)
- [Benchmark](#benchmark)
  - [Single long read input](#single-long-read-input)
  - [Paired-end input](#paired-end-input)
- [Contributing](#contributing)
- [Citing](#citing)
  - [Bibtex](#bibtex)

## Motivation

I couldn't find a tool for subsampling reads that met my requirements. All the
strategies I could find fell short as they either just wanted a number or percentage of
reads to subsample to or, if they did subsample to a coverage, they assume all reads are
the same size (i.e Illumina). As I mostly work with long-read data this posed a problem
if I wanted to subsample a file to certain coverage, as length of reads was never taken
into account. `rasusa` addresses this shortcoming.

A workaround I had been using for a while was using [`filtlong`][filtlong]. It was
simple enough, I just figure out the number of bases I need to achieve a (theoretical)
coverage for my sample. Say I have a fastq from an _E. coli_ sample with 5 million reads
and I want to subset it to 50x coverage. I just need to multiply the expected size of
the sample's genome, 4.6 million base pairs, by the coverage I want and I have my target
bases - 230 million base pairs. In `filtlong`, I can do the following

```sh
target=230000000
filtlong --target_bases "$target" reads.fq > reads.50x.fq
```

However, this is technically not the intended function of `filtlong`; it's a quality
filtering tool. What you get in the end is a subset of the ["highest scoring"][score]
reads at a (theoretical) coverage of 50x. Depending on your circumstances, this might be
what you want. However, you bias yourself towards the best/longest reads in the dataset
\- not a fair representation of your dataset as a whole. There is also the possibility
of favouring regions of the genome that produce longer/higher quality reads. [De Maio
_et al._][mgen-ref] even found that by randomly subsampling nanopore reads you achieve
_better_ genome assemblies than if you had filtered.

So, depending on your circumstances, an unbiased subsample of your reads might be what
you need. And if this is the case, `rasusa` has you covered.

## Install

Some of these installation options require the [`rust` toolchain][rust], which is
extremely easy to set up. However, if you do not wish to install `rust` then there are a
number of options available.

### `cargo`

[![Crates.io](https://img.shields.io/crates/v/rasusa.svg)](https://crates.io/crates/rasusa)

Prerequisite: [`rust` toolchain][rust] (min. v1.53.0)

```sh
cargo install rasusa
```

### `conda`

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/rasusa)](https://anaconda.org/bioconda/rasusa)
[![bioconda version](https://anaconda.org/bioconda/rasusa/badges/platforms.svg)](https://anaconda.org/bioconda/rasusa)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```sh
conda install rasusa
```

Thank you to Devon Ryan ([@dpryan79][dpryan79]) for [help debugging the bioconda
recipe][pr-help].

### Container

Docker images are hosted at [quay.io]. For versions 0.3.0 and earlier, the images were
hosted on [Dockerhub][dockerhub].

#### `singularity`

Prerequisite: [`singularity`][singularity]

```sh
URI="docker://quay.io/mbhall88/rasusa"
singularity exec "$URI" rasusa --help
```

The above will use the latest version. If you want to specify a version then use a
[tag][quay.io] (or commit) like so.

```sh
VERSION="0.6.1"
URI="docker://quay.io/mbhall88/rasusa:${VERSION}"
```

#### `docker`

[![Docker Repository on Quay](https://quay.io/repository/mbhall88/rasusa/status "Docker Repository on Quay")](https://quay.io/repository/mbhall88/rasusa)

Prerequisite: [`docker`][docker]

```sh
docker pull quay.io/mbhall88/rasusa
docker run quay.io/mbhall88/rasusa --help
```

You can find all the available tags on the [quay.io repository][quay.io]. Note: versions
prior to 0.4.0 were housed on [Docker Hub](https://hub.docker.com/r/mbhall88/rasusa).


### `homebrew`

Prerequisite: [`homebrew`][homebrew]

The `homebrew` installation is done via the [homebrew-bio tap][brew-tap].

```sh
brew tap brewsci/bio
brew install rasusa
```

or

```sh
brew install brewsci/bio/rasusa
```

### Release binaries

**tl;dr**: Run the following snippet to download the binary for your system to the
current directory and show the help menu.

```shell
VERSION="0.6.1"  # change accordingly
OS=$(uname -s)                                                                                                       
if [ "$OS" = "Linux" ]; then                                                                                         
    triple="x86_64-unknown-linux-musl"                                                                              
elif [ "$OS" = "Darwin" ]; then                                                                                        
    triple="x86_64-apple-darwin"                                                         
else                                                      
    echo "ERROR: $OS not a recognised operating system"
fi              
if [ -n "$triple" ]; then   
    URL="https://github.com/mbhall88/rasusa/releases/download/${VERSION}/rasusa-${VERSION}-${triple}.tar.gz"
    wget "$URL" -O - | tar -xzf -
    ./rasusa --help             
fi
```

These binaries _do not_ require that you have the `rust` toolchain installed.

Currently, there are two pre-compiled binaries available:
- Linux kernel `x86_64-unknown-linux-musl` (works on most Linux distributions I tested)
- OSX kernel `x86_64-apple-darwin` (works for any post-2007 Mac)

If these binaries do not work on your system please raise an issue and I will
potentially add some additional [target triples][triples].

### Build locally

Prerequisite: [`rust` toolchain][rust]

```sh
git clone https://github.com/mbhall88/rasusa.git
cd rasusa
cargo build --release
target/release/rasusa --help
# if you want to check everything is working ok
cargo test --all
```


## Usage

### Basic usage

```
rasusa --input in.fq --coverage 30 --genome-size 4.6mb
```

The above command will output the subsampled file to `stdout`.

Or, if you have paired Illumina

```
rasusa -i r1.fq -i r2.fq --coverage 30 --genome-size 4g -o out.r1.fq -o out.r2.fq
```

For more details on the above options, and additional options, see below.

### Required parameters

There are three required options to run `rasusa`.

#### Input

##### `-i`, `--input`

This option specifies the file(s) containing the reads you would like to subsample. The
file(s) must be valid fasta or fastq format and can be compressed (with a tool such as
`gzip`).  
Illumina paired files can be passed in two ways.
1. Using `--input` twice `-i r1.fq -i r2.fq`
2. Using `--input` once, but passing both files immediately after `-i r1.fq r2.fq`

> Bash wizard tip 🧙: Let globs do the work for you `-i r*.fq`

#### Coverage

##### `-c`, `--coverage`

> Not required if [`--bases`](#target-number-of-bases) is present

This option is used to determine the minimum coverage to subsample the reads to. It can
be specified as an integer (100), a decimal/float (100.7), or either of the previous
suffixed with an 'x' (100x).

_**Note**: Due to the method for determining how many bases are required to achieve the
desired coverage, the actual coverage, in the end, could be slightly higher than
requested. For example, if the last included read is very long. The log messages should
inform you of the actual coverage in the end._

#### Genome size

##### `-g`, `--genome-size`

> Not required if [`--bases`](#target-number-of-bases) is present

The genome size of the input is also required. It is used to determine how many bases
are necessary to achieve the desired coverage. This can, of course, be as precise or
rough as you like.  
Genome size can be passed in many ways. As a plain old integer (1600), or with a metric
suffix (1.6kb). All metric suffixes can have an optional 'b' suffix and be lower, upper,
or mixed case. So 'Kb', 'kb' and 'k' would all be inferred as 'kilo'. Valid metric
suffixes include:

- Base (b) - multiplies by 1
- Kilo (k) - multiplies by 1,000
- Mega (m) - multiplies by 1,000,000
- Giga (g) - multiplies by 1,000,000,000
- Tera (t) - multiplies by 1,000,000,000,000

Alternatively, a [FASTA/Q index file][faidx] can be given and the genome size will be
set to the sum of all reference sequences in it.

[faidx]: https://www.htslib.org/doc/faidx.html

### Optional parameters

#### Output

##### `-o`, `--output`

NOTE: This parameter is required if passing paired Illumina data.

By default, `rasusa` will output the subsampled file to `stdout` (if one file is given).
If you would prefer to specify an output file path, then use this option.

Output for Illumina paired files can be specified in the same manner as
[`--input`](#input)
1. Using `--output` twice `-o out.r1.fq -o out.r2.fq`
2. Using `--output` once, but passing both files immediately after `-o out.r1.fq
   out.r2.fq`

The ordering of the output files is assumed to be the same as the input.  
_Note: The output will always be in the same format as the input. You cannot pass fastq
as input and ask for fasta as output._

`rasusa` will also attempt to automatically infer whether comression of the output
file(s) is required. It does this by detecting any of the supported extensions:
- `.gz`: will compress the output with [`gzip`][gzip]
- `.bz` or `.bz2`: will compress the output with [`bzip2`][bzip]
- `.lzma`: will compress the output with the [`xz`][xz] LZMA algorithm

[gzip]: http://www.gzip.org/
[bzip]: https://sourceware.org/bzip2/
[xz]: https://tukaani.org/xz/

#### Output compression format

##### `-O`, `--output-type`

Use this option to manually set the compression algoritm to use for the output file(s).
It will override any format automatically detected from the output path.

Valid options are:
- `g`: [`gzip`][gzip]
- `b`: [`bzip2`][bzip]
- `l`: [`xz`][xz] LZMA algorithm
- `u`: no compression

*Note: these options are case insensitive.*

#### Compresion level

##### `-l`, `--compress-level`

Compression level to use if compressing the output. 1 is for fastest/least compression
and 9 is for slowest/best. By default this is set to 6, which is also the default for
most compression programs.

#### Target number of bases

##### `-b`, `--bases`

Explicitly set the number of bases required in the subsample. This option takes the
number in the same format as [genome size](#genome-size).

*Note: if this option is given, genome size and coverage are not required, or ignored if
they are provided.*

#### Random seed

##### `-s`, `--seed`

This option allows you to specify the [random seed][seed] used by the random subsampler.
By explicitly setting this parameter, you make the subsample for the input reproducible.
The seed is an integer, and by default it is not set, meaning the operating system will
seed the random subsampler. You should only pass this parameter if you are likely to
want to subsample the same input file again in the future and want the same subset of
reads.

#### Verbosity

##### `-v`

Adding this optional flag will make the logging more verbose. By default, logging will
produce messages considered "info" or above (see [here][log-lvl] for more details). If
verbosity is switched on, you will additionally get "debug" level logging messages.

### Full usage

```text
$ rasusa --help

rasusa 0.6.1
Randomly subsample reads to a specified coverage.

USAGE:
    rasusa [FLAGS] [OPTIONS] --bases <bases> --coverage <FLOAT> --genome-size <size|faidx> --input <input>...

FLAGS:
    -h, --help
            Prints help information

    -V, --version
            Prints version information

    -v
            Switch on verbosity.


OPTIONS:
    -b, --bases <bases>
            Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB

            If this option is given, --coverage and --genome-size are ignored
    -l, --compress-level <1-9>
            Compression level to use if compressing output [default: 6]

    -c, --coverage <FLOAT>
            The desired coverage to sub-sample the reads to

            If --bases is not provided, this option and --genome-size are required
    -g, --genome-size <size|faidx>
            Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB

            Alternatively, a path to a FASTA/Q index file can be provided and the genome size will be set to the sum of
            all reference sequences.

            If --bases is not provided, this option and --coverage are required
    -i, --input <input>...
            The fast{a,q} file(s) to subsample.

            For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two files consecutively
            `-i r1.fq r2.fq`.
    -o, --output <output>...
            Output filepath(s); stdout if not present.

            For paired Illumina you may either pass this flag twice `-o o1.fq -o o2.fq` or give two files consecutively
            `-o o1.fq o2.fq`. NOTE: The order of the pairs is assumed to be the same as that given for --input. This
            option is required for paired input.
    -O, --output-type <u|b|g|l>
            u: uncompressed; b: Bzip2; g: Gzip; l: Lzma

            Rasusa will attempt to infer the output compression format automatically from the filename extension. This
            option is used to override that. If writing to stdout, the default is uncompressed
    -s, --seed <INT>
            Random seed to use.
```

### Snakemake

If you want to use `rasusa` in a [`snakemake`][snakemake] pipeline, it is advised to use
the [wrapper][wrapper].

```py
rule subsample:
    input:
        r1="{sample}.r1.fq",
        r2="{sample}.r2.fq",
    output:
        r1="{sample}.subsampled.r1.fq",
        r2="{sample}.subsampled.r2.fq",
    params:
        options="--seed 15",  # optional
        genome_size="3mb",  # required
        coverage=20,  # required
    log:
        "logs/subsample/{sample}.log",
    wrapper:
        "0.70.0/bio/rasusa"
```

*See the [latest wrapper][wrapper] documentation for the most up-to-date version
number.*

## Benchmark

> “Time flies like an arrow; fruit flies like a banana.”  
> ― Anthony G. Oettinger

The real question is: will `rasusa` just needlessly eat away at your precious time on
earth?

To do this benchmark, I am going to use [hyperfine][hyperfine].

The data I used comes from

> [Bainomugisa, Arnold, et al. "A complete high-quality MinION nanopore assembly of an
> extensively drug-resistant Mycobacterium tuberculosis Beijing lineage strain
> identifies novel variation in repetitive PE/PPE gene regions." Microbial genomics 4.7
> (2018).][1]

### Single long read input

Download and rename the fastq

```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR649/008/SRR6490088/SRR6490088_1.fastq.gz"
wget "$URL" -O - | gzip -d -c > tb.fq
```

The file size is 2.9G, and it has 379,547 reads.  
We benchmark against `filtlong` using the same strategy outlined in
[Motivation](#motivation).

```shell
TB_GENOME_SIZE=4411532
COVG=50
TARGET_BASES=$(( TB_GENOME_SIZE * COVG ))
FILTLONG_CMD="filtlong --target_bases $TARGET_BASES tb.fq"
RASUSA_CMD="rasusa -i tb.fq -c $COVG -g $TB_GENOME_SIZE -s 1"
hyperfine --warmup 3 --runs 10 --export-markdown results-single.md \
     "$FILTLONG_CMD" "$RASUSA_CMD" 
```

#### Results

| Command                                   |       Mean [s] | Min [s] | Max [s] |     Relative |
|:------------------------------------------|---------------:|--------:|--------:|-------------:|
| `filtlong --target_bases 220576600 tb.fq` | 21.685 ± 0.055 |  21.622 |  21.787 | 21.77 ± 0.29 |
| `rasusa -i tb.fq -c 50 -g 4411532 -s 1`   | 0.996 ±  0.013 |   0.983 |   1.023 |         1.00 |

**Summary**: `rasusa` ran 21.77 ± 0.29 times faster than `filtlong`.

### Paired-end input

Download and then deinterleave the fastq with [`pyfastaq`][pyfastaq]

```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6488968/SRR6488968.fastq.gz"
wget "$URL" -O - | gzip -d -c - | fastaq deinterleave - r1.fq r2.fq
```

Each file's size is 179M and has 283,590 reads.  
For this benchmark, we will use [`seqtk`][seqtk]. As `seqtk` requires a fixed number of
reads to subsample to, I ran `rasusa` initially and got the number of reads it was using
for its subsample. We will also test `seqtk`'s 2-pass mode as this is analogous to
`rasusa`.  
We use a lower coverage for the Illumina as the samples only have ~38x coverage.

```shell
TB_GENOME_SIZE=4411532
COVG=20
NUM_READS=147052
SEQTK_CMD_1="seqtk sample -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
SEQTK_CMD_2="seqtk sample -2 -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
RASUSA_CMD="rasusa -i r1.fq r2.fq -c $COVG -g $TB_GENOME_SIZE -s 1 -o /tmp/r1.fq -o /tmp/r2.fq"
hyperfine --warmup 5 --runs 20 --export-markdown results-paired.md \
     "$SEQTK_CMD_1" "$SEQTK_CMD_2" "$RASUSA_CMD"
```

#### Results

| Command                                                                                           |    Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------------------------------------------------------------------------------------|-------------:|---------:|---------:|------------:|
| `seqtk sample -s 1 r1.fq 147052 > /tmp/r1.fq; seqtk sample -s 1 r2.fq 147052 > /tmp/r2.fq;`       | 693.8 ± 17.4 |    678.8 |    780.5 | 1.21 ± 0.19 |
| `seqtk sample -2 -s 1 r1.fq 147052 > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq 147052 > /tmp/r2.fq;` | 739.5 ± 28.4 |    695.6 |    811.8 | 1.29 ± 0.20 |
| `./rasusa -i r1.fq r2.fq -c 20 -g 4411532 -s 1 -o /tmp/r1.fq -o /tmp/r2.fq`                       | 574.7 ± 88.1 |    432.3 |    784.4 |        1.00 |

**Summary**: `rasusa` ran 1.21 times faster than `seqtk` (1-pass) and 1.29 times faster
than `seqtk` (2-pass)

So, `rasusa` is faster than `seqtk` but doesn't require a fixed number of reads -
allowing you to avoid doing maths to determine how many reads you need to downsample to
a specific coverage. 🤓

## Contributing

If you would like to help improve `rasusa` you are very welcome!

For changes to be accepted, they must pass the CI and coverage checks. These include:

- Code is formatted with `rustfmt`. This can be done by running `cargo fmt` in the
  project directory.
- There are no compiler errors/warnings. You can check this by running `cargo clippy
  --all-features --all-targets -- -D warnings`
- Code coverage has not reduced. If you want to check coverage before pushing changes, I
  use [`kcov`][kcov].

## Citing

If you use `rasusa` in your research, it would be very much appreciated if you could
cite it.

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03941/status.svg)](https://doi.org/10.21105/joss.03941)

> Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, https://doi.org/10.21105/joss.03941

### Bibtex

```Bibtex
@article{Hall2022,
  doi = {10.21105/joss.03941},
  url = {https://doi.org/10.21105/joss.03941},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {3941},
  author = {Michael B. Hall},
  title = {Rasusa: Randomly subsample sequencing reads to a specified coverage},
  journal = {Journal of Open Source Software}
}

```

[1]: https://doi.org/10.1099/mgen.0.000188
[brew-tap]: https://github.com/brewsci/homebrew-bio
[channels]: https://bioconda.github.io/user/install.html#set-up-channels
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[docker]: https://docs.docker.com/v17.12/install/
[dockerhub]: https://hub.docker.com/r/mbhall88/rasusa
[dpryan79]: https://github.com/dpryan79
[filtlong]: https://github.com/rrwick/Filtlong
[homebrew]: https://docs.brew.sh/Installation
[hyperfine]: https://github.com/sharkdp/hyperfine
[kcov]: https://github.com/SimonKagstrom/kcov
[log-lvl]: https://docs.rs/log/0.4.6/log/enum.Level.html#variants
[mgen-ref]: https://doi.org/10.1099/mgen.0.000294
[pr-help]: https://github.com/bioconda/bioconda-recipes/pull/18690
[pyfastaq]: https://github.com/sanger-pathogens/Fastaq
[quay.io]: https://quay.io/repository/mbhall88/rasusa
[rust]: https://www.rust-lang.org/tools/install
[score]: https://github.com/rrwick/Filtlong#read-scoring
[seed]: https://en.wikipedia.org/wiki/Random_seed
[seqtk]: https://github.com/lh3/seqtk
[singularity]: https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[triples]: https://clang.llvm.org/docs/CrossCompilation.html#target-triple
[wrapper]: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/rasusa.html


---
title: 'Rasusa: Randomly subsample sequencing reads to a specified coverage'
tags:
  - Rust
  - bioinformatics
  - genomics
  - fastq
  - fasta
  - subsampling
  - random
authors:
  - name: Michael B. Hall
    orcid: 0000-0003-3683-6208
    affiliation: 1
affiliations:
 - name: European Molecular Biology Laboratory, European Bioinformatics Institute EMBL-EBI, Hinxton, UK
   index: 1
date: 18 October 2021
bibliography: paper.bib
---

# Summary

A fundamental requirement for many applications in genomics is the sequencing of genetic
material (DNA/RNA). Different sequencing technologies exist, but all aim to accurately
reproduce the sequence of nucleotides (the individual units of DNA and RNA) in the
genetic material under investigation. The result of such efforts is a text file
containing the individual fragments of genetic material - termed "reads" - represented
as strings of letters (A, C, G, and T/U).

The amount of data in one of these read files depends on how much genetic material was
present and how long the sequencing device was operated. Read depth (coverage) is a
measure of the volume of genetic data contained in a read file. For example, coverage of
5x indicates that, on average, each nucleotide in the original genetic material is represented five
times in the read file.

Many of the computational methods employed in genomics are affected by coverage;
counterintuitively, more is not always better. For example, because sequencing devices
are not perfect, reads inevitably contain errors. As such, higher coverage increases the
number of errors and potentially makes them look like alternative sequences.
Furthermore, for some applications, too much coverage can cause a degradation in
computational performance via increased runtimes or memory usage.

We present Rasusa, a software program that randomly subsamples a given read file to a
specified coverage. Rasusa is written in the Rust programming language and is much
faster than current solutions for subsampling read files. In addition, it provides an
ergonomic command-line interface and allows users to specify a desired coverage or a
target number of nucleotides.


# Statement of need

Read subsampling is a useful mechanism for creating artificial datasets, allowing
exploration of a computational method's performance as data becomes more scarce. In
addition, the coverage of a sample can have a significant impact on a variety of
computational methods, such as RNA-seq [@Baccarella2018], taxonomic classification
[@Gweon2019], antimicrobial resistance detection [@Gweon2019], and genome assembly
[@Maio2019] - to name a few.

There is limited available software for subsampling read files. Assumably, most
researchers use custom scripts for this purpose. However, two existing programs for
subsampling are Filtlong [@filtlong] and Seqtk [@seqtk]. Unfortunately, neither of these
tools provides subsampling to a specified coverage "out of the box".

Filtlong is technically a filtering tool, not a subsampling one. It scores each read
based on its length and quality and outputs the highest-scoring subset. Additionally,
minimum and maximum read lengths can be specified, along with the size of the subset
required. Ultimately, the subset produced by Filtlong is not necessarily representative
of the original reads but is biased towards those with the greatest length or quality.
While this may sound like a good thing, in some applications, such as genome assembly,
it has been shown that a random subsample produces superior results to a filtered subset
[@Maio2019].

Seqtk does do random subsampling via the sample subcommand. However, the only option
available is to specify the number of reads required. Thus, it is up to the user to
determine the number of reads required to reach the desired coverage. While this serves
for Illumina sequencing data, which generally have uniform(ish) read lengths, it does
not work for other modalities like PacBio and Nanopore, where read lengths vary
significantly.

Rasusa provides a random subsample of a read file (FASTA or FASTQ), with two ways of
specifying the size of the subset. One method takes a genome size and the desired
coverage, while the other takes a target number of bases (nucleotides). In the genome
size and coverage option, we multiply the genome size by the coverage to obtain the
target number of bases for the subset. As such, the resulting read file will have, on 
average, the amount of coverage requested. In addition, Rasusa allows setting a random seed
to allow reproducible subsampling. Other features include user control over whether the
output is compressed and specifying the compression algorithm and level.

Rasusa is 21 and 1.2 times faster than Filtlong and Seqtk, respectively.

# Availability

Rasusa is open-source and available under an MIT license at
https://github.com/mbhall88/rasusa.

# Acknowledgements

We acknowledge contributions from Pierre Marijon and suggestions from Zamin Iqbal. In
addition, MBH is funded by the EMBL International PhD Programme.

# References

