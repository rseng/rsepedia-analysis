<!-- markdownlint-disable MD013 MD041 -->

![Logo](https://cbg-ethz.github.io/V-pipe/img/logo.svg)

[![bio.tools](https://img.shields.io/badge/bio-tools-blue.svg)](https://bio.tools/V-Pipe)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.8.1-blue.svg)](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe)
[![Deploy Docker image](https://github.com/cbg-ethz/V-pipe/actions/workflows/deploy-docker.yaml/badge.svg)](https://github.com/cbg-ethz/V-pipe/pkgs/container/v-pipe)
[![Tests](https://github.com/cbg-ethz/V-pipe/actions/workflows/run_regression_tests.yaml/badge.svg)](https://github.com/cbg-ethz/V-pipe/actions/workflows/run_regression_tests.yaml)
[![Mega-Linter](https://github.com/cbg-ethz/V-pipe/actions/workflows/mega-linter.yml/badge.svg)](https://github.com/cbg-ethz/V-pipe/actions/workflows/mega-linter.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

V-pipe is a workflow designed for the analysis of next generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format (e.g., consensus sequences, SNV calls, local/global haplotypes).
V-pipe is written using the Snakemake workflow management system.

## Usage

Different ways of initializing V-pipe are presented below. We strongly encourage you to deploy it [using the quick install script](#using-quick-install-script), as this is our preferred method.

To configure V-pipe refer to the documentation present in [config/README.md](config/README.md).

V-pipe expects the input samples to be organized in a [two-level](config/README.md#samples) directory hierarchy,
and the sequencing reads must be provided in a sub-folder named `raw_data`. Further details can be found on the [website](https://cbg-ethz.github.io/V-pipe/usage/).
Check the utils subdirectory for [mass-importers tools](utils/README.md#samples-mass-importers) that can assist you in generating this hierarchy.

We provide [virus-specific base configuration files](config/README.md#virus-base-config) which contain handy defaults for, e.g., HIV and SARS-CoV-2. Set the virus in the general section of the configuration file:

```yaml
general:
  virus_base_config: hiv
```

Also see [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) to learn more about the command-line options available when executing the workflow.

### Using quick install script

To deploy V-pipe, use the [installation script](utils/README.md#quick-installer) with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
./quick_install.sh -w work
```

This script will download and install miniconda, checkout the V-pipe git repository (use `-b` to specify which branch/tag) and setup a work directory (specified with `-w`) with an executable script that will execute the workflow:

```bash
cd work
# edit config.yaml and provide samples/ directory
./vpipe --jobs 4 --printshellcmds --dry-run
```

### Using Docker

Note: the [docker image](https://github.com/cbg-ethz/V-pipe/pkgs/container/v-pipe) is only setup with components to run the workflow for HIV and SARS-CoV-2 virus base configurations.
Using V-pipe with other viruses or configurations might require internet connectivity for additional software components.

Create `config.yaml` or `vpipe.config` and then populate the `samples/` directory.
For example, the following config file could be used:

```yaml
general:
  virus_base_config: hiv

output:
  snv: true
  local: true
  global: false
  visualization: true
  QA: true
```

Then execute:

```bash
docker run --rm -it -v $PWD:/work ghcr.io/cbg-ethz/v-pipe:master --jobs 4 --printshellcmds --dry-run
```

### Using Snakedeploy

First install [mamba](https://github.com/conda-forge/miniforge#mambaforge), then create and activate an environment with Snakemake and Snakedeploy:

```bash
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
conda activate snakemake
```

Snakemake's [official workflow installer Snakedeploy](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe) can now be used:

```bash
snakedeploy deploy-workflow https://github.com/cbg-ethz/V-pipe --tag master .
# edit config/config.yaml and provide samples/ directory
snakemake --use-conda --jobs 4 --printshellcmds --dry-run
```

## Dependencies

- **[Conda](https://conda.io/docs/index.html)**

  Conda is a cross-platform package management system and an environment manager application. Snakemake uses mamba as a package manager.

- **[Snakemake](https://snakemake.readthedocs.io/)**

  Snakemake is the central workflow and dependency manager of V-pipe. It determines the order in which individual tools are invoked and checks that programs do not exit unexpectedly.

- **[VICUNA](https://www.broadinstitute.org/viral-genomics/vicuna)**

  VICUNA is a *de novo* assembly software designed for populations with high mutation rates. It is used to build an initial reference for mapping reads with ngshmmalign aligner when a `references/cohort_consensus.fasta` file is not provided. Further details can be found in the [wiki](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#input-files) pages.

### Computational tools

Other dependencies are managed by using isolated conda environments per rule, and below we list some of the computational tools integrated in V-pipe:

- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**

  FastQC gives an overview of the raw sequencing data. Flowcells that have been overloaded or otherwise fail during sequencing can easily be determined with FastQC.

- **[PRINSEQ](http://prinseq.sourceforge.net/)**

  Trimming and clipping of reads is performed by PRINSEQ. It is currently the most versatile raw read processor with many customization options.

- **[ngshmmalign](https://github.com/cbg-ethz/ngshmmalign)**

  We perform the alignment of the curated NGS data using our custom ngshmmalign that takes structural variants into account. It produces multiple consensus sequences that include either majority bases or ambiguous bases.

- **[bwa](https://github.com/lh3/bwa)**

  In order to detect specific cross-contaminations with other probes, the Burrows-Wheeler aligner is used. It quickly yields estimates for foreign genomic material in an experiment.
  Additionally, It can be used as an alternative aligner to ngshmmalign.

- **[MAFFT](http://mafft.cbrc.jp/alignment/software/)**

  To standardise multiple samples to the same reference genome (say HXB2 for HIV-1), the multiple sequence aligner MAFFT is employed. The multiple sequence alignment helps in determining regions of low conservation and thus makes standardisation of alignments more robust.

- **[Samtools and bcftools](https://www.htslib.org/)**

  The Swiss Army knife of alignment postprocessing and diagnostics. bcftools is also used to generate consensus sequence with indels.

- **[SmallGenomeUtilities](https://github.com/cbg-ethz/smallgenomeutilities)**

  We perform genomic liftovers to standardised reference genomes using our in-house developed python library of utilities for rewriting alignments.

- **[ShoRAH](https://github.com/cbg-ethz/shorah)**

  ShoRAh performs SNV calling and local haplotype reconstruction by using bayesian clustering.

- **[LoFreq](https://csb5.github.io/lofreq/)**

  LoFreq (version 2) is SNVs and indels caller from next-generation sequencing data, and can be used as an alternative engine for SNV calling.

- **[SAVAGE](https://bitbucket.org/jbaaijens/savage) and [Haploclique](https://github.com/cbg-ethz/haploclique)**

  We use HaploClique or SAVAGE to perform global haplotype reconstruction for heterogeneous viral populations by using an overlap graph.

## Citation

If you use this software in your research, please cite:

Posada-Céspedes S., Seifert D., Topolsky I., Jablonski K.P., Metzner K.J., and Beerenwinkel N. 2021.
"V-pipe: a computational pipeline for assessing viral genetic diversity from high-throughput sequencing data."
_Bioinformatics_, January. doi:[10.1093/bioinformatics/btab015](https://doi.org/10.1093/bioinformatics/btab015).

## Contributions

- [Ivan Topolsky\* ![orcid]](https://orcid.org/0000-0002-7561-0810), [![github]](https://github.com/dryak)
- [Kim Philipp Jablonski ![orcid]](https://orcid.org/0000-0002-4166-4343), [![github]](https://github.com/kpj)
- [Lara Fuhrmann ![orcid]](https://orcid.org/0000-0001-6405-0654), [![github]](https://github.com/LaraFuhrmann)
- [Uwe Schmitt ![orcid]](https://orcid.org/0000-0002-4658-0616), [![github]](https://github.com/uweschmitt)
- [Monica Dragan ![orcid]](https://orcid.org/0000-0002-7719-5892), [![github]](https://github.com/monicadragan)
- [Susana Posada Céspedes ![orcid]](https://orcid.org/0000-0002-7459-8186), [![github]](https://github.com/sposadac)
- [David Seifert ![orcid]](https://orcid.org/0000-0003-4739-5110), [![github]](https://github.com/SoapZA)
- Tobias Marschall
- [Niko Beerenwinkel\*\* ![orcid]](https://orcid.org/0000-0002-0573-6119)

\* software maintainer ;
\** group leader

[github]: https://cbg-ethz.github.io/V-pipe/img/mark-github.svg
[orcid]: https://cbg-ethz.github.io/V-pipe/img/ORCIDiD_iconvector.svg

## Contact

We encourage users to use the [issue tracker](https://github.com/cbg-ethz/V-pipe/issues). For further enquiries, you can also contact the V-pipe Dev Team <v-pipe@bsse.ethz.ch>.
<!-- markdownlint-disable MD013 -->

# Configuring V-pipe

In order to start using V-pipe, you need to provide three things:

 1. Samples in a specific directory structure
 2. _(optional)_ TSV file listing the samples
 3. Configuration file

The utils subdirectory provides [tools](../utils/README.md#samples-mass-importers) that can assist in importing samples files and structuring them.


## Configuration file

The V-pipe workflow is customized using a structured configuration file called `config.yaml`, `config.json` or, for backward compatibility, `vpipe.config` (INI-like format).

This configuration file is a text file written using a basic structure composed of sections, properties and values. When using [YAML](https://yaml.org/spec/1.0/#id2564813) or [JSON](https://www.json.org/json-en.html) format use these languages associative array/dictionaries in two levels for sections and properties. When using the older [INI format](https://docs.python.org/3/library/configparser.html), sections are expected in squared brackets, and properties are followed by corresponding values.

Further more, it is possible to specify additional options on the command line using Snakemake's `--configfile` to pass additional YAML/JSON configuration files, and/or using Snakemake's `--config` to pass sections and properties in a [YAML Flow style](https://yaml.org/spec/1.2.0/#Flow)/JSON syntax.

Here is an **example** of `config.yaml`:

```yaml
general:
  virus_base_config: hiv

input:
  datadir: samples
  samples_file: config/samples.tsv

output:
  datadir: results
  snv: true
  local: true
  global: false
  visualization: true
  QA: true
```

At minimum, a valid configuration **MUST** provide a reference sequence against which to align the short reads from the raw data. This can be done in several ways:

- by using a [_virus base config_](#virus-base-config) that will provide default presets for specific viruses
- by directly passing a reference .fasta file in the section _input_ -> property _reference_ that will override the default

### virus base config

We provide virus-specific base configuration files which contain handy defaults for some viruses.

Currently, the following _virus base config_ are available:

- [hiv](hiv.yaml): provides HXB2 as a reference sequence for HIV, and sets the default aligner to _ngshmmalign_.
- [sars-cov-2](sars-cov-2.yaml): provides NC\_045512.2 as a reference sequence for SARS-CoV-2, sets the default aligner to _bwa_ and sets the variant calling to be done against the reference instead of the cohort's consensus.

### configuration manual

More information about all the available configuration options and an exhaustive list can be found in [config.html](config.html)
or [online](https://htmlpreview.github.io/?https://github.com/cbg-ethz/V-pipe/blob/master/config/config.html).

### legacy V-pipe 1.xx/2.xx users

If you want to re-use your old configuration
from a [legacy V-pipe v1.x/2.x installation](https://github.com/cbg-ethz/V-pipe/wiki/options)
or [sars-cov2 branch](https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/#running-v-pipe)
it is possible, if you keep in mind the following caveats:

- The older INI-like syntax is still supported for a `vpipe.config` configuration file.
  - This configuration will be overridden by `config.yaml` or `config.json`,
    you might want to delete those files from your working directory if you are not using them.
- V-pipe starting from version 2.99.1 follows the [Standardized usage](https://snakemake.github.io/snakemake-workflow-catalog/?rules=true) rules of the
  [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe)
  - This defines a newer [directory structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility)
    - samples TSV table is now expected to be in `config/samples.tsv`
      (use the section _input_ ->  property _samples_file_ to override).
    - the per sample output isn't written in the same `samples/` directory as the input anymore, but in a separate directory called `results/`
      (use the section _output_ -> property _datadir_ to override).
    - the cohort-wide output isn't written in a different `variants/` directory anymore, but at at the base of the _output datadir_ - i.e by default in `results/`
      (use the section _output_ -> property _cohortdir_ to specify a different path **relative to the output datadir**).
  - Add the following sections and properties to your `vpipe.config` configuration file to **bring back the legacy behaviour**:

```ini
[input]
datadir=samples
samples_file=samples.tsv

[output]
datadir=samples
cohortdir=../variants
```

As of version 2.99.1, only the analysis of viral sequencing data has been
[extensively tested](https://github.com/cbg-ethz/V-pipe/actions/workflows/run_regression_tests.yaml)
and is guaranteed stable.
For other more advanced functionality you might want to wait until a future release.

## samples tsv

File containing sample unique identifiers and dates as tab-separated values.

**Example:** here, we have two samples from patient 1 and one sample from patient 2:

```tsv
patient1	20100113
patient1	20110202
patient2	20081130
```

By default, V-pipe searches for a file named `config/samples.tsv`, if this file does not exist, a list of samples is built by searching the contents of the input datadir.

### read-lenght

The samples' read-length is used for critical steps of the pipeline (e.g.: quality filtering). Different possibilities are available to set its value:

- by default, V-pipe expects a read-length of 250bp
- this default can be globally overridden in the configuration file in section _input_ -> property _read_length_
  ```yaml
  input:
    read_length: 150
  ```
- the samples TSV file can contain an optional third column specifying the read length.
  This is particularly useful when samples are sequenced using protocols with different read lengths.
  ```tsv
  patient1	20100113	150
  patient1	20110202	200
  patient2	20081130	150
  ```

  The utils subdirectory contain [mass-importers tools](../utils/README.md#samples-mass-importers) that can generate this third column while importing samples.

## samples

V-pipe expects the input samples to be organized in a two-level directory hierarchy.

- The first level can be, e.g., patient samples or biological replicates of an experiment.
- The second level can be, e.g., different sampling dates or different sequencing runs of the same sample.
- Inside that directory, the sub-directory `raw_data/` holds the sequencing data in FASTQ format (optionally compressed with GZip).

**For example:**

```lang-none
samples
├── patient1
│   ├── 20100113
│   │   └──raw_data
│   │      ├──patient1_20100113_R1.fastq
│   │      └──patient1_20100113_R2.fastq
│   └── 20110202
│       └──raw_data
│          ├──patient1_20100202_R1.fastq
│          └──patient1_20100202_R2.fastq
└── patient2
    └── 20081130
        └──raw_data
           ├──patient2_20081130_R1.fastq.gz
           └──patient2_20081130_R2.fastq.gz
```

The utils subdirectory contain [mass-importers tools](../utils/README.md#samples-mass-importers) to assist you in generating this hierarchy.
<!-- markdownlint-disable MD013 MD010 -->

# Utils

This directory contains several command-line utilities that can help in some ancillary tasks to get V-pipe running.

## Quick installer

`quick_install.sh` is a script that can assist deploying V-pipe: It can automatically download and install bioconda, snakemake and fetch V-pipe from the repository.

It is possible to directly run it from the web with the following commands:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -w working
cd ./working/
```

- This will download and install bioconda, snakemake and V-pipe in the current directory (use option `-p` to specify another directory).
- The installed version by default is the `master` git branch
  - It is possible to specify another git branch or a tag using the `-b` options, e.g. `-b v2.99.1`.
  - Alternatively, the `-r` option will download and install the `.tar.gz` tar-ball package of a specific version.
  - See the [release page](https://github.com/cbg-ethz/V-pipe/releases) for available tags and releases.
- using `-w` will create a working directory and populate it.
  - It will copy over a default `config.yaml` that you can [edit to your liking](../config/README.md)
  - And it will create a handy `vpipe` short-cut script to invoke `snakemake`:
```bash
cd working
# edit config.yaml and provide samples/ directory
./vpipe --jobs 4 --printshellcmds --dry-run
```

> **Tips**:
> To create and populate other new working directories, you can call `init_project.sh` from within the new directory:
>
> ```bash
> mkdir -p working_2
> cd working_2
> ../V-pipe/init_project.sh
> ```

Available options:

```console
# ./quick_install.sh -h
usage: ./quick_install.sh [options]
options:
-f           force overwriting directories
-p PREFIX    prefix directory under which to install V-pipe
             [default: current directory]
-b BRANCH    install specified branch of V-pipe
             [default: master]
-r RELEASE   install specified release package of V-pipe
-w WORKDIR   create and populate working directory
-m           only minimal working directory
-h           print this help message and exit"
```


# Samples mass-importers

Here are also some tools that can help mass-importing samples into V-pipe.
They will search for `.fastq.gz` files, put them in the [two-level hierarchy that V-pipe expects in  `samples/`](../config/README.md#samples) and generate the corresponding `samples.tsv`. By default, hard-links will be used to save space and avoid duplicating these large sequencing files.

Whichever of these tools is most suitable for you depends on how you receive the files from the sequencing lab.

- [`sort_samples_dumb`](#sort_samples_dumb) is for loose collection of FASTQ files
- [`sort_samples_demultiplexstats`](#sort_samples_demultiplexstats) and [`sort_samples_jobinfo`](#sort_samples_jobinfo) are better suited when additional information has been provided by the base calling and demultiplexing software.
  - these two also support _patch maps_ TSV files, helping to rename the samples from the name used in the original sample sheet in the LIMS (laboratory's information management system) to something more flexible.
    For example, to remap simple sequential numbers to longer names, use the following patch map TSV:

```tsv
1	05_2021_02_01
2	05_2021_02_02
3	05_2021_02_03
4	05_2021_02_04
5	05_2021_02_05
268	05_2021_02_06
269	05_2021_02_07
270	05_2021_02_08
271	05_2021_02_09
272	05_2021_02_10
273	05_2021_02_11
274	05_2021_02_12
275	05_2021_02_13
```

## sort_samples_dumb

`sort_samples_dumb` is useful when labs are providing the sequencing data simply as a loose collection of FASTQ files.

Example of usage:

```console
V-pipe/utils/sort_samples_dumb -f download/ -t working/samples.tsv -o working/samples -b 20210110
```

Where:

- `-f` : specifies the main directory containing the downloaded the `.fastq.gz` files.
  - all its subdirectories will be searched recursively
- `-t` : specifies the `samples.tsv` file to create
- `-o` : specifies the output directory where to store the files
- `-b` : is the directory to use for the second level (see V-pipe tutorials), e.g.: dates
  - the first level is usually the patient or sample name and `sort_samples_dumb` will attempt to guess it from the file names.
  - the second level is usually the sampling date or sequencing batch, but `sort_samples_dumb` has no simple way to guess that.

Running this command will immediately copy over the `.fastq.gz` files (using hard links) in the working directory and generate the `samples.tsv`. After checking the content, you should be able to run V-pipe.

Available options:

```console
./sort_samples_dumb -h
Usage: ./sort_samples_dumb -f <DIR> -b <BATCH> [-l <LEN>] [-L {''|--link|--symbolic-link|--reflink}] -o <DIR> -m <MODE>
	-f : directory containing .fastq.gz files
	-b : batch name to use for 2nd level (e.g.: date)
	-l : read lenght (default: autodetect)
	-L : link parameter to pass to cp when copying (default: --link)
	-t : tsv file (default: samples.<BATCH>.tsv)
	-T : do not truncate (empty) the file before starting
	-D : sample have duplicates (e.g.: across lanes)
	-p : prefix to prepend to fastq files (e.g.: for fusing runs)
	-s : suffix to append to fastq files (e.g.: for fusing runs)
	-o : output directory
	-m : POSIX mode parameter to pass to mkdir (e.g.: 0770)
```

Of special interest:

- `-l`: sets the read length instead of trying to autodetect it with an _awk_ script.
- `-D`: when multiple files have the same name, the importer will group them for merging as the same sample (e.g.: lane-duplicates).


## sort_samples_demultiplexstats

`sort_samples_demultiplexstats` can help if the lab provides the DemultiPlex Stats files generated by Illumina's `bcl2fastq` version 2.x software. This tool will use the information provided in the JSON file to match files to the respective samples they came from.

The general syntax is:

```console
V-pipe/utils/sort_samples_demultiplexstats --statsdir downloads/Demultiplex --fastqdir downloads/RawReads --qcdir downloads/FastQC --outdir working/samples
```

Where:

- `--statsdir` : is the directory containing the file `Stats/Stats.json` produced by `bcl2fastq`.
- `--outdir` : directory where to create the output
optionally:
- `--fastqdir`: is the directory containing the `.fastq.gz` files if they are not in the same directory as stats.
- `--qcdir`: is the directory with the FastQC's quality checks if provided by the lab.
- `--noempty`: if the lab deleted the empty (0 reads) `.fastq.gz` files.

This command is the **first step** (analysing `Stats.json`):

- It will generate a samples TSV file in the output directory: `working/samples/samples.`_{some date}_`.tsv`.
  You can copy the content of this file into your `samples.tsv`
- It will generate a file `working/samples/movedata.sh`, you can run it with:
  ```console
  bash working/samples/movedata.sh
  ```
  and that file will in turn perform the **second step**: hard-linking all the files from `download/RawReads` into `working/samples/`_{sample name}_`/`_{some date}_.
- once the second step has been performed, you should be able to run V-pipe.

Available options:

```console
./sort_samples_demultiplexstats  -h
usage: sort_samples_demultiplexstats [-h] -S DIR [-f DIR] [-q DIR] [-o DIR] [-m MODE] [-L CPLINK] [-s] [-a] [-n] [-p TSV]

Uses bcl2fastq's demultiplexing stats as metadata to organise samples

optional arguments:
  -h, --help            show this help message and exit
  -S DIR, --statsdir DIR
                        directory containing 'Stats/Stats.json'
  -f DIR, --fastqdir DIR
                        directory containing .fastq.gz files if different from above
  -q DIR, --qcdir DIR   if set, import FastQC's _fastqc.html files from there
  -o DIR, --outdir DIR  output directory
  -m MODE, --mode MODE  POSIX file access mode to be passed to mkdir
  -L CPLINK, --linking CPLINK
                        parameter to pass to `cp` for linking files instead of copying their data
  --force               Force overwriting any existing file when moving
  -s, --summary         Only display a summary of datasets, not an exhaustive list of all samples
  -a, --append          Append to the end of movedatafiles.sh, instead of overwritting (use when calling from an external combiner wrapper)
  -n, --noempty         skip fastq.gz files with bad yield (0 reads)
  -p TSV, --patchmap TSV
                        patchmap file to rename samples
```


## sort_samples_jobinfo

`sort_samples_jobinfo` is similar in concept to [`sort_samples_demultiplexstats`](#sort_samples_demultiplexstats), but it relies on the `CompletedJobInfo.xml` and `SampleSheetUsed.csv` files generated by the Illumina Analysis Software on Windows, if those are files are provided by the lab, and will try to match samples listed therein to  `.fastq.gz` files.

The general syntax is:

```console
V-pipe/utils/sort_samples_jobinfo --sourcedir=downloads/20210528_061936 --outdir=working/samples
```

Where:

- `--sourcedir` : is the directory, usually named _{yymmdd}_`_`_{hhmmss}_ and found inside the subdirectory `Alignment_1` of the Illumina's run folder (e.g.: `E:\210527_NDX550487_RUO_0008_AHTWG2AFX2\Alignment_1\20210528_061936`).
  It contains the two files `CompletedJobInfo.xml` and `SampleSheetUsed.csv` and a subdirectory named `Fastq` containing all the`.fastq.gz` sequencing files.
- `--outdir` : directory where to create the output

This command is the **first step** (analyzing `CompletedJobInfo.xml` and `SampleSheetUsed.csv` ):

- It will generate a samples TSV file in the output directory: `working/samples/samples.`_{some date}_`.tsv`.
  You can copy the content of this file into your `samples.tsv`
- It will generate a file `working/samples/movedata.sh`, you can run it with:
  ```console
  bash working/samples/movedata.sh
  ```
  and that file will in turn perform the **second step**: hard-linking all the files from `downloads/20210528_061936/Fastq` into `working/samples/`_{sample name}_`/`_{some date}_.
- once the second step has been performed, you should be able to run V-pipe.
- if the option `--batch` is provided, it will also generate an additional file `working/samples/batch.`_{some date}_`.yaml` including extra information gathered from the files (e.g. the _Library preparation kit_ listed in the input CSV). The parameter of `--batch` is used to provide the name of the lab for the `lab:` field in this file.

Available options:

```console
./sort_samples_jobinfo  -h
usage: sort_samples_jobinfo [-h] -S DIR [-f DIR] [-o DIR] [-m MODE] [-L CPLINK] [-b LAB] [-s] [-a] [-l] [-p TSV]

Uses CompletedJobInfo.xml and SampleSheetUsed.csv from Illumina Analysis Software

optional arguments:
  -h, --help            show this help message and exit
  -S DIR, --sourcedir DIR
                        directory containing CompletedJobInfo.xml and SampleSheetUsed.csv
  -f DIR, --fastqdir DIR
                        directory containing .fastq.gz files if not in 'Fastq' subdirectory
  -o DIR, --outdir DIR  output directory
  -m MODE, --mode MODE  POSIX file access mode to be passed to mkdir
  -L CPLINK, --linking CPLINK
                        parameter to pass to `cp` for linking files instead of copying their data
  --force               Force overwriting any existing file when moving
  -b LAB, --batch LAB   generate batch description
  -s, --summary         Only display a summary of datasets, not an exhaustive list of all samples
  -a, --append          Append to the end of movedatafiles.sh, instead of overwritting (use when calling from an external combiner wrapper)
  -l, --forcelanes      Explicitly look for sample in each lane (for replicates across lanes)
  -p TSV, --patchmap TSV
                        patchmap file to rename samples
```
<!-- markdownlint-disable MD013 -->

# SARS-CoV-2 test data

This data is based on the two synthetic positive controls used as part of the
[SARS-CoV-2 sequencing project done at the BSSE](https://ethz.ch/en/news-and-events/eth-news/news/2020/04/analyses-for-getting-to-grips-with-the-pandemic.html),
and was sequenced at the
[Genomics Facility Basel](https://bsse.ethz.ch/genomicsbasel)
led by Dr. Christian Beisel.

The
[Twist Synthetic SARS-CoV-2 RNA Controls](https://www.twistbioscience.com/resources/twist-synthetic-sars-cov-2-rna-controls)
were diluted 10000 on the plate and used in the initial RT step.

Library preparation adjusted from
[nCoV-2019 sequencing protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w).

A [quick guide to tiling amplicon sequencing and bioinformatics](https://artic.network/quick-guide-to-tiling-amplicon-sequencing-bioinformatics.html).

## Short

- Step 1: reverse transcription, random priming with random hexamers
- Step 2: PCR amplification of viral genome

  [nCoV-2019 PrimalSeq sequencing primers](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3)
- Step 3: Illumina adapter ligation; canonical Illumina TruSeq adapter sequences

  For [adapter trimming](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html?langsel=/us/)

  TruSeq LT and TruSeq HT-based kits:

  - Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  - Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

  Trimming should not really be necessary. Amplicons 400ish bp. Read length 251 bases

## Sequencing

[MiSeq](https://emea.illumina.com/systems/sequencing-platforms/miseq.html)
PE 2x251

(v3 600 cycles kit)

## Testdata packaging

The REF_aln.bam files generated by [V-pipe](https://cbg-ethz.github.io/V-pipe/sars-cov-2/)
(as part of the
[SARS-CoV-2 sequencing project](https://bsse.ethz.ch/cevo/cevo-press/2020/05/first-data-for-genomic-surveillance-of-sars-cov-2-in-switzerland-made-available.html))
have been:

- [capped to 512](https://www.biostars.org/p/154220/)
  using [jvarkit](https://lindenb.github.io/jvarkit/Biostar154220.html),
- [converted back to .fastq raw_reads](http://www.metagenomics.wiki/tools/samtools/converting-bam-to-fastq)
  using [`samtools bam2fq`](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/samtools/bam2fq/separate.html)
- and finally GZ-compressed using
  [`advdef --shrink-insane`](https://www.advancemame.it/comp-readme.html).
---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---


# SARS-CoV-2 Tutorial

This tutorial shows the basics of how to interact with V-pipe. A recording of our webinar covering the subject is available at the bottom of the current page.

For the purpose of this Tutorial, we will work with the sars-cov2 branch which is adapted for the SARS-CoV-2 virus.


## Organizing Data

V-pipe expects the input samples to be organized in a two-level hierarchy:

At the first level, input files grouped by samples (e.g.: patients or biological replicates of an experiment).
A second level for distinction of datasets belonging to the same sample (e.g.: sample dates).
Inside that directory, the sub-directory raw_data holds the sequencing data in FASTQ format (optionally compressed with GZip).
Paired-ended reads need to be in split files with _R1 and _R2 suffixes.


## Preparing a small dataset

You can run the first test on your workstation or a good laptop.

First, you need to prepare the data:

* For that test, you need to download the following runs from SRA: SRR10903401 and SRR10903402

If you have difficulties, check this shared directory. You can obtain there a copy of the .fastq files. More information on the steps necessary to generate the .fastq files from SRA can be found in the README.md file.

```bash
mkdir -p samples/SRR10903401/20200102/raw_data
cd samples/SRR10903401/20200102/raw_data
fasterq-dump --progress SRR10903401
cd -
```

```bash
mkdir -p samples/SRR10903402/20200102/raw_data
cd samples/SRR10903402/20200102/raw_data
fasterq-dump --progress SRR10903402
cd -
```

You then have to rename the files so that they have `_R1` and `_R2` suffixes:

```bash
mv samples/SRR10903401/20200102/raw_data/SRR10903401_1.fastq samples/SRR10903401/20200102/raw_data/SRR10903401_R1.fastq
mv samples/SRR10903401/20200102/raw_data/SRR10903401_2.fastq samples/SRR10903401/20200102/raw_data/SRR10903401_R2.fastq

mv samples/SRR10903402/20200102/raw_data/SRR10903402_1.fastq samples/SRR10903402/20200102/raw_data/SRR10903402_R1.fastq
mv samples/SRR10903402/20200102/raw_data/SRR10903402_2.fastq samples/SRR10903402/20200102/raw_data/SRR10903402_R2.fastq
```

The downloaded files will have the following structure:

```bash
tree samples
```


## Install V-pipe

V-pipe uses the [Bioconda](https://bioconda.github.io/) bioinformatics software repository for all its pipeline components. The pipeline itself is written using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

For advanced users: If your are fluent with these tools, you can:

* directly download and install [bioconda](https://bioconda.github.io/user/install.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda),
* make sure to configure V-pipe to use the `sars-cov-2` virus-config
* and start using V-pipe with them, using the --use-conda to [automatically download and install](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) any further pipeline dependencies.
* please refer to the documentation for additional instructions.

The present tutorial will show simplified commands that automate much of this process.

To deploy V-pipe, you can use the installation script with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p testing -w work
```

* using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
* using `-w` will create a working directory and populate it. It will copy over the references and the default `config/config.yaml`, and create a handy `vpipe` short-cut script to invoke `snakemake`.

Tip: To create and populate other new working directories, you can call init_project.sh from within the new directory:

```console
mkdir -p working_2
cd working_2
../V-pipe/init_project.sh
```


## Running V-pipe

Copy the samples directory you created in the step Preparing a small dataset to this working directory. You can display the directory structure with tree samples or find samples3.

```bash
mv ./samples ./testing/work/
```

Prepare V-pipe's configuration:

```bash
cat <<EOT >> ./testing/work/config.yaml
general:
    virus_base_config: 'sars-cov-2'

output:
    snv: false
    local: false
    global: false
    visualization: false
    QA: false
EOT
```

Check what will be executed:

```bash
cd ./testing/work/
./vpipe --dryrun
```

As it is your first run of V-pipe, this will also generate the sample collection table. Check `samples.tsv` in your editor.

Note that the demo files you downloaded have reads of length 150 only. V-pipe’s default parameters are optimized for reads of length 250; add the third column in the tab-separated file:

```bash
cat <<EOT > ./testing/work/samples.tsv
SRR10903401	20200102	150
SRR10903402	20200102	150
EOT
```

Tip: Always check the content of the `samples.tsv` file.

If you didn’t use the correct structure, this file might end up empty or some entries might be missing.
You can safely delete it and re-run the `--dryrun` to regenerate it.

Run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
cd ./testing/work/
./vpipe -p --cores 2
```


## Output

The Wiki contains an overview of the output files. The output of the SNV calling is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in `samples/​{hierarchy}​/variants/SNVs/snvs.vcf`, you can open it with your favorite VCF tools for visualisation or downstream processing. It is also available in a tabular format in `samples/​{hierarchy}​/variants/SNVs/snvs.csv`.

### Expected output

The small dataset that we used in this tutorial section has been analyzed by [doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036). The results of the original analysis (using bwa, samtools mpileup, and bcftools) are displayed in Table 2 in the article:

Using either the VCF or CSV files, compare with the results given out by V-pipe (with bwa and ShoRAH).

* For positions 19164 and 24323 of SRR10903401 and position 11563 of SRR10903402, we expect to see similar results in V-pipe.
* For the remaining positions (1821, 26314 and 26590 of SRR10903401), we expect that ShoRAH will consider the variants of poor quality and reject them because there is very little support ( <= than 5 reads supporting the alt).


## Swapping component

The default configuration uses ShoRAH to call the SNVs and to reconstruct the local (windowed) haplotypes.

Components can be swapped simply by changing the `config.yaml` file. For example to call SNVs using lofreq:

```yaml
general:
  snv_caller: lofreq
```


## Cluster deployment

It is possible to ask snakemake to submit jobs on a cluster using the batch submission command-line interface of your cluster.

The platform LSF by IBM is one of the popular systems you might find (Others include SLURM, Grid Engine).


...TODO...
