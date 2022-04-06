# ![nf-core/mag](docs/images/nf-core-mag_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/mag/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mag/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/mag/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/mag/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/mag/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3589527-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3589527)
[![Cite Preprint](https://img.shields.io/badge/Cite%20Us!-Cite%20Preprint-orange)](https://doi.org/10.1101/2021.08.29.458094)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mag-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/mag)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/mag** is a bioinformatics best-practise analysis pipeline for assembly, binning and annotation of metagenomes.

<p align="center">
    <img src="docs/images/mag_workflow.png" alt="nf-core/mag workflow overview" width="90%">
</p>

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/mag/results).

## Pipeline summary

By default, the pipeline currently performs the following: it supports both short and long reads, quality trims the reads and adapters with [fastp](https://github.com/OpenGene/fastp) and [Porechop](https://github.com/rrwick/Porechop), and performs basic QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
The pipeline then:

* assigns taxonomy to reads using [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) and/or [Kraken2](https://github.com/DerrickWood/kraken2/wiki)
* performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks their quality using [Quast](http://quast.sourceforge.net/quast)
* predicts protein-coding genes for the assemblies using [Prodigal](https://github.com/hyattpd/Prodigal)
* performs metagenome binning using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), and checks the quality of the genome bins using [Busco](https://busco.ezlab.org/)
* assigns taxonomy to bins using [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) and/or [CAT](https://github.com/dutilh/CAT)

Furthermore, the pipeline creates various reports in the results directory specified, including a [MultiQC](https://multiqc.info/) report summarizing some of the findings and software versions.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/mag -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run nf-core/mag -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '*_R{1,2}.fastq.gz'
    ```

    or

    ```console
    nextflow run nf-core/mag -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv
    ```

See [usage docs](https://nf-co.re/mag/usage) and [parameter docs](https://nf-co.re/mag/parameters) for all of the available options when running the pipeline.

## Documentation

The nf-core/mag pipeline comes with documentation about the pipeline [usage](https://nf-co.re/mag/usage), [parameters](https://nf-co.re/mag/parameters) and [output](https://nf-co.re/mag/output). Detailed information about how to specify the input can be found under [input specifications](https://nf-co.re/mag/usage#input_specifications).

### Group-wise co-assembly and co-abundance computation

Each sample has an associated group ID (see [input specifications](https://nf-co.re/mag/usage#input_specifications)). This group information can be used for group-wise co-assembly with `MEGAHIT` or `SPAdes` and/or to compute co-abundances for the binning step with `MetaBAT2`. By default, group-wise co-assembly is disabled, while the computation of group-wise co-abundances is enabled. For more information about how this group information can be used see the documentation for the parameters [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode).

When group-wise co-assembly is enabled, `SPAdes` is run on accordingly pooled read files, since `metaSPAdes` does not yet allow the input of multiple samples or libraries. In contrast, `MEGAHIT` is run for each group while supplying lists of the individual readfiles.

## Credits

nf-core/mag was written by [Hadrien Gourlé](https://hadriengourle.com) at [SLU](https://slu.se), [Daniel Straub](https://github.com/d4straub) and [Sabrina Krakau](https://github.com/skrakau) at the [Quantitative Biology Center (QBiC)](http://qbic.life).

Long read processing was inspired by [caspargross/HybridAssembly](https://github.com/caspargross/HybridAssembly) written by Caspar Gross [@caspargross](https://github.com/caspargross)

We thank the following people for their extensive assistance in the development of this pipeline:

* [Alexander Peltzer](https://github.com/apeltzer)
* [Antonia Schuster](https://github.com/antoniaschuster)
* [Phil Ewels](https://github.com/ewels)
* [Gisela Gabernet](https://github.com/ggabernet)
* [Harshil Patel](https://github.com/drpatelh)
* [Johannes Alneberg](https://github.com/alneberg)
* [Maxime Borry](https://github.com/maxibor)
* [Maxime Garcia](https://github.com/MaxUlysse)
* [Michael L Heuer](https://github.com/heuermh)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mag` channel](https://nfcore.slack.com/channels/mag) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/mag for your analysis, please cite the preprint as follows:

> **nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning**
>
> Sabrina Krakau, Daniel Straub, Hadrien Gourlé, Gisela Gabernet, Sven Nahnsen.
>
> bioRxiv 2021.08.29.458094. doi: [10.1101/2021.08.29.458094](https://doi.org/10.1101/2021.08.29.458094).

additionally you can cite the pipeline directly with the following doi: [10.5281/zenodo.3589527](https://doi.org/10.5281/zenodo.3589527)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
# nf-core/mag: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.1 - 2021/11/25

### `Added`

- [#240](https://github.com/nf-core/mag/pull/240) - Add prodigal to predict protein-coding genes for assemblies.
- [#241](https://github.com/nf-core/mag/pull/241) - Add parameter `--skip_prodigal`.
- [#244](https://github.com/nf-core/mag/pull/244) - Add pipeline preprint information.
- [#245](https://github.com/nf-core/mag/pull/245) - Add Prokka to annotate binned genomes.

### `Changed`

- [#249](https://github.com/nf-core/mag/pull/249) - Update workflow overview figure.
- [#258](https://github.com/nf-core/mag/pull/258) - Updated MultiQC 1.9 to 1.11.
- [#260](https://github.com/nf-core/mag/pull/260) - Updated SPAdes 3.13.1 -> 3.15.3, MEGAHIT 1.2.7 -> 1.2.7

### `Fixed`

- [#256](https://github.com/nf-core/mag/pull/256) - Fix `--skip_busco`.
- [#236](https://github.com/nf-core/mag/pull/236) - Fix large assemblies (> 4 billion nucleotides in length).
- [#254](https://github.com/nf-core/mag/pull/254) - Fix MetaBAT2 error with nextflow version 21.10.x (21.04.03 is the latest functional version for nf-core/mag 2.1.0).
- [#255](https://github.com/nf-core/mag/pull/255) - Update gtdbtk conda channel.
- [#258](https://github.com/nf-core/mag/pull/258) - FastP results are now in MultiQC.

## v2.1.0 - 2021/07/29

### `Added`

- [#212](https://github.com/nf-core/mag/pull/212), [#214](https://github.com/nf-core/mag/pull/214) - Add bin abundance estimation based on median sequencing depths of corresponding contigs (results are written to `results/GenomeBinning/bin_depths_summary.tsv` and `results/GenomeBinning/bin_summary.tsv`) [#197](https://github.com/nf-core/mag/issues/197).
- [#214](https://github.com/nf-core/mag/pull/214) - Add generation of (clustered) heat maps with bin abundances across samples (using centered log-ratios)
- [#217](https://github.com/nf-core/mag/pull/217) - Publish genes predicted with Prodigal within BUSCO run (written to `results/GenomeBinning/QC/BUSCO/[assembler]-[bin]_prodigal.gff`).

### `Changed`

- [#218](https://github.com/nf-core/mag/pull/218) - Update to nf-core 2.0.1 `TEMPLATE` (DSL2)

### `Fixed`

- [#226](https://github.com/nf-core/mag/pull/226) - Fix handling of `BUSCO` output when run in auto lineage selection mode and selected specific lineage is the same as the generic one.

## v2.0.0 - 2021/06/01

### `Added`

- [#179](https://github.com/nf-core/mag/pull/179) - Add BUSCO automated lineage selection functionality (new default). The pameter `--busco_auto_lineage_prok` can be used to only consider prokaryotes and the parameter `--busco_download_path` to run BUSCO in `offline` mode.
- [#178](https://github.com/nf-core/mag/pull/178) - Add taxonomic bin classification with `GTDB-Tk` `v1.5.0` (for bins filtered based on `BUSCO` QC metrics).
- [#196](https://github.com/nf-core/mag/pull/196) - Add process for CAT database creation as an alternative to using pre-built databases.

### `Changed`

- [#162](https://github.com/nf-core/mag/pull/162) - Switch to DSL2
- [#162](https://github.com/nf-core/mag/pull/162) - Changed `--input` file format from `TSV` to `CSV` format, requires header now
- [#162](https://github.com/nf-core/mag/pull/162) - Update `README.md`, `docs/usage.md` and `docs/output.md`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastP` from version `0.20.0` to `0.20.1`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `Bowtie2` from version `2.3.5` to `2.4.2`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastQC` from version `0.11.8` to `0.11.9`
- [#172](https://github.com/nf-core/mag/pull/172) - Compressed discarded MetaBAT2 output files
- [#176](https://github.com/nf-core/mag/pull/176) - Update CAT DB link
- [#179](https://github.com/nf-core/mag/pull/179) - Update `BUSCO` from version `4.1.4` to `5.1.0`
- [#179](https://github.com/nf-core/mag/pull/179) - By default BUSCO now performs automated lineage selection instead of using the bacteria_odb10 lineage as reference. Specific lineage datasets can still be provided via `--busco_reference`.
- [#178](https://github.com/nf-core/mag/pull/178) - Change output file: `results/GenomeBinning/QC/quast_and_busco_summary.tsv` -> `results/GenomeBinning/bin_summary.tsv`, contains GTDB-Tk results as well.
- [#191](https://github.com/nf-core/mag/pull/191) - Update to nf-core 1.14 `TEMPLATE`
- [#193](https://github.com/nf-core/mag/pull/193) - Compress CAT output files [#180](https://github.com/nf-core/mag/issues/180)
- [#198](https://github.com/nf-core/mag/pull/198) - Requires nextflow version `>= 21.04.0`
- [#200](https://github.com/nf-core/mag/pull/200) - Small changes in GitHub Actions tests
- [#203](https://github.com/nf-core/mag/pull/203) - Renamed `fastp` params and improved description in documentation: `--mean_quality` -> `--fastp_qualified_quality`, `--trimming_quality` -> `--fastp_cut_mean_quality`

### `Fixed`

- [#175](https://github.com/nf-core/mag/pull/175) - Fix bug in retrieving the `--max_unbinned_contigs` longest unbinned sequences that are longer than `--min_length_unbinned_contigs` (`split_fasta.py`)
- [#175](https://github.com/nf-core/mag/pull/175) - Improved runtime of `split_fasta.py` in `METABAT2` process (important for large assemblies, e.g. when computing co-assemblies)
- [#194](https://github.com/nf-core/mag/pull/194) - Allow different folder structures for Kraken2 databases containing `*.k2d` files [#187](https://github.com/nf-core/mag/issues/187)
- [#195](https://github.com/nf-core/mag/pull/195) - Fix documentation regarding required compression of input FastQ files [#160](https://github.com/nf-core/mag/issues/160)
- [#196](https://github.com/nf-core/mag/pull/196) - Add process for CAT database creation as solution for problem caused by incompatible `DIAMOND` version used for pre-built `CAT database` and `CAT classification` [#90](https://github.com/nf-core/mag/issues/90), [#188](https://github.com/nf-core/mag/issues/188)

## v1.2.0 - 2021/02/10

### `Added`

- [#146](https://github.com/nf-core/mag/pull/146) - Add `--coassemble_group` parameter to allow group-wise co-assembly
- [#146](https://github.com/nf-core/mag/pull/146) - Add `--binning_map_mode` parameter allowing different mapping strategies to compute co-abundances used for binning (`all`, `group`, `own`)
- [#149](https://github.com/nf-core/mag/pull/149) - Add two new parameters to allow custom SPAdes and MEGAHIT options (`--spades_options` and `--megahit_options`)

### `Changed`

- [#141](https://github.com/nf-core/mag/pull/141) - Update to nf-core 1.12.1 `TEMPLATE`
- [#143](https://github.com/nf-core/mag/pull/143) - Manifest file has to be handed over via `--input` parameter now
- [#143](https://github.com/nf-core/mag/pull/143) - Changed format of manifest input file: requires a '.tsv' suffix and additionally contains group ID
- [#143](https://github.com/nf-core/mag/pull/143) - TSV `--input` file allows now also entries containing only short reads
- [#145](https://github.com/nf-core/mag/pull/145) - When using TSV input files, uses sample IDs now for `FastQC` instead of basenames of original read files. Allows non-unique file basenames.

### `Removed`

- [#143](https://github.com/nf-core/mag/pull/143) - Change parameter: `--manifest` -> `--input`

## v1.1.2 - 2020/11/24

### `Changed`

- [#135](https://github.com/nf-core/mag/pull/135) - Update to nf-core 1.12 `TEMPLATE`

### `Fixed`

- [#133](https://github.com/nf-core/mag/pull/133) - Fixed processing of `--input` parameter [#131](https://github.com/nf-core/mag/issues/131)

## v1.1.1 - 2020/11/10

### `Added`

- [#121](https://github.com/nf-core/mag/pull/121) - Add full-size test
- [#124](https://github.com/nf-core/mag/pull/124) - Add worfklow overview figure to `README`

### `Changed`

- [#123](https://github.com/nf-core/mag/pull/123) - Update to new nf-core 1.11 `TEMPLATE`

### `Fixed`

- [#118](https://github.com/nf-core/mag/pull/118) - Fix `seaborn` to `v0.10.1` to avoid `nanoplot` error
- [#120](https://github.com/nf-core/mag/pull/120) - Fix link to CAT database in help message
- [#124](https://github.com/nf-core/mag/pull/124) - Fix description of `CAT` process in `output.md`

## v1.1.0 - 2020/10/06

### `Added`

- [#35](https://github.com/nf-core/mag/pull/35) - Add social preview image
- [#49](https://github.com/nf-core/mag/pull/49) - Add host read removal with `Bowtie 2` and according custom section to `MultiQC`
- [#49](https://github.com/nf-core/mag/pull/49) - Add separate `MultiQC` section for `FastQC` after preprocessing
- [#65](https://github.com/nf-core/mag/pull/65) - Add `MetaBAT2` RNG seed parameter `--metabat_rng_seed` and set the default to 1 which ensures reproducible binning results
- [#65](https://github.com/nf-core/mag/pull/65) - Add parameters `--megahit_fix_cpu_1`, `--spades_fix_cpus` and `--spadeshybrid_fix_cpus` to ensure reproducible results from assembly tools
- [#66](https://github.com/nf-core/mag/pull/66) - Export `depth.txt.gz` into result folder
- [#67](https://github.com/nf-core/mag/pull/67) - Compress assembly files
- [#82](https://github.com/nf-core/mag/pull/82) - Add `nextflow_schema.json`
- [#104](https://github.com/nf-core/mag/pull/104) - Add parameter `--save_busco_reference`

### `Changed`

- [#56](https://github.com/nf-core/mag/pull/56) - Update `MetaBAT2` from `v2.13` to `v2.15`
- [#46](https://github.com/nf-core/mag/pull/46) - Update `MultiQC` from `v1.7` to `v1.9`
- [#88](https://github.com/nf-core/mag/pull/88) - Update to new nf-core 1.10.2 `TEMPLATE`
- [#88](https://github.com/nf-core/mag/pull/88) - `--reads` is now removed, use `--input` instead
- [#101](https://github.com/nf-core/mag/pull/101) - Prevented PhiX alignments from being stored in work directory [#97](https://github.com/nf-core/mag/issues/97)
- [#104](https://github.com/nf-core/mag/pull/104), [#111](https://github.com/nf-core/mag/pull/111) - Update `BUSCO` from `v3.0.2` to `v4.1.4`

### `Fixed`

- [#29](https://github.com/nf-core/mag/pull/29) - Fix `MetaBAT2` binning discards unbinned contigs [#27](https://github.com/nf-core/mag/issues/27)
- [#31](https://github.com/nf-core/mag/pull/31), [#36](https://github.com/nf-core/mag/pull/36), [#76](https://github.com/nf-core/mag/pull/76), [#107](https://github.com/nf-core/mag/pull/107) - Fix links in README
- [#47](https://github.com/nf-core/mag/pull/47) - Fix missing `MultiQC` when `--skip_quast` or `--skip_busco` was specified
- [#49](https://github.com/nf-core/mag/pull/49), [#89](https://github.com/nf-core/mag/pull/89) - Added missing parameters to summary
- [#50](https://github.com/nf-core/mag/pull/50) - Fix missing channels when `--keep_phix` is specified
- [#54](https://github.com/nf-core/mag/pull/54) - Updated links to `minikraken db`
- [#54](https://github.com/nf-core/mag/pull/54) - Fixed `Kraken2` dp preparation: allow different names for compressed archive file and contained folder as for some minikraken dbs
- [#55](https://github.com/nf-core/mag/pull/55) - Fixed channel joining for multiple samples causing `MetaBAT2` error [#32](https://github.com/nf-core/mag/issues/32)
- [#57](https://github.com/nf-core/mag/pull/57) - Fix number of threads used by `MetaBAT2` program `jgi_summarize_bam_contig_depths`
- [#70](https://github.com/nf-core/mag/pull/70) - Fix `SPAdes` memory conversion issue [#61](https://github.com/nf-core/mag/issues/61)
- [#71](https://github.com/nf-core/mag/pull/71) - No more ignoring errors in `SPAdes` assembly
- [#72](https://github.com/nf-core/mag/pull/72) - No more ignoring of `BUSCO` errors
- [#73](https://github.com/nf-core/mag/pull/73), [#75](https://github.com/nf-core/mag/pull/75) - Improved output documentation
- [#96](https://github.com/nf-core/mag/pull/96) - Fix missing bin names in `MultiQC` BUSCO section [#78](https://github.com/nf-core/mag/issues/78)
- [#104](https://github.com/nf-core/mag/pull/104) - Fix `BUSCO` errors causing missing summary output [#77](https://github.com/nf-core/mag/issues/77)

### `Deprecated`

- [#29](https://github.com/nf-core/mag/pull/29) - Change depreciated parameters: `--singleEnd` -> `--single_end`, `--igenomesIgnore` -> `--igenomes_ignore`

## v1.0.0 - 2019/12/20

Initial release of nf-core/mag, created with the [nf-core](http://nf-co.re/) template.

### `Added`

- short and long reads QC (fastp, porechop, filtlong, fastqc)
- Lambda and PhiX detection and filtering (bowtie2, nanolyse)
- Taxonomic classification of reads (centrifuge, kraken2)
- Short read and hybrid assembly (megahit, metaspades)
- metagenome binning (metabat2)
- QC of bins (busco, quast)
- annotation (cat/bat)
# Code of Conduct at nf-core (v1.0)

## Our Pledge

In the interest of fostering an open, collaborative, and welcoming environment, we as contributors and maintainers of nf-core, pledge to making participation in our projects and community a harassment-free experience for everyone, regardless of:

- Age
- Body size
- Familial status
- Gender identity and expression
- Geographical location
- Level of experience
- Nationality and national origins
- Native language
- Physical and neurological ability
- Race or ethnicity
- Religion
- Sexual identity and orientation
- Socioeconomic status

Please note that the list above is alphabetised and is therefore not ranked in any order of preference or importance.

## Preamble

> Note: This Code of Conduct (CoC) has been drafted by the nf-core Safety Officer and been edited after input from members of the nf-core team and others. "We", in this document, refers to the Safety Officer and members of the nf-core core team, both of whom are deemed to be members of the nf-core community and are therefore required to abide by this Code of Conduct. This document will amended periodically to keep it up-to-date, and in case of any dispute, the most current version will apply.

An up-to-date list of members of the nf-core core team can be found [here](https://nf-co.re/about). Our current safety officer is Renuka Kudva.

nf-core is a young and growing community that welcomes contributions from anyone with a shared vision for [Open Science Policies](https://www.fosteropenscience.eu/taxonomy/term/8). Open science policies encompass inclusive behaviours and we strive to build and maintain a safe and inclusive environment for all individuals.

We have therefore adopted this code of conduct (CoC), which we require all members of our community and attendees in nf-core events to adhere to in all our workspaces at all times. Workspaces include but are not limited to Slack, meetings on Zoom, Jitsi, YouTube live etc.

Our CoC will be strictly enforced and the nf-core team reserve the right to exclude participants who do not comply with our guidelines from our workspaces and future nf-core activities.

We ask all members of our community to help maintain a supportive and productive workspace and to avoid behaviours that can make individuals feel unsafe or unwelcome. Please help us maintain and uphold this CoC.

Questions, concerns or ideas on what we can include? Contact safety [at] nf-co [dot] re

## Our Responsibilities

The safety officer is responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behaviour.

The safety officer in consultation with the nf-core core team have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

Members of the core team or the safety officer who violate the CoC will be required to recuse themselves pending investigation. They will not have access to any reports of the violations and be subject to the same actions as others in violation of the CoC.

## When are where does this Code of Conduct apply?

Participation in the nf-core community is contingent on following these guidelines in all our workspaces and events. This includes but is not limited to the following listed alphabetically and therefore in no order of preference:

- Communicating with an official project email address.
- Communicating with community members within the nf-core Slack channel.
- Participating in hackathons organised by nf-core (both online and in-person events).
- Participating in collaborative work on GitHub, Google Suite, community calls, mentorship meetings, email correspondence.
- Participating in workshops, training, and seminar series organised by nf-core (both online and in-person events). This applies to events hosted on web-based platforms such as Zoom, Jitsi, YouTube live etc.
- Representing nf-core on social media. This includes both official and personal accounts.

## nf-core cares 😊

nf-core's CoC and expectations of respectful behaviours for all participants (including organisers and the nf-core team) include but are not limited to the following (listed in alphabetical order):

- Ask for consent before sharing another community member’s personal information (including photographs) on social media.
- Be respectful of differing viewpoints and experiences. We are all here to learn from one another and a difference in opinion can present a good learning opportunity.
- Celebrate your accomplishments at events! (Get creative with your use of emojis 🎉 🥳 💯 🙌 !)
- Demonstrate empathy towards other community members. (We don’t all have the same amount of time to dedicate to nf-core. If tasks are pending, don’t hesitate to gently remind members of your team. If you are leading a task, ask for help if you feel overwhelmed.)
- Engage with and enquire after others. (This is especially important given the geographically remote nature of the nf-core community, so let’s do this the best we can)
- Focus on what is best for the team and the community. (When in doubt, ask)
- Graciously accept constructive criticism, yet be unafraid to question, deliberate, and learn.
- Introduce yourself to members of the community. (We’ve all been outsiders and we know that talking to strangers can be hard for some, but remember we’re interested in getting to know you and your visions for open science!)
- Show appreciation and **provide clear feedback**. (This is especially important because we don’t see each other in person and it can be harder to interpret subtleties. Also remember that not everyone understands a certain language to the same extent as you do, so **be clear in your communications to be kind.**)
- Take breaks when you feel like you need them.
- Using welcoming and inclusive language. (Participants are encouraged to display their chosen pronouns on Zoom or in communication on Slack.)

## nf-core frowns on 😕

The following behaviours from any participants within the nf-core community (including the organisers) will be considered unacceptable under this code of conduct. Engaging or advocating for any of the following could result in expulsion from nf-core workspaces.

- Deliberate intimidation, stalking or following and sustained disruption of communication among participants of the community. This includes hijacking shared screens through actions such as using the annotate tool in conferencing software such as Zoom.
- “Doxing” i.e. posting (or threatening to post) another person’s personal identifying information online.
- Spamming or trolling of individuals on social media.
- Use of sexual or discriminatory imagery, comments, or jokes and unwelcome sexual attention.
- Verbal and text comments that reinforce social structures of domination related to gender, gender identity and expression, sexual orientation, ability, physical appearance, body size, race, age, religion or work experience.

### Online Trolling

The majority of nf-core interactions and events are held online. Unfortunately, holding events online comes with the added issue of online trolling. This is unacceptable, reports of such behaviour will be taken very seriously, and perpetrators will be excluded from activities immediately.

All community members are required to ask members of the group they are working within for explicit consent prior to taking screenshots of individuals during video calls.

## Procedures for Reporting CoC violations

If someone makes you feel uncomfortable through their behaviours or actions, report it as soon as possible.

You can reach out to members of the [nf-core core team](https://nf-co.re/about) and they will forward your concerns to the safety officer(s).

Issues directly concerning members of the core team will be dealt with by other members of the core team and the safety manager, and possible conflicts of interest will be taken into account. nf-core is also in discussions about having an ombudsperson, and details will be shared in due course.

All reports will be handled with utmost discretion and confidentially.

## Attribution and Acknowledgements

- The [Contributor Covenant, version 1.4](http://contributor-covenant.org/version/1/4)
- The [OpenCon 2017 Code of Conduct](http://www.opencon2017.org/code_of_conduct) (CC BY 4.0 OpenCon organisers, SPARC and Right to Research Coalition)
- The [eLife innovation sprint 2020 Code of Conduct](https://sprint.elifesciences.org/code-of-conduct/)
- The [Mozilla Community Participation Guidelines v3.1](https://www.mozilla.org/en-US/about/governance/policies/participation/) (version 3.1, CC BY-SA 3.0 Mozilla)

## Changelog

### v1.0 - March 12th, 2021

- Complete rewrite from original [Contributor Covenant](http://contributor-covenant.org/) CoC.
# nf-core/mag: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

* [Bowtie2](https:/dx.doi.org/10.1038/nmeth.1923)
    > Langmead, B. and Salzberg, S. L. 2012 Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), p. 357–359. doi: 10.1038/nmeth.1923.

* [Busco](https://doi.org/10.1007/978-1-4939-9173-0_14)
    > Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: assessing genome assembly and annotation completeness. In Gene prediction (pp. 227-245). Humana, New York, NY. doi: 10.1007/978-1-4939-9173-0_14.

* [CAT](https://doi.org/10.1186/s13059-019-1817-x)
    > von Meijenfeldt, F. B., Arkhipova, K., Cambuy, D. D., Coutinho, F. H., & Dutilh, B. E. (2019). Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT. Genome biology, 20(1), 1-14. doi: 10.1186/s13059-019-1817-x.

* [Centrifuge](https://doi.org/10.1101/gr.210641.116)
    > Kim, D., Song, L., Breitwieser, F. P., & Salzberg, S. L. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome research, 26(12), 1721-1729. doi: 10.1101/gr.210641.116.

* [FastP](https://doi.org/10.1093/bioinformatics/bty560)
    > Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884–i890. doi: 10.1093/bioinformatics/bty560.

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* [Filtlong](https://github.com/rrwick/Filtlong)

* [GTDB-Tk](https://doi.org/10.1093/bioinformatics/btz848)
    > Chaumeil, P. A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics , 36(6), 1925–1927. doi: 10.1093/bioinformatics/btz848.

* [Kraken2](https://doi.org/10.1186/s13059-019-1891-0)
    > Wood, D et al., 2019. Improved metagenomic analysis with Kraken 2. Genome Biology volume 20, Article number: 257. doi: 10.1186/s13059-019-1891-0.

* [Krona](https://doi.org/10.1186/1471-2105-12-385)
    > Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC bioinformatics, 12(1), 1-10. doi: 10.1186/1471-2105-12-385.

* [MEGAHIT](https://doi.org/10.1016/j.ymeth.2016.02.020)
    > Li, D., Luo, R., Liu, C. M., Leung, C. M., Ting, H. F., Sadakane, K., ... & Lam, T. W. (2016). MEGAHIT v1. 0: a fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods, 102, 3-11. doi: 10.1016/j.ymeth.2016.02.020.

* [MetaBAT2](https://doi.org/10.7717/peerj.7359)
    > Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. doi: 10.7717/peerj.7359.

* [MultiQC](https://doi.org/10.1093/bioinformatics/btw354)
    > Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: doi.org/10.1093/bioinformatics/btw354.

* [NanoLyse](https://doi.org/10.1093/bioinformatics/bty149)
    > De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666-2669. doi: 10.1093/bioinformatics/bty149.

* [NanoPlot](https://doi.org/10.1093/bioinformatics/bty149)
    > De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666-2669. doi: 10.1093/bioinformatics/bty149.

* [Porechop](https://github.com/rrwick/Porechop)

* [Prodigal](https://pubmed.ncbi.nlm.nih.gov/20211023/)
    > Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

* [Prokka](https://pubmed.ncbi.nlm.nih.gov/24642063/)
    > Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014 Jul 15;30(14):2068-9. doi: 10.1093/bioinformatics/btu153. Epub 2014 Mar 18. PMID: 24642063.

* [SAMtools](https://doi.org/10.1093/bioinformatics/btp352)
    > Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics , 25(16), 2078–2079. doi: 10.1093/bioinformatics/btp352.

* [SPAdes](https://doi.org/10.1101/gr.213959.116)
    > Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834. doi: 10.1101/gr.213959.116.

## Data

* [Full-size test data](https://doi.org/10.1038/s41587-019-0191-2)
    > Bertrand, D., Shaw, J., Kalathiyappan, M., Ng, A. H. Q., Kumar, M. S., Li, C., ... & Nagarajan, N. (2019). Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes. Nature biotechnology, 37(8), 937-944. doi: 10.1038/s41587-019-0191-2.

## Software packaging/containerisation tools

* [Anaconda](https://anaconda.com)
    > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)
    > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

* [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)
    > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

* [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

* [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)
    > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
# nf-core/mag: Documentation

The nf-core/mag documentation is split into the following pages:

* [Usage](usage.md)
    * An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
* [Output](output.md)
    * An overview of the different results produced by the pipeline and how to interpret them.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)
# nf-core/mag: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mag/usage](https://nf-co.re/mag/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Input specifications

The input data can be passed to nf-core/mag in two possible ways using the `--input` parameter.

### Direct FASTQ input (short reads only)

The easiest way is to specify directly the path (with wildcards) to your input FASTQ files. For example:

```console
--input 'path/to/data/sample_*_R{1,2}.fastq.gz'
```

This input method only works with short read data and will assign all files to the same group. By default, this group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

Please note the following additional requirements:

* Files names must be unique
* Valid file extensions: `.fastq.gz`, `.fq.gz` (files must be compressed)
* The path must be enclosed in quotes
* The path must have at least one `*` wildcard character
* When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs
* To run single-end data you must additionally specify `--single_end`
* If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### Samplesheet input file

Alternatively, to assign different groups or to include long reads for hybrid assembly with metaSPAdes, you can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata.

This CSV file should contain the following columns:

`sample,group,short_reads_1,short_reads_2,long_reads`

The path to `long_reads` and `short_reads_2` is optional. Valid examples could look like the following:

```console
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
```

or

```console
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1.fastq.gz,,
sample2,0,data/sample2.fastq.gz,,
```

Please note the following requirements:

* 5 comma-seperated columns
* Valid file extension: `.csv`
* Must contain the header `sample,group,short_reads_1,short_reads_2,long_reads`
* Sample IDs must be unique
* FastQ files must be compressed (`.fastq.gz`, `.fq.gz`)
* `long_reads` can only be provided in combination with paired-end short read data
* Within one samplesheet either only single-end or only paired-end reads can be specified
* If single-end reads are specified, the command line parameter `--single_end` must be specified as well

Again, by default, the group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/mag --input samplesheet.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

See the [nf-core/mag website documentation](https://nf-co.re/mag/usage#usage) for more information about pipeline specific parameters.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/mag
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

Additionally, to enable also reproducible results from the individual assembly tools this pipeline provides extra parameters. SPAdes is designed to be deterministic for a given number of threads. To generate reproducible results set the number of cpus with `--spades_fix_cpus` or `--spadeshybrid_fix_cpus`. This will overwrite the number of cpus specified in the `base.config` file and additionally ensure that it is not increased in case of retries for individual samples. MEGAHIT only generates reproducible results when run single-threaded.
You can fix this by using the prameter `--megahit_fix_cpu_1`. In both cases, do not specify the number of cpus for these processes in additional custom config files, this would result in an error.

MetaBAT2 is run by default with a fixed seed within this pipeline, thus producing reproducible results.

To allow also reproducible bin QC with BUSCO, run BUSCO providing already downloaded lineage datasets with `--busco_download_path` (BUSCO will be run using automated lineage selection in offline mode) or provide a specific lineage dataset via `--busco_reference` and use the parameter `--save_busco_reference`. This may be useful since BUSCO datasets are frequently updated and old versions do not always remain (easily) accessible.

For the taxonomic bin classification with [CAT](https://github.com/dutilh/CAT), when running the pipeline with `--cat_db_generate` the parameter `--save_cat_db` can be used to also save the generated database to allow reproducibility in future runs. Note that when specifying a pre-built database with `--cat_db`, currently the database can not be saved.

The taxonomic classification of bins with GTDB-Tk is not guaranteed to be reproducible, since the placement of bins in the reference tree is non-deterministic. However, the authors of the GTDB-Tk article examined the reproducibility on a set of 100 genomes across 50 trials and did not observe any difference (see [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)).

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`, `test_hybrid`, `test_host_rm`, `test_hybrid_host_rm`, `test_busco_auto`
    * Profiles with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

Note, do not change number of CPUs with custom config files for the processes `spades`, `spadeshybrid` or `megahit` when specifying the parameters `--spades_fix_cpus`, `--spadeshybrid_fix_cpus` and `--megahit_fix_cpu_1` respectively.

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Tool-specific options

For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change tool-specific command-line arguments (e.g. providing an additional command-line argument to the `STAR_ALIGN` process) as well as publishing options (e.g. saving files produced by the `STAR_ALIGN` process that aren't saved by default by the pipeline). In the majority of instances, as a user you won't have to change the default options set by the pipeline developer(s), however, there may be edge cases where creating a simple custom config file can improve the behaviour of the pipeline if for example it is failing due to a weird error that requires setting a tool-specific parameter to deal with smaller / larger genomes.

The command-line arguments passed to STAR in the `STAR_ALIGN` module are a combination of:

* Mandatory arguments or those that need to be evaluated within the scope of the module, as supplied in the [`script`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L49-L55) section of the module file.

* An [`options.args`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L56) string of non-mandatory parameters that is set to be empty by default in the module but can be overwritten when including the module in the sub-workflow / workflow context via the `addParams` Nextflow option.

The nf-core/rnaseq pipeline has a sub-workflow (see [terminology](https://github.com/nf-core/modules#terminology)) specifically to align reads with STAR and to sort, index and generate some basic stats on the resulting BAM files using SAMtools. At the top of this file we import the `STAR_ALIGN` module via the Nextflow [`include`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L10) keyword and by default the options passed to the module via the `addParams` option are set as an empty Groovy map [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L5); this in turn means `options.args` will be set to empty by default in the module file too. This is an intentional design choice and allows us to implement well-written sub-workflows composed of a chain of tools that by default run with the bare minimum parameter set for any given tool in order to make it much easier to share across pipelines and to provide the flexibility for users and developers to customise any non-mandatory arguments.

When including the sub-workflow above in the main pipeline workflow we use the same `include` statement, however, we now have the ability to overwrite options for each of the tools in the sub-workflow including the [`align_options`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L225) variable that will be used specifically to overwrite the optional arguments passed to the `STAR_ALIGN` module. In this case, the options to be provided to `STAR_ALIGN` have been assigned sensible defaults by the developer(s) in the pipeline's [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L70-L74) and can be accessed and customised in the [workflow context](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L201-L204) too before eventually passing them to the sub-workflow as a Groovy map called `star_align_options`. These options will then be propagated from `workflow -> sub-workflow -> module`.

As mentioned at the beginning of this section it may also be necessary for users to overwrite the options passed to modules to be able to customise specific aspects of the way in which a particular tool is executed by the pipeline. Given that all of the default module options are stored in the pipeline's `modules.config` as a [`params` variable](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L24-L25) it is also possible to overwrite any of these options via a custom config file.

Say for example we want to append an additional, non-mandatory parameter (i.e. `--outFilterMismatchNmax 16`) to the arguments passed to the `STAR_ALIGN` module. Firstly, we need to copy across the default `args` specified in the [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L71) and create a custom config file that is a composite of the default `args` as well as the additional options you would like to provide. This is very important because Nextflow will overwrite the default value of `args` that you provide via the custom config.

As you will see in the example below, we have:

* appended `--outFilterMismatchNmax 16` to the default `args` used by the module.
* changed the default `publish_dir` value to where the files will eventually be published in the main results directory.
* appended `'bam':''` to the default value of `publish_files` so that the BAM files generated by the process will also be saved in the top-level results directory for the module. Note: `'out':'log'` means any file/directory ending in `out` will now be saved in a separate directory called `my_star_directory/log/`.

```nextflow
params {
    modules {
        'star_align' {
            args          = "--quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outFilterMismatchNmax 16"
            publish_dir   = "my_star_directory"
            publish_files = ['out':'log', 'tab':'log', 'bam':'']
        }
    }
}
```

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
# nf-core/mag: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Quality control](#quality-control) of input reads - trimming and contaminant removal
* [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
* [Assembly](#assembly) of trimmed reads
* [Protein-coding gene prediction](#gene-prediction) of assemblies
* [Binning](#binning) of assembled contigs
* [Taxonomic classification of binned genomes](#taxonomic-classification-of-binned-genomes)
* [Genome annotation of binned genomes](#genome-annotation-of-binned-genomes)
* [Additional summary for binned genomes](#additional-summary-for-binned-genomes)
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Note that when specifying the parameter `--coassemble_group`, for the corresponding output filenames/directories of the assembly or downsteam processes the group ID, or more precisely the term `group-[group_id]`, will be used instead of the sample ID.

## Quality control

These steps trim away the adapter sequences present in input reads, trims away bad quality bases and sicard reads that are too short.
It also removes host contaminants and sequencing controls, such as PhiX or the Lambda phage.
FastQC is run for visualising the general quality metrics of the sequencing runs before and after trimming.

<!-- TODO: Add example MultiQC plots generated for this pipeline -->

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `QC_shortreads/fastqc/`
    * `[sample]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
    * `[sample].trimmed_[1/2]_fastqc.html`: FastQC report, containing quality metrics for trimmed and, if specified, filtered read files

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the results folder and part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `QC_shortreads/fastp/[sample]/`
    * `fastp.html`: Interactive report
    * `fastp.json`: Report in json format

</details>

### Remove PhiX sequences from short reads

The pipeline uses bowtie2 to map the reads against PhiX and removes mapped reads.

<details markdown="1">
<summary>Output files</summary>

* `QC_shortreads/remove_phix/`
    * `[sample].phix_removed.bowtie2.log`: Contains a brief log file indicating how many reads have been retained.

</details>

### Host read removal

The pipeline uses bowtie2 to map short reads against the host reference genome specified with `--host_genome` or `--host_fasta` and removes mapped reads. The information about discarded and retained reads is also included in the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

* `QC_shortreads/remove_host/`
    * `[sample].host_removed.bowtie2.log`: Contains the bowtie2 log file indicating how many reads have been mapped as well as a file listing the read ids of discarded reads.

</details>

### Remove Phage Lambda sequences from long reads

The pipeline uses Nanolyse to map the reads against the Lambda phage and removes mapped reads.

<details markdown="1">
<summary>Output files</summary>

* `QC_longreads/NanoLyse/`
    * `[sample]_nanolyse.log`: Contains a brief log file indicating how many reads have been retained.

</details>

### Filtlong and porechop

The pipeline uses filtlong and porechop to perform quality control of the long reads that are eventually provided with the TSV input file.

No direct host read removal is performed for long reads.
However, since within this pipeline filtlong uses a read quality based on k-mer matches to the already filtered short reads, reads not overlapping those short reads might be discarded.
The lower the parameter `--longreads_length_weight`, the higher the impact of the read qualities for filtering.
For further documentation see the [filtlong online documentation](https://github.com/rrwick/Filtlong).

### Quality visualisation for long reads

NanoPlot is used to calculate various metrics and plots about the quality and length distribution of long reads. For more information about NanoPlot see the [online documentation](https://github.com/wdecoster/NanoPlot).

<details markdown="1">
<summary>Output files</summary>

* `QC_longreads/NanoPlot/[sample]/`
    * `raw_*.[png/html/txt]`: Plots and reports for raw data
    * `filtered_*.[png/html/txt]`: Plots and reports for filtered data

</details>

## Taxonomic classification of trimmed reads

### Kraken

Kraken2 classifies reads using a k-mer based approach as well as assigns taxonomy using a Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

* `Taxonomy/kraken2/[sample]/`
    * `kraken2.report`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
    * `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

</details>

### Centrifuge

Centrifuge is commonly used for the classification of DNA sequences from microbial samples. It uses an indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.

More information on the [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) website

<details markdown="1">
<summary>Output files</summary>

* `Taxonomy/centrifuge/[sample]/`
    * `report.txt`: Tab-delimited result file. See the [centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output) for information about the fields
    * `kreport.txt`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
    * `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

</details>

## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes. Hybrid assembly is only supported by SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

<details markdown="1">
<summary>Output files</summary>

* `Assembly/MEGAHIT/`
    * `[sample/group].contigs.fa.gz`: Compressed metagenome assembly in fasta format
    * `[sample/group].log`: Log file
    * `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
        * `MEGAHIT-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
        * `MEGAHIT-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

</details>

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

<details markdown="1">
<summary>Output files</summary>

* `Assembly/SPAdes/`
    * `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
    * `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
    * `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
    * `[sample/group].log`: Log file
    * `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
        * `SPAdes-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
        * `SPAdes-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

</details>

### SPAdesHybrid

SPAdesHybrid is a part of the [SPAdes](http://cab.spbu.ru/software/spades/) software and is used when the user provides both long and short reads.

<details markdown="1">
<summary>Output files</summary>

* `Assembly/SPAdesHybrid/`
    * `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
    * `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
    * `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
    * `[sample/group].log`: Log file
    * `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
        * `SPAdesHybrid-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
        * `SPAdesHybrid-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

</details>

### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

<details markdown="1">
<summary>Output files</summary>

* `Assembly/[assembler]/QC/[sample/group]/`
    * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
    * `quast.log`: QUAST log file
    * `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format

</details>

## Gene prediction

Protein-coding genes are predicted for each assembly.

<details markdown="1">
<summary>Output files</summary>

* `Prodigal/`
    * `[sample/group].gff`: Gene Coordinates in GFF format
    * `[sample/group].faa`: The protein translation file consists of all the proteins from all the sequences in multiple FASTA format.
    * `[sample/group].fna`: Nucleotide sequences of the predicted proteins using the DNA alphabet, not mRNA (so you will see 'T' in the output and not 'U').
    * `[sample/group]_all.txt`: Information about start positions of genes.

</details>

## Binning

### Contig sequencing depth

Sequencing depth per contig and sample is generated by `jgi_summarize_bam_contig_depths --outputDepth`. The values correspond to `(sum of exactely aligned bases) / ((contig length)-2*75)`. For example, for two reads aligned exactly with `10` and `9` bases on a 1000 bp long contig the depth is calculated by `(10+9)/(1000-2*75)` (1000bp length of contig minus 75bp from each end, which is excluded).

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/`
    * `[assembler]-[sample/group]-depth.txt.gz`: Sequencing depth for each contig and sample or group, only for short reads.

</details>

### MetaBAT2

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/MetaBAT2/`
    * `[assembler]-[sample/group].*.fa`: Genome bins retrieved from input assembly
    * `[assembler]-[sample/group].unbinned.*.fa`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported

</details>

All the files and contigs in this folder will be assessed by QUAST and BUSCO.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/MetaBAT2/discarded/`
    * `*.lowDepth.fa.gz`: Low depth contigs that are filtered by MetaBat2
    * `*.tooShort.fa.gz`: Too short contigs that are filtered by MetaBat2
    * `*.unbinned.pooled.fa.gz`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
    * `*.unbinned.remaining.fa.gz`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

</details>

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### Bin sequencing depth

For each genome bin the median sequencing depth is computed based on the corresponding contig depths given in `GenomeBinning/[assembler]-[sample/group]-depth.txt.gz`.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/`
    * `bin_depths_summary.tsv`: Summary of bin sequencing depths for all samples. Depths are available for samples mapped against the corresponding assembly, i.e. according to the mapping strategy specified with `--binning_map_mode`. Only for short reads.
    * `[assembler]-[sample/group]-binDepths.heatmap.png`: Clustered heatmap showing bin abundances of the assembly across samples. Bin depths are transformed to centered log-ratios and bins as well as samples are clustered by Euclidean distance. Again, sample depths are available according to the mapping strategy specified with `--binning_map_mode`.

</details>

### QC for metagenome assembled genomes with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/QC/QUAST/[assembler]-[bin]/`
    * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
    * `quast.log`: QUAST log file
    * `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format
* `GenomeBinning/QC/`
    * `quast_summary.tsv`: QUAST output for all bins summarized

</details>

### QC for metagenome assembled genomes with BUSCO

[BUSCO](https://busco.ezlab.org/) is a tool used to assess the completeness of a genome assembly. It is run on all the genome bins and high quality contigs obtained by MetaBAT2. By default, BUSCO is run in automated lineage selection mode in which it first tries to select the domain and then a more specific lineage based on phylogenetic placement. If available, result files for both the selected domain lineage and the selected more specific lineage are placed in the output directory. If a lineage dataset is specified already with `--busco_reference`, only results for this specific lineage will be generated.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/QC/BUSCO/`
    * `[assembler]-[bin]_busco.log`: Log file containing the standard output of BUSCO.
    * `[assembler]-[bin]_busco.err`: File containing potential error messages returned from BUSCO.
    * `short_summary.domain.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results for the selected domain when run in automated lineage selection mode. Not available for bins for which a viral lineage was selected.
    * `short_summary.specific_lineage.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results in case a more specific lineage than the domain could be selected or for the lineage provided via `--busco_reference`.
    * `[assembler]-[bin]_buscos.[lineage].fna.gz`: Nucleotide sequence of all identified BUSCOs for used lineages (domain or specific).
    * `[assembler]-[bin]_buscos.[lineage].faa.gz`: Aminoacid sequence of all identified BUSCOs for used lineages (domain or specific).
    * `[assembler]-[bin]_prodigal.gff`: Genes predicted with Prodigal.

</details>

If the parameter `--save_busco_reference` is set, additionally the used BUSCO lineage datasets are stored in the output directy.

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/QC/BUSCO/`
    * `busco_downloads/`: All files and lineage datasets downloaded by BUSCO when run in automated lineage selection mode. (Can currently not be used to reproduce analysis, see the [nf-core/mag website documentation](https://nf-co.re/mag/usage#reproducibility) how to achieve reproducible BUSCO results).
    * `reference/*.tar.gz`: BUSCO reference lineage dataset that was provided via `--busco_reference`.

</details>

Besides the reference files or output files created by BUSCO, the following summary files will be generated:

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/QC/`
    * `busco_summary.tsv`: A summary table of the BUSCO results, with % of marker genes found. If run in automated lineage selection mode, both the results for the selected domain and for the selected more specific lineage will be given, if available.

</details>

## Taxonomic classification of binned genomes

### CAT

[CAT](https://github.com/dutilh/CAT) is a toolkit for annotating contigs and bins from metagenome-assembled-genomes. The MAG pipeline uses CAT to assign taxonomy to genome bins based on the taxnomy of the contigs.

<details markdown="1">
<summary>Output files</summary>

* `Taxonomy/CAT/[assembler]/`
    * `[assembler]-[sample/group].ORF2LCA.names.txt.gz`: Tab-delimited files containing the lineage of each contig, with full lineage names
    * `[assembler]-[sample/group].bin2classification.names.txt.gz`: Taxonomy classification of the genome bins, with full lineage names
* `Taxonomy/CAT/[assembler]/raw/`
    * `[assembler]-[sample/group].concatenated.predicted_proteins.faa.gz`: Predicted protein sequences for each genome bin, in fasta format
    * `[assembler]-[sample/group].concatenated.predicted_proteins.gff.gz`: Predicted protein features for each genome bin, in gff format
    * `[assembler]-[sample/group].ORF2LCA.txt.gz`: Tab-delimited files containing the lineage of each contig
    * `[assembler]-[sample/group].bin2classification.txt.gz`: Taxonomy classification of the genome bins
    * `[assembler]-[sample/group].log`: Log files

</details>

If the parameters `--cat_db_generate` and `--save_cat_db` are set, additionally the generated CAT database is stored:

<details markdown="1">
<summary>Output files</summary>

* `Taxonomy/CAT/CAT_prepare_*.tar.gz`: Generated and used CAT database.

</details>

### GTDB-Tk

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is a toolkit for assigning taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). nf-core/mag uses GTDB-Tk to classify binned genomes which satisfy certain quality criteria (i.e. completeness and contamination assessed with the BUSCO analysis).

<details markdown="1">
<summary>Output files</summary>

* `Taxonomy/GTDB-Tk/[assembler]/[sample/group]/`
    * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.summary.tsv`: Classifications for bacterial and archaeal genomes (see the [GTDB-Tk documentation for details](https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html).
    * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.classify.tree.gz`: Reference tree in Newick format containing query genomes placed with pplacer.
    * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.markers_summary.tsv`: A summary of unique, duplicated, and missing markers within the 120 bacterial marker set, or the 122 archaeal marker set for each submitted genome.
    * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.msa.fasta.gz`: FASTA file containing MSA of submitted and reference genomes.
    * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.filtered.tsv`: A list of genomes with an insufficient number of amino acids in MSA.
    * `gtdbtk.[assembler]-[sample/group].*.log`: Log files.
    * `gtdbtk.[assembler]-[sample/group].failed_genomes.tsv`: A list of genomes for which the GTDB-Tk analysis failed, e.g. because Prodigal could not detect any genes.
* `Taxonomy/GTDB-Tk/gtdbtk_summary.tsv`: A summary table of the GTDB-Tk classification results for all bins, also containing bins which were discarded based on the BUSCO QC, which were filtered out by GTDB-Tk ((listed in `*.filtered.tsv`) or for which the analysis failed (listed in `*.failed_genomes.tsv`).

</details>

## Genome annotation of binned genomes

### Prokka

Whole genome annotation is the process of identifying features of interest in a set of genomic DNA sequences, and labelling them with useful information. [Prokka](https://github.com/tseemann/prokka) is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files.

<details markdown="1">
<summary>Output files</summary>

* `Prokka/[assembler]/[bin]/`
    * `[bin].gff`: annotation in GFF3 format, containing both sequences and annotations
    * `[bin].gbk`: annotation in GenBank format, containing both sequences and annotations
    * `[bin].fna`: nucleotide FASTA file of the input contig sequences
    * `[bin].faa`: protein FASTA file of the translated CDS sequences
    * `[bin].ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    * `[bin].sqn`: an ASN1 format "Sequin" file for submission to Genbank
    * `[bin].fsa`: nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file
    * `[bin].tbl`: feature Table file, used by "tbl2asn" to create the .sqn file
    * `[bin].err`: unacceptable annotations - the NCBI discrepancy report.
    * `[bin].log`: contains all the output that Prokka produced during its run
    * `[bin].txt`: statistics relating to the annotated features found
    * `[bin].tsv`: tab-separated file of all features (locus_tag, ftype, len_bp, gene, EC_number, COG, product)

</details>

## Additional summary for binned genomes

<details markdown="1">
<summary>Output files</summary>

* `GenomeBinning/bin_summary.tsv`: Summary of bin sequencing depths together with BUSCO, QUAST and GTDB-Tk results, if at least one of the later was generated.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.tsv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
<!--
# nf-core/mag pull request

Many thanks for contributing to nf-core/mag!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the dev branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/nf-core/mag/tree/master/.github/CONTRIBUTING.md)
-->
<!-- markdownlint-disable ul-indent -->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
    - [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/nf-core/mag/tree/master/.github/CONTRIBUTING.md)
    - [ ] If necessary, also make a PR on the nf-core/mag _branch_ on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
- [ ] Make sure your code lints (`nf-core lint`).
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
# nf-core/mag: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf-core/mag.

We try to manage the required tasks for nf-core/mag using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

> If you need help using or modifying nf-core/mag then the best place to ask is on the nf-core Slack [#mag](https://nfcore.slack.com/channels/mag) channel ([join our Slack here](https://nf-co.re/join/slack)).

## Contribution workflow

If you'd like to write some code for nf-core/mag, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [nf-core/mag issues](https://github.com/nf-core/mag/issues) to avoid duplicating work
    * If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/mag repository](https://github.com/nf-core/mag) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core schema build` and add any new parameters to the pipeline JSON schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

### Pipeline tests

Each `nf-core` pipeline should be set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

* On your own fork, make a new branch `patch` based on `upstream/master`.
* Fix the bug, and bump version (X.Y.Z+1).
* A PR should be made on `master` from patch to directly this particular bug.

## Getting help

For further information/help, please consult the [nf-core/mag documentation](https://nf-co.re/mag/usage) and don't hesitate to get in touch on the nf-core Slack [#mag](https://nfcore.slack.com/channels/mag) channel ([join our Slack here](https://nf-co.re/join/slack)).

## Pipeline contribution conventions

To make the nf-core/mag code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new flags/options to `nextflow.config` with a default (see below).
5. Add any new flags/options to `nextflow_schema.json` with help text (with `nf-core schema build`).
6. Add any new flags/options to the help message (for integer/text parameters, print to help the corresponding `nextflow.config` parameter).
7. Add sanity checks for all relevant parameters.
8. Add any new software to the `scrape_software_versions.py` script in `bin/` and the version command to the `scrape_software_versions` process in `main.nf`.
9. Do local tests that the new code works properly and as expected.
10. Add a new test command in `.github/workflow/ci.yml`.
11. If applicable add a [MultiQC](https://https://multiqc.info/) module.
12. Update MultiQC config `assets/multiqc_config.yaml` so relevant suffixes, name clean up, General Statistics Table column order, and module figures are in the right order.
13. Optional: Add any descriptions of MultiQC report sections and output files to `docs/output.md`.

### Default values

Parameters should be initialised / defined with default values in `nextflow.config` under the `params` scope.

Once there, use `nf-core schema build` to add to `nextflow_schema.json`.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. A nf-core standard set of labels that should be followed where possible can be seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config), which has the default process as a single core-process, and then different levels of multi-core configurations for increasingly large memory requirements defined with standardised labels.

The process resources can be passed on to the tool dynamically within the process with the `${task.cpu}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

* initial process channel: `ch_output_from_<process>`
* intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of nextflow in the pipeline with: `nf-core bump-version --nextflow . [min-nf-version]`

### Software version reporting

If you add a new tool to the pipeline, please ensure you add the information of the tool to the `get_software_version` process.

Add to the script block of the process, something like the following:

```bash
<YOUR_TOOL> --version &> v_<YOUR_TOOL>.txt 2>&1 || true
```

or

```bash
<YOUR_TOOL> --help | head -n 1 &> v_<YOUR_TOOL>.txt 2>&1 || true
```

You then need to edit the script `bin/scrape_software_versions.py` to:

1. Add a Python regex for your tool's `--version` output (as in stored in the `v_<YOUR_TOOL>.txt` file), to ensure the version is reported as a `v` and the version number e.g. `v2.1.1`
2. Add a HTML entry to the `OrderedDict` for formatting in MultiQC.

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).
---
name: Feature request
about: Suggest an idea for the nf-core/mag pipeline
labels: enhancement
---

<!--
# nf-core/mag feature request

Hi there!

Thanks for suggesting a new feature for the pipeline!
Please delete this text and anything that's not relevant from the template below:
-->

## Is your feature request related to a problem? Please describe

<!-- A clear and concise description of what the problem is. -->

<!-- e.g. [I'm always frustrated when ...] -->

## Describe the solution you'd like

<!-- A clear and concise description of what you want to happen. -->

## Describe alternatives you've considered

<!-- A clear and concise description of any alternative solutions or features you've considered. -->

## Additional context

<!-- Add any other context about the feature request here. -->
---
name: Bug report
about: Report something that is broken or incorrect
labels: bug
---

<!--
# nf-core/mag bug report

Hi there!

Thanks for telling us about a problem with the pipeline.
Please delete this text and anything that's not relevant from the template below:
-->

## Check Documentation

I have checked the following places for your error:

- [ ] [nf-core website: troubleshooting](https://nf-co.re/usage/troubleshooting)
- [ ] [nf-core/mag pipeline documentation](https://nf-co.re/mag/usage)

## Description of the bug

<!-- A clear and concise description of what the bug is. -->

## Steps to reproduce

Steps to reproduce the behaviour:

1. Command line: <!-- [e.g. `nextflow run ...`] -->
2. See error: <!-- [Please provide your error message] -->

## Expected behaviour

<!-- A clear and concise description of what you expected to happen. -->

## Log files

Have you provided the following extra information/files:

- [ ] The command used to run the pipeline
- [ ] The `.nextflow.log` file <!-- this is a hidden file in the directory where you launched the pipeline -->

## System

- Hardware: <!-- [e.g. HPC, Desktop, Cloud...] -->
- Executor: <!-- [e.g. slurm, local, awsbatch...] -->
- OS: <!-- [e.g. CentOS Linux, macOS, Linux Mint...] -->
- Version <!-- [e.g. 7, 10.13.6, 18.3...] -->

## Nextflow Installation

- Version: <!-- [e.g. 21.04.0] -->

## Container engine

- Engine: <!-- [e.g. Conda, Docker, Singularity, Podman, Shifter or Charliecloud] -->
- version: <!-- [e.g. 1.0.0] -->

## Additional context

<!-- Add any other context about the problem here. -->
