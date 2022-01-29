For contributions
see
[contributors](https://github.com/vcflib/vcflib/graphs/contributors)
and
[commits](https://github.com/vcflib/vcflib/commits/master).

## ChangeLog v1.0.3 (20220122)

This is a maintenance release of vcflib.

+ Merge intervaltree changes (thanks @jnunn and @timmassingham)
+ Built with gcc-11
+ Fix issue #251 hapLrt: fix segfault when accessing genotype field. (thanks @mphschmitt)
+ Fix vcfflatten: fix segfault when no 'AF' field is present (#47, thanks @mphschmitt)
+ Fixes on vcfnulldotslashdot #310 (thanks @WinterFor)
+ Fix issue #301: Replace raw pointer usage with std::unique_ptr #306 (thanks @Glebanister)
+ Fix man page installation #321 (thanks @alexreg)
+ Use `guix shell` instead of `guix environment` for development
+ Regenerated online docs
+ README: add matrix badge (removed gitter badge)

## ChangeLog v1.0.2 (20210104)

This is a maintenance release of vcflib, mostly improving the build
system, CI and generating markdown docs as well as man pages.

+ Removed tabixpp and htslib source dependencies, i.e., we are now using
  the distro provided libraries and include files through pkg-config.
  See also the [README](README.md#build-from-source)
+ Removed the tabixpp+htslib git submodules
+ Generalise and document the cmake build system
+ Added tests to the cmake build system and build instructions to README
+ Added support for ARM64 and PowerPC, see #292 (thanks @genisysram and @mr-c)
+ Added github actions for the issue tracker
+ Added githum CI
+ Updated header files in src with copyright/license info, see #16
+ Created markdown [docs](./doc/vcflib.md) and [man pages](./man/) for
  all utilities. Created a script bin2md for markdown generation and
  use pandoc for the man page generation.

## Older changes

For older changes view the git [log](https://github.com/vcflib/vcflib/commits/master).
# vcflib

### A C++ library for parsing and manipulating VCF files.

[![Github-CI](https://github.com/vcflib/vcflib/workflows/CI/badge.svg)](https://github.com/vcflib/vcflib/actions?query=workflow%3ACI) [![Travis-CI](https://travis-ci.com/vcflib/vcflib.svg?branch=master)](https://travis-ci.com/github/vcflib/vcflib) [![AnacondaBadge](https://anaconda.org/bioconda/vcflib/badges/installer/conda.svg)](https://anaconda.org/bioconda/vcflib) [![DL](https://anaconda.org/bioconda/vcflib/badges/downloads.svg)](https://anaconda.org/bioconda/vcflib) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-vcflib-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-vcflib-brightgreen.svg)](https://www.gnu.org/software/guix/packages/V/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/libvcflib-dev/version.svg)](https://packages.debian.org/testing/libvcflib-dev) [![C++0x](https://img.shields.io/badge/Language-C++0x-steelblue.svg)](https://www.cprogramming.com/c++11/what-is-c++0x.html) [![Chat on Matrix](https://matrix.to/img/matrix-badge.svg)](https://matrix.to/#/#vcflib:matrix.org)

Vcflib and related tools are the workhorses in bioinformatics for processing the VCF variant calling format. See

Vcflib and tools for processing the VCF variant call format;
Erik Garrison, Zev N. Kronenberg, Eric T. Dawson, Brent S. Pedersen, Pjotr Prins;
doi: https://doi.org/10.1101/2021.05.21.445151

## overview

The [Variant Call Format
(VCF)](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)
is a flat-file, tab-delimited textual format that
describes reference-indexed variations between individuals.  VCF
provides a common interchange format for the description of variation
in individuals and populations of samples, and has become the
*de facto* standard reporting format for a wide array of genomic
variant detectors.

vcflib provides methods to manipulate and interpret sequence variation
described by VCF.  It is both:

 * an API for parsing and operating on records of genomic variation as
   it can be described by the VCF format
 * a collection of command-line utilities for executing complex
   manipulations on VCF files

vclib is both a library (with an API) and a collection of useful
tools. The API provides a quick and extremely permissive method to
read and write VCF files.  Extensions and applications of the library
provided in the included utilities (*.cpp) comprise the vast bulk of
the library's utility.

---

Short index:

- [Install](#INSTALL)
- [Usage](#USAGE)
- [TOOLS](#TOOLS)
  * [Filter](#filter)
  * [Metrics](#metrics)
  * [Phenotype](#phenotype)
  * [Genotype](#genotype)
  * [Transformation](#transformation)
  * [Statistics](#statistics)
  * [Scripts](#scripts)
- [Link library](#link-library)
- [Build from source](#build-from-source)
- [Development](#Development)
- [LICENSE](#LICENSE)
- [Credit work](#Credit)

---

## INSTALL

For latest updates see [RELEASE NOTES](./RELEASE_NOTES.md).

### [Bioconda](https://bioconda.github.io/user/install.html)

Conda installs in user land without root access

```sh
conda install -c bioconda vcflib
```

### [Homebrew](https://brew.sh)

Homebrew installs on Linux and Mac OSX

```sh
brew install brewsci/bio/vcflib
```

### [Debian](https://debian.org/)

For Debian and Ubuntu

```sh
apt-get install libvcflib-tools libvcflib-dev
```

### [GNU Guix](https://guix.gnu.org/)

We develop against guix and vcflib is packaged as

```sh
guix package -i vcflib
```

See also the Guix shell below.

## USAGE

Users are encouraged to drive the utilities in the library in a
streaming fashion, using Unix pipes to fully utilize resources on
multi-core systems.  Piping provides a convenient method to interface
with other libraries (vcf-tools, BedTools, GATK, htslib,
[bio-vcf](https://github.com/vcflib/bio-vcf), bcftools,
[freebayes](https://github.com/freebayes)) which interface via VCF
files, allowing the composition of an immense variety of processing
functions. Examples can be found in the scripts,
e.g. [script](./scripts/vcfgtcompare.sh).


# TOOLS

<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->


<!--
  Created with ./scripts/bin2md.rb --index
-->


## filter

| filter command | description |
| :-------------- | :---------- |
 | [vcfuniq](./doc/vcfuniq.md) | List unique genotypes. Like GNU uniq, but for VCF records. Remove records which have the same position, ref, and alt as the previous record. |
 | [vcfuniqalleles](./doc/vcfuniqalleles.md) | List unique alleles For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files. |
 | [vcffilter](./doc/vcffilter.md) | VCF filter the specified vcf file using the set of filters |

## metrics

| metrics command | description |
| :-------------- | :---------- |
 | [vcfcheck](./doc/vcfcheck.md) | Validate integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file. |
 | [vcfhethomratio](./doc/vcfhethomratio.md) | Generates the het/hom ratio for each individual in the file |
 | [vcfhetcount](./doc/vcfhetcount.md) | Calculate the heterozygosity rate: count the number of alternate alleles in heterozygous genotypes in all records in the vcf file |
 | [vcfdistance](./doc/vcfdistance.md) | Adds a tag to each variant record which indicates the distance to the nearest variant. (defaults to BasesToClosestVariant if no custom tag name is given. |
 | [vcfentropy](./doc/vcfentropy.md) | Annotate VCF records with the Shannon entropy of flanking sequence. Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight, EntropyCenter, which are the entropies of the sequence of the given window size to the left, right, and center of the record. Also adds EntropyRef and EntropyAlt for each alt. |

## phenotype

| phenotype command | description |
| :-------------- | :---------- |
 | [permuteGPAT++](./doc/permuteGPAT++.md) | **permuteGPAT++** is a method for adding empirical p-values to a GPAT++ score. |

## genotype

| genotype command | description |
| :-------------- | :---------- |
 | [normalize-iHS](./doc/normalize-iHS.md) | normalizes iHS or XP-EHH scores. |
 | [hapLrt](./doc/hapLrt.md) | HapLRT is a likelihood ratio test for haplotype lengths. The lengths are modeled with an exponential distribution. The sign denotes if the target has longer haplotypes (1) or the background (-1). |
 | [abba-baba](./doc/abba-baba.md) | **abba-baba** calculates the tree pattern for four indviduals. This tool assumes reference is ancestral and ignores non **abba-baba** sites. The output is a boolian value: 1 = true , 0 = false for abba and baba. the tree argument should be specified from the most basal taxa to the most derived. |

## transformation

| transformation command | description |
| :-------------- | :---------- |
 | [vcfinfo2qual](./doc/vcfinfo2qual.md) | Sets QUAL from info field tag keyed by [key]. The VCF file may be omitted and read from stdin. The average of the field is used if it contains multiple values. |
 | [vcfsamplediff](./doc/vcfsamplediff.md) | Establish putative somatic variants using reported differences between germline and somatic samples. Tags each record where the listed sample genotypes differ with <tag>. The first sample is assumed to be germline, the second somatic. Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples. |
 | [vcfaddinfo](./doc/vcfaddinfo.md) | Adds info fields from the second file which are not present in the first vcf file. |
 | [vcfremoveaberrantgenotypes](./doc/vcfremoveaberrantgenotypes.md) | strips samples which are homozygous but have observations implying heterozygosity. Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO). |
 | [vcfglxgt](./doc/vcfglxgt.md) | Set genotypes using the maximum genotype likelihood for each sample. |
 | [dumpContigsFromHeader](./doc/dumpContigsFromHeader.md) | Dump contigs from header |
 | [vcfevenregions](./doc/vcfevenregions.md) | Generates a list of regions, e.g. chr20:10..30 using the variant density information provided in the VCF file to ensure that the regions have even numbers of variants. This can be use to reduce the variance in runtime when dividing variant detection or genotyping by genomic coordinates. |
 | [vcfcat](./doc/vcfcat.md) | Concatenates VCF files |
 | [vcfannotategenotypes](./doc/vcfannotategenotypes.md) | Examine genotype correspondence. Annotate genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant. |
 | [vcfafpath](./doc/vcfafpath.md) | Display genotype paths |
 | [vcfclassify](./doc/vcfclassify.md) | Creates a new VCF where each variant is tagged by allele class: snp, ts/tv, indel, mnp |
 | [vcfallelicprimitives](./doc/vcfallelicprimitives.md) | If multiple allelic primitives (gaps or mismatches) are specified in a single VCF record, split the record into multiple lines, but drop all INFO fields. Does not handle genotypes (yet). MNPs are split into multiple SNPs unless the -m flag is provided. Records generated by splits have th |
 | [vcfqual2info](./doc/vcfqual2info.md) | Puts QUAL into an info field tag keyed by [key]. |
 | [vcfcreatemulti](./doc/vcfcreatemulti.md) | If overlapping alleles are represented across multiple records, merge them into a single record. Currently only for indels. |
 | [vcfgeno2alleles](./doc/vcfgeno2alleles.md) | modifies the genotypes field to provide the literal alleles rather than indexes |
 | [vcfsample2info](./doc/vcfsample2info.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcfld](./doc/vcfld.md) | Compute LD |
 | [vcfnumalt](./doc/vcfnumalt.md) | outputs a VCF stream where NUMALT has been generated for each record using sample genotypes |
 | [vcfstreamsort](./doc/vcfstreamsort.md) | Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart. |
 | [vcfinfosummarize](./doc/vcfinfosummarize.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcflength](./doc/vcflength.md) | Add length info field |
 | [vcfkeepgeno](./doc/vcfkeepgeno.md) | Reduce file size by removing FORMAT fields not listed on the command line from sample specifications in the output |
 | [vcfcombine](./doc/vcfcombine.md) | Combine VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted. |
 | [vcfprimers](./doc/vcfprimers.md) | For each VCF record, extract the flanking sequences, and write them to stdout as FASTA records suitable for alignment. |
 | [vcfflatten](./doc/vcfflatten.md) | Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin. |
 | [vcf2dag](./doc/vcf2dag.md) | Modify VCF to be able to build a directed acyclic graph (DAG) |
 | [vcfcleancomplex](./doc/vcfcleancomplex.md) | Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change. |
 | [vcfbreakmulti](./doc/vcfbreakmulti.md) | If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields. |
 | [vcfindex](./doc/vcfindex.md) | Adds an index number to the INFO field (id=position) |
 | [vcfkeepinfo](./doc/vcfkeepinfo.md) | To decrease file size remove INFO fields not listed on the command line |
 | [vcfgeno2haplo](./doc/vcfgeno2haplo.md) | Convert genotype-based phased alleles within --window-size into haplotype alleles. Will break haplotype construction when encountering non-phased genotypes on input. |
 | [vcfintersect](./doc/vcfintersect.md) | VCF set analysis |
 | [vcfannotate](./doc/vcfannotate.md) | Intersect the records in the VCF file with targets provided in a BED file. Intersections are done on the reference sequences in the VCF file. If no VCF filename is specified on the command line (last argument) the VCF read from stdin. |
 | [smoother](./doc/smoother.md) | smoothes is a method for window smoothing many of the GPAT++ formats. |
 | [vcf2fasta](./doc/vcf2fasta.md) | Generates sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy]. Each sequence in the fasta file is named using the same pattern used for the file name, allowing them to be combined. |
 | [vcfsamplenames](./doc/vcfsamplenames.md) | List sample names |
 | [vcfleftalign](./doc/vcfleftalign.md) | Left-align indels and complex variants in the input using a pairwise ref/alt alignment followed by a heuristic, iterative left realignment process that shifts indel representations to their absolute leftmost (5') extent. |
 | [vcfglbound](./doc/vcfglbound.md) | Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max. |
 | [vcfcommonsamples](./doc/vcfcommonsamples.md) | Generates each record in the first file, removing samples not present in the second |
 | [vcfecho](./doc/vcfecho.md) | Echo VCF to stdout (simple demo) |
 | [vcfkeepsamples](./doc/vcfkeepsamples.md) | outputs each record in the vcf file, removing samples not listed on the command line |
 | [vcf2tsv](./doc/vcf2tsv.md) | Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table. Specifying -g will output one line per sample with genotype information. When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index |
 | [vcfoverlay](./doc/vcfoverlay.md) | Overlay records in the input vcf files with order as precedence. |
 | [vcfgenosamplenames](./doc/vcfgenosamplenames.md) | Get samplenames |
 | [vcfremovesamples](./doc/vcfremovesamples.md) | outputs each record in the vcf file, removing samples listed on the command line |
 | [vcfremap](./doc/vcfremap.md) | For each alternate allele, attempt to realign against the reference with lowered gap open penalty. If realignment is possible, adjust the cigar and reference/alternate alleles. Observe how different alignment parameters, including context and entropy-dependent ones, influence variant classification and interpretation. |
 | [vcffixup](./doc/vcffixup.md) | Generates a VCF stream where AC and NS have been generated for each record using sample genotypes |

## statistics

| statistics command | description |
| :-------------- | :---------- |
 | [vcfgenosummarize](./doc/vcfgenosummarize.md) | Adds summary statistics to each record summarizing qualities reported in called genotypes. Uses: RO (reference observation count), QR (quality sum reference observations) AO (alternate observation count), QA (quality sum alternate observations) |
 | [vcfcountalleles](./doc/vcfcountalleles.md) | Count alleles |
 | [meltEHH](./doc/meltEHH.md) |  |
 | [genotypeSummary](./doc/genotypeSummary.md) | Generates a table of genotype counts. Summarizes genotype counts for bi-allelic SNVs and indel |
 | [vcfrandomsample](./doc/vcfrandomsample.md) | Randomly sample sites from an input VCF file, which may be provided as stdin. Scale the sampling probability by the field specified in KEY. This may be used to provide uniform sampling across allele frequencies, for instance. |
 | [pVst](./doc/pVst.md) | **pVst** calculates vst, a measure of CNV stratification. |
 | [vcfrandom](./doc/vcfrandom.md) | Generate a random VCF file |
 | [segmentFst](./doc/segmentFst.md) | **segmentFst** creates genomic segments (bed file) for regions with high wcFst |
 | [sequenceDiversity](./doc/sequenceDiversity.md) | The **sequenceDiversity** program calculates two popular metrics of haplotype diversity: pi and extended haplotype homozygoisty (eHH). Pi is calculated using the Nei and Li 1979 formulation. eHH a convenient way to think about haplotype diversity. When eHH = 0 all haplotypes in the window are unique and when eHH = 1 all haplotypes in the window are identical. |
 | [segmentIhs](./doc/segmentIhs.md) | Creates genomic segments (bed file) for regions with high wcFst |
 | [vcfgenotypes](./doc/vcfgenotypes.md) | Report the genotypes for each sample, for each variant in the VCF. Convert the numerical represenation of genotypes provided by the GT field to a human-readable genotype format. |
 | [vcfaltcount](./doc/vcfaltcount.md) | count the number of alternate alleles in all records in the vcf file |
 | [plotHaps](./doc/plotHaps.md) | **plotHaps** provides the formatted output that can be used with 'bin/plotHaplotypes.R'. |
 | [vcfsitesummarize](./doc/vcfsitesummarize.md) | Summarize by site |
 | [vcfgenotypecompare](./doc/vcfgenotypecompare.md) | adds statistics to the INFO field of the vcf file describing the amount of discrepancy between the genotypes (GT) in the vcf file and the genotypes reported in the <other-genotype-tag>. use this after vcfannotategenotypes to get correspondence statistics for two vcfs. |
 | [vcfstats](./doc/vcfstats.md) | Prints statistics about variants in the input VCF file. |
 | [wcFst](./doc/wcFst.md) | **wcFst** is Weir & Cockerham's Fst for two populations. Negative values are VALID, they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. Specifically **wcFst** uses equations 1,2,3,4. |
 | [permuteSmooth](./doc/permuteSmooth.md) | **permuteSmooth** is a method for adding empirical p-values smoothed wcFst scores. |
 | [bFst](./doc/bFst.md) | **bFst** is a Bayesian approach to Fst. Importantly **bFst** accounts for genotype uncertainty in the model using genotype likelihoods. For a more detailed description see: `A Bayesian approach to inferring population structure from dominant markers' by Holsinger et al. Molecular Ecology Vol 11, issue 7 2002. The likelihood function has been modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population's allele frequency, and Fst. |
 | [vcfroc](./doc/vcfroc.md) | Generates a pseudo-ROC curve using sensitivity and specificity estimated against a putative truth set. Thresholding is provided by successive QUAL cutoffs. |
 | [vcfparsealts](./doc/vcfparsealts.md) | Alternate allele parsing method. This method uses pairwise alignment of REF and ALTs to determine component allelic primitives for each alternate allele. |
 | [pFst](./doc/pFst.md) | **pFst** is a probabilistic approach for detecting differences in allele frequencies between two populations. |
 | [iHS](./doc/iHS.md) | **iHS** calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014). |
 | [popStats](./doc/popStats.md) | General population genetic statistics for each SNP |

See also [vcflib.md](./doc/vcflib.md).

## scripts

The vcflib source repository contains a number of additional scripts.
Click on the link to see the source code.

| script | description |
| :-------------- | :---------- |
| [vcfclearinfo](./scripts/vcfclearinfo) | clear INFO field |
| [vcfqualfilter](./scripts/vcfqualfilter) | quality filter |
| [vcfnulldotslashdot](./scripts/vcfnulldotslashdot) | rewrite null genotypes to ./. |
| [vcfprintaltdiscrepancy.r](./scripts/vcfprintaltdiscrepancy.r) | show ALT discrepancies in a table |
| [vcfremovenonATGC](./scripts/vcfremovenonATGC) | remove non-nucleotides in REF or ALT |
| [plotSmoothed.R](./scripts/plotSmoothed.R) | smooth plot of wcFst, pFst or abba-baba |
| [vcf_strip_extra_headers](./scripts/vcf_strip_extra_headers) | strip headers |
| [plotHapLrt.R](./scripts/plotHapLrt.R) | plot results of pFst |
| [vcfbiallelic](./scripts/vcfbiallelic) | remove anything that is not biallelic |
| [vcfsort](./scripts/vcfsort) | sort VCF using shell script |
| [vcfnosnps](./scripts/vcfnosnps) | remove SNPs |
| [vcfmultiwayscripts](./scripts/vcfmultiwayscripts) | more multiway comparisons |
| [vcfgtcompare.sh](./scripts/vcfgtcompare.sh) | annotates records in the first file with genotypes and sites from the second |
| [plotPfst.R](./scripts/plotPfst.R) | plot pFst |
| [vcfregionreduce_and_cut](./scripts/vcfregionreduce_and_cut) | reduce, gzip, and tabix |
| [plotBfst.R](./scripts/plotBfst.R) | plot results of pFst |
| [vcfnobiallelicsnps](./scripts/vcfnobiallelicsnps) | remove biallelic SNPs |
| [vcfindels](./scripts/vcfindels) | show INDELS |
| [vcfmultiway](./scripts/vcfmultiway) | multiway comparison |
| [vcfregionreduce](./scripts/vcfregionreduce) | reduce VCFs using a BED File, gzip them up and create tabix index |
| [vcfprintaltdiscrepancy.sh](./scripts/vcfprintaltdiscrepancy.sh) | runner |
| [vcfclearid](./scripts/vcfclearid) | clear ID field |
| [vcfcomplex](./scripts/vcfcomplex) | remove all SNPs but keep SVs |
| [vcffirstheader](./scripts/vcffirstheader) | show first header |
| [plotXPEHH.R](./scripts/plotXPEHH.R) | plot XPEHH |
| [vcfregionreduce_pipe](./scripts/vcfregionreduce_pipe) | reduce, gzip and tabix in a pipe |
| [vcfplotaltdiscrepancy.sh](./scripts/vcfplotaltdiscrepancy.sh) | plot ALT discrepancy runner |
| [vcfplottstv.sh](./scripts/vcfplottstv.sh) | runner |
| [vcfnoindels](./scripts/vcfnoindels) | remove INDELs |
| [bgziptabix](./scripts/bgziptabix) | runs bgzip on the input and tabix indexes the result |
| [plotHaplotypes.R](./scripts/plotHaplotypes.R) | plot results |
| [vcfplotsitediscrepancy.r](./scripts/vcfplotsitediscrepancy.r) | plot site discrepancy |
| [vcfindelproximity](./scripts/vcfindelproximity) | show SNPs around an INDEL |
| [bed2region](./scripts/bed2region) | convert VCF CHROM column in VCF file to region |
| [vcfplotaltdiscrepancy.r](./scripts/vcfplotaltdiscrepancy.r) | plot ALT discrepancies |
| [plot_roc.r](./scripts/plot_roc.r) | plot ROC |
| [vcfmultiallelic](./scripts/vcfmultiallelic) | remove anything that is not multiallelic |
| [vcfsnps](./scripts/vcfsnps) | show SNPs |
| [vcfvarstats](./scripts/vcfvarstats) | use fastahack to get stats |
| [vcfregionreduce_uncompressed](./scripts/vcfregionreduce_uncompressed) | reduce, gzip and tabix |
| [plotWCfst.R](./scripts/plotWCfst.R) | plot wcFst |
| [vcf2bed.py](./scripts/vcf2bed.py) | transform VCF to BED file |
| [vcfjoincalls](./scripts/vcfjoincalls) | overlay files using QUAL and GT from a second VCF |
| [vcf2sqlite.py](./scripts/vcf2sqlite.py) | push VCF file into SQLite3 database using dbname |

# Development

## build from source

VCFLIB uses the cmake build system, after a recursive checkout
of the sources make the files in the ./build directory with:

```sh
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
mkdir -p build && cd build
cmake ..
cmake --build .
cmake --install .
```

and to run the tests

```sh
ctest --verbose
```

Executables are built into the `./build` directory in the repository.

Note, if you have an existing repo update submodules with

```sh
git submodule update --init --recursive --progress
cd build
cmake --build . --target clean
```

Build dependencies can be viewed in the Travis-CI and github-CI
scripts (see badges above), as well as [guix.scm](./guix.scm) used by
us to create the build environment (for instructions see the header of
guix.scm). Essentially:

- C++ compiler
- htslib
- tabixpp

For include files add

- libhts-dev
- libtabixpp-dev
- libtabixpp0

And for some of the VCF executables

- python
- perl

### Using a different htslib

Check out htslib in tabixpp (recursively) and

    cmake -DHTSLIB_LOCAL:STRING=./tabixpp/htslib/ ..
    cmake --build .

## link library

The standard build creates `build/vcflib.a`. Take a hint from the
[cmake](./CMakeLists.txt) file that builds all the vcflib tools.

## source code

See [vcfecho.cpp](./src/vcfecho.cpp) for basic usage.
[Variant.h](./src/Variant.h) and [Variant.cpp](./src/Variant.cpp)
describe methods available in the API.  vcflib is incorporated into
several projects, such as
[freebayes](https://github.com/freebayes/freebayes), which may provide
a point of reference for prospective developers.  Note vcflib contains
submodules (git repositories) comprising some dependencies. A full
Guix development environment we use is defined [here](./guix.scm).

# adding tests

vcflib uses different test systems. The most important one is the
[doctest](https://docs.python.org/3/library/doctest.html) because it
doubles as documentation. For an example see
[vcf2tsv.md](./test/pytest/vcf2tsv.md) which can be run from the
command line with

```sh
cd test
python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF pytest/vcf2tsv.md
```

# Contributing

To contribute code to vcflib send a github pull request. We may ask
you to add a working test case as described in 'adding tests'.

## LICENSE

This software is distributed under the free software [MIT LICENSE](./LICENSE).

## CREDIT

Citations are the bread and butter of Science.  If you are using this
software in your research and want to support our future work, please
cite the following publication:

Vcflib and tools for processing the VCF variant call format;
Erik Garrison, Zev N. Kronenberg, Eric T. Dawson, Brent S. Pedersen, Pjotr Prins;
doi: https://doi.org/10.1101/2021.05.21.445151

## Bibtex reference

```bibtex
@article {Garrison2021.05.21.445151,
	author = {Garrison, Erik and Kronenberg, Zev N. and Dawson, Eric T. and Pedersen, Brent S. and Prins, Pjotr},
	title = {Vcflib and tools for processing the VCF variant call format},
	elocation-id = {2021.05.21.445151},
	year = {2021},
	doi = {10.1101/2021.05.21.445151},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/05/23/2021.05.21.445151},
	eprint = {https://www.biorxiv.org/content/early/2021/05/23/2021.05.21.445151.full.pdf},
	journal = {bioRxiv}
}
```
---
name: Bug report 🐞
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''
---
**Only bug reports!**

The C++ version of VCFLIB is in *maintenance* mode. Use the github issue
tracker to report bugs *only*. For comments, questions and features,
please use the google group mailing list as stated on the
[README](https://github.com/vcflib/vcflib/blob/master/README.md)!

**Describe the bug**

A clear and concise description of what the bug is.

**To Reproduce**

Include all steps to reproduce the behavior and paste any complete
errors from the terminal.

**Expected behavior**

A clear and concise description of what you expected to happen.

**Screenshots**

If applicable, add screenshots to help explain your problem.

**Additional context**

Add any other context about the problem here.

Include a set of VCF files to reproduce the issue

+ bonus points if you try to minimize the test case yourself, as issues are often localized:
  - try to use sambamba or samtools slice to first extract the reference where the error occurs
  - if that succeeds (the error is still reproducible), continue to crop the file in binary-search fashion

**Finally**

Please check the README and docs carefully. Everyone working on vcflib is doing that for free. Please respect our time (too).
% VCFNULLDOTSLASHDOT(1) vcfnulldotslashdot | Convert VCF . to ./.
% Erik Garrison and other vcflib contributors

# NAME

vcfnulldotslashdot converts single dots to ./.

This is useful for some downstream analysis tools.

# SYNOPSIS

**vcf2tsv** < *file*

# DESCRIPTION

vcfnulldotslashdot converts single dots to ./.

This is useful for some downstream analysis tools.

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES


<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

vcf2tsv converts a VCF to VCF

```python

>>> cat("../scripts/vcfnulldotslashdot data/issue_307_vcfnulldotslashdot.vcf")
##fileformat=VCFv4.1
##fileDate=20210209
##source=freeBayes v0.9.21
##reference=ahy.fa
##phasing=none
##filter="TYPE = snp & QUAL > 30 & AF > 0.05 & AF < 0.95 genotypes filtered with: GQ > 20"
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	rab-field.AH13
contig71493_103113	50	.	G	A	79.0603	.	AC=5;AF=0.078125;LEN=1;TYPE=snp	GT	./.:0,0,0

```

## Source code

[vcfnulldotslashdot](../scripts/vcfnulldotslashdot)

# LICENSE

Copyright 2021 (C) Erik Garrison and vcflib contributors. MIT licensed.
% VCF2TSV(1) vcf2tsv | Convert VCF to TSV
% Erik Garrison and other vcflib contributors

# NAME

vcf2tsv - Converts stdin or given VCF file to tab-delimited format,
using null string to replace empty values in the table.

# SYNOPSIS

**vcf2tsv** \[-n null_string] \[-g] \[*file*]

# DESCRIPTION

**vcf2tsv** converts stdin or given VCF file to tab-delimited format,
using null string to replace empty values in the table.

Specifying *-g* will output one line per sample with genotype
information.  When there is more than one alt allele there will be
multiple rows, one for each allele and, the info will match the 'A'
index

## Options

-h, --help

: shows help message and exits.

-g

: Output one line per sample with genotype information.

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES


<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

```

>>> head("vcf2tsv -h",1)
usage: vcf2tsv [-n null_string] [-g] [vcf file]


```

vcf2tsv converts a CSV to a tabulated test file, e.g.

```python

>>> head("vcf2tsv ../samples/sample.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2      NS
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .
19      112     .       A       G       10      .       .       .       .       .       .       .       .       .
20      14370   rs6054257       G       A       29      PASS    .       .       0.5     .       .       14      .       3

```

Use the `-g` switch to show genotypes

```python

>>> head("vcf2tsv -g ../samples/sample.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2      NS      SAMPLE  DP      GQ      GT      HQ
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00001 .       .       0|0     10,10
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00002 .       .       0|0     10,10
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .       NA00003 .       .       0/1     3,3

```

## Source code

[vcf2tsv.cpp](../../src/vcf2tsv.cpp)

## Regression tests

The following commands run full regression tests:

>>> run_stdout("vcf2tsv ../samples/sample.vcf", ext="tsv")
output in <a href="../data/regression/vcf2tsv_4.tsv">vcf2tsv_4.tsv</a>

>>> run_stdout("vcf2tsv -g ../samples/sample.vcf", ext="tsv")
output in <a href="../data/regression/vcf2tsv_5.tsv">vcf2tsv_5.tsv</a>


# LICENSE

Copyright 2020 (C) Erik Garrison and vcflib contributors. MIT licensed.
% WCFST(1) wcFst (vcflib) | wcFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**wcFst**

# SYNOPSIS

**wcFst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

# DESCRIPTION

**wcFst** is Weir & Cockerham's Fst for two populations. Negative values are VALID, they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. Specifically **wcFst** uses equations 1,2,3,4.



# OPTIONS

```


Output : 3 columns :                 
     1. seqid                        
     2. position                     
     3. target allele frequency      
     4. background allele frequency  
     5. **wcFst**                        

required: t,target     -- argument: a zero based comma separated list of target individuals corrisponding to VCF columns        
required: b,background -- argument: a zero based comma separated list of background individuals corrisponding to VCF columns    
required: f,file       -- argument: proper formatted VCF                                                                        
required, y,type       -- argument: genotype likelihood format; genotype : GT,GL,PL,GP                                             
optional: r,region     -- argument: a tabix compliant genomic range: seqid or seqid:start-end                                   
optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero 

Type: statistics


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[wcFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/wcFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFANNOTATE(1) vcfannotate (vcflib) | vcfannotate (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfannotate**

# SYNOPSIS

**vcfannotate** [options] [<vcf file>]

# DESCRIPTION

Intersect the records in the VCF file with targets provided in a BED file. Intersections are done on the reference sequences in the VCF file. If no VCF filename is specified on the command line (last argument) the VCF read from stdin.



# OPTIONS

```


options:
    -b, --bed   use annotations provided by this BED file
    -k, --key   use this INFO field key for the annotations
    -d, --default  use this INFO field key for records without annotations

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfannotate.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotate.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCHECK(1) vcfcheck (vcflib) | vcfcheck (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfcheck**

# SYNOPSIS

**vcfcheck** [options] <vcf file>

# DESCRIPTION

Validate integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file.



# OPTIONS

```

options:
    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
    -x, --exclude-failures If a record fails, don't print it.  Otherwise do.
    -k, --keep-failures    Print if the record fails, otherwise not.
    -h, --help       Print this message.
    -v, --version    Print version.


Type: metrics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcheck.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcheck.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFHETCOUNT(1) vcfhetcount (vcflib) | vcfhetcount (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfhetcount**

# SYNOPSIS

**vcfhetcount** <vcf file>

# DESCRIPTION

Calculate the heterozygosity rate: count the number of alternate alleles in heterozygous genotypes in all records in the vcf file



# OPTIONS

```

outputs a count for each individual in the file

Type: metrics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfhetcount.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfhetcount.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% SMOOTHER(1) smoother (vcflib) | smoother (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**smoother**

# SYNOPSIS

**smoother** --format pFst --file GPA.output.txt

# DESCRIPTION

smoothes is a method for window smoothing many of the GPAT++ formats.



# OPTIONS

```


      **smoother** averages a set of scores over a sliding genomic window.            
      **smoother** slides over genomic positions not the SNP indices. In other words  
      the number of scores within a window will not be constant. The last         
      window for each seqid can be smaller than the defined window size.          
      **smoother** automatically analyses different seqids separately.                
Output : 4 columns :     
     1. seqid            
     2. window start     
     2. window end       
     3. averaged score   

required: f,file     -- argument: a file created by GPAT++                           
required: o,format   -- argument: format of input file, case sensitive               
                              available format options:                                    
                                wcFst, pFst, bFst, iHS, xpEHH, abba-baba, col3             
optional: w,window   -- argument: size of genomic window in base pairs (default 5000)
optional: s,step     -- argument: window step size in base pairs (default 1000)      
optional: t,truncate -- flag    : end last window at last position (zero based)      

Type: transformation


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[smoother.cpp](https://github.com/vcflib/vcflib/blob/master/src/smoother.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFINFOSUMMARIZE(1) vcfinfosummarize (vcflib) | vcfinfosummarize (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfinfosummarize**

# SYNOPSIS

**vcfinfosummarize** [options] <vcf file>

# DESCRIPTION

Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO.



# OPTIONS

```

options:
    -f, --field         Summarize this field in the INFO column
    -i, --info          Store the computed statistic in this info field
    -a, --average       Take the mean for field (default)
    -m, --median        Use the median
    -n, --min           Use the min
    -x, --max           Use the max
    -h, --help          Print this message
    -v, --version       Print version


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfinfosummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfinfosummarize.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% PVST(1) pVst (vcflib) | pVst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**pVst**

# SYNOPSIS

**pVst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --type CN

# DESCRIPTION

**pVst** calculates vst, a measure of CNV stratification.



# OPTIONS

```




The statistic Vst is used to test the difference in copy numbers at
each SV between two groups: Vst = (Vt-Vs)/Vt, where Vt is the overall
variance of copy number and Vs the average variance within
populations.

Output : 4 columns :     
     1. seqid            
     2. position         
     3. end              
     3. vst              
     4. probability      

required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       
required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   
required: f,file       -- argument: a properly formatted VCF.                                                                  
required: y,type       -- argument: the genotype field with the copy number: e.g. CN|CNF                           
optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 
optional: x,cpu        -- argument: number of CPUs [1] 
optional: n,per        -- argument: number of permutations [1000] 

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[pVst.cpp](https://github.com/vcflib/vcflib/blob/master/src/pVst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% IHS(1) iHS (vcflib) | iHS (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**iHS**

# SYNOPSIS

**iHS** --target 0,1,2,3,4,5,6,7 --file my.phased.vcf \ --region chr1:1-1000 > STDOUT 2> STDERR

# DESCRIPTION

**iHS** calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014).



# OPTIONS

```


Our code is highly concordant with both implementations mentioned. However, we do not set an upper limit to the allele frequency.  **iHS** can be run without a genetic map, in which case the change in EHH is integrated over a constant.  Human genetic maps for GRCh36 and GRCh37 (hg18 & hg19) can be found at: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ . **iHS** by default interpolates SNV positions to genetic position (you don't need a genetic position for every VCF entry in the map file).

**iHS** analyses requires normalization by allele frequency.  It is important that **iHS** is calculated over large regions so that the normalization does not down weight real signals.  For genome-wide runs it is recommended to run slightly overlapping windows and throwing out values that fail integration (columns 7 & 8 in the output) and then removing duplicates by using the 'sort' and 'uniq' linux commands.  Normalization of the output is as simple as running 'normalize-**iHS**'.



     **iHS** calculates the integrated ratio of haplotype decay between the reference and non-reference allele.
Output : 4 columns :
     1. seqid
     2. position
     3. target allele frequency
     4. integrated EHH (alternative)
     5. integrated EHH (reference)
     6. **iHS** ln(iEHHalt/iEHHref)
     7. != 0 integration failure
     8. != 0 integration failure

Params:
       required: t,target  <STRING>  A zero base comma separated list of target
                                     individuals corresponding to VCF columns
       required: r,region  <STRING>  A tabix compliant genomic range
                                     format: "seqid:start-end" or "seqid"
       required: f,file    <STRING>  Proper formatted and phased VCF.
       required: y,type    <STRING>  Genotype likelihood format: GT,PL,GL,GP
       optional: a,af      <DOUBLE>  Alternative alleles with frquences less
                                     than [0.05] are skipped.
       optional: x,threads <INT>     Number of CPUS [1].
       recommended: g,gen <STRING>   A PLINK formatted map file.



Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/iHS.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% POPSTATS(1) popStats (vcflib) | popStats (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**popStats**

# SYNOPSIS

popStat --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf

# DESCRIPTION

General population genetic statistics for each SNP



# OPTIONS

```



  Calculates basic population statistics at bi-allelic sites. The allele frequency is the number of non-reference alleles divided by the total number of alleles.  The expected hetrozygosity is 2*p*q, where p is the non-reference allele frequency and q is 1-p.  The observed heterozgosity is the fraction of 0/1 genotypes out of all genotypes.  The inbreeding coefficient, Fis, is the relative heterozygosity of each individual vs. compared to the target group. 

Output : 9 columns :                 
     1. seqid                        
     2. position                     
     3. target allele frequency      
     4. expected heterozygosity      
     5. observed heterozygosity      
     6. number of hets               
     7. number of homozygous ref     
     8. number of homozygous alt     
     9. target Fis                   
required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns        
required: f,file       -- proper formatted VCF                                                                        
required, y,type       -- genotype likelihood format; genotype : GL,PL,GP                                             
optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1                                              

Type: statistics


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[popStats.cpp](https://github.com/vcflib/vcflib/blob/master/src/popStats.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSAMPLENAMES(1) vcfsamplenames (vcflib) | vcfsamplenames (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfsamplenames**

# SYNOPSIS

**vcfsamplenames**

# DESCRIPTION

List sample names



# OPTIONS

```


Type: transformation

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfsamplenames.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsamplenames.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFALTCOUNT(1) vcfaltcount (vcflib) | vcfaltcount (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfaltcount**

# SYNOPSIS

**vcfaltcount** <vcf file>

# DESCRIPTION

count the number of alternate alleles in all records in the vcf file



# OPTIONS

```


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfaltcount.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfaltcount.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCLASSIFY(1) vcfclassify (vcflib) | vcfclassify (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfclassify**

# SYNOPSIS

**vcfclassify** <vcf file>

# DESCRIPTION

Creates a new VCF where each variant is tagged by allele class: snp, ts/tv, indel, mnp



# OPTIONS

```


Type: transformation

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfclassify.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfclassify.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFPRIMERS(1) vcfprimers (vcflib) | vcfprimers (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfprimers**

# SYNOPSIS

**vcfprimers** [options] <vcf file>

# DESCRIPTION

For each VCF record, extract the flanking sequences, and write them to stdout as FASTA records suitable for alignment.



# OPTIONS

```

options:
    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
    -l, --primer-length    The length of the primer sequences on each side of the variant

This tool is intended for use in designing validation
experiments.  Primers extracted which would flank all of the alleles at multi-allelic
sites.  The name of the FASTA "reads" indicates the VCF record which they apply to.
The form is >CHROM_POS_LEFT for the 3' primer and >CHROM_POS_RIGHT for the 5' primer,
for example:

>20_233255_LEFT
CCATTGTATATATAGACCATAATTTCTTTATCCAATCATCTGTTGATGGA
>20_233255_RIGHT
ACTCAGTTGATTCCATACCTTTGCCATCATGAATCATGTTGTAATAAACA


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfprimers.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfprimers.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% vcflib(1) vcflib | vcflib (index)
% Erik Garrison and vcflib contributors

# NAME

**vcflib** index

# DESCRIPTION

vcflib contains tools and libraries for dealing with the Variant Call
Format (VCF) which is a flat-file, tab-delimited textual format
intended to describe reference-indexed variations between
individuals.

VCF provides a common interchange format for the description of
variation in individuals and populations of samples, and has become
the defacto standard reporting format for a wide array of genomic
variant detectors.

vcflib provides methods to manipulate and interpret sequence variation
as it can be described by VCF. It is both:

* an API for parsing and operating on records of genomic variation as it can be described by the VCF format,
* and a collection of command-line utilities for executing complex manipulations on VCF files.

The API itself provides a quick and extremely permissive method to
read and write VCF files. Extensions and applications of the library
provided in the included utilities (*.cpp) comprise the vast bulk of
the library's utility for most users.

<!--
  Created with ./scripts/bin2md.rb --index
-->


## filter

| filter command | description |
| :-------------- | :---------- |

## metrics

| metrics command | description |
| :-------------- | :---------- |

## phenotype

| phenotype command | description |
| :-------------- | :---------- |

## genotype

| genotype command | description |
| :-------------- | :---------- |

## transformation

| transformation command | description |
| :-------------- | :---------- |

## statistics

| statistics command | description |
| :-------------- | :---------- |

# SOURCE CODE

See the source code repository at https://github.com/vcflib/vcflib

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

% VCFRANDOMSAMPLE(1) vcfrandomsample (vcflib) | vcfrandomsample (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfrandomsample**

# SYNOPSIS

**vcfrandomsample** [options] [<vcf file>]

# DESCRIPTION

Randomly sample sites from an input VCF file, which may be provided as stdin. Scale the sampling probability by the field specified in KEY. This may be used to provide uniform sampling across allele frequencies, for instance.



# OPTIONS

```

options:
    -r, --rate RATE          base sampling probability per locus
    -s, --scale-by KEY       scale sampling likelihood by this Float info field
    -p, --random-seed N      use this random seed (by default read from /dev/random)
    -q, --pseudorandom-seed  use a pseudorandom seed (by default read from /dev/random)


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfrandomsample.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfrandomsample.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFREMOVESAMPLES(1) vcfremovesamples (vcflib) | vcfremovesamples (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfremovesamples**

# SYNOPSIS

**vcfremovesamples** <vcf file> [SAMPLE1] [SAMPLE2] ...

# DESCRIPTION

outputs each record in the vcf file, removing samples listed on the command line





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfremovesamples.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfremovesamples.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGLXGT(1) vcfglxgt (vcflib) | vcfglxgt (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfglxgt**

# SYNOPSIS

**vcfglxgt** [options] <vcf file>

# DESCRIPTION

Set genotypes using the maximum genotype likelihood for each sample.



# OPTIONS

```

options:
    -n, --fix-null-genotypes   only apply to null and partly-null genotypes



Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfglxgt.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglxgt.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCOMMONSAMPLES(1) vcfcommonsamples (vcflib) | vcfcommonsamples (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcommonsamples**

# SYNOPSIS

**vcfcommonsamples** <vcf file> <vcf file>

# DESCRIPTION

Generates each record in the first file, removing samples not present in the second



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcommonsamples.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcommonsamples.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% PLOTHAPS(1) plotHaps (vcflib) | plotHaps (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**plotHaps**

# SYNOPSIS



# DESCRIPTION

**plotHaps** provides the formatted output that can be used with 'bin/plotHaplotypes.R'.



# OPTIONS

```


Output : haplotype matrix and positions

**plotHaps**  --target 0,1,2,3,4,5,6,7  --file my.phased.vcf.gz                                                           

required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF column s        
required: r,region     -- argument: a tabix compliant genomic range : "seqid:start-end" or "seqid"                          
required: f,file       -- argument: proper formatted phased VCF file                                                            
required: y,type       -- argument: genotype likelihood format: PL,GP,GP                                                        

Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[plotHaps.cpp](https://github.com/vcflib/vcflib/blob/master/src/plotHaps.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFFLATTEN(1) vcfflatten (vcflib) | vcfflatten (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfflatten**

# SYNOPSIS

**vcfflatten** [options][file] -h --help display this help message and exit. -i --ignore-errors do not flatten locus if 'AF' is not specified.

# DESCRIPTION

Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfflatten.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfflatten.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFBREAKMULTI(1) vcfbreakmulti (vcflib) | vcfbreakmulti (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfbreakmulti**

# SYNOPSIS

**vcfbreakmulti** [options] [file]

# DESCRIPTION

If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfbreakmulti.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfbreakmulti.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% PERMUTESMOOTH(1) permuteSmooth (vcflib) | permuteSmooth (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**permuteSmooth**

# SYNOPSIS

**permuteSmooth** -s wcFst.smooth.txt -f wcFst.txt -n 5 -s 1

# DESCRIPTION

**permuteSmooth** is a method for adding empirical p-values smoothed wcFst scores.



# OPTIONS

```


Required:
      file:     f   -- argument: original wcFst data     
      smoothed: s   -- argument: smoothed wcFst data     
      format:   y   -- argument: [swcFst, segwcFst]      
Optional:
      number:   n   -- argument: the number of permutations to run for each value [1000]
      success:  u   -- argument: stop permutations after 's' successes [1]
      success:  x   -- argument: number of threads [1]

OUTPUT: **permuteSmooth** will append three additional columns:
        1. The number of successes                            
        2. The number of trials                               
        3. The empirical p-value                              


Type: statistics


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[permuteSmooth.cpp](https://github.com/vcflib/vcflib/blob/master/src/permuteSmooth.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFLD(1) vcfld (vcflib) | vcfld (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfld**

# SYNOPSIS

**vcfld** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf -e -d -r

# DESCRIPTION

Compute LD



# OPTIONS

```


required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        
required: b,background -- argument: a zero base comma separated list of background individuals corresponding to VCF columns    
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: w,window     -- argument: window size to average LD; default is 1000                                                 
optional: e,external   -- switch: population to calculate LD expectation; default is target                                    
optional: d,derived    -- switch: which haplotype to count "00" vs "11"; default "00",                                   


Type: transformation


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfld.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfld.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFPARSEALTS(1) vcfparsealts (vcflib) | vcfparsealts (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfparsealts**

# SYNOPSIS

**vcfparsealts** <vcf file>

# DESCRIPTION

Alternate allele parsing method. This method uses pairwise alignment of REF and ALTs to determine component allelic primitives for each alternate allele.





# EXAMPLES

```

Example:

**vcfparsealts** samples/sample.vcf
##fileformat=VCFv4.0
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
19      111     .       A       C       9.6     .       .       GT:HQ   0|0:10,10       0|0:10,10       0/1:3,3
 ( A :: 111 A -> A;  )  ( C :: 111 A -> C;  )
19      112     .       A       G       10      .       .       GT:HQ   0|0:10,10       0|0:10,10       0/1:3,3
 ( A :: 112 A -> A;  )  ( G :: 112 A -> G;  )
20      14370   rs6054257       G       A       29      PASS    AF=0.5;DP=14;NS=3;DB;H2 GT:GQ:DP:HQ     0|0:48:1:51,51     1|0:48:8:51,51  1/1:43:5:.,.
 ( A :: 14370 G -> A;  )  ( G :: 14370 G -> G;  )
20      17330   .       T       A       3       q10     AF=0.017;DP=11;NS=3     GT:GQ:DP:HQ     0|0:49:3:58,50     0|1:3:5:65,3    0/0:41:3:.,.
 ( A :: 17330 T -> A;  )  ( T :: 17330 T -> T;  )
20      1110696 rs6040355       A       G,T     67      PASS    AA=T;AF=0.333,0.667;DP=10;NS=2;DB       GT:GQ:DP:HQ        1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4:.,.
 ( A :: 1110696 A -> A;  )  ( G :: 1110696 A -> G;  )  ( T :: 1110696 A -> T;  )
20      1230237 .       T       .       47      PASS    AA=T;DP=13;NS=3 GT:GQ:DP:HQ     0|0:54:.:56,60  0|0:48:4:51,51     0/0:61:2:.,.
 ( . :: 1230237 T -> .;  )  ( T :: 1230237 T -> T;  )
20      1234567 microsat1       G       GA,GAC  50      PASS    AA=G;AC=3,1;AN=6;DP=9;NS=3      GT:GQ:DP  0/1:.:4  0/2:17:2        1/1:40:3
 ( G :: 1234567 G -> G;  )  ( GA :: 1234567 G -> G; 1234568  -> A;  )  ( GAC :: 1234567 G -> G; 1234568  -> AC;  )
20      1235237 .       T       .       0       .       .       GT      0/0     0|0     ./.
 ( . :: 1235237 T -> .;  )  ( T :: 1235237 T -> T;  )
X       10      rsTest  AC      A,ATG   10      PASS    .       GT      0       0/1     0|2
 ( A :: 10 A -> A; 11 C -> ;  )  ( AC :: 10 AC -> AC;  )  ( ATG :: 10 A -> A; 11  -> T; 11 C -> G;  )


Type: statistics

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfparsealts.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfparsealts.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFDISTANCE(1) vcfdistance (vcflib) | vcfdistance (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfdistance**

# SYNOPSIS

**vcfdistance** [customtagname] < [vcf file]

# DESCRIPTION

Adds a tag to each variant record which indicates the distance to the nearest variant. (defaults to BasesToClosestVariant if no custom tag name is given.



# OPTIONS

```


Type: metrics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfdistance.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfdistance.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFROC(1) vcfroc (vcflib) | vcfroc (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfroc**

# SYNOPSIS

**vcfroc** [options] [<vcf file>]

# DESCRIPTION

Generates a pseudo-ROC curve using sensitivity and specificity estimated against a putative truth set. Thresholding is provided by successive QUAL cutoffs.



# OPTIONS

```

options:
    -t, --truth-vcf FILE      use this VCF as ground truth for ROC generation
    -w, --window-size N       compare records up to this many bp away (default 30)
    -c, --complex             directly compare complex alleles, don't parse into primitives
    -r, --reference FILE      FASTA reference file


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfroc.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfroc.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFKEEPSAMPLES(1) vcfkeepsamples (vcflib) | vcfkeepsamples (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfkeepsamples**

# SYNOPSIS

**vcfkeepsamples** <vcf file> [SAMPLE1] [SAMPLE2] ...

# DESCRIPTION

outputs each record in the vcf file, removing samples not listed on the command line





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfkeepsamples.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfkeepsamples.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% HAPLRT(1) hapLrt (vcflib) | hapLrt (VCF genotype)
% Erik Garrison and vcflib contributors

# NAME

**hapLrt**

# SYNOPSIS

hapLRT --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --type GP --file my.vcf

# DESCRIPTION

HapLRT is a likelihood ratio test for haplotype lengths. The lengths are modeled with an exponential distribution. The sign denotes if the target has longer haplotypes (1) or the background (-1).



# OPTIONS

```


Output : 4 columns :                             
     1. seqid                                    
     2. position                                 
     3. mean target haplotype length             
     4. mean background haplotype length         
     5. p-value from LRT                         
     6. sign                                     

required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        
required: b,background -- argument: a zero base comma separated list of background individuals corresponding to VCF columns    
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: r,region     -- argument: a genomic range to calculate **hapLrt** on in the format : "seqid:start-end" or "seqid" 


Type: genotype

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[hapLrt.cpp](https://github.com/vcflib/vcflib/blob/master/src/hapLrt.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENO2ALLELES(1) vcfgeno2alleles (vcflib) | vcfgeno2alleles (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfgeno2alleles**

# SYNOPSIS

**vcfgeno2alleles** <[vcf file]

# DESCRIPTION

modifies the genotypes field to provide the literal alleles rather than indexes



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgeno2alleles.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgeno2alleles.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFUNIQALLELES(1) vcfuniqalleles (vcflib) | vcfuniqalleles (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcfuniqalleles**

# SYNOPSIS

**vcfuniqalleles** <vcf file>

# DESCRIPTION

List unique alleles For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files.



# OPTIONS

```


Type: filter

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfuniqalleles.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfuniqalleles.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFUNIQ(1) vcfuniq (vcflib) | vcfuniq (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcfuniq**

# SYNOPSIS

**vcfuniq** <vcf file>

# DESCRIPTION

List unique genotypes. Like GNU uniq, but for VCF records. Remove records which have the same position, ref, and alt as the previous record.



# OPTIONS

```


Type: filter

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfuniq.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfuniq.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFINTERSECT(1) vcfintersect (vcflib) | vcfintersect (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfintersect**

# SYNOPSIS

**vcfintersect** [options] [<vcf file>]

# DESCRIPTION

VCF set analysis



# OPTIONS

```


options:
    -b, --bed FILE            use intervals provided by this BED file
    -R, --region REGION       use 1-based tabix-style region (e.g. chrZ:10-20), multiples allowed
    -S, --start-only          don't use the reference length information in the record to determine
                              overlap status, just use the start posiion
    -v, --invert              invert the selection, printing only records which would
                                not have been printed out
    -i, --intersect-vcf FILE  use this VCF for set intersection generation
    -u, --union-vcf FILE      use this VCF for set union generation
    -w, --window-size N       compare records up to this many bp away (default 30)
    -r, --reference FILE      FASTA reference file, required with -i and -u
    -l, --loci                output whole loci when one alternate allele matches
    -m, --ref-match           intersect on the basis of record REF string
    -t, --tag TAG             attach TAG to each record's info field if it would intersect
    -V, --tag-value VAL       use this value to indicate that the allele is passing
                              '.' will be used otherwise.  default: 'PASS'
    -M, --merge-from FROM-TAG
    -T, --merge-to   TO-TAG   merge from FROM-TAG used in the -i file, setting TO-TAG
                              in the current file.

For bed-vcf intersection, alleles which fall into the targets are retained.

Haplotype aware intersection, union and complement. Use for intersection and union of VCF files: unify on equivalent alleles within window-size bp
as determined by haplotype comparison alleles.

type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfintersect.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfintersect.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENOTYPECOMPARE(1) vcfgenotypecompare (vcflib) | vcfgenotypecompare (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenotypecompare**

# SYNOPSIS

**vcfgenotypecompare** <other-genotype-tag> <vcf file>

# DESCRIPTION

adds statistics to the INFO field of the vcf file describing the amount of discrepancy between the genotypes (GT) in the vcf file and the genotypes reported in the <other-genotype-tag>. use this after vcfannotategenotypes to get correspondence statistics for two vcfs.



# OPTIONS

```


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgenotypecompare.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenotypecompare.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFADDINFO(1) vcfaddinfo (vcflib) | vcfaddinfo (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfaddinfo**

# SYNOPSIS

**vcfaddinfo** <vcf file> <vcf file>

# DESCRIPTION

Adds info fields from the second file which are not present in the first vcf file.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfaddinfo.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfaddinfo.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSTATS(1) vcfstats (vcflib) | vcfstats (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfstats**

# SYNOPSIS

**vcfstats** [options] <vcf file>

# DESCRIPTION

Prints statistics about variants in the input VCF file.



# OPTIONS

```


    -r, --region          specify a region on which to target the stats, requires a BGZF
                          compressed file which has been indexed with tabix.  any number of
                          regions may be specified.
    -a, --add-info        add the statistics intermediate information to the VCF file,
                          writing out VCF records instead of summary statistics
    -t, --add-type        only add the type= field to the info (faster than -a)
    -l, --no-length-frequency    don't out the indel and mnp length-frequency spectra
    -m, --match-score N          match score for SW algorithm
    -x, --mismatch-score N       mismatch score for SW algorithm
    -o, --gap-open-penalty N     gap open penalty for SW algorithm
    -e, --gap-extend-penalty N   gap extension penalty for SW algorithm


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfstats.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfstats.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCOMBINE(1) vcfcombine (vcflib) | vcfcombine (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcombine**

# SYNOPSIS

**vcfcombine** [vcf file] [vcf file] ...

# DESCRIPTION

Combine VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted.



# OPTIONS

```


options:
    -h --help           This text.
    -v --version        Print version.
    -r --region REGION  A region specifier of the form chrN:x-y to bound the merge

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcombine.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcombine.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFEVENREGIONS(1) vcfevenregions (vcflib) | vcfevenregions (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfevenregions**

# SYNOPSIS

**vcfevenregions** [options] <vcf file>

# DESCRIPTION

Generates a list of regions, e.g. chr20:10..30 using the variant density information provided in the VCF file to ensure that the regions have even numbers of variants. This can be use to reduce the variance in runtime when dividing variant detection or genotyping by genomic coordinates.



# OPTIONS

```

options:
    -f, --fasta-reference REF    FASTA reference file to use to obtain primer sequences.
    -n, --number-of-regions N    The number of desired regions.
    -p, --number-of-positions N  The number of positions per region.
    -o, --offset N               Add an offset to region positioning, to avoid boundary
                                 related artifacts in downstream processing.
    -l, --overlap N              The number of sites to overlap between regions.  Default 0.
    -s, --separator SEQ          Specify string to use to separate region output.  Default '-'

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfevenregions.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfevenregions.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFNUMALT(1) vcfnumalt (vcflib) | vcfnumalt (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfnumalt**

# SYNOPSIS

**vcfnumalt** <vcf file>

# DESCRIPTION

outputs a VCF stream where NUMALT has been generated for each record using sample genotypes





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfnumalt.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfnumalt.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% SEQUENCEDIVERSITY(1) sequenceDiversity (vcflib) | sequenceDiversity (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**sequenceDiversity**

# SYNOPSIS

**sequenceDiversity** --target 0,1,2,3,4,5,6,7 --file my.vcf

# DESCRIPTION

The **sequenceDiversity** program calculates two popular metrics of haplotype diversity: pi and extended haplotype homozygoisty (eHH). Pi is calculated using the Nei and Li 1979 formulation. eHH a convenient way to think about haplotype diversity. When eHH = 0 all haplotypes in the window are unique and when eHH = 1 all haplotypes in the window are identical.



# OPTIONS

```


Output : 5 columns:
         1.  seqid
         2.  start of window
         3.  end of window  
         4.  pi             
         5.  eHH            


required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: a,af         -- sites less than af  are filtered out; default is 0                                          
optional: r,region     -- argument: a tabix compliant region : "seqid:0-100" or "seqid"                                    
optional: w,window     -- argument: the number of SNPs per window; default is 20                                               

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[sequenceDiversity.cpp](https://github.com/vcflib/vcflib/blob/master/src/sequenceDiversity.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFALLELICPRIMITIVES(1) vcfallelicprimitives (vcflib) | vcfallelicprimitives (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfallelicprimitives**

# SYNOPSIS

**vcfallelicprimitives** [options] [file]

# DESCRIPTION

If multiple allelic primitives (gaps or mismatches) are specified in a single VCF record, split the record into multiple lines, but drop all INFO fields. Does not handle genotypes (yet). MNPs are split into multiple SNPs unless the -m flag is provided. Records generated by splits have th



# OPTIONS

```

options:
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Tag records which are split apart of a complex allele with this flag.
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: 200).
    -k, --keep-info         Maintain site and allele-level annotations when decomposing.
                            Note that in many cases, such as multisample VCFs, these won't
                            be valid post-decomposition.  For biallelic loci in single-sample
                            VCFs, they should be usable with caution.
    -g, --keep-geno         Maintain genotype-level annotations when decomposing.  Similar
                            caution should be used for this as for --keep-info.

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfallelicprimitives.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfallelicprimitives.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGLBOUND(1) vcfglbound (vcflib) | vcfglbound (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfglbound**

# SYNOPSIS

**vcfglbound** [options] <vcf file>

# DESCRIPTION

Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max.



# OPTIONS

```


Then cap (bound) at N (e.g. -10).options:
    -b, --bound N          Bound GLs to this limit.
    -x, --exclude-broken   If GLs are > 0, remove site.


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfglbound.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglbound.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFLENGTH(1) vcflength (vcflib) | vcflength (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcflength**

# SYNOPSIS

**vcflength**

# DESCRIPTION

Add length info field





# EXAMPLES

```

Example:

**vcflength** samples/sample.vcf
##fileformat=VCFv4.0
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
19      111     .       A       C       9.6     .       length=0;length.alt=1;length.ref=1      GT:HQ   0|0:10,10  0|0:10,10       0/1:3,3
19      112     .       A       G       10      .       length=0;length.alt=1;length.ref=1      GT:HQ   0|0:10,10  0|0:10,10       0/1:3,3
20      14370   rs6054257       G       A       29      PASS    AF=0.5;DP=14;NS=3;length=0;length.alt=1;length.ref=1;DB;H2 GT:GQ:DP:HQ     0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
20      17330   .       T       A       3       q10     AF=0.017;DP=11;NS=3;length=0;length.alt=1;length.ref=1     GT:GQ:DP:HQ     0|0:49:3:58,50  0|1:3:5:65,3    0/0:41:3:.,.
20      1110696 rs6040355       A       G,T     67      PASS    AA=T;AF=0.333,0.667;DP=10;NS=2;length=0,0;length.alt=1,1;length.ref=1;DB   GT:GQ:DP:HQ     1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4:.,.
20      1230237 .       T       .       47      PASS    AA=T;DP=13;NS=3;length=0;length.alt=1;length.ref=1GT:GQ:DP:HQ      0|0:54:.:56,60  0|0:48:4:51,51  0/0:61:2:.,.
20      1234567 microsat1       G       GA,GAC  50      PASS    AA=G;AC=3,1;AN=6;DP=9;NS=3;length=1,2;length.alt=2,3;length.ref=1  GT:GQ:DP        0/1:.:4 0/2:17:2        1/1:40:3
20      1235237 .       T       .       0       .       length=0;length.alt=1;length.ref=1      GT      0/00|0     ./.
X       10      rsTest  AC      A,ATG   10      PASS    length=-1,1;length.alt=1,3;length.ref=2 GT      0 0/1      0|2

Type: transformation

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcflength.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcflength.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFECHO(1) vcfecho (vcflib) | vcfecho (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfecho**

# SYNOPSIS

**vcfecho** <vcf file>

# DESCRIPTION

Echo VCF to stdout (simple demo)



# OPTIONS

```


Type: transformation

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfecho.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfecho.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCF2DAG(1) vcf2dag (vcflib) | vcf2dag (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcf2dag**

# SYNOPSIS

**vcf2dag** [options] [<vcf file>]

# DESCRIPTION

Modify VCF to be able to build a directed acyclic graph (DAG)



# OPTIONS

```

options:
    -r, --reference FILE         FASTA reference file.

Modify the VCF file so that homozygous regions are included as REF/. calls.
For each ref and alt allele, assign an index.  These steps are sufficient to
enable use of the VCF as a DAG (specifically a partially-ordered graph).

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcf2dag.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2dag.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFREMOVEABERRANTGENOTYPES(1) vcfremoveaberrantgenotypes (vcflib) | vcfremoveaberrantgenotypes (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfremoveaberrantgenotypes**

# SYNOPSIS

**vcfremoveaberrantgenotypes** <vcf file>

# DESCRIPTION

strips samples which are homozygous but have observations implying heterozygosity. Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO).



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfremoveaberrantgenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfremoveaberrantgenotypes.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFFILTER(1) vcffilter (vcflib) | vcffilter (VCF filter)
% Erik Garrison and vcflib contributors

# NAME

**vcffilter**

# SYNOPSIS

**vcffilter** [options] <vcf file>

# DESCRIPTION

VCF filter the specified vcf file using the set of filters



# OPTIONS

```


options:
    -f, --info-filter     specifies a filter to apply to the info fields of records,
                          removes alleles which do not pass the filter
    -g, --genotype-filter specifies a filter to apply to the genotype fields of records
    -k, --keep-info       used in conjunction with '-g', keeps variant info, but removes genotype
    -s, --filter-sites    filter entire records, not just alleles
    -t, --tag-pass        tag vcf records as positively filtered with this tag, print all records
    -F, --tag-fail        tag vcf records as negatively filtered with this tag, print all records
    -A, --append-filter   append the existing filter tag, don't just replace it
    -a, --allele-tag      apply -t on a per-allele basis.  adds or sets the corresponding INFO field tag
    -v, --invert          inverts the filter, e.g. grep -v
    -o, --or              use logical OR instead of AND to combine filters
    -r, --region          specify a region on which to target the filtering, requires a BGZF
                          compressed file which has been indexed with tabix.  any number of
                          regions may be specified.

Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields

Operators can be any of: =, !, <, >, |, &

Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"

For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.

type: filter

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcffilter.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcffilter.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% NORMALIZE-IHS(1) normalize-iHS (vcflib) | normalize-iHS (VCF genotype)
% Erik Garrison and vcflib contributors

# NAME

**normalize-iHS**

# SYNOPSIS

normalizeHS -s 0.01 -f input.txt

# DESCRIPTION

normalizes iHS or XP-EHH scores.



# OPTIONS

```




A cross-population extended haplotype homozygosity (XP-EHH) score is
directional: a positive score suggests selection is likely to have
happened in population A, whereas a negative score suggests the same
about population B. See for example
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687721/


Output : **normalize-iHS** adds one additional column to input (normalized score).
required: -f            -- Output from iHS or XPEHH 
optional: -s            -- Max AF diff for window [0.01]

Type: genotype



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[normalize-iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/normalize-iHS.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENOSAMPLENAMES(1) vcfgenosamplenames (vcflib) | vcfgenosamplenames (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenosamplenames**

# SYNOPSIS

**vcfgenosamplenames**

# DESCRIPTION

Get samplenames





# EXAMPLES

```

Example:

vcfsamplenames samples/sample.vcf

NA00001
NA00002
NA00003


Type: transformation

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgenosamplenames.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenosamplenames.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFENTROPY(1) vcfentropy (vcflib) | vcfentropy (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfentropy**

# SYNOPSIS

**vcfentropy** [options] <vcf file>

# DESCRIPTION

Annotate VCF records with the Shannon entropy of flanking sequence. Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight, EntropyCenter, which are the entropies of the sequence of the given window size to the left, right, and center of the record. Also adds EntropyRef and EntropyAlt for each alt.



# OPTIONS

```

options:
    -f, --fasta-reference  FASTA reference file to use to obtain flanking sequences
    -w, --window-size      Size of the window over which to calculate entropy



Type: metrics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfentropy.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfentropy.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSITESUMMARIZE(1) vcfsitesummarize (vcflib) | vcfsitesummarize (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfsitesummarize**

# SYNOPSIS

**vcfsitesummarize** <vcf file>

# DESCRIPTION

Summarize by site





# EXAMPLES

```

Example:

**vcfsitesummarize** samples/sample.vcf

CHROM   POS     ID      REF     QUAL    FILTER  AA      AC      AF      AN      DP      NS      DB      H2
19      111     .       A       9.6     .                                                       0       0
19      112     .       A       10      .                                                       0       0
20      14370   rs6054257       G       29      PASS                    0.5             14      3       1 1
20      17330   .       T       3       q10                     0.017           11      3       0       0
20      1110696 rs6040355       A       67      PASS    T                               10      2       1 0
20      1230237 .       T       47      PASS    T                               13      3       0       0
20      1234567 microsat1       G       50      PASS    G                       6       9       3       0 0
20      1235237 .       T       0       .                                                       0       0
X       10      rsTest  AC      10      PASS


Type: statistics

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfsitesummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsitesummarize.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFAFPATH(1) vcfafpath (vcflib) | vcfafpath (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfafpath**

# SYNOPSIS

**vcfafpath** <vcf file>

# DESCRIPTION

Display genotype paths





# EXAMPLES

```

Example:

    **vcfafpath** samples/scaffold612.vcf

```

T -> A
A -> G
T -> C
C -> A
C -> T
A -> G
T -> C
G -> C
C -> CAGA
A -> G
```


Type: transformation

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfafpath.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfafpath.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% SEGMENTIHS(1) segmentIhs (vcflib) | segmentIhs (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**segmentIhs**

# SYNOPSIS

segmentFst -s 2 -f iHS.normalized.output.txt

# DESCRIPTION

Creates genomic segments (bed file) for regions with high wcFst



# OPTIONS

```

Output : 8 columns :                 
     1. Seqid                        
     2. Start (zero based)           
     3. End   (zero based)           
     4. Average iHS                  
     5. Average high Fst (iHS > -s)  
     6. N iHS values in segment      
     7. N high iHS values in segment 
     8. Segment length               
required: -f            -- Output from normalizeIHS     
optional: -s            -- High absolute iHS cutoff [2] 

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[segmentIhs.cpp](https://github.com/vcflib/vcflib/blob/master/src/segmentIhs.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFHETHOMRATIO(1) vcfhethomratio (vcflib) | vcfhethomratio (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfhethomratio**

# SYNOPSIS

**vcfhethomratio** <vcf file>

# DESCRIPTION

Generates the het/hom ratio for each individual in the file



# OPTIONS

```


Type: metrics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfhethomratio.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfhethomratio.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFINFO2QUAL(1) vcfinfo2qual (vcflib) | vcfinfo2qual (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfinfo2qual**

# SYNOPSIS

**vcfinfo2qual** [key] [vcf_file]

# DESCRIPTION

Sets QUAL from info field tag keyed by [key]. The VCF file may be omitted and read from stdin. The average of the field is used if it contains multiple values.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfinfo2qual.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfinfo2qual.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSTREAMSORT(1) vcfstreamsort (vcflib) | vcfstreamsort (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfstreamsort**

# SYNOPSIS

**vcfstreamsort** [options] [vcf file]

# DESCRIPTION

Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart.



# OPTIONS

```

options:

    -h, --help             this dialog
    -w, --window N         number of sites to sort (default 10000)
    -a, --all              load all sites and then sort in memory

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfstreamsort.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfstreamsort.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFFIXUP(1) vcffixup (vcflib) | vcffixup (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcffixup**

# SYNOPSIS

**vcffixup** <vcf file>

# DESCRIPTION

Generates a VCF stream where AC and NS have been generated for each record using sample genotypes



# OPTIONS

```




Count the allele frequencies across alleles present in each record in the VCF file. (Similar to vcftools --freq.)

Uses genotypes from the VCF file to correct AC (alternate allele count), AF
(alternate allele frequency), NS (number of called), in the VCF records.  For
example:

    % vcfkeepsamples file.vcf NA12878 | **vcffixup** - | vcffilter -f "AC > 0"

Would downsample file.vcf to only NA12878, removing sites for which the sample
was not called as polymorphic.

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcffixup.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcffixup.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCLEANCOMPLEX(1) vcfcleancomplex (vcflib) | vcfcleancomplex (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcleancomplex**

# SYNOPSIS

**vcfcleancomplex** <vcf file>

# DESCRIPTION

Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change.



# OPTIONS

```


Generate a VCF stream in which 'long' non-complexalleles have their position corrected.
assumes that VCF records can't overlap 5'->3'

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcleancomplex.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcleancomplex.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% ABBA-BABA(1) abba-baba (vcflib) | abba-baba (VCF genotype)
% Erik Garrison and vcflib contributors

# NAME

**abba-baba**

# SYNOPSIS

**abba-baba** --tree 0,1,2,3 --file my.vcf --type PL

# DESCRIPTION

**abba-baba** calculates the tree pattern for four indviduals. This tool assumes reference is ancestral and ignores non **abba-baba** sites. The output is a boolian value: 1 = true , 0 = false for abba and baba. the tree argument should be specified from the most basal taxa to the most derived.



# OPTIONS

```


     Example:
     D   C  B   A 
     \ /  /    /  
      \  /    /   
       \    /    
        \  /     
         /        
        /         
 --tree A,B,C,D

Output : 4 columns :     
     1. seqid            
     2. position         
     3. abba             
     4. baba             
required: t,tree       -- a zero based comma separated list of target individuals corrisponding to VCF columns
required: f,file       -- a properly formatted VCF.                                                           
required: y,type       -- genotype likelihood format ; genotypes: GP,GL or PL;                                


type: genotype

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[abba-baba.cpp](https://github.com/vcflib/vcflib/blob/master/src/abba-baba.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENOSUMMARIZE(1) vcfgenosummarize (vcflib) | vcfgenosummarize (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenosummarize**

# SYNOPSIS

**vcfgenosummarize** <[input file] >[output vcf]

# DESCRIPTION

Adds summary statistics to each record summarizing qualities reported in called genotypes. Uses: RO (reference observation count), QR (quality sum reference observations) AO (alternate observation count), QA (quality sum alternate observations)



# OPTIONS

```


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgenosummarize.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenosummarize.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% BFST(1) bFst (vcflib) | bFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**bFst**

# SYNOPSIS

**bFst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1

# DESCRIPTION

**bFst** is a Bayesian approach to Fst. Importantly **bFst** accounts for genotype uncertainty in the model using genotype likelihoods. For a more detailed description see: `A Bayesian approach to inferring population structure from dominant markers' by Holsinger et al. Molecular Ecology Vol 11, issue 7 2002. The likelihood function has been modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population's allele frequency, and Fst.



# OPTIONS

```


Output : 11 columns :                          
     1.  Seqid                                     
     2.  Position				     
     3.  Observed allele frequency in target.	     
     4.  Estimated allele frequency in target.     
     5.  Observed allele frequency in background.  
     6.  Estimated allele frequency in background. 
     7.  Observed allele frequency combined. 	     
     8.  Estimated allele frequency in combined.   
     9.  ML estimate of Fst (mean)		     
     10. Lower bound of the 95% credible interval  
     11. Upper bound of the 95% credible interval  

required: t,target     -- a zero bases comma separated list of target individuals corrisponding to VCF columns
required: b,background -- a zero bases comma separated list of background individuals corrisponding to VCF columns
required: f,file a     -- a proper formatted VCF file.  the FORMAT field MUST contain "PL"
required: d,deltaaf    -- skip sites were the difference in allele frequency is less than deltaaf


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[bFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/bFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% SEGMENTFST(1) segmentFst (vcflib) | segmentFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**segmentFst**

# SYNOPSIS

**segmentFst** -s 0.7 -f wcFst.output.txt

# DESCRIPTION

**segmentFst** creates genomic segments (bed file) for regions with high wcFst



# OPTIONS

```


**segmentFst** provides a way to find continious regions with high Fst values.  It takes the output of wcFst and produces a BED file.  These high Fst region can be permutated with 'permuteGPATwindow'
Output : 8 columns :                 
     1. Seqid                        
     2. Start (zero based)           
     3. End   (zero based)           
     4. Average Fst                  
     5. Average high Fst (Fst > -s)  
     6. N Fst values in segment      
     7. N high fst values in segment 
     8. Segment length               
required: -f            -- Output from wcFst     
optional: -s            -- High Fst cutoff [0.8] 

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[segmentFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/segmentFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFQUAL2INFO(1) vcfqual2info (vcflib) | vcfqual2info (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfqual2info**

# SYNOPSIS

**vcfqual2info** [key] [vcf_file]

# DESCRIPTION

Puts QUAL into an info field tag keyed by [key].



# OPTIONS

```



Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfqual2info.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfqual2info.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENO2HAPLO(1) vcfgeno2haplo (vcflib) | vcfgeno2haplo (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfgeno2haplo**

# SYNOPSIS

**vcfgeno2haplo** [options] [<vcf file>]

# DESCRIPTION

Convert genotype-based phased alleles within --window-size into haplotype alleles. Will break haplotype construction when encountering non-phased genotypes on input.



# OPTIONS

```

options:
    -h, --help              Print this message
    -v, --version           Print version
    -r, --reference FILE    FASTA reference file
    -w, --window-size N     Merge variants at most this many bp apart (default 30)
    -o, --only-variants     Don't output the entire haplotype, just concatenate
                            REF/ALT strings (delimited by ":")



Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgeno2haplo.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgeno2haplo.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFRANDOM(1) vcfrandom (vcflib) | vcfrandom (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfrandom**

# SYNOPSIS

**vcfrandom**

# DESCRIPTION

Generate a random VCF file





# EXAMPLES

```

Example:

    **vcfrandom**

##fileformat=VCFv4.0
##source=**vcfrandom**
##reference=/d2/data/references/build_37/human_reference_v37.fa
##phasing=none
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bill
one     1       .       G       G,A     100     .       DP=83   GT:DP   0/1:1
one     2       .       G       G,A     100     .       DP=3    GT:DP   0/1:49
one     3       .       G       C,T     100     .       DP=5    GT:DP   0/1:12
one     4       .       C       G,T     100     .       DP=51   GT:DP   0/1:60
one     5       .       A       T,A     100     .       DP=31   GT:DP   0/1:89
one     6       .       T       T,A     100     .       DP=56   GT:DP   0/1:60
one     7       .       T       A,C     100     .       DP=78   GT:DP   0/1:75
one     8       .       T       G,A     100     .       DP=73   GT:DP   0/1:78
one     9       .       C       C,G     100     .       DP=42   GT:DP   0/1:67


Type: statistics

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfrandom.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfrandom.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFLEFTALIGN(1) vcfleftalign (vcflib) | vcfleftalign (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfleftalign**

# SYNOPSIS

**vcfleftalign** [options] [file]

# DESCRIPTION

Left-align indels and complex variants in the input using a pairwise ref/alt alignment followed by a heuristic, iterative left realignment process that shifts indel representations to their absolute leftmost (5') extent.



# OPTIONS

```


This is the same procedure used in the internal left alignment in
freebayes, and can be used when preparing VCF files for input to
freebayes to decrease positional representation differences between
the input alleles and left-realigned alignments.

options:

        -r, --reference FILE  Use this reference as a basis for realignment.
        -w, --window N        Use a window of this many bp when left aligning (150).

Left-aligns variants in the specified input file or stdin.  Window
size is determined dynamically according to the entropy of the regions
flanking the indel.  These must have entropy > 1 bit/bp, or be shorter
than ~5kb.


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfleftalign.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfleftalign.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFKEEPINFO(1) vcfkeepinfo (vcflib) | vcfkeepinfo (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfkeepinfo**

# SYNOPSIS

**vcfkeepinfo** <vcf file> [FIELD1] [FIELD2] ...

# DESCRIPTION

To decrease file size remove INFO fields not listed on the command line



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfkeepinfo.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfkeepinfo.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSAMPLE2INFO(1) vcfsample2info (vcflib) | vcfsample2info (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfsample2info**

# SYNOPSIS

**vcfsample2info** [options] <vcf file>

# DESCRIPTION

Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO.



# OPTIONS

```

options:
    -f, --field         Add information about this field in samples to INFO column
    -i, --info          Store the computed statistic in this info field
    -a, --average       Take the mean of samples for field (default)
    -m, --median        Use the median
    -n, --min           Use the min
    -x, --max           Use the max


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfsample2info.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsample2info.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% GENOTYPESUMMARY(1) genotypeSummary (vcflib) | genotypeSummary (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**genotypeSummary**

# SYNOPSIS

genotypeSummmary --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf --snp

# DESCRIPTION

Generates a table of genotype counts. Summarizes genotype counts for bi-allelic SNVs and indel



# OPTIONS

```


output: table of genotype counts for each individual.
required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns        
required: f,file       -- proper formatted VCF                                                                        
required, y,type       -- genotype likelihood format; genotype : GL,PL,GP                                             
optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1                                              
optional, s,snp        -- Only count SNPs                                              
optional, a,ancestral  -- describe counts relative to the ancestral allele defined as AA in INFO

Type: statistics


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[genotypeSummary.cpp](https://github.com/vcflib/vcflib/blob/master/src/genotypeSummary.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% DUMPCONTIGSFROMHEADER(1) dumpContigsFromHeader (vcflib) | dumpContigsFromHeader (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**dumpContigsFromHeader**

# SYNOPSIS

**dumpContigsFromHeader** file

# DESCRIPTION

Dump contigs from header





# EXAMPLES

```

Example:

    **dumpContigsFromHeader** samples/scaffold612.vcf

    ##contig=<ID=scaffold4,length=1524>
    ##contig=<ID=scaffold12,length=56895>
    (...)

    output

    scaffold4       1524
    scaffold12      56895
    (...)

Type: transformation
      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[dumpContigsFromHeader.cpp](https://github.com/vcflib/vcflib/blob/master/src/dumpContigsFromHeader.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFREMAP(1) vcfremap (vcflib) | vcfremap (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfremap**

# SYNOPSIS

**vcfremap** [options] [<vcf file>]

# DESCRIPTION

For each alternate allele, attempt to realign against the reference with lowered gap open penalty. If realignment is possible, adjust the cigar and reference/alternate alleles. Observe how different alignment parameters, including context and entropy-dependent ones, influence variant classification and interpretation.



# OPTIONS

```

options:
    -w, --ref-window-size N      align using this many bases flanking each side of the reference allele
    -s, --alt-window-size N      align using this many flanking bases from the reference around each alternate allele
    -r, --reference FILE         FASTA reference file, required with -i and -u
    -m, --match-score N          match score for SW algorithm
    -x, --mismatch-score N       mismatch score for SW algorithm
    -o, --gap-open-penalty N     gap open penalty for SW algorithm
    -e, --gap-extend-penalty N   gap extension penalty for SW algorithm
    -z, --entropy-gap-open       use entropy scaling for the gap open penalty
    -R, --repeat-gap-extend N    penalize non-repeat-unit gaps in repeat sequence
    -a, --adjust-vcf TAG         supply a new cigar as TAG in the output VCF


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfremap.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfremap.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% MELTEHH(1) meltEHH (vcflib) | meltEHH (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**meltEHH**

# SYNOPSIS

**meltEHH** --target 0,1,2,3,4,5,6,7 --pos 10 --file my.phased.vcf \ --region chr1:1-1000 > STDOUT 2> STDERR

# DESCRIPTION





# OPTIONS

```


 **meltEHH** provides the data to plot extended haplotype homozygosity
(EHH) curves and produces the data to generate the following plot:
<img src="https://github.com/vcflib/vcflib/blob/master/examples/example-ehh.png?raw=true" alt="" width=400>



     **meltEHH** provides the data to plot EHH curves.
Output : 4 columns :
     1. seqid
     2. position
     3. EHH
     4. ref or alt [0 == ref]
Params:
       required: t,target   <STRING>  A zero base comma separated list of target
                                     individuals corresponding to VCF columns
       required: r,region   <STRING>  A tabix compliant genomic range
                                     format: "seqid:start-end" or "seqid"
       required: f,file     <STRING>  Proper formatted and phased VCF.
       required: y,type     <STRING>  Genotype likelihood format: GT,PL,GL,GP
       required: p,position <INT>     Variant position to melt.
       optional: a,af       <DOUBLE>  Alternative alleles with frequencies less
                                     than [0.05] are skipped.



Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[meltEHH.cpp](https://github.com/vcflib/vcflib/blob/master/src/meltEHH.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCREATEMULTI(1) vcfcreatemulti (vcflib) | vcfcreatemulti (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcreatemulti**

# SYNOPSIS

**vcfcreatemulti** [options] [file]

# DESCRIPTION

If overlapping alleles are represented across multiple records, merge them into a single record. Currently only for indels.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcreatemulti.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcreatemulti.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFKEEPGENO(1) vcfkeepgeno (vcflib) | vcfkeepgeno (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfkeepgeno**

# SYNOPSIS

**vcfkeepgeno** <vcf file> [FIELD1] [FIELD2] ...

# DESCRIPTION

Reduce file size by removing FORMAT fields not listed on the command line from sample specifications in the output



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfkeepgeno.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfkeepgeno.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFOVERLAY(1) vcfoverlay (vcflib) | vcfoverlay (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfoverlay**

# SYNOPSIS

**vcfoverlay** [options] [<vcf file> ...]

# DESCRIPTION

Overlay records in the input vcf files with order as precedence.



# OPTIONS

```

options:
    -h, --help       this dialog
    -v, --version    prints version


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfoverlay.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfoverlay.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCF2TSV(1) vcf2tsv (vcflib) | vcf2tsv (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcf2tsv**

# SYNOPSIS

**vcf2tsv** [-n null_string] [-g] [vcf file]

# DESCRIPTION

Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table. Specifying -g will output one line per sample with genotype information. When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcf2tsv.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2tsv.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFINDEX(1) vcfindex (vcflib) | vcfindex (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfindex**

# SYNOPSIS

**vcfindex** <vcf file>

# DESCRIPTION

Adds an index number to the INFO field (id=position)



# OPTIONS

```


Type: transformation

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfindex.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfindex.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCOUNTALLELES(1) vcfcountalleles (vcflib) | vcfcountalleles (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfcountalleles**

# SYNOPSIS

**vcfcountalleles** <vcf file>

# DESCRIPTION

Count alleles





# EXAMPLES

```

Example:

**vcfcountalleles** samples/scaffold612.vcf
42603

Type: statistics

      

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcountalleles.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcountalleles.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFSAMPLEDIFF(1) vcfsamplediff (vcflib) | vcfsamplediff (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfsamplediff**

# SYNOPSIS

**vcfsamplediff** [options] <tag> <sample> <sample> [ <sample> ... ] <vcf file>

# DESCRIPTION

Establish putative somatic variants using reported differences between germline and somatic samples. Tags each record where the listed sample genotypes differ with <tag>. The first sample is assumed to be germline, the second somatic. Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples.



# OPTIONS

```


options:
    -s --strict     Require that no observations in the germline support the somatic alternate.


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfsamplediff.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfsamplediff.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFANNOTATEGENOTYPES(1) vcfannotategenotypes (vcflib) | vcfannotategenotypes (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfannotategenotypes**

# SYNOPSIS

**vcfannotategenotypes** <annotation-tag> <vcf file> <vcf file>

# DESCRIPTION

Examine genotype correspondence. Annotate genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant.



# OPTIONS

```


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfannotategenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfannotategenotypes.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCF2FASTA(1) vcf2fasta (vcflib) | vcf2fasta (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcf2fasta**

# SYNOPSIS

**vcf2fasta** [options] [file]

# DESCRIPTION

Generates sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy]. Each sequence in the fasta file is named using the same pattern used for the file name, allowing them to be combined.



# OPTIONS

```

options:
    -f, --reference REF     Use this reference when decomposing samples.
    -p, --prefix PREFIX     Affix this output prefix to each file, none by default
    -P, --default-ploidy N  Set a default ploidy for samples which do not have information in the first record (2).


Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcf2fasta.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2fasta.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% PFST(1) pFst (vcflib) | pFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**pFst**

# SYNOPSIS

**pFst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

# DESCRIPTION

**pFst** is a probabilistic approach for detecting differences in allele frequencies between two populations.



# OPTIONS

```




**pFst** is a likelihood ratio test (LRT) quantifying allele frequency
differences between populations.  The LRT by default uses the binomial
distribution.  If Genotype likelihoods are provided it uses a modified
binomial that weights each allele count by its certainty.  If type is
set to 'PO' the LRT uses a beta distribution to fit the allele
frequency spectrum of the target and background.  PO requires the AD
and DP genotype fields and requires at least two pools for the target
and background.  The p-value calculated in **pFst** is based on the
chi-squared distribution with one degree of freedom.


Output : 3 columns :     
     1. seqid            
     2. position         
     3. **pFst** probability 

required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       
required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   
required: f,file       -- argument: a properly formatted VCF.                                                                  
required: y,type       -- argument: genotype likelihood format ; genotypes: GP, GL or PL; pooled: PO                           
optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero
optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 
optional: c,counts     -- switch  : use genotype counts rather than genotype likelihoods to estimate parameters, default false 

Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[pFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/pFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFGENOTYPES(1) vcfgenotypes (vcflib) | vcfgenotypes (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenotypes**

# SYNOPSIS

**vcfgenotypes** <vcf file>

# DESCRIPTION

Report the genotypes for each sample, for each variant in the VCF. Convert the numerical represenation of genotypes provided by the GT field to a human-readable genotype format.



# OPTIONS

```




```





# EXAMPLES

```

Example:

      **vcfgenotypes** samples/sample.vcf

19      111     A       C       A,C     NA00001:A/A     NA00002:A/A     NA00003:A/C
19      112     A       G       A,G     NA00001:A/A     NA00002:A/A     NA00003:A/G
20      14370   G       A       G,A     NA00001:G/G     NA00002:G/A     NA00003:A/A
20      17330   T       A       T,A     NA00001:T/T     NA00002:T/A     NA00003:T/T
20      1110696 A       G,T     A,G,T   NA00001:G/T     NA00002:G/T     NA00003:T/T
20      1230237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:T/T
20      1234567 G       GA,GAC  G,GA,GAC        NA00001:G/GA    NA00002:G/GAC   NA00003:GA/GA
20      1235237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:./.
X       10      AC      A,ATG   AC,A,ATG        NA00001:AC      NA00002:AC/A    NA00003:AC/ATG

Type: statistics

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenotypes.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% VCFCAT(1) vcfcat (vcflib) | vcfcat (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfcat**

# SYNOPSIS

**vcfcat** [file1] [file2] ... [fileN]

# DESCRIPTION

Concatenates VCF files



# OPTIONS

```


Type: transformation

      

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfcat.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcat.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
% PERMUTEGPAT++(1) permuteGPAT++ (vcflib) | permuteGPAT++ (VCF phenotype)
% Erik Garrison and vcflib contributors

# NAME

**permuteGPAT++**

# SYNOPSIS

**permuteGPAT++** -f gpat.txt -n 5 -s 1

# DESCRIPTION

**permuteGPAT++** is a method for adding empirical p-values to a GPAT++ score.



# OPTIONS

```


     Currently **permuteGPAT++** only supports wcFst, but will be extended.    

OUTPUT: **permuteGPAT++** will append three additional columns:
        1. The number of successes                         
        2. The number of trials                            
        3. The empirical p-value                           

file:    f   -- argument: the input file     
number:  n   -- argument: the number of permutations to run for each value [1000]
success: s   -- argument: stop permutations after 's' successes [1]

Type: phenotype

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[permuteGPAT++.cpp](https://github.com/vcflib/vcflib/blob/master/src/permuteGPAT++.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
