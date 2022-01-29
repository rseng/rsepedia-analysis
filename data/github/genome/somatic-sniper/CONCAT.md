# SomaticSniper

   The purpose of this program is to identify single nucleotide positions that are different between tumor and normal (or in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. Complete documentation is available at the project [web site](http://gmt.genome.wustl.edu/somatic-sniper/) or in the included [pdf](https://github.com/genome/somatic-sniper/blob/master/docs/sniper_manual.pdf).

## License
SomaticSniper is licensed under the [MIT license](docs/copyright). 

## Documentation
* [Installation](gmt/install.md)
* [Usage](gmt/documentation.md)
* [FAQ](gmt/install.md#faq)

## User Support

Please first search [Biostar](http://www.biostars.org) and then ask a question there if needed.  We automatically monitor [Biostar](http://www.biostars.org) for questions related to our tools.

## SomaticSniper User Manual

__David E. Larson, Travis E. Abbott and Christopher C. Harris__

__October 26, 2011__


The purpose of this program is to identify single nucleotide positions that are different between tumor and normal (or, in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. It outputs a file in a format very similar to Samtools consensus format. It uses the genotype likelihood model of MAQ (as implemented in Samtools) and then calculates the probability that the tumor and normal genotypes are different. This probability is reported as a somatic score. The somatic score is the Phred-scaled probability (between 0 to 255) that the Tumor and Normal genotypes are not different where 0 means there is no probability that the genotypes are different and 255 means there is a probability of 1 − 10^(255∕−10) that the genotypes are different between tumor and normal. This is consistent with how the SAM format reports such probabilities.

There are two modes, the joint genotyping mode (-J) takes into account the fact that the tumor and normal samples are not entirely independent and also takes into account the prior probability of a somatic mutation. This probability can be scaled to control the sensitivity of the algorithm. An accurate value for this prior would be 0.000001, but this may result in a severe lack of sensitivity at lower depths. A less realistic prior probability will generate more sensitive results at the expense of an increase in the number of false positives. To get a similar sensitivity to the default mode, we recommend using a prior of 0.01. The default mode treats the two samples as if they came from two different individuals. This mode uses a less accurate mathematical model, but yields good results, especially if the normal may contain some tumor cells or the tumor is quite impure.

## Usage

`bam-somaticsniper [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snv_output_file>`

__Required Option:__

`-f   FILE REQUIRED reference sequence in the FASTA format`

__Options:__

`-q  INT filtering reads with mapping quality less than INT [0]`

`-Q  INT filtering somatic snv output with somatic quality less than INT [15]`

`-L FLAG do not report LOH variants as determined by genotypes`

`-G FLAG do not report Gain of Referene variants as determined by genotypes`

`-p  FLAG disable priors in the somatic calculation. Increases sensitivity for solid tumors.`

`-J  FLAG Use prior probabilities accounting for the somatic mutation rate`

`-s  FLOAT prior probability of a somatic mutation (implies -J) [0.01]`

`-T  FLOAT theta in maq consensus calling model (for -c/-g) [0.850000]`

`-N  INT number of haplotypes in the sample (for -c/-g) [2]`

`-r  FLOAT prior of a diﬀerence between two haplotypes (for -c/-g) [0.001000]`

`-F  STRING select output format (vcf or classic) [classic]`

### Notes on running SomaticSniper

Minimally, you must provide the program the reference fasta the bams were aligned against (passed with the -f option), a tumor bam, a normal bam, and the filename of the resulting output file. We recommend filtering out reads with a mapping quality of 0 (i.e. use -q 1) as they are typically randomly placed in the genome. We have also found that few variants with a somatic score less than 15 validate, but you may decrease the minimum score or increase it to a higher threshold (eg -Q 40). To obtain high confidence sites, we recommend also thresholding the minimum average mapping quality for the variant base to 40 for reads aligned with BWA or 70 for reads aligned with MAQ. We have not tested other aligners at this time. Disabling priors is not recommended, but may increase sensitivity at the cost of a decrease in specificity.

#### Current Recommended Settings

We recommend that you utilize both the `-G` and `-L` options when running in order to reduce likely false positives with little impact on sensitivity. An example command-line is below:

```
bam-somaticsniper -Q 40 -G -L -f reference.fa tumor.bam normal.bam output.txt
```

#### Basic filtering with provided Perl scripts

A small number of basic Perl scripts are included in the SomaticSniper package (located in src/scripts of the source code release) to aid in filtering out likely false positives. In order to get the recommended filtering you should do the following. Defaults are set assuming that BWA short is the aligner used. Other aligners have not been tested and recommendations are not available. Before proceeding you will need to obtain and compile bam-readcount (https://github.com/genome/bam-readcount). You will also need to generate a samtools pileup (not mpileup) indel file. Handling of indel containing VCFs is not implemented.

1. Filter on standard filters using the indel file. This will also remove LOH calls e.g. `perl snpfilter.pl –snp-file your_sniper_file –indel-file your_indel_pileup`

2. Adapt the remainder for use with bam-readcount e.g. `perl prepare_for_readcount.pl –snp-file your_sniper_file.SNPfilter`

3. Run bam-readcount (I’d recommend using the same mapping quality -q setting as you ran SomaticSniper with) e.g. `bam-readcount -b 15 -f your_ref.fasta -l your_sniper_file.SNPfilter.pos your_tumor.bam > your_readcounts.rc`
    Run the false positive filter e.g. perl fpfilter.pl –snp-file your_sniper_file.SNPfilter –readcount-file your_readcounts.rc

5. Lastly, run the "high confidence" filter which filters based on the Somatic Score and mapping quality e.g. `perl highconfidence.pl –snp-file your_sniper_file.SNPfilter.fp_pass`

Your final set of high confidence and highly filtered indels is now in the file `your_sniper_file.SNPfilter.fp_pass.hc`

### File Formats

The output by SomaticSniper consists of line for all sites whose consensus differs from the reference base. Each of the three available output formats is described below

__Classic:__

Each line contains the following tab-separated values:

 1. Chromosome
 2. Position
 3. Reference base
 4. IUB genotype of tumor
 5. IUB genotype of normal
 6. Somatic Score
 7. Tumor Consensus quality
 8. Tumor variant allele quality
 9. Tumor mean mapping quality
 10. Normal Consensus quality
 11. Normal variant allele quality
 12. Normal mean mapping quality
 13. Depth in tumor (# of reads crossing the position)
 14. Depth in normal (# of reads crossing the position)
 15. Mean base quality of reads supporting reference in tumor
 16. Mean mapping quality of reads supporting reference in tumor
 17. Depth of reads supporting reference in tumor
 18. Mean base quality of reads supporting variant(s) in tumor
 19. Mean mapping quality of reads supporting variant(s) in tumor
 20. Depth of reads supporting variant(s) in tumor
 21. Mean base quality of reads supporting reference in normal
 22. Mean mapping quality of reads supporting reference in normal
 23. Depth of reads supporting reference in normal
 24. Mean base quality of reads supporting variant(s) in normal
 25. Mean mapping quality of reads supporting variant(s) in normal
 26. Depth of reads supporting variant(s) in normal

__VCF__

VCF output from SomaticSniper conforms to version 4.1 of the VCF specification. Hence, each non-header output line contains the following fields:

 1. Chromosome
 2. Position
 3. ID (unused)
 4. Reference base
 5. Alternate bases (comma separated)
 6. Quality (unused)
 7. Filters (unused)
 8. INFO (unused)
 9. FORMAT specification for each sample
 10. NORMAL sample data
 11. TUMOR sample data

The following FORMAT fields will be populated for each of NORMAL and TUMOR.
 
| ID | Number | Type | Description |
| --- | ------ | ---- | -----------|
| GT | 1 | String | Genotype |
|IGT | 1 | String |Genotype when called independently (only filled if called in joint prior mode)|
| DP | 1 | Integer| Total read depth |
| DP4| 4 | Integer | Number of high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases |
| BCOUNT | 4 | Integer | Occurrence count for each base at this site (A,C,G,T) |
| GQ |1 | Integer | Genotype quality |
| JGQ | 1 | Integer | Joint genotype quality (only filled if called in joint prior mode) |
| VAQ | 1 | Integer | Variant quality |
| BQ  | . | Integer | Average base quality of each base in the call, reported in alphabetical order (A,C,G,T) | 
| MQ | 1 | Integer  | Average mapping quality across all reads. |
| AMQ | . | Integer | Average mapping quality of each base in the call, reported in alphabetical order (A,C,G,T) |
| SS | 1 | Integer  | Variant status relative to non-adjacent normal: 0=wildtype, 1=germline, 2=somatic, 3=LOH, 4=unknown |
| SSC | 1 | Integer | Somatic Score |

### User Support

Please first search [Biostar](http://www.biostars.org) and then ask a question there if needed. We automatically monitor Biostar for questions related to our tools. 
## Build Instructions

### Build dependencies

* For APT-based systems (Debian, Ubuntu), install the following packages:

```
sudo apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev
```

* For RPM-based systems (Fedora, CentOS, RHEL), install the following packages instead:

```
sudo yum groupinstall "Development tools" 
sudo yum install zlib-devel ncurses-devel cmake
```

Note that for some RPM based systems (like RHEL), you will need to install cmake 2.8 or greater yourself as the packaged version is much older.

### Clone the SomaticSniper repository

Clone the git repository.

```
git clone git://github.com/genome/somatic-sniper.git
```

### Build SomaticSniper

SomaticSniper does not support in-source builds. So create a subdirectory, enter it, build, and run tests:

```
mkdir somatic-sniper/build
cd somatic-sniper/build
cmake ../
make deps
make -j
make test
```

The binary `bam-somaticsniper` can then be found under `somatic-sniper/build/bin`. If you have administrative rights, then run `sudo make install` to install the tool for all users under `/usr/bin`.

## FAQ

### I get lots of compile errors indicating that files are missing. How do I fix this?

SomaticSniper requires that it be linked to an old version of samtools (v0.1.6). This typically happens because you have attempted to link to a newer version. As of version [1.0.3](https://github.com/genome/somatic-sniper/releases/tag/v1.0.3), SomaticSniper includes samtools as part of its build process and you do not need to download samtools yourself.

### I get errors from cmake about missing modules. How do I fix this?

As of commit [09ef624](https://github.com/genome/somatic-sniper/commit/09ef624e5bb275e0fd62396a14a878711e746cb9) or version [1.0.4](https://github.com/genome/somatic-sniper/releases/tag/v1.0.4), this should no longer be an issue and tarballs from github should function as intended. In earlier versions, SomaticSniper contained a git submodule called build-common. This submodule contains helper modules for cmake. If you downloaded the source as a tarball from github or forgot to do a recursive clone using git, then you will not have this submodule and will see cmake errors. If you are using git, we recommend you go back and use the `--recursive` option when cloning the SomaticSniper repository. If you cannot use git, follow the instructions below to remedy the situation.

1. Download the build-common module separately [here](https://github.com/genome/build-common/tarball/master).
2. Extract that tarball and rename the directory it creates to `build-common`.
3. Replace the empty build-common subdirectory in the sniper directory with directory you just created.
4. Resume following the build instructions.

### I get a floating point exception on running `bam-somaticsniper`. What's going on?

This has been reported when using reference fasta indexes available via the GATK resource bundle. Please try reindexing your fasta with samtools and rerunning `bam-somaticsniper`.
