cerebra
================================
<a href="https://pypi.org/project/cerebra/"><img alt="PyPI" src="https://badge.fury.io/py/cerebra.svg"></a>
[![Docker Build](https://img.shields.io/docker/cloud/build/lincolnharris/cerebra)](https://hub.docker.com/r/lincolnharris/cerebra)    
[![Build Status](https://travis-ci.org/czbiohub/cerebra.svg?branch=master)](https://travis-ci.org/czbiohub/cerebra) 
[![codecov](https://codecov.io/gh/czbiohub/cerebra/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/cerebra)              
[![DOI](https://zenodo.org/badge/180649463.svg)](https://zenodo.org/badge/latestdoi/180649463)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02432/status.svg)](https://doi.org/10.21105/joss.02432)

What is `cerebra`?
-------------------------------------
This tool allows you to quickly summarize and interpret VCF files.      

If you're interested in learning the full complement of genomic variants present in your DNA/RNA samples, a typical workflow might involve sequencing, followed by variant calling with a tool such as GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php), which generates [variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (VCF) files.       

A VCF file looks like this:

```##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
chr1 631391 . C T 72.28 . AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;QD=25.36;SOR=2.303 GT:AD:DP:GQ:PL 1/1:0,2:2:6:84,6,0
```
Note that only a single VCF record is displayed here.
A sequencing run can generate on the order of 10<sup>8</sup> unique VCF records. 
Only a small portion of these VCF records are likely to be relevant to the researcher.
Thus drawing conclusions from VCF files remains a substantial challenge.


`cerebra` provides a fast and intuitive framework for summarizing VCF records across samples.
It comprises three modules that do the following:      

        1) remove germline variants from samples of interest        
        2) count the total number of variants in a given sample, on a per-gene basis           
        3) report protein variants for each sample                 
        
`cerebra` gets its name from the eponymous X-men [character](https://en.wikipedia.org/wiki/Cerebra), who had the ability to detect mutant individuals among the general public. 

If you're working with tumor data, it might be a good idea to limit the variant search space to only known, and potentially actionable, cancer variants. 
Therefore `cerebra` implements an optional method for restricting to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database.  

This tool was developed for, but is certainly not limited to, single-cell RNA sequencing data. 
`cerebra` is free software available under the MIT license.


What makes `cerebra` different from traditional VCF parsers? 
-------------------------------------
Python libraries exist (_i.e._ [PyVCF](https://pyvcf.readthedocs.io/en/latest/) and [vcfpy](https://vcfpy.readthedocs.io/en/stable/index.html)) for extracting information from VCF files.
Another is [GATK VariantsToTable](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php), which produces a file that looks like this:
 
    CHROM    POS ID      QUAL    AC
     1        10  .       50      1
     1        20  rs10    99      10

This table contains only genomic (_i.e._ DNA-level) coordinates. 
Often the next questions are: what are the genes associated with these variants, and what are the peptide-level effects of these variants?
`cerebra` queries a reference genome (.fa) and annotation (.gtf) to match each DNA-level variant with its associated gene, and its probable protein variant.
`cerebra` produces the following outfile: 

```
$ python
> import json
> f = open(/path/to/cerebra/output.json)
> data = json.load(f)
> print(json.dumps(data, indent=4))

{
    "CCN1": {
        "A16_B000563": [],
        "A1_B001546": [],
        "A1_B002531": [
            "ENSP00000398736.2:p.(Glu189=)"
        ],
        "A1_B002570": [],
        "A2_B002558": [],
        "A3_B000561": [
            "ENSP00000398736.2:p.(Arg209Trp)",
            "ENSP00000398736.2:p.(Ile90Val)"
        ],
        "A3_B000568": [],
        "A3_B001544": [
            "ENSP00000398736.2:p.(Ala82Thr)"
        ],
        "A3_B002090": [],
        "A3_B002531": [
            "ENSP00000398736.2:p.(Pro217Ser)"
        ]
    }
}
```

Here _CCN1_ is a gene name while _A16_B000563_, _A1_B001546_, _A1_B002531,_... are RNA-seq sample IDs.
`cerebra` reports variants to every gene in the genome, for every sample in a given experiment. 
The _ENSP*_ numbers refer to [Ensembl](https://www.ensembl.org/index.html) translation IDs -- unique identifiers that correspond to exactly one polypeptide in the Ensembl database. 
The strings enclosed in parentheses refer to the amino-acid level variants reported in that particular sample. 
For example the string `Arg209Trp` indicates that position 209 of this particular polypeptide should contain an _Arg_, but the experimental sample instead contains a _Trp_. 
`cerebra` adheres to HGVS sequence variant [nomenclature](https://varnomen.hgvs.org/) in reporting amino-acid variants.


Features
--------
### `germline-filter`

Variant calling is often applied to the study of cancer. 
If the research project is centered around a “tumor vs. normal” question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between tumor and normal samples, and thus excludes variants unlikely to be pathogenic for the cancer under study.
The user provides a very simple metadata file (see [USAGE.md](https://github.com/czbiohub/cerebra/blob/master/docs/USAGE.md)) that indicates which tumor samples correspond to which normal samples.
The output of germline-filter is a set of trimmed-down VCF files, which will be used for the next two steps.
If you do not have access to “normal” samples then proceed directly to `count-variants` or `find-peptide-variants`.

### `count-variants`
This module reports the raw variant counts for every gene across every sample.
The output is a CSV file that contains counts for each sample versus every gene in the genome. 
If working with cancer variants, the user has the option of limiting the search space to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database. 

### `find-peptide-variants`
This module reports the protein variations associated with each genomic variant.
VCF records are converted to peptide-level variants, and then [ENSEMBL](https://uswest.ensembl.org/index.html) protein IDs, 
in accordance to the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature. 
The output is a hierarchically ordered text file (CSV or JSON) that reports the Ensemble protein ID and the gene associated with each variant, for each experimental sample. 
The user again has the option of limiting the search space to variants also found in the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database. 
This module also has an option to report the number of variant vs. wildtype reads found at each locus. 


Dependencies
------------
`cerebra` depends on some (fairly standard) packages and libraries. 
Before installing, it might be a good idea to make sure all of the requisite packages are installed on your system (_note:_ if installing with Docker you can skip this step). 

__MacOS Dependencies:__
```
sudo pip install setuptools
brew update
brew install openssl
brew install zlib
```

__Linux Dependencies:__
```
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

As of present `cerebra` is not installable on Windows. 
`cerebra` depends on the [`pysam`](https://pysam.readthedocs.io/en/latest/index.html) library (or rather, `pysam` is a dependency-of-a-dependency) and currently this library is only available on Unix-like systems. 
Windows solutions like [WSL](https://docs.microsoft.com/en-us/windows/wsl/) exist for overcoming precisely this challenge, however, `cerebra` has not been tested on WSL or any other Unix-like subsystem for Windows.    


Installation (for users)
------------
There are four different methods available to install `cerebra`.
Choose one of the following:

__With [Docker](https://hub.docker.com/r/lincolnharris/cerebra) (recommended)__          
```
docker pull lincolnharris/cerebra
docker run -it lincolnharris/cerebra
```      
_warning: this image will take up ~1Gb of space._               

__With the python standard library [`venv`](https://docs.python.org/3/library/venv.html) module__
```
python3 -m venv cerebra
source cerebra/bin/activate
pip3 install cerebra 
```

__With [conda](https://docs.conda.io/en/latest/)__
``` 
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install cerebra
```

__From [PyPi](https://pypi.org/project/cerebra/) (system-wide installation, NOT RECOMMENDED)__    
For novice users, it's generally a better idea to install packages within virtual environments. 
However, `cerebra` can be installed system-wide, if you so choose. 
```
pip3 install cerebra

# OR, if you dont have root privileges
pip3 install --user cerebra
```


Installation (for developers)
------------
Here's how to set up cerebra for local development. 
After installing the requisite [dependencies](https://github.com/czbiohub/cerebra/blob/master/README.md#dependencies):

1.  Fork the `cerebra` repo on GitHub: https://github.com/czbiohub/cerebra
2.  Clone your fork locally:

        $ git clone https://github.com/your-name/cerebra.git

3.  Install your local copy into a virtualenv. Using the standard library [`venv`](https://docs.python.org/3/library/venv.html) module: 

        $ cd cerebra
        $ python3 -m venv cerebra-dev
        $ source cerebra-dev/bin/activate
        $ pip3 install -r requirements.txt -r test_requirements.txt -e .

4.  Create a branch for local development:

        $ git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5.  When you're done making changes, check that your changes pass `flake8` and `pytest`:

        $ make test
        $ make coverage
        $ make lint

6.  Commit your changes and push your branch to GitHub:

        $ git add .
        $ git commit -m "Your detailed description of your changes."
        $ git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.
See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for more. 


(Quickstart) Usage
-----
`$ cerebra` should return usage information

```
Usage: cerebra  <command>

  a tool for fast and accurate summarizing of variant calling format (VCF)
  files

Options:
  -h, --help  Show this message and exit.

Commands:
  germline-filter    filter out common SNPs/indels between tumor and normal samples
  count-variants    count total number of variants in each sample, and report on a per-gene basis
  find-peptide-variants  report peptide-level SNPs and indels in each sample, and associated coverage
```

Note that the `-h` command will display usage information for each of the three commands. 

An example workflow might look like this:   

**Step 1:**     
```
cerebra germline-filter --processes 2 --normal_path /path/to/normal/vcfs --tumor_path /path/to/tumor/vcfs --metadata /path/to/metadata/file --outdir /path/to/filtered/vcfs
```

**Step 2:**     
```
cerebra count-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --refgenome /path/to/genome/annotation --outfile /path/to/output/file /path/to/filtered/vcfs/*
```

**Step 3:**          
```
cerebra find-peptide-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --annotation /path/to/genome/annotation --genomefa /path/to/genome/fasta --report_coverage 1 --output /path/to/output/file /path/to/filtered/vcfs/*
```

For advanced usage information, see [USAGE.md](https://github.com/czbiohub/cerebra/blob/master/docs/USAGE.md). 


Authors
--------
This work was produced by [Lincoln Harris](https://github.com/lincoln-harris), [Rohan Vanheusden](https://github.com/rvanheusden), [Olga Botvinnik](https://github.com/olgabot) and [Spyros Darmanis](https://spyrosdarmanis.wixsite.com/mylab) of the Chan Zuckerberg Biohub. 
For questions please contact ljharris018@gmail.com


Contributing
--------
We welcome any bug reports, feature requests or other contributions. 
Please submit a well documented report on our [issue tracker](https://github.com/czbiohub/cerebra/issues). 
For substantial changes please fork this repo and submit a pull request for review. 

See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for additional details. 

You can find official releases [here](https://github.com/czbiohub/cerebra/releases). 
---
title: '`cerebra`: A tool for fast and accurate summarizing of variant calling format (VCF) files'
tags:
    - python
    - genomics
    - variant calling
    - VCF
    - single-cell
    - cancer
authors:
 - name: Lincoln Harris
   orcid: 0000-0003-2544-4225
   affiliation: 1
 - name: Rohan Vanheusden
   affiliation: 1 
 - name: Olga Botvinnik
   orcid: 0000-0003-4412-7970
   affiliation: 1 
 - name: Spyros Darmanis
   orcid: 0000-0003-4002-8158
   affiliation: 1 
affiliations: 
- name: Chan Zuckerberg Biohub, San Francisco, CA, United States
  index: 1
date: 20, August, 2020 
bibliography: paper.bib
---

## Motivation

A single "typo" in the genome can have profound consequences on an organism's biology.
Identifying the protein changes that accompany these genomic typos is a fundamental challenge in bioinformatics. 
Tools exist for identifying genomic variants and predicting their associated polypeptide-level changes, however, wrangling genomic variant calls and protein variant predictions across thousands of samples remains a substantial challenge. 
`cerebra` addresses this need by offering a fast and accurate framework for summarizing genomic variant calls and protein variant predictions across many samples. 

To find variants in the genome, researchers often begin with a [DNA-sequencing](https://en.wikipedia.org/wiki/DNA_sequencing) (DNA-seq) or [RNA-sequencing](https://en.wikipedia.org/wiki/RNA-Seq) (RNA-seq) experiment on their samples of interest.
Sequencing is followed by alignment of reads to the reference genome with tools like [STAR](https://github.com/alexdobin/STAR) or [BWA](http://bio-bwa.sourceforge.net/), followed by variant calling with tools like [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) 
or [freebayes](https://github.com/ekg/freebayes) [@star; @bwa; @haplocaller; @freebayes]. 
Variant callers produce tab delimited text files in the [variant calling format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (VCF) for each processed sample.
VCF files encode: i) the genomic position, ii) reference vs. observed DNA sequence, and iii) quality associated with each observed variant. 
Shown below are the first 4 lines of a sample VCF file. 
Note that only a single record is displayed, and that the record line has been artificially wrapped.

```
##fileformat=VCFv4.2
##source=HaplotypeCaller
#CHROM  POS ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
chr1	631391	.	C	T	72.28	.	AC=2;AF=1.00;AN=2;DP=2;
        ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=NaN;
        QD=25.36;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:84,6,0
```

Current methods for variant calling are incredibly powerful and robust, however, a single sequencing run can generate as many as $10^{8}$ unique VCF records.
Only a small portion of these VCF records are likely to be relevant to the researcher.
In addition, variant callers report only the genomic location of the variant, and not the effect the variant has on the translated protein sequence.
To address the unmet need for high-throughput VCF summary tools, we introduce `cerebra`, a python package that provides fast and accurate protein variant summarizing of VCF files.


## Functionality

`cerebra` comprises three modules: i) `germline-filter` is intended for working with cancer data and removes variants that are common between tumor and normal samples, ii) `count-variants` reports total number of variants in each sample, and iii) `find-peptide-variants` reports the likely protein variants in each sample. 
Here we use _variant_ to refer to single nucleotide polymorphisms (SNPs) and short insertions/deletions. 
`cerebra` is not capable of reporting larger structural variants such as copy number variations and chromosomal rearrangements.

A data structure crucial to `cerebra` is the *genome interval tree*, which matches RNA transcripts
and polypeptides to each feature in the genome (\autoref{workflow}). 
[Interval trees](https://en.wikipedia.org/wiki/Interval_tree) are self-balancing binary search trees that store numeric intervals and can quickly retrieve every such interval that overlaps a given query interval (_[see also](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_). 
Given _n_ nodes, interval trees have theoretical average-case O(log*n*) and worst-case O(*n*) time complexity for search operations, making them tractable for genome-scale operations [@Cormen:2009, _[see also](https://www.coursera.org/lecture/algorithms-part1/interval-search-trees-ot9vw)_].
Tree construction proceeds at O(*n*log*n*) time complexity, making construction rather than search the bottleneck for most VCF sets [@Alekseyenko:2007]. 
The _genome interval tree_ is constructed with a reference genome sequence ([FASTA format](https://en.wikipedia.org/wiki/FASTA_format), often with a `.fa` extension), and a genome annotation 
([gene transfer format, GTF](https://www.gencodegenes.org/pages/data_format.html), `.gtf` extension).
We rely on the [ncls](https://github.com/biocore-ntnu/ncls) python library for fast interval tree construction and lookup operations [@Alekseyenko:2007].

In order to analyze multiple VCF records at once, we use [multiprocessing](https://en.wikipedia.org/wiki/Multiprocessing) with the Python [pathos](https://pypi.org/project/pathos/) library module [@pathos].
We extract relevant information -- including genomic interval, observed base, and read coverage -- from each variant record. 
In the `germline-filter` module variants are compared to one another and filtered out if found to be identical.
In `count-variants` variants are simply matched to whichever gene they came from. 
In `find-peptide-variants` variants are queried against our _genome interval tree_ -- if a matching interval is found we convert the DNA-level variant to a protein variant. 
Finally, protein variants across all VCF files are reported in tabular format. 

![Workflow describing the `find-peptide-variants` module.
We construct a genome interval tree from a genome annotation (.gtf) and a reference genome sequence (.fa), then process VCF files in parallel to create a single tabular output file (CSV or JSON).\label{workflow}](fig1.jpg)

## `germline-filter`

Variant calling is often applied to the study of cancer. 
If the research project is centered around a "tumor vs. normal" question, then `germline-filter` is the proper starting point. 
This module removes germline variants that are common between tumor and normal samples, and thus excludes variants unlikely to be pathogenic for the cancer under study.
The user provides a very simple metadata file (see [USAGE.md](https://github.com/czbiohub/cerebra/blob/master/docs/USAGE.md)) that indicates which tumor samples correspond to which normal samples.
Using the [vcfpy](https://pypi.org/project/vcfpy/) library we quickly identify shared variants across tumor/normal matched VCF files, then write new VCFs that contain only the unique variants [@vcfpy].
These steps are performed by a [subprocess pool](https://pypi.org/project/pathos/) so that we can process multiple discrete chunks of input at the same time. 

The output of `germline-filter` is a set of trimmed-down VCF files, which will be used for the next two steps. 
If you do not have access to "normal" samples then proceed directly to `count-variants` or `find-peptide-variants`. 

## `count-variants`

The `count-variants` module reports the raw variant counts for every gene across every sample.
We first create a _genome interval tree_ from the reference GTF, then read in a VCF file and convert it to a [vcfpy](https://pypi.org/project/vcfpy/) object, then processes VCF records in [parallel](https://en.wikipedia.org/wiki/Multiprocessing). 
Each variant is matched to its corresponding gene, and gene-wise counts are stored in [shared memory](https://en.wikipedia.org/wiki/Shared_memory).

If working with cancer samples, the user has the option to limit the reported variants to those also found in Wellcome Sanger Institute's [COSMIC](https://cancer.sanger.ac.uk/cosmic) database [@cosmic]. 
While certainly not exhaustive, this database contains an extensive list of known human variants. 
This option is designed to limit the search space to known and potentially actionable targets. 

`count-variants` produces two output files, one containing raw variant counts and one containing COSMIC filtered variant counts for every gene in the genome. 

## `find-peptide-variants`

The `find-peptide-variants` module reports the protein variants potentially associated with each genomic variant. 
First we load the reference GTF, then construct an index (.fai) of the genome fasta file with [pyfaidx](https://pypi.org/project/pyfaidx/) to enable fast random memory access [@pyfaidx].
We then create a _genome interval tree_ that will be used to quickly match genomic coordinates from VCF records to protein variants. 
The user again has the option to limit the search space to variants found in the COSMIC database. 
VCF records are read in simultaneously; individual records are converted to _GenomePosition_ objects to keep track of their genomic intervals and observed DNA bases.
_GenomePositions_ are then queried against the _genome interval tree_. 
If an overlapping interval is found, we retrieve the protein variant from this node of the _genome interval tree_. 
Protein variants are converted to [ENSEMBL](https://www.ensembl.org/index.html) protein IDs, in accordance with the [HGVS](https://varnomen.hgvs.org/) sequence variant nomenclature [@ensembl; @hgvs].
The output is a hierarchically ordered text file (CSV or JSON) that reports the the ENSEMBL protein ID and the gene associated with each variant, for each experimental sample.    

Variant callers are known to produce a great deal of false positives, especially when applied to single-cell RNA-seq data [@Enge:2017].
To address this concern, we have included the `coverage` option. 
If indicated this option will report counts for both variant and wildtype reads at all variant loci. 
We reasoned that variants with a high degree of read support are less likely to be false positives.
This option is designed to give the user more confidence in individual variant calls.        

We should emphasize that `find-peptide-variants` does not *definitively* report protein variants but rather the *likely* set of protein variants. 
Definitively reporting protein variants from RNA-seq requires knowledge of alternate splicing -- this represents an open problem in the field [@Huang:2017]. 
For example, if a read picks up a variant in exon 2 of a given gene, we can report each of the potential isoforms of that gene that contain exon 2, but we **cannot** infer which of those particular isoforms are actually present in our sample (see \autoref{splice}). 
For the example shown in \autoref{splice} we would translate and report _t1_ and _t3_ as both of these contain exon 2. 
It is possible the sample does not actually express both of these isoforms, however, determining the isoform landscape of a sample from RNA-seq is outside the scope of this project. 

![For a given mutational event, `cerebra` reports ALL potentially affected isoforms.\label{splice}](fig2.jpg)

To assess performance of `find-peptide-variants` we obtained VCFs from a single-cell RNA-seq study conducted on lung adenocarcinoma patient samples [@Maynard:2020]. 
These VCFs were produced with STAR (alignment) and GATK HaplotypeCaller (variant calling), and are on the order of megabytes, typical of a single-cell RNA-seq experiment. 
`cerebra` was run on standard hardware (MacBook Pro, 2.5GHz quad-core processor, 16 GB RAM).
As show in \autoref{runtime} `cerebra` processed a set of 100 VCF files in approximately 34 minutes. 

![`cerebra` processes 100 VCF files (~400 Mb in total) in ~34 minutes.\label{runtime}](fig3.jpg)

The first 10 or so minutes of `cerebra find-peptide-variants` do not involve any VCF processing, instead, this time is attributed to the _genome interval tree_ construction phase.
After the tree is built, files are processed in a near-linear manner. 
Also of note is that `cerebra`'s search operations take advantage of multiprocessing.
`cerebra` should scale better to high-memory machines than single-threaded tools, though it has been designed to run on standard hardware.

## Conclusions

RNA/DNA sequencing paired with fast and accurate summarizing of variants is often crucial to understanding the biology of an experimental system. 
We present a tool that can be used to quickly summarize the variant calls contained within a large set of VCF files.
As sequencing costs continue to drop, large-scale variant calling will become more accessible, and summary tools like `cerebra` will become essential for drawing meaningful conclusions in a reasonable time frame. 
Our tool offers the advantages of parallel processing and a single, easy-to-interpret output file (CSV or JSON).

`cerebra` is already enabling research, see [@Maynard:2020], a study that examines the tumor microenvironment of late-stage drug-resistant carcinomas with single-cell RNA-sequencing. 
Understanding the mutational landscape of individual tumors was essential to this study, and given the sheer volume of VCF records, would not have been possible without `cerebra`. 
We hope that `cerebra` can provide an easy-to-use framework for future studies in the same vein. 

## Acknowledgments

Funding for this work was provided by the [Chan Zuckerberg Biohub](https://www.czbiohub.org/).
The authors would like to thank Ashley Maynard, Angela Pisco and Daniel Le for helpful discussions and support.

## Correspondence

Please contact `ljharris018@gmail.com`

## Code

`cerebra` is written in Python 3. 
Code and detailed installation instructions can be found at https://github.com/czbiohub/cerebra. 
In addition, `cerebra` can be found on [PyPi](https://pypi.org/project/cerebra/) and [Dockerhub](https://hub.docker.com/r/lincolnharris/cerebra).

## References

History
=======

v1.1.3 (2020-09-15)
---------------------
some pkg requirement versions pinned -- should work better with PyPi now.   
various updates to READMEs and instruction files.


v1.1.2 (2020-07-31)
---------------------         
the majority of reviewer comments have now been addressed.

v1.1.1 (2020-07-20)
---------------------      
added docker support, updated requirement files, added multi-platform CI testing, added Docker continuous build feature, fixed refs in text, added adnl refs, updated CONTRIBUTING.md, updated MAKEFILE, updated ORCID ids, removed docs/ folder.


v1.1.0 (2020-07-13)
---------------------       
starting to address JOSS reviewer comments.     
A lot of work still to be done.


v1.0.0 (2020-01-31)
---------------------

_initial release_             
basic functionality achieved.





Credits
=======

Development Lead
----------------

-   Lincoln Harris <[@lincoln-harris](https://github.com/lincoln-harris)> &lt;ljharris018@gmail.com&gt;

Contributors
------------

-   Olga Botvinnik <[@olgabot](https://github.com/olgabot)>
-   Rohan Vanheusden <[@rvanheusden](https://github.com/rvanheusden)>
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

### Report Bugs

Report bugs at https://github.com/czbiohub/cerebra/issues.

If you are reporting a bug, please include:

-   Your operating system name and version.
-   The python version you are using. 
-   Any other details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.
-   Screenshots or code snippets, when applicable. 

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features/enhancement requests. Anything tagged with "enhancement" is open to whoever wants to implement it.

### Write Documentation

`cerebra` could always use more documentation, whether as
part of the official `cerebra` docs, in docstrings, or
even on the web in blog posts, articles, and such.

If there is a function/module that you think could use additional documentation, or is simply unclear to you, please let us know. 

### Submit Feedback

The best way to send feedback is to file an issue at https://github.com/czbiohub/cerebra/issues.

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :)

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1.  The pull request should include tests, ideally with >60% code coverage. 
2.  Functionality should be encapsulated within a function(s) that includes a docstring. 
3.  The pull request should work for Python 3.6, 3.7 and potentially 3.8. Check
    https://travis-ci.org/github/czbiohub/cerebra/pull_requests and make sure that the tests pass
    for all supported Python versions.

Tips
----

To run a subset of tests:

    $ pytest test_you_want_to_run.py
Usage (advanced)
============

## General Note

If working with extraordinarily large VCF files (>= 1Mb), the following processing steps are likely to be very slow. 
Before starting it might be a good idea to trim down your VCF files to only the region/chromosome of interest, using something like [bcftools](http://samtools.github.io/bcftools/bcftools-man.html#view). 

## Filtering germline variants

If you have access to germline or "normal" tissue, then `germline-filter` is a good place to start. 
All you have to do is provide a simple metadata file that links each tumor sample to its corresponding normal sample. 
For example: 
```
tumor_sample_id,normal_sample_id
sample1,gl_sample1
sample2,gl_sample1
sample3,gl_sample2
sample4,gl_sample2
sample5,gl_sample2
```
Once you have made this metadata file you're ready to run `germline-filter.` 
An example command line:   
```
cerebra germline-filter --processes 2 --normal_path /path/to/normal/vcfs --tumor_path /path/to/tumor/vcfs --metadata /path/to/metadata/file --outdir /path/to/filtered/vcfs
```

This will create a new directory (`/path/to/filtered/vcfs/`) that contains a set of entirely new VCFs. 

## Counting variants 

The module `count-variants` module can be run after `germline-filter`, on the new vcfs contained in `/path/to/filtered/vcfs/`.
However, `germline-filter` is entirely optional -- if you dont have access to germline or "normal" samples, `count-variants` is the place to start. 
An example command line:   
```
cerebra count-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --refgenome /path/to/genome/annotation --outfile /path/to/output/file /path/to/filtered/vcfs/*
```

NOTE that the cosmic database is also optional. If you'd like you can download one of the database files from [here](https://cancer.sanger.ac.uk/cosmic/download), however, you can also run `count-variants` without this option. 

## Finding peptide variants

Like `count-variants`, `find-peptide-variants` is a standalone module. 
You can run it on the VCFs generated by `germline-filter` or on unfiltered VCFs. 
Also like `count-variants`, this module gives you the option of filtering through a cosmic database. 
An example command line:   
```
cerebra find-peptide-variants --processes 2 --cosmicdb /optional/path/to/cosmic/database --annotation /path/to/genome/annotation --genomefa /path/to/genome/fasta --report_coverage 1 --output /path/to/output/file /path/to/filtered/vcfs/*
```

`report_coverage` is a BOOLEAN option will report counts for both variant and wildtype reads at all variant loci.
`--report_coverage 1` turns this option on, while `--report_coverage 0` turns it off.
We reasoned that variants with a high degree of read support are less likely to be false positives. 
This option is designed to give the user more confidence in individual variant calls.

For example, when run on the pre-packaged test VCF set (_cerebra/tests/data/test_find_peptide_variants/vcf_),    
`cerebra find-peptide-variants --report_coverage 1` should yield the following (partial) entry:

```
A1,['ENSP00000395243.3:p.(Leu813delinsArgTrp),[2:0]', 'ENSP00000415559.1:p.(Leu813delinsArgTrp),[2:0]'],
```
This tells us that the sample _A1_ contains likely variants in the Ensembl peptide IDs _ENSP00000395243.3_ and _ENSP00000415559.1_. 
Both variants are insertions of _ArgTrp_ in place of _Leu_ at the 813th amino acid.

The [_x_:_y_] string represents the absolute number of variant and wildtype reads at that loci. 
Thus [2:0] means 2 variant reads and 0 wildtype reads were found at each of these loci. 
A coverage string in the format of [_x_:_y_:_z_] would indicate there are two variant alleles at a given loci, _x_ and _y_, in addition to wildtype, _z_. 

## Testing

First install the packages specified in [test_requirements.txt](https://github.com/czbiohub/cerebra/blob/master/test_requirements.txt). 
Now you should be able to run:

` $ make test `

If you've installed `cerebra` in a virtual environment make sure the environment is active. 
Confirm that all tests have passed.
If otherwise, feel free to submit an [issue report](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md). 
Installation Instructions
============

Dependencies
------------
`cerebra` depends on some (fairly standard) packages and libraries. 
Before installing it might be a good idea to make sure all of the requisite packages are installed on your system (_note:_ if installing with Docker you can skip this step). 

__MacOS Dependencies:__
```
sudo pip install setuptools
brew update
brew install openssl
brew install zlib
```

__Linux Dependencies:__
```
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

As of present `cerebra` is not installable on Windows. 
`cerebra` depends on the [`pysam`](https://pysam.readthedocs.io/en/latest/index.html) library (or rather, `pysam` is a dependency-of-a-dependency) and currently this library is only available on Unix-like systems. 
Windows solutions like [WSL](https://docs.microsoft.com/en-us/windows/wsl/) exist for overcoming precisely this challange, however, `cerebra` has not been tested on WSL or any other Unix-like subsystem for Windows.    


Installation (for users)
------------
There are four different methods available to install `cerebra`.
Choose one of the following:

__With [Docker](https://hub.docker.com/r/lincolnharris/cerebra) (recommended)__          
```
docker pull lincolnharris/cerebra
docker run -it lincolnharris/cerebra
```      
_warning: this image will take up ~1Gb of space._               

__With traditional git clone and the python standard library [`venv`](https://docs.python.org/3/library/venv.html) module__
```
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
python3 -m venv cerebra-dev
source cerebra-dev/bin/activate
pip3 install [--user] . 
```

__With traditional git clone and [conda](https://docs.conda.io/en/latest/)__
``` 
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install [--user] . 
```

__From [PyPi](https://pypi.org/project/cerebra/)(system-wide installation, NOT RECOMMENDED)__           
For novice users, its generally a better idea to install packages within virtual environments. 
However, `cerebra` can be installed system-wide, if you so choose. 
```
pip install cerebra

# OR, if you dont have root privileges
pip install --user cerebra
```


Installation (for developers)
------------
Here's how to set up cerebra for local development. 
After installing the requisite [dependencies](https://github.com/czbiohub/cerebra/blob/master/docs/INSTALL.md#dependencies):

1.  Fork the `cerebra` repo on GitHub: https://github.com/czbiohub/cerebra
2.  Clone your fork locally:

        $ git clone https://github.com/your-name/cerebra.git

3.  Install your local copy into a virtualenv. Using the standard library [`venv`](https://docs.python.org/3/library/venv.html) module: 

        $ cd cerebra
        $ python3 -m venv cerebra-dev
        $ source cerebra-dev/bin/activate
        $ pip install -r requirements.txt -r test_requirements.txt -e . 

4.  Create a branch for local development:

        $ git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5.  When you're done making changes, check that your changes pass flake8 and the tests:

        $ make test
        $ make coverage
        $ make lint

6.  Commit your changes and push your branch to GitHub:

        $ git add .
        $ git commit -m "Your detailed description of your changes."
        $ git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.
See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for more. 
