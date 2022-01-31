# WHdenovo
A cost-effective approach to diploid assembly for single samples and trios. It includes the following steps: construct sequence graph from illumina, align long reads to the graph and partition these long reads to two haplotypes.

### Installation

Conda installation
```
conda install -c shilpagarg whdenovo
```

### Data Simulation

For simulating Illumina data:
```
whdenovo simulate illumina N-bp.fasta <het> out/illumina/
```
And it will show you the file you'll need for pacbio data simulation and WHdenovo test

For simulating PacBio data:
```
whdenovo simulate pacbio sample.fastq <coverage> out/pacbio/ \
                       <mom_hap1.fasta> <mom_hap2.fasta> \
                       <dad_hap1.fasta> <dad_hap2.fasta> \
                       <child_hap1.fasta> <child_hap2.fasta>
```
### Assembly

Trio case:
```
whdenovo partition --illumina1 <illumina_child_1.fq> --illumina2 <illumina_child_2.fq> \
                       --pacbio <pacbio_mom.fasta> <pacbio_dad.fasta> <pacbio_child.fasta>
                       -p ped [-t <thread>] [-o out/path] [--lowc INT] [--highc INT]
```

e.g.

```
whdenovo partition --illumina1 16513_illumina_10k_1.5/child.het1.5.cov30_1.fq --illumina2 16513_illumina_10k_1.5/child.het1.5.cov30_2.fq \
                       --pacbio 16513_pacbio_10k_1.5_20/pacbio_mom.fasta 16513_pacbio_10k_1.5_20/pacbio_dad.fasta 16513_pacbio_10k_1.5_20/pacbio_child.fasta \
                       -p ped -t 24 -o test.simu --lowc 5 --highc 60
```

Individual case
```
whdenovo partition --illumina1 <illumina_who_1.fq> --illumina2 <illumina_who_2.fq> \ 
                       --pacbio <pacbio_who.fasta> [-t <thread>] [-o out/path]
```
For assembling the genome from partitioned reads given by ```partition.py```:

```
conda activate flye
whdenovo assemble -f son.inputreads.fa \
                      -0 path/to/output/HP0.reads -1 path/to/output/HP1.reads \
                      --assemble -s 15k -t 40 <--pacbio|--nano>
```
If you wish to use other assemblers with partitioned reads, just remove the ```--assemble``` argument and you may not need to activate the virtual environment.

### Result Validation

For validating the partitioning of simulated the data.
```
whdenovo validate -p out/path -f <pacbio_who.fasta>
```
For validating the partitioning of real data when you have ground truth classification:
```
whdenovo validate -p out/path -t <tagged.reads.txt>
```
An example for tagged.reads.txt is at test/haplotagged.reads.txt, which should include the HP tag and PS tag from haplotagged BAM file.
***
We acknowledge the support of dependencies such as [bfc](https://github.com/lh3/bfc), [SPAdes](http://cab.spbu.ru/software/spades/), [vg](https://github.com/vgteam/vg) and [GraphAligner](https://github.com/maickrau/GraphAligner).

### Citations
1. A graph-based approach to diploid genome assembly. [Link](https://doi.org/10.1093/bioinformatics/bty279)
2. A haplotype-aware de novo assembly of related individuals using pedigree sequence graph. [Link](https://doi.org/10.1093/bioinformatics/btz942)
[![PyPI](https://img.shields.io/pypi/v/whatshap.svg)](https://pypi.python.org/pypi/whatshap)
[![Build Status](https://semaphoreci.com/api/v1/whatshap/whatshap/branches/master/shields_badge.svg)](https://semaphoreci.com/whatshap/whatshap)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/whatshap/README.html)

![WhatsHap logo](https://bitbucket.org/repo/8AjxBd/images/543323972-whatshap_logo.png)

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.

For documentation and information on how to cite WhatsHap, please visit the [WhatsHap Homepage](https://whatshap.readthedocs.io/)
whatshap_logo_text.svg: Logo with editable text (Syntax font needs to be installed)
whatshap_logo.svg: Logo with text converted to paths
whatshap_logo_small.svg: Top part of logo only, without the text
=======
Changes
=======

v0.13 (2016-10-27)
------------------

* Use ``PS`` tag instead of ``HP`` tag by default to store phasing information.
  This applies to the ``phase`` and ``hapcut2vcf`` subcommands. ``PS`` is also
  used by other tools and standard according to the VCF specification.
* Incorporated genotype likelihoods into our phasing framework. On request
  (by using option ``--distrust-genotypes``), genotypes can now be changed at a cost
  corresponding to their input genotype likelihoods. The changed genotypes are
  written to the output VCF. The behavior of ``--distrust-genotypes`` can be
  fine-tuned by the added options ``--include-homozygous``, ``--default-gq``,
  ``--gl-regularizer``, and ``--changed-genotype-list``.
* Correctly handle cases when processing VCFs with two or more disjoint
  families.

v0.12 (2016-07-01)
------------------

* Speed up allele detection
* Add an ``unphase`` subcommand which removes all phasing from a VCF file
  (``HP`` and ``PS`` tags, pipe notation).
* Add option ``--tag=`` to the ``phase`` subcommand, which allows to choose
  whether ReadBackedPhasing-compatible ``HP`` tags or standard ``PS`` tags are
  used to describe phasing in the output VCF.
* Manage versions with `versioneer <https://github.com/warner/python-versioneer>`_.
  This means that ``whatshap --version`` and the program version in the VCF header
  will include the Git commit hash, such as ``whatshap 0.11+50.g1b7af7a``.
* Add subcommand "haplotag" to tag reads in a BAM file with their haplotype.
* Fix a bug where re-alignment around variants at the very end of a chromosome
  would lead to an AssertionError.

v0.11 (2016-06-09)
------------------

* When phasing a pedigree, blocks that are not connected by reads but
  can be phased based on genotypes will be connected per default. This
  behavior can be turned off using option ``--no-genetic-haplotyping``.
* Implemented allele detection through re-alignment: To detect which allele of a
  variant is seen in a read, the query is aligned to the two haplotypes at that
  position. This results in better quality phasing, especially for
  low-quality reads (PacBio). Enabled if ``--reference`` is provided. Current
  limitation: No score for the allele is computed.
* As a side-effect of the new allele detection, we can now also phase
  insertions, deletions, MNPs and "complex" variants.
* Added option ``--chromosome`` to only work on specifed chromosomes.
* Use constant recombination rate per default, allows to use ``--ped``
  without using ``--genmap``.
* ``whatshap`` has become a command with subcommands. From now on, you need
  to run ``whatshap phase`` to phase VCFs.
* Add a ``stats`` subcommand that prints statistics about phased VCFs.

v0.10 (2016-04-27)
------------------

* Use ``--ped`` to phase pedigrees with the PedMEC algorithm
* Phase all samples in a multi-sample VCF
* Drop support for Python 3.2 - we require at least Python 3.3 now

v0.9 (2016-01-05)
-----------------

* This is the first release available via PyPI (and that can therefore be
  installed via ``pip install whatshap``)

January 2016
------------

* Trio phasing implemented in a branch

September 2015
--------------

* pWhatsHap implemented (in a branch)

April 2015
----------

* Create haplotype-specific BAM files

February 2015
-------------

* Smart read selection

January 2015
------------

* Ability to read multiple BAM files and merge them on the fly

December 2014
-------------

* Logo
* Unit tests

November 2014
-------------

* Cython wrapper for C++ code done
* Ability to write a phased VCF (using HP tags).

June 2014
---------

* Repository for WhatsHap refactoring created

April 2014
----------

* The WhatsHap algorithm is introduced at RECOMB
.. image:: https://img.shields.io/pypi/v/whatshap.svg?branch=master
    :target: https://pypi.python.org/pypi/whatshap

.. image:: https://semaphoreci.com/api/v1/whatshap/whatshap/branches/master/shields_badge.svg
    :target: https://semaphoreci.com/whatshap/whatshap

|

.. image:: https://bitbucket.org/repo/8AjxBd/images/3378940113-whatshap_logo.png
    :scale: 50%

|

WhatsHap
========

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *haplotype assembly*. It is especially suitable for long
reads, but works also well with short reads.

If you use WhatsHap, please cite:

    Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel,
    Leen Stougie, Gunnar W. Klau, Alexander Schönhuth.
    `WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads
    <http://dx.doi.org/10.1089/cmb.2014.0157>`_.
    Journal of Computational Biology, 22(6), pp. 498-509, 2015.
    (`Get a self-archived version <https://bioinf.mpi-inf.mpg.de/homepage/publications.php?&account=marschal>`_)

The version of WhatsHap you find here is the result of further development
focused on making the software easy and straightforward to use. WhatsHap is now
Open Source software under the MIT license and we welcome contributions.


Pedigree phasing
----------------

WhatsHap is capable of using pedigree information
about individuals to further improve phasing results, and to drastically reduce
the needed coverage. A preprint is available on bioRxiv:

    Read-Based Phasing of Related Individuals.
    Shilpa Garg, Marcel Martin, Tobias Marschall.
    `doi:10.1101/037101 <http://dx.doi.org/10.1101/037101>`_


Parallel Version: pWhatsHap
---------------------------
A parallelization of the core dynamic programming algorithm has been described in 

    M. Aldinucci, A. Bracciali, T. Marschall, M. Patterson, N. Pisanti, M. Torquati. 
    `High-Performance Haplotype Assembly <http://dx.doi.org/10.1007/978-3-319-24462-4_21>`_. Proceedings of the 11th International
    Meeting on Computational Intelligence Methods for Bioinformatics and
    Biostatistics (CIBB), 245-258, 2015.

The current implementation can be found in `branch parallel <https://bitbucket.org/whatshap/whatshap/branch/parallel>`_.


Documentation
-------------

* `Bitbucket page <https://bitbucket.org/whatshap/whatshap/>`_
* `Read the documentation online <https://whatshap.readthedocs.io/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the
  repository and in the downloaded tar distribution.


Mailing list
------------
We run a `public mailing list <https://lists.cwi.nl/mailman/listinfo/whatshap>`_. Please
don't hesitate to post questions and comments.
Various notes
=============

* There is a step in which variants are re-discovered in the BAM file. This may
  fail when the variant caller has used some type of re-alignment (as
  freebayes does). Would be better to integrate this into the variant caller or
  to get the information out of it. This applies only to indels, which are not
  supported right now anyway.
* Input format for HapCompass: http://www.brown.edu/Research/Istrail_Lab/resources/hapcompass_manual.html#sec11


Allele detection with re-alignment
----------------------------------

WhatsHap can detect which allele a read contains at a variant position by
aligining a part of the read to the two possible haplotypes. The haplotype
for which the alignment is better wins.

Allele detection through re-alignment is enabled when the ``--reference``
parameter is used on the command-line.

Re-alignment in this version detects slightly *fewer* alleles than the old
algorithm, but this is typically justified because the old algorithm gave
wrong results. Re-alignment however correctly detects that both haplotypes are
equally good and then refuses to choose.

The alignment algorithm uses edit distance at the moment, which allows us to
detect alleles correctly most of the time, but does not allow us to make use
of base qualities (in fact, the weighted algorithm degenerates into an
unweighted one). To fix this, we need a better alignment algorithm.

Here are some examples for how re-alignment works.

Insertion next to a SNP
~~~~~~~~~~~~~~~~~~~~~~~

Haplotypes::

    ref:   CCTTAGT
    alt:   CCTCAGT

Alignment as reported in BAM file::

    ref:   CCT-TAGT
    query: CCTCAAGT

The second ``T`` is aligned to an ``A``, which is not one of the expected bases.
Thus, no variant would be detected here.

Re-aligning the query to the "alt" haplotype, we get::

    alt:   CCTCA-GT
    query: CCTCAAGT

This alignment has lower cost and we therefore detect that the allele in this
read is probably the alternative one.


Ambiguous
~~~~~~~~~

This was previously detected incorrectly::

    ref:   TGCTTTAAGG
    alt:   TGCTTTCAGG
    query: TGCCTTCAAGG

Two possible alignments are ::

    ref:   TGC-TTTAAGG
    query: TGCCTTCAAGG

and ::

    alt:   TGCTTTCA-GG
    query: TGCCTTCAAGG

Both have cost two and therefore the correct allele is unclear.
.. include:: ../CHANGES.rst
Developing
==========

This is the developer documentation for WhatsHap.


Development installation
------------------------

For development, make sure that you install Cython and tox. We also recommend
using a virtualenv. This sequence of commands should work::

	git clone https://bitbucket.org/whatshap/whatshap
	cd whatshap
	virtualenv -p python3 venv
	venv/bin/pip3 install Cython nose tox
	venv/bin/pip3 install -e .

Then you can run WhatsHap like this (or activate the virtualenv and omit the
``venv/bin`` part)::

	venv/bin/whatshap --help

The tests can be run like this::

	venv/bin/tox


Development installation (alternative)
--------------------------------------

Alternatively, if you do not want to use virtualenv, you can do the following::

	git clone https://bitbucket.org/whatshap/whatshap.git
	cd whatshap
	python3 setup.py build_ext -i
	bin/whatshap

This requires Cython, pysam, and pyvcf to be installed.


Installing other Python versions in Ubuntu
------------------------------------------

Ubuntu comes with one default Python 3 version, and in order to test WhatsHap
with other Python versions (3.2, 3.3 and 3.4), use the “deadsnakes” repository.
Ensure you have the following packages::

	sudo apt-get install build-essential python-software-properties

Then get and install the desired Python versions. For example, for Python 3.2::

	sudo add-apt-repository ppa:fkrull/deadsnakes
	sudo apt-get update
	sudo apt-get install python3.2-dev python3-setuptools

If pip and virtualenv are not available, install them (Since they are so essential,
we use sudo to install them system-wide, but you can also install them into
your $HOME by omitting the sudo and adding the ``--user`` option instead)::

	sudo easy_install3 pip
	sudo pip3 install virtualenv


Debugging
---------

Here is one way to get a backtrace from gdb (assuming the bug occurs while
running the tests)::

	$ gdb python3
	(gdb) run -m nose

After you get a SIGSEGV, let gdb print a backtrace:

	(gdb) bt


Making a release
----------------

If this is the first time you attempt to upload a distribution to PyPI, create a
configuration file named ``.pypirc`` in your home directory with the following
contents::

	[distutils]
	index-servers =
	    pypi

	[pypi]
	username=my-user-name
	password=my-password

See also `this blog post about getting started with
PyPI <http://peterdowns.com/posts/first-time-with-pypi.html>`_. In particular,
note that a ``%`` in your password needs to be doubled and that the password
must *not* be put between quotation marks even if it contains spaces.

#. Set the correct version number in the changelog. Ensure that the list of changes is up-to-date.

#. Ensure you have no uncommitted changes in the working copy.

#. Run ``tox``, ensuring all tests pass.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag v0.1

#. Create a distribution (``.tar.gz`` file), ensuring that the auto-generated version number in
   the tarball is as you expect it::

       python3 setup.py sdist

#. Upload the distribution to PyPI (the tarball must be regenerated since ``upload`` requires a preceding ``sdist``)::

       python3 setup.py sdist upload

#. Push the tag::

       git push --tags

#. Update the `bioconda recipe <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/whatshap/meta.yaml>`_.
   It is probly easiest to edit the recipe via the web interface and send in a
   pull request. Ensure that the list of dependencies (the ``requirements:``
   section in the recipe) is in sync with the ``setup.py`` file.

   Since this is just a version bump, the pull request does not need a
   review by other bioconda developers. As soon as the tests pass and if you
   have the proper permissions, it can be merged directly.

If something went wrong, fix the problem and follow the above instructions again,
but with an incremented revision in the version number. That is, go from version
x.y to x.y.1. Do not change a version that has already been uploaded.


Adding a new subcommand
-----------------------

Follow the instructions in ``whatshap/example.py``.
==========
User guide
==========

Run WhatsHap like this::

    whatshap phase -o phased.vcf input.vcf input.bam

Phasing information is added to the VCF file in a way that is compatible with
GATK’s ReadBackedPhasing. That is, the HP tag denotes which set of phased
variants a variant belongs to. The VCF file can also be gzip-compressed.


Features and limitations
========================

WhatsHap supports phasing of variants in diploid genomes.

Supported variant types are SNVs (single-nucleotide variants), insertions,
deletions, MNPs (multiple adjacent SNVs) and “complex” variants. Complex
variants are those that do not fall in any of the other categories, but
are not structural variants. An example is the variant TGCA → AAC.
Structural variants are not phased.

If no reference sequence is provided (using ``--reference``), only
SNVs, insertions and deletions can be phased.

All variants in the input VCF that are marked as being heterozygous
(genotype 0/1) and that have appropriate coverage are used as input for the core
phasing algorithm. If the algorithm could determine how the variant should be
phased, that information will be added to the variant in the output VCF.

Variants can be left unphased for two reasons: Either the variant type is
not supported or the phasing algorithm could not make a phasing decision.
In both cases, the information from the input VCF is simply copied to output
VCF unchanged.


Recommended workflow
====================

Best phasing results are obtained if you sequence your sample(s) on both PacBio
and Illumina: Illumina for high-quality variant calls and PacBio for its long
reads.

1. Map your reads to the reference, making sure that you assign each read to a
read group (the "@RG" header line in the BAM file). WhatsHap supports VCF files
with multiple samples and in order to determine which reads belong to which
sample, it uses the 'sample name' (SM) of the read group. If you have a single
sample only and no or incorrect read group headers, you can run WhatsHap with
``--ignore-read-groups`` instead.

2. Call variants in your sample(s) using the most accurate reads you have. These
will typically be Illumina reads, resulting in a a set of variant calls you can
be reasonably confident in. If you do not know which variant caller to use, yet,
we recommend FreeBayes, which is fast, Open Source and easy to use. In any case,
you will need a standard VCF file as input for WhatsHap in the next step.

3. Run WhatsHap with the VCF file of high-confidence variant calls (obtained in
the previous step) and with the *longest* reads you have. These will typically
be PacBio reads. Phasing works best with long reads, but WhatsHap can use any
read that covers at least two heterozygous variant calls, so even paired-end or
mate-pair reads are somewhat helpful. If you have multiple sets of reads, you
can combine them by providing multiple BAM files on the command line.


Representation of phasing information in VCFs
=============================================

WhatsHap supports two ways in which it can store phasing information in a VCF
file: The standards-compliant ``PS`` tag and the ``HP`` tag used by GATK’s
ReadBackedPhasing tool. When you run ``whatshap phase``, you can select which
format is used by setting ``--tag=PS`` or ``--tag=HP``.

We will use a small VCF file as an example in the following. Unphased, it
looks like this::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0/1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      0/1      1/1
    chr1    300  .   G    T    50.0   .       .    GT      0/1      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0/1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0/1      1/1

Note that sample 1 is heterozygous at all shown loci (expressed with
``0/1`` in the ``GT`` field).


Phasing represented by pipe (``|``) notation
--------------------------------------------

The ``GT`` fields can be phased by ordering the alleles by haplotype and
separating them with a pipe symbol (``|``) instead of a slash (``/``)::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0|1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      1|0      0/1
    chr1    300  .   G    T    50.0   .       .    GT      1|0      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0|1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0|1      1/1

The alleles on one of the haplotypes of sample1 are: A, G, T, T, A.
On the other haplotype, they are: T, C, G, A, G.

Swapping ones and zeros in the ``GT`` fields would result in a VCF file with
the equivalent information.


Phasing represented by PS ("phase set") tag
-------------------------------------------

The pipe notation has problems when not all variants in the VCF file can be
phased. The `VCF specification <https://github.com/samtools/hts-specs>`_
introduces the ``PS`` tag to solve some of them. The ``PS`` is a
unique identifier for a "phase set", which is a set of variants that were
be phased relative to each other. There are usually multiple phase sets in
the file, and variants that belong to the same phase set do not need to
be consecutive in the file::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1      sample2
    chr1    100  .   A    T    50.0   .       .    GT:PS:PQ   0|1:100:22   0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:PS:PQ   1|0:100:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:PS:PQ   1|0:300:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:PS:PQ   0|1:300:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:PS:PQ   0|1:100:12   0/1:.:.

This VCF contains two phase sets named ``100`` and ``300``. The names are
arbitrary, but WhatsHap will choose the position of the leftmost variant
of the phase set as its name. The variants at 100, 150 and 500 are in the same
phase set, while the variants at 300 and 350 are in a different phase set.
Such a configuration is typically seen when paired-end or mate-pair reads are
used for phasing.

In the case of WhatsHap, the phase sets are identical to the connected components
of the variant connectivity graph. Two variants in that graph are connected if a
read exists that covers them.

The above example also shows usage of the ``PQ`` tag for "phasing quality".
WhatsHap currently does not add this tag.


Phasing represented by HP tag
-----------------------------

GATK’s ReadBackedPhasing tool uses a different way to represent phased variants.
It is in principle the same as the combination of pipe notation with the ``PS``
tag, but the ``GT`` field is left unchanged and all information is added to a
separate ``HP`` tag ("haplotype identifier") instead. This file encodes the same
information as the example above::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1         sample2
    chr1    100  .   A    T    50.0   .       .    GT:HP      0/1:100-1,100-2      0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:HP:PQ   0/1:100-2,100-1:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:HP:PQ   0/1:300-2,300-1:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:HP:PQ   0/1:300-1,300-2:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:HP:PQ   0/1:100-1,100-2:12   0/1:.:.

A few notes:

* ReadBackedPhasing does not add the ``PQ`` to the first variant in a phase set/haplotype
  group. This probably means that the phasing quality is to be interpreted as relative to
  the previous or first variant in the set.
* ReadBackedPhasing does not phase indels
* Discussions on the GATK forum on this topic:
   - https://gatkforums.broadinstitute.org/discussion/4226
   - https://gatkforums.broadinstitute.org/discussion/4038/


Trusting the variant caller
===========================

WhatsHap will trust the variant caller to have made the right decision of
whether a variant is heterozygous or homozygous. If you use the option
``--distrust-genotypes``, then this assumption is softened: An optimal solution
could involve switching a variant from being heterozygous to homozygous.
Currently, if that option is enabled and such a switch occurs, the variant
will simply appear as being unphased. No change of the genotype in the VCF is
done.

If you use this option, fewer variants will be phased.

Note that switching homozygous variants to heterozygous is never possible since
only heterozygous variants are considered for phasing.


.. _phasing-pedigrees:

Phasing pedigrees
=================

WhatsHap can take advantage of pedigree information to obtain a much
better phasing. To turn on pedigree mode, run WhatsHap like this::

	whatshap phase --ped pedigree.ped -o phased.vcf input.vcf input.bam

where ``pedigree.ped`` is a plink-compatible PED file to describe the
relationships between samples and ``input.vcf`` is a multi-sample VCF
with all individuals that should be phased. The reads for all individuals
can be in one or more BAM files. WhatsHap will match them based on sample
names provided in the read groups (just like for the default single-individual
mode).

PED file format
---------------

WhatsHap recognizes `PLINK-compatible PED
files <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>`_.
A PED file is a white-space (space or tab) delimited file with at least six
columns. WhatsHap checks the column count, but uses only

  * column 2: individual ID
  * column 3: paternal ID
  * column 4: maternal ID

The other columns are ignored. Lines starting with ``#`` are considered
comments and are ignored. Empty lines are also ignored.

To define a single trio, it is sufficient to have a single row in the PED file
with the child, mother and father. It is *not* necessary to include "dummy" rows
for individuals whose parents are unknown. (You will currently get a warning if
you do, but this will be changed.)

Here is an example defining a trio::

    # Fields: family, individual_id, paternal_id, maternal_id, sex, phenotype
    FAMILY01 the_child father mother 0 1

A quartet (note how multiple consecutive spaces are fine)::

    # Fields: family, individual_id, paternal_id, maternal_id, sex, phenotype
    FAMILY01 one_child   father mother 0 1
    FAMILY01 other_child father mother 0 1

*Important*: The names in the PED file *must* match the sample names in your VCF
and BAM files!

Pedigree phasing parameters
---------------------------

Phasing in pedigree mode requires costs for recombination events. Per
default, WhatsHap will assume a constant recombination rate across the
chromosome to be phased. The recombination rate (in cM/Mb) can be
changed by providing option ``--recombrate``. The default value of
1.26 cM/Mb is suitable for human genomes.

In order to use region-specific recombination rates, a genetic map file
can be provided via option ``--genmap``. WhatsHap expects a three-column
text file like this::

	position COMBINED_rate(cM/Mb) Genetic_Map(cM)
	55550 0 0
	568322 0 0
	568527 0 0
	721290 2.685807669 0.410292036939447
	723819 2.8222713027 0.417429561063975
	723891 2.9813105581 0.417644215424158
	...

The first (header) line is ignored and the three columns are expected to
give the pysical position (in bp), the local recombination rate between the
given position and the position given in the previous row (in cM/Mb), and
the cumulative genetic distance from the start of the chromosome (in cM).
The above example was taken from the 1000 Genomes genetic map `provided by
SHAPEIT
<https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap>`_.
Since genetic map files provide information for only one chromosome, the
``--genmap`` option has to be combined with ``--chromosome``.
============
Installation
============


Requirements
------------

WhatsHap is implemented in C++ and Python. You need to have a C++ compiler,
Python 3.3 (or later) and the corresponding Python header files. Python 3.5
is slightly faster than 3.4. In Ubuntu, make sure the packages
``build-essential`` and ``python3-dev`` are installed.


Installation
------------

WhatsHap can be installed with pip::

	pip3 install --user whatshap

This installs WhatsHap into ``$HOME/.local/bin``.  Then add
``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

	export PATH=$HOME/.local/bin:$PATH
	whatshap --help


Installation with Conda
-----------------------

WhatsHap is also available as a conda package from the `bioconda
channel <https://bioconda.github.io/>`_::

    conda install -c bioconda whatshap
=====================
Questions and Answers
=====================

 * **Can WhatsHap use a reference panel?** Reference panels are used by population-based phasers (like Beagle or ShapeIt). Although we are considering integrating this, WhatsHap cannot take advantage of reference panels right now. In case you have population data, we suggest to produce a population-based phasing (using `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_, `ShapeIt <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_, etc.) and a read-based phasing using WhatsHap separately and then compare/integrate results in a postprocessing step.
 * **Will Illumina data lead to a good read-based phasing?** Illumina paired-end data is not ideal for read-based phasing, since most pairs of heterozygous SNPs will not be bridged by a read pair. However, WhatsHap will attempt to produce as long haplotype blocks as possible. Running WhatsHap will hence tell you how phase-informative your input data is. Just take a look at the number (and sizes) of produced haplotype blocks.
 * **How large can/should a pedigree be for pedigree-aware read-based phasing (i.e. using option** :code:`--ped` **)?** The pedigree mode in WhatsHap is intended for intermediate-size pedigrees. The runtime of the core phasing step will be linear in :math:`2^{2t}`, where :math:`t` is the number of trio relationships (= number of children) in your pedigree. We do not recommend to use pedigrees with :math:`t>5`. For such pedigrees, read-data is unnecessary in most cases anyway and a very high-quality phasing can be obtained by genetic haplotyping methods (like `MERLIN <http://csg.sph.umich.edu/abecasis/Merlin/tour/haplotyping.html>`_).
.. include:: ../README.rst

=================
Table of contents
=================

.. toctree::
   :maxdepth: 2

   installation
   guide
   faq
   develop
   notes
   changes


..
   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
Parts of WhatsHap have been described in different articles. Please choose
an appropriate citation depending on your use case.

If you use WhatsHap as a tool:

    | Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O. Fischer,
      Nadia Pisanti, Gunnar W. Klau, Alexander Schoenhuth, Tobias Marschall.
    | *WhatsHap: fast and accurate read-based phasing*
    | bioRxiv 085050
    | doi: `10.1101/085050 <https://doi.org/10.1101/085050>`_

To refer to the core WhatsHap phasing algorithm:

    | Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel,
      Leen Stougie, Gunnar W. Klau, Alexander Schönhuth.
    | *WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads*
    | Journal of Computational Biology, 22(6), pp. 498-509, 2015.
    | doi: `10.1089/cmb.2014.0157 <http://dx.doi.org/10.1089/cmb.2014.0157>`_
      (`Get self-archived PDF <https://bioinf.mpi-inf.mpg.de/homepage/publications.php?&account=marschal>`_)

To refer to the pedigree-phasing algorithm and the PedMEC problem:

    | Shilpa Garg, Marcel Martin, Tobias Marschall.
    | *Read-based phasing of related individuals*
    | Bioinformatics 2016; 32 (12): i234-i242.
    | doi: `10.1093/bioinformatics/btw276 <https://doi.org/10.1093/bioinformatics/btw276>`_

WhatsHap's genotyping algorithm is described here:

    | Jana Ebler, Marina Haukness, Trevor Pesout, Tobias Marschall, Benedict Paten.
    | Haplotype-aware genotyping from noisy long reads
    | bioRxiv
    | doi: `10.1101/293944 <https://doi.org/10.1101/293944>`_

The HapChat algorithm is an alternative MEC solver able to handle higher coverages. It can be used
through "whatshap phase --algorithm=hapchat". It has been described in this paper:

    | Stefano Beretta, Murray Patterson, Simone Zaccaria, Gianluca Della Vedova, Paola Bonizzoni.
    | *HapCHAT: adaptive haplotype assembly for efficiently leveraging high coverage in long reads*.
    | BMC Bioinformatics, 19:252, 2018.
    | doi: `10.1186/s12859-018-2253-8 <https://doi.org/10.1186/s12859-018-2253-8>`_
    
A parallelization of the core dynamic programming algorithm (“pWhatsHap”)
has been described in

    | M. Aldinucci, A. Bracciali, T. Marschall, M. Patterson, N. Pisanti, M. Torquati.
    | *High-Performance Haplotype Assembly*
    | Proceedings of the 11th International Meeting on Computational Intelligence
      Methods for Bioinformatics and Biostatistics (CIBB), 245-258, 2015.
    | doi: `10.1007/978-3-319-24462-4_21 <http://dx.doi.org/10.1007/978-3-319-24462-4_21>`_

pWhatsHap is currently not integrated into the main WhatsHap source code. It
is available in
`branch parallel <https://bitbucket.org/whatshap/whatshap/branch/parallel>`_
in the Git repository.
=======
Changes
=======

master branch
-------------
* Integration of the HapChat algorithm as an alternative MEC solver, available
  through ``whatshap phase --algorithm=hapchat``. Contributed by the HapChat
  team, see https://doi.org/10.1186/s12859-018-2253-8.

v0.17 (2018-07-20)
------------------
* :issue:`140`: Haplotagging now works when chromosomes are missing in the VCF.
* Added option ``--merge-reads``, which is helpful for high coverage data.
* When phasing pedigrees, ensure that haplotypes are ordered as
  paternal_allele|maternal_allele in the output VCF. This seems to be a common
  convention and also used by 1000G.
* Test cases now use pytest instead of nose (which is discontinued).

v0.16 (2018-05-22)
------------------

* :issue:`167`: Fix the ``haplotag`` command. It would tag reads incorrectly.
* :issue:`154`: Use barcode information in BX tags when running ``haplotag``
  on 10x Genomics linked read data.
* :issue:`153`: Allow combination of ``--ped`` and ``--samples`` to only work
  on a subset of samples in a pedigree. Added ``--use-ped-samples`` to only
  phase samples mentioned in PED file (while ignoring other samples in input VCF).

v0.15 (2018-04-07)
------------------

* New subcommand ``genotype`` for haplotype-aware genotyping 
  (see https://doi.org/10.1101/293944 for details on the method).
* Support CRAM files in addition to BAM.
* :issue:`133`:
  No longer create BAM/CRAM index if it does not exist. This is safer when running multiple
  WhatsHap instances in parallel. From now on, you need to create the index yourself
  (for example with ``samtools index``) before running WhatsHap.
* :issue:`152`: Reads marked as “duplicate” in the input BAM/CRAM file are now ignored.
* :issue:`157`: Adapt to changed interface in Pysam 0.14.
* :issue:`158`: Handle read groups with missing sample (SM) tag correctly.

v0.14.1 (2017-07-07)
--------------------

* Fix compilation problem by distinguishing gcc and clang.

v0.14 (2017-07-06)
------------------

* Added ``--full-genotyping`` to (re-)genotype the given variants based on the reads
* Added option ``whatshap compare --switch-error-bed`` to write BED file with switch
  error positions
* Added ``whatshap compare --plot-blocksizes`` to plot histogroms of block sizes
* Added option ``--longest-block-tsv`` to output position-wise stats on longest joint
  haplotype block
* Added option ``whatshap compare --tsv-multiway`` to write results of multi-way
  comparison to tab-separated file
* Added option --chromosome to whatshap stats
* ``whatshap compare`` can now compute the block-wise Hamming distance
* ``whatshap stats`` can now compute an N50 for the phased blocks
* Fixed compilation issues on OS X (clang)
* Detect unsorted VCFs and chromosome name mismatches between BAM and VCF
* Fix crash when whatshap compare encounteres unphased VCFs
* Expanded documentation.

v0.13 (2016-10-27)
------------------

* Use ``PS`` tag instead of ``HP`` tag by default to store phasing information.
  This applies to the ``phase`` and ``hapcut2vcf`` subcommands. ``PS`` is also
  used by other tools and standard according to the VCF specification.
* Incorporated genotype likelihoods into our phasing framework. On request
  (by using option ``--distrust-genotypes``), genotypes can now be changed at a cost
  corresponding to their input genotype likelihoods. The changed genotypes are
  written to the output VCF. The behavior of ``--distrust-genotypes`` can be
  fine-tuned by the added options ``--include-homozygous``, ``--default-gq``,
  ``--gl-regularizer``, and ``--changed-genotype-list``.
* Correctly handle cases when processing VCFs with two or more disjoint
  families.

v0.12 (2016-07-01)
------------------

* Speed up allele detection
* Add an ``unphase`` subcommand which removes all phasing from a VCF file
  (``HP`` and ``PS`` tags, pipe notation).
* Add option ``--tag=`` to the ``phase`` subcommand, which allows to choose
  whether ReadBackedPhasing-compatible ``HP`` tags or standard ``PS`` tags are
  used to describe phasing in the output VCF.
* Manage versions with `versioneer <https://github.com/warner/python-versioneer>`_.
  This means that ``whatshap --version`` and the program version in the VCF header
  will include the Git commit hash, such as ``whatshap 0.11+50.g1b7af7a``.
* Add subcommand "haplotag" to tag reads in a BAM file with their haplotype.
* Fix a bug where re-alignment around variants at the very end of a chromosome
  would lead to an AssertionError.

v0.11 (2016-06-09)
------------------

* When phasing a pedigree, blocks that are not connected by reads but
  can be phased based on genotypes will be connected per default. This
  behavior can be turned off using option ``--no-genetic-haplotyping``.
* Implemented allele detection through re-alignment: To detect which allele of a
  variant is seen in a read, the query is aligned to the two haplotypes at that
  position. This results in better quality phasing, especially for
  low-quality reads (PacBio). Enabled if ``--reference`` is provided. Current
  limitation: No score for the allele is computed.
* As a side-effect of the new allele detection, we can now also phase
  insertions, deletions, MNPs and "complex" variants.
* Added option ``--chromosome`` to only work on specifed chromosomes.
* Use constant recombination rate per default, allows to use ``--ped``
  without using ``--genmap``.
* ``whatshap`` has become a command with subcommands. From now on, you need
  to run ``whatshap phase`` to phase VCFs.
* Add a ``stats`` subcommand that prints statistics about phased VCFs.

v0.10 (2016-04-27)
------------------

* Use ``--ped`` to phase pedigrees with the PedMEC algorithm
* Phase all samples in a multi-sample VCF
* Drop support for Python 3.2 - we require at least Python 3.3 now

v0.9 (2016-01-05)
-----------------

* This is the first release available via PyPI (and that can therefore be
  installed via ``pip install whatshap``)

January 2016
------------

* Trio phasing implemented in a branch

September 2015
--------------

* pWhatsHap implemented (in a branch)

April 2015
----------

* Create haplotype-specific BAM files

February 2015
-------------

* Smart read selection

January 2015
------------

* Ability to read multiple BAM files and merge them on the fly

December 2014
-------------

* Logo
* Unit tests

November 2014
-------------

* Cython wrapper for C++ code done
* Ability to write a phased VCF (using HP tags).

June 2014
---------

* Repository for WhatsHap refactoring created

April 2014
----------

* The WhatsHap algorithm is introduced at RECOMB
Various notes
=============

* There is a step in which variants are re-discovered in the BAM file. This may
  fail when the variant caller has used some type of re-alignment (as
  freebayes does). Would be better to integrate this into the variant caller or
  to get the information out of it. This applies only to indels, which are not
  supported right now anyway.
* Input format for HapCompass: http://www.brown.edu/Research/Istrail_Lab/resources/hapcompass_manual.html#sec11


Allele detection with re-alignment
----------------------------------

WhatsHap can detect which allele a read contains at a variant position by
aligining a part of the read to the two possible haplotypes. The haplotype
for which the alignment is better wins.

Allele detection through re-alignment is enabled when the ``--reference``
parameter is used on the command-line.

Re-alignment in this version detects slightly *fewer* alleles than the old
algorithm, but this is typically justified because the old algorithm gave
wrong results. Re-alignment however correctly detects that both haplotypes are
equally good and then refuses to choose.

The alignment algorithm uses edit distance at the moment, which allows us to
detect alleles correctly most of the time, but does not allow us to make use
of base qualities (in fact, the weighted algorithm degenerates into an
unweighted one). To fix this, we need a better alignment algorithm.

Here are some examples for how re-alignment works.

Insertion next to a SNP
~~~~~~~~~~~~~~~~~~~~~~~

Haplotypes::

    ref:   CCTTAGT
    alt:   CCTCAGT

Alignment as reported in BAM file::

    ref:   CCT-TAGT
    query: CCTCAAGT

The second ``T`` is aligned to an ``A``, which is not one of the expected bases.
Thus, no variant would be detected here.

Re-aligning the query to the "alt" haplotype, we get::

    alt:   CCTCA-GT
    query: CCTCAAGT

This alignment has lower cost and we therefore detect that the allele in this
read is probably the alternative one.


Ambiguous
~~~~~~~~~

This was previously detected incorrectly::

    ref:   TGCTTTAAGG
    alt:   TGCTTTCAGG
    query: TGCCTTCAAGG

Two possible alignments are ::

    ref:   TGC-TTTAAGG
    query: TGCCTTCAAGG

and ::

    alt:   TGCTTTCA-GG
    query: TGCCTTCAAGG

Both have cost two and therefore the correct allele is unclear.
.. include:: ../CHANGES.rst
Developing
==========

The `WhatsHap source code is on Bitbucket <https://bitbucket.org/whatshap/whatshap/>`_.
WhatsHap is developed in Python 3, Cython and C++.


Development installation
------------------------

For development, make sure that you install Cython and tox. We also recommend
using a virtualenv. This sequence of commands should work::

	git clone https://bitbucket.org/whatshap/whatshap
	cd whatshap
	python3 -m venv venv
	source venv/bin/activate
	pip install -e .[dev]

The last command installs also all the development dependencies, such as Cython. Use only
``pip install -e .`` to omit those.

Next, you can run WhatsHap like this::

	whatshap --help


Development installation when using Conda
-----------------------------------------

If you are using `Bioconda <https://bioconda.github.io/>`_, it is convenient to develop WhatsHap in a separate environment::

	conda create -n whatshap-dev python=3.6 pysam PyVCF pyfaidx xopen Cython pytest sphinx-issues
	source activate whatshap-dev
	git clone https://bitbucket.org/whatshap/whatshap
	cd whatshap
	pip install -e .

The last command installs WhatsHap into your Conda environment named ``whatshap-dev``. So when
executing ``whatshap`` you will run the latest version you just cloned.

Running tests
-------------

While in the virtual environment, you can run the tests for the current Python version like this::

	pytest

Whenever you change any Cython code (``.pyx`` files), you need to re-run the
``pip install -e .`` step in order to compile it.

Optionally, to run tests for *all* supported Python versions, you can run
`tox <https://tox.readthedocs.io/>`_. It creates separate virtual environments for each Python
version, installs WhatsHap into it and then runs the tests. It also tests documentation generation
with ``sphinx``. Run it like this::

    tox

If ``tox`` is installed on the system, you do not need to be inside a virtual environment for this.
However, you need to have all tested Python versions installed on the system! See the instructions
below for how to do this on Ubuntu.


Installing other Python versions in Ubuntu
------------------------------------------

Ubuntu comes with one default Python 3 version, and in order to test WhatsHap
with older or newer Python versions, follow the instructions for enabling the
`“deadsnakes” repository <https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa>`_.
After you have done so, ensure you have the following packages::

	sudo apt install build-essential python-software-properties

Then get and install the desired Python versions. Make sure you install the ``-dev`` package.
For example, for Python 3.4::

	sudo apt update
	sudo apt install python3.4-dev


Debugging
---------

Here is one way to get a backtrace from gdb (assuming the bug occurs while
running the tests)::

	$ gdb python3
	(gdb) run -m pytest

After you get a SIGSEGV, let gdb print a backtrace:

	(gdb) bt


Wrapping C++ classes
--------------------

The WhatsHap phasing algorithm is written in C++, as are many of the core
data structures such as the “Read” class. To make the C++ classes usable from
Python, we use Cython to wrap the classes. All these definitions are spread
across multiple files. To add new attributes or methods to an existing class
or to add a new class, changes need to be made in different places.

Let us look at the “Read” class. The following places in the code may need to
be changed if the Read class is changed or extended:

* ``src/read.cpp``: Implementation of the class (C++).
* ``src/read.h``: Header with the class declaration (also normal C++).
* ``whatshap/cpp.pxd``: Cython declarations of the class. This repeats – using
  the Cython syntax this time – a subset of the information from the
  ``src/read.h`` file. This duplication is required because Cython
  cannot read ``.h`` files (it would need a full C++ parser for that).

  Note that the ``cpp.pxd`` file contains definitions for *all* the ``.h``
  headers. (It would be cleaner to have them in separate ``.pxd`` files, but
  this leads to problems when linking the compiled files.)
* ``whatshap/core.pxd``: This contains declarations of all *Cython* classes
  wrapping C++ classes. Note that the class ``Read`` in this file has the
  same name as the C++ class, but that it is not the same as the C++ one!
  The distinction is made by prefixing the C++ class with ``cpp.``, which is
  the name of the module in which it is declared in (that is, the C++ class
  ``Read`` is declared in ``cpp.pxd``). The wrapping (Cython) class ``Read``
  stores the C++ class in an attribute named ``thisptr``. If you add a new
  class, it needs to be added to this file. If you only modify an existing one,
  you probably do not need to change this file.
* ``whatshap/core.pyx``: The Cython implementation of the wrapper classes.
  Again, the name ``Read`` by itself is the Python wrapper class and
  ``cpp.Read`` is the name for the C++ class.

Before adding yet more C++ code, which then requires extra code for wrapping it,
consider writing an implementation in Cython instead. See ``readselect.pyx``,
for example, which started out as a Python module and was then transferred to
Cython to make it faster. Here, the Cython code is not merely a wrapper, but
contains the implementation itself.


Writing documentation
---------------------

Documentation is located in the ``doc/`` subdirectory. It is written in
`reStructuredText format <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
and is translated by `Sphinx <http://www.sphinx-doc.org/>`_ into HTML format.

Documentation is hosted on `Read the Docs <https://readthedocs.org/>`_.
In theory, it is built automatically whenever a commit is made. The documentation in the
``master`` branch should be visible at `https://whatshap.readthedocs.io/en/latest/ <https://whatshap.readthedocs.io/en/latest/>`_
and documentation for the most recent released version should be visible at `https://whatshap.readthedocs.io/en/stable/ <https://whatshap.readthedocs.io/en/stable/>`_.
However, the connection between Bitbucket and Read the Docs has never worked
well in this particular project, so builds actually need to be triggered manually
until we have solved this problem.

To generate documentation locally, ensure that you installed sphinx and add-ons
necessary to build documantation (running ``pip install -e .[dev]`` will take
care of this). Then go into the ``doc/`` directory and run ``make``. You can
then open ``doc/_build/html/index.html`` in your browser. The theme that is
used is a bit different from the one the Read the Docs uses.


Making a release
----------------

If this is the first time you attempt to upload a distribution to PyPI, create a
configuration file named ``.pypirc`` in your home directory with the following
contents::

	[distutils]
	index-servers =
	    pypi

	[pypi]
	username=my-user-name
	password=my-password

See also `this blog post about getting started with
PyPI <http://peterdowns.com/posts/first-time-with-pypi.html>`_. In particular,
note that a ``%`` in your password needs to be doubled and that the password
must *not* be put between quotation marks even if it contains spaces.

#. Set the correct version number in the changelog. Ensure that the list of changes is up-to-date.

#. Ensure you have no uncommitted changes in the working copy.

#. Run ``tox``, ensuring all tests pass.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag v0.1

#. Create a distribution (``.tar.gz`` file), ensuring that the auto-generated version number in
   the tarball is as you expect it::

       python3 setup.py sdist

#. Upload the distribution to PyPI (the tarball must be regenerated since ``upload`` requires a preceding ``sdist``)::

       twine upload dist/whatshap-x.yz.tar.gz

   You may need to install the ``twine`` tool to run this command.
#. Push the tag::

       git push --tags

#. Update the `bioconda recipe <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/whatshap/meta.yaml>`_.
   It is probly easiest to edit the recipe via the web interface and send in a
   pull request. Ensure that the list of dependencies (the ``requirements:``
   section in the recipe) is in sync with the ``setup.py`` file.

   Since this is just a version bump, the pull request does not need a
   review by other bioconda developers. As soon as the tests pass and if you
   have the proper permissions, it can be merged directly.

If something went wrong, fix the problem and follow the above instructions again,
but with an incremented revision in the version number. That is, go from version
x.y to x.y.1. Do not change a version that has already been uploaded.


Adding a new subcommand
-----------------------

Follow the instructions in ``whatshap/example.py``.
.. _user-guide:

==========
User guide
==========

WhatsHap is a read-based phasing tool. In the typical case, it expects
1) a VCF file with variants of an individual and 2) a BAM or CRAM file with
sequencing reads from that same individual. WhatsHap uses the sequencing reads
to reconstruct the haplotypes and then writes out the input VCF augmented with
phasing information.

The basic command-line for running WhatsHap is this::

    whatshap phase -o phased.vcf input.vcf input.bam

The reads used for variant calling (to create the input VCF) do not
need to be the same as the ones that are used for phasing. We
recommend that high-quality short reads are used for variant calling and
that the phasing is then done with long reads, see :ref:`the recommended
workflow <recommended-workflow>`.

If the input VCF is a multi-sample VCF, WhatsHap will haplotype all
samples individually. For this, the input must contain reads from all
samples.

:ref:`Multiple BAM/CRAM files can be provided <multiple-bam-files>`,
even from different technologies.

If you want to phase samples of individuals that are related, you can use
:ref:`pedigree phasing <phasing-pedigrees>` mode to improve results.
In this mode, WhatsHap is no longer purely a read-based phasing tool.

With error-prone reads (PacBio, Nanopore), we strongly recommend that you
enable re-alignment mode by providing a reference in FASTA format::

    whatshap phase --reference ref.fasta -o phased.vcf input.vcf input.bam

You can also phase indels by adding the option ``--indels``.

WhatsHap adds the phasing information to the input VCF file and writes it to
the output VCF file. :ref:`See below to understand how phasing information
is represented <phasing_in_vcfs>`.
The VCF file can also be gzip-compressed.


Features and limitations
========================

WhatsHap supports phasing of variants in diploid genomes.

Supported variant types are SNVs (single-nucleotide variants), insertions,
deletions, MNPs (multiple adjacent SNVs) and “complex” variants. Complex
variants are those that do not fall in any of the other categories, but
are not structural variants. An example is the variant TGCA → AAC.
Structural variants are not phased.

If no reference sequence is provided (using ``--reference``), only
SNVs, insertions and deletions can be phased.

All variants in the input VCF that are marked as being heterozygous
(genotype 0/1) and that have appropriate coverage are used as input for the core
phasing algorithm. If the algorithm could determine how the variant should be
phased, that information will be added to the variant in the output VCF.

Variants can be left unphased for two reasons: Either the variant type is
not supported or the phasing algorithm could not make a phasing decision.
In both cases, the information from the input VCF is simply copied to output
VCF unchanged.


Subcommands
===========

WhatsHap comes with the following subcommands.

========== ===================================================
Subcommand Description
========== ===================================================
phase      Phase variants in a VCF with the WhatsHap algorithm
stats      Print phasing statistics
compare    Compare two or more phasings
hapcut2vcf Convert hapCUT output format to VCF
unphase    Remove phasing information from a VCF file
haplotag   Tag reads by haplotype
haplofasta Write haplotypes in FASTA format from a phased VCF
genotype   Genotype variants
========== ===================================================

Not all are fully documented in this manual, yet. To get help for a
subcommand named ``SUBCOMMAND``, run ::

    whatshap SUBCOMMAND --help


.. _recommended-workflow:

Recommended workflow
====================

Best phasing results are obtained if you sequence your sample(s) on both PacBio
and Illumina: Illumina for high-quality variant calls and PacBio for its long
reads.

1. Map your reads to the reference, making sure that you assign each read to a
read group (the ``@RG`` header line in the BAM/CRAM file). WhatsHap supports VCF
files with multiple samples and in order to determine which reads belong to which
sample, it uses the 'sample name' (SM) of the read group. If you have a single
sample only and no or incorrect read group headers, you can run WhatsHap with
``--ignore-read-groups`` instead.

2. Call variants in your sample(s) using the most accurate reads you have. These
will typically be Illumina reads, resulting in a a set of variant calls you can
be reasonably confident in. If you do not know which variant caller to use, yet,
we recommend FreeBayes, which is fast, Open Source and easy to use. In any case,
you will need a standard VCF file as input for WhatsHap in the next step.

3. Run WhatsHap with the VCF file of high-confidence variant calls (obtained in
the previous step) and with the *longest* reads you have. These will typically
be PacBio reads. Phasing works best with long reads, but WhatsHap can use any
read that covers at least two heterozygous variant calls, so even paired-end or
mate-pair reads are somewhat helpful. If you have multiple sets of reads, you
can combine them by providing multiple BAM/CRAM files on the command line.


.. _input-data-requirements:

Input data requirements
=======================

WhatsHap needs correct metadata in the VCF and the BAM/CRAM input files so that
it can figure out which read belongs to which sample. As an example, assume you
give WhatsHap a VCF file that starts like this::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  SampleA  SampleB
    chr1    100  .   A    T    50.0   .       .    GT      0/1      0/1
    ...

WhatsHap sees that there are two samples in it named “SampleA” and “SampleB”
and expects to find the reads for these samples somewhere in the BAM/CRAM file
(or files) that you provide. For that to happen, all reads belonging to a sample
must have the ``RG`` tag, and at the same time, the read group must occur in the
header of the BAM/CRAM file and have the correct sample name. In this example, a
header might look like this::

	@HD     VN:1.4  SO:coordinate
	@SQ     SN:...  LN:...
	...
	@RG   ID:1  SM:SampleA
	@RG   ID:2  SM:SampleB

The ``@RG`` header line will often contain more fields, such as ``PL`` for
the platform and ``LB`` for the library name. WhatsHap only uses the ``SM``
attribute.

With the above header, the individual alignments in the file will be tagged with
a read group of ``1`` or ``2``. For example, an alignment in the BAM/CRAM file
that comes from SampleA would be tagged with ``RG:Z:1``. This is also described
in the `SAM/BAM specification <https://samtools.github.io/hts-specs/>`_.

It is perfectly fine to have multiple read groups for a single sample::

	@RG   ID:1a  SM:SampleA
	@RG   ID:1b  SM:SampleA
	@RG   ID:2   SM:SampleB


What to do when the metadata is not correct
-------------------------------------------

If WhatsHap complains that it cannot find the reads for a sample, then chances
are that the metadata in the BAM/CRAM and/or VCF file are incorrect. You have the
following options:

* Edit the sample names in the VCF header.
* Set the correct read group info in the BAM/CRAM file, for example with the Picard
  tool AddOrReplaceReadGroups.
* Re-map the reads and pass the correct metadata-setting options to your mapping
  tool.
* Use the ``--ignore-read-groups`` option of WhatsHap. In this case, WhatsHap
  ignores all read group metadata in the BAM/CRAM input file(s) and assumes that all
  reads come from the sample that you want to phase. In this mode, you can
  only phase a single sample at a time. If the input VCF file contains more than
  one sample, you need to specify which one to phase by using
  ``--sample=The_Sample_Name``.


.. _multiple-bam-files:

Using multiple input BAM/CRAM files
-----------------------------------

WhatsHap supports reading from multiple BAM or CRAM files. Just provide all BAM
and CRAM files you want to use on the command-line. All the reads across all
those files that to a specific sample are used to phase that sample. This can be
used to combine reads from multiple technologies. For example, if you have
Nanopore reads in one BAM file and PacBio reads in another CRAM file, you can
run the phasing like this::

    whatshap phase -o phased.vcf input.vcf nanopore.bam pacbio.cram

You need to make sure that read group information
:ref:`is accurate in all files <input-data-requirements>`.


.. _vcfs-as-reads:

Using a phased VCF instead of a BAM/CRAM file
---------------------------------------------

It is possible to provide a phased VCF file instead of a BAM/CRAM file. WhatsHap
will then treat the haplotype blocks (:ref:`phase sets <phase-sets>`) it
describes as "reads". For example, if the phased VCF contains only
chromosome-sized haplotypes, then each chromosome would give rise to two such
"reads". These reads are then used as any other read in the phasing algorithm,
that is, they are combined with the normal sequencing reads and the best
solution taking all reads into account is computed.


.. _selection-and-merging:

Read selection and merging
--------------------------

Whatshap has multiple ways to reduce the coverage of the input ---
allowing faster runtimes --- in a way that attempts to minimize the
amount of information lost in this process.  The default behaviour is
to ensure a maximum coverage via read selection: a heuristic that
extracts a subset of the reads that is most informative for phasing.
An optional step which can be done before selection is to merge
subsets of reads together to form superreads according to a
probabilistic model of how likely subsets of reads are to appear
together on the same haplotype (p_s) or different haplotypes (p_d).
By default, this feature is not activated, however it can be activated
by specifying the ``--merge-reads`` flag when running ``whatshap
phase``.  This model is parameterized by the following four parameters

====================== ======================================================
Parameter              Description
====================== ======================================================
error-rate             Probability that a nucleotide is wrong
maximum-error-rate     Maximum error any edge of the merging graph can have
threshold              Threshold ratio of p_s/p_d to merge two sets
negative-threshold     Threshold ratio of p_d/p_s to not merge two sets
====================== ======================================================

which can be specified by the respective flags ``--error-rate=0.15``,
``--maximum-error-rate=0.25``, ``--threshold=100000`` and
``--negative-threshold=1000`` (note that defaults are shown here for
example) when running ``whatshap phase``.


.. _phasing_in_vcfs:

Representation of phasing information in VCFs
=============================================

WhatsHap supports two ways in which it can store phasing information in a VCF
file: The standards-compliant ``PS`` tag and the ``HP`` tag used by GATK’s
ReadBackedPhasing tool. When you run ``whatshap phase``, you can select which
format is used by setting ``--tag=PS`` or ``--tag=HP``.

We will use a small VCF file as an example in the following. Unphased, it
looks like this::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0/1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      0/1      1/1
    chr1    300  .   G    T    50.0   .       .    GT      0/1      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0/1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0/1      1/1

Note that sample 1 is heterozygous at all shown loci (expressed with
``0/1`` in the ``GT`` field).


Phasing represented by pipe (``|``) notation
--------------------------------------------

The ``GT`` fields can be phased by ordering the alleles by haplotype and
separating them with a pipe symbol (``|``) instead of a slash (``/``)::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0|1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      1|0      0/1
    chr1    300  .   G    T    50.0   .       .    GT      1|0      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0|1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0|1      1/1

The alleles on one of the haplotypes of sample1 are: A, G, T, T, A.
On the other haplotype, they are: T, C, G, A, G.

Swapping ones and zeros in the ``GT`` fields would result in a VCF file with
the equivalent information.


.. _phase-sets:

Phasing represented by PS ("phase set") tag
-------------------------------------------

The pipe notation has problems when not all variants in the VCF file can be
phased. The `VCF specification <https://github.com/samtools/hts-specs>`_
introduces the ``PS`` tag to solve some of them. The ``PS`` is a
unique identifier for a "phase set", which is a set of variants that were
be phased relative to each other. There are usually multiple phase sets in
the file, and variants that belong to the same phase set do not need to
be consecutive in the file::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1      sample2
    chr1    100  .   A    T    50.0   .       .    GT:PS:PQ   0|1:100:22   0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:PS:PQ   1|0:100:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:PS:PQ   1|0:300:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:PS:PQ   0|1:300:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:PS:PQ   0|1:100:12   0/1:.:.

This VCF contains two phase sets named ``100`` and ``300``. The names are
arbitrary, but WhatsHap will choose the position of the leftmost variant
of the phase set as its name. The variants at 100, 150 and 500 are in the same
phase set, while the variants at 300 and 350 are in a different phase set.
Such a configuration is typically seen when paired-end or mate-pair reads are
used for phasing.

In the case of WhatsHap, the phase sets are identical to the connected components
of the variant connectivity graph. Two variants in that graph are connected if a
read exists that covers them.

The above example also shows usage of the ``PQ`` tag for "phasing quality".
WhatsHap currently does not add this tag.


Phasing represented by HP tag
-----------------------------

GATK’s ReadBackedPhasing tool uses a different way to represent phased variants.
It is in principle the same as the combination of pipe notation with the ``PS``
tag, but the ``GT`` field is left unchanged and all information is added to a
separate ``HP`` tag ("haplotype identifier") instead. This file encodes the same
information as the example above::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1         sample2
    chr1    100  .   A    T    50.0   .       .    GT:HP      0/1:100-1,100-2      0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:HP:PQ   0/1:100-2,100-1:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:HP:PQ   0/1:300-2,300-1:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:HP:PQ   0/1:300-1,300-2:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:HP:PQ   0/1:100-1,100-2:12   0/1:.:.

A few notes:

* ReadBackedPhasing does not add the ``PQ`` to the first variant in a phase set/haplotype
  group. This probably means that the phasing quality is to be interpreted as relative to
  the previous or first variant in the set.
* ReadBackedPhasing does not phase indels
* Discussions on the GATK forum on this topic:
   - https://gatkforums.broadinstitute.org/discussion/4226
   - https://gatkforums.broadinstitute.org/discussion/4038/


Trusting the variant caller
===========================

WhatsHap will trust the variant caller to have made the right decision of
whether a variant is heterozygous or homozygous. If you use the option
``--distrust-genotypes``, then this assumption is softened: An optimal solution
could involve switching a variant from being heterozygous to homozygous.
Currently, if that option is enabled and such a switch occurs, the variant
will simply appear as being unphased. No change of the genotype in the VCF is
done.

If you use this option, fewer variants will be phased.

Note that switching homozygous variants to heterozygous is never possible since
only heterozygous variants are considered for phasing.


.. _phasing-pedigrees:

Phasing pedigrees
=================

When phasing multiple samples from individuals that are related (such as
parent/child or a trio), then it is possible to provide WhatsHap with
a ``.ped`` file that describes the pedigree. WhatsHap will use the
pedigree *and* the reads to infer a combined, much better phasing.

To turn on pedigree mode, run WhatsHap like this::

	whatshap phase --ped pedigree.ped -o phased.vcf input.vcf input.bam

where ``pedigree.ped`` is a plink-compatible PED file to describe the
relationships between samples and ``input.vcf`` is a multi-sample VCF
with all individuals that should be phased. The reads for all individuals
can be in one or more BAM/CRAM files. WhatsHap will match them based on sample
names provided in the read groups (just like for the default single-individual
mode). 
In the resulting VCF file (``phased.vcf``), 
haplotype alleles of a child are given as paternal|maternal, i.e.
the first allele is the one inherited from the father and the second one
the allele inherited from the mother.

PED file format
---------------

WhatsHap recognizes `PLINK-compatible PED
files <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>`_.
A PED file is a white-space (space or tab) delimited file with at least six
columns. WhatsHap checks the column count, but uses only

  * column 2: individual ID
  * column 3: paternal ID
  * column 4: maternal ID

The other columns are ignored. Lines starting with ``#`` are considered
comments and are ignored. Empty lines are also ignored.

To define a single trio, it is sufficient to have a single row in the PED file
with the child, mother and father. It is *not* necessary to include "dummy" rows
for individuals whose parents are unknown. (You will currently get a warning if
you do, but this will be changed.)

Here is an example defining a trio::

    # Fields: family, individual_id, paternal_id, maternal_id, sex, phenotype
    FAMILY01 the_child father mother 0 1

A quartet (note how multiple consecutive spaces are fine)::

    # Fields: family, individual_id, paternal_id, maternal_id, sex, phenotype
    FAMILY01 one_child   father mother 0 1
    FAMILY01 other_child father mother 0 1

*Important*: The names in the PED file *must* match the sample names in your VCF
and BAM/CRAM files!

Pedigree phasing parameters
---------------------------

Phasing in pedigree mode requires costs for recombination events. Per
default, WhatsHap will assume a constant recombination rate across the
chromosome to be phased. The recombination rate (in cM/Mb) can be
changed by providing option ``--recombrate``. The default value of
1.26 cM/Mb is suitable for human genomes.

In order to use region-specific recombination rates, a genetic map file
can be provided via option ``--genmap``. WhatsHap expects a three-column
text file like this::

	position COMBINED_rate(cM/Mb) Genetic_Map(cM)
	55550 0 0
	568322 0 0
	568527 0 0
	721290 2.685807669 0.410292036939447
	723819 2.8222713027 0.417429561063975
	723891 2.9813105581 0.417644215424158
	...

The first (header) line is ignored and the three columns are expected to
give the pysical position (in bp), the local recombination rate between the
given position and the position given in the previous row (in cM/Mb), and
the cumulative genetic distance from the start of the chromosome (in cM).
The above example was taken from the 1000 Genomes genetic map `provided by
SHAPEIT
<https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap>`_.
Since genetic map files provide information for only one chromosome, the
``--genmap`` option has to be combined with ``--chromosome``.


Creating phased references in FASTA format
==========================================

To reconstruct the two haplotypes that a phased VCF describes, the
``bcftools consensus`` command can be used. It is part of
`bcftools <http://www.htslib.org/>`_. As input, it expects a reference
FASTA file and either an indexed BCF or a compressed and indexed VCF file.
To work with the uncompressed VCF output that WhatsHap produces, proceed
as follows::

    bgzip phased.vcf
    tabix phased.vcf.gz
    bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
    bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta

Here, ``reference.fasta`` is the reference in FASTA format and ``phased.vcf``
is the phased VCF. Afterwards, ``haplotype1.fasta`` and ``haplotype2.fasta``
will contain the two haplotypes.

.. note:
    If there are problems in the input VCF, bcftools (as of version 1.3) may
    not give an error message and instead create files that are identical to
    the input ``reference.fasta``. As a precaution, you may want to make sure
    that the two haplotype FASTA files are indeed different from the input
    reference FASTA.


Visualizing phasing results
===========================

Sometimes it is helpful to visually inspect phasing results by looking at them
in a genome browser. The steps here assume that you use the Integrative Genomics
Viewer (IGV).


GTF with haplotype blocks
-------------------------

WhatsHap can create a GTF file from a phased VCF file that describes the
haplotype blocks. With phasing results in ``phased.vcf``, run ::

    whatshap stats --gtf=phased.gtf phased.vcf

WhatsHap will print some statistics about the phasing in the VCF, and it
will also create the file ``phased.gtf``.

Now open both ``phased.vcf`` and ``phased.gtf`` in IGV in order to inspect the
haplotype block structure. In this example, there are four haplotype blocks and
it is clear which variants they connect:

|

.. image:: _static/gtf.png

|

Haplotype blocks can be interleaved or nested if mate-pair or paired-end reads
are used for phasing. In the GTF track, you will note this because the blocks
appear as “exons” connected with a horizontal line (not shown in the screenshot).


Coloring reads
--------------

It is often a lot more interesting to also show the reads along with the
variants.

For that, run the ``whatshap haplotag`` subcommand on your phased VCF file. It
tags each read in a BAM file with ``HP:i:1`` or ``HP:i:2`` depending on which
haplotype it belongs to, and also adds a ``PS`` tag that describes in which
haplotype block the read is. With your aligned reads in ``alignments.bam``,
run ::

    whatshap haplotag -o haplotagged.bam --reference reference.fasta phased.vcf alignments.bam

The ``haplotag`` commands re-detects the alleles in the reads in the same way
the main ``phase`` command does it. Since availability of a reference influences
how this is done, if you used ``--reference`` with your ``phase`` command, you
should alse use ``--reference`` here.

When using 10X Genomics BAM files, ``haplotag`` reads the BX tags and per default
assigns reads that belong to the same read cloud to the same haplotype. 
This feature can be switched off using the ``--ignore-linked-read`` flag.

The input VCF may have been phased by any program, not only WhatsHap, as long as
the phasing info is recorded with a ``PS`` or ``HP`` tag.

Also, the reads in the input BAM file do not have to be the ones that were used
for phasing. That is, you can even phase using one set of reads and then assign
haplotypes to an entirely different set of reads (but from the same sample).

The command above creates a BAM file ``haplotagged.bam`` with the tagged reads,
which you can open in IGV.

To visualize the haplotype blocks, right click on the BAM track and choose
*Color Alignments by* → *tag*. Then type in ``PS`` and click “Ok”. Here is an
example of how this can look like. From the colors of the reads alone,
it is easy to see that there are four haplotype blocks.

|

.. image:: _static/haplotagged-PS.png

|

You can also visualize the haplotype assignment. For that, choose
*Color Alignments by* → *tag* and type in ``HP``. Additionally, you may want to
also sort the alignments by the ``HP`` tag using the option *Sort Alignments by*
in the right-click context menu.

Here is an impression of how this can look like. The reads colored in red belong
to one haplotype, while the ones in blue belong to the other. Gray reads are
those that could not be tagged, usually because they don’t cover any
heterozygous variants.

|

.. image:: _static/haplotagged-HP.png

|

Genotyping Variants
===================

Besides phasing them, WhatsHap can also re-genotype variants. Given a VCF file
containing variant positions, it computes genotype likelihoods for all three
genotypes (0/0, 0/1, 1/1) and outputs them in a VCF file together with a
genotype prediction. Genotyping can be run using the following command::

    whatshap genotype -o genotyped.vcf input.vcf input.bam

The predicted genotype is stored in the output VCF using the ``GT`` tag and the ``GL`` tag
provides (log10-scaled) likelihoods computed by the genotyping algorithm.
As for phasing, providing a reference sequence is strongly recommended in order to
enable re-alignment mode::

    whatshap genotype --reference ref.fasta -o genotyped.vcf input.vcf input.bam

.. _installation:

============
Installation
============

Installation of WhatsHap is easiest if you use Conda.


Installation with Conda
-----------------------

First, ensure you have Conda (miniconda or Anaconda) installed and made the
proper settings to enable the “bioconda” channel. For that, follow
`the bioconda installation instructions <https://bioconda.github.io/#install-conda>`_.

Then WhatsHap can be installed with this command::

    conda install whatshap

If you have Conda, but not enabled bioconda, use this command::

    conda install -c bioconda -c conda-forge whatshap


Installation with pip
---------------------

Before you can `pip install`, you need to install dependencies that pip cannot
install for you. WhatsHap is implemented in C++ and Python. You need to have a
C++ compiler, Python 3.3 (or later) and the corresponding Python header files.
Python 3.5 is slightly faster than 3.4. In Ubuntu, installing the packages
``build-essential`` and ``python3-dev`` will take care of all required
dependencies.

WhatsHap can then be installed with pip::

	pip3 install --user whatshap

This installs WhatsHap into ``$HOME/.local/bin``.  Then add
``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

    export PATH=$HOME/.local/bin:$PATH
    whatshap --help

Alternatively, you can also install WhatsHap into a virtual environment if you
are familiar with that.


Installing an unreleased development version
--------------------------------------------

If you for some reason want to use the most recent development version of
WhatsHap, you can install it in the following way. These instructions will
create a virtual environment in the directory ``whatshap-env`` that contains
WhatsHap. Simply delete that directory to uninstall the software. Other WhatsHap
versions you may have installed in other locations remain unaffected. ::

	python3 -m venv whatshap-env
	whatshap-env/bin/pip install Cython
	whatshap-env/bin/pip install git+https://bitbucket.org/whatshap/whatshap

You can then run WhatsHap like this::

	whatshap-env/bin/whatshap --version

You should see a version number like ``0.17+103.g71e5b3c``, which means that
this version is 103 Git commits ahead of version 0.17.
=====================
Questions and Answers
=====================

 * **Can WhatsHap use a reference panel?** Reference panels are used by population-based phasers (like Beagle or ShapeIt). Although we are considering integrating this, WhatsHap cannot take advantage of reference panels right now. In case you have population data, we suggest to produce a population-based phasing (using `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_, `ShapeIt <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_, etc.) and a read-based phasing using WhatsHap separately and then compare/integrate results in a postprocessing step.
 * **Will Illumina data lead to a good read-based phasing?** Illumina paired-end data is not ideal for read-based phasing, since most pairs of heterozygous SNPs will not be bridged by a read pair. However, WhatsHap will attempt to produce as long haplotype blocks as possible. Running WhatsHap will hence tell you how phase-informative your input data is. Just take a look at the number (and sizes) of produced haplotype blocks.
 * **How large can/should a pedigree be for pedigree-aware read-based phasing (i.e. using option** :code:`--ped` **)?** The pedigree mode in WhatsHap is intended for intermediate-size pedigrees. The runtime of the core phasing step will be linear in :math:`2^{2t}`, where :math:`t` is the number of trio relationships (= number of children) in your pedigree. We do not recommend to use pedigrees with :math:`t>5`. For such pedigrees, read-data is unnecessary in most cases anyway and a very high-quality phasing can be obtained by genetic haplotyping methods (like `MERLIN <http://csg.sph.umich.edu/abecasis/Merlin/tour/haplotyping.html>`_).
.. include:: README.rst

=================
Table of contents
=================

.. toctree::
   :maxdepth: 2

   installation
   guide
   faq
   develop
   notes
   howtocite
   changes


..
   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
.. image:: https://img.shields.io/pypi/v/whatshap.svg?branch=master
    :target: https://pypi.python.org/pypi/whatshap

.. image:: https://semaphoreci.com/api/v1/whatshap/whatshap/branches/master/shields_badge.svg
    :target: https://semaphoreci.com/whatshap/whatshap

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io/recipes/whatshap/README.html

|

.. image:: https://bitbucket.org/repo/8AjxBd/images/543323972-whatshap_logo.png

|

WhatsHap
========

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.

Please :ref:`cite us if you use WhatsHap <howtocite>`.


Features
========

  * Very accurate results (Martin et al.,
    `WhatsHap: fast and accurate read-based phasing <https://doi.org/10.1101/085050>`_)
  * Works well with Illumina, PacBio, Oxford Nanopore and other types of reads
  * It phases SNVs, indels and even “complex” variants (such as ``TCG`` → ``AGAA``)
  * Pedigree phasing mode uses reads from related individuals (such as trios)
    to improve results and to reduce coverage requirements
    (Garg et al., `Read-Based Phasing of Related Individuals <https://doi.org/10.1093/bioinformatics/btw276>`_).
  * WhatsHap is :ref:`easy to install <installation>`
  * It is :ref:`easy to use <user-guide>`: Pass in a VCF and one or more BAM files, get out a phased VCF.
    Supports multi-sample VCFs.
  * It produces standard-compliant VCF output by default
  * If desired, get output that is compatible with ReadBackedPhasing
  * Open Source (MIT license)


Documentation
-------------

* `Bitbucket page <https://bitbucket.org/whatshap/whatshap/>`_
* `Read the documentation online <https://whatshap.readthedocs.io/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the
  repository and in the downloaded tar distribution.


Mailing list
------------
We run a `public mailing list <https://lists.cwi.nl/mailman/listinfo/whatshap>`_. Please
don't hesitate to post questions and comments.
.. _howtocite:

How to cite
===========

.. include:: ../CITATION.rst
