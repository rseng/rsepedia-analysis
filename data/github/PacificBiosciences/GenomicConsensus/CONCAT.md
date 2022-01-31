<h1 align="center"><img src="http://www.pacb.com/wp-content/themes/pacific-biosciences/img/pacific-biosciences-logo-mobile.svg"/></h1>
<h1 align="center">GenomicConsensus</h1>
<p align="center">Genome polishing and variant calling.</p>

***

The ``GenomicConsensus`` package provides the ``variantCaller`` tool,
which allows you to apply the Quiver or Arrow algorithm to mapped
PacBio reads to get consensus and variant calls.

## Availability
Latest version can be installed via bioconda package `genomicconsensus`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Background on Quiver and Arrow

*Quiver* is the legacy consensus model based on a conditional random
field approach.  Quiver enables consensus accuracies on genome
assemblies at accuracies approaching or even exceeding Q60 (one error
per million bases).  If you use the HGAP assembly protocol in
SMRTportal 2.0 or later, Quiver runs automatically as the final
"assembly polishing" step.

Over the years Quiver has proven difficult to train and develop, so we are
phasing it out in favor of the new model, Arrow.  *Arrow* is an
improved consensus model based on a more straightforward hidden Markov
model approach.

Quiver is supported for PacBio RS data.  Arrow is supported for PacBio
Sequel data and RS data with the P6-C4 chemistry.


## Running
-------
Basic usage is as follows:

```sh
% quiver aligned_reads{.cmp.h5, .bam, .fofn, or .xml}    \
>     -r reference{.fasta or .xml} -o variants.gff       \
>     -o consensus.fasta -o consensus.fastq
```

``quiver`` is a shortcut for ``variantCaller --algorithm=quiver``.
Naturally, to use arrow you could use the ``arrow`` shortcut or
``variantCaller --algorithm=arrow``.

in this example we perform haploid consensus and variant calling on
the mapped reads in the ``aligned_reads.bam`` which was aligned to
``reference.fasta``.  The ``reference.fasta`` is only used for
designating variant calls, not for computing the consensus.  The
consensus quality score for every position can be found in the output
FASTQ file.

*Note that 2.3 SMRTanalysis does not support "dataset" input (FOFN
 or XML files); those who need this feature should wait for the forthcoming
 release of SMRTanalysis 3.0 or build from GitHub sources.*


## More documentation

- [More detailed installation and running instructions](./doc/HowTo.rst)
- [FAQ](./doc/FAQ.rst)
- [variants.gff spec](./doc/VariantsGffSpecification.rst)
- [CHANGELOG](./CHANGELOG)

DISCLAIMER
----------
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
Quiver FAQ
==========

What are EviCons? GenomicConsensus? Quiver? Plurality?  
------------------------------------------------------------
**GenomicConsensus** is the current PacBio consensus and variant
calling suite.  It contains a main driver program, ``variantCaller``,
which provides two consensus/variant calling algorithms: **Arrow** and
**Quiver**.  These algorithms can be run by calling
``variantCaller.py --algorithm=[arrow|quiver|plurality]`` or by going
through the convenience wrapper scripts ``quiver`` and ``arrow``.

**EviCons** was the previous generation PacBio variant caller (removed
 in software release v1.3.1).

Separate packages called **ConsensusCore** and **ConsensusCore2** are
C++ libraries where all the computation behind Quiver and Arrow are
done, respectively.  This is transparent to the user after
installation.


What is Plurality?
------------------
**Plurality** is a very simple variant calling algorithm: it stacks up the
aligned reads (alignment as produced by BLASR, or alternate mapping
tool), and for each column under a reference base, calls the most
abundant (i.e., the plurality) read base (or bases, or deletion) as
the consensus at that reference position.


Why is Plurality a weak algorithm?
----------------------------------
Plurality does not perform any local realignment.  This means it is
heavily biased by the alignment produced by the mapper (BLASR,
typically).  It also means that it is insensitive at detecting indels.
Consider this example::

    Reference    AAAA
                 ----
      Aligned    A-AA
        reads    AA-A
                 -AAA
                 ----
    Plurality    AAAA
    consensus

Note here that every read has a deletion and the correct consensus
call would be "AAA", but due to the mapper's freedom in gap-placement
at the single-read level, the plurality sequence is "AAAA"---so the
deletion is missed.  Local realignment, which plurality does not do,
but which could be considered as implicit in the Quiver algorithm,
essentially pushes the gaps here to the same column, thus identifying
the deletion.  While plurality could be adjusted to use a simple "gap
normalizing" realignment, in practice noncognate extras (spurious
non-homopolymer base calls) in the midst of homopolymer runs pose
challenges.

What is Quiver?
---------------
**Quiver** is a more sophisticated algorithm that finds the maximum
quasi-likelihood template sequence given PacBio reads of the template. 
PacBio reads are modeled using a conditional random field approach that
scores the quasi-likelihood of a read given a template sequence.  In
addition to the base sequence of each read, Quiver uses several
additional *QV* covariates that the basecaller provides.  Using these
covariates provides additional information about each read, allowing
more accurate consensus calls.

Quiver does not use the alignment provided by the mapper (BLASR,
typically), except for determining how to group reads together at a
macro level.  It implicitly performs its own realignment, so it is
highly sensitive to all variant types, including indels---for example,
it resolves the example above with ease.

The name **Quiver** reflects a consensus-calling algorithm that is
`QV-aware`.

We use the lowercase "quiver" to denote the quiver *tool* in GenomicConsensus,
which applies the Quiver algorithm to mapped reads to derive sequence 
consensus and variants.

Quiver is described in detail in the supplementary material to the
`HGAP paper`_.


What is Arrow?
--------------
Arrow is a newer model intended to supercede Quiver in the near future. 
The key differences from Quiver are that it uses an HMM model instead 
of a CRF, it computes true likelihoods, and it uses a smaller set of
covariates.  We expect a whitepaper on Arrow to be available soon.

We use the lowercase "arrow" to denote the arrow *tool*, which applies
the Arrow algorithm to mapped reads to derive sequence 
consensus and variants.


What data works with Arrow/Quiver?
----------------------------------

Quiver supports all PacBio *RSII* data since the C2 chemistry release
(late 2012).

Arrow supports all data from the PacBio Sequel instrument, as well as
data from the PacBio RSII using the P6-C4 chemistry.


How do I run `quiver`/`arrow`?
------------------------------
For general instructions on installing and running, see the
HowTo_ document.



What is the output from `quiver`/`arrow`?
-----------------------------------------
There are three output files from the GenomicConsensus tools:

1. A consensus *FASTA* file containing the consensus sequence
2. A consensus *FASTQ* file containing the consensus sequence with quality annotations
3. A variants *GFF* file containing a filtered, annotated list of variants identified

It is important to note that the variants included in the output
variants GFF file are *filtered* by coverage and quality, so not all
variants that are apparent in comparing the reference to the consensus
FASTA output will correspond to variants in the output variants GFF
file.

To enable all output files, the following can be run (for example)::

    % quiver -j16 aligned_reads.cmp.h5 -r ref.fa \
     -o consensus.fa                             \
     -o consensus.fq                             \
     -o variants.gff

The extension is used to determine the output file format.


What does it mean that `quiver` consensus is *de novo*?
-------------------------------------------------------
Quiver's consensus is *de novo* in the sense that the reference and the reference
alignment are not used to inform the consensus output.  Only the reads
factor into the determination of the consensus.

The only time the reference sequence is used to make consensus calls -
when the ``--noEvidenceConsensusCall`` flag is set to ``reference`` or
``lowercasereference`` (the default)- is when there is no effective
coverage in a genomic window, so Quiver has no evidence for computing
consensus.  One can set ``--noEvidenceConsensusCall=nocall`` to
avoid using the reference even in zero coverage regions.


What is the expected `quiver` accuracy?
---------------------------------------
Quiver's expected accuracy is a function of coverage and chemistry.
The C2 chemistry (no longer available), P6-C4 and P4-C2 chemistries
provide the most accuracy.  Nominal consensus accuracy levels are as
follows:

+----------+-------------------------------+
|Coverage  |Expected consensus accuracy    |
|          +------------------+------------+
|          | C2, P4-C2, P6-C4 | P5-C3      |
+==========+==================+============+
|10x       | > Q30            | > Q30      |
+----------+------------------+------------+
|20x       | > Q40            | > Q40      |
+----------+------------------+------------+
|40x       | > Q50            | > Q45      |
+----------+------------------+------------+
|60-80x    | ~ Q60            | > Q55      |
+----------+------------------+------------+

The "Q" values referred to are Phred-scaled
quality values:

.. math::
   q = -10 \log_{10} p_{error}

for instance, Q50 corresponds to a p_error of 0.00001---an accuracy
of 99.999%.  These accuracy expectations are based on routine
validations performed on multiple bacterial genomes before each
chemistry release.


What is the expected accuracy from `arrow`
------------------------------------------
`arrow` achieves similar accuracy to `quiver`.  Numbers will be published soon.


What are the residual errors after applying `quiver`?
-----------------------------------------------------

If there are errors remaining applying Quiver, they will almost
invariably be homopolymer run-length errors (insertions or deletions).



Does `quiver`/`arrow` need to know what sequencing chemistry was used?
----------------------------------------------------------------------

At present, the Quiver model is trained per-chemistry, so it is very
important that Quiver knows the sequencing chemistries used.

If SMRT Analysis software was used to build the `cmp.h5` or BAM input file, the
`cmp.h5` will be loaded with information about the sequencing
chemistry used for each SMRT Cell, and GenomicConsensus will automatically
identify the right parameters to use.

If custom software was used to build the `cmp.h5`, or an
override of Quiver's autodetection is desired,  then the
chemistry or model must be explicity entered. For example::

  % quiver -p P4-C2 ...
  % quiver -p P4-C2.AllQVsMergingByChannelModel ...



Can a mix of chemistries be used in a cmp.h5 file for quiver/arrow?
-------------------------------------------------------------------

Yes!  GenomicConsensus tools automatically see the chemistry *per-SMRT Cell*, so it
can figure out the right parameters for each read and model them
appropriately.


What chemistries and chemistry mixes are supported?
---------------------------------------------------


For Quiver: all PacBio RS chemistries are supported.  Chemistry
mixtures of P6-C4, P4-C2, P5-C3, and C2 are supported.

For Arrow: the RS chemistry P6-C4, and all PacBio Sequel chemistries
are supported.  Mixes of these chemistries are supported.



What are the QVs that the Quiver model uses?
--------------------------------------------
Quiver uses additional QV tracks provided by the basecaller.  
These QVs may be looked at as little breadcrumbs that are left behind by
the basecaller to help identify positions where it was likely that
errors of a given type occurred.  Formally, the QVs for a given read are
vectors of the same length as the number of bases called; the QVs
used are as follows:

  - DeletionQV
  - InsertionQV
  - MergeQV
  - SubstitutionQV
  - DeletionTag

To find out if your cmp.h5 file is loaded with these QV tracks, run the command
::

    % h5ls -rv aligned_reads.cmp.h5

and look for the QV track names in the output.  If your cmp.h5 file is
lacking some of these tracks, Quiver will still run, though it will
issue a warning that its performance will be suboptimal.


Why is `quiver`/`arrow` making errors in some region?
-----------------------------------------------------
The most likely cause for *true* errors made by these tools is that the
coverage in the region was low.  If there is 5x coverage over a
1000-base region, then 10 errors in that region can be expected.

It is important to understand that the effective coverage available to
`quiver`/`arrow` is not the full coverage apparent in plots---the tools
filter out ambiguously mapped reads by default.  The
remaining coverage after filtering is called the /effective coverage/.
See the next section for discussion of `MapQV`.

If you have verified that there is high effective coverage in the region
in question, it is highly possible---given the high accuracy quiver and arrow
can achieve---that the apparent errors actually
reflect true sequence variants.  Inspect the FASTQ output file to
ensure that the region was called at high confidence; if an erroneous
sequence variant is being called at high confidence, please report a
bug to us.


What does Quiver do for genomic regions with no effective coverage?
-------------------------------------------------------------------
For regions with no effective coverage, no variants are outputted, and
the FASTQ confidence is 0.

The output in the FASTA and FASTQ consensus sequence tracks is
dependent on the setting of the ``--noEvidenceConsensusCall`` flag.
Assuming the reference in the window is "ACGT", the options are:

+---------------------------------------------+---------+
|``--noEvidenceConsensusCall=...``            |Consensus|
|                                             |output   |
+=============================================+=========+
|``nocall`` (default in 1.4)                  |NNNN     |
+---------------------------------------------+---------+
|``reference``                                |ACGT     |
+---------------------------------------------+---------+
|``lowercasereference`` (new post 1.4, and the|         |
|default)                                     |acgt     |
+---------------------------------------------+---------+




What is `MapQV` and why is it important?
----------------------------------------
`MapQV` is a single scalar Phred-scaled QV per aligned read that
reflects the mapper's degree of certainty that the read aligned to
*this* part of the reference and not some other.  Unambigously mapped
reads will have a high `MapQV` (typically 255), while a read that was
equally likely to have come from two parts of the reference would have
a `MapQV` of 3.

`MapQV` is pretty important when you want highly accurate variant
calls.  Quiver and Plurality both filter out aligned reads with a
MapQV below 20 (by default), so as not to call a variant using data of
uncertain genomic origin.

This can be problematic if using quiver/arrow to get a consensus
sequence.  If the genome of interest contains long (relative to the library
insert size) highly-similar repeats, the effective coverage (after
`MapQV` filtering) may be reduced in the repeat regions---this is termed
these `MapQV` dropouts.  If the coverage is sufficiently reduced in
these regions, quiver/arrow will not call consensus in these regions---see
`What do quiver/arrow do for genomic regions with no effective coverage?`_.

If you want to use ambiguously mapped reads in computing a consensus
for a denovo assembly, the `MapQV` filter can be turned off entirely.
In this case, the consensus for each instance of a genomic repeat will
be calculated using reads that may actually be from other instances of
the repeat, so the exact trustworthiness of the consensus in that
region may be suspect.  The next section describes how to disable the
`MapQV` filter.


How can the `MapQV` filter be turned off and when should it be?
--------------------------------------------------------------
The `MapQV` filter can be disabled using the flag
``--mapQvThreshold=0`` (shorthand: ``-m=0``).  If running a
quiver/arrow job via SMRT Portal, this can be done by unchecking the "Use
only unambiguously mapped reads" option. Consider this in
de novo assembly projects, but it is not recommended for variant
calling applications.


How can variant calls made by quiver/arrow be inspected or validated?
---------------------------------------------------------------------
When in doubt, it is easiest to inspect the region in a tool like
SMRT View, which enables you to view the reads aligned to the region.
Deletions and substitutions should be fairly easy to spot; to view
insertions, right-click on the reference base and select "View
Insertions Before...".


What are the filtering parameters that quiver/arrow use?
--------------------------------------------------------

The available options limit read coverage, filters reads by `MapQV`, and filters
variants by quality and coverage.

- The overall read coverage used to call consensus in every window is
  100x by default, but can be changed using ``-X=value``.
- The `MapQV` filter, by default, removes reads with MapQV < 20.  This
  is configured using ``--mapQvThreshold=value`` / ``-m=value``
- Variants are only called if the read coverage of the site exceeds
  5x, by default---this is configurable using ``-x=value``.
  Further, they will not be called if the confidence (Phred-scaled)
  does not exceed 40---configurable using ``-q=value``.


What happens when the sample is a mixture, or diploid?
-----------------------------------------------------
At present, quiver/arrow assume a haploid sample, and the behavior of
on sample mixtures or diploid/polyploid samples is
*undefined*.  The program will not crash, but the output results are
not guaranteed to accord with any one of the haplotypes in the sample,
as opposed to a potential patchwork.  


Why would I want to *iterate* the mapping+(quiver/arrow) process?
-----------------------------------------------------------------
Some customers using quiver for polishing highly repetitive genomes
have found that if they take the consensus FASTA output of quiver, use
it as a new reference, and then perform mapping and Quiver again to
get a new consensus, they get improved results from the second round
of quiver.

This can be explained by noting that the output of the first round of
quiver is more accurate than the initial draft consensus output by the
assembler, so the second round's mapping to the quiver consensus can
be more sensitive in mapping reads from repetitive regions.  This can
then result in improved consensus in those repetitive regions, because
the reads have been assigned more correctly to their true genomic
loci.  However there is also a possibility that the potential shifting
of reads around from one rounds' mapping to the next might alter
borderline (low confidence) consensus calls even away from repetitive
regions.

We recommend the (mapping+quiver) iteration for customers polishing
repetitive genomes, and it could also prove useful for resequencing
applications.  However we caution that this is very much an
*exploratory* procedure and we make no guarantees about its
performance.  In particular, borderline consensus calls can change
when the procedure is iterated, and the procedure is *not* guaranteed
to be convergent.


Is iterating the (mapping+quiver/arrow) process a convergent procedure?
-----------------------------------------------------------------------
We have seen many examples where (mapping+quiver), repeated many
times, is evidently *not* a convergent procedure.  For example, a
variant call may be present in iteration n, absent in n+1, and then
present again in n+2.  It is possible for subtle changes in mapping to
change the set of reads examined upon inspecting a genomic window, and
therefore result in a different consensus sequence there.  We expect
this to be the case primarily for "borderline" (low confidence) base
calls.



.. _HowTo: ./HowTo.rst
.. _`HGAP paper`: http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2474.html

``variants.gff`` File Format (Version 2.1)
============================================

As of this version, ``variants.gff`` is our primary variant call file
format.  The ``variants.gff`` file is based on the `GFFv3 standard`_.
The GFFv3 standard describes a tab-delimited plain-text file
meta-format for describing genomic "features."  Each gff file consists
of some initial "header" lines supplying metadata, and then a number
of "feature" lines providing information about each identified
variant.

The GFF Coordinate System
-------------------------

All coordinates in GFF files are 1-based, and all intervals ``start,
end`` are understood as including both endpoints.

Headers
-------

The ``variants.gff`` file begins with a block of metadata headers,
which looks like the following:

::

    ##gff-version 3
    ##pacbio-variant-version 2.1
    ##date Tue Feb 28 17:44:18 2012
    ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
    ##source GenomicConsensus v0.1.0
    ##source-commandline callVariants.py --algorithm=plurality aligned_reads.cmp.h5 -r spinach.fasta -o variants.gff
    ##source-alignment-file /home/popeye/data/aligned_reads.cmp.h5
    ##source-reference-file /home/popeye/data/spinach.fasta
    ##sequence-region EGFR_Exon_23 1 189
    ##sequence-header EGFR_Exon_24 1 200

The ``source`` and ``source-commandline`` describe the name and
version of the software generating the file.
``pacbio-variant-version`` reflects the specification version that the
file contents should adhere to.

  The ``sequence-region`` headers describe the names and extents of
the reference groups (i.e. reference contigs) that will be refered to
in the file.  The names are the same as the full FASTA header.

``source-alignment-file`` and ``source-reference-file`` record
absolute paths to the primary input files.


Feature lines
-------------

After the headers, each line in the file describes a genomic
*feature*; in this file, all the features are potential variants
flagged by the variant caller.  The general format of a variant line
is a 9-column (tab-delimited) record, where the first 8 columns
correspond to fixed, predefined entities in the GFF standard, while
the 9th column is a flexible semicolon-delimited list of mappings
``key=value``.

The 8 predefined columns are as follows:

+------+-------+--------------------------------+------------------+
|Column|Name   |Description                     |Example           |
|Number|       |                                |                  |
+------+-------+--------------------------------+------------------+
|1     |seqId  |The full FASTA header for the   |``lambda_NEB3011``|
|      |       |reference contig.               |                  |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|2     |source |(unused; always populated with  |``.``             |
|      |       |``.``)                          |                  |
+------+-------+--------------------------------+------------------+
|3     |type   |the type of variant.  One of    |``substitution``  |
|      |       |``insertion``, ``deletion``, or |                  |
|      |       |``substitution``.               |                  |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|4     |start  |1-based start coordinate for the|200               |
|      |       |variant.                        |                  |
+------+-------+--------------------------------+------------------+
|5     |end    |1-based end coordinate for the  |215               |
|      |       |variant.  start<=end always     |                  |
|      |       |obtains, regardless of strand.  |                  |
+------+-------+--------------------------------+------------------+
|6     |score  |unused; populated with ``.``    |``.``             |
+------+-------+--------------------------------+------------------+
|7     |strand |unused; populated with ``.``    |``.``             |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|8     |phase  |unused; populated with ``.``    |``.``             |
+------+-------+--------------------------------+------------------+


The attributes in the 9th (final) column are as follows:

+--------------+----------------------------+-----------------+
|Key           |Description                 |Example          |
|              |                            |value            |
+--------------+----------------------------+-----------------+
|``coverage``  |the read coverage of the    |``42``           |
|              |variant site (not the       |                 |
|              |variant itself)             |                 |
+--------------+----------------------------+-----------------+
|``confidence``|the phred-scaled probability|``37``           |
|              |that the variant is real,   |                 |
|              |rounded to the nearest      |                 |
|              |integer and truncated at 93 |                 |
+--------------+----------------------------+-----------------+
|``reference`` |the reference base or bases |``T``, ``.``     |
|              |for the variant site.  May  |                 |
|              |be ``.`` to represent a     |                 |
|              |zero-length substring (for  |                 |
|              |insertion events)           |                 |
+--------------+----------------------------+-----------------+
|``variantSeq``|the read base or bases      |``T``            |
|              |corresponding to the        | (haploid);      |
|              |variant. ``.`` encodes a    |``T/C``, ``T/.`` |
|              |zer-length string, as for a | (heterozygous)  |
|              |deletion.                   |                 |
+--------------+----------------------------+-----------------+
|``frequency`` |the read coverage of the    |``13``           |
|              |variant itself; for         | (haploid)       |
|              |heterozygous variants, the  |                 |
|              |frequency of both observed  |``15/12``        |
|              |alleles.  This is an        | (heterozygous)  |
|              |optional field.             |                 |
+--------------+----------------------------+-----------------+


The attributes may be present in any order.

The four types of variant we support are as follows. *(Recall that the
field separator is a tab, not a space.)*

1. Insertion.  Examples::

    ref00001 . insertion 8 8 . . . reference=.;variantSeq=G;confidence=22;coverage=18;frequency=10
    ref00001 . insertion 19 19 . . . reference=.;variantSeq=G/.;confidence=22;coverage=18;frequency=7/5

  For insertions, start==end, and the insertion event is understood as
  taking place *following* the reference position `start`.

2. Deletion.  Examples::

    ref00001 . deletion 348 349 . . . reference=G;variantSeq=.;confidence=39;coverage=25;frequency=20
    ref00001 . deletion 441 443 . . . reference=GG;variantSeq=GG/.;confidence=39;coverage=25;frequency=8/8

3. Substitution.  Examples::

    ref000001 . substitution 100 102 . . . reference=GGG;variantSeq=CCC;confidence=50;coverage=20;frequency=16
    ref000001 . substitution 200 201 . . . reference=G;variantSeq=G/C;confidence=50;coverage=20;frequency=10/6



Compression
-----------

The gff metaformat is verbose, so for practical purposes we will gzip
encode ``variants.gff`` files as ``variants.gff.gz``.  Consumers of
the variant file should be able to read it in either form.


Other file formats
------------------

The VCF and BED standards describe variant-call specific file formats.
We can currently translate `variants.gff` files to these formats, but
they are not the primary output of the variant callers.


.. _GFFv3 standard: http://www.sequenceontology.org/gff3.shtml


Variant Caller Functional Specification
=======================================

Version 3.3


Introduction
------------

This document describes the interface, input/output, and performance
characteristics of ``variantCaller``, a variant calling tool
provided by the ``GenomicConsensus`` package.


Software Overview
-----------------

The ``GenomicConsensus`` package provides a command-line tool,
``variantCaller``, which provides several consensus and variant-calling  algorithms for
PacBio sequencing data.


Functional Requirements
-----------------------

Command-line interface
``````````````````````

``variantCaller`` is invoked from the command line.  For example, a simple
invocation is::

        variantCaller -j8 --algorithm=arrow  \
                         -r lambdaNEB.fa     \
                         -o variants.gff     \
                         aligned_subreads.bam

which requests that variant calling proceed,
- using 8 worker processes,
- employing the **arrow** algorithm,
- taking input from the file ``aligned_subreads.bam``,
- using the FASTA file ``lambdaNEB.fa`` as the reference,
- and writing output to ``variants.gff``.

A particularly useful option is ``--referenceWindow/-w``: this option
allows the user to direct the tool to perform variant calling
exclusively on a *window* of the reference genome.

Invoking

::

    variantCaller --help

will provide a help message explaining all available options; they will be
documented here shortly.



Input and output
````````````````
``variantCaller`` requires two input files:

- A sorted file of reference-aligned reads in `PacBio's standard BAM format`_;
- A FASTA file adhering to `PacBio's FASTA file conventions`_

The input file is the main argument to ``variantCaller``, while the
output files are provided as arguments to the ``-o`` flag.  For
example,

::

        variantCaller aligned_subreads.bam -r lambda.fa  -o myVariants.gff -o myConsensus.fasta

will read input from ``aligned_subreads.bam``, using the reference
``lambda.fa``, and send variant call output to the file
``myVariants.gff``, and consensus output to ``myConsensus.fasta``.
The extension of the filename provided to the ``-o`` flag is
meaningful, as it determines the output file format.  The file formats
presently supported, by extension, are

``.gff``
        PacBio GFFv3 variants format; convertable to VCF or BED.

``.fasta``
        FASTA file recording the consensus sequence calculated for each reference contig

``.fastq``
        FASTQ file recording the consensus sequence calculated for
        each reference contig, as well as per-base confidence scores


.. note::

   The *quiver* and *arrow* algorithms require that certain metrics
   are in place in the input BAM file.

   *quiver*, which operates on RSII data only, requires the
   basecaller-computed "pulse features" ``InsertionQV``,
   ``SubstitutionQV``, ``DeletionQV``, and ``DeletionTag``.  These
   features are populated in BAM tags by the ``bax2bam`` conversion
   program.

   *arrow*, which operates on RSII P6-C4 data and all Sequel data,
   requires per-read SNR metrics, and the per-base ``PulseWidth``
   metric for Sequel data (but not for RSII P6-C4).  These metrics are
   populated by Sequel instrument software or the ``bax2bam``
   converter (for RSII data).

   The selected algorithm will halt with an error message if features
   it requires are unavailable.


Available algorithms
````````````````````

At this time there are two algorithms available for variant calling:
**plurality**, **quiver**, and **arrow**.

**Plurality** is a simple and very fast procedure that merely tallies
the most frequent read base or bases found in alignment with each
reference base, and reports deviations from the reference as potential
variants.  This is a very insensitive and flawed approach for PacBio
sequence data, which is prone to insertion and deletion errors.

**Quiver** is a more complex procedure based on algorithms originally
developed for CCS.  Quiver leverages the quality values (QVs) provided by
upstream processing tools, which provide insight into whether
insertions/deletions/substitutions were deemed likely at a given read
position.  Use of **quiver** requires the ``ConsensusCore``
library.

**Arrow** is the successor to Quiver; it uses a more principled HMM
model approach.  It does not require basecaller quality value metrics;
rather, it uses the per-read SNR metric and the per-pulse pulsewidth
metric as part of its likelihood model.  Beyond the model specifics,
other aspects of the Arrow algorithm are similar to Quiver.  Use of
**arrow** requires the ``ConsensusCore2`` library, which is provided
by the ``unanimity`` codebase.


Software interfaces
```````````````````
The ``GenomicConsensus`` module has two essential dependencies:

1. **pbcore**, the PacBio Python bioinformatics library
2. **ConsensusCore**, a C++ library with SWIG bindings that provides
   access to the Quiver algorithm.
3. **ConsensusCore2**, a C++ library with SWIG bindings that provides access to
   the Arrow algorithm.

These modules are easily installed using their ``setup`` scripts,
which is the canonical means of installing Python packages.


Confidence values
-----------------

The arrow*, *quiver*, and *plurality* algorithms make a confidence
metric available for every position of the consensus sequence.  The
confidence should be interpreted as a phred-transformed posterior
probability that the consensus call is incorrect; i.e.

.. math::

    QV = -10 \log_{10}(p_{err})

``variantCaller`` clips reported QV values at 93---larger values
cannot be encoded in a standard FASTQ file.



Chemistry specificity
---------------------

The Quiver and Arrow algorithm parameters are trained per-chemistry.
Quiver and Arrow identify the sequencing chemistry used for each run
by looking at metadata contained in the data file (the input BAM or
cmp.h5 file).  This behavior can be overriden by a command line flag.

When multiple chemistries are represented in the reads in the input
file, Quiver/Arrow will model each read appropriately using the
parameter set for its chemistry, thus yielding optimal results.


Performance Requirements
------------------------

``variantCaller`` performs variant calling in parallel using multiple
processes.  Work splitting and inter-process communication are handled using
the Python ``multiprocessing`` module.  Work can be split among an arbitrary
number of processes (using the ``-j`` command-line flag), but for best
performance one should use no more worker processes than there are CPUs in the
host computer.

The running time of the *plurality* algorithm should not exceed the
runtime of the BLASR process that produced the cmp.h5. The running
time of the *quiver* algorithm should not exceed 4x the runtime of
BLASR.

The amount of core memory (RAM) used by a ``variantCaller`` run should
not exceed 2GB per active CPU core (as selected using the ``-j`` flag).


.. _PacBio's standard BAM format: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
.. _PacBio's FASTA file conventions: http://pacbiofileformats.readthedocs.io/en/3.0/FASTA.html

How to install and use GenomicConsensus
=======================================

**We recommend that you obtain GenomicConsensus by installing the most
recent version of SMRTanalysis.  Other means of installation are not
officially supported.**


Basic running instructions
--------------------------

Basic usage---using 8 CPUs to compute consensus of mapped reads and
variants relative to a reference---is as follows::

    % quiver -j8 aligned_reads{.cmp.h5, .bam, .fofn, or .xml} \
    >     -r reference{.fasta or .xml} -o variants.gff        \
    >     -o consensus.fasta -o consensus.fastq

``quiver`` is a shortcut for ``variantCaller --algorithm=quiver``.
Naturally, to use arrow you could use the ``arrow`` shortcut or
``variantCaller --algorithm=arrow``.

in this example we perform haploid consensus and variant calling on
the mapped reads in the ``aligned_reads.bam`` which was aligned to
``reference.fasta``.  The ``reference.fasta`` is only used for
designating variant calls, not for computing the consensus.  The
consensus quality score for every position can be found in the output
FASTQ file.

*Note that 2.3 SMRTanalysis does not support "dataset" input (FOFN
or XML files); those who need this feature should wait for the forthcoming
release of SMRTanalysis 3.0 or build from GitHub sources.*


Running a large-scale resequencing/polishing job in SMRTanalysis 2.3
--------------------------------------------------------------------

We do not recommend attempting  to construct a single giant cmp.h5 file and
then processing it on a single node.  This is inefficient and users attempting to do this
have run into many problems with the instability of the HDF5 library (which PacBio is
moving away from, in favor of BAM_.)

To run a large-scale resequencing job (>50 megabase genome @ 50x
coverage,nominally), you want to spread the computation load across
multiple nodes in your computing cluster.  

The `smrtpipe` workflow engine in SMRTanalysis 2.3 provides a
convenient workflow automating this---it will automatically spread the
load for both mapping and quiver jobs among your available cluster
nodes.  This is accessible via the SMRTportal UI; the simplest way to 
set up and run thse workflows is via tha UI.  Nonetheless, we include 
command-line instructions for completeness.

If you have to run the `smrtpipe` workflow manually from the command
line, a recipe is as folows::

1. Make sure the reference you will align and compare against is
   present in a SMRTportal "reference repository".  Even if you
   don't want to use SMRTportal, you need to build/import the
   reference appropriately, and the simplest way to do that is
   via SMRTportal.  If you don't have a SMRTportal instance, 
   you can use the ``referenceUploader`` command to prepare your
   reference repository.

2. Prepare an "input.fofn" file listing, one-per-line, each "bax.h5"
   file in your input data set.

3. Convert the "input.fofn" to an "input.xml" file that SMRTpipe can
   understand::

   $ fofnToSmrtpipeInput.py input.fofn > input.xml

4. Prepare your "params.xml" file.  Here is a `params.xml template`_
   you can use; you should just need to edit the reference path.

5. Activate your SMRTanalysis environment, and invoke smrtpipe::

   $ source <SMRT Analysis>/etc/setup.sh
   $ smrtpipe.py --distribute --params=params.xml xml:input.xml

6. After successful execution is complete, the results should be
   available as `data/consensus.fast[aq].gz` and
   `data/variants.gff.gz`, etc.

Please consult the `SMRTpipe reference manual`_ for further information.

*Note that resequencing (mapping reads against a reference genome and
then calling consensus and identifying variants) and polishing
(mapping reads against a draft assembly and then taking the consensus
output as the final, polished, assembly) are the same algorithmic
operation, the only effective difference is that the "variants.gff"
output is not biologically meaningful in the polishing case---it just
records the edits that were made to the draft to produce the polished
assembly.*

Running a large-scale quiver/arrow job in SMRTanalysis 3.0+
-----------------------------------------------------------

(Forthcoming)


Building bleeding-edge code (unsupported)
----------------------------------------

If you need to access the the latest code for some reason, a
convenient way to build it is to use PacBio's pitchfork_ build
system, which will take care of all third party dependencies for you.
Here's a recipe::

  git clone git@github.com:PacificBiosciences/pitchfork.git
  cd pitchfork
  make GenomicConsensus   # may take some time, as it builds dependencies...

Now, with GenomicConsensus built, you can use it via::

  bash --init-file deployment/setup-env.sh  # Puts you in a subshell where your build is available
  quiver --help                             # now you have quiver, arrow, etc. available

If you encounter build issues using `pitchfork`, please report the
issues there.  Note that you can deploy PacBio software to a location
of your choice using pitchfork.


Further questions?
------------------

Please consult the `FAQ document`_.

.. _`FAQ document`: ./FAQ.rst
.. _pitchfork : https://github.com/PacificBiosciences/pitchfork
.. _`params.xml template`: ./params-template.xml
.. _`SMRTpipe reference manual`: http://www.pacb.com/wp-content/uploads/2015/09/SMRT-Pipe-Reference-Guide.pdf
.. _`BAM`: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html

Known Issues
============

Python 2.6 multiprocessing is susceptible to a bug where exceptions
are occasionally thrown at shutdown because the daemon processes are
allowed to continue executing while the interpreter is shutting down.
(See: http://bugs.python.org/issue4106, http://bugs.python.org/issue9207)
The bug is fixed in 2.7 but not in 2.6.  I haven't been able to find a
workaround.

.. GenomicConsensus documentation master file, created by
   sphinx-quickstart on Sat Jan 28 18:28:19 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GenomicConsensus
================

Contents:

.. toctree::
   :maxdepth: 2

   VariantCallerFunctionalSpecification
   VariantsGffSpecification
   VariantCallerKnownIssues
   HowToQuiver
   QuiverFAQ
Variant Caller Validation Specification
=======================================

Created: 2012/02/20
Author: jdrake

Synopsis
--------
There are several algorithms implemented for detecting SNPs in the data using alignments
against a known reference.  These include Evicons (PacBio), GATK (Broad) and, most recently,
GenomicConsensus (PacBio).  The first was built back in the days of 70% accurate reads
and is quickly becoming deprecated, though currently used only as our haploid caller (though it
does have diploid calling functionality).  The second is part of a comprehensive tool kit
that provides better diploid calling and is currently used as such in the secondary pipeline.
The third is the most recent incarnation and the heir apparent going forward.

There are no metrics built to measure, for example, the sensitivy and specificity of these
algorithms, and thus are difficult to evaluate against eachother.  Ostensibly, since we're
closer to the data, we should be able to better tune an algorithm to maximize true +/-
variant calls.  This exercise will create datasets to generate ROC curves and potentially
other user metrics to properly evaluate algorithms.

Workflow
--------
A dataset will be generated using a set of mutated reference genomes to align real reads
against.  Each mutated reference will have a list of 'ground-truth' mutations associated
with them.  The alignments will then be processed by each of the candidate variant caller
algorithms and their results evaluated against the ground truth for true +/-.

1. Generate mutated reference(s)
2. Align reads to each mutated reference
3. Run variant callers using alignments
4. Evaluate calls vs ground truth
5. Generate metrics
6. Repeat steps 3 - 5 as necessary

*NOTE*: The mutated references could be generated on the fly using the mutation information,
thus, obviating the need to save all the mutated references.

The automation of the workflow should eventually be packaged and deployed into the Millhouse
framework. Initially, it can configured to run from the command line on the secondary cluster
and the smrtpipe infrastructure.

Mutation Sets
-------------
Starting with lambda (well represented amongst currently available runs), generate a set, M, of
mutated lambda references with n randomly generated point mutations within each m mutated genome.
Each point mutation p will be one of P = {Substitution(s), Insertion(i), Deletion(d)}. Locations will
be associated with each mutation as a 1-based offset using the wild-type genome (w) coordinates.

Offsets are stored in 0-based coordinates, but displayed and manipulated using a 1-based coordinate
system because that's what GFF, SAM and Tablet uses.

Mutation sets will be stored in a file per reference mutated.  The file will be used as input to
mutate genomes on the fly just prior to alignment as well as input to the validation procedure.
GFF files could be generated from them fairly easily if, though probably not, necessary.  Multiple
versions of this file could co-exist for the same reference.

The format of the mutations file is extremely simple and compact making it suitable for source
control.  For simplicity, we'll use the python pickle protocol vers 2.

mutation = {
    offs, # offset in the wild-type genome (1-based)
    typ,  # type of mutation
    wild, # base(s), wild-type strain
    mut   # base(s), the mutated strain
}

Comparison
~~~~~~~~~~

Alignments play a big role in this.  Homopolymer regions are treated slightly differently when
QV's are involved.  Without them, affine gaps are used which push gaps to the left, e.g.,

GAATGAAGCCG
GAATG-AGCCG

These degenerate cases will be handled by collapsing the homopolymer aligment gaps to the left before
comparing.  See BLASR subsection for more detail.

Some alignments against a substitution generate this type of call (some field removed for brevity):
deletion 10201 length=1;confidence=22;coverage=16;reference=G
insertion 10201 length=1;variantSeq=T;confidence=23;coverage=16

This happens when the substitution forms a homopolymer.  These will be collapsed and labelled substitutions
during comparisons.

Does it use CCS reads by default, QV values?  The production version does not use CCS reads but does use
QV values.


Data Sets
---------
Key things to pay attention to in datasets:
- coverage level
- quality
- location of mutation (e.g., homopolymer)
- nature of mutation (e.g., insertion followed by deletion)

Start with a positive control using an unmutated reference.  Zero mutations should be found.

Metrics
-------
Confusion matrix
ROC Curves (using quality scores)
QQ Plot

Notes
-----
http://www.sequenceontology.org/gff3.shtml
http://smrtwiki/wiki/SMRTpipe
http://web/~mhsieh/dag/index.html

Evicons
~~~~~~~
Top level module src: //depot/software/assembly/pbpy/pbpy/smrtpipe/modules
Evicons smrtpipe modules: P_Consensus, P_ConsensusAlgorithm
    wraps runChunkConnectPost.py (same dir) which ...
    wraps eviConsWrapper.py (../../../bin/) which ...
    wraps jConnectPost[.sh] (../../../../seymour/dist2/analysis/bin/) which ...
    wraps a call an evicons jar file
*Un*-wrapping this may be more cumbersome than generating the appropriate inputs to the module.

GATK
~~~~
P_GATKVC

Uses the UnifiedGenotyper, TableRecalibration, CountCovariates components

Uses BAM inputs, generated after alignment (blasr)
//depot/software/bioinformatics/tools/pbsamtools

BLASR
~~~~~
Running blasr to get a cmp.h5 file (super basic, with crappy alignments)::

> compareSequences.py --algorithm=blasr --h5fn=aligned.cmp.h5 input.fofn refdir

More productiony way::

> compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=6  --noXML --h5mode=w \
    --h5fn=control.cmp.h5 --minAccuracy=0.75 --minLength=50  -x -minMatch 12 -x -bestn 1 -x -minPctIdentity 70.0 \
    --regionTable=/mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/jobs/037/037285/data/filtered_regions.chunk001of002.fofn \
    input.fofn /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/references/lambda


`refdir` is a directory containing a set of information related to a reference sequence.
The key files appear to be <reference>.fa and reference.info.xml.  It can work with just
a fasta file, but will produce a cmp.h5 that breaks evicons (reference length is 0). There
is a utility to generate these ref dirs:

/mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_LastSuccessfulBuild/analysis/bin/referenceUploader


Validation tests could be source controlled under the siv tree, given they're likely to
transition into that group eventually (//depot/software/assembly/siv-test/...)

Using what we've already got:
//depot/software/assembly/siv-test/module-test/bin/
- mutateRef.py (?)
- evalVariantCalls.py (?)

1.3.3 Enhancements
==================

Bug 20100
---------
**Genomic Consensus to support rare variant calling**

Adds the ability to to call rare variants.  Rare variants are defined here as
mutations detected at a 1% < frequency < 50%.  There is an initial minimum
coverage requirement set at 500x.  The information provided by this feature
will be limited to deviations from the reference.

This will limited to SNPs only.  Indels will be ignored.

Codon-aware filtering could easily be applied to the output as a post-processing
step.  It could also be used *in situ* as an additional filtering mechanism to
reduce potentially noisy output.  We can start with the post-processing option
(easy) then evolve towards being codon aware as necessary.

This functionality will be optional.

*Inputs*:

    A column-oriented set of short (~1 - 3 bases) sequence calls and their
    corresponding frequencies.

*Outputs*:

    A GFF-style record per rare variant call. See VariantsGffSpecification for
    standard format. This feature will augment the standard record with a new
    key/value pair indicating frequency. Example: freq=0.10.  There may be more
    than one variant per reference position.

    Please note that no consensus file(s) will be generated for rare variants,
    though enough information is provided to build one in a separate
    tools/module.

Bug 20628
---------
**Add support for detecting and reporting correlated mutations**

Provides support for determing whether or not sets of mutations are co-located.
This only includes SNPs, not indels.  Correlations may only be established
using overrlapping reads, i.e., gaps not supportable. Correlations will have
some confidence metric associated with them (TBD).  This functionality may also
be combined with rare variant calling output.

The guiding use-case for this feature is the BCR-ABL project, which targeted an
863 bp region of the human genome (11,000x coverage). `BCR-ABL Project Details`_.

This functionality will be optional.

*Inputs*:

    CCS-based variant calls at each position including read information: ID,
    start, stop.  'start' and 'stop' are in (+) strand genomic coordinates.

*Outputs*:

    A table (possibly) of correlated mutations that could look like:

    =====  =======  =====  ===================
    Freq   # Reads  Conf   Mutations
    =====  =======  =====  ===================
    40.4%  4,321    40     123a, 140c
    30.3%  3,210    30     50t, 350a
    20.2%  2,500    20     1400g, 1500c, 1550t
    =====  =======  =====  ===================

    We may also choose to include an output of read IDs associated with reported
    sets of co-located mutations.  Formats TBD.

.. _BCR-ABL Project Details: http://web/~mbrown/workspace2011Q4/bcrAblASHRuns/
