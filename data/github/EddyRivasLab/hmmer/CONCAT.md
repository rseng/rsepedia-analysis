## HMMER - biological sequence analysis using profile HMMs

[![](https://travis-ci.org/EddyRivasLab/hmmer.svg?branch=develop)](https://travis-ci.org/EddyRivasLab/hmmer)
![](http://img.shields.io/badge/license-BSD-brightgreen.svg)

[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models" (profile HMMs). HMMER is used by many
protein family domain databases and large-scale annotation pipelines,
including [Pfam](http://pfam.xfam.org) and other members of the
[InterPro Consortium](http://www.ebi.ac.uk/interpro/).

To obtain HMMER releases, please visit [hmmer.org](http://hmmer.org).

To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).


### to download and build the current source code release:

```
   % wget http://eddylab.org/software/hmmer/hmmer.tar.gz
   % tar zxf hmmer.tar.gz
   % cd hmmer-3.3.2
   % ./configure --prefix /your/install/path
   % make
   % make check                 # optional: run automated tests
   % make install               # optional: install HMMER programs, man pages
   % (cd easel; make install)   # optional: install Easel tools
``` 

Executable programs will be installed in `/your/install/path/bin`. If
you leave this optional `./configure` argument off, the default prefix
is `/usr/local`.

Files to read in the source directory:

   * INSTALL - brief installation instructions.
   * Userguide.pdf - the HMMER User's Guide.
 
To get started after installation, see the Tutorial section in the
HMMER User's Guide (Userguide.pdf).



### to clone a copy of HMMER3 source from github:

The tarball way, above, is a better way to install HMMER (it includes
a precompiled Userguide.pdf, for example), but you can also clone our
github repo. You need to clone both the HMMER and Easel repositories,
as follows:

```
   % git clone https://github.com/EddyRivasLab/hmmer
   % cd hmmer
   % git clone https://github.com/EddyRivasLab/easel
   % autoconf
```

and to build:

```bash
   % ./configure
   % make
```

Our [git workflow](https://github.com/EddyRivasLab/hmmer/wiki/Git-workflow)
includes three main branches:

 * **master** is the stable branch for HMMER3 releases (including when
   H3 is released as a library inside Infernal)
 * **develop** is the HMMER3 development branch
 * **h4-develop** is the HMMER4 development branch.

To build the most recent official release, leave both HMMER and Easel
on their default **master** branch.  To contribute to HMMER3
development, you want to be on the **develop** branches. If you want
to send us a pull request on GitHub, please base your changes on our
**develop** branches.


### to report a problem:

Visit our
[issues tracking page at github](https://github.com/EddyRivasLab/hmmer/issues).

# HMMER 3.3.2 release notes (Nov 2020)


## bug fixes:

* Fixed a recently introduced bug that could cause hmmsearch (and
  presumably hmmscan) to segfault on rare comparisons involving highly
  biased sequences. In domain postprocessing
  (p7_domaindef.c::rescore_isolated_domain()), when p7_Decoding()
  returns an eslERANGE error on garbage sequences in long_target mode,
  nhmmer needs to reset its background model, which it modifies on the
  fly. This was in the fix for iss #198, which we added in the 3.3.1
  release. However, that fix failed to check for long_target mode,
  which introduced this new bug for hmmsearch/hmmscan.

* ./configure --enable-PIC wasn't setting -fPIC option in impl-sse,
  impl-vmx Makefiles.  
  (Thanks to Martin Larralde.)
  
* fixed an uninitialized ptr in makehmmerdb, in the fm_data
  structure, which could cause makehmmerdb to crash.  
  (Thanks to Augustin Zidek.)



## new thingies:

* added a `make install-strip` target.  
  (Thanks to Sebastien Jaenicke.)



For more information, you can peruse the
[git log for our master (stable release) branch](https://github.com/EddyRivasLab/hmmer/commits/master).


# HMMER 3.3.1 release notes (Jul 2020)


## most important changes:

* Default sequence weighting behavior was slightly changed, making it
  independent of whether an RF reference annotation line is present on
  the input alignment or not. Previously, we used the RF line (if
  present) to define consensus columns. This would mean that
  reformatting the same alignment from Stockholm format (with
  reference column annotation) to aligned FASTA format (without it)
  would result in different profile HMM parameters, because of small
  differences in sequence weighting. (iss #180)

* Revised how nhmmer guesses whether an input query file consists of
  single sequence(s), multiple sequence alignment(s), or profile
  HMM(s). `--qformat` allows user to specify input sequence file format
  (unaligned or aligned formats) and override the guesser;
  `--qsingle_seqs` allows an MSA format to be read as set(s) of single
  query sequences instead of as MSA(s). Also, fixed a problem where
  nhmmer would run on a protein target sequence file without complaint.
  (iss #171, iss #196, PR #194)

## bug fixes:

* The comparison engine used by hmmsearch, hmmscan, phmmer, and
  jackhmmer has a design limit of 100K residues for target sequence
  length (usually protein sequences, sometimes RNA transcripts).  Only
  nhmmer/nhmmscan are designed for searching arbitrary length genome
  sequences. That limit was not being enforced, and it was possible
  for a user to run hmmsearch inadvertently instead of nhmmer, which
  can cause numerical overflows and give infinite scores, among other
  problems. The comparison engine now exits with an error if the
  target sequence length is >100K.

* nhmmer now gives a less confusing error message when it tries to
  open a sequence file as FMindex format and fails. (iss #195, PR
  #200)

* fixed a problem where nhmmer failed to reset its background model
  composition after p7_Decoding() fails with a range error on highly
  repetitive sequence, causing remainder of target sequence
  comparisons to be performed with a corrupted bg. (iss #198, PR
  #199).

* removed a stray printf() from nhmmer.

* fixed a problem where if p7_Decoding() failed with a range error on
  a highly repetitive sequence, a trace structure wasn't being reused
  properly, and it could throw an exception with `trace not empty;
  needs to be Reuse()'d?`. (PR #191)

* `nhmmscan --aliscoreout` was segfaulting. The `--aliscoreout` option
  was only present in nhmmer for some experiment anyway, and was never 
  hooked up in hmmscan. Fixed by removing the option. (iss #190)

* fixed p7_gmx_utest unit test failures (segfaults and malloc
  failures) for large L,M inputs. (iss #176)

* fixed a problem where nhmmer's sequence windows, when scanning a
  long target sequence, depended on previous sequences in the file,
  thus meaning that subseq windows depended on target database
  order. In some cases, this could affect whether hits were or were
  not found. The fix introduces a flag in esl_sqio_ReadBlock() that
  makes sequence window choices reproducible. (PR #174)

* fixed possible buffer overflow in
  p7_tophits_GetMaxPositionLength(). (Did not occur in practice;
  required a terabase of target sequence to exercise.)

* fixed p7_hmmd_search_stats unit test problems with memory corruption
  and leaks, detected on OS/X. 
  
* fixed problem with `make check` failing on hmmc2 on FreeBSD.


For even more information, you can peruse the
[git log for our develop branch](https://github.com/EddyRivasLab/hmmer/commits/develop).


# HMMER 3.3 release notes (Nov 2019)


## most important changes:

* We improved the `hmmpgmd` search daemon. (`hmmpgmd` is the compute
  server used at EBI to support their HMMER web servers.)  
  Now `hmmpgmd` handles large target sequence databases more
  efficiently by using "sharding". Instead of loading the entire
  target sequence database into memory on every node, a target
  database can be sharded across multiple nodes. The `hmmpgmd` daemon
  also now has its own user guide.

* We improved how we calculate our default sequence weights (Henikoff
  position-based weights), especially on deep alignments of 10-100K+
  sequences. Now we calculate PB weights only on consensus columns,
  not all columns. This avoids some cases of artifactually extreme
  weights on very gappy alignments. We also changed to a better rule
  for defining sequence fragments.  These changes give a small
  improvement in sensitivity/specificity benchmarking of HMMER/Pfam
  searches, and substantial speed improvements in building profiles on
  deep alignments. Because of these changes, profile HMMs built with
  version 3.3 give slightly different scores compared to previous
  HMMER3 versions.

* Fixed a bug in the `hmmstat` "compKL" (composition KL divergence)
  calculation, which was off by a constant. Now it is reported as a
  standard KL divergence in bits.


## bug fixes:

* fixed a bug where in some rare and complicated situations, `nhmmer`
  could report overlapping envelopes on reverse strand. [iss#159]

* Several bugs were fixed in MPI mode, including iss#157 (`hmmsim
  --mpi` was always segfaulting, even on simple examples) and iss#154
  (problems in `hmmscan` that were corrupting output data).

* Fixed some 32-bit integer overflow bugs in `hmmpgmd`. 

* Fixed a bug in the `hmmstat` "compKL" (composition KL divergence)
  calculation, which was off by a constant. Now it is reported as a
  standard KL divergence in bits.

* Fixed a bug where `hmmconvert --outfmt` wasn't recognizing `3/f` as
  a format.


## Smaller changes

* Our fasta format parser now detects aligned FASTA format (.afa
  files) more robustly, and will not attempt to read a .afa file as an
  unaligned sequence file. [iss#153]

* Our `make check` tests depend on Python >= 3.5. Added checks in
  `./configure` and `make` to fail gracefully if python3 isn't available.

* `./configure` now always calls `AC_PROG_CC_STDC` to add compiler flags
  for C99 (even with icc).

* Removed undocumented `--{tmm,tmi,tmd,tim,tii,tdm,tdd}` options from
  hmmbuild.

* Data serialization routines (used in `hmmpgmd`) were rewritten and
  improved.


For even more information, you can peruse the
[git log for our develop branch](https://github.com/EddyRivasLab/hmmer/commits/develop).


