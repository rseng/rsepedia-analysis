<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/img/logo.png width=150px>

[![PyPI version shields.io](https://img.shields.io/pypi/v/epitopepredict.svg)](https://pypi.python.org/pypi/epitopepredict/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/epitopepredict/badge/?version=latest)](https://epitopepredict.readthedocs.io/en/latest/?badge=latest)

### Background

**epitopepredict** provides a standardized programmatic interface and command line tool for executing multiple epitope prediction methods. Currently this largely consists of interfaces to several MHC binding prediction, the results of which can then be processed and visualized in a consistent manner. There is a built-in method for MHC-class I prediction and the TEPITOPEPan method is provided as a 'built in' method for MHC-class II. The IEDB tools and netMHCpan, netMHCIIpan and MHCFlurry are also supported. Those tools are free for academic use but have to be installed separately. This software runs on most linux systems.

Documentation is at http://epitopepredict.readthedocs.io

### Installation

current release:

`pip install epitopepredict`

or for latest version on github:

`pip install -e git+https://github.com/dmnfarrell/epitopepredict.git#egg=epitopepredict`
<img src=https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/img/logo.png width=150px>

[![PyPI version shields.io](https://img.shields.io/pypi/v/epitopepredict.svg)](https://pypi.python.org/pypi/epitopepredict/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/epitopepredict/badge/?version=latest)](https://epitopepredict.readthedocs.io/en/latest/?badge=latest)

### Background

**epitopepredict** provides a standardized programmatic interface and command line tool for executing multiple epitope prediction methods. Currently this largely consists of interfaces to several MHC binding prediction, the results of which can then be processed and visualized in a consistent manner. There is a built-in method for MHC-class I prediction and the TEPITOPEPan method is provided as a 'built in' method for MHC-class II. The IEDB tools and netMHCpan, netMHCIIpan and MHCFlurry are also supported. Those tools are free for academic use but have to be installed separately. This software runs on most linux systems. 

Documentation is at http://epitopepredict.readthedocs.io

### Installation

`pip install epitopepredict`

or for latest version on github:

`pip install -e git+https://github.com/dmnfarrell/epitopepredict.git#egg=epitopepredict`
Code Examples
=============

This page is for those using the Python API. For those wanting to use the command line application see the
Command line interface page. General usage of this package is to provide convenient access to binding prediction methods
and perform analysis on the results. There are multiple potential applications.

Methodology
-----------

MHC binding and other prediction methods are implemented by inheriting from a `Predictor` object. All such classes should at minimum override the `predict` method for scoring a single sequence. This may wrap methods from other python packages or call command line predictors. For example the `TepitopePredictor` uses the `epitopepredict.tepitope` module provided with this package.

The predict method should return a Pandas DataFrame. The `predict_sequences` method is used for multiple protein sequences contained in a dataframe of sequences in a standard format. This is created from a genbank or fasta file (see examples below). For large numbers of sequences `predict_sequences` should be called with save=True so that the results are saved as each protein is completed to avoid memory issues, since many alleles might be called for each protein. Results are saved with one file per protein/sequence in csv format.

The results are of the following form and are returned sorted by the score column::

        peptide       core      pos  score      name         allele  rank
   198  VIFRLMRTNFL  FRLMRTNFL  198    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   199  IFRLMRTNFLI  FRLMRTNFL  199    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   200  FRLMRTNFLIK  FRLMRTNFL  200    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   709  NRFVTLDGQQF  FVTLDGQQF  709    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   710  RFVTLDGQQFY  FVTLDGQQF  710    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   711  FVTLDGQQFYW  FVTLDGQQF  711    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   70   DSFLLMLCLHH  FLLMLCLHH   70    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   71   SFLLMLCLHHA  FLLMLCLHH   71    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   72   FLLMLCLHHAY  FLLMLCLHH   72    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   32   QGIVRQRVIPV  IVRQRVIPV   32    1.7  ZEBOVgp1  HLA-DRB1*0101    10

where name is the protein identifier from the input file (a locus tag for example) and a score column which will differ between methods. MHC-II methods can be run for varying lengths, with the core usually being the highest scoring in that peptide/n-mer (but not always).

Basics
------

imports::

    import epitopepredict as ep
    from epitopepredict import base, sequtils, analysis, plotting

create a Predictor object::

    #get list of predictors
    print base.predictors
    ['tepitope', 'netmhciipan', 'iedbmhc1', 'iedbmhc2', 'mhcflurry', 'mhcnuggets', 'iedbbcell']
    p = base.get_predictor('tepitope')

get sequence data::

    #get data in genbank format into a dataframe
    df = sequtils.genbank2Dataframe(genbankfile, cds=True)
    #get sequences from fasta file
    df = sequtils.fasta2Dataframe(fastafile)

run predictions for a protein sequence::

    seq = ep.testsequence
    label = 'myprot' #optional label for your sequence
    p = base.get_predictor('tepitope')
    p.predict(sequence=seq, allele='HLA-DRB1*01:01', length=11, name=label)

run predictions for multiple proteins::

    #run for 2 alleles and save results to savepath
    alleles = ["HLA-DRB1*01:01", "HLA-DRB1*03:05"]
    p = base.get_predictor('tepitope')
    p.predict_proteins(df, length=11, alleles=alleles, save=True, path=savepath)

run predictions for a list of peptides::

    from epitopepredict import peptutils
    seqs = peptutils.create_random_sequences(5000)
    p = ep.get_predictor('tepitope')
    x = p.predict_peptides(seqs, alleles=alleles)

run with multiple threads::

    x = p.predict_peptides(seqs, alleles=alleles, threads=4)

load previous results into a predictor::

    p.load(path=path) #where path stores csv files for multiple proteins
    p.load(filename=file) # where file is a csv formatted file of prediction results (can be 1 or more proteins)

Analysis
--------

get all the binders using the current data loaded into the predictor::

    #default is to use percentile cutoff per allele, returns a dataframe
    p.get_binders(cutoff=.95)

get binders for only one protein by top median rank::

    p.get_binders(name=name, cutoff=10, cutoff_method='rank')

get all promiscuous binders, returns a dataframe::

    pb = p.promiscuous_binders(n=2, cutoff=.95)
    #same using score cutoff
    pb = p.promiscuous_binders(n=2, cutoff_method='score', cutoff=500)

find clusters of binders in these results::

    cl = analysis.find_clusters(b, method, dist=9, minsize=3)
Introduction
============

epitopepredict provides a standardized programmatic and command line interface for executing multiple MHC binding prediction methods.
The results from each method can then be processed and visualized in a consistent manner.

Installation
============

This software should be run on a Linux operating system. Ubuntu is recommended but most major distributions will be fine. Windows is not supported. macOS (OS X) may work but has not been tested. If you use windows or a mac you can simply install a linux virtual machine and run from there. You can then run the command line interface or use in python. Install with pip using::

    pip install epitopepredict

Or for latest version on github::

    pip install -e git+https://github.com/dmnfarrell/epitopepredict.git#egg=epitopepredict

**Python dependencies**

* numpy
* pandas
* matplotlib
* biopython
* tornado
* bokeh
* wtforms

Prediction algorithms
---------------------

There are now multiple MHC binding prediction algorithms available freely online. Often the problem is determining how to use them and which alleles they support. The 'state of the art' algorithms are probably those based on neural networks such as netMHC class I and II routines. These are packaged as external tools and can be installed freely on your system.

**Supported algorithms**

+---------------------+-------------------------------------------------------------+
| name                | description                                                 |
+=====================+=============================================================+
| basicmhc1           | built-in MHC-class I predictor                              |
+---------------------+-------------------------------------------------------------+
| tepitope            | implements the TEPITOPEPan method, built in (MHC-II)        |
+---------------------+-------------------------------------------------------------+
| netMHCpan           | http://www.cbs.dtu.dk/services/NetMHCpan/  (MHC-I)          |
+---------------------+-------------------------------------------------------------+
| netMHCIIpan         | http://www.cbs.dtu.dk/services/NetMHCIIpan/ (MHC-II)        |
+---------------------+-------------------------------------------------------------+
| mhcflurry           | https://github.com/openvax/mhcflurry (MHC-I)                |
+---------------------+-------------------------------------------------------------+
| IEDB MHC-I tools    | http://tools.immuneepitope.org/mhci/download/               |
+---------------------+-------------------------------------------------------------+

MHCFlurry can be easily installed easily with pip.

Installing netMHCpan and netMHCIIpan
------------------------------------

Due to license restrictions these programs must be installed separately. You can go to these respective links to fill in forms that will give you access to the install file:

 * http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan
 * http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan

The install instructions can then be found in the readme files when you untar the downloaded file e.g. netMHCpan-4.0.readme. Remember to test the software is working before you use it in epitopepredict.

Installing IEDB MHC-I tools
---------------------------

Note that if using the netMHCpan programs above you probably **DO NOT** need to use the IEDB tools unless you have specific requirements to do so. The distributions 'IEDB_MHC*.tar.gz' contain a collection of peptide binding prediction tools for Major Histocompatibility Complex (MHC) class I and II molecules. The collection is a mixture of pythons scripts and linux 32-bit environment specific binaries. Linux environment is required. Under ubuntu you should also install tcsh and gawk::

    sudo apt install tcsh gawk

Download from http://tools.iedb.org/mhci/download/. Unpack the tar.gz files. Run the 'configure' script to set up path variables for trained models. This has been tested to work with version 2.17.

::

    tar -zxvf IEDB_MHC_I-*.*.*.tar.gz
    cd mhc_i
    ./configure.py

*MHC-II tools are not currently supported.*

Submit Bugs
===========

This software is under active development particularly with a view to improve the command line tools. Please use the github project page to submit bugs or suggestions: http://dmnfarrell.github.io/epitopepredict

References
==========

* Zhang, L., Chen, Y., Wong, H.-S., Zhou, S., Mamitsuka, H., & Zhu, S. (2012). TEPITOPEpan: extending TEPITOPE for peptide binding prediction covering over 700 HLA-DR molecules. PloS One, 7(2), e30483. http://doi.org/10.1371/journal.pone.0030483

* Nielsen, M., Lund, O., Buus, S., & Lundegaard, C. (2010). MHC class II epitope predictive algorithms. Immunology, 130(3), 319–28. http://doi.org/10.1111/j.1365-2567.2010.03268.x

* Karosiene, E., Rasmussen, M., Blicher, T., Lund, O., Buus, S., & Nielsen, M. (2013). NetMHCIIpan-3.0, a common pan-specific MHC class II prediction method including all three human MHC class II isotypes, HLA-DR, HLA-DP and HLA-DQ. Immunogenetics, 65(10), 711–24. http://doi.org/10.1007/s00251-013-0720-y

* Chaves, F. a, Lee, A. H., Nayak, J. L., Richards, K. a, & Sant, A. J. (2012). The utility and limitations of current Web-available algorithms to predict peptides recognized by CD4 T cells in response to pathogen infection. Journal of Immunology (Baltimore, Md. : 1950), 188(9), 4235–48. http://doi.org/10.4049/jimmunol.1103640
MHC Allele Nomenclature
=======================

Human or animal MHC allele names have a unique number corresponding to up to four sets of digits separated by colons. The length of the allele designation is dependent on the sequence of the allele and that of its nearest relative. All alleles receive at least a four digit name, which corresponds to the first two sets of digits, longer names are only assigned when necessary. See the full explanation at http://hla.alleles.org/nomenclature/naming.html.

For epitopemap, you can see in the job submission form that we only use the four digit names e.g. HLA-A*68:02. The binding prediction algorithms are not designed to distinguish on a finer level. The four digit names are usually enough for most purposes.

**Examples**

MHC-I allele (using the IEDB tools)::

    HLA-A*68:02

For MHC-II alleles we use the following format::

    HLA-DRB1*01:01

Some methods will work if you leave out the colon separator but not all, so it's best to use the standard naming scheme.

**References**

* http://hla.alleles.org/nomenclature
* https://www.ebi.ac.uk/ipd/mhc/Command Line Interface
======================

Installing the package provides the command `epitopepredict` in your path. This is a command line interface
to the library without the need for any Python coding. It provides pre-defined functionality with settings
specified in a text configuration file. Using this you can make MHC predictions with your chosen alleles and
predictors. If you are using the IEDB prediction tools they should be installed locally and you can specify
the path in the [iedbtools] section. Otherwise ignore those settings. Note that if settings are left out
generally defaults will be used so you can have a minimal file as in the examples.

Usage
-----

Usage largely involves setting up the config file and having your input files prepared.
Running the command `epitopepredict -c <yourfilename>.conf` will create a new config file for you to work from if it doesn't exist.
Just edit this with a text editor and then to execute::

    epitopepredict -c <yourfilename>.conf -r

You can also test the pipeline after installing by running::

    epitopepredict -t

This will generate predictions using a set of sample HIV-1 sequences and save the results to a folder called hiv1_test which you can open in the web app to view (see below). This should work 'out of the box' as it only uses the built in prediction algorithm, tepitope.

Configuration file settings
---------------------------

The advantage of configuration files is in avoiding long commands that have to be remembered or are prone to mistakes. Also the config files can be kept to recall what setting we used or to copy them for another set of files. The current options available in the file are shown below::

    [base]
    predictors = tepitope
    mhc2_alleles = HLA-DRB1*01:01,HLA-DRB1*04:01
    mhc1_alleles = HLA-A*01:01
    mhc1_length = 11
    mhc2_length = 15
    n = 2
    cutoff_method = default
    cutoffs = .95
    sequence_file =
    path = results
    overwrite = no
    verbose = no
    names =
    cpus = 1

    [iedbtools]
    iedbmhc1_path =
    iedb_mhc1_method = IEDB_recommended

Settings explained
------------------

+------------------+-----------------------------+------------------------------------------------------------------------------+
| name             | example value               | meaning                                                                      |
+==================+=============================+==============================================================================+
| predictors       | tepitope                    | name of predictor: e.g. tepitope, iedbmhc1, netmhciipan, mhcflurry           |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc1_alleles     | HLA-A*01:01,HLA-A*03:01     | list of MHC-I alleles or preset name                                         |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc2_alleles     | HLA-DRB1*0101,HLA-DRB1*0103 | list of MHC-II alleles or preset name                                        |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc1_length      | 11                          | length of n-mers for MHC-I prediction                                        |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc2_length      | 11                          | length of n-mers for MHC-II prediction                                       |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| n                | 3                           | minimum number of alleles for counting as promiscuous binders                |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cutoff_method    | score                       | cutoff method: default, score or rank used for getting binders (see below)   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cutoffs          | .95                         | percentile/score/rank cutoff for counting binders                            |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| sequence_file    | zaire-ebolavirus.gb         | set of protein sequences in genbank or fasta format                          |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| peptide_file     | peptides.txt                | set of peptides in a plain text file, one per row                            |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| path             | results                     | folder to save results to, can be empty for current folder                   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| overwrite        | no                          | overwrite the previous results                                               |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| names            | Rv0011c,Rv0019c             | subset of protein/sequence names to predict from your input file, optional   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| verbose          | no                          | displays more information while running                                      |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cpus             | 1                           | number of processors to use, use 0 for all available                         |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedbmhc1_path    |                             | folder where the IEDB MHC-I tools are installed, not required unless used    |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedb_mhc1_method | IEDB_recommended            | predictor to use within the IEDB MHC-I tools (see below)                     |
+------------------+-----------------------------+------------------------------------------------------------------------------+

Cutoff methods
--------------

Methods for achieving an appropriate cutoff for considering a peptide to be a binder are somewhat arbitrary. They vary with the application. There are three methods provided to select binders:

* **default** - allele specific global cutoffs, this uses a percentile cutoff to select peptides using pre-calculated quantile scores for each allele. This may avoid an issue where certain alleles will dominate if using a single score cutoff. Though there is limited evidence to suggest this is more appropriate. Typical value would be .95 i.e. top 95% in each allele.
* **rank** - Select top ranking peptides in each sequence above the cutoff. This would be useful for small numbers of sequence but for a lot of proteins might produce too many false positives.
* **score** - Use a single score cutoff for all peptides/alleles. This is probably the standard method. Typical binding predictors produce an affinity score and a cutoff of 500 is used. However this might also produce a lot of false positives.

Binding promiscuity
-------------------

Promiscuous binders are those above the cutoffs in more than n alleles. The rationale for this is that a peptide is more likely to be immunogenic in your target population if it is a binder in multiple alleles. This may not be the case in reality of course. By default the command line tool will calculate the promiscuous binders to give you a unique list of peptides and include the number of alleles in which it is a binder. The table is ranked by this value and the maximum score over the alleles tested.

Preset allele lists
-------------------

For convenience there are some lists of common alleles that you can use without having to type allele names into the config file. These have been taken from various sources and are only a rough guide. Use `epitopepredict -p` to see the available presets. The format of allele names is discussed on the MHC Allele Nomenclature page.

The current selection is:

+---------------------+--------------------------------------------------------+
| name                | description                                            |
+---------------------+--------------------------------------------------------+
| mhc1_supertypes     | 6 MHC-I supertypes                                     |
+---------------------+--------------------------------------------------------+
| mhc2_supertypes     | 7 MHC-II supertypes                                    |
+---------------------+--------------------------------------------------------+
| us_caucasion_mhc1   | 30 most common US caucasion MHC-I                      |
+---------------------+--------------------------------------------------------+
| us_african_mhc1     | 30 most common US african MHC-I                        |
+---------------------+--------------------------------------------------------+
| human_common_mhc2   | 11 most prevalent HLA-DR alleles worldwide             |
+---------------------+--------------------------------------------------------+
| broad_coverage_mhc1 | 26 alleles providing broad coverage                    |
+---------------------+--------------------------------------------------------+
| bovine_like_mhc2    | 8 HLA-DR alleles chosen to approximate bovine response |
+---------------------+--------------------------------------------------------+

IEDB tool methods
-----------------

The IEDB combines multiple prediction methods into its tools. Generally it's recommended to use their IEDB_recommended method but individual methods may be preferred. You can specify these using the iedb_mhc1_method option. Remember they do not all support all alleles. See Installing IEDB MHC-I tools.

::

    ann
    comblib_sidney2008
    consensus
    IEDB_recommended
    netmhcpan
    smm
    smmpmbec

Examples
--------

**MHC-II binding predictions for preset alleles of proteins in a genbank file**

Using preset allele lists saves you the trouble of writing the alleles out. You can get the built-in presets by using -p at the command line. If you provide MHC-I alleles for a class II predictor like tepitope the program will give an error. More cpus means speed improvements::

    [base]
    predictors = tepitope
    mhc2_alleles = human_common_mhc2
    n = 2
    cutoffs = .95
    sequence_file = zaire-ebolavirus.gb
    path = results
    cpus = 2

**A small set of peptides**

Say we want to predict for small list of peptides with multiple prediction methods and select the top 10 ranking in at least 3 alleles. Here input.txt is just simple text file with all the individual peptides. They should be of an appropriate length::

    [base]
    predictors = tepitope,mhcflurry
    mhc1_alleles = human_common_mhc2
    mhc2_alleles = human_common_mhc2
    cutoff_method = rank
    cutoffs = 10
    n=3
    path = results
    peptide_file = input.txt

**Strict cutoffs**

For selection you can use very strict score cutoff level or high global percentile. In this example we use a score cutoff so must provide a cutoff value for each method::

    [base]
    predictors = tepitope,netmhciipan
    mhc1_alleles = human_common_mhc2
    cutoff_method = score
    cutoffs = 6,50
    n=3
    path = results
    peptide_file = input.txt

Outputs
-------

In each results folder you will find a sub-folder for each method. This has csv files with the predictions for each sequence, if using multiple protein sequences. This is the primary raw output. These folders can be re-used as input in the analysis section without re-running predictions and read by the web interface for presentation if needed. There are also files of the form final_method_n.csv which contain the promiscuous binders for each method.
File Formats
============

The command line interface currently accepts protein sequences as fasta or genbank files. Users will probably be familiar with these formats anyway if they are using them. A few useful tips are provided here:

* try to use sensible names in fasta file headers as they are used as identifiers for the results. When you get fasta files from some sources they can have very long headers like this::

    >lcl|NC_001802.1_prot_NP_057849.4_1 [gene=gag-pol] [locus_tag=HIV1gp1] [db_xref=GeneID:155348] [protein=Gag-Pol] [exception=ribosomal slippage] [protein_id=NP_057849.4] [location=join(336..1637,1637..4642)] [gbkey=CDS]


By default the text before the first space is used as the identifier for each protein so that should be unique. In this case it will be `lcl|NC_001802.1_prot_NP_057849.4_1`. You can also include an option called `fasta_header_sep` in the configuration file that will split the fasta name with another symbol as well, in this way you can shorten the names further, but they should still be unique.

* only CDS (coding sequence) features are used from genbank files, as obviously these are the ones with sequences

* make sure the /translation qualifier is present in the features of the genbank file. Some files might not have it and therefore no sequence is present. A typical genbank feature looks like this::

     CDS             360172..360507
                     /locus_tag="lmo0332"
                     /experiment="EXISTENCE:[PMID:19448609]"
                     /note="lmo0332"
                     /codon_start=1
                     /transl_table=11
                     /product="hypothetical protein"
                     /protein_id="NP_463862.1"
                     /db_xref="GeneID:987567"
                     /translation="MIYYICALYTFISALVSFGFSLDALLKSRKVNGDALINAKYAVS
                     RSLSLLIVALGLFIFKSDAFLVALSLVMIGAQLFDGIIGIKISTFKTVGPLLTAVGNV
                     IMLILFLTI"
Neoepitope Prediction
=====================

About
-----

It is now known that tumors elicit adaptive immune responses and that the antigens driving effective T cell response are generated from somatically mutated genes. Cancer vaccines and adoptive T cell therapy therefore depends on identification of these patient-specific potential neo-epitopes that might be targeted. This is achieved by applying whole exome sequencing of matched cancer and normal tissues, usually along with RNA-seq quantification of tumour gene expression. Variant calling is then performed, the new mutated peptides extracted and potential epitopes predicted. The computational aspect is broadly outlined below. There are several other software pipelines designed for this task, notably pvactool_ and MuPeXI. Users are encouraged to try at least one of these also.

.. _pvactool: http://pvactools.readthedocs.io/en/latest/pvacseq/

.. image:: neoepitope_workflow.png

Method
------

The program currently accepts vcf or maf files that have been output from a variant calling program. In future this step will be added so that the user can provide raw reads. The vcf/maf file is processed using the `varcode` python library for variant effect prediction to estimate potential mutated coding sequence. The resulting mutated sequence regions are broken up into peptides which are then filtered for similarity to the self proteome. MHC binding predictions are then performed on the remainder.

Usage and Configuration
------------------------

The pipeline is run via the same command line tool using `epitopepredict`. This uses the same text confguration file to provide the inputs and settings. Settings specific to neoepitope prediction are in the [neoepitope] section. These are explained below. Running the tool is then as simple as calling this command::

    epitopepredict -c <yourfilename>.conf -n

You can test the neoepitope pipeline after installing by running::

    epitopepredict -n -t

This will check that the human reference genome is available. (Note for users running the snap package: these files will be placed in your home directory usually under */home/user/snap/epitopepredict/x1/.cache/pyensembl/*.  If you uninstall the package you can also delete this folder to clear space).

neopredict accepts one or more vcf or maf files which have been created from a variant calling program. If more than one file, they should be comma separated. These file names specified in the [neoepitope] of the configuration file::

    [neoepitope]
    vcf_files =
    maf_files =

The remaining options not in this section are in the other sections covered in the command line interface section of the docs. For example you can specify the prediction algorithm, alleles and length of peptides in those sections.

References
----------

* Y. C. Lu and P. F. Robbins, “Cancer immunotherapy targeting neoantigens,” Semin. Immunol., vol. 28, no. 1, pp. 22–27, 2016.
* M. Efremova, F. Finotello, D. Rieder, and Z. Trajanoski, “Neoantigens generated by individual mutations and their role in cancer immunity and immunotherapy,” Front. Immunol., vol. 8, no. November, pp. 1–8, 2017.
* N. P. Restifo, M. E. Dudley, and S. A. Rosenberg, “Adoptive immunotherapy for cancer: Harnessing the T cell response,” Nat. Rev. Immunol., vol. 12, no. 4, pp. 269–281, 2012.
* A. M. Bjerregaard, M. Nielsen, S. R. Hadrup, Z. Szallasi, and A. C. Eklund, “MuPeXI: prediction of neo-epitopes from tumor sequencing data,” Cancer Immunol. Immunother., vol. 66, no. 9, pp. 1123–1130, 2017.
* J. Hundal et al., “pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens,” Genome Med., vol. 8, no. 1, p. 11, 2016.
* A. Rubinsteyn et al., “Computational pipeline for the PGV-001 neoantigen vaccine trial,” Front. Immunol., vol. 8, no. January, p. 1807, 2017.epitopepredict
==============

.. toctree::
   :maxdepth: 4

   epitopepredict
Web Application
===============

A web app that is launched from the command line can be used to view and analyse results from a set of predictions that you have made. This is an improved and much easier to use form of a previous web interface called epitopemap and replaces it. Note: this app is still under development, suggestions for additional functionality are welcome. In the future it will probably be possible to launch predictions inside the app. For now you run predictions on the command line and view them in the browser.

Usage and Interface
-------------------

After you have made some predictions you can run the app (usually from the same folder where you ran your predictions) using::

    epitopepredict -s

The default port is 8888. You can use a different port by specifying it with the -x option. This should open a new web page in your browser. To view all results from a run you then enter the path (folder) where your results were saved in the form and press submit. This should refresh your form with a drop down list of all the available sequences/proteins.

There are several ways to view a set of binding predictions, all of which allow views for whichever predictors you have used. There is currently

**Summary view**

Summarizes the results for multiple sequences in one page. You can choose to view the table of all predicted binders, promiscuous binders or a summary over each sequence. These tables can be downloaded to csv files.

**Sequence view**

For viewing the detailed results for a single sequence, often representing a protein coding sequence. Graphical views of the prediction scoring across the sequence are designed to provide a quick look at the pattern of peptide binding prediction in multiple alleles. By default a track view of each allele/predictor is shown as below:

.. image:: web_app_scr1.png

**Track plots** are useful for overall results over protein-length and longer sequences. The plot can be zoomed and panned using the mouse. A hover tooltip shows the particular peptide details.

**Grid plots** are another way to view peptide scores across a protein sequence for each allele. A hover tooltip shows the particular peptide details.

.. image:: web_app_grid_view.png

**Config**

Allows you to generate a configuration file from a form for running a set of predictions. In future this could be used to submit jobs directly.

Future features
---------------

* Improved graphical features for genome based prediction.
* Location of clusters of binders in sequences.
* Export of peptide lists/n-mers for experimental use.
* Mutation/conservation analysis.
* Edit config and run predictions from web page.
Welcome to the epitopepredict documentation.
============================================

Contents:

.. toctree::
   :maxdepth: 2

   description
   cli
   webapp
   examples
   nomenclature
   formats
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
epitopepredict package
======================

Submodules
----------

epitopepredict.analysis module
------------------------------

.. automodule:: epitopepredict.analysis
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.app module
-------------------------

.. automodule:: epitopepredict.app
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.base module
--------------------------

.. automodule:: epitopepredict.base
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.cluster module
-----------------------------

.. automodule:: epitopepredict.cluster
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.config module
----------------------------

.. automodule:: epitopepredict.config
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.neo module
-------------------------

.. automodule:: epitopepredict.neo
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.peptutils module
-------------------------------

.. automodule:: epitopepredict.peptutils
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.plotting module
------------------------------

.. automodule:: epitopepredict.plotting
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.sequtils module
------------------------------

.. automodule:: epitopepredict.sequtils
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.tepitope module
------------------------------

.. automodule:: epitopepredict.tepitope
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.tests module
---------------------------

.. automodule:: epitopepredict.tests
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.utilities module
-------------------------------

.. automodule:: epitopepredict.utilities
    :members:
    :undoc-members:
    :show-inheritance:

epitopepredict.web module
-------------------------

.. automodule:: epitopepredict.web
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: epitopepredict
    :members:
    :undoc-members:
    :show-inheritance:
