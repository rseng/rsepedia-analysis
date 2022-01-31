[![Codacy Badge](https://api.codacy.com/project/badge/Grade/887cb2a956394c9ea9b35855fc046663)](https://www.codacy.com/manual/marieBvr/virAnnot?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=marieBvr/virAnnot&amp;utm_campaign=Badge_Grade)

# Funding

<img src="docs/source/INRA_logo.jpg" width="30%"/> &emsp;<img src="docs/source/ubx-logo.png" width="30%"/>

# License

This work is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0).

# VirAnnot

VirAnnot was build to ease the assembly, blast search and taxonomic annotation of metagenomic NGS data. It is used in Virologie team of [UMR1332 BFP](http://www6.bordeaux-aquitaine.inra.fr/bfp) laboratory at INRA.
VirAnnot also takes part of the Euphresco project "[Plant Health Bioinformatics Network](https://doi.org/10.5281/zenodo.3245830)". See [more](https://gitlab.com/ahaegeman/phbn-wp2-training).

It is designed to identify viruses in plant metagenomic data but it can be used to assemble and annotate any sequences with the NCBI taxonomy.

Link to the [article](https://doi.org/10.1094/PBIOMES-07-19-0037-A)

# Documentation

See documentation:
<https://virannot-docs.readthedocs.io/en/latest/>

Pipeline general scheme:

![scheme](docs/source/dia-intro.png)
Example execution
=================

Create a directory for your experiment::

  mkdir test_virAnnot
  cd test_virAnnot

Copy example read files, Illumina adapters fasta file, the sample id mapping file, the step and parameter file::

  cp /path/to/virAnnot/examples/reads/*.fq .
  cp /path/to/virAnnot/examples/adapters.fa .
  cp /path/to/virAnnot/examples/map.txt .
  cp /path/to/virAnnot/examples/step.yaml .
  cp /path/to/virAnnot/examples/parameters.yaml .

You have to modify this file to fit your configuration.

This example contains all modules and options available and must be used as a template for your own analysis.

Step **ReadSoustraction**
*************************

This module uses bowtie2 to map reads against nucleotide sequence and Samtools to extract unmapped pairs.

Corresponding ``step.yaml`` section:

.. literalinclude:: ../../examples/step.yaml
  :lines: 1-9

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n ReadSoustraction_phiX

Step **Demultiplex**
********************

This module uses cutadapt demultiplex reads from library and also trim reads from adapters and chimeric reads.
Each demultiplexing step are described in the module section.
Corresponding ``step.yaml`` section:

.. literalinclude:: ../../examples/step.yaml
  :lines: 10-19

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Demultiplex

Step **DemultiplexHtml**
************************

This module gather all ``*_demultiplex.stats.csv`` files and create and html report.
Corresponding ``step.yaml`` section:

.. literalinclude:: ../../examples/step.yaml
  :lines: 20-25

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n DemultiplexHtml

Output example:

.. raw:: html

  <iframe src="_static/stat_demultiplex/index.html" height="400px" width="100%"></iframe>


Step **Diamond**
****************

This module launch Diamond similarity search for reads and produce an XML file per sample.

.. literalinclude:: ../../examples/step.yaml
  :lines: 63-72

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Diamond


Step **Assembly**
*****************

This module simply launch udba_ud, newbler and metaspades assemblers in each sample folder, rename scaffolds id and move the resulting fasta file.

.. literalinclude:: ../../examples/step.yaml
  :lines: 33-46

Test both idba and spades:

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Assembly_idba
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Assembly_spades


Example of idba assembly::

  >ds2015-149_0
  GTGTAAGGTGGTGAAGG...
  >ds2015-149_1
  CCTGCGAATTGGGCCAA...

Step **Map**
************

This module uses bowtie2 to map reads back on assemblies and samtools will count reads per contig.

Step file:

.. literalinclude:: ../../examples/step.yaml
  :lines: 47-62

Command:

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Map_idba
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Map_spades

Output a two column tabular file, column 1 sequence ID, column 2 number of reads.
Example of ``.rn`` file produce::

  ds2015-149_0    1179
  ds2015-149_1    444
  ds2015-149_10   26
  ds2015-149_11   44
  ds2015-149_12   14
  ds2015-149_13   4
  ds2015-149_14   6

Step **Blast**
**************

This module is able to launch Blast(s) against provided databases localy or remotely.
The script blast_launch.py must be present on distant servers and ``parameter.yaml`` modified to fit your servers.

Step file:

.. literalinclude:: ../../examples/step.yaml
  :lines: 94-136

Commands:

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast_nr
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast_refvirTX
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast_allvirTX
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast_RPS

Step **Blast2ecsv**
*******************

  This module uses the XML file generated by the corresponding Blast module and the taxonomy contained in the SQLITE database to create a csv file containing match options, taxonomy string and sequences.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast2ecsv_nr
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast2ecsv_refvirTX
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Blast2ecsv_allvirTX

Example output of ds2015-149_idba.scaffold.tbltx.all_vir.csv::

  #algo   query_id        nb_reads        query_length    accession       description     organism        percentIdentity nb_hsps queryOverlap    hitOverlap      evalue  score   tax_id  taxonomy        sequence
  "TBLASTX"       "ds2015-149_52" "6"     "117"   "KX274275.1"    "Grapevine rupestris stem pitting associated virus isolate SK704 B, complete genome"    "Grapevine rupestris stem pitting-associated virus"     "95.8"  "3"     "100"   "3"     "7.55823333338424e-05"  "222.2257"      "196400"        "Viruses;ssRNA viruses;Betaflexiviridae;Foveavirus;Grapevine rupestris stem pitting-associated virus"   "GAACACTATGAACGACAACTGGAAATCTGAGCACGCTATAAACACCCACAAACTCAAACGTAGACAAAGCTTTAACTAAGTTATTCATAATAATCACACCATGCCAAACACTTAAGG"

Step **Rps2ecsv**
*****************

  This module uses the rpstblastn XML file and the PFAM taxonomy to annotate query sequences and produce a readable CSV file.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Rps2ecsv

.. literalinclude:: ../../examples/step.yaml
  :lines: 184-188

Example output of ds2015-149_idba.scaffold.rps.pfam.csv::

  #query_id       query_length    cdd_id  hit_id  evalue  startQ  endQ    frame   description     superkingdom    no rank family  genus
  "ds2015-149_0"  "1428"  "pfam01443"     "gnl|CDD|279750"        "1.33194e-06"   "29"    "223"   "2"     "pfam01443, Viral_helicase1, Viral (Superfamily 1) RNA helicase.  Helicase activity for this family has been demonstrated and NTPase activity. This helicase has multiple roles at different stages of viral RNA replication, as dissected by mutational analysis."     "Viruses(1.00);"        "ssRNA viruses(0.99);unclassified viruses(0.01);"       "Alphaflexiviridae(0.36);Virgaviridae(0.24);Betaflexiviridae(0.15);Tymoviridae(0.10);Bromoviridae(0.07);"       "Potexvirus(0.26);Allexivirus(0.10);Carlavirus(0.08);Tymovirus(0.08);Tobamovirus(0.08);"


Step **Ecsv2excel**
*******************

  This module takes all csv files and create a compile them in a single Excel file.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Ecsv2excel

.. literalinclude:: ../../examples/step.yaml
  :lines: 189-195

Example output of ds2015-149_idba.scaffold.xlsx:

.. image:: ecsv2excel.png

Step **Ecsv2krona**
*******************

  This module uses CSV files from Blast2ecsv module to create Krona html file.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Ecsv2krona
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Ecsv2krona_dmd

Example output of Krona:

.. raw:: html

  <iframe src="_static/krona_blast/blast.global.krona.html" height="400px" width="100%"></iframe>
  <iframe src="_static/krona_diamond/global_krona_dmd.html" height="400px" width="100%"></iframe>


Step **Automapper**
*******************

  This module uses Blast CSV annotation file to select reference sequences, map query sequences and produce png of identity plot and alignment file in fasta format.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Automapper_nr
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Automapper_allvirTX
  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Automapper_refseqTX

Example output of ds2015-149/ds2015-149_autoMapper_nr:

.. raw:: html

  <iframe src="_static/results/index.html" height="400px" width="100%"></iframe>


Step **Rps2tree**
*****************

  This module use Rps2ecsv results of all sample to cut and group sequences based on identified domains and create OTUs, identity matrix, tree nexus files and png for each domains colored by SampleID.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Rps2tree

.. raw:: html

  <iframe src="_static/rps2tree_global/index.html" height="400px" width="100%"></iframe>

Step **Getresults**
*******************

  This module simply copy important results file to a result directory.

.. code-block:: bash

  virAnnot.py -m map.txt -s step.yaml -p parameters.yaml -n Getresults
Prerequisite
============
External programs
-----------------
* NCBI Blast+ suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST)
* SQLite (https://www.sqlite.org/)
* Mummer3 (http://mummer.sourceforge.net/)
* Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Cutadapt (https://github.com/marcelm/cutadapt)
* ETE tree (http://etetoolkit.org/)
* IDBA-UD (https://github.com/loneknightpy/idba)
* drVM (https://sourceforge.net/projects/sb2nhri/files/drVM/)
* Open Grid Scheduler (http://gridscheduler.sourceforge.net/)
* Diamond (https://github.com/bbuchfink/diamond)
* SortMeRNA (https://bioinfo.lifl.fr/RNA/sortmerna/)
* PrintSeq-lite (http://prinseq.sourceforge.net/)
* Samtools (http://samtools.sourceforge.net/)
* Bcftools (https://samtools.github.io/bcftools/)
* Seqtk (https://github.com/lh3/seqtk)
* NCBI utils (https://www.ncbi.nlm.nih.gov/books/NBK179288/)


External databases
------------------
* NCBI nr, nt (ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
* NCBI Taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy)
* PFAM (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) (rpsblast files, fasta files, and smp files)


Perl external libraries
-----------------------
* Getopt::Long
* \File::Basename
* DBI
* \Data::Dumper
* Bioperl
* Color::Rgb
* List::Util
* Spreadsheet::WriteExcel
* Log::Log4perl
* DBD::SQLite
* SQL::SplitStatement
* Math::Round
* String::Random
* Bio::SearchIO:blastxml
* Bio::SeqIO


Perl included libraries
-----------------------
* Tools::Fasta
* Tools::Fastq
* Tools::Blast
* Tools::Taxonomy
* Logger::Logger


Python library
--------------
* os
* call
* logging
* random
* string
* argparse
* re
* sys
* Bio
* time
* glob
* shutil
* yaml
* csv
* importlib
* matplotlib


Install
=======
Add tools and launchers folders to your $PATH.

.. code-block:: bash

  export PATH=/path/to/tools:/path/to/launchers:$PATH

Add lib folder to your $PERL5LIB.

.. code-block:: bash

 export PERL5LIB=/path/to/lib:$PERL5LIB


Database
========

NCBI taxonomy and the homemade per domain Pfam taxonomy are stored in a simple SQLite database.

Schema:

.. image:: taxo_sql.svg


NCBI Taxonomy
-------------

- Download and extract NCBI taxonomy files.

.. code-block:: bash

 wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz ; gunzip taxdump.tar.gz; tar -xf taxdump.tar;
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ; gunzip prot.accession2taxid.gz;
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz ; gunzip nucl_gb.accession2taxid.gz;
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz ; gunzip dead_prot.accession2taxid.gz;
 cat prot.accession2taxid dead_prot.accession2taxid > acc2taxid.prot
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz ; gunzip nucl_wgs.accession2taxid.gz;
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz ; gunzip dead_wgs.accession2taxid.gz
 cat nucl_wgs.accession2taxid nucl_gb.accession2taxid dead_wgs.accession2taxid > acc2taxid.nucl
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz; gunzip dead_nucl.accession2taxid.gz;
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz; gunzip gi_taxid_prot.dmp.gz;

Optionally you can combine multiple accession2taxid file with a simple cat. But keep separated nucl and prot accessions as they will be loaded in two different tables.

Launch the loadTaxonomy.pl script that will create the sqlite database. The script needs two provided sqlite files: ``taxonomyIndex.sql`` and ``taxonomyStructure.sql`` that describe the database struture. All these files are in virAnnot/db/.

.. code-block:: bash

 ./loadTaxonomy.pl -struct taxonomyStructure.sql -index taxonomyIndex.sql -acc_prot acc2taxid.prot -acc_nucl acc2taxid.nucl -names names.dmp -nodes nodes.dmp -gi_prot gi_taxid_prot.dmp


PFAM taxonomy
-------------

The pipeline modules ``rps2ecsv`` and ``rps2tree`` need taxonomic information of the PFAM domains to work.
You need to extract these informations and load it into the sqlite database.
Be carefull, the files you will download have a size of ~900Mo.

- Download and extract Pfam FASTA files:

.. code-block:: bash

 ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz
 tar -xzf fasta.tar.gz;
 mkdir pfam
 mv pfam*.FASTA fasta/
 cd pfam/

- Extract taxonomic information for each sequence of each PFAM domain and store it in ``*.tax.txt`` files:

.. code-block:: bash

 ls -1 pfam*.FASTA | sed 's,^\(.*\)\.FASTA,./gi2taxonomy.pl -i & -o \1.tax.txt -db taxonomy.tmp.sqlite -r,' | bash

- Create a file of file for the ``*.tax.txt`` files:

.. code-block:: bash

 listPath.pl -d . | grep 'tax.txt' > idx

- Compute taxonomy statistic for each domain and create a sql file to load into the database:

.. code-block:: bash

 taxo_profile_to_sql.pl -i idx > taxo_profile.sql

- Load into the database:

.. code-block:: bash

 sqlite3 taxonomy.tmp.sqlite < taxo_profile.sql

- Modify path to the database by editing the following scripts:

.. code-block:: bash

 ./tools/rps2ecsv.pl:my $db = '/path/to/taxonomy.tmp.sqlite';
 ./tools/2krona_new.pl:my $db = '/path/to/taxonomy.tmp.sqlite';
 ./tools/ecsv2krona.pl:my $db = '/path/to/taxonomy.tmp.sqlite';
 ./tools/blast2ecsv.pl:my $db = '/path/to/taxonomy.tmp.sqlite';
 ./tools/rps2tree.pl:my $db = '/path/to/taxonomy.tmp.sqlite';
 ./tools/autoMapper.pl:  'taxonomyDatabase'  => '/path/to/taxonomy.tmp.sqlite'


NCBI Blast database
-------------------

NCBI non redundant databases are very large and similarity search using Blast is an intensive task. I recommand to use those databases on computer clusters.

- Download NCBI nr et nt Blast files.

.. code-block:: bash

  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz


Modify the parameters.yaml to fit your configuration.

::

  servers:
  genotoul:
    adress: 'genotoul.toulouse.inra.fr'
    username: 'stheil'
    db:
      nr: '/bank/blastdb/nr'
      nt: '/bank/blastdb/nt'


Reduced databases are a good choice for limited computer ressources and drastically faster similarity search. Here are some example commands using NCBI tools to download sequences.

- Reduced NCBI databases:

Get all viroids nucleotide sequence from genbank::

 esearch -db "nucleotide" -query "txid12884[Organism]" | efetch -format fasta > viroids_nucl.fna

Get all viruses nucleotide sequences from genbank::

 esearch -db "nucleotide" -query "txid10239[Organism]" | efetch -format fasta > viruses_nucl.fna

Create Blast DB example::

 makeblastdb -in viruses_nucl.fna -parse_seqids -dbtype nucl


- Download PFAM files for RPSBLAST.

.. code-block:: bash

  wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz
  wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz
  wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz

Here I use only PFAM domains but ``fasta.tar.gz`` and ``cdd.tar.gz`` contains files for the entire CDD database. You can either delete files that are not from PFAM database or use the complete CDD.

- Delete file that are not from PFAM:

.. code-block:: bash

  \ls -1 | grep -v 'pfam' | sed 's,^.*$,rm &,'

Add '| bash' if correct.

- Download entire CDD database:

.. code-block:: bash

  wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/CDD_LE.tar.gz
Parameters files
================

parameters.yaml
---------------
Defines paths for both local and remote binaries and databases. A template is provided in the ``examples`` directory.

.. literalinclude:: ../../examples/parameters.yaml


step.yaml
---------
Defines the steps that the pipeline will execute. A template is provided in the ``/examples`` directory.

Step names correspond to a python module that will launch the step. Step names are split based on the '_' character so you can launch multiple instance. For example you might want to launch blastx and blastn, so step names could be 'Blast_N' and 'Blast_X'. What is after the underscore do not matters, it is just used to differanciate the two steps.

Special words in bracket are used as substitution string.
- (file), (file1) and (file2)
- (SampleID)
- (library)

.. literalinclude:: ../../examples/step.yaml


.. _map.txt:

map.txt
-------
The map file describe the experiment. It is a tabulated file with the first line containing headers starting with '#'. It must contain at least two column: SampleID and file.
A template is provided in the ``examples`` directory.
This is a minimum map.txt file:

.. literalinclude:: ../../examples/map.txt

You can add categories for each sample so they can be used when coloring sequences in trees from the Rps2tree module.
One library can be attributed to multiple samples, as shown in the example. Thus the demultiplexing step will be able to differentiate each sample and separate them.
Modules description
===================
ReadSoustraction
----------------

ReadSoustraction module use bowtie2 to map reads against a reference. The output is directly piped to samtools for bam convertion and again piped to samtools to select unmapped paired reads (-f 12 -F 256). The pipe continues with bamtofastq tool to create two fastq files (r1 and r2).

.. code-block:: python3

  def _create_cmd (self):
    cmd = ''
    cmd += self.params['bin']['bowtie'] + ' -p ' + str(self.n_cpu)
    cmd += ' -x ' + self.db
    cmd += ' -1 ' + self.i1 + ' -2 ' + self.i2 + ' | '
    cmd += self.params['bin']['samtools'] + ' view -bS - '
    cmd += ' | ' + self.params['bin']['samtools'] + ' view -u -f 12 -F 256 - | ' + self.params['bin']['bedtools'] + ' bamtofastq -i - -fq ' + self.o1 + ' -fq2 ' + self.o2
    log.debug(cmd)
    return cmd

Options
*******
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``db``: Bowtie database.
- ``o1``: R1 output fastq file.
- ``o2``: R2 output fastq file.
- ``sge``: [BOOL]
- ``n_cpu``: [INT] Number of CPU to use.
- ``iter``: [SampleID]


Demultiplex
-----------

.. image:: demultiplex_0.png
    :scale: 15 %

Warning: This demultiplex procedure is specifc to our sequencing methods.

The demultiplex.pl script uses `cutadapt
<https://github.com/marcelm/cutadapt/>`_ to demultiplex and trim step by step sequences used by the sequencing technology and the dsRNA extraction protocol.

Our extraction and sequencing protocol:


Informations about indexes, common sequences, and sample ids are stored in the map.txt file which must be tab delimited file.
The script produce lots of temporary file (\*\_step\_\*) that can be deleted at the end of the execution.
Each step produce 3 type of files:

.info:
  the matching information.

.log:
  the execution information.

.out:
  the fastq file.

Options
*******
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``adapters``: Fasta file of adapters.
- ``middle``: Check for MIDs in middle of the reads and 1 : trim the reads or 2: exclude the read.
- ``min_qual``: Trim the read if the quality is below this threshold.
- ``polyA``: Trim poly-A tail.
- ``min_len``: Exclude reads under this size threshold.
- ``iter``: [SampleID]

Step 01: 5' index search
************************

.. image:: demultiplex_1.png
    :scale: 20 %



.. code-block:: perl

  _launchCutAdapt($self,$files->{1}, $self->{index}, $tmp_file_prefix . "_step.01_R1",'d','k','-g','0','1','0.8');

- ``'d','k'``: This step discard untrimmed reads, so reads that do not contain indexes are excluded, and keep trimmed reads.
- ``'-g'``: search indexes in 5'.
- ``'0'``: no errors allowed in the index sequence.
- ``'1'``: search for one index in each read.
- ``'0.8'``: 80% of the index length to be considered as a match.


Step 02: 5' common sequence
***************************

.. image:: demultiplex_2.png
    :scale: 20 %


.. code-block:: perl

  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.02_R1",'k','k','-g','0.1',scalar(keys(%{$self->{_common}})),'0.7');

- ``'k','k'``: Keep reads that contains or not the common part.
- ``'-g'``: search in 5' part.
- ``'0.1'``: 10% of sequencing errors.
- ``scalar(keys(%{$self->{_common}}))``: will search as many common part as provided.
- ``'0.7'``: 70% of the common part length to be considered as a match.


Step 03: 5' common sequence fragments
**************************************

.. image:: demultiplex_3.png
    :scale: 20 %


.. code-block:: perl

  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.021_R1",'k','k','-g','0.2',scalar(keys(%{$self->{_common}})),'0.5');

- ``'k','k'``: Keep reads that contains or not the common part.
- ``'-g'``: search in 5' part.
- ``'0.2'``: 20% of sequencing errors.
- ``scalar(keys(%{$self->{_common}}))``: will search as many common part as provided.
- ``'0.5'``: 50% of the common part length to be considered as a match.


Step 04: Trimming sequencing adapters
*************************************

.. image:: demultiplex_4.png
    :scale: 20 %


.. code-block:: perl

  _launchCutAdapt($self,$files->{1},$self->{illuminaAdapter},$tmp_file_prefix . "_step.03_R1",'k','k','-b','0.2',scalar(keys(%{$self->{illuminaAdapter}})),'0.6');

- ``'k','k'``: Keep reads that contains or not the common part.
- ``'-b'``: search adapters anywhere in the read.
- ``'0.2'``: 20% of sequencing errors.
- ``scalar(keys(%{$self->{illuminaAdapter}}))``: will search as many adapters as provided.
- ``'0.6'``: 60% of the adapters length to be considered as a match.

Step 05: Search for hybrid reads
********************************

.. image:: demultiplex_5.png
    :scale: 20 %



This step is really specific to our extraction method since very short DNA fragment can be link together during the aspecific adapters ligation step of the Illumina kits. This creating reads composed of two different PCR product. Thus our program search for index sequence in the middle of the read and trim it to keep the 5' part or exlude the read.
The research is done both on provided indexes sequences and reverse complement of thoses sequences.

``-middle [1|2]``   Search for common tag in the middle of the read. 1: trim the read. 2: exclude the read.


.. code-block:: perl

  _launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.04_R1",'k','k','-b','0.1','1','0.5');

- ``'k','k'``: Keep reads that contains or not the index , or ``'k','d'`` if the ``-middle`` option is provided.
- ``'-b'``: search adapters anywhere in the read.
- ``'0.1'``: 20% of sequencing errors.
- ``'1'``: search for one index in each read.
- ``'0.5'``: 50% of the adapters length to be considered as a match.

Step 06: Search for polyA (optional)
************************************

.. image:: demultiplex_6.png
    :scale: 15 %



In Illumina technology, if the sequencing matrix is too short compared to the sequencing length, the sequencing machine adds a bunch of A's and then radom sequence.

.. code-block:: perl

  _launchCutAdapt($self,$files->{1}, $h, $tmp_file_prefix . "_step.05_R1",'k','k','-a','0','1','0.8');

- ``'k','k'``: Keep reads that contains or not the index , or ``'k','d'`` if the ``-middle`` option is provided.
- ``'-a'``: search in the 3' end.
- ``'0'``: no sequencing errors.
- ``'1'``: search for one index in each read.
- ``'0.8'``: 80% of the polyA length to be considered as a match.

Assembly
--------

This module can launch three assemblers, `IDBA
<https://github.com/loneknightpy/idba/>`_, `MetaSpades
<http://bioinf.spbau.ru/en/spades3.7/>`_ and `Newbler 
<https://wikis.utexas.edu/display/bioiteam/GS+De+novo+assembler>`_ for single-end data.

Foreach assembler, the module convert reads files to the proper format, launch the assembly in a separate directory, rename scaffolds identifier and move results file to the sample root directory.

Options
*******


Map
---

This module uses bowtie2, samtools and readPerContig.pl script to map reads back on the assembly and count for each scaffold the number of reads aligned resulting a simple two column file scaffoldID and nb_reads used by other modules.

Options
*******
- ``contigs``: fasta file of contigs to map reads on. [mandatory]
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``ising``: singletons fastq file
- ``n_cpu``: [INT] number of CPU to use.
- ``sge``: [BOOL] use SGE scheduler.
- ``bam``: BAM file name.
- ``rn``: output file name.


Normalization
-------------
This module randomly select NUM reads from paired-files.

Options
*******
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``o1``: Output R1 normalized file. [mandatory]
- ``o2``: Output R2 normalized file. [mandatory]
- ``num``: [INT] Number of reads to randomly select. [mandatory]
- ``iter``: Iteration on [sample, library].
- ``n_cpu``: [INT] number of CPU to use.
- ``sge``: [BOOL] use SGE scheduler.


Diamond
-------

This module launches Diamond similarity search on reads and produce an XML file simalar to what Blast does so it can be treated by the Blast2ecsv module and so on.

Options
*******
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``db``: Values are defined in the parameters.yaml file. [mandatory]
- ``ising``: singletons fastq file
- ``n_cpu``: [INT] number of CPU to use.
- ``sge``: [BOOL] use SGE scheduler.
- ``sensitive``: [BOOL]
- ``more_sensitive``: [BOOL]
- ``out``: XML output file
- ``score``: Report matches above this score.
- ``max_target_seqs``: Maximum match per query sequences.
- ``evalue``: Min e-value.
- ``identity``: Report matches above this identity percent. 0 > X > 100.
- ``qov``: Query overlap.
- ``hov``: Hit overlap.

Diamond2Blast
-------------

Options
*******
- ``i``: CSV file with DIAMOND results. [mandatory]
- ``contigs``: Fasta file. [mandatory]
- ``out``: XML output file.
- ``type``: Blast type. ['tblastx','blastx','blastn','blastp','rpstblastn']. [mandatory]
- ``db``: Values are defined in the parameters.yaml file. [mandatory]
- ``evalue``: Min e-value.
- ``server``: ['enki','genologin','avakas'] Values are defined in the parameters.yaml file.
- ``n_cpu``: [INT] number of CPU to use.
- ``tc``: Number of task launched at the same time on SGE.
- ``num_chunk``: Number of chunks to split the original fasta file for parallel execution.
- ``max_target_seqs``: Maximum match per query sequences.
- ``sge``: [BOOL] use SGE scheduler.


Blast
-----

This module launches all type of Blast on local machine or distant servers. This module has been developped for our own local machines and servers, but it can be easly modified to fit your needs.

This module mainly depends on the parameters.yaml file and the blast_launch.py script which has to be present on the server you want to use and modified to fit your server configuration.

Options
*******
- ``contigs``: Fasta file. [mandatory]
- ``db``: Values are defined in the parameters.yaml file. [mandatory]
- ``type``: Blast type. ['tblastx','blastx','blastn','blastp','rpstblastn']. [mandatory]
- ``n_cpu``: [INT] number of CPU to use.
- ``tc``: Number of task launched at the same time on SGE. (Experimental, works on Genotoul)
- ``max_target_seqs``: Maximum match per query sequences.
- ``num_chunk``: Number of chunks to split the original fasta file for parallel execution.
- ``out``: Output file name.
- ``server``: ['enki','genologin','avakas', 'curta'] Values are defined in the parameters.yaml file.
- ``sge``: [BOOL] use SGE scheduler.

This module is able to launch Blast instance on distant servers if the database and the blast_launch.py script is present on the server. Then you have to edit the parameters.yaml file to fit your configuration. The script has been developped to use two computer cluster, Avakas (PBS + Torque) and Genotoul (SGE) but each cluster has its own configuration so you may have to modify this script to adapt it to your configuration.

Blast2ecsv
----------

This module parse Blast xml outputs, filter matches on different criteria and link Accession number to NCBI taxonomy.

Options
*******
- ``b``: Blast file.
- ``if``: Input format ['xml','m8']
- ``out``: Output file name.
- ``evalue``: Min e-value.
- ``fhit``: Only report first hit.
- ``fhsp``: Only report first hsp.
- ``pm``:
- ``r``: Reduced taxonomy. Report only 5 consistent rank.
- ``vs``: Only report sequences when match is virus or viroids.
- ``rn``: Read number. File created by the Map module.
- ``type``: Blast type. ['tblastx','blastx','blastn','blastp','rpstblastn']
- ``score``: Report matches above this score.
- ``identity``: Report matches above this identity percent. 0 > X > 100.
- ``qov``: Query overlap.
- ``hov``: Hit overlap.
- ``pd``: Parse description. Useful when the query ID is stored in the dscription field in the XML file.
- ``sge``: [BOOL] use SGE scheduler.

Ecsv2excel
----------

This module aggregates multiple ecsv file to create a colored XLSX file.
It launches the ecsv2krona.pl script.

Options
*******
- ``b[INT]``: CSV Blast file from 1 to 10.
- ``out``: Outpuy file name.
- ``r``: RPSBLAST csv file.
- ``sge``: [BOOL] use SGE scheduler.


Ecsv2krona
----------

This module launch the ecsv2krona.pl script. It will aggregate multiple ecsv file into one Krona html file.

Options
*******
- ``b``: [INT] CSV Blast file.
- ``id``: [INT] ID wanted corresponding to the Blast file.
- ``x``: [INT] XML Blast file. If used, this file will be split by species and link in the Krona file.
- ``out``: Output file name.
- ``data``: ['both','reads','contigs','none']
- ``r``: Use reduced taxonomy.
- ``c``: ['identity','taxid','none']
- ``iter``: ['global']
- ``sge``: [BOOL] use SGE scheduler.

Rps2ecsv
--------

This module launch rps2ecsv.pl script for each sample.
This module parse XML files from rpsblast and create csv file as a result.

Options
*******
- ``b``: RPSBLAST XML file.
- ``contigs``: Fasta file.
- ``sge``: [BOOL] use SGE scheduler.
- ``out``: Output file name.
- ``evalue``: e-value threshold.

Rps2tree
--------

This module launch rps2tree.pl script for all sample provided.
This module generates Operational Taxonomic Unit (OTU) based on RPS-Blast results.
For each CDD motifs, contigs are clustered together based on matrix distance.
The tree is generated thanks to `ete3 toolkit <http://etetoolkit.org/>`_.

Options
*******
- ``pfam``: CSV file from Rps2ecsv.
- ``contigs``: Fasta file.
- ``ecsv``: CSV file from Blast2ecsv.
- ``out``: Output file name.
- ``sge``: [BOOL]
- ``viral_portion``: Minimum percentage of viral sequence in a domain to be selected.
- ``min_prot``: Minimum protein length to be included in a tree.
- ``perc``: Percentage of identity. Threshold set to define OTU.
- ``iter``: ['global']

Rps2merge
---------
This module launches rps2merge script.
It generates a summary file, merging OTU results, blast results and rps results.
For each OTU, if multiple contigs corresponding, one is randomly selected.

Options
*******
- ``pfam``: CSV file from Rps2ecsv.
- ``blastx``: CSV file from Blast2ecsv.
- ``rps_folder``: Name of output folder of Rps2tree.
- ``id``: Sample or library ID
- ``out``: CSV output file
- ``iter``: ['global']
- ``sge``: [BOOL]

AutoMapper
----------

This module launch autoMapper.pl script for every sample.

Options
*******
- ``contigs``: Fasta file.
- ``ecsv``: CSV file from Blast2ecsv.
- ``i1``: R1 fastq file. [mandatory]
- ``i2``: R2 fastq file. [mandatory]
- ``out``: Output folder.
- ``sge``: [BOOL]
- ``ref``: Blast database name.


Blast2hits
----------

This takes multiple CSV Blast file from Blast2ecsv and draw histograms for by taxonomy.

Options
*******

- ``b[INT]`` CSV Blast file from Blast2ecsv
- ``id[INT]`` ID associated with the Blast file.
- ``iter`` [global]
- ``sge`` [BOOL]
- ``out`` Output file name.
.. image:: INRA_logo.jpg
  :align: left
  :scale: 3 %

.. image:: ubx-logo.png
  :align: right
  :scale: 3 %

|

Welcome to VirAnnot's documentation!
====================================

VirAnnot was build to ease the assembly, blast search and taxonomic annotation of metagenomic multi-sample datasets. It is used in the Virologie team of `UMR1332 BFP <http://www6.bordeaux-aquitaine.inra.fr/bfp>`_ laboratory at INRA.

It was designed to identify viruses in plants but it can be used to assemble and then annotate any sequences with the NCBI taxonomy.

NR and NT must be present localy and/or on distant servers and NCBI taxonomy is loaded in SQLITE database with a provided script.

Blast step is the most time consumming step and the use of large computer cluster is clearly an advantage.
Here we used two clusters :

`CURTA
<https://redmine.mcia.univ-bordeaux.fr/>`_ at Bordeaux University.

`GENOTOUL
<http://www.genotoul.fr//>`_ at Toulouse INRA datacenter.

Pipeline general scheme:
------------------------

.. image:: dia-intro.png


Guide
=====
.. toctree::
  :maxdepth: 2

  prerequisite
  parameter
  modules
  execution
