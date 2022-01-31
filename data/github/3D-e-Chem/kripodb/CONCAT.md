# Change log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/.

## Unreleased

## [3.0.0] - 2018-03-28

### Changed

* Fingerprints changed format from intbitset to Roaring bitmap (#52)

Fingerprint data files created with KripoDB < 3 will not work anymore.
They must be exported to makebits with old version of KripoDB and imported in new version.

### Fixed

* When generating similarities of the same fingerprints force ignore upper triangle

## [2.4.0] - 2018-03-28

### Added

- Commands to
    - Merge fingerprint/fragment/pharmacophore database files
    - Merge fragment database files
    - Convert phar formatted file from/to pharmacophore database file
- Methods required by kripo (https://github.com/3D-e-Chem/kripo)

### Fixed

- GET /fragments should accept PDB codes in uppercase (#51)
- Donor and acceptor features are swapped (#50)

## [2.3.1] - 2017-07-21

### Fixed

- clients for *.phar and *.svg endpoints have no response (#49) 

## [2.3.0] - 2017-05-17

### Added

- Similarities 
    - Histogram can output raw scores
    - Histogram can read frozen matrices using either lower or upper triangle
    - Export can be filtered by frag1 and/or pdb codes
    - Filter can be filtered by skip list or keep db
- Use scripts in update steps
- Pharmacophores
  - Store pharmocophore points in pytables table (#29)
  - Export pharmacophore points in [\*.phar format](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html#format)
  - Sub command to add points to table from a directory
  - Sub command to filter the pharmacophores points based on a fragments database
  - Sub command and webservice endpoint to fetch the pharmacophore points of a single fragment identifier (#30)
  - Canned method to fetch pharmacophores of a list fragment identifiers
- Dive
  - Tag pdb in file by filename
  - Scripts to run it

### Fixed

- Connexion internal change broke web service server 

## [2.2.1] - 2017-03-07

### Fixed

- Web service has internal server error when fragment has no molblock (#45)

## [2.2.0] - 2017-02-23

### Changed

- Canned methods can now raise exception with ids which could not be found and data for ids which could

### Fixed

- Fetch fragment with no molblock throws error (#41)
- Not found response of web service should be JSON (#42)

## [2.1.0] - 2017-01-17

### Added

- Webservice endpoint to render 2D fragment in SVG
- Published documentation on http://kripodb.readthedocs.io
- Documented update pipeline
- Documented command line interface in Sphinx
- Retrieve fragments from webservice based on fragment id and pdb code (#35)

### Changed

- merge similarity pairs files in chunks instead of loading whole source file in memory
- canned fragments_by_* methods can use local file or webservice
- error when duplicate fragment insert is performed
- Renamed `kripodb similarities serve` to `kripodb serve`, as it now also serves the fragments
- Switched from nosetest to py.test (#36)

### Removed

- no longer create indices for similarity pairs file, querying is done on dense matrix

## [2.0.0] - 2016-07-14

### Changed

- Renamed distance to similarity (#21)
- Flag to ignore upper triangle when calculating distances, instead of always ignore (#20)

## [1.4.2] - 2016-06-03

### Changed

- Lower webservice cutoff to 0.45 (#18)

## [1.4.1] - 2016-05-31

### Added

- Webservice online at http://3d-e-chem.vu-compmedchem.nl/kripodb/ui/
- Ignore_upper triangle option in distance import sub command

## [1.4.0] - 2016-05-03

### Changed

- Using nested sub-commands instead of long sub-command. For example `kripodb distmatrix_import` now is `kripodb distances import`

### Added

- Faster distance matrix storage format
- Python3 support (#12)
- Automated build to docker hub.

### Removed

- CLI argument `--precision`

## [1.3.0] - 2016-04-23

### Added

- webservice server/client for distance matrix (#16). The CLI and canned commands can now take a local file or a url.

### Fixed

- het_seq_nr contains non-numbers (#15)

## [1.2.5] - 2016-03-24

### Fixed

- fpneigh2tsv not available as sub command

## [1.2.4] - 2016-03-24

### Added

- Sub command to convert fpneight distance file to tsv.

## [1.2.3] - 2016-03-01

### Changed

- Converting distances matrix will load id2label lookup into memory to speed up conversion

## [1.2.2] - 2016-02-22

### Added

- Added sub command to read fpneigh formatted distance matrix file (#14)

## [1.2.1] - 2016-02-12

### Added

- Added sub commands to read/write distance matrix in tab delimited format (#13)
- Created repo for Knime example and plugin at https://github.com/3D-e-Chem/knime-kripodb (#8)

## [1.2.0] - 2016-02-11

### Added

- Prefix to canned fragments lookups (#11)
- PDB meta data to fragments db (#6)
- Limit to distance matrix searches (#5)

### Changed

- Merging of distance matrix files more robust (#10)
- Tanimoto coefficient is rounded up (#7)

## [1.0.0] - 2016-02-05

### Added

- Convert fragments shelve to sqlite
- Convert SDF molecules file to sqlite
- Convert Makebits formated file to sqlite
- Create distance matrix using modified tanimoto coefficient in hdf5 format
# Kripo DB

[![Build Status](https://travis-ci.org/3D-e-Chem/kripodb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/kripodb)
[![Build status](https://ci.appveyor.com/api/projects/status/diign2fenvai0dst?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/kripodb)
[![Codacy Grade Badge ](https://api.codacy.com/project/badge/Grade/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/3D-e-Chem/kripodb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/kripodb&amp;utm_campaign=Badge_Grade)
[![Codacy Coverage Badge](https://api.codacy.com/project/badge/Coverage/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/3D-e-Chem/kripodb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/kripodb&amp;utm_campaign=Badge_Coverage)
[![DOI](https://zenodo.org/badge/19641/3D-e-Chem/kripodb.svg)](https://zenodo.org/badge/latestdoi/19641/3D-e-Chem/kripodb)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/3dechem/kripodb/)
[![Documentation Status](https://readthedocs.org/projects/kripodb/badge/?version=latest)](http://kripodb.readthedocs.io/en/latest/?badge=latest)

Library to interact with Kripo fragment, fingerprint, pharmacophore and similarity data files.

Use [kripo](https://github.com/3D-e-Chem/kripo) to generate fragments, pharmacophores and fingerprints from PDB files.

KRIPO stands for Key Representation of Interaction in POckets, see [reference](http://dx.doi.org/10.1186/1758-2946-6-S1-O26) for more information.

# Glossary

* Pocket, binding site of the ligand in the protein of a crystal structure
* Fragment, part of the ligand
* Subpocket, part of the protein pocket which binds with the fragment
* Fingerprint, fingerprint of structure-based pharmacophore of subpocket
* Similarity matrix, similarities between all fingerprint pairs calculated using the modified tanimoto similarity index
* Kripo fragment identifier, used as identifier for fragment, subpocket and fingerprint

# Install

Requirements:

* RDKit, http://www.rdkit.org/docs/Install.html, to read SDF files and generate smile strings from molecules
* pip, version 8.0.0 or greater, for wheel support
* git, to clone kripodb repository during installation

```
pip install -U setuptools
pip install numpy
pip install git+https://github.com/3D-e-Chem/kripodb.git
```

# Usage

To see available commands
```
kripodb --help
```

## Create all

Commands to create all data files see [update documentation](docs/data-update.rst).

## Search for most similar fragments

Command to find fragments most similar to `3kxm_K74_frag1` fragment.
```
kripodb similar sim_all.h5 3kxm_K74_frag1 --cutoff 0.45
```

## Create similarity matrix from text files

Commands to create similarity matrix see [update documentation](docs/data-update.rst).

# Data sets

## Example

An example data set included in the [data/](data/) directory of this repo. See [data/README.md](data/README.md) for more information.

## GPCR

All fragments based on GPCR proteins compared with all proteins in PDB.

* kripo.gpcrandhits.sqlite - Fragments sqlite database
* kripo.gpcr.h5 - HDF5 file with similarity matrix

The data set has been published at [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.50835.svg)](http://dx.doi.org/10.5281/zenodo.50835)

## Protein Data Bank

All fragments form all proteins-ligand complexes in PDB compared with all.

* Fragments sqlite database - Download from http://3d-e-chem.vu-compmedchem.nl/kripodb/fragments.sqlite
* Pharmacophores database - Download from http://3d-e-chem.vu-compmedchem.nl/kripodb/pharmacophores.h5
* Similarity matrix - Can be queried on webservice at http://3d-e-chem.vu-compmedchem.nl/kripodb. For build instructions see http://kripodb.readthedocs.io/en/latest/data-update.html
* Fragment fingerprints - See http://kripodb.readthedocs.io/en/latest/data-update.html for instructions how to convert to a similarity matrix

Date at which the data of the 3d-e-chem.vu-compmedchem.nl webservice was last updated can found at http://3d-e-chem.vu-compmedchem.nl/kripodb/version.txt

A data set with PDB entries till 23 December 2015 has been published at [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.55254.svg)](http://dx.doi.org/10.5281/zenodo.55254)

# KNIME

The [Knime-KripoDB-example.zip](https://github.com/3D-e-Chem/knime-kripodb/blob/master/examples/Knime-KripoDB-example.zip) file is an example workflow showing how to use KripoDB python package inside KNIME (http://www.knime.org).
It can be run by importing it into KNIME.
Make sure the Python used by KNIME is the same as the Python with kripodb package installed.

The https://github.com/3D-e-Chem/knime-kripodb repo adds KripoDB code templates to KNIME.

# Development of KripoDB

Install the development deps with:
```
pip install -r requirements.txt
```

# Docker

## Create image

```
docker build -t 3dechem/kripodb .
```

## Run container

Show the kripodb help with
```
docker run --rm 3dechem/kripodb kripodb --help
```

To calculate the mean bit density of the fingerprints in the `fingerprints.sqlite` file in the current working directory use following command.
```
docker run --rm -u $UID -v $PWD:/data 3dechem/kripodb kripodb meanbitdensity /data/fingerprints.sqlite
```

# Web service

The Kripo data files can be queried using a web service.

Start webservice with:
```
kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5
```
It will print the urls for the swagger spec and UI.

Note! The webservice returns a limited amount of results. To get all results use local files.

On http://3d-e-chem.vu-compmedchem.nl/kripodb/ui/ there is a KripoDB webservice with the full PDB fragment all vs all matrix.
The date of the latest PDB record included in the webservice can be found in http://3d-e-chem.vu-compmedchem.nl/kripodb/version.txt

# Documentation

API and data update pipeline documentation can be found at http://kripodb.readthedocs.io/en/latest/.

# Reference

KRIPO – a structure-based pharmacophores approach explains polypharmacological effects;
Tina Ritschel, Tom JJ Schirris, and Frans GM Russel; J Cheminform. 2014; 6(Suppl 1): O26;
Published online 2014 Mar 11; http://dx.doi.org/10.1186/1758-2946-6-S1-O26
# Example data set

* fragments.sqlite - Fragments sqlite database containing a small number of fragments with their smiles string and molblock.
* fingerprints.sqlite - Fingerprints sqlite database with fingerprint stored as [fastdumped intbitset](http://intbitset.readthedocs.org/en/latest/index.html#intbitset.intbitset.fastdump)
* similarities.h5 - HDF5 file with similarities matrix of fingerprints using modified tanimoto similarity index 

## Creating tiny data set

1. Create fingerprints db with 1000 fingerprints
```
gunzip -c fingerprint01.fp.gz | head -1001 | kripodb fingerprints import - fingerprints.sqlite
```

2. Shrink fragments db to only contain fragments which have a fingerprint
```
cat | sqlite3 fragments.sqlite <<EOF
ATTACH DATABASE 'fingerprints.sqlite' AS fp;
DELETE FROM molecules WHERE frag_id NOT IN (SELECT frag_id FROM fp.bitsets);
DELETE FROM fragments WHERE frag_id NOT IN (SELECT frag_id FROM fp.bitsets);
DELETE FROM pdbs WHERE pdb_code || prot_chain NOT IN (SELECT pdb_code || prot_chain FROM fragments);
VACUUM;
EOF

```

3. Create similarity matrix

```
kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite fingerprints.sqlite fingerprints.sqlite similarities.h5
```

kripodb.script
==============

.. automodule:: kripodb.script
    :members:

.. automodule:: kripodb.script.fragments
    :members:

.. automodule:: kripodb.script.fingerprints
    :members:

.. automodule:: kripodb.script.similarities
    :members:

.. automodule:: kripodb.script.dive
    :members:
Data update
===========

The Kripo data can be updated in 2 ways:

.. toctree::
    :maxdepth: 1

    Baseline update <baseline-update.rst>
    Incremental update <incremental-update.rst>

Steps
-----

Overview of steps involved in updating Kripo:

1. Create staging directory
2. Create sub-pocket pharmacophore fingerprints
3. Create fragment information
4. Add new fragment information to fragment sqlite db
5. Populate PDB metadata in fragments database
6. Check no fragments are duplicated
7. Calculate similarity scores between fingerprints
8. Convert pairs file into dense similarity matrix
9. Switch staging to current
10. Update web service

.. note:: Steps 2 through 3 require undisclosed scripts or https://github.com/3D-e-Chem/kripo
.. note:: Steps 4 and 6 through 7 can be done using the KripoDB Python library.

.. todo:: Remove Kripo fragment/fingerprints of obsolete PDBs (ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat)

Disk layout
-----------

Directories for Kripo:

* **current/**, directory which holds current dataset
* **staging/**, which is used to compute new items and combine new and old items.
* **old/**, which is used as a backup containing the previous update.

Files and directories for a data set (inside current, staging and old directories):

* **pharmacophores.h5**, pharmacophores database file
* **out.fp.sqlite**, fingerprints file
* **fragments.sqlite**, fragment information database file
* **similarities.h5**, similarities as pairs table
* **similarities.packedfrozen.h5**, similarities as dense matrix

Input directories:

* **$PDBS_ADDED_DIR**, directory containing new PDB files to be processed

Requirements
------------

* Slurm batch scheduler
* KripoDB and it's dependencies installed and in path
* Posix filesystem, NFS of Virtualbox share do not accept writing of hdf5 or sqlite files
kripodb.webservice
==================

.. automodule:: kripodb.webservice
    :members:

.. automodule:: kripodb.webservice.client
    :members:

.. automodule:: kripodb.webservice.server
    :members:
kripodb.db
==========

.. automodule:: kripodb.db
    :members:
Command line interface
======================

.. argparse::
   :module: kripodb.script.__init__
   :func: make_parser
   :prog: kripodb

kripodb.hdf5
============

.. automodule:: kripodb.hdf5
    :members:
kripodb.pdb
===========

.. automodule:: kripodb.pdb
    :members:kripodb.canned
==============

.. automodule:: kripodb.canned
    :members:
kripodb.dive
============

.. automodule:: kripodb.dive
    :members:kripodb.makebits
================

.. automodule:: kripodb.makebits
    :members:
kripodb.modifiedtanimoto
========================

.. automodule:: kripodb.modifiedtanimoto
    :members:kripodb.pharmacophores
======================

.. automodule:: kripodb.pharmacophores
    :members:Incremental update
==================

.. contents::

The Kripo data set can be incrementally updated with new PDB entries.

1. Create staging directory
---------------------------

Setup path with update scripts using::

    export SCRIPTS=$PWD/../kripodb/update_scripts

Create a new directory::

  mkdir staging
  cd ..

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

The `ids.txt` file must contain a list of PDB identifiers which have not been processed before.
It can be fetched from https://www.rcsb.org/.

Adjust the PDB save location in the `singleprocess.py` script to the staging directory.

Run the following command to generate fragments/pharmacophores/fingerprints for each PDB listed in `ids.txt`::

  python singleprocess.py


3. Create fragment information
------------------------------

1. Fragment shelve
^^^^^^^^^^^^^^^^^^

Where the fragment came from is stored in a Python shelve file.
It can be generated from the pharmacophore files using::

  compiledDatabase.py

2. Fragment sdf
^^^^^^^^^^^^^^^

The data generated thus far contains the molblocks of the ligands and atom nrs of each fragment.
The fragment molblocks can be generated into a fragment sdf file with::

  fragid2sd.py fragments.shelve > fragments.sd


3. Pharmacophores
^^^^^^^^^^^^^^^^^

The raw pharmacophores are stored in the FRAGMENT_PPHORES sub-directory.
Each pocket has a \*_pphore.sd.gz file which contains the pharmacophore points of the whole pocket and
a \*_pphores.txt file which contains the indexes of pharmacophore points for each sub pocket or fragment.
The raw pharmacophores of the update can be added to the existing pharmacophores datafile with::

    cp ../current/pharmacophores.h5 .
    kripodb pharmacophores add FRAGMENT_PPHORES pharmacophores.h5

4. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

Step 4 and 5 can be submitted to scheduler with::

   jid_db=$(sbatch --parsable -n 1 -J db_append $SCRIPTS/db_append.sh)

5. Populate PDB metadata in fragments database
----------------------------------------------
The following command will updated the PDB metadata to fragments database::

    kripodb fragments pdb fragments.sqlite

6. Check no fragments are duplicated
------------------------------------

The similarity matrix can not handle duplicates. It will result in addition of scores::

    jid_dups=$(sbatch --parsable -n 1 -J check_dups --dependency=afterok:$jid_db $SCRIPTS/incremental_duplicates.sh)

7. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between the new and existing fingerprints and between new fingerprints themselves can be calculated with::

    current_chunks=$(ls ../current/*fp.gz |wc -l)
    all_chunks=$(($current_chunks + 1))
    jid_fpneigh=$(sbatch --parsable -n $all_chunks -J fpneigh --dependency=afterok:$jid_dups $SCRIPTS/incremental_similarities.sh)
    jid_merge_matrices=$(sbatch --parsable -n 1 -J merge_matrices --dependency=afterok:$jid_fpneigh $SCRIPTS/incremental_merge_similarities.sh)

8. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. note:: Converting the pairs file into a dense matrix goes quicker with more memory.

    The frame size (-f) should be as big as possible, 100000000 requires 6Gb RAM.

The following commands converts the pairs into a compressed dense matrix::

    jid_compress_matrix=$(sbatch --parsable -n 1 -J compress_matrix --dependency=afterok:$jid_merge_matrices $SCRIPTS/freeze_similarities.sh)

The output of this step is ready used to find similar fragments,
using either the webservice with the `kripodb serve` command or with the `kripodb similarities similar` command directly.

9. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The current and new pharmacophores need to be combined::

    mv staging/FRAGMENT_PPHORES staging/FRAGMENT_PPHORES.new
    rsync -a current/FRAGMENT_PPHORES staging/FRAGMENT_PPHORES
    rm -r staging/FRAGMENT_PPHORES.new

.. todo:: rsync of current/FRAGMENT_PPHORES to destination, maybe too slow due large number of files.
    Switch to move old pharmacohores and rsync new pharmacophores into it when needed.

The current and new fingerprints need to be combined::

    cp -n current/*.fp.gz staging/

The staging can be made current with the following commands::

    mv current old && mv staging current

9.1 Merge fingerprint files (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To keep the number of files to a minimum it is advised to merge the fingerprint files from incremental updates of a year.

The incremental fingerprint files are named like `out.<year><week>.fp.gz`, to generate kripo_fingerprints_<year>_fp.gz run::

    sbatch --parsable -n 1 -J merge_fp $SCRIPTS/incremental_merge_fp.sh <year>

10.0 Update web service
-----------------------

The webservice running at http://3d-e-chem.vu-compmedchem.nl/kripodb must be updated with the new datafiles.

The following files must copied to the server

* fragments.sqlite
* pharmacophores.h5
* similarities.packedfrozen.h5

The webservice must be restarted.

To show how up to date the webservice is the release date of the latest PDB is stored in `version.txt` which can be reached at http://3d-e-chem.vu-compmedchem.nl/kripodb/version.txt
The content `version.txt` must be updated.
DiVE visualization
==================

DiVE homepage at https://github.com/NLeSC/DiVE

The Kripo similarity matrix can be embedded to 2D or 3D using largevis and then visualized using DiVE.

Steps

1. LargeVis input file from Kripo similarity matrix
2. Perform embedding using LargeVis
3. Generate DiVE metadata datafiles
4. Create DiVE input file

Input datasets

1. only fragment1 or whole unfragmented ligands
2. all fragments
3. only gpcr frag1
4. only kinase frag1
5. only gpcr and kinase frag1

Output datasets

1. 2D
2. 3D

1. LargeVis input file from Kripo similarity matrix
---------------------------------------------------

Dump the similarity matrix to csv of \*frag1 fragments::

    kripodb similarities export --no_header --frag1 similarities.h5 similarities.frag1.txt

Similarities between GPCR pdb entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `GPCRDB web service <http://gpcrdb.org/services/reference/#!/structure/structure_list>`_ to fetch a list of PDB codes which contain GPCR proteins::

    curl -X GET --header 'Accept: application/json' 'http://gpcrdb.org/services/structure/' | jq  -r '.[] | .pdb_code' > pdb.gpcr.txt

Dump the similarity matrix to csv::

    kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.txt similarities.h5 similarities.frag1.gpcr.txt

Similarities between GPCR and Kinase pdb entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `KLIFS KNIME nodes <https://github.com/3D-e-Chem/knime-klifs>`_ to create a file with of PDB codes of Kinases called `pdb.kinase.txt`.

Dump the similarity matrix to csv::

    cat pdb.gpcr.txt pdb.kinase.txt > pdb.gpcr.kinase.txt
    kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.kinase.txt similarities.h5 similarities.frag1.gpcr.kinase.txt

2. Perform embedding using LargeVis
-----------------------------------

Get or compile LargeVis binaries from https://github.com/lferry007/LargeVis

Compile using miniconda::

    conda install gsl gcc
    cd LargeVis/Linux
    c++ LargeVis.cpp main.cpp -o LargeVis -lm -pthread -lgsl -lgslcblas -Ofast -Wl,-rpath,$CONDA_PREFIX/lib -march=native -ffast-math
    cp LargeVis $CONDA_PREFIX/bin/

Then embed frag1 similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.3d.txt

Then embed frag1 similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.2d.txt

Then embed similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.txt -output largevis.3d.txt

Then embed similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.txt -output largevis.2d.txt


The `kripo export` in step 1 and the LargeVis command can be submitted to scheduler with::

   sbatch -n 1 $SCRIPTS/dive_frag1.sh
   sbatch -n 1 $SCRIPTS/dive_frag1_gpcr_kinase.sh

3. Generate DiVE metadata datafiles
-----------------------------------

Command to generate properties files::

    wget -O uniprot.txt 'http://www.uniprot.org/uniprot/?query=database:pdb&format=tab&columns=id,genes(PREFERRED),families,database(PDB)'
    kripodb dive export --pdbtags pdb.gpcr.txt --pdbtags pdb.kinase.txt fragments.sqlite uniprot.txt

Will generate in current working directory the following files:

* kripo.props.txt
* kripo.propnames.txt

4. Create DiVE input file
-------------------------

DiVE has a script which can combine the LargeVis coordinates together with metadata. 
Download the MakeVizDataWithProperMetadata.py script from https://github.com/NLeSC/DiVE/blob/master/scripts_prepareData/MakeVizDataWithProperMetadata.py

For more information about the script see https://github.com/NLeSC/DiVE#from-output-of-largevis-to-input-of-dive .

Example command to generate new DiVE input file::

    python MakeVizDataWithProperMetadata.py -coord largevis2.similarities.frag1.gpcr.kinase.txt -metadata kripo.props.txt -np kripo.propnames.txt -json largevis2.similarities.frag1.gpcr.kinase.json -dir .

The generated file (largevis2.similarities.frag1.gpcr.kinase.json) can be uploaded at https://nlesc.github.io/DiVE/ to visualize.Baseline update
===============

.. contents::

The Kripo data set is generated from scratch every year or when algorithms change.

1. Create staging directory
---------------------------

Setup path with update scripts using::

    export SCRIPTS=$PWD/../kripodb/update_scripts

Create a new directory::

  mkdir staging
  cd ..

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

Use directory listing of new pdb files as input::

  ls $PDBS_ADDED_DIR | pdblist2fps_final_local.py

.. todo:: Too slow when run on single cpu.
    Chunkify input, run in parallel and merge results

.. _create-fragment-information:

3. Create fragment information
------------------------------

1. Fragment shelve
^^^^^^^^^^^^^^^^^^

Where the fragment came from is stored in a Python shelve file.
It can be generated from the pharmacophore files using::

  compiledDatabase.py

2. Fragment sdf
^^^^^^^^^^^^^^^

The data generated thus far contains the molblocks of the ligands and atom nrs of each fragment.
The fragment molblocks can be generated into a fragment sdf file with::

  fragid2sd.py > fragments.sd

3. Pharmacophores
^^^^^^^^^^^^^^^^^

The raw pharmacophores are stored in the FRAGMENT_PPHORES sub-directory.
Each pocket has a \*_pphore.sd.gz file which contains the pharmacophore points of the whole pocket and
a \*_pphores.txt file which contains the indexes of pharmacophore points for each sub pocket or fragment.
The raw pharmacophores need to be added to the pharmacophores datafile with::

    kripodb pharmacophores add FRAGMENT_PPHORES pharmacophores.h5

4. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

Step 4 and 5 can be submitted to scheduler with::

   jid_db=$(sbatch --parsable -n 1 -J db_append $SCRIPTS/db_append.sh)


5. Populate PDB metadata in fragments database
----------------------------------------------
The following command will updated the PDB metadata to fragments database::

    kripodb fragments pdb fragments.sqlite


6. Check no fragments are duplicated
------------------------------------

The similarity matrix can not handle duplicates. It will result in addition of scores::

    jid_dups=$(sbatch --parsable -n 1 -J check_dups --dependency=afterok:$jid_db $SCRIPTS/baseline_duplicates.sh)

7. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between fingerprints can be calculated with::

    all_chunks=$(ls *fp.gz |wc -l)
    jid_fpunzip=$(sbatch --parsable -n $all_chunks -J fpunzip --dependency=afterok:$jid_dups $SCRIPTS/baseline_fpunzip.sh)
    nr_chunks="$(($all_chunks * $all_chunks / 2 - $all_chunks))"
    jid_fpneigh=$(sbatch --parsable -n $nr_chunks -J fpneigh --dependency=afterok:$jid_fpunzip $SCRIPTS/baseline_similarities.sh)
    jid_fpzip=$(sbatch --parsable -n $all_chunks -J fpzip --dependency=afterok:$jid_fpneigh $SCRIPTS/baseline_fpzip.sh)
    jid_merge_matrices=$(sbatch --parsable -n 1 -J merge_matrices --dependency=afterok:$jid_fpneigh $SCRIPTS/baseline_merge_similarities.sh)

To prevent duplicates similarities of a chunk against itself should ignore the upper triangle.

.. todo:: Don't fpneigh run sequentially but submit to batch queue system and run in parallel

8. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. tip:: Converting the pairs file into a dense matrix goes quicker with more memory.

The following commands converts the pairs into a compressed dense matrix::

    jid_compress_matrix=$(sbatch --parsable -n 1 -J compress_matrix --dependency=afterok:$jid_merge_matrices $SCRIPTS/freeze_similarities.sh)

The output of this step is ready to be served as a webservice using the `kripodb serve` command.

9. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The staging can be made current with the following commands::

    mv current old
    mv staging current

10.0 Update web service
-----------------------

The webservice running at http://3d-e-chem.vu-compmedchem.nl/kripodb must be updated with the new datafiles.

The following files must copied to the server

* fragments.sqlite
* pharmacophores.h5
* similarities.packedfrozen.h5

The webservice must be restarted.

To show how up to date the webservice is the release date of the latest PDB is stored in `version.txt` which can be reached at http://3d-e-chem.vu-compmedchem.nl/kripodb/version.txt
The content `version.txt` must be updated.
.. KripoDB documentation master file, created by
   sphinx-quickstart on Thu Feb  4 21:13:05 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to KripoDB's documentation!
===================================

For installation and usage see https://github.com/3D-e-Chem/kripodb/blob/master/README.md

.. toctree::
    :maxdepth: 2

    Data update <data-update.rst>
    DiVe visualization <dive-visualization.rst>
    API <api.rst>
    CLI <cli.rst>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
kripodb.frozen
==============

.. automodule:: kripodb.frozen
    :members:
kripodb.pairs
=============

.. automodule:: kripodb.pairs
    :members:API
===

.. toctree::
    :maxdepth: 2

    kripodb.canned <canned.rst>
    kripodb.db <db.rst>
    kripodb.dive <dive.rst>
    kripodb.frozen <frozen.rst>
    kripodb.hdf5 <hdf5.rst>
    kripodb.makebits <makebits.rst>
    kripodb.modifiedtanimoto <modifiedtanimoto.rst>
    kripodb.pairs <pairs.rst>
    kripodb.pharmacophores <pharmacophores.rst>
    kripodb.pdb <pdb.rst>
    kripodb.script <script.rst>
    kripodb.webservice <webservice.rst>

