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

