# Changelog

### 2.0.0

* Add Zenodo integration for uploading the created CSV file


### 1.1.0

* Add executable (frbcatdb-image) to add images to the database 
* Add functionality to decode_VOEvent to dump database to a CSV file

# frbcatdb
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build Status](https://travis-ci.org/TRASAL/frbcatdb.svg?branch=master)](https://travis-ci.org/TRASAL/frbcatdb)[![codecov](https://codecov.io/gh/TRASAL/frbcatdb/branch/master/graph/badge.svg)](https://codecov.io/gh/TRASAL/frbcatdb)[![Codacy Badge](https://api.codacy.com/project/badge/Grade/de13488f778e4843a8922ee2417a3416)](https://www.codacy.com/app/omrubi/frbcatdb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=TRASAL/frbcatdb&amp;utm_campaign=Badge_Grade)[![Readthedocs badge](https://readthedocs.org/projects/frbcatdb/badge/)](http://frbcatdb.readthedocs.io/en/latest/?badge=latest)[![PyPI version](https://badge.fury.io/py/pyfrbcatdb.svg)](https://badge.fury.io/py/pyfrbcatdb)[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1326399.svg)](https://doi.org/10.5281/zenodo.1326399)

The frbcatdb is a database to store a catalog of Fast Radio Bursts (FRBs).
The DB is intended to contain old FRB events as well as new FRBs detected by the
AA-ALERT FRB detection pipeline from Apertif observations and also possible follow-up observations or others FRBs detected by other telescopes.
The frbcatdb is attached to the VOEvent backbone and uses this infrastructure as its source.

The `db` folder contains scripts to create an empty frbcat DB (`create_db.csh`),
to import it from an existing dump file (`import_db.sh`) and
to dump an existing DB to a dump file (`dump_db.csh`).
It also contains the model (Entity-Relationship diagram) to be opened with mysql-workbench. ![frbcatdb ER diagram](db/relationships.real.compact.png)

The `pyfrbcatdb` is Python package for manipulating the frbcatdb and its linking
with the VOEvent backbone.

# pyfrbcatdb usage
A default configuration file is installed in /etc/pyfrbcatdb/dbase.config. In this file the FRBCat database configuration can be defined. Alternatively, a user may supply their own configuration file with a command line argument of the executable, or define the database configuration via argument switches, or, alternatively via environment variables.

For inserting a VOEvent XML file into the FRBCat database, the decode_VOEvent executable is used:
```
usage: decode_VOEvent [-h] [-c MY_CONFIG] --dbName DBNAME [--dbHost DBHOST]
                      [--dbPort DBPORT] --dbUser DBUSER
                      [--dbPassword DBPASSWORD] [--CSV CSV] [--log LOG]
                      [--zenodo ZENODO]
                      [VOEvent [VOEvent ...]]

Process VOEvent XML file and add it to FRB database Args that start with '--'
(eg. --dbName) can also be set in a config file
(/etc/pyfrbcatdb/dbase.config or specified via -c). Config file syntax
allows: key=value, flag=true, stuff=[a,b,c] (for details, see syntax at
https://goo.gl/R74nmi). If an arg is specified in more than one place, then
commandline values override environment variables which override config file
values which override defaults.

positional arguments:
  VOEvent               List of VOEvent XML files

optional arguments:
  -h, --help            show this help message and exit
  -c MY_CONFIG, --my-config MY_CONFIG
                        config file path
  --dbName DBNAME       name postgres database [env var: dbNameFRBCat]
  --dbHost DBHOST       name postgres database [env var: dbHostFRBCat]
  --dbPort DBPORT       name postgres database [env var: dbPortFRBCat]
  --dbUser DBUSER       user postgres database [env var: dbUserFRBCat]
  --dbPassword DBPASSWORD
                        user postgres database password [env var:
                        dbPasswordFRBCat]
  --CSV CSV             CSV filename to dump database to [env var: CSVFRBCat]
  --log LOG             log file, default=[HOME]/pyfrbcatdb_decode.log
  --zenodo ZENODO       upload CSV to Zenodo, access token [env var: zenodoFRBCat]
```
For inserting an image into the database, the frbcatdb-image executable is used. Apart from the database configuration, the tool takes two positional arguments. The first is the filename of the image to be added, the second is the 'id' in the 'radio measurement params' table that the image should be connected to:
```
usage: frbcatdb-image [-h] [-c MY_CONFIG] --dbName DBNAME [--dbHost DBHOST]
                      [--dbPort DBPORT] --dbUser DBUSER
                      [--dbPassword DBPASSWORD] [--caption CAPTION]
                      [--title TITLE]
                      filename rmpid

Create VOEvent XML file from FRB database Args that start with '--' (eg.
--dbName) can also be set in a config file (/data/github/venv-
aa/etc/pyfrbcatdb/dbase.config or specified via -c). Config file syntax
allows: key=value, flag=true, stuff=[a,b,c] (for details, see syntax at
https://goo.gl/R74nmi). If an arg is specified in more than one place, then
commandline values override environment variables which override config file
values which override defaults.

positional arguments:
  filename              Name of file to fetch from
  rmpid                 rmp_id

optional arguments:
  -h, --help            show this help message and exit
  -c MY_CONFIG, --my-config MY_CONFIG
                        config file path
  --dbName DBNAME       name postgres database [env var: dbNameFRBCat]
  --dbHost DBHOST       name postgres database [env var: dbHostFRBCat]
  --dbPort DBPORT       name postgres database [env var: dbPortFRBCat]
  --dbUser DBUSER       user postgres database [env var: dbUserFRBCat]
  --dbPassword DBPASSWORD
                        user postgres database password [env var:
                        dbPasswordFRBCat]
  --caption CAPTION     figure caption
  --title TITLE         figure title
```

For extracting a VOEvent from the FRBCat database, the create_VOEvent executable is used. Note that some features might still be missing for the current release from this utility.
```
usage: create_VOEvent [-h] [-c MY_CONFIG] --dbName DBNAME [--dbHost DBHOST]
                      [--dbPort DBPORT] --dbUser DBUSER
                      [--dbPassword DBPASSWORD] [--log LOG]
                      frb_ids [frb_ids ...]

Create VOEvent XML file from FRB database Args that start with '--' (eg.
--dbName) can also be set in a config file
(/etc/pyfrbcatdb/dbase.config or specified via -c). Config
file syntax allows: key=value, flag=true, stuff=[a,b,c] (for details, see
syntax at https://goo.gl/R74nmi). If an arg is specified in more than one
place, then commandline values override environment variables which override
config file values which override defaults.

positional arguments:
  frb_ids               List of frbs ids

optional arguments:
  -h, --help            show this help message and exit
  -c MY_CONFIG, --my-config MY_CONFIG
                        config file path
  --dbName DBNAME       name postgres database [env var: dbNameFRBCat]
  --dbHost DBHOST       name postgres database [env var: dbHostFRBCat]
  --dbPort DBPORT       name postgres database [env var: dbPortFRBCat]
  --dbUser DBUSER       user postgres database [env var: dbUserFRBCat]
  --dbPassword DBPASSWORD
                        user postgres database password [env var:
                        dbPasswordFRBCat]
  --log LOG             log file, default=[HOME]/pyfrbcatdb_create.log
```
