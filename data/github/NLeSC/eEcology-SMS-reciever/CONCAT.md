eEcology-SMS-reciever
=====================

[![Build Status](https://travis-ci.org/NLeSC/eEcology-SMS-reciever.svg?branch=master)](https://travis-ci.org/NLeSC/eEcology-SMS-reciever)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/NLeSC/eEcology-SMS-reciever/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/NLeSC/eEcology-SMS-reciever/?branch=master)
[![Code Coverage](https://scrutinizer-ci.com/g/NLeSC/eEcology-SMS-reciever/badges/coverage.png?b=master)](https://scrutinizer-ci.com/g/NLeSC/eEcology-SMS-reciever/?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.47325.svg)](https://doi.org/10.5281/zenodo.47325)

Webservice to store SMS messages into a database table.

Accepts SMS messages from SMSSync Android app see https://play.google.com/store/apps/details?id=org.addhen.smssync, http://smssync.ushahidi.com/ and https://github.com/ushahidi/SMSSync .

Getting Started
---------------

    cd <directory containing this file>
    pipenv --two
    pipenv shell
    pip install -r requirements.txt
    cp development.ini-dist development.ini

Edit development.ini to configure db connection, etc.

Create sms database schema:

    psql -h db.e-ecology.sara.nl eecology < sms.sql

Grant `<someone>` user rights to perform inserts on table, see comments in `sms.sql` for required grants.

Start service:

    pserve development.ini

Service running at http://localhost:6566/sms/

SMSSync configuration
---------------------

Configure Sync URL with:

* Secret key = Same as value of 'secret_key' key in *ini file.
* Keyword = '^ID', all tracker messages start with 'ID'
* URL = Depends on where you run this server and if it is reversed proxied.
* HTTP Method = POST
* Data Format = URLEncoded

Docker build
------------

### Construct image

1. `sudo docker build -t sverhoeven/smsreciever:1.0.0 .`
2. Export or push to registry

### Run container

1. Import or pull from registry
2. `sudo docker run -p 6566:6566 --env DB_URL="postgresql://*******:********@db.e-ecology.sara.nl/eecology?sslmode=require" --env SECRET_KEY=supersecretkey -d --name smsreciever sverhoeven/smsreciever:1.0.0`

Error log is available with `sudo docker logs smsreciever`.

Web application will run on http://localhost:6566/sms/

Database upgrades
-----------------

The latest schema is specified in `sms.sql`.

Using alembic (http://pythonhosted.org/alembic/) for database migrations.

To create a new migration step run:

    alembic revision -m "Timestamp without time zone"
    # edit alembic/versions/*py files to specify changes

To upgrade an existing schema run:

    alembic upgrade head --sql | psql ...


Tests
-----

### Unit tests

The unit tests can be run with
```
nosetests
```

### Functional tests

The functional tests needs

* Postgresql database with PostGIS extension
* Postgresql user with schema and user creation permission

The tests can be run with

```
DB_URL=postgresql://postgres:mysecretpassword@172.17.0.2/postgres nosetests -a functional
```

Postgis in a docker container can be used to test.
```
docker run --name some-postgis -e POSTGRES_PASSWORD=mysecretpassword -d mdillon/postgis
```

1.0.10
------

### Fixed

- Catch InvalidRequestError and perform rollback

### Added

- /status endpoint will give error when no positions have been received recently

1.0.9
-----

### Added

- Example apache config
- Example nagios check
- Functional tests

1.0.8
-----

### Fixed

- Message with no positions, caused a permission denied which caused transaction is aborted error

1.0.7
-----

- Force UTC timezone during insert

1.0.6
-----

- When raw message already exists return OK, instead of error

1.0.5
-----

- Inserts raw message before parsing message, instead of parsing message and then inserting raw message and parsed message.

1.0.4
-----

- Use timestamp without time zone, as it is used in other eEcology tables.

1.0.3
-----

- Debug string changed, detect it by number of columns
- Treat each insert as a seperate transaction. Position duplicates are skipped, but new message is stored.

1.055
-----

- Date time has changed in firmware 1.055 to

    * date = 2 digits year (base 2000) + 3 digits day of year
    * time = seconds of day

1.0.0
-----

-  Initial version
