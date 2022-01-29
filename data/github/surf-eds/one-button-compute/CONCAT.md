# One button compute

[![Build Status](https://travis-ci.org/surf-eds/one-button-compute.svg?branch=master)](https://travis-ci.org/surf-eds/one-button-compute)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1033817.svg)](https://doi.org/10.5281/zenodo.1033817)

One button compute is aweb site that runs a workflow.

# Feature/Limitations

* Workflow is a single file in [Common Workflow format](http://www.commonwl.org/)
* Workflow must take single input file (--input option) and generates a single output file (--output option)
* Web application runs workflow on directory of input files
* The workflow, the directory with input files is downloaded from a remote storage server
* The directory with output files is uploaded to a remote storage server

The remote storage server can be WebDAV or S3 or Swift.

* The WebDAV server used for production is BeeHub (https://www.beehub.nl)
* The WebDAV server used for local development can be the Docker container `nlesc/xenon-webdav` (https://hub.docker.com/r/nlesc/xenon-webdav/)
* The S3 server used for local development is Minio (https://minio.io) instance, Minio is a single user S3 compliant server
* The Swift server used for production is a Openstack Swift (http://swift.openstack.org) instance

# Requirements

* Python2
* Docker
* Read/write access to a remote storage server. Can be a WebDAV or S3 or Swift server/account.

# Install

## 1. Redis server

Redis server is used to perform computation asynchronous from http request.

Use Docker to start a redis server

```
docker run -d -p 6379:6379 redis
```

Note!: When Celery workers are going to be run on different machines make sure they can connect to the redis server.

## 2. Install dependencies

Install the Python dependencies with
```
pip install -r requirements.txt
```

## 3. Configure the application

```
cp settings.cfg-dist settings.cfg
```

Configure remote storage type, location and credentials in settings.cfg.

## S3 development server (optional)

A Minio server can be started with
```
mkdir -p minio/export
docker run -d --name obc-minio -p 9000:9000 -v $PWD/minio:/root minio/minio /root/export
docker logs obc-minio
```

The log output contains the credentials, urls and access instructions.

To use known credentials from settings.cfg start it with
```
docker run -d --name obc-minio -p 9000:9000 -v $PWD/minio:/root -e "MINIO_ACCESS_KEY=AKIAIOSFODNN7EXAMPLE" \
  -e "MINIO_SECRET_KEY=wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY" minio/minio /root/export
```

Use `mc` (https://docs.minio.io/docs/minio-client-quickstart-guide) as CLI client.

## WebDAV development server (optional)

A WebDAV server can be started with
```
docker run -d nlesc/xenon-webdav
```

Read/write can be done in `~/xenon/uploads` path with xenon:javagat credentials.

Use `cadaver` (http://www.webdav.org/cadaver/) as CLI client.

## Reverse proxy (optional)

Configure Nginx as reverse proxy for the flask app port 5000.

```
location / {
  proxy_pass http://localhost:5000;
}
```

## Auto start (optional)

Automatically start one-button-compute on boot with upstart file

```
cat /etc/init/onebuttoncompute.conf
# Running on port 5000

description "One button compute"

start on filesystem or runlevel [2345]
stop on runlevel [!2345]

script
  cd /opt/one-button-compute
  python onebuttoncompute.py
end script
```

# Run

Start Celery worker and web server with
```
celery worker -A onebuttoncompute.celery &
python onebuttoncompute.py
```

# Usage

Add a CWL workflow and input files to remote storage.
See `example/` sub-directory for an example workflow.

Go to http://localhost:5000/ (or http://&lt;server-name&gt;/ when reverse proxy is setup) to submit a computation.
Examples which can be used in One Button Compute application.

For more examples see https://github.com/surf-eds/cwl-examples

# Prerequisites

Install CWL runner.
```
pip install cwl-runner
```

# Word count example

Simple docker container that performs a word count.

## Run using Docker

It takes 2 arguments:
1. Input text file
2. Output text file

Create a input file and run with input and output volumes mounted.
```
echo Lorem ipsum dolor sit amet > input
docker run -ti --rm -u $UID -v $PWD:/input -v $PWD:/output wca /input/input /output/output
cat output
```

## Run using cwl-runner

```
./example/cwa.tool.cwl --input README.md --output README.wc
```
This will execute cwa script inside Docker container using cwl-runner.

## Run from web app

Must use version >= v2.0.0 and < 3.0.0 of one-button-compute repo.

1. Upload a text file and workflow file (cwa.tool.cwl) to the remote storage server configured in the settings.cfg file.
2. Create a output directory on remote storage server.
3. Submit in web application

* CWL workflow file: cwa.tool.cwl
* Input file:
* Output directory:

## Run using cwl-runner with multiple files


The job order file (cwa-files.job.yml) contains the list of input files and output filenames.

```
./example/cwa-files.tool.cwl cwa-files.job.yml
```

### Run from web app with multiple files and using Minio server

Must use version >= v3.0.0 of one-button-compute repo.

Requires a Minio server, see "S3 development server" section ../README.md for instructions.

With S3_ROOT of 'http://localhost:9000/mybucket/obc'
```
mc config host add myminio http://localhost:9000 *** ***
mc mb myminio/mybucket
mc cp cwa.tool.cwl myminio/mybucket/obc/run1/cwa.tool.cwl
mc cp README.md myminio/mybucket/obc/run1/input/file1.txt
mc cp cwa.tool.cwl myminio/mybucket/obc/run1/input/file2.txt
```

In the One Button Compute web interface fill form with

* Input directory = run1/input
* CWL workflow file = run1/cwa.tool.cwl
* Output directory = run1/output

The output can be cleared with
```
mc rm -r --force myminio/mybucket/obc/run1/output
```

### Use Swift storage

To run a cwl from Openstack Swift storage.

```
# Create a container for input
swift post -w $OS_USERNAME -r $OS_USERNAME eds
# Upload input file
swift upload --object-name pointcloud2images/input/rock-section.zip eds rock-section.zip
# Upload cwl
swift upload --object-name pointcloud2images/images2pointcloud.workflow.packed.cwl eds images2pointcloud.workflow.packed.cwl
```

In the One Button Compute web interface fill form with

* Input directory = pointcloud2images/input
* CWL workflow file = pointcloud2images/images2pointcloud.workflow.packed.cwl
* Output directory = pointcloud2images/output

The output can be cleared with
```
swift delete -p pointcloud2images/output eds
```

Use s3curl with s3sara alias in ~/.s3curl
```
# List bucket
s3curl.pl --id s3sara -- https://s3.swift.surfsara.nl/eds |xml_pp
# Download
s3curl.pl --id s3sara -- https://s3.swift.surfsara.nl/eds/pointcloud2images/input/rock-section.zip > rock-section.zip
# Upload
s3curl.pl --id s3sara --put=rr.zip -- https://s3.swift.surfsara.nl/eds/pointcloud2images/input/rr.zip
s3curl.pl --id s3sara --delete -- https://s3.swift.surfsara.nl/eds/pointcloud2images/input/rr.zip
```
