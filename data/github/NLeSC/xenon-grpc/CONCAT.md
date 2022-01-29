# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

## [3.0.2] - 2020-03-23

### Changed

* Depends on Xenon 3.1.0

## [3.0.1] - 2019-09-11

### Changed

* Depends on Xenon 3.0.4
* Depends on Xenon cloud adaptors 3.0.2

## [3.0.0] - 2019-06-14

### Added

* getDefaultRuntime rpc to SchedulerService ([#44](https://github.com/xenon-middleware/xenon-grpc/issues/44))
* start_time and temp_space fields to JobDescription message ([#43](https://github.com/xenon-middleware/xenon-grpc/issues/43))
* [at](https://linux.die.net/man/1/at) scheduler

### Changed

* Replaced tasks+cores+nodes fields in JobDescription message with nodes+processes+thread fields ([#625](https://github.com/xenon-middleware/xenon/issues/625)).
* Require Java 11 or greater, as xenon package has same compatibility ([#42](https://github.com/xenon-middleware/xenon-grpc/issues/42https://github.com/xenon-middleware/xenon-grpc/issues/42))
* Upgraded to Xenon 3.0.0 ([#40](https://github.com/xenon-middleware/xenon-grpc/issues/40))

### Removed

* hdfs filesystem
* options field in JobDescription message ([#630](https://github.com/xenon-middleware/xenon/issues/630))

## [2018-03-14] 2.3.0

### Added

* scheduler argument for job description (#38)

### Changed

* Depends on Xenon 2.6.0

## [2018-03-06] 2.2.1

### Fixed

* hadoop/grpc netty version conflict (#37)

## [2018-03-05] 2.2.0

### Added

* support for KeytabCredential (#33)
* supportedCredentials (#35)

### Changed

* Depends on Xenon 2.5.0

### Fixed

* FileSystemAdaptorDescription fields synced (#36)

## [2018-02-26] 2.1.0

### Added

* Name to JobDescription and JobStatus
* Max memory to JobDescription

### Changed

* Use latest dependencies and plugins
* Depends on Xenon 2.4.0

## [2018-01-04] 2.0.1

### Changed

* Use latest dependencies and plugins

### Fixed

* Class of exception lost in translation [#32]
* Interactive job: sometimes output to stdout is repeated, sometimes skipped [#34]

## [2017-11-07] 2.0.0

Initial release
gRPC (http://www.grpc.io/) server for Xenon (https://xenon-middleware.github.io/xenon/).

Can be used to use Xenon in a non-java based language.
For example pyxenon (https://github.com/NLeSC/pyxenon) uses the Xenon gRPC server.

The server tries to mimic the Xenon library API as much as possible, differences are described in the [proto file](src/main/proto/xenon.proto) .

[![Build Status](https://travis-ci.org/xenon-middleware/xenon-grpc.svg?branch=master)](https://travis-ci.org/xenon-middleware/xenon-grpc)
[![Build status](https://ci.appveyor.com/api/projects/status/tep8bad05e76a69w/branch/master?svg=true)](https://ci.appveyor.com/project/NLeSC/xenon-grpc/branch/master)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=xenon-middleware_xenon-grpc&metric=alert_status)](https://sonarcloud.io/dashboard?id=xenon-middleware_xenon-grpc)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=xenon-middleware_xenon-grpc&metric=coverage)](https://sonarcloud.io/dashboard?id=xenon-middleware_xenon-grpc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1043481.svg)](https://doi.org/10.5281/zenodo.1043481)

# Install

On [releases page](https://github.com/xenon-middleware/xenon-grpc/releases) download a tarball (or zipfile).

The tarball can be installed with:

```bash
tar -xf xenon-grpc-shadow*.tar
```

Add `xenon-grpc*/bin` to your PATH environment variable for easy usage.

# Usage

To start the grpc server with default arguments run

```bash
./xenon-grpc*/bin/xenon-grpc
```

To get help run

```bash
./xenon-grpc*/bin/xenon-grpc --help
```

Or call the jar directly with

```bash
java -jar xenon-grpc-*/lib/xenon-grpc-*-all.jar
```

# Development

## Run server

```bash
./gradlew installDist
./build/install/xenon-grpc/bin/xenon-grpc
```

## Run client

For use [polyglot](https://github.com/grpc-ecosystem/polyglot)

```bash
wget https://github.com/grpc-ecosystem/polyglot/releases/download/v2.0.0/polyglot.jar
java -jar polyglot.jar --proto_discovery_root=src/main/proto list_services
echo {} | java -jar polyglot.jar call --endpoint=localhost:50051 --full_method=xenon.SchedulerService/getAdaptorDescriptions
```

## Python client

Compile proto into python stubs

```sh
pip install grpcio grpcio-tools
xenon-grpc --proto > xenon.proto
python -m grpc_tools.protoc -I. --python_out=. --grpc_python_out=. xenon.proto
```

Now use the generated stubs, see https://grpc.io/docs/tutorials/basic/python.html#creating-the-client

## Mutual TLS

Create self-signed certificate and use for server and client on same machine.
Make sure `Common Name` field is filled with hostname of machine.
See http://httpd.apache.org/docs/2.4/ssl/ssl_faq.html#selfcert

```bash
openssl req -new -x509 -nodes -out server.crt -keyout server.key
./build/install/xenon-grpc/bin/xenon-grpc --server-cert-chain server.crt --server-private-key server.key --client-cert-chain server.crt
```

Test with polyglot

```bash
echo {} | java -jar polyglot.jar call --endpoint=<hostname as used in certificate>:50051 --full_method=xenon.FileSystemService/getAdaptorNames --use_tls=true --tls_client_cert_path=$PWD/server.crt --tls_client_key_path=$PWD/server.key --tls_ca_cert_path=$PWD/server.crt
```

In a ipython shell with generated stubs in working directory:

```python
import grpc
import xenon_pb2
import xenon_pb2_grpc
import socket

creds = grpc.ssl_channel_credentials(
    root_certificates=open('../../../server.crt').read(),
    private_key=open('../../../server.key', 'rb').read(),
    certificate_chain=open('../../../server.crt', 'rb').read()
)
channel = grpc.secure_channel(socket.gethostname() + ':50051', creds)
stub = xenon_pb2_grpc.XenonJobsStub(channel)
response = stub.getAdaptorDescriptions(xenon_pb2.Empty())
print(response)
```

## New release

```sh
./gradlew build
```

Generates application tar/zip in `build/distributions/` directory.

1. Bump version in `build.gradle`, `CITATION.cff`, add version to `CHANGELOG.md` and commit/push
2. Create a new GitHub release
3. Upload the files in `build/distributions/` directory to that release
4. Publish release
5. Edit Zenodo entry, add [Xenon doi](https://doi.org/10.5281/zenodo.597993) as `is referenced by this upload`.
