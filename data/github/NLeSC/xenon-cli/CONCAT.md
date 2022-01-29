# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

## [3.0.5] - 2020-03-23

### Changed

* Conda package now noarch instead of just linux-64 ([#73](https://github.com/xenon-middleware/xenon-cli/issues/73))
* Upgraded to Xenon library 3.1.0

### Fixed

* Slurm 19 compliant ([#72](https://github.com/xenon-middleware/xenon-cli/issues/72))

## [3.0.4] - 2019-09-11

### Changed

* Upgraded to Xenon library 3.0.4

## [3.0.3] - 2019-09-09

### Changed

* Upgraded to Xenon library 3.0.3

## [3.0.2] - 2019-08-07

### Changed

* Upgraded to Xenon cloud adaptors 3.0.2

## [3.0.1] - 2019-08-01

### Changed

* Upgraded to Xenon 3.0.1
* Test against xenonmiddleware Docker images ([#42](https://github.com/xenon-middleware/xenon-docker-images/issues/42))

## [3.0.0] - 2019-06-14

### Added

* --temp-space argument ([#61](https://github.com/xenon-middleware/xenon-cli/issues/61))
* [at](https://linux.die.net/man/1/at) scheduler

### Changed

* submit and exec sub command use a tasks+cores+nodes arguments, instead of nodes+processes+thread ([#625](https://github.com/xenon-middleware/xenon/issues/625)).
* Require Java 11 or greater, as xenon package has same compatibility
* Upgraded to Xenon 3.0.0
* Switched to [testcontainers](https://www.testcontainers.org/) for testing against Docker containers

### Removed

* hdfs filesystem

## [2.4.1] - 2019-02-26

### Fixed

* s3 adaptor gives wrong error due to gson conflict

## [2.4.0] - 2018-03-14

### Added

* scheduler arguments for exec and submit

### Changed

* Upgraded to Xenon 2.6.0

## [2.3.0] - 2018-03-05

### Added

* KeytabCredential support (#46)

### Changed

* Upgraded to Xenon 2.5.0

### Fixed

* Only show credentials flags supported by adaptor (#32)

## [2.2.1] - 2018-02-28

### Changed

* Upgraded to Xenon 2.4.1

### Fixed

* slf4j multiple bindings warning
* Slurm maxtime for interactive job does not appear functional (#29)
* On InvalidLocationException return supported locations (#31)
* On UnknownPropertyException return supported props (#34)

## [2.2.0] - 2018-02-26

### Added

* --name to to submit sub commands (#59)
* --max-memory to exec and submit sub command (#58)

### Changed

* Upgraded to Xenon 2.4.0

## [2.1.0] - 2018-02-16

### Added

* --start-single-process to exec and submit sub commands (#56)
* --inherit-env to exec and submit sub commands (#55)

### Changed

* Upgraded to Xenon 2.3.2-nohadoop

## [2.0.1] - 2018-01-30

### Changed

* Upgraded to Xenon 2.3.0-nohadoop

## [2.0.0] - 2017-11-07

### Added

* Subcommands
  * mkdir
  * rename
  * wait
* Status details to jobs list
* --long format for files list (#16)
* --verbose and --stacktrace arguments
* In `xenon --help`, added type column to list of adaptors

### Changed

* Upgraded to Xenon 2.2.0
* Renamed `--format cwljson` argument to `--json`
* Xenon CLI now has same major version as Xenon

## [1.0.3] - 2017-07-20

### Fixed

* Filter scheme properties based on XenonPropertyDescription.Component.SCHEDULER or XenonPropertyDescription.Component.FILESYSTEM (#12)

## [1.0.2] - 2017-05-08

### Changed

* Upgraded to Xenon 1.2.2

## Fixed

* Weird behavior with sftp (#10)

## [1.0.1] - 2017-03-02

### Changed

* Upgraded to Xenon 1.2.1

## [1.0.0] - 2017-02-23

Initial release
# Xenon Command Line Interface

[![Build Status](https://travis-ci.org/xenon-middleware/xenon-cli.svg?branch=master)](https://travis-ci.org/xenon-middleware/xenon-cli)
[![Build status](https://ci.appveyor.com/api/projects/status/vki0xma8y7glpt09/branch/master?svg=true)](https://ci.appveyor.com/project/NLeSC/xenon-cli/branch/master)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=xenon-middleware_xenon-cli&metric=alert_status)](https://sonarcloud.io/dashboard?id=xenon-middleware_xenon-cli)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=xenon-middleware_xenon-cli&metric=coverage)](https://sonarcloud.io/dashboard?id=xenon-middleware_xenon-cli)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597603.svg)](https://doi.org/10.5281/zenodo.597603)
[![Anaconda-Server Badge](https://anaconda.org/nlesc/xenon-cli/badges/installer/conda.svg)](https://anaconda.org/NLeSC/xenon-cli)

Command line interface which uses the [Xenon library](https://xenon-middleware.github.io/xenon) to perform job and file operations.

## Install

Dependencies:

* Java runtime version 11 or greater

Goto [releases](https://github.com/xenon-middleware/xenon-cli/releases) and download a tarball (or zipfile).
The tarball can be installed with:

```bash
tar -xf build/distributions/xenon*.tar
xenon*/bin/xenon --help
```

Add `xenon*/bin` to your PATH environment variable for easy usage.

Or install with [ananconda](https://anaconda.org/nlesc/xenon-cli):

```bash
conda install -c conda-forge -c nlesc xenon-cli
```

## Usage

```bash
# List files on local filesystem
xenon filesystem file list /etc
# List files on remote filesystem using sftp
xenon filesystem sftp --location localhost list /etc
# Copy local file to remote filesystem
xenon filesystem sftp --location localhost upload /etc/passwd /tmp/copy-of-passwd
# Execute a program remotely using ssh
xenon scheduler ssh --location localhost exec /bin/hostname
# Pipe to a remote file
echo "sleep 30;echo Hello" | xenon sftp --location localhost upload - /tmp/myjob.sh
# Submit to a remote Slurm batch scheduler
xenon scheduler slurm --location ssh://localhost submit /bin/sh /tmp/myjob.sh
```

The above commands use your current username and keys from ~/.ssh.

To keep password or passphrase invisible in process list put the password in a text file (eg. 'password.txt') and then use '@password.txt' as argument.
For example:

```sh
xenon filesystem sftp --location localhost --username $USER --password @password.txt list $PWD/src
```

## Build

```sh
./gradlew build
```

Generates application tar/zip in `build/distributions/` directory.

## Tests

Requirements for the integration tests:

* [docker](https://docs.docker.com/engine/installation/), v1.13 or greater
* [docker-compose](https://docs.docker.com/compose/), v1.10 or greater

The unit and integration tests can be run with:

```sh
./gradlew check
```

## Release

1. Bump version in `build.gradle`, `conda/xenon-cli/meta.yaml` files, add version to `CHANGELOG.md` and commit/push
2. Run `./gradlew build` to build distributions
3. Create a new GitHub release
4. Upload the files in `build/distributions/` directory to that release
5. Publish release
6. Edit [Zenodo entry](https://doi.org/10.5281/zenodo.597603), add [Xenon doi](https://doi.org/10.5281/zenodo.597993) as `is referenced by this upload`.
7. Create conda release, see [conda/README.md](conda/README.md)

## Docker

Run Xenon CLI as a Docker container.

The Docker image can be build with

```sh
./gradlew docker
```

Generates a `xenonmiddleware/xenon-cli` Docker image.

To use local files use volume mounting (watch out as the path should be relative to mount point):

```sh
docker run -ti --rm xenonmiddleware/xenon-cli --user $USER -v $PWD:/work --adaptor ssh upload --source /work/somefile.txt --location localhost --path /tmp/copy-of-somefile.txt
```

## Common Workflow Language

Run Xenon CLI using a cwl-runner or as a tool in a [Common Workflow Language](http://www.commonwl.org/) workflow.

Requires `xenonmiddleware/xenon-cli` Docker image to be available locally.

Example to list contents of `/etc` directory via a ssh to localhost connection with cwl-runner:

```sh
./xenon-ls.cwl --adaptor sftp --location $USER@172.17.0.1 --certfile ~/.ssh/id_rsa --path /etc
# Copy file from localhost to working directory inside Docker container
./xenon-upload.cwl --adaptor sftp --certfile ~/.ssh/id_rsa --location $USER@172.17.0.1 --source $PWD/README.md --target /tmp/copy-of-README.md
# Copy file inside Docker container to localhost
./xenon-download.cwl --adaptor sftp --certfile ~/.ssh/id_rsa --location $USER@172.17.0.1 --source /etc/passwd --target $PWD/copy-of-passwd
```

(Replace `<user>@<host>` with actual username and hostname + expects docker with default network range)
Conda recipe or xenon command line interface.

# Build

Update version number in xenon-cli/meta.yaml.

```sh
./gradlew installShadowDist
cd conda
conda install conda-build conda-verify
conda-build -c conda-forge xenon-cli
```

# Upload

```sh
conda install anaconda-client
anaconda upload --user nlesc <path to xenon-cli*.tar.bz2>
```

# Install

```sh
conda install -c conda-forge -c nlesc xenon-cli
```
