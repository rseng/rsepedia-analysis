The Xenon library has different test suites:

* unit tests
* integration tests - Run the integration tests against docker containers
* fixed client environment tests - Run the integration tests within and against docker containers
* live tests - Runs the integration tests against live systems

# Run unit tests

The unit tests can be run using:
```bash
./gradlew test
```

# Run integration tests against Docker containers

Requirements:
* [docker](https://docs.docker.com/engine/installation/), v1.13 or greater
* [docker-compose](https://docs.docker.com/compose/), v1.10 or greater

The integration tests of Xenon is run against the Docker images.
 
The schedulers and file systems supported by Xenon all have Docker images on [Docker Hub](https://hub.docker.com/r/nlesc/).
See [https://github.com/NLeSC/xenon-docker-images](https://github.com/NLeSC/xenon-docker-images) for the source code of the images. The integration tests will start the Docker containers they need using [Docker compose junit rules](https://github.com/palantir/docker-compose-rule).

The integration tests against Docker containers can be run with:

```bash
./gradlew integrationTest
```

To filter the tests use [Gradle's filtering mechanism](https://docs.gradle.org/3.3/userguide/java_plugin.html#test_filtering).
For example to only run the tests with `ftp` in their package name or class name use:

```bash
./gradlew integrationTest --tests '*ftp*'
```

The tests can be slow, to see which test is running set the `CI` environment variable.
```
# To print test started/passed/etc. event to stdout
CI=1 ./gradlew integrationTest
```

# Run fixed client environment tests

Requirements:
* [docker](https://docs.docker.com/engine/installation/), v1.13 or greater
* [docker-compose](https://docs.docker.com/compose/), v1.10 or greater
* docker server running on localhost, which excludes MacOS and Windows as they run the Docker server in a VM.

Run tests which expect the client environment (like credentials and ssh agent) to be in a fixed state.

The [xenonmiddleware/fixed-client](https://hub.docker.com/r/nlesc/xenon-fixed-client/) Docker image is used as the fixed state client environment.

The fixed client environment tests can be run with:
```bash
# change to root directory of Xenon repository
./src/fixedClientEnvironmentTest/resources/run-fixed-client-environment-test.sh
```

The script will startup the [xenonmiddleware/fixed-client](https://hub.docker.com/r/nlesc/xenon-fixed-client/) Docker container which:
* has the current code mounted
* can communicate with Docker server, so tests can start filesystem/scheduler containers
* can connect to host local ports, so tests can connect to filesystem/scheduler containers
* has the gradle cache mounted, so it does not need to download all plugins and dependencies again
* runs as the current user, so it writes test results and coverage in the Gradle build directory as the current user
* runs the tests by running `./gradlew fixedClientEnvironmentTest` command

# Run live tests

Run tests against a (remote) system like a cluster with Slurm or an sftp server. 

It is the user's responsibility to manage the live system. For example, if you want to test whether certain files 
exist, you need to create them yourself. If you want to run the tests against the live system, you can run the `src/liveTest/resources/scripts/create_symlinks` script.

Supported arguments (`-D<name>=<value>`):
* `xenon.scheduler`, name of scheduler
* `xenon.filesystem`, name of file system
* `xenon.scheduler.location`, location of scheduler
* `xenon.scheduler.isEmbedded`, whether scheduler is embedded, when set is true
* `xenon.scheduler.supportsInteractive`, whether scheduler supports interactive jobs, when set is true
* `xenon.filesystem.location`, location of filesystem
* `xenon.username`, username for location
* `xenon.password`, password for username
* `xenon.certfile`, path to certificate file
* `xenon.passphrase`, passphrase for certificate file
* `xenon.basedir`, path at location where `create_symlinks` script was run, will get combined with filesystem.getEntryPath() to form absolute path
* `xenon.separator`, overwrite the default OS path separator

Run examples
```bash
# slurm on das5 with default credentials
./gradlew liveTest -Dxenon.scheduler=slurm -Dxenon.scheduler.location=das5.vu.nl
# slurm on das5 with username/password
./gradlew liveTest -Dxenon.scheduler=slurm -Dxenon.scheduler.location=das5.vu.nl -Dxenon.username=username -Dxenon.password=password
# slurm on das5 with certificate file
./gradlew liveTest -Dxenon.scheduler=slurm -Dxenon.scheduler.location=das5.vu.nl -Dxenon.username=username -Dxenon.certfile=pathtocertfile [ -Dxenon.passphrase=passphrase ]
# slurm on das5 with default credentials and a custom property
./gradlew liveTest -Dxenon.scheduler=slurm -Dxenon.scheduler.location=das5.vu.nl -Dxenon.adaptors.slurm.strictHostKeyChecking=false
# sftp on localhost:10022
./gradlew liveTest -Dxenon.filesystem=sftp -Dxenon.filesystem.location=localhost:10022  -Dxenon.username=xenon -Dxenon.password=javagat -Dxenon.adaptors.file.sftp.strictHostKeyChecking=false -Dxenon.adaptors.file.sftp.loadKnownHosts=false
# local filesystem and scheduler
./gradlew liveTest -Dxenon.scheduler=local -Dxenon.filesystem=file -Dxenon.filesystem.location=$PWD -Dxenon.username=$USERNAME -Dxenon.basedir=$PWD -Dxenon.scheduler.supportsInteractive=1 -Dxenon.scheduler.isEmbedded=1 -Dxenon.scheduler.location=$PWD -Dxenon.scheduler.workdir=$PWD
```
1. Write documentation on how to start writing a new adaptor
1. Add more adaptors, for example for:
    - SWIFT
    - [Hadoop](https://en.wikipedia.org/wiki/Apache_Hadoop) HDFS (almost done), YARN, see [hadoop repo](https://github.com/xenon-middleware/xenon-adaptors-hadoop)
    - [Azure-Batch](https://azure.microsoft.com/en-us/services/batch/) [#393](https://github.com/xenon-middleware/xenon/issues/393)
    - [Amazon-Batch](https://aws.amazon.com/batch/) [#415](https://github.com/xenon-middleware/xenon/issues/415)
    - GridFTP, see [grid repo](https://github.com/xenon-middleware/xenon-adaptors-grid)
    - glite
    - [iRods](https://irods.org/) [#597](https://github.com/xenon-middleware/xenon/issues/597)
Xenon 3.1.0
-----------

This is release 3.1.0 of Xenon.

Notable changes compared to v3.0.4:
-----------------------------------

- improved error state returned by getJobStatusses to show the difference between no job and connection loss.
- added support for slurm 19
- added heartbeat in the SSH adaptor quickly detect lost connections.
- moved to version 2.4.0 of SSHD 

Notable changes compared to v3.0.3:
-----------------------------------

- fixed bug in FTP adaptor when listing the root directory
- fixed bug in SFTP adaptor when listing the root directory

Notable changes compared to v3.0.2:
-----------------------------------

- fixed the FTP adaptor which lost bytes due to being in ASCII mode 
- fixed the numbering in this changelog

Notable changes compared to v3.0.1:
-----------------------------------

- updated jaxb dependencies to prevent illegal access warning from JVM
- removed leftover debug print in webdav 

Notable changes compared to v3.0.0:
-----------------------------------

- fixed minute delay on shutdown when SSH adaptor was used (discussed as part of #654).

Notable changes compared to v2.6.2:
-----------------------------------

- Moved adaptors with large dependencies (such as S3 and HDFS) into separate libraries, resulting in a much smaller "core" distribution.
- Changed JobDescription to a tasks+cores model, instead of nodes+processes+thread (#625).   
- Remove the JOB_OPTIONS hack from JobDescription (#629 and #628)
- Added support for memory requirements and job name in JobDescription (#562 and #609)
- Added an adaptor for the at scheduler (#381)
- Dropped offline support (#649)
- Require Java 11 or greater (#647)
- Many smaller bugfixes and updates of dependencies. 

Notable changes compared to v2.6.1:
-----------------------------------

- added support for temp space in JobDescription. 
- added support stdout, stderr and stdin to Torque.
- fixed several unit tests that failed on OSX

Notable changes compared to v2.6.0:
-----------------------------------

- fixed hashCode and equals of JobDescription

Notable changes compared to v2.5.0:
-----------------------------------

- added support for scheduler specific arguments to JobDescription
- fixed specification of runtime limit in gridengine adaptor 

Notable changes compared to v2.4.1:
-----------------------------------

- added equals to KeytabCredential (#615)
- added getSupportedCrenentials to AdaptorDescription (#595)
- clarified description of JobState.getState (#596)

Notable changes compared to v2.4.0:
-----------------------------------

- fixed JobDescription equals, hashCode and toString (#612)
- fixed slurm adaptors status retrieval of finished jobs (#613)
- fixed slurm adaptors parsing of scontrol output on pre 17 slurm versions

Notable changes compared to v2.3.0:
-----------------------------------

- added name to job description and job status (#609)
- added max memory to job description (#562)
- added threads per process to job description

Notable changes compared to v2.2.0:
-----------------------------------

- added an HDFS filesystem adaptor 
- fixed bug in GridEngineSchedulers for complex configurations of number of slots per node
- various code cleanups, etc.

Notable changes compared to v2.1.0:
-----------------------------------

- extended CredentialMap to retrieve all keys 
- removed logback config from jar 
- fixed bug in handling workdir of Local and TorqueSchedulers
- many small bugfixes, additional tests, etc.


Notable changes compared to v2.0.0:
-----------------------------------

- added getCredential to Scheduler and FileSystem 
- fixed a bug in equals of CredentialMap
- added proper check of supported credential types in adaptors 
- many small bugfixes, additional tests, etc.

Notable changes compared to v1.2.3:
-----------------------------------

- complete overhaul of public API, which should increase ease-of-use significantly. 
- complete overhaul of integration test framework, which should improve performance and make it easier to test against different versions of the same middleware.
- complete overhaul of implementation, which should make implementing adaptors much more straightforward.
- replaced Jsch with Apache SSHD in the SSH and SFTP adaptors
- replaced Apache Jackrabbit with Sardine in the Webdav adaptor. 
- added an S3 filesystem adaptor. 

Notable changes compared to v1.2.2:
-----------------------------------

- fixed various issues flagged by sonarqube

Notable changes compared to v1.2.1:
-----------------------------------

- fixed bug in the copy engine that would ignore a copy if source and destination had exactly the same path (even when on different machines).
- added timeout overflow detection in Jobs.waitUntilDone and Jobs.waitUntilRunning.
- added SonarQube code for quality analysis and coverage
- we have a new logo!

Notable changes compared to v1.2.0:
-----------------------------------

- fixed nasty inconsistency in adaptor implementations of waiting for jobs to start or finish.

Notable changes between v1.2.0 and v1.1.0:
------------------------------------------

- added support for WebDAV file access.
- added OSX testing in Travis
- added support for Slurm version 15.08.6
- fixed several bugs related to Windows local file system semantics
- many small bugfixes, additional tests, etc. 


Notable changes between v1.1.0 and v1.0.0:
------------------------------------------
 
- added support for FTP file access.
- added support for Torque resource manager.
- added support for Slurm versions 2.6.0, 14.03.0 and 14.11.9-Bull.1.0.
- added option to specify resources in sge adaptor.
- added support for SSH agent and agent proxies
- added a -lot- of unit and integration tests
- javadoc is java 8 compliant.
- the adaptor documentation is now part of the javadoc.
- Xenon releases are now available in jCenter.
- switched from and ant to gradle based build system, this is reflected in the directory structure 
- split unit and integration tests
- added docker support for integration tests
- now using travis-ci for continous integration 
- now using PMD and codacy for code quality
- now using codecov for unit and integration test coverage
- moved examples and tutorial to a separate repository https://github.com/NLeSC/Xenon-examples


What's missing:
---------------
	
The GridFTP adaptor is not considered stable yet. It is not part of this release.

There is no adaptor writing documentation at the moment, nor is the Javadoc complete for the internals methods of the adaptor implementations.

It should be made easier to inspect at runtime which adaptors are available and what properties they support.

We can always use more adaptors, e.g, for SWIFT, HDFS, YARN, Azure-Batch, etc. These are planned for 2.1.0 or later.

We can always use more interfaces, e.g. for clouds (to start and stop VMs). This is planned for 3.0.0.




# Xenon

[![Build Status](https://travis-ci.org/xenon-middleware/xenon.svg?branch=master)](https://travis-ci.org/xenon-middleware/xenon)
[![Build status Windows](https://ci.appveyor.com/api/projects/status/h4l4wn158db23kuf?svg=true)](https://ci.appveyor.com/project/NLeSC/xenon/branch/master)
[![codecov](https://codecov.io/gh/xenon-middleware/xenon/branch/master/graph/badge.svg)](https://codecov.io/gh/xenon-middleware/xenon)
[![SonarCloud](https://sonarcloud.io/api/project_badges/measure?project=nlesc%3AXenon&metric=alert_status)](https://sonarcloud.io/dashboard?id=nlesc%3AXenon)
[![GitHub license](https://img.shields.io/badge/license-Apache--2.0%20-blue.svg)](https://github.com/xenon-middleware/xenon/blob/master/LICENSE)
[![Download](https://api.bintray.com/packages/nlesc/xenon/xenon/images/download.svg)](https://bintray.com/nlesc/xenon/xenon/_latestVersion)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597993.svg)](https://doi.org/10.5281/zenodo.597993)
[![Research Software Directory](https://img.shields.io/badge/rsd-xenon-00a3e3.svg)](https://www.research-software.nl/software/xenon)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3692/badge)](https://bestpractices.coreinfrastructure.org/projects/3692)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)

Copyright 2013-2021 The Netherlands eScience Center

## What problem does Xenon solve?

Many applications use remote storage and compute resources. To do so, they need
to include code to interact with the scheduling systems and file transfer
protocols used on those remote machines.

Unfortunately, many different scheduler systems and file transfer protocols
exist, often with completely different programming interfaces. This makes it
hard for applications to switch to a different system or support multiple
remote systems simultaneously.

Xenon solves this problem by providing a single programming interface to many
different types of remote resources, allowing applications to switch without
changing a single line of code.

![Xenon abstraction](/docs/images/readme-xenon-abstraction.svg.png "Xenon abstraction")

## How does Xenon work?

Xenon is an abstraction layer that sits between your application and the (usually remote)
resource it uses. Xenon is written in Java, but is also accessible from other
languages (e.g. Python) through its gRPC interface and via the command line.

![Xenon API](/docs/images/readme-xenon-api.svg.png "Xenon API")

### Overview of the Xenon ecosystem of tools

| component | repository |
|---|---|
| Xenon library | https://github.com/xenon-middleware/Xenon |
| Xenon cloud adaptors like s3 | https://github.com/xenon-middleware/xenon-adaptors-cloud |
| Xenon grid adaptors like gridftp| https://github.com/xenon-middleware/xenon-adaptors-grid |
| Xenon hadoop adaptors like hdfs | https://github.com/xenon-middleware/xenon-adaptors-hadoop |
| gRPC extension for Xenon | https://github.com/xenon-middleware/xenon-grpc |
| command line interface to Xenon | https://github.com/xenon-middleware/xenon-cli |
| Python API for Xenon | https://github.com/xenon-middleware/pyxenon |
| Docker images | https://github.com/xenon-middleware/xenon-docker-images |

## Supported middleware

Xenon currently supports the following file access mechanisms:

- ``file`` (local file manipulation)
- ``ftp``
- ``sftp``
- ``webdav``
- ``s3``
- ``hdfs``

Xenon currently supports the following job submission mechanisms:

- ``local``
- ``ssh``
- ``at``
- ``gridengine``
- ``slurm``
- ``torque``  

See the [roadmap](/ROADMAP.md) for the planned extensions.

## Adding Xenon as a dependency to your project

Follow the instructions from [bintray.com](https://bintray.com/nlesc/xenon/xenon) to include Xenon as a 
dependency for Gradle, Maven, SBT, or Leiningen projects, e.g. Gradle:

```gradle
	allprojects {
		repositories {
			...
			jcenter()
		}
	}
```

and 

```gradle
	dependencies {
	        compile group: 'nl.esciencecenter.xenon', name: 'xenon', version: '3.1.0'
	}

```

This will give the core adaptors to get cloud, grid and hadoop adaptors add the following dependencies:
```gradle
    compile group: 'nl.esciencecenter.xenon.adaptors', name: 'xenon-adaptors-cloud', version: '3.0.2'
    compile group: 'nl.esciencecenter.xenon.adaptors', name: 'xenon-adaptors-grid', version: '3.0.0'
    compile group: 'nl.esciencecenter.xenon.adaptors', name: 'xenon-adaptors-hadoop', version: '3.0.0'
```

## Simple examples

Here are some examples of basic operations you can perform with Xenon: 

### Copying a file from a local filesystem to a remote filesystem

```java
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.credentials.PasswordCredential;
import nl.esciencecenter.xenon.filesystems.CopyMode;
import nl.esciencecenter.xenon.filesystems.CopyStatus;
import nl.esciencecenter.xenon.filesystems.FileSystem;
import nl.esciencecenter.xenon.filesystems.Path;

public class CopyFileLocalToSftpAbsolutePaths {

    public static void main(String[] args) throws Exception {

        // Use the file system adaptors to create file system representations; the remote file system
        // requires credentials, so we need to create those too.
        //
        // Assume the remote system is actually just a Docker container (e.g.
        // https://hub.docker.com/r/xenonmiddleware/ssh/), accessible via
        // port 10022 on localhost
        String location = "localhost:10022";
        String username = "xenon";
        char[] password = "javagat".toCharArray();
        PasswordCredential credential = new PasswordCredential(username, password);
        FileSystem localFileSystem = FileSystem.create("file");
        FileSystem remoteFileSystem = FileSystem.create("sftp", location, credential);

        // create Paths for the source and destination files, using absolute paths
        Path sourceFile = new Path("/etc/passwd");
        Path destFile = new Path("/tmp/password");

        // create the destination file only if the destination path doesn't exist yet
        CopyMode mode = CopyMode.CREATE;
        boolean recursive = false;

        // perform the copy and wait 1000 ms for the successful or otherwise
        // completion of the operation
        String copyId = localFileSystem.copy(sourceFile, remoteFileSystem, destFile, mode, recursive);
        long timeoutMilliSecs = 1000;
        CopyStatus copyStatus = localFileSystem.waitUntilDone(copyId, timeoutMilliSecs);

        // throw any exceptions
        XenonException copyException = copyStatus.getException();
        if (copyException != null) {
          throw copyException;
        }
    }
}
```

### Submitting a job

The following code performs a wordcount of a file residing on a remote machine: 

```java 
import nl.esciencecenter.xenon.credentials.PasswordCredential;
import nl.esciencecenter.xenon.schedulers.JobDescription;
import nl.esciencecenter.xenon.schedulers.JobStatus;
import nl.esciencecenter.xenon.schedulers.Scheduler;

public class SlurmSubmitWordCountJob {

    public static void main(String[] args) throws Exception {

        // Assume the remote system is actually just a Docker container (e.g.
        // https://hub.docker.com/r/xenonmiddleware/slurm/), accessible to user 'xenon' via
        // port 10022 on localhost, using password 'javagat'
        String location = "localhost:10022";
        String username = "xenon";
        char[] password = "javagat".toCharArray();
        PasswordCredential credential = new PasswordCredential(username, password);

        // create the SLURM scheduler representation
        Scheduler scheduler = Scheduler.create("slurm", location, credential);

        JobDescription description = new JobDescription();
        description.setExecutable("/usr/bin/wc");
        description.setArguments("-l", "/etc/passwd");
        description.setStdout("/tmp/wc.stdout.txt");

        // submit the job
        String jobId = scheduler.submitBatchJob(description);

        long WAIT_INDEFINITELY = 0;
        JobStatus jobStatus = scheduler.waitUntilDone(jobId, WAIT_INDEFINITELY);

        // print any exceptions
        Exception jobException = jobStatus.getException();
        if (jobException != null)  {
          throw jobException;
        }

    }
}
```

The output of the job will be written to ``/tmp/wc.stdout.txt`` file in the ``xenonmiddleware/slurm`` Docker container.

For more examples, see the tutorial at [Read The Docs](http://xenonrse2017.readthedocs.io/).

## Documentation

Xenon's JavaDoc is available online at <http://xenon-middleware.github.io/xenon/>.

## Documentation for maintainers

For developers of Xenon itself 

* see [RELEASE.md](RELEASE.md) how to perform a release.
* see [ADAPTOR_DEVELOPMENT.md](ADAPTOR_DEVELOPMENT.md) how to write a new adaptor.

## Legal

The Xenon library is copyrighted by the Netherlands eScience Center and released
under the Apache License, Version 2.0. A copy of the license may be obtained
from [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

Xenon uses several third-party libraries that have their own (permissive, open 
source) licenses. See the file [legal/README.md](legal/README.md) for an overview.

# Adding an adaptor

> This documentation is out of date, for now use https://github.com/xenon-middleware/xenon-adaptors-cloud as an example.

To add an adaptor the following steps must be performed:
* source code
* dependencies (optional)
* unit tests
* integration tests
* Docker based server for adaptor to test against in integration tests
* registration in Xenon engine
* registration in build system

For a new file adaptor use `webdav` as an example. 
For a new job adaptor use `slurm` as an example. 

Adding an adaptor can be completed by adding/changing files in the following locations:
1. Source code in `src/main/java/nl/esciencecenter/xenon/adaptors/<adaptor name>`.
2. Specify dependencies of adaptor in `build.gradle`.
3. Unit tests in `src/test/java/nl/esciencecenter/xenon/adaptors/<adaptor name>`.
4. Register adaptor in `src/main/java/nl/esciencecenter/xenon/engine/XenonEngine.java:loadAdaptors()`
5. Integration tests
  1. Create a Dockerfile in `src/integrationTest/docker/xenon-<adaptor-name>` for a server of the adaptor to test against
  2. Register the Docker image in `src/integrationTest/docker/docker-compose.yml` and  `gradle/docker.gradle`
  3. Add the Docker container credentials/location/configuration to `src/integrationTest/docker/xenon.test.properties.docker`
  4. Create an integration test in `src/integrationTest/java/esciencecenter/xenon/adaptors/<adaptor name>/`
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [xenon@esciencecenter.nl](mailto:xenon@esciencecenter.nl). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/NLeSC/Xenon/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/NLeSC/Xenon/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of concensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing unit tests still work by running ``./gradlew test``;
1. make sure that the existing integration tests still work by running ``./gradlew check``; and ``./src/fixedClientEnvironmentTest/resources/run-fixed-client-environment-test.sh``;
1. add your own unit tests and integration tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the Xenon repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request (have a look at some of our old pull requests to see how this works, for example [#294](https://github.com/NLeSC/Xenon/pull/294)).


# Creating a release

Step-by-step guide on creating a Xenon release

## 1. Up the version number and add changelog

To create a release, version numbers should be updated in:

- README.md
- gradle/common.gradle
- CHANGELOG.md

    Also, add a section to the CHANGELOG.md with notable changes in this version.

- CITATION.cff (see section below)

## 2. Update the site

To update Xenon version in `docs/_data/version.yml` and update Javadocs inside `docs/` directory run:

```bash
./gradlew publishSite
```

The site is a [Jekyll](https://jekyllrb.com/) powered site and hosted by GitHub pages at http://nlesc.github.io/Xenon/

## 2.1 Update the CITATION.cff and generate the Zenodo metadata

The ``CITATION.cff`` contains the information on how to cite Xenon in scientific
applications. Update the authors list if needed, and change the
``date-released`` and ``version`` fields to reflect the release that you are
about to make.

NOTE: the CITATION.cff contains the concept DOI of xenon (refering to all versions), so this does not need to be changed. 

After updating the citation metadata, use [``cffconvert``](https://pypi.org/project/cffconvert/)
to generate the metadata as used by Zenodo. Install ``cffconvert`` as follows:

```bash
# install cffconvert in user space from PyPI
pip3 install --user cffconvert

# change directory to xenon root dir (if needed)
cd <directory where the CITATION.cff lives>

# check if the CITATION.cff file is valid (if there is no output, that 
# means it's all good)
cffconvert --validate

# generate the zenodo file by (note the dot in the file name):
cffconvert --outputformat zenodo --ignore-suspect-keys > .zenodo.json
```


## 3. Commit the changes

Next, commit all changes you have made to the master branch. If you check with  

    git status

you should get something like this:

    On branch master
    Your branch is up-to-date with 'origin/master'.
    Changes not staged for commit:

	modified:   .zenodo.json
	modified:   CITATION.cff
	modified:   CHANGELOG.md
	modified:   README.md
	modified:   gradle/common.gradle
    modified:   docs/_config.yml

    Untracked files:
    (use "git add <file>..." to include in what will be committed)

	docs/versions/3.0.1/

Add and commit these files using `git add` and `git commit` and `git push`.

## 4. Create a GitHub release

On GitHub, go to the releases tab and draft a new release.

The tag and title should be set to the version.

The release description should contain a small text to explain what Xenon is and the part of the CHANGELOG.md which pertains to this version.

## 5. Check if DOI is created in Zenodo

Zenodo is linked to the Xenon github, so when a release is created, a DOI 
will be created automatically. Click the DOI badge on the github page to check 
this.

Check that the authors and license of the Zenodo entry are correct.
If not, then correct, save and publish the Zenodo entry.

### 6. Add jars to bintray

To add the release to bintray, do the following: 

```bash
export BINTRAY_USER=<your bintray username>
export BINTRAY_KEY=<your bintray API key>
./gradlew publishToBintray
```

This should create the new release on bintray and upload all necessary data, jars, etc.

On https://bintray.com/nlesc/xenon check that xenon* packages have been updated.

The latest version tag usually takes a few minutes to update.

### Alternative manual bintray step

Note: this step should not be needed if the automated bintray publishing works! It is only here for reference.

To add the necessary jar and pom files to bintray manually, it is easiest to ensure you 
have them locally first. By calling: 

    ./gradlew publishToMavenLocal

gradle should put the required files in:

    ~/.m2/repository/nl/esciencecenter/xenon/xenon/<version>/

Next goto: 

    https://bintray.com/nlesc/xenon/xenon

Click on 'new version' and insert <version> as name. Next, in the overview page 
click on the version number and on `UI` (in the version files box). Bintray will 
ask for the 'Target Repository Path', which should be: 

    /nl/esciencecenter/xenon/xenon/<version>

Click on the `click to add files` button and add the files that where generated in 
the .m2 repository. Click on `save` and then on `publish` to publish the version.
 
### 7. Update applications using Xenon.

Update related repos such as Xenon-examples, pyXenon, xenon-cli, etc

And finally celebrate.
# Serve site

Serve site on [http://localhost:4000/xenon/](http://localhost:4000/xenon/) with the following command
```
docker run --rm --label=jekyll --volume=$(pwd):/srv/jekyll \
-i -t -p 127.0.0.1:4000:4000 jekyll/jekyll:pages jekyll serve
```

# Updating current version

1. In `_data/version.yml` update the value of the `current` key to new current version, do manually or by using `./gradlew versionSite`
2. Add version specific artifacts like javadoc to `versions/<version>/`, do manually or by using `./gradlew copyJavadoc`
3. Commit and push
![logo](images/NLeSC_Xenon_logo.png "Xenon Logo")

Xenon
=====

Copyright 2013 The Netherlands eScience Center

Author: Jason Maassen (<J.Maassen@esciencecenter.nl>)

Version: Userguide v1.0, Xenon v1.0

Last modified: 11 September 2013


Copyrights & Disclaimers
------------------------

Xenon is copyrighted by the Netherlands eScience Center and 
releases under the Apache License, Version 2.0.

See the "LICENSE" and "NOTICE" files in the Xenon distribution for
more information. 

For more information on the Netherlands eScience Center see:

<http://www.esciencecenter.nl> 

The Xenon project web site can be found at:

<https://github.com/NLeSC/Xenon>.


Third party libraries
---------------------

This product includes the SLF4J library, which is Copyright 
(c) 2004-2013 QOS.ch See "notices/LICENSE.slf4j.txt" for the licence
information of the SLF4J library.

This product includes the JSch library, which is Copyright 
(c) 2002-2012 Atsuhiko Yamanaka, JCraft,Inc. 
See "notices/LICENSE.jsch.txt" for the licence information of the 
JSch library.

This product includes the Logback library, which is Copyright 
(c) 1999-2012, QOS.ch. See "notices/LICENSE.logback.txt" for the 
licence information of the Logback library.

This product includes the JaCoCo library, which is Copyright
(c) 2009, 2013 Mountainminds GmbH & Co. KG and Contributors. See
"notices/LICENSE.jacoco.txt" for the licence information of the 
JaCoCo library.

This project includes the JUnit library. 
See "notices/LICENSE.junit.txt" for the licence information of the 
JUnit library.

This project includes the Mockito library, which is Copyright 
(c) 2007 Mockito contributors. See "notices/LICENSE.mockito.txt" for
the licence information of the Mockito library.


What is it?
-----------

Xenon is a middleware abstraction library. It provides a simple 
Java programming interface to various pieces of software that can be
used to access distributed compute and storage resources. 


Why Xenon?
------------

Xenon is developed by the Netherlands eScience Center as a support
library for our projects. Several projects develop end-user
applications that require access to distributed compute and storage
resources. Xenon provides a simple API to access those resources,
allowing those applications to be developed more rapidly. The
experience gained during the development of these end-user 
applications is used to improve the Xenon API and implementation. 


Installation
------------

The installation procedure and dependencies of the Xenon library
can be found in the file "INSTALL.md" in the Xenon distribution. 


Design
------

Xenon is designed with extensibility in mind. It uses a modular
and layered design as shown in the figure below:

![Xenon design](images/xenon-design.png "Xenon design.")
	
Xenon consists of three layers, an *interface layer*, an 
*engine layer* and an *adaptor layer*. 

The *interface layer* is used by the application using Xenon. It
contains several specialized interfaces:

- Xenon: this is the main entry point used to retrieve the other
  interfaces. 
- Files: contains functionality related to files, e.g., creation,
  deletion, copying, reading, writing, obtaining directory listings,
  etc. 
- Jobs: contains functionality related to job submission, e.g.,
  submitting, polling status, cancelling, etc. 
- Credentials: contains functionality related to credentials.
  Credentials (such as a username password combination) are often
  needed to gain access to files or to submit jobs. 

The modular design of Xenon allows us to add additional interfaces
in later versions, e.g., a Clouds interface to manage virtual
machines, or a Networks interface to manage bandwidth-on-demand
networks. 

The *adaptor layer* contains the adaptors for the each of the
middlewares that Xenon supports. An *adaptor* offers a middleware
specific implementation for the functionality offered by one of the
interfaces in Xenon.
 
For example, an adaptor may provide an *sftp* specific implementation
of the functions in the Xenon *file interface* (such as *copy* or
*delete*) by translating each of these functions to *sftp* specific
code and commands.

For each interface in Xenon there may be multiple adaptors
translating its functionality to different middlewares. To
distinguises between these adaptors Xenon uses the *scheme* they
support, such as "sftp", "http" or "ssh". There can be only one
adaptor for each scheme. 

The *engine layer* of Xenon contains the "glue" that connects each
interface to the adaptors that implement its functionality. When a
function of the interface layer is invoked, the call will be
forwarded to the engine layer. It is then the responsibility of the
engine layer to forward this call to the right adaptor. 

To perform this selection, the engine layer matches the *scheme* of
the object on which the operation needs to be performed, to the
*schemes* supported by each of the adaptors. When the schemes match,
the adaptor is selected. 


Interfaces and datatypes
------------------------

This section will briefly explain each of the interfaces and related
datatypes. Detailed information about Xenon can be found in the
online Javadoc at: 

<http://nlesc.github.io/xenon/javadoc/>

### Package Structure ##

The Xenon API uses the following package structure:

- `nl.esciencecenter.xenon` Entry point into Xenon.
- `nl.esciencecenter.xenon.credentials` Credential interface.
- `nl.esciencecenter.xenon.files`  Files interface.
- `nl.esciencecenter.xenon.jobs`  Jobs interface.
- `nl.esciencecenter.xenon.util`  Various utilty classes.

We will now briefly describe the most important classes and
interfaces of these packages.

### Getting started ###

The [`nl.esciencecenter.xenon`][1] package contains the entry
point into the Xenon library: [__XenonFactory__][2]

    public class XenonFactory {
       public static Xenon newXenon(Map<String,String> properties) 
       public static void endXenon(Xenon xenon) 
       public static void endAll()
    }

The __newXenon__ method can be used to create a new __Xenon__
instance, while the __endXenon__ method can be used to release the
Xenon instance once it is no longer needed. It is important to end
the Xenon when it is no longer needed, as this allows it to release 
any resources it has obtained. 

When creating an Xenon using __newXenon__, the _properties_
parameter can be used to configure the Xenon instance. If no
configuration is necessary, `null` can be used. 

Properties consist of a set of key-value pairs. In Xenon all keys
__must__ start with "xenon.". To configure the adaptors,
properties of the form "xenon.adaptors.(name).(property)" can be
used, where "(name)" is the name of the adaptor (for example "local"
or "ssh") and "(property)" is the name of the property to be 
configured. Note that this name can be futher qualified, for example 
"xenon.adaptors.local.a.b.c". The available properties can be found
in the documentation of the individual adaptors (see Appendix A). 

A call to __newXenon__ will return an [__Xenon__][3]:

    public interface Xenon {
        Files files()
        Jobs jobs()
        Credentials credentials()
        Map<String,String> getProperties()
        AdaptorStatus getAdaptorStatus(String adaptorName)
        AdaptorStatus[] getAdaptorStatuses()
    }

The __files__, __jobs__ and __credentials__ methods in this
interface can be used to retrieve various interfaces that the
Xenon library offers. They will be described in more detail below. 

The __getProperties__ method can be used to retrieve the properties
used when the Xenon was created. Most objects created by Xenon
contain such a __getProperties__ method. For brevity, we will not
explain these further.

The __getAdaptorStatus__ method can be used to retrieve information
about the adaptors. This information is returned in an 
[__AdaptorStatus__][4]:

    public interface AdaptorStatus {
        String getName()
        String getDescription()
        String[] getSupportedSchemes()
        XenonPropertyDescription[] getSupportedProperties()
        Map<String, String> getAdaptorSpecificInformation()
    }
    
An __AdaptorStatus__ contains __getName__ to retrieve the name of an
adaptor,  __getDescription__ to get a human readable description of
what functionality it has to offer and __getSupportedSchemes__ to
retrieve a list of the schemes it supports.

The __getSupportedProperties__ method can be used to retrieve a list
of configuration options the adaptor supports. Each returned
[__XenonPropertyDescription__][5] gives a full description of a
single property, including its name (of the form 
"xenon.adaptors.(name).(property)"), the expected type of its 
value, a human readable description of its purpose, etc. More
information on the supported properties can be found in Appendix A.

Finally, __getAdaptorSpecificInformation__ can be used to retrieve
status information from the adaptor. Each key contains a property of
the form described above. 

### Credentials interface ###

The [`nl.esciencecenter.xenon.credentials`][6] package contains the
[__Credentials__][7] interface of Xenon:

    public interface Credentials {
        Credential newCertificateCredential(String scheme, String keyfile, String certfile, String username, 
            char [] password, Map<String,String> properties) 

        Credential newPasswordCredential(String scheme, String username, char [] password, Map<String,String> properties)
        Credential getDefaultCredential(String scheme)
        void close(Credential credential)
    }

The __Credentials__ interface contains various methods for creating
credentials, based on certificates or passwords. For each method,
the desired _scheme_ needs to be provided as a parameter (for example,
"ssh" or "sftp"). This allows Xenon to forward the call to the
correct adaptor. Note that some types of credentials may not be
supported by all adaptors. An exception will be thrown when an
unsupported __new**Credential__ methods is invoked. 

Additional configuration can also be provides using the _properties_
parameter, which use the same form as described in the
_Xenon factory and interface_ section above. If no additional
configuration is needed, `null` can be used. The 
__getDefaultCredential__ method returns the default credential for
 the given scheme. All adaptors are guarenteed to support this
method. 

All __new**Credential__ methods return a [__Credential__][13] that contains the following 
methods: 

    public interface Credential {
       String getAdaptorName()
       Map<String,String> getProperties()
    }

The __getAdaptorName__ method can be used to retrieve the name of
the adaptor that created the credential. Many adaptor specific
objects returned by Xenon contain this method. For brevity we will
not explain this further.

When a __Credential__ is no longer used, it __must__ be closed using
__close__. This releases any resources held by the __Credential__.
The __isOpen__ method can be used to check if a __Credential__ is
open or closed. 

### Files interface ###

The [`nl.esciencecenter.xenon.files`][8] package contains the
[__Files__][9] interface of Xenon. For readability we will split
the explanation of __Files__ into several parts:

    public interface Files {
       FileSystem newFileSystem(String scheme, String location, Credential credential, Map<String,String> properties)
       void close(FileSystem filesystem) 
       boolean isOpen(FileSystem filesystem)
       // ... more follows
    }

The __Files__ interface contains several method for creating and
closing a [__FileSystem__][10]. A __FileSystem__ provides an
abstraction for a (possibly remote) file system. 

To create a __FileSystem__ the __newFileSystem__ method can be used.
As before, the desired __scheme__ must be provided as a parameter.
In addition, the _location_ parameter provides information on the 
location of the file system using an adaptor specific string. For
local file systems, the location must contain the root of the file
system to access, such as "/" on Linux or "C:" on Windows. For remote 
file systems, the location typically contains the host name of the
machine to connect to. The exact format of accepted location strings
can be found in the adaptor documentation. 

The following are all valid combinations of file system schemes and
locations:

    "file","/"                   connect to local Linux file system
    "file","C:"                  connect to local Windows C: drive
    "sftp","example.com"         connect to example.com using sftp 
    "sftp","test@example.com:44" connect to example.com using sftp 
                                 on port 44 with "test" as user name.

The __newFileSystem__ method also has a _credential_ parameter to
provide the credential needed to access the file system. If this
parameter is set to `null` the default credentials will be used for
the scheme. The _properties_ parameter can be used to provide
additional configuration properties. Again, `null` can be used if no
additional configuration is required. The returned __FileSystem__
contains the following methods:

    public interface FileSystem {
        String getScheme()
        String getLocation()
        Path getEntryPath()
    }

The __getScheme__ and __getLocation__ methods returns the scheme and
location strings used to create the __FileSystem__. The
__getEntryPath__ method returns the _path at which the file system
was entered_. For example, when accessing a file system using "sftp"
it is customary (but not manditory) to enter the file system at the
users' home directory. Therefore, the entry path of the 
__FileSystem__ will be similar to "/home/(username)". For local file
systems the entry path is typically set to the root of the file
system (such as "/" or "C:").

When a __FileSystem__ is no longer used, it __must__ be closed using
__close__. This releases any resources held by the __FileSystem__.
The __isOpen__ method can be used to check if a __FileSystem__ is
open or closed. Once a __FileSystem__ is created, it can be used to
access files: 

    public interface Files {
       Path newPath(FileSystem filesystem, RelativePath location) 
       void createFile(Path path)
       void createDirectories(Path dir)
       void createDirectory(Path dir)
       boolean exists(Path path)
       void delete(Path path)
       FileAttributes getAttributes(Path path)
       // ... more follows
    }

The __newPath__ method can be used to create a new [__Path__][11].
A __Path__ represents a path on a specific __FileSystem__. This path
does not necessarily exists. To create an __Path__, both the target
__FileSystem__ and a [__RelativePath__][12] are needed. A
__RelativePath__ contains a sequence of strings separated using a
special _separator_ character, and is used to identify a location
on a file system (for example "/tmp/dir"). __RelativePath__ contains
many utility methods for manipulating these string sequences. The
details can be found in the Javadoc of [__RelativePath__][12].

__Files__ contains several methods to create and delete files and
directories. When creating files and directories Xenon checks if
the target already exists. If so, an exception will be thrown.
Similary, an exception is thrown when attempting to delete
non-existing file or a directory that is not empty. The __exists__
method can be used to check if a path exists.

Using the __getAttributes__ method the attributes of a file can be
retrieved. The returned [__FileAttributes__][14] contains information
on the type of file (regular file, directory, link, etc), it size,
creation time, access rights, etc. 

To list directories, the following methods are available:

    public interface Files {
       DirectoryStream<Path> newDirectoryStream(Path dir)
       DirectoryStream<PathAttributesPair> newAttributesDirectoryStream(Path dir)  
       // ... more follows
    }

Both __newDirectoryStream__ and __newAttributesDirectoryStream__
return a [__DirectoryStream__][15] which can be used to iterate over
the contents of a directory. For the latter, the __FileAttributes__ 
for each of the files are also included. alternatively, these methods
are also available with an extra _filter_ parameter, which can be
used to filter the stream in advance.

To read or write files, the following methods are available:

    public interface Files {
       InputStream newInputStream(Path path)
       OutputStream newOutputStream(Path path, OpenOption... options)
    }

Using these methods, an __InputStream__ can be created to read a
file, and an __OutputStream__ can be created to write a file. The
__newOutputStream__ method requires a _options_ parameter to specify
how the file should be opened for writing (for example, should the
data be append or should the file be truncated first). These options
are describe in more detail in the Javadoc.

To copy files, the following methods are available:

    public interface Files {
       Copy copy(Path source, Path target, CopyOption... options)
       CopyStatus getCopyStatus(Copy copy)
       CopyStatus cancelCopy(Copy copy)
    }

The __copy__ method supports various copy operations such as a
regular copy, a resume or an append. The _options_ parameter can be
used to specify the desired operation. Normally, __copy__ performs
its operation _synchronously_, that is, the call blocks until the
copy is completed. However, _asynchronous_ operations are also
supported by providing the option [__CopyOption.ASYNCHRONOUS__][17].
In that case a [__Copy__][16] object is returned that can be used 
to retrieve the status of the copy (using __getCopyStatus__) or
cancel it (using __cancelCopy__). The details of the available copy
operations can be found in the Javadoc of [__CopyOption__][17].

### Jobs interface ###

The [`nl.esciencecenter.xenon.job`][18] package contains the
[__Jobs__][19] interface of Xenon. For readability we will split
the explanation of __Jobs__ into several parts:

    public interface Jobs {
        Scheduler newScheduler(String scheme, String location, Credential credential, Map<String,String> properties)
        void close(Scheduler scheduler)
        boolean isOpen(Scheduler scheduler)
        // ... more follows
    }

The __Jobs__ interface contains the __newScheduler__ method that can
be used to create a [__Scheduler__][20]. A __Scheduler__ provides an
abstraction for a (possibly remote) scheduler that can be used to
run jobs. The __newScheduler__ method has __scheme__ and __location__
parameters that specify how to access the scheduler. As with
__newFileSystem__ the __location__ is adaptor specific. To access
the local scheduler, passing `null` or an empty string is sufficient.
To access remote schedulers, the location typically contains the
host name of the machine to connect to. The exact format of accepted
location strings can be found in the adaptor documentation.

The following are valid examples of scheduler schemes and locations:

    "local",""                     the local scheduler 
    "ssh","example.com"            connect to a remote scheduler at example.com using SSH
    "slurm",""                     connect to a local slurm scheduler
    "slurm","test@example.com:44"  connect to a remote slurm scheduler at example.com via SSH on port 44 with user "test".

When a __Scheduler__ is no longer used, is __must__ be closed using
the __close__ method. The __isOpen__ method can be use to check if a
__Scheduler__ is open or closed. A __Scheduler__ contains the
following:

    public interface Scheduler {
        String[] getQueueNames()
        boolean isOnline()
        boolean supportsInteractive()
        boolean supportsBatch()
    }

Each __Scheduler__ contains one or more queues to which jobs can be
submitted. Each queue has a name that is unique to the __Scheduler__.
The __getQueueNames__ method can be used to retrieve all queue names. 

The __isOnline__ method can be used to determine if the __Scheduler__
is an _online scheduler_ or an _offline scheduler_. Online schedulers
need to remain active for their jobs to run. Closing an online 
scheduler will kill all jobs that were submitted to it. Offline
schedulers do not need to remains active for their jobs to run. A
submitted job will typically be handed over to some external server
that will manage the job for the rest of its lifetime.

The __supportsInteractive__ and __supportsBatch__ method can be use
to check if the __Scheduler__ supports interactive and/or batch jobs.
Interactive jobs are jobs where the user gets direct control over
the standard streams of the job (the _stdin_, _stdout_ and _stderr_
streams). The user __must__ retrieve these streams using the
__getStreams__ method in __Jobs__ and then provide input and output,
or close the streams. Failing to do so may cause the job to block
indefinately. Batch jobs are jobs where the standard streams are
redirected from and to files. The location of these files must be
set before the job is started, as will be explained below.

Once a __Scheduler__ is created, __Jobs__ contains several methods
to retrieve information about the __Scheduler__:

    public interface Jobs {
        String getDefaultQueueName(Scheduler scheduler)
        QueueStatus getQueueStatus(Scheduler scheduler, String queueName)
        QueueStatus[] getQueueStatuses(Scheduler scheduler, String... queueNames).
        Job[] getJobs(Scheduler scheduler, String... queueNames)
        // ... more follows
    }

The __getQueueStatuses__ method can be used to retrieve information
about a queue. If no queue names are provided as a parameter,
information on all queues in the scheduler will be returned. Using
the __getDefaultQueueName__ the default queue can be retrieved for
the __Scheduler__. The __getJobs__ method can be used to retrieve
information on all jobs in a queue. Note that this may also include
jobs from other users.

To submit and manage jobs, the __Jobs__ interface contains the
following methods:

    public interface Jobs {
        Job submitJob(Scheduler scheduler, JobDescription description)
        Streams getStreams(Job job)
        JobStatus getJobStatus(Job job)
        JobStatus[] getJobStatuses(Job... jobs)
        JobStatus waitUntilRunning(Job job, long timeout)
        JobStatus waitUntilDone(Job job, long timeout)
        JobStatus cancelJob(Job job)
    }    

The __submitJob__ method can be used to submit a job to a
__Scheduler__. A [__JobDescription__][21] must be provided as
parameter. A __JobDescription__ contains all necessary information
on how to start the job, for example, the location of the executable,
any command line arguments that are required, the working directory,
if the job is an interactive of batch job, the location of the files
for stream redirection (in case of a batch job), etc. See the Javadoc
for details of the __JobDescription__.

Once a job is submitted, a [__Job__][22] object is returned that can
be used with __getJobStatus__ to retrieve the status of the job, and
with __cancelJob__ to cancel it. This __Job__ contains the following:

    public interface Job {
        JobDescription getJobDescription()
        Scheduler getScheduler()
        String getIdentifier()
        boolean isInteractive()
        boolean isOnline()
    } 

Besides methods for retrieveing the __JobDescription__ and
__Scheduler__ that created it, each __Job__ also contains the
__isInteractive__ method to determine if the __Job__ is interactive,
and the __isOnline__ method to determine if the job is running on an
_online scheduler_ (explained above).
 
After submitting a job, __waitUntilRunning__ can be used to wait
until a job is no longer waiting in the queue and __waitUntilDone__
can be used to wait until the job has finished.  

For all methods returning a [__JobStatus__][23], the following rule
applies: after a job has finished, the status is only guarenteed to
be returned _once_. Any subsequent calls to a method that returns a 
__JobStatus__ _may_ throw an exception stating that the job does not
exist. Some adaptors may return a result however.  


### Utilities classes ###

The [`nl.esciencecenter.xenon.util`][25] package contains various
utility classes. The main entry points are __Utils__, __Sandbox__
and __JavaJobDescription__.

In [__Utils__][42] various utility methods can be found that make it
easier to use Xenon. Many methods provide simple shortcuts to
often used code constructs. Some examples are shown below:

    public class Utils {
        // Create a new local Scheduler.
        public static Scheduler getLocalScheduler(Jobs jobs)             

        // Create a new Scheduler without Credentials or properties.
        public static Scheduler newScheduler(Jobs jobs, String scheme)   

        // Create a Path that represents the home directory of the current user.
        public static Path getLocalHome(Files files)                     

        // Create a Path that represents the current working directory.
        public static Path getLocalCWD(Files files)                      

        // Convert a String containing a local path into a Path.
        public static Path fromLocalPath(Files files, String path)       

        // Retrieve all local file systems.
        public static FileSystem [] getLocalFileSystems(Files files)     
 
        //  Are we running on a Linux machine ?
        public static boolean isLinux()   

        // Are we running on a Windows machine ?
        public static boolean isWindows()                     

        // Are we running on a OSX machine ?
        public static boolean isOSX()                                   
    }

In addition many methods are provided for reading data from files or
streams to various output targets, writing data to files or streams
from various input sources, recursive copying, recursive deletion,
etc. See the Javadoc of [__Utils__][42] for details.

A [__Sandbox__][43] is a utility class that makes is it easier to
create a (possibly remote) temporary directory and transfer files to
and from this directory. A Sandbox is often used in when submitting 
jobs that require input files and / or produce output files. Sandbox
contains the following methods:

    public class Sandbox {
       Sandbox(Files files, Path root, String sandboxName)
       void addUploadFile(Path src, String dest)
       void addDownloadFile(String src, Path dest)
       void upload(CopyOption... options)
       void download(CopyOption... options)
       void delete()
    }

Creating a Sandbox requires an Xenon __Files__ interface and a
__root__ directory. The Sandbox will then create a temporary
directory __sandboxName__ in __root__. If __sandboxName__ is `null`,
a random name will be generated. Using __addUploadFile__ files can
be added to the upload queue. These files will be transferred to the
Sandbox directory when __upload__ is invoked. Similarly, using
__addDownloadFile__, files can be added to the download queue. They
will be downloaded from the Sandbox directory when __download__ is
invoked. Finally, __delete__ can be used to delete the Sandbox
directory.

A [__JavaJobDescription__][44] is a utility class that makes is it
easier to create a __JobDescription__ for running a Java application.
In addition to the command line arguments used by the application,
Java applications typically require a number of _special_ command
line argument for the Java Virtual Machine (JVM), such as a
_class path_, _system properties_, and _JVM options_. 

The JavaJobDescription class extends the regular JobDescription with
support for these additional arguments. When a Job a submitted to an
Xenon Scheduler that uses a JavaJobDescription, the various types
of command line arguments will be merged automatically into a single
arguments list. See the Javadoc of [__JavaJobDescription__][44] for
details.

Examples
--------

Examples of how to use Xenon can be found in the [examples][26]
directory. We will list the examples here in order of increasing
complexity, and with a short description of each example.

### Initializing Xenon ###

Creating an __Xenon__: 
[CreatingXenon.java][27]

Creating an __Xenon__ with configuration properties: 
[CreatingXenonWithProperties.java][28]

### Creating Credentials ###

Creating a password and default __Credential__:
[CreatingCredential.java][29]

### File Access ###

Creating a local __FileSystem__:
[CreateLocalFileSystem.java][30]

Checking if a local file exists:
[LocalFileExists.java][32]

Creating a __FileSystem__ based on a URI. 
[CreateFileSystem.java][31]

Checking if a (possibly remote) file exists:
[FileExists.java][32]

Listing a directory:
[DirectoryListing.java][33]

Listing the attributes of a file:
[ShowFileAttributes.java][45]

Copying a file:
[CopyFile.java][34]

### Job Submission ###

Creating a __Scheduler__ and retrieving the status of its queues:
[ListQueueStatus.java][36]

Creating a __Scheduler__ and retrieving the jobs:
[ListJobs.java][37]

Listing the status of a Job:
[ListJobStatus.java][38]

Submitting a batch job without output: 
[SubmitSimpleBatchJob.java][39]

Submitting a batch job with output: 
[SubmitBatchJobWithOutput.java][40]

Submitting an interactive job with output: 
[SubmitInteractiveJobWithOutput.java][41]

[1]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/package-summary.html
[2]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/XenonFactory.html
[3]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/Xenon.html
[4]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/AdaptorStatus.html
[5]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/XenonPropertyDescription.html
[6]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/credentials/package-summary.html 
[7]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/credentials/Credentials.html
[8]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/package-summary.html
[9]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/Files.html 
[10]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/FileSystem.html
[11]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/Path.html
[12]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/RelativePath.html
[13]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/credentials/Credential.html
[14]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/FileAttributes.html
[15]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/DirectoryStream.html
[16]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/Copy.html
[17]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/files/CopyOption.html
[18]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/package-summary.html
[19]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/Jobs.html
[20]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/Scheduler.html
[21]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/JobDescription.html
[22]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/Job.html
[23]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/jobs/JobStatus.html
[25]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/utils/package-summary.html
[26]: https://github.com/NLeSC/xenon/tree/develop/examples/src/nl/esciencecenter/xenon/examples
[27]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/CreatingXenon.java
[28]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/CreatingXenonWithProperties.java
[29]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/credentials/CreatingCredential.java
[30]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/CreateLocalFileSystem.java
[31]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/CreateFileSystem.java
[32]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/LocalFileExists.java  
[33]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/FileExists.java
[34]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/DirectoryListing.java
[35]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/CopyFile.java 
[36]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/ListQueueStatus.java
[37]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/ListJobs.java
[38]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/ListJobStatus.java
[39]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/SubmitSimpleBatchJob.java
[40]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/SubmitBatchJobWithOutput.java
[41]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/jobs/SubmitInteractiveJobWithOutput.java
[42]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/utils/Utils.html
[43]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/utils/Sandbox.html
[44]: http://nlesc.github.io/xenon/javadoc/nl/esciencecenter/xenon/utils/JavaJobDescription.html
[45]: https://github.com/NLeSC/xenon/blob/develop/examples/src/nl/esciencecenter/xenon/examples/files/ShowFileAttributes.java

Appendix A: Adaptor Documentation
---------------------------------

This section contains the adaptor documentation which is generated from the information provided by the adaptors themselves.

Xenon currently supports 4 adaptors: local, ssh, gridengine, slurm.

Adaptor: local
--------

The local adaptor implements all functionality with  standard java classes such as java.lang.Process and java.nio.file.Files.

#### Supported schemes: ####
local, file

#### Supported locations: ####
(null), (empty string), /

#### Supported properties: ####


__`xenon.adaptors.local.queue.pollingDelay`__

The polling delay for monitoring running jobs (in milliseconds).

- Expected type: INTEGER

- Default value: 1000

- Valid for: [XENON]


__`xenon.adaptors.local.queue.multi.maxConcurrentJobs`__

The maximum number of concurrent jobs in the multiq..

- Expected type: INTEGER

- Default value: 4

- Valid for: [XENON]



Adaptor: ssh
--------

The SSH adaptor implements all functionality with remove ssh servers.

#### Supported schemes: ####
ssh, sftp

#### Supported locations: ####
[user@]host[:port]

#### Supported properties: ####


__`xenon.adaptors.ssh.autoAddHostKey`__

Automatically add unknown host keys to known_hosts.

- Expected type: BOOLEAN

- Default value: true

- Valid for: [SCHEDULER, FILESYSTEM]


__`xenon.adaptors.ssh.strictHostKeyChecking`__

Enable strict host key checking.

- Expected type: BOOLEAN

- Default value: true

- Valid for: [SCHEDULER, FILESYSTEM]


__`xenon.adaptors.ssh.loadKnownHosts`__

Load the standard known_hosts file.

- Expected type: BOOLEAN

- Default value: true

- Valid for: [XENON]


__`xenon.adaptors.ssh.queue.pollingDelay`__

The polling delay for monitoring running jobs (in milliseconds).

- Expected type: LONG

- Default value: 1000

- Valid for: [SCHEDULER]


__`xenon.adaptors.ssh.queue.multi.maxConcurrentJobs`__

The maximum number of concurrent jobs in the multiq..

- Expected type: INTEGER

- Default value: 4

- Valid for: [SCHEDULER]


__`xenon.adaptors.ssh.gateway`__

The gateway machine used to create an SSH tunnel to the target.

- Expected type: STRING

- Default value: null

- Valid for: [SCHEDULER, FILESYSTEM]



Adaptor: gridengine
--------

The SGE Adaptor submits jobs to a (Sun/Ocacle/Univa) Grid Engine scheduler. This adaptor uses either the local or the ssh adaptor to gain access to the scheduler machine.

#### Supported schemes: ####
ge, sge

#### Supported locations: ####
(locations supported by local), (locations supported by ssh)

#### Supported properties: ####


__`xenon.adaptors.gridengine.ignore.version`__

Skip version check is skipped when connecting to remote machines. WARNING: it is not recommended to use this setting in production environments!

- Expected type: BOOLEAN

- Default value: false

- Valid for: [SCHEDULER]


__`xenon.adaptors.gridengine.accounting.grace.time`__

Number of milliseconds a job is allowed to take going from the queue to the qacct output.

- Expected type: LONG

- Default value: 60000

- Valid for: [SCHEDULER]


__`xenon.adaptors.gridengine.poll.delay`__

Number of milliseconds between polling the status of a job.

- Expected type: LONG

- Default value: 1000

- Valid for: [SCHEDULER]



Adaptor: slurm
--------

The Slurm Adaptor submits jobs to a Slurm scheduler. This adaptor uses either the local or the ssh adaptor to gain access to the scheduler machine.

#### Supported schemes: ####
slurm

#### Supported locations: ####
(locations supported by local), (locations supported by ssh)

#### Supported properties: ####


__`xenon.adaptors.slurm.ignore.version`__

Skip version check is skipped when connecting to remote machines. WARNING: it is not recommended to use this setting in production environments!

- Expected type: BOOLEAN

- Default value: false

- Valid for: [SCHEDULER]


__`xenon.adaptors.slurm.disable.accounting.usage`__

Do not used accounting info of slurm, even when available. Mostly for testing purposes

- Expected type: BOOLEAN

- Default value: false

- Valid for: [SCHEDULER]


__`xenon.adaptors.slurm.poll.delay`__

Number of milliseconds between polling the status of a job.

- Expected type: LONG

- Default value: 1000

- Valid for: [SCHEDULER]



This directory contains documentation for Xenon developers that could not be placed elsewhere.

Torque setup with docker (OS X)
===============================

Dependencies:

```
brew install boot2docker docker
```

Setup:

```
boot2docker up
docker pull agaveapi/torque
docker run -d -h test.torque.xenon.esciencecenter.nl -p 10022:22 \
    --privileged --name test.torque.xenon agaveapi/torque
VBoxManage controlvm \
    "boot2docker-vm" natpf1 "tcp-port10022,tcp,,10022,,10022";
echo -n "\nHost torque\nUser testuser\nHostName localhost\nPort 10022\n" >> \
    ~/.ssh/config
```

Log in:

```
ssh torque # password testuser
```##GridEngine adaptor for Xenon

###XML Parsing

This adaptors parses qstat info from xml. See this website for the technique used:

http://www.gridengine.info/2005/11/17/building-xml-bindings-for-qstat/

####Schemas

The schemas used for generating the XML binding are kept in docs/schemas

# Legal

The Xenon library is copyrighted by the Netherlands eScience Center and released
under the Apache License, Version 2.0. A copy of the license may be obtained from [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

Xenon uses several third-party libraries that have their own (permissive, open 
source) licenses. This document uses the short names from https://spdx.org/licenses/ 
to unambiguously identify licenses. There are 2 sets of licenses to consider.
These are:

1. runtime libraries
1. development libraries, for example for testing and checking code style



## runtime libraries

For now, I limited the runtime libraries to the list of jars from 
./build/install/xenon/lib/ that you get after issuing this command:

```bash
./gradlew installDist
```

- items labeled _N_ have a copy of notice file in legal/
- items labeled _L_ have a copy of license file in legal/

[*Apache-2.0*](https://spdx.org/licenses/Apache-2.0.html#licenseText)

1. aws-s3-2.0.2.jar _N, L_
1. cglib-2.2.1-v20090111.jar _N, L_
1. commons-codec-1.9.jar _N, L_
1. commons-logging-1.2.jar _N, L_
1. commons-net-3.3.jar _N, L_
1. gson-2.5.jar _L_
1. guava-16.0.1.jar _L_
1. guice-3.0.jar _N, L_
1. guice-assistedinject-3.0.jar _N, L_
1. httpclient-4.5.1.jar _N, L_
1. httpcore-4.4.3.jar _N, L_
1. javax.inject-1.jar _L_
1. java-xmlbuilder-1.1.jar _L_
1. jclouds-blobstore-2.0.2.jar _N, L_
1. jclouds-core-2.0.2.jar _N, L_
1. joda-time-2.8.1.jar _N, L_
1. s3-2.0.2.jar _N, L_
1. sardine-5.7.jar _L_
1. sshd-core-1.4.0.jar _N, L_
1. sts-2.0.2.jar _N, L_

[*BSD-3-Clause*](https://spdx.org/licenses/BSD-3-Clause.html#licenseText)

1. asm-3.1.jar _L_

[*CDDL-1.0*](https://spdx.org/licenses/CDDL-1.0.html#licenseText)

1. jsr250-api-1.0.jar _L_
1. jsr311-api-1.1.1.jar _L_

[*LGPL-2.1*](https://spdx.org/licenses/LGPL-2.1.html#licenseText) or [*EPL-1.0*](https://spdx.org/licenses/EPL-1.0.html#licenseText) dual-license

1. logback-classic-1.0.11.jar _L_
1. logback-core-1.0.11.jar _L_

[*MIT*](https://spdx.org/licenses/MIT.html#licenseText)

1. slf4j-api-1.7.22.jar _L_

*public domain*

1. aopalliance-1.0.jar
1. base64-2.3.8.jar

## development libraries

[*Apache-2.0*](https://spdx.org/licenses/Apache-2.0.html#licenseText)

1. gradle-wrapper.jar _L_

[*BSD-3-Clause*](https://spdx.org/licenses/BSD-3-Clause.html#licenseText)

1. hamcrest-core-1.3.jar _L_
1. hamcrest-library-1.3.jar _L_

[*CPL-1.0*](https://spdx.org/licenses/CPL-1.0.html#licenseText)

1. system-rules-1.16.0.jar _L_

[*EPL-1.0*](https://spdx.org/licenses/EPL-1.0.html#licenseText)

1. junit-4.12.jar _L_

[*MIT*](https://spdx.org/licenses/MIT.html#licenseText)

1. mockito-all-1.9.5.jar _L_


