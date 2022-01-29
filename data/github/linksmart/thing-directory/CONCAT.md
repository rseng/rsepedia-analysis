# LinkSmart Thing Directory
[![Docker Pulls](https://img.shields.io/docker/pulls/linksmart/td.svg)](https://hub.docker.com/r/linksmart/td/tags)
[![GitHub tag (latest pre-release)](https://img.shields.io/github/tag-pre/linksmart/thing-directory.svg?label=pre-release)](https://github.com/linksmart/thing-directory/tags)
[![CICD](https://github.com/linksmart/thing-directory/workflows/CICD/badge.svg)](https://github.com/linksmart/thing-directory/actions?query=workflow:CICD)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03075/status.svg)](https://doi.org/10.21105/joss.03075)
  
This is an implementation of the [W3C WoT Thing Description Directory (TDD)](https://w3c.github.io/wot-discovery/), a registry of [Thing Descriptions](https://www.w3.org/TR/wot-thing-description/).

## Getting Started
Visit the following pages to get started:
* [Deployment](https://github.com/linksmart/thing-directory/wiki/Deployment): How to deploy the software, as Docker container, Debian package, or platform-specific binary distributions
* [Configuration](https://github.com/linksmart/thing-directory/wiki/Configuration): How to configure the server software with JSON files and environment variables
* [API Documentation](https://linksmart.github.io/swagger-ui/dist/?url=https://raw.githubusercontent.com/linksmart/thing-directory/master/apidoc/openapi-spec.yml): How to interact with the networking APIs

**Further documentation are available in the [wiki](https://github.com/linksmart/thing-directory/wiki)**.

## Features
* Service Discovery
  * [DNS-SD registration](https://github.com/linksmart/thing-directory/wiki/Discovery-with-DNS-SD)
  * [LinkSmart Service Catalog](https://github.com/linksmart/service-catalog) registration
* RESTful API
  * [HTTP API](https://linksmart.github.io/swagger-ui/dist/?url=https://raw.githubusercontent.com/linksmart/thing-directory/master/apidoc/openapi-spec.yml)
    * Thing Description (TD) CRUD, catalog, and validation
    * XPath 3.0 and JSONPath [query languages](https://github.com/linksmart/thing-directory/wiki/Query-Language)
    * TD validation with JSON Schema ([default](https://github.com/linksmart/thing-directory/blob/master/wot/wot_td_schema.json))
    * Request [authentication](https://github.com/linksmart/go-sec/wiki/Authentication) and [authorization](https://github.com/linksmart/go-sec/wiki/Authorization)
    * JSON-LD response format
* Persistent Storage
  * LevelDB
* CI/CD ([Github Actions](https://github.com/linksmart/thing-directory/actions?query=workflow:CICD))
  * Automated testing
  * Automated builds and releases ([Docker images](https://hub.docker.com/r/linksmart/td/tags?page=1&ordering=last_updated), [binaries](https://github.com/linksmart/thing-directory/releases))

## Development
The dependencies of this package are managed by [Go Modules](https://github.com/golang/go/wiki/Modules).

Clone this repo:
```bash
git clone https://github.com/linksmart/thing-directory.git
cd thing-directory
```

Compile from source:
```bash
go build
```
This will result in an executable named `thing-directory` (linux/macOS) or `thing-directory.exe` (windows).

Get the CLI arguments help (linux/macOS):
```bash
$ ./thing-directory -help
Usage of ./thing-directory:
  -conf string
        Configuration file path (default "conf/thing-directory.json")
  -version
        Print the API version
```

Run (linux/macOS):
```bash
$ ./thing-directory --conf=sample_conf/thing-directory.json
```

To build and run together:
```bash
go run . --conf=sample_conf/thing-directory.json
```

Test all packages (add `-v` flag for verbose results):
```bash
go test ./...
```


## Contributing
Contributions are welcome. 

Please fork, make your changes, and submit a pull request. For major changes, please open an issue first and discuss it with the other authors.
---
title: 'Thing Directory: Simple and lightweight registry of IoT device metadata'
tags:
  - internet of things
  - web of things
  - wireless sensor networks
  - discovery
  - catalog
authors:
  - name: Farshid Tavakolizadeh
    affiliation: 1
  - name: Shreekantha Devasya
    affiliation: 1
affiliations:
 - name: Fraunhofer Institute for Applied Information Technology, Sankt Augustin, Germany
   index: 1
date: 16 December 2020
bibliography: paper.bib

---

# Statement of Need  
<!-- A clear Statement of Need that illustrates the research purpose of the software. -->

The fast emergence of IoT (Internet of Things) technologies has influenced scientific communities to embrace novel sources of information and their potential use in various domains. While the vast amount of sensory data is beneficial, the lack of uniform access interfaces hinders researchers from efficient exploitation. A structured, yet flexible registry is needed to model device metadata and allow exploration and interaction with such devices. In the IoT context, the Things are physical devices (e.g., sensors, actuators, gateways) or virtual ones (e.g., digital twins, proxies, aggregated data sources). An ideal metadata registry for such Things should have a flat learning curve and be easily deployable. This would allow researchers to focus less on interoperability and more on fast prototyping and data application. Registries that are based on established standards are preferred since they can incorporate metadata with many existing Things out-of-the-box. Moreover, the registry software should be lightweight, easily deployable, and free of rigid requirements to support fast prototyping across a wide range of use cases from the edge to the cloud. 

# Thing Directory 
<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

Thing Directory is a searchable registry of metadata for Things. The API is based on W3C 
Web of Things (WoT) Discovery [@WoTDiscovery], a specification for secure discovery of Things. The registry uses JSON-LD (JSON for Linked Data) encoding by default. The JSON format is human-readable and portable; the linked data support makes the data machine-interpretable. The architecture of the registry is modular with a pluggable storage backend, allowing connection to various database systems using drivers. The current implementation, backed by an embedded LevelDB storage, can be deployed on highly constrained single-board computers such as the Raspberry Pi Zero series. It has very minimal idle processing and memory footprints and can scale on demand to utilize all locally available resources. More powerful storage backends can be added to create a horizontally scalable directory in cloud environments. 

The application is packaged as binary distributions, Debian packages, as well as Docker images for easy deployment on a variety of platforms. The data model of the metadata is based on W3C WoT Thing Description (TD) [@WoTTD] which is extensible by design. Thing Directory validates all inputs using a JSON-Schema, describing the data model. This makes it possible to extend the server’s structured data model and validation mechanism with no programming. The various modules of the system are provided as re-usable Go libraries which can be imported by other applications to build functionalities around Thing Descriptions.

<!-- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. -->


## Use case: Assessment of Building Energy Efficiency 
Construction companies often deal with the challenge of delivering target energy-efficient buildings, given specific budget and time constraints. Energy efficiency, as one of the key factors for renovation investments, depends on the availability of various data sources to support the renovation design and planning. These include climate data and building material along with residential comfort and energy consumption patterns. 

As part of the pre-renovation activities, the construction planners deploy various sensors to collect relevant data over a period. Such sensors become part of a wireless sensor network (WSN) and expose data endpoints with the help of one or more gateways. Depending on the protocols, the endpoints require different interaction flows to securely deliver current and historical measurements. The renovation applications need to discover the sensors, their endpoints, and how to interact with them based on search criteria such as the physical location, mapping to the building model, or measurement type. 

The Thing Directory supports engineers in the assessment of building energy efficiency by providing the means to collect and discover the metadata of deployed sensors in an easy and standardized way. Instances of Thing Directory have been deployed in four renovation sites (two apartments, two buildings) across Europe as registries of wireless sensors which are locally accessible over Z-Wave or WiFi. The API has been integrated into applications for profiling of resident usage of building systems, building information modeling, and process modeling and automation. Such applications query sensor metadata based on zoning and sensor types. Once discovered, the metadata provides these applications with necessary details on how to authenticate and query data. Since the Thing Directory is based on an open standard, being integrated with it adds interoperability with WSNs beyond the scope of this use case. The applications will be able to interact with the growing number of compliant WoT devices.  

# Related Work 
<!-- A list of key references, including to other software addressing related needs. -->

There are multiple attempts to modeling and discovering the connected Things and their interfaces. OGC SensorThings API [@sensorthingsSensing;@sensorthingsTasking] has been a successful model for the representation of Things. Sensorthings API consists of two parts: sensing and tasking. The popular implementations of the OGC SensorThings are FROST [@frost] and GOST [@gost] which support mainly the sensing part. FROST has preliminary tasking support. These solutions focus on centralized deployments and are not suitable for a federated scenario and gateways with limited computational power. On top of that, the specification which they are based upon, enforces both metadata and observation modeling. That is not practical in IoT scenarios with heterogenous data formats and interfaces. The survey articles by @DiscoveryCategorization2016 and @DiscoveryReview2020 discuss and categorize several other technologies related to discovery in the IoT field and evaluate them.

The WoT TD [@WoTTD] covers a wide range of Things by providing the semantics to describe the textual metadata, interaction affordances, data schemas, and relations. Thingweb Directory [@Thingweb] has implemented the discovery of WoT TDs using a proprietary API which does not comply with W3C WoT Discovery [@WoTDiscovery] standard. The Thing Directory complies with both W3C WoT Discovery and W3C WoT TD.

# Acknowledgement 
<!-- Acknowledgement of any financial support. -->

This work was conducted as part of the BIMERR project, a European Commission’s Horizon 2020 research and innovation program under grant agreement No 820621. The resulting software is released as part of LinkSmart, an open-source IoT prototyping platform.  

# References
# Building `sys/unix`

The sys/unix package provides access to the raw system call interface of the
underlying operating system. See: https://godoc.org/golang.org/x/sys/unix

Porting Go to a new architecture/OS combination or adding syscalls, types, or
constants to an existing architecture/OS pair requires some manual effort;
however, there are tools that automate much of the process.

## Build Systems

There are currently two ways we generate the necessary files. We are currently
migrating the build system to use containers so the builds are reproducible.
This is being done on an OS-by-OS basis. Please update this documentation as
components of the build system change.

### Old Build System (currently for `GOOS != "linux"`)

The old build system generates the Go files based on the C header files
present on your system. This means that files
for a given GOOS/GOARCH pair must be generated on a system with that OS and
architecture. This also means that the generated code can differ from system
to system, based on differences in the header files.

To avoid this, if you are using the old build system, only generate the Go
files on an installation with unmodified header files. It is also important to
keep track of which version of the OS the files were generated from (ex.
Darwin 14 vs Darwin 15). This makes it easier to track the progress of changes
and have each OS upgrade correspond to a single change.

To build the files for your current OS and architecture, make sure GOOS and
GOARCH are set correctly and run `mkall.sh`. This will generate the files for
your specific system. Running `mkall.sh -n` shows the commands that will be run.

Requirements: bash, go

### New Build System (currently for `GOOS == "linux"`)

The new build system uses a Docker container to generate the go files directly
from source checkouts of the kernel and various system libraries. This means
that on any platform that supports Docker, all the files using the new build
system can be generated at once, and generated files will not change based on
what the person running the scripts has installed on their computer.

The OS specific files for the new build system are located in the `${GOOS}`
directory, and the build is coordinated by the `${GOOS}/mkall.go` program. When
the kernel or system library updates, modify the Dockerfile at
`${GOOS}/Dockerfile` to checkout the new release of the source.

To build all the files under the new build system, you must be on an amd64/Linux
system and have your GOOS and GOARCH set accordingly. Running `mkall.sh` will
then generate all of the files for all of the GOOS/GOARCH pairs in the new build
system. Running `mkall.sh -n` shows the commands that will be run.

Requirements: bash, go, docker

## Component files

This section describes the various files used in the code generation process.
It also contains instructions on how to modify these files to add a new
architecture/OS or to add additional syscalls, types, or constants. Note that
if you are using the new build system, the scripts/programs cannot be called normally.
They must be called from within the docker container.

### asm files

The hand-written assembly file at `asm_${GOOS}_${GOARCH}.s` implements system
call dispatch. There are three entry points:
```
  func Syscall(trap, a1, a2, a3 uintptr) (r1, r2, err uintptr)
  func Syscall6(trap, a1, a2, a3, a4, a5, a6 uintptr) (r1, r2, err uintptr)
  func RawSyscall(trap, a1, a2, a3 uintptr) (r1, r2, err uintptr)
```
The first and second are the standard ones; they differ only in how many
arguments can be passed to the kernel. The third is for low-level use by the
ForkExec wrapper. Unlike the first two, it does not call into the scheduler to
let it know that a system call is running.

When porting Go to an new architecture/OS, this file must be implemented for
each GOOS/GOARCH pair.

### mksysnum

Mksysnum is a Go program located at `${GOOS}/mksysnum.go` (or `mksysnum_${GOOS}.go`
for the old system). This program takes in a list of header files containing the
syscall number declarations and parses them to produce the corresponding list of
Go numeric constants. See `zsysnum_${GOOS}_${GOARCH}.go` for the generated
constants.

Adding new syscall numbers is mostly done by running the build on a sufficiently
new installation of the target OS (or updating the source checkouts for the
new build system). However, depending on the OS, you may need to update the
parsing in mksysnum.

### mksyscall.go

The `syscall.go`, `syscall_${GOOS}.go`, `syscall_${GOOS}_${GOARCH}.go` are
hand-written Go files which implement system calls (for unix, the specific OS,
or the specific OS/Architecture pair respectively) that need special handling
and list `//sys` comments giving prototypes for ones that can be generated.

The mksyscall.go program takes the `//sys` and `//sysnb` comments and converts
them into syscalls. This requires the name of the prototype in the comment to
match a syscall number in the `zsysnum_${GOOS}_${GOARCH}.go` file. The function
prototype can be exported (capitalized) or not.

Adding a new syscall often just requires adding a new `//sys` function prototype
with the desired arguments and a capitalized name so it is exported. However, if
you want the interface to the syscall to be different, often one will make an
unexported `//sys` prototype, an then write a custom wrapper in
`syscall_${GOOS}.go`.

### types files

For each OS, there is a hand-written Go file at `${GOOS}/types.go` (or
`types_${GOOS}.go` on the old system). This file includes standard C headers and
creates Go type aliases to the corresponding C types. The file is then fed
through godef to get the Go compatible definitions. Finally, the generated code
is fed though mkpost.go to format the code correctly and remove any hidden or
private identifiers. This cleaned-up code is written to
`ztypes_${GOOS}_${GOARCH}.go`.

The hardest part about preparing this file is figuring out which headers to
include and which symbols need to be `#define`d to get the actual data
structures that pass through to the kernel system calls. Some C libraries
preset alternate versions for binary compatibility and translate them on the
way in and out of system calls, but there is almost always a `#define` that can
get the real ones.
See `types_darwin.go` and `linux/types.go` for examples.

To add a new type, add in the necessary include statement at the top of the
file (if it is not already there) and add in a type alias line. Note that if
your type is significantly different on different architectures, you may need
some `#if/#elif` macros in your include statements.

### mkerrors.sh

This script is used to generate the system's various constants. This doesn't
just include the error numbers and error strings, but also the signal numbers
an a wide variety of miscellaneous constants. The constants come from the list
of include files in the `includes_${uname}` variable. A regex then picks out
the desired `#define` statements, and generates the corresponding Go constants.
The error numbers and strings are generated from `#include <errno.h>`, and the
signal numbers and strings are generated from `#include <signal.h>`. All of
these constants are written to `zerrors_${GOOS}_${GOARCH}.go` via a C program,
`_errors.c`, which prints out all the constants.

To add a constant, add the header that includes it to the appropriate variable.
Then, edit the regex (if necessary) to match the desired constant. Avoid making
the regex too broad to avoid matching unintended constants.

### mkmerge.go

This program is used to extract duplicate const, func, and type declarations
from the generated architecture-specific files listed below, and merge these
into a common file for each OS.

The merge is performed in the following steps:
1. Construct the set of common code that is idential in all architecture-specific files.
2. Write this common code to the merged file.
3. Remove the common code from all architecture-specific files.


## Generated files

### `zerrors_${GOOS}_${GOARCH}.go`

A file containing all of the system's generated error numbers, error strings,
signal numbers, and constants. Generated by `mkerrors.sh` (see above).

### `zsyscall_${GOOS}_${GOARCH}.go`

A file containing all the generated syscalls for a specific GOOS and GOARCH.
Generated by `mksyscall.go` (see above).

### `zsysnum_${GOOS}_${GOARCH}.go`

A list of numeric constants for all the syscall number of the specific GOOS
and GOARCH. Generated by mksysnum (see above).

### `ztypes_${GOOS}_${GOARCH}.go`

A file containing Go types for passing into (or returning from) syscalls.
Generated by godefs and the types file (see above).
# Exponential Backoff [![GoDoc][godoc image]][godoc] [![Build Status][travis image]][travis] [![Coverage Status][coveralls image]][coveralls]

This is a Go port of the exponential backoff algorithm from [Google's HTTP Client Library for Java][google-http-java-client].

[Exponential backoff][exponential backoff wiki]
is an algorithm that uses feedback to multiplicatively decrease the rate of some process,
in order to gradually find an acceptable rate.
The retries exponentially increase and stop increasing when a certain threshold is met.

## Usage

See https://godoc.org/github.com/cenkalti/backoff#pkg-examples

## Contributing

* I would like to keep this library as small as possible.
* Please don't send a PR without opening an issue and discussing it first.
* If proposed change is not a common use case, I will probably not accept it.

[godoc]: https://godoc.org/github.com/cenkalti/backoff
[godoc image]: https://godoc.org/github.com/cenkalti/backoff?status.png
[travis]: https://travis-ci.org/cenkalti/backoff
[travis image]: https://travis-ci.org/cenkalti/backoff.png?branch=master
[coveralls]: https://coveralls.io/github/cenkalti/backoff?branch=master
[coveralls image]: https://coveralls.io/repos/github/cenkalti/backoff/badge.svg?branch=master

[google-http-java-client]: https://github.com/google/google-http-java-client/blob/da1aa993e90285ec18579f1553339b00e19b3ab5/google-http-client/src/main/java/com/google/api/client/util/ExponentialBackOff.java
[exponential backoff wiki]: http://en.wikipedia.org/wiki/Exponential_backoff

[advanced example]: https://godoc.org/github.com/cenkalti/backoff#example_
ZeroConf: Service Discovery with mDNS
=====================================
ZeroConf is a pure Golang library that employs Multicast DNS-SD for

* browsing and resolving services in your network
* registering own services

in the local network.

It basically implements aspects of the standards
[RFC 6762](https://tools.ietf.org/html/rfc6762) (mDNS) and
[RFC 6763](https://tools.ietf.org/html/rfc6763) (DNS-SD).
Though it does not support all requirements yet, the aim is to provide a compliant solution in the long-term with the community.

By now, it should be compatible to [Avahi](http://avahi.org/) (tested) and Apple's Bonjour (untested).
Target environments: private LAN/Wifi, small or isolated networks.

[![GoDoc](https://godoc.org/github.com/grandcat/zeroconf?status.svg)](https://godoc.org/github.com/grandcat/zeroconf)
[![Go Report Card](https://goreportcard.com/badge/github.com/grandcat/zeroconf)](https://goreportcard.com/report/github.com/grandcat/zeroconf)
[![Build Status](https://travis-ci.com/grandcat/zeroconf.svg?branch=master)](https://travis-ci.com/grandcat/zeroconf)

## Install
Nothing is as easy as that:
```bash
$ go get -u github.com/grandcat/zeroconf
```
This package requires **Go 1.7** (context in std lib) or later.

## Browse for services in your local network

```go
// Discover all services on the network (e.g. _workstation._tcp)
resolver, err := zeroconf.NewResolver(nil)
if err != nil {
    log.Fatalln("Failed to initialize resolver:", err.Error())
}

entries := make(chan *zeroconf.ServiceEntry)
go func(results <-chan *zeroconf.ServiceEntry) {
    for entry := range results {
        log.Println(entry)
    }
    log.Println("No more entries.")
}(entries)

ctx, cancel := context.WithTimeout(context.Background(), time.Second*15)
defer cancel()
err = resolver.Browse(ctx, "_workstation._tcp", "local.", entries)
if err != nil {
    log.Fatalln("Failed to browse:", err.Error())
}

<-ctx.Done()
```
A subtype may added to service name to narrow the set of results. E.g. to browse `_workstation._tcp` with subtype `_windows`, use`_workstation._tcp,_windows`.

See https://github.com/grandcat/zeroconf/blob/master/examples/resolv/client.go.

## Lookup a specific service instance

```go
// Example filled soon.
```

## Register a service

```go
server, err := zeroconf.Register("GoZeroconf", "_workstation._tcp", "local.", 42424, []string{"txtv=0", "lo=1", "la=2"}, nil)
if err != nil {
    panic(err)
}
defer server.Shutdown()

// Clean exit.
sig := make(chan os.Signal, 1)
signal.Notify(sig, os.Interrupt, syscall.SIGTERM)
select {
case <-sig:
    // Exit by user
case <-time.After(time.Second * 120):
    // Exit by timeout
}

log.Println("Shutting down.")
```
Multiple subtypes may be added to service name, separated by commas. E.g `_workstation._tcp,_windows` has subtype `_windows`.

See https://github.com/grandcat/zeroconf/blob/master/examples/register/server.go.

## Features and ToDo's
This list gives a quick impression about the state of this library.
See what needs to be done and submit a pull request :)

* [x] Browse / Lookup / Register services
* [x] Multiple IPv6 / IPv4 addresses support
* [x] Send multiple probes (exp. back-off) if no service answers (*)
* [ ] Timestamp entries for TTL checks
* [ ] Compare new multicasts with already received services

_Notes:_

(*) The denoted features might not be perfectly standards compliant, but shouldn't cause any problems.
    Some tests showed improvements in overall robustness and performance with the features enabled.

## Credits
Great thanks to [hashicorp](https://github.com/hashicorp/mdns) and to [oleksandr](https://github.com/oleksandr/bonjour) and all contributing authors for the code this projects bases upon.
Large parts of the code are still the same.

However, there are several reasons why I decided to create a fork of the original project:
The previous project seems to be unmaintained. There are several useful pull requests waiting. I merged most of them in this project.
Still, the implementation has some bugs and lacks some other features that make it quite unreliable in real LAN environments when running continously.
Last but not least, the aim for this project is to build a solution that targets standard conformance in the long term with the support of the community.
Though, resiliency should remain a top goal.
# UUID package for Go language

[![Build Status](https://travis-ci.org/satori/go.uuid.png?branch=master)](https://travis-ci.org/satori/go.uuid)
[![Coverage Status](https://coveralls.io/repos/github/satori/go.uuid/badge.svg?branch=master)](https://coveralls.io/github/satori/go.uuid)
[![GoDoc](http://godoc.org/github.com/satori/go.uuid?status.png)](http://godoc.org/github.com/satori/go.uuid)

This package provides pure Go implementation of Universally Unique Identifier (UUID). Supported both creation and parsing of UUIDs.

With 100% test coverage and benchmarks out of box.

Supported versions:
* Version 1, based on timestamp and MAC address (RFC 4122)
* Version 2, based on timestamp, MAC address and POSIX UID/GID (DCE 1.1)
* Version 3, based on MD5 hashing (RFC 4122)
* Version 4, based on random numbers (RFC 4122)
* Version 5, based on SHA-1 hashing (RFC 4122)

## Installation

Use the `go` command:

	$ go get github.com/satori/go.uuid

## Requirements

UUID package requires Go >= 1.2.

## Example

```go
package main

import (
	"fmt"
	"github.com/satori/go.uuid"
)

func main() {
	// Creating UUID Version 4
	u1 := uuid.NewV4()
	fmt.Printf("UUIDv4: %s\n", u1)

	// Parsing UUID from string input
	u2, err := uuid.FromString("6ba7b810-9dad-11d1-80b4-00c04fd430c8")
	if err != nil {
		fmt.Printf("Something gone wrong: %s", err)
	}
	fmt.Printf("Successfully parsed: %s", u2)
}
```

## Documentation

[Documentation](http://godoc.org/github.com/satori/go.uuid) is hosted at GoDoc project.

## Links
* [RFC 4122](http://tools.ietf.org/html/rfc4122)
* [DCE 1.1: Authentication and Security Services](http://pubs.opengroup.org/onlinepubs/9696989899/chap5.htm#tagcjh_08_02_01_01)

## Copyright

Copyright (C) 2013-2018 by Maxim Bublis <b@codemonkey.ru>.

UUID package released under MIT License.
See [LICENSE](https://github.com/satori/go.uuid/blob/master/LICENSE) for details.
# gojsonreference
An implementation of JSON Reference - Go language

## Dependencies
https://github.com/xeipuuv/gojsonpointer

## References
http://tools.ietf.org/html/draft-ietf-appsawg-json-pointer-07

http://tools.ietf.org/html/draft-pbryan-zyp-json-ref-03
[![GoDoc](https://godoc.org/github.com/xeipuuv/gojsonschema?status.svg)](https://godoc.org/github.com/xeipuuv/gojsonschema)
[![Build Status](https://travis-ci.org/xeipuuv/gojsonschema.svg)](https://travis-ci.org/xeipuuv/gojsonschema)
[![Go Report Card](https://goreportcard.com/badge/github.com/xeipuuv/gojsonschema)](https://goreportcard.com/report/github.com/xeipuuv/gojsonschema)

# gojsonschema

## Description

An implementation of JSON Schema for the Go  programming language. Supports draft-04, draft-06 and draft-07.

References :

* http://json-schema.org
* http://json-schema.org/latest/json-schema-core.html
* http://json-schema.org/latest/json-schema-validation.html

## Installation

```
go get github.com/xeipuuv/gojsonschema
```

Dependencies :
* [github.com/xeipuuv/gojsonpointer](https://github.com/xeipuuv/gojsonpointer)
* [github.com/xeipuuv/gojsonreference](https://github.com/xeipuuv/gojsonreference)
* [github.com/stretchr/testify/assert](https://github.com/stretchr/testify#assert-package)

## Usage

### Example

```go

package main

import (
    "fmt"
    "github.com/xeipuuv/gojsonschema"
)

func main() {

    schemaLoader := gojsonschema.NewReferenceLoader("file:///home/me/schema.json")
    documentLoader := gojsonschema.NewReferenceLoader("file:///home/me/document.json")

    result, err := gojsonschema.Validate(schemaLoader, documentLoader)
    if err != nil {
        panic(err.Error())
    }

    if result.Valid() {
        fmt.Printf("The document is valid\n")
    } else {
        fmt.Printf("The document is not valid. see errors :\n")
        for _, desc := range result.Errors() {
            fmt.Printf("- %s\n", desc)
        }
    }
}


```

#### Loaders

There are various ways to load your JSON data.
In order to load your schemas and documents,
first declare an appropriate loader :

* Web / HTTP, using a reference :

```go
loader := gojsonschema.NewReferenceLoader("http://www.some_host.com/schema.json")
```

* Local file, using a reference :

```go
loader := gojsonschema.NewReferenceLoader("file:///home/me/schema.json")
```

References use the URI scheme, the prefix (file://) and a full path to the file are required.

* JSON strings :

```go
loader := gojsonschema.NewStringLoader(`{"type": "string"}`)
```

* Custom Go types :

```go
m := map[string]interface{}{"type": "string"}
loader := gojsonschema.NewGoLoader(m)
```

And

```go
type Root struct {
	Users []User `json:"users"`
}

type User struct {
	Name string `json:"name"`
}

...

data := Root{}
data.Users = append(data.Users, User{"John"})
data.Users = append(data.Users, User{"Sophia"})
data.Users = append(data.Users, User{"Bill"})

loader := gojsonschema.NewGoLoader(data)
```

#### Validation

Once the loaders are set, validation is easy :

```go
result, err := gojsonschema.Validate(schemaLoader, documentLoader)
```

Alternatively, you might want to load a schema only once and process to multiple validations :

```go
schema, err := gojsonschema.NewSchema(schemaLoader)
...
result1, err := schema.Validate(documentLoader1)
...
result2, err := schema.Validate(documentLoader2)
...
// etc ...
```

To check the result :

```go
    if result.Valid() {
    	fmt.Printf("The document is valid\n")
    } else {
        fmt.Printf("The document is not valid. see errors :\n")
        for _, err := range result.Errors() {
        	// Err implements the ResultError interface
            fmt.Printf("- %s\n", err)
        }
    }
```


## Loading local schemas

By default `file` and `http(s)` references to external schemas are loaded automatically via the file system or via http(s). An external schema can also be loaded using a `SchemaLoader`.

```go
	sl := gojsonschema.NewSchemaLoader()
	loader1 := gojsonschema.NewStringLoader(`{ "type" : "string" }`)
	err := sl.AddSchema("http://some_host.com/string.json", loader1)
```

Alternatively if your schema already has an `$id` you can use the `AddSchemas` function
```go
	loader2 := gojsonschema.NewStringLoader(`{
			"$id" : "http://some_host.com/maxlength.json",
			"maxLength" : 5
		}`)
	err = sl.AddSchemas(loader2)
```

The main schema should be passed to the `Compile` function. This main schema can then directly reference the added schemas without needing to download them.
```go
	loader3 := gojsonschema.NewStringLoader(`{
		"$id" : "http://some_host.com/main.json",
		"allOf" : [
			{ "$ref" : "http://some_host.com/string.json" },
			{ "$ref" : "http://some_host.com/maxlength.json" }
		]
	}`)

	schema, err := sl.Compile(loader3)

	documentLoader := gojsonschema.NewStringLoader(`"hello world"`)

	result, err := schema.Validate(documentLoader)
```

It's also possible to pass a `ReferenceLoader` to the `Compile` function that references a loaded schema.

```go
err = sl.AddSchemas(loader3)
schema, err := sl.Compile(gojsonschema.NewReferenceLoader("http://some_host.com/main.json"))
``` 

Schemas added by `AddSchema` and `AddSchemas` are only validated when the entire schema is compiled, unless meta-schema validation is used.

## Using a specific draft
By default `gojsonschema` will try to detect the draft of a schema by using the `$schema` keyword and parse it in a strict draft-04, draft-06 or draft-07 mode. If `$schema` is missing, or the draft version is not explicitely set, a hybrid mode is used which merges together functionality of all drafts into one mode.

Autodectection can be turned off with the `AutoDetect` property. Specific draft versions can be specified with the `Draft` property.

```go
sl := gojsonschema.NewSchemaLoader()
sl.Draft = gojsonschema.Draft7
sl.AutoDetect = false
```

If autodetection is on (default), a draft-07 schema can savely reference draft-04 schemas and vice-versa, as long as `$schema` is specified in all schemas.

## Meta-schema validation
Schemas that are added using the `AddSchema`, `AddSchemas` and `Compile` can be validated against their meta-schema by setting the `Validate` property.

The following example will produce an error as `multipleOf` must be a number. If `Validate` is off (default), this error is only returned at the `Compile` step. 

```go
sl := gojsonschema.NewSchemaLoader()
sl.Validate = true
err := sl.AddSchemas(gojsonschema.NewStringLoader(`{
     $id" : "http://some_host.com/invalid.json",
    "$schema": "http://json-schema.org/draft-07/schema#",
    "multipleOf" : true
}`))
 ```
``` 
 ```

Errors returned by meta-schema validation are more readable and contain more information, which helps significantly if you are developing a schema.

Meta-schema validation also works with a custom `$schema`. In case `$schema` is missing, or `AutoDetect` is set to `false`, the meta-schema of the used draft is used.


## Working with Errors

The library handles string error codes which you can customize by creating your own gojsonschema.locale and setting it
```go
gojsonschema.Locale = YourCustomLocale{}
```

However, each error contains additional contextual information. 

Newer versions of `gojsonschema` may have new additional errors, so code that uses a custom locale will need to be updated when this happens.

**err.Type()**: *string* Returns the "type" of error that occurred. Note you can also type check. See below

Note: An error of RequiredType has an err.Type() return value of "required"

    "required": RequiredError
    "invalid_type": InvalidTypeError
    "number_any_of": NumberAnyOfError
    "number_one_of": NumberOneOfError
    "number_all_of": NumberAllOfError
    "number_not": NumberNotError
    "missing_dependency": MissingDependencyError
    "internal": InternalError
    "const": ConstEror
    "enum": EnumError
    "array_no_additional_items": ArrayNoAdditionalItemsError
    "array_min_items": ArrayMinItemsError
    "array_max_items": ArrayMaxItemsError
    "unique": ItemsMustBeUniqueError
    "contains" : ArrayContainsError
    "array_min_properties": ArrayMinPropertiesError
    "array_max_properties": ArrayMaxPropertiesError
    "additional_property_not_allowed": AdditionalPropertyNotAllowedError
    "invalid_property_pattern": InvalidPropertyPatternError
    "invalid_property_name":  InvalidPropertyNameError
    "string_gte": StringLengthGTEError
    "string_lte": StringLengthLTEError
    "pattern": DoesNotMatchPatternError
    "multiple_of": MultipleOfError
    "number_gte": NumberGTEError
    "number_gt": NumberGTError
    "number_lte": NumberLTEError
    "number_lt": NumberLTError
    "condition_then" : ConditionThenError
    "condition_else" : ConditionElseError

**err.Value()**: *interface{}* Returns the value given

**err.Context()**: *gojsonschema.JsonContext* Returns the context. This has a String() method that will print something like this: (root).firstName

**err.Field()**: *string* Returns the fieldname in the format firstName, or for embedded properties, person.firstName. This returns the same as the String() method on *err.Context()* but removes the (root). prefix.

**err.Description()**: *string* The error description. This is based on the locale you are using. See the beginning of this section for overwriting the locale with a custom implementation.

**err.DescriptionFormat()**: *string* The error description format. This is relevant if you are adding custom validation errors afterwards to the result.

**err.Details()**: *gojsonschema.ErrorDetails* Returns a map[string]interface{} of additional error details specific to the error. For example, GTE errors will have a "min" value, LTE will have a "max" value. See errors.go for a full description of all the error details. Every error always contains a "field" key that holds the value of *err.Field()*

Note in most cases, the err.Details() will be used to generate replacement strings in your locales, and not used directly. These strings follow the text/template format i.e.
```
{{.field}} must be greater than or equal to {{.min}}
```

The library allows you to specify custom template functions, should you require more complex error message handling.
```go
gojsonschema.ErrorTemplateFuncs = map[string]interface{}{
	"allcaps": func(s string) string {
		return strings.ToUpper(s)
	},
}
```

Given the above definition, you can use the custom function `"allcaps"` in your localization templates:
```
{{allcaps .field}} must be greater than or equal to {{.min}}
```

The above error message would then be rendered with the `field` value in capital letters. For example:
```
"PASSWORD must be greater than or equal to 8"
```

Learn more about what types of template functions you can use in `ErrorTemplateFuncs` by referring to Go's [text/template FuncMap](https://golang.org/pkg/text/template/#FuncMap) type.

## Formats
JSON Schema allows for optional "format" property to validate instances against well-known formats. gojsonschema ships with all of the formats defined in the spec that you can use like this:

````json
{"type": "string", "format": "email"}
````

Not all formats defined in draft-07 are available. Implemented formats are:

* `date`
* `time`
* `date-time`
* `hostname`. Subdomains that start with a number are also supported, but this means that it doesn't strictly follow [RFC1034](http://tools.ietf.org/html/rfc1034#section-3.5) and has the implication that ipv4 addresses are also recognized as valid hostnames.
* `email`. Go's email parser deviates slightly from [RFC5322](https://tools.ietf.org/html/rfc5322). Includes unicode support.
* `idn-email`. Same caveat as `email`.
* `ipv4`
* `ipv6`
* `uri`. Includes unicode support.
* `uri-reference`. Includes unicode support.
* `iri`
* `iri-reference`
* `uri-template`
* `uuid`
* `regex`. Go uses the [RE2](https://github.com/google/re2/wiki/Syntax) engine and is not [ECMA262](http://www.ecma-international.org/publications/files/ECMA-ST/Ecma-262.pdf) compatible.
* `json-pointer`
* `relative-json-pointer`

`email`, `uri` and `uri-reference` use the same validation code as their unicode counterparts `idn-email`, `iri` and `iri-reference`. If you rely on unicode support you should use the specific 
unicode enabled formats for the sake of interoperability as other implementations might not support unicode in the regular formats.

The validation code for `uri`, `idn-email` and their relatives use mostly standard library code.

For repetitive or more complex formats, you can create custom format checkers and add them to gojsonschema like this:

```go
// Define the format checker
type RoleFormatChecker struct {}

// Ensure it meets the gojsonschema.FormatChecker interface
func (f RoleFormatChecker) IsFormat(input interface{}) bool {

    asString, ok := input.(string)
    if ok == false {
        return false
    }

    return strings.HasPrefix("ROLE_", asString)
}

// Add it to the library
gojsonschema.FormatCheckers.Add("role", RoleFormatChecker{})
````

Now to use in your json schema:
````json
{"type": "string", "format": "role"}
````

Another example would be to check if the provided integer matches an id on database:

JSON schema:
```json
{"type": "integer", "format": "ValidUserId"}
```

```go
// Define the format checker
type ValidUserIdFormatChecker struct {}

// Ensure it meets the gojsonschema.FormatChecker interface
func (f ValidUserIdFormatChecker) IsFormat(input interface{}) bool {

    asFloat64, ok := input.(float64) // Numbers are always float64 here
    if ok == false {
        return false
    }

    // XXX
    // do the magic on the database looking for the int(asFloat64)

    return true
}

// Add it to the library
gojsonschema.FormatCheckers.Add("ValidUserId", ValidUserIdFormatChecker{})
````

Formats can also be removed, for example if you want to override one of the formats that is defined by default.

```go
gojsonschema.FormatCheckers.Remove("hostname")
```


## Additional custom validation
After the validation has run and you have the results, you may add additional
errors using `Result.AddError`. This is useful to maintain the same format within the resultset instead
of having to add special exceptions for your own errors. Below is an example.

```go
type AnswerInvalidError struct {
    gojsonschema.ResultErrorFields
}

func newAnswerInvalidError(context *gojsonschema.JsonContext, value interface{}, details gojsonschema.ErrorDetails) *AnswerInvalidError {
    err := AnswerInvalidError{}
    err.SetContext(context)
    err.SetType("custom_invalid_error")
    // it is important to use SetDescriptionFormat() as this is used to call SetDescription() after it has been parsed
    // using the description of err will be overridden by this.
    err.SetDescriptionFormat("Answer to the Ultimate Question of Life, the Universe, and Everything is {{.answer}}")
    err.SetValue(value)
    err.SetDetails(details)

    return &err
}

func main() {
    // ...
    schema, err := gojsonschema.NewSchema(schemaLoader)
    result, err := gojsonschema.Validate(schemaLoader, documentLoader)

    if true { // some validation
        jsonContext := gojsonschema.NewJsonContext("question", nil)
        errDetail := gojsonschema.ErrorDetails{
            "answer": 42,
        }
        result.AddError(
            newAnswerInvalidError(
                gojsonschema.NewJsonContext("answer", jsonContext),
                52,
                errDetail,
            ),
            errDetail,
        )
    }

    return result, err

}
```

This is especially useful if you want to add validation beyond what the
json schema drafts can provide such business specific logic.

## Uses

gojsonschema uses the following test suite :

https://github.com/json-schema/JSON-Schema-Test-Suite
# gojsonpointer
An implementation of JSON Pointer - Go language

## Usage
	jsonText := `{
		"name": "Bobby B",
		"occupation": {
			"title" : "King",
			"years" : 15,
			"heir" : "Joffrey B"			
		}
	}`
	
    var jsonDocument map[string]interface{}
    json.Unmarshal([]byte(jsonText), &jsonDocument)
    
    //create a JSON pointer
    pointerString := "/occupation/title"
    pointer, _ := NewJsonPointer(pointerString)
    
    //SET a new value for the "title" in the document     
    pointer.Set(jsonDocument, "Supreme Leader of Westeros")
    
    //GET the new "title" from the document
    title, _, _ := pointer.Get(jsonDocument)
    fmt.Println(title) //outputs "Supreme Leader of Westeros"
    
    //DELETE the "heir" from the document
    deletePointer := NewJsonPointer("/occupation/heir")
    deletePointer.Delete(jsonDocument)
    
    b, _ := json.Marshal(jsonDocument)
    fmt.Println(string(b))
    //outputs `{"name":"Bobby B","occupation":{"title":"Supreme Leader of Westeros","years":15}}`


## References
http://tools.ietf.org/html/draft-ietf-appsawg-json-pointer-07

### Note
The 4.Evaluation part of the previous reference, starting with 'If the currently referenced value is a JSON array, the reference token MUST contain either...' is not implemented.

[![GoDoc](https://godoc.org/github.com/eclipse/paho.mqtt.golang?status.svg)](https://godoc.org/github.com/eclipse/paho.mqtt.golang)
[![Go Report Card](https://goreportcard.com/badge/github.com/eclipse/paho.mqtt.golang)](https://goreportcard.com/report/github.com/eclipse/paho.mqtt.golang)

Eclipse Paho MQTT Go client
===========================


This repository contains the source code for the [Eclipse Paho](http://eclipse.org/paho) MQTT Go client library. 

This code builds a library which enable applications to connect to an [MQTT](http://mqtt.org) broker to publish messages, and to subscribe to topics and receive published messages.

This library supports a fully asynchronous mode of operation.


Installation and Build
----------------------

This client is designed to work with the standard Go tools, so installation is as easy as:

```
go get github.com/eclipse/paho.mqtt.golang
```

The client depends on Google's [websockets](https://godoc.org/golang.org/x/net/websocket) and [proxy](https://godoc.org/golang.org/x/net/proxy) package, 
also easily installed with the commands:

```
go get golang.org/x/net/websocket
go get golang.org/x/net/proxy
```


Usage and API
-------------

Detailed API documentation is available by using to godoc tool, or can be browsed online
using the [godoc.org](http://godoc.org/github.com/eclipse/paho.mqtt.golang) service.

Make use of the library by importing it in your Go client source code. For example,
```
import "github.com/eclipse/paho.mqtt.golang"
```

Samples are available in the `cmd` directory for reference.


Runtime tracing
---------------

Tracing is enabled by assigning logs (from the Go log package) to the logging endpoints, ERROR, CRITICAL, WARN and DEBUG


Reporting bugs
--------------

Please report bugs by raising issues for this project in github https://github.com/eclipse/paho.mqtt.golang/issues 


More information
----------------

Discussion of the Paho clients takes place on the [Eclipse paho-dev mailing list](https://dev.eclipse.org/mailman/listinfo/paho-dev).

General questions about the MQTT protocol are discussed in the [MQTT Google Group](https://groups.google.com/forum/?hl=en-US&fromgroups#!forum/mqtt).

There is much more information available via the [MQTT community site](http://mqtt.org).
Contributing to Paho
====================

Thanks for your interest in this project.

Project description:
--------------------

The Paho project has been created to provide scalable open-source implementations of open and standard messaging protocols aimed at new, existing, and emerging applications for Machine-to-Machine (M2M) and Internet of Things (IoT).
Paho reflects the inherent physical and cost constraints of device connectivity. Its objectives include effective levels of decoupling between devices and applications, designed to keep markets open and encourage the rapid growth of scalable Web and Enterprise middleware and applications. Paho is being kicked off with MQTT publish/subscribe client implementations for use on embedded platforms, along with corresponding server support as determined by the community.

- https://projects.eclipse.org/projects/technology.paho

Developer resources:
--------------------

Information regarding source code management, builds, coding standards, and more.

- https://projects.eclipse.org/projects/technology.paho/developer

Contributor License Agreement:
------------------------------

Before your contribution can be accepted by the project, you need to create and electronically sign the Eclipse Foundation Contributor License Agreement (CLA).

- http://www.eclipse.org/legal/CLA.php

Contributing Code:
------------------

The Go client is developed in Github, see their documentation on the process of forking and pull requests; https://help.github.com/categories/collaborating-on-projects-using-pull-requests/

Git commit messages should follow the style described here;

http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html

Contact:
--------

Contact the project developers via the project's "dev" list.

- https://dev.eclipse.org/mailman/listinfo/paho-dev

Search for bugs:
----------------

This project uses Github issues to track ongoing development and issues.

- https://github.com/eclipse/paho.mqtt.golang/issues

Create a new bug:
-----------------

Be sure to search for existing bugs before you create another one. Remember that contributions are always welcome!

- https://github.com/eclipse/paho.mqtt.golang/issues
[![Build Status](https://travis-ci.org/miekg/dns.svg?branch=master)](https://travis-ci.org/miekg/dns)
[![Code Coverage](https://img.shields.io/codecov/c/github/miekg/dns/master.svg)](https://codecov.io/github/miekg/dns?branch=master)
[![Go Report Card](https://goreportcard.com/badge/github.com/miekg/dns)](https://goreportcard.com/report/miekg/dns)
[![](https://godoc.org/github.com/miekg/dns?status.svg)](https://godoc.org/github.com/miekg/dns)

# Alternative (more granular) approach to a DNS library

> Less is more.

Complete and usable DNS library. All Resource Records are supported, including the DNSSEC types.
It follows a lean and mean philosophy. If there is stuff you should know as a DNS programmer there
isn't a convenience function for it. Server side and client side programming is supported, i.e. you
can build servers and resolvers with it.

We try to keep the "master" branch as sane as possible and at the bleeding edge of standards,
avoiding breaking changes wherever reasonable. We support the last two versions of Go.

# Goals

* KISS;
* Fast;
* Small API. If it's easy to code in Go, don't make a function for it.

# Users

A not-so-up-to-date-list-that-may-be-actually-current:

* https://github.com/coredns/coredns
* https://cloudflare.com
* https://github.com/abh/geodns
* http://www.statdns.com/
* http://www.dnsinspect.com/
* https://github.com/chuangbo/jianbing-dictionary-dns
* http://www.dns-lg.com/
* https://github.com/fcambus/rrda
* https://github.com/kenshinx/godns
* https://github.com/skynetservices/skydns
* https://github.com/hashicorp/consul
* https://github.com/DevelopersPL/godnsagent
* https://github.com/duedil-ltd/discodns
* https://github.com/StalkR/dns-reverse-proxy
* https://github.com/tianon/rawdns
* https://mesosphere.github.io/mesos-dns/
* https://pulse.turbobytes.com/
* https://github.com/fcambus/statzone
* https://github.com/benschw/dns-clb-go
* https://github.com/corny/dnscheck for <http://public-dns.info/>
* https://namesmith.io
* https://github.com/miekg/unbound
* https://github.com/miekg/exdns
* https://dnslookup.org
* https://github.com/looterz/grimd
* https://github.com/phamhongviet/serf-dns
* https://github.com/mehrdadrad/mylg
* https://github.com/bamarni/dockness
* https://github.com/fffaraz/microdns
* http://kelda.io
* https://github.com/ipdcode/hades <https://jd.com>
* https://github.com/StackExchange/dnscontrol/
* https://www.dnsperf.com/
* https://dnssectest.net/
* https://dns.apebits.com
* https://github.com/oif/apex
* https://github.com/jedisct1/dnscrypt-proxy
* https://github.com/jedisct1/rpdns
* https://github.com/xor-gate/sshfp
* https://github.com/rs/dnstrace
* https://blitiri.com.ar/p/dnss ([github mirror](https://github.com/albertito/dnss))
* https://github.com/semihalev/sdns
* https://render.com
* https://github.com/peterzen/goresolver
* https://github.com/folbricht/routedns

Send pull request if you want to be listed here.

# Features

* UDP/TCP queries, IPv4 and IPv6
* RFC 1035 zone file parsing ($INCLUDE, $ORIGIN, $TTL and $GENERATE (for all record types) are supported
* Fast
* Server side programming (mimicking the net/http package)
* Client side programming
* DNSSEC: signing, validating and key generation for DSA, RSA, ECDSA and Ed25519
* EDNS0, NSID, Cookies
* AXFR/IXFR
* TSIG, SIG(0)
* DNS over TLS (DoT): encrypted connection between client and server over TCP
* DNS name compression

Have fun!

Miek Gieben  -  2010-2012  -  <miek@miek.nl>
DNS Authors 2012-

# Building

This library uses Go modules and uses semantic versioning. Building is done with the `go` tool, so
the following should work:

    go get github.com/miekg/dns
    go build github.com/miekg/dns

## Examples

A short "how to use the API" is at the beginning of doc.go (this also will show when you call `godoc
github.com/miekg/dns`).

Example programs can be found in the `github.com/miekg/exdns` repository.

## Supported RFCs

*all of them*

* 103{4,5} - DNS standard
* 1348 - NSAP record (removed the record)
* 1982 - Serial Arithmetic
* 1876 - LOC record
* 1995 - IXFR
* 1996 - DNS notify
* 2136 - DNS Update (dynamic updates)
* 2181 - RRset definition - there is no RRset type though, just []RR
* 2537 - RSAMD5 DNS keys
* 2065 - DNSSEC (updated in later RFCs)
* 2671 - EDNS record
* 2782 - SRV record
* 2845 - TSIG record
* 2915 - NAPTR record
* 2929 - DNS IANA Considerations
* 3110 - RSASHA1 DNS keys
* 3123 - APL record
* 3225 - DO bit (DNSSEC OK)
* 340{1,2,3} - NAPTR record
* 3445 - Limiting the scope of (DNS)KEY
* 3597 - Unknown RRs
* 403{3,4,5} - DNSSEC + validation functions
* 4255 - SSHFP record
* 4343 - Case insensitivity
* 4408 - SPF record
* 4509 - SHA256 Hash in DS
* 4592 - Wildcards in the DNS
* 4635 - HMAC SHA TSIG
* 4701 - DHCID
* 4892 - id.server
* 5001 - NSID
* 5155 - NSEC3 record
* 5205 - HIP record
* 5702 - SHA2 in the DNS
* 5936 - AXFR
* 5966 - TCP implementation recommendations
* 6605 - ECDSA
* 6725 - IANA Registry Update
* 6742 - ILNP DNS
* 6840 - Clarifications and Implementation Notes for DNS Security
* 6844 - CAA record
* 6891 - EDNS0 update
* 6895 - DNS IANA considerations
* 6944 - DNSSEC DNSKEY Algorithm Status
* 6975 - Algorithm Understanding in DNSSEC
* 7043 - EUI48/EUI64 records
* 7314 - DNS (EDNS) EXPIRE Option
* 7477 - CSYNC RR
* 7828 - edns-tcp-keepalive EDNS0 Option
* 7553 - URI record
* 7858 - DNS over TLS: Initiation and Performance Considerations
* 7871 - EDNS0 Client Subnet
* 7873 - Domain Name System (DNS) Cookies
* 8080 - EdDSA for DNSSEC
* 8499 - DNS Terminology

## Loosely Based Upon

* ldns - <https://nlnetlabs.nl/projects/ldns/about/>
* NSD - <https://nlnetlabs.nl/projects/nsd/about/>
* Net::DNS - <http://www.net-dns.org/>
* GRONG - <https://github.com/bortzmeyer/grong>
# elog
[![GoDoc](https://godoc.org/github.com/farshidtz/elog?status.svg)](https://godoc.org/github.com/farshidtz/elog)
[![Build Status](https://travis-ci.org/farshidtz/elog.svg?branch=master)](https://travis-ci.org/farshidtz/elog)

elog extends Go's built-in [log](https://golang.org/pkg/log) package to enable simple levelled logging and to modify the formatting. 

## Debugging mode
The debugging mode can be enabled by [setting the configured environment variable](https://github.com/farshidtz/elog#setting-an-environment-variable) ("DEBUG" by default) to 1. Alternatively, the debugging can be enabled using a [flag](https://github.com/farshidtz/elog#enable-debugging-with-a-flag).

## Usage
Get the package

    go get github.com/farshidtz/elog
Import
```go
import "github.com/farshidtz/elog"
```

Examples
```go
10	// Initialize with the default configuration
11	logger := elog.New("[main] ", nil)
12	
13	logger.Println("Hello world!")
14	logger.Debugln("Hello world!")
15	logger.Fatalln("Hello world!")
```
Debugging not enabled
```
2016/07/14 16:47:04 [main] Hello world!
2016/07/14 16:47:04 [main] Hello world!
exit with status 1
```
Debugging enabled
```
2016/07/14 16:47:04 [main] main.go:13: Hello world!
2016/07/14 16:47:04 [debug] main.go:14 Hello world!
2016/07/14 16:47:04 [main] main.go:15 Hello world!
exit with status 1
```
### Error logging
```go
01	package main
02
03	import "github.com/farshidtz/elog"
04
05	var logger *elog.Logger
06
07	func main() {
08		logger = elog.New("[main] ", nil)
09
10		v, err := divide(10, 0)
11		if err != nil {
12			logger.Fatalln(err)
13		}
14		logger.Println(v)
15	}
16
17	func divide(a, b int) (float64, error) {
18		if b == 0 {
19			return 0, logger.Errorf("Cannot divide by zero")
20			// The error is logged in debugging mode
21		}
22		return float64(a) / float64(b), nil
23	}
```
Debugging not enabled
```
2016/07/20 16:38:31 [main] main.go:12: Cannot divide by zero
```
Debugging enabled
```
2016/07/20 16:38:31 [debug] main.go:19: Cannot divide by zero
2016/07/20 16:38:31 [main] main.go:12: Cannot divide by zero
```
### Configuration
```go
10	logger := elog.New("[I] ", &elog.Config{
11	  TimeFormat: time.RFC3339, 
12	  DebugPrefix: "[D] ", 
13	})
14	
15	logger.Println("Hello world!")
16	logger.Debugln("Hello world!")
17	logger.Fatalln("Hello world!")
```
Debugging not enabled
```
2016-07-14T16:57:15Z [I] Hello world!
2016-07-14T16:57:15Z [I] Hello world!
exit with status 1
```
Debugging enabled
```
2016-07-14T16:57:15Z [I] main.go:15 Hello world!
2016-07-14T16:57:15Z [D] main.go:16 Hello world!
2016-07-14T16:57:15Z [I] main.go:17 Hello world!
exit with status 1
```
### Initialization with init()
Alternatively, a global logger can be instantiated in the init() function and used throughout the package.  
```go
import "github.com/farshidtz/elog"

var logger *elog.Logger

func init() {
	logger = elog.New("[main] ", nil)
}
```

### Enable debugging with a flag
```go
package main

import (
	"github.com/farshidtz/elog"
	"flag"
)

var debugFlag = flag.Bool("d", false, "Enable debugging")
func main() {
	flag.Parse()
	logger := elog.New("[main] ", &elog.Config{
	  DebugEnabled: debugFlag,
	})
	
	logger.Println("Hello World!")
	logger.Debugln("Hello World!")
}
```
Example (Windows PowerShell):
```
PS C:\logging> .\main.exe
2016/07/20 15:55:31 [main] Hello World!

PS C:\logging> .\main.exe -d
2016/07/20 15:55:32 [main] main.go:15: Hello World!
2016/07/20 15:55:32 [debug] main.go:16: Hello World!
```

### Setting an environment variable
```
Unix:
export DEBUG=1

Command Prompt:
set DEBUG=1

PowerShell:
$env:DEBUG=1
```
## Documentation
For Go log's build-in methods: [golang.org/pkg/log](https://golang.org/pkg/log)

For extended elog methods: [godoc.org/github.com/farshidtz/elog](https://godoc.org/github.com/farshidtz/elog)# mqtt-match [![GoDoc](https://godoc.org/github.com/farshidtz/mqtt-match?status.svg)](https://godoc.org/github.com/farshidtz/mqtt-match) [![Build Status](https://travis-ci.org/farshidtz/mqtt-match.svg?branch=master)](https://travis-ci.org/farshidtz/mqtt-match)

Match mqtt formatted topic strings to strings, e.g. `foo/+` should match `foo/bar`.

### Usage

```go
package main

import (
    "fmt"
    mqtttopic "github.com/farshidtz/mqtt-match"
)

func main(){
    fmt.Println(mqtttopic.Match("foo/+", "foo/bar"))
    // true
}
```
### Copyrights Notice
This package is a translation of [mqtt-match](https://github.com/ralphtheninja/mqtt-match) for Go.
[![Build Status](https://travis-ci.org/ancientlore/go-avltree.svg?branch=master)](https://travis-ci.org/ancientlore/go-avltree)
[![Coverage Status](https://coveralls.io/repos/ancientlore/go-avltree/badge.svg)](https://coveralls.io/r/ancientlore/go-avltree)
[![GoDoc](https://godoc.org/github.com/ancientlore/go-avltree?status.png)](https://godoc.org/github.com/ancientlore/go-avltree)
[gocover](http://gocover.io/github.com/ancientlore/go-avltree)

An [AVL tree](http://en.wikipedia.org/wiki/AVL_tree) (Adel'son-Vel'skii & Landis) is a binary search tree in which the heights of the left and right subtrees of the root differ by at most one and in which the left and right subtrees are again AVL trees.

With each node of an AVL tree is associated a balance factor that is Left High, Equal, or Right High according, respectively, as the left subtree has height greater than, equal to, or less than that of the right subtree.

The AVL tree is, in practice, balanced quite well. It can (at the worst case) become skewed to the left or right, but never so much that it becomes inefficient. The balancing is done as items are added or deleted.

This version is enhanced to allow "indexing" of values in the tree; however, the indexes are not stable as the tree could be resorted as items are added or removed.

It is safe to iterate or search a tree from multiple threads provided that no threads are modifying the tree.

The tree works on interface{} types and there is also a specialization for strings, pairs, and objects. Additionally, the tree supports iteration and a channel iterator.

	t.Do(func(z interface{}) { if z.(int) % 3333 == 0 { fmt.Printf("%d ", z); } })

	for v := range t.Iter() {
        	if v.(int) % 3333 == 0 {
                	fmt.Printf("%d ", v);
        	}
	}

To install, you can use:

	go get github.com/ancientlore/go-avltree

See some sample code at https://gist.github.com/ancientlore/8855122
context
=======
[![Build Status](https://travis-ci.org/gorilla/context.png?branch=master)](https://travis-ci.org/gorilla/context)

gorilla/context is a general purpose registry for global request variables.

> Note: gorilla/context, having been born well before `context.Context` existed, does not play well
> with the shallow copying of the request that [`http.Request.WithContext`](https://golang.org/pkg/net/http/#Request.WithContext) (added to net/http Go 1.7 onwards) performs. You should either use *just* gorilla/context, or moving forward, the new `http.Request.Context()`.

Read the full documentation here: http://www.gorillatoolkit.org/pkg/context
# gorilla/mux

[![GoDoc](https://godoc.org/github.com/gorilla/mux?status.svg)](https://godoc.org/github.com/gorilla/mux)
[![Build Status](https://travis-ci.org/gorilla/mux.svg?branch=master)](https://travis-ci.org/gorilla/mux)
[![CircleCI](https://circleci.com/gh/gorilla/mux.svg?style=svg)](https://circleci.com/gh/gorilla/mux)
[![Sourcegraph](https://sourcegraph.com/github.com/gorilla/mux/-/badge.svg)](https://sourcegraph.com/github.com/gorilla/mux?badge)

![Gorilla Logo](http://www.gorillatoolkit.org/static/images/gorilla-icon-64.png)

https://www.gorillatoolkit.org/pkg/mux

Package `gorilla/mux` implements a request router and dispatcher for matching incoming requests to
their respective handler.

The name mux stands for "HTTP request multiplexer". Like the standard `http.ServeMux`, `mux.Router` matches incoming requests against a list of registered routes and calls a handler for the route that matches the URL or other conditions. The main features are:

* It implements the `http.Handler` interface so it is compatible with the standard `http.ServeMux`.
* Requests can be matched based on URL host, path, path prefix, schemes, header and query values, HTTP methods or using custom matchers.
* URL hosts, paths and query values can have variables with an optional regular expression.
* Registered URLs can be built, or "reversed", which helps maintaining references to resources.
* Routes can be used as subrouters: nested routes are only tested if the parent route matches. This is useful to define groups of routes that share common conditions like a host, a path prefix or other repeated attributes. As a bonus, this optimizes request matching.

---

* [Install](#install)
* [Examples](#examples)
* [Matching Routes](#matching-routes)
* [Static Files](#static-files)
* [Registered URLs](#registered-urls)
* [Walking Routes](#walking-routes)
* [Graceful Shutdown](#graceful-shutdown)
* [Middleware](#middleware)
* [Handling CORS Requests](#handling-cors-requests)
* [Testing Handlers](#testing-handlers)
* [Full Example](#full-example)

---

## Install

With a [correctly configured](https://golang.org/doc/install#testing) Go toolchain:

```sh
go get -u github.com/gorilla/mux
```

## Examples

Let's start registering a couple of URL paths and handlers:

```go
func main() {
    r := mux.NewRouter()
    r.HandleFunc("/", HomeHandler)
    r.HandleFunc("/products", ProductsHandler)
    r.HandleFunc("/articles", ArticlesHandler)
    http.Handle("/", r)
}
```

Here we register three routes mapping URL paths to handlers. This is equivalent to how `http.HandleFunc()` works: if an incoming request URL matches one of the paths, the corresponding handler is called passing (`http.ResponseWriter`, `*http.Request`) as parameters.

Paths can have variables. They are defined using the format `{name}` or `{name:pattern}`. If a regular expression pattern is not defined, the matched variable will be anything until the next slash. For example:

```go
r := mux.NewRouter()
r.HandleFunc("/products/{key}", ProductHandler)
r.HandleFunc("/articles/{category}/", ArticlesCategoryHandler)
r.HandleFunc("/articles/{category}/{id:[0-9]+}", ArticleHandler)
```

The names are used to create a map of route variables which can be retrieved calling `mux.Vars()`:

```go
func ArticlesCategoryHandler(w http.ResponseWriter, r *http.Request) {
    vars := mux.Vars(r)
    w.WriteHeader(http.StatusOK)
    fmt.Fprintf(w, "Category: %v\n", vars["category"])
}
```

And this is all you need to know about the basic usage. More advanced options are explained below.

### Matching Routes

Routes can also be restricted to a domain or subdomain. Just define a host pattern to be matched. They can also have variables:

```go
r := mux.NewRouter()
// Only matches if domain is "www.example.com".
r.Host("www.example.com")
// Matches a dynamic subdomain.
r.Host("{subdomain:[a-z]+}.example.com")
```

There are several other matchers that can be added. To match path prefixes:

```go
r.PathPrefix("/products/")
```

...or HTTP methods:

```go
r.Methods("GET", "POST")
```

...or URL schemes:

```go
r.Schemes("https")
```

...or header values:

```go
r.Headers("X-Requested-With", "XMLHttpRequest")
```

...or query values:

```go
r.Queries("key", "value")
```

...or to use a custom matcher function:

```go
r.MatcherFunc(func(r *http.Request, rm *RouteMatch) bool {
    return r.ProtoMajor == 0
})
```

...and finally, it is possible to combine several matchers in a single route:

```go
r.HandleFunc("/products", ProductsHandler).
  Host("www.example.com").
  Methods("GET").
  Schemes("http")
```

Routes are tested in the order they were added to the router. If two routes match, the first one wins:

```go
r := mux.NewRouter()
r.HandleFunc("/specific", specificHandler)
r.PathPrefix("/").Handler(catchAllHandler)
```

Setting the same matching conditions again and again can be boring, so we have a way to group several routes that share the same requirements. We call it "subrouting".

For example, let's say we have several URLs that should only match when the host is `www.example.com`. Create a route for that host and get a "subrouter" from it:

```go
r := mux.NewRouter()
s := r.Host("www.example.com").Subrouter()
```

Then register routes in the subrouter:

```go
s.HandleFunc("/products/", ProductsHandler)
s.HandleFunc("/products/{key}", ProductHandler)
s.HandleFunc("/articles/{category}/{id:[0-9]+}", ArticleHandler)
```

The three URL paths we registered above will only be tested if the domain is `www.example.com`, because the subrouter is tested first. This is not only convenient, but also optimizes request matching. You can create subrouters combining any attribute matchers accepted by a route.

Subrouters can be used to create domain or path "namespaces": you define subrouters in a central place and then parts of the app can register its paths relatively to a given subrouter.

There's one more thing about subroutes. When a subrouter has a path prefix, the inner routes use it as base for their paths:

```go
r := mux.NewRouter()
s := r.PathPrefix("/products").Subrouter()
// "/products/"
s.HandleFunc("/", ProductsHandler)
// "/products/{key}/"
s.HandleFunc("/{key}/", ProductHandler)
// "/products/{key}/details"
s.HandleFunc("/{key}/details", ProductDetailsHandler)
```


### Static Files

Note that the path provided to `PathPrefix()` represents a "wildcard": calling
`PathPrefix("/static/").Handler(...)` means that the handler will be passed any
request that matches "/static/\*". This makes it easy to serve static files with mux:

```go
func main() {
    var dir string

    flag.StringVar(&dir, "dir", ".", "the directory to serve files from. Defaults to the current dir")
    flag.Parse()
    r := mux.NewRouter()

    // This will serve files under http://localhost:8000/static/<filename>
    r.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(dir))))

    srv := &http.Server{
        Handler:      r,
        Addr:         "127.0.0.1:8000",
        // Good practice: enforce timeouts for servers you create!
        WriteTimeout: 15 * time.Second,
        ReadTimeout:  15 * time.Second,
    }

    log.Fatal(srv.ListenAndServe())
}
```

### Registered URLs

Now let's see how to build registered URLs.

Routes can be named. All routes that define a name can have their URLs built, or "reversed". We define a name calling `Name()` on a route. For example:

```go
r := mux.NewRouter()
r.HandleFunc("/articles/{category}/{id:[0-9]+}", ArticleHandler).
  Name("article")
```

To build a URL, get the route and call the `URL()` method, passing a sequence of key/value pairs for the route variables. For the previous route, we would do:

```go
url, err := r.Get("article").URL("category", "technology", "id", "42")
```

...and the result will be a `url.URL` with the following path:

```
"/articles/technology/42"
```

This also works for host and query value variables:

```go
r := mux.NewRouter()
r.Host("{subdomain}.example.com").
  Path("/articles/{category}/{id:[0-9]+}").
  Queries("filter", "{filter}").
  HandlerFunc(ArticleHandler).
  Name("article")

// url.String() will be "http://news.example.com/articles/technology/42?filter=gorilla"
url, err := r.Get("article").URL("subdomain", "news",
                                 "category", "technology",
                                 "id", "42",
                                 "filter", "gorilla")
```

All variables defined in the route are required, and their values must conform to the corresponding patterns. These requirements guarantee that a generated URL will always match a registered route -- the only exception is for explicitly defined "build-only" routes which never match.

Regex support also exists for matching Headers within a route. For example, we could do:

```go
r.HeadersRegexp("Content-Type", "application/(text|json)")
```

...and the route will match both requests with a Content-Type of `application/json` as well as `application/text`

There's also a way to build only the URL host or path for a route: use the methods `URLHost()` or `URLPath()` instead. For the previous route, we would do:

```go
// "http://news.example.com/"
host, err := r.Get("article").URLHost("subdomain", "news")

// "/articles/technology/42"
path, err := r.Get("article").URLPath("category", "technology", "id", "42")
```

And if you use subrouters, host and path defined separately can be built as well:

```go
r := mux.NewRouter()
s := r.Host("{subdomain}.example.com").Subrouter()
s.Path("/articles/{category}/{id:[0-9]+}").
  HandlerFunc(ArticleHandler).
  Name("article")

// "http://news.example.com/articles/technology/42"
url, err := r.Get("article").URL("subdomain", "news",
                                 "category", "technology",
                                 "id", "42")
```

### Walking Routes

The `Walk` function on `mux.Router` can be used to visit all of the routes that are registered on a router. For example,
the following prints all of the registered routes:

```go
package main

import (
	"fmt"
	"net/http"
	"strings"

	"github.com/gorilla/mux"
)

func handler(w http.ResponseWriter, r *http.Request) {
	return
}

func main() {
	r := mux.NewRouter()
	r.HandleFunc("/", handler)
	r.HandleFunc("/products", handler).Methods("POST")
	r.HandleFunc("/articles", handler).Methods("GET")
	r.HandleFunc("/articles/{id}", handler).Methods("GET", "PUT")
	r.HandleFunc("/authors", handler).Queries("surname", "{surname}")
	err := r.Walk(func(route *mux.Route, router *mux.Router, ancestors []*mux.Route) error {
		pathTemplate, err := route.GetPathTemplate()
		if err == nil {
			fmt.Println("ROUTE:", pathTemplate)
		}
		pathRegexp, err := route.GetPathRegexp()
		if err == nil {
			fmt.Println("Path regexp:", pathRegexp)
		}
		queriesTemplates, err := route.GetQueriesTemplates()
		if err == nil {
			fmt.Println("Queries templates:", strings.Join(queriesTemplates, ","))
		}
		queriesRegexps, err := route.GetQueriesRegexp()
		if err == nil {
			fmt.Println("Queries regexps:", strings.Join(queriesRegexps, ","))
		}
		methods, err := route.GetMethods()
		if err == nil {
			fmt.Println("Methods:", strings.Join(methods, ","))
		}
		fmt.Println()
		return nil
	})

	if err != nil {
		fmt.Println(err)
	}

	http.Handle("/", r)
}
```

### Graceful Shutdown

Go 1.8 introduced the ability to [gracefully shutdown](https://golang.org/doc/go1.8#http_shutdown) a `*http.Server`. Here's how to do that alongside `mux`:

```go
package main

import (
    "context"
    "flag"
    "log"
    "net/http"
    "os"
    "os/signal"
    "time"

    "github.com/gorilla/mux"
)

func main() {
    var wait time.Duration
    flag.DurationVar(&wait, "graceful-timeout", time.Second * 15, "the duration for which the server gracefully wait for existing connections to finish - e.g. 15s or 1m")
    flag.Parse()

    r := mux.NewRouter()
    // Add your routes as needed

    srv := &http.Server{
        Addr:         "0.0.0.0:8080",
        // Good practice to set timeouts to avoid Slowloris attacks.
        WriteTimeout: time.Second * 15,
        ReadTimeout:  time.Second * 15,
        IdleTimeout:  time.Second * 60,
        Handler: r, // Pass our instance of gorilla/mux in.
    }

    // Run our server in a goroutine so that it doesn't block.
    go func() {
        if err := srv.ListenAndServe(); err != nil {
            log.Println(err)
        }
    }()

    c := make(chan os.Signal, 1)
    // We'll accept graceful shutdowns when quit via SIGINT (Ctrl+C)
    // SIGKILL, SIGQUIT or SIGTERM (Ctrl+/) will not be caught.
    signal.Notify(c, os.Interrupt)

    // Block until we receive our signal.
    <-c

    // Create a deadline to wait for.
    ctx, cancel := context.WithTimeout(context.Background(), wait)
    defer cancel()
    // Doesn't block if no connections, but will otherwise wait
    // until the timeout deadline.
    srv.Shutdown(ctx)
    // Optionally, you could run srv.Shutdown in a goroutine and block on
    // <-ctx.Done() if your application should wait for other services
    // to finalize based on context cancellation.
    log.Println("shutting down")
    os.Exit(0)
}
```

### Middleware

Mux supports the addition of middlewares to a [Router](https://godoc.org/github.com/gorilla/mux#Router), which are executed in the order they are added if a match is found, including its subrouters.
Middlewares are (typically) small pieces of code which take one request, do something with it, and pass it down to another middleware or the final handler. Some common use cases for middleware are request logging, header manipulation, or `ResponseWriter` hijacking.

Mux middlewares are defined using the de facto standard type:

```go
type MiddlewareFunc func(http.Handler) http.Handler
```

Typically, the returned handler is a closure which does something with the http.ResponseWriter and http.Request passed to it, and then calls the handler passed as parameter to the MiddlewareFunc. This takes advantage of closures being able access variables from the context where they are created, while retaining the signature enforced by the receivers.

A very basic middleware which logs the URI of the request being handled could be written as:

```go
func loggingMiddleware(next http.Handler) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        // Do stuff here
        log.Println(r.RequestURI)
        // Call the next handler, which can be another middleware in the chain, or the final handler.
        next.ServeHTTP(w, r)
    })
}
```

Middlewares can be added to a router using `Router.Use()`:

```go
r := mux.NewRouter()
r.HandleFunc("/", handler)
r.Use(loggingMiddleware)
```

A more complex authentication middleware, which maps session token to users, could be written as:

```go
// Define our struct
type authenticationMiddleware struct {
	tokenUsers map[string]string
}

// Initialize it somewhere
func (amw *authenticationMiddleware) Populate() {
	amw.tokenUsers["00000000"] = "user0"
	amw.tokenUsers["aaaaaaaa"] = "userA"
	amw.tokenUsers["05f717e5"] = "randomUser"
	amw.tokenUsers["deadbeef"] = "user0"
}

// Middleware function, which will be called for each request
func (amw *authenticationMiddleware) Middleware(next http.Handler) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        token := r.Header.Get("X-Session-Token")

        if user, found := amw.tokenUsers[token]; found {
        	// We found the token in our map
        	log.Printf("Authenticated user %s\n", user)
        	// Pass down the request to the next middleware (or final handler)
        	next.ServeHTTP(w, r)
        } else {
        	// Write an error and stop the handler chain
        	http.Error(w, "Forbidden", http.StatusForbidden)
        }
    })
}
```

```go
r := mux.NewRouter()
r.HandleFunc("/", handler)

amw := authenticationMiddleware{}
amw.Populate()

r.Use(amw.Middleware)
```

Note: The handler chain will be stopped if your middleware doesn't call `next.ServeHTTP()` with the corresponding parameters. This can be used to abort a request if the middleware writer wants to. Middlewares _should_ write to `ResponseWriter` if they _are_ going to terminate the request, and they _should not_ write to `ResponseWriter` if they _are not_ going to terminate it.

### Handling CORS Requests

[CORSMethodMiddleware](https://godoc.org/github.com/gorilla/mux#CORSMethodMiddleware) intends to make it easier to strictly set the `Access-Control-Allow-Methods` response header.

* You will still need to use your own CORS handler to set the other CORS headers such as `Access-Control-Allow-Origin`
* The middleware will set the `Access-Control-Allow-Methods` header to all the method matchers (e.g. `r.Methods(http.MethodGet, http.MethodPut, http.MethodOptions)` -> `Access-Control-Allow-Methods: GET,PUT,OPTIONS`) on a route
* If you do not specify any methods, then:
> _Important_: there must be an `OPTIONS` method matcher for the middleware to set the headers.

Here is an example of using `CORSMethodMiddleware` along with a custom `OPTIONS` handler to set all the required CORS headers:

```go
package main

import (
	"net/http"
	"github.com/gorilla/mux"
)

func main() {
    r := mux.NewRouter()

    // IMPORTANT: you must specify an OPTIONS method matcher for the middleware to set CORS headers
    r.HandleFunc("/foo", fooHandler).Methods(http.MethodGet, http.MethodPut, http.MethodPatch, http.MethodOptions)
    r.Use(mux.CORSMethodMiddleware(r))
    
    http.ListenAndServe(":8080", r)
}

func fooHandler(w http.ResponseWriter, r *http.Request) {
    w.Header().Set("Access-Control-Allow-Origin", "*")
    if r.Method == http.MethodOptions {
        return
    }

    w.Write([]byte("foo"))
}
```

And an request to `/foo` using something like:

```bash
curl localhost:8080/foo -v
```

Would look like:

```bash
*   Trying ::1...
* TCP_NODELAY set
* Connected to localhost (::1) port 8080 (#0)
> GET /foo HTTP/1.1
> Host: localhost:8080
> User-Agent: curl/7.59.0
> Accept: */*
> 
< HTTP/1.1 200 OK
< Access-Control-Allow-Methods: GET,PUT,PATCH,OPTIONS
< Access-Control-Allow-Origin: *
< Date: Fri, 28 Jun 2019 20:13:30 GMT
< Content-Length: 3
< Content-Type: text/plain; charset=utf-8
< 
* Connection #0 to host localhost left intact
foo
```

### Testing Handlers

Testing handlers in a Go web application is straightforward, and _mux_ doesn't complicate this any further. Given two files: `endpoints.go` and `endpoints_test.go`, here's how we'd test an application using _mux_.

First, our simple HTTP handler:

```go
// endpoints.go
package main

func HealthCheckHandler(w http.ResponseWriter, r *http.Request) {
    // A very simple health check.
    w.Header().Set("Content-Type", "application/json")
    w.WriteHeader(http.StatusOK)

    // In the future we could report back on the status of our DB, or our cache
    // (e.g. Redis) by performing a simple PING, and include them in the response.
    io.WriteString(w, `{"alive": true}`)
}

func main() {
    r := mux.NewRouter()
    r.HandleFunc("/health", HealthCheckHandler)

    log.Fatal(http.ListenAndServe("localhost:8080", r))
}
```

Our test code:

```go
// endpoints_test.go
package main

import (
    "net/http"
    "net/http/httptest"
    "testing"
)

func TestHealthCheckHandler(t *testing.T) {
    // Create a request to pass to our handler. We don't have any query parameters for now, so we'll
    // pass 'nil' as the third parameter.
    req, err := http.NewRequest("GET", "/health", nil)
    if err != nil {
        t.Fatal(err)
    }

    // We create a ResponseRecorder (which satisfies http.ResponseWriter) to record the response.
    rr := httptest.NewRecorder()
    handler := http.HandlerFunc(HealthCheckHandler)

    // Our handlers satisfy http.Handler, so we can call their ServeHTTP method
    // directly and pass in our Request and ResponseRecorder.
    handler.ServeHTTP(rr, req)

    // Check the status code is what we expect.
    if status := rr.Code; status != http.StatusOK {
        t.Errorf("handler returned wrong status code: got %v want %v",
            status, http.StatusOK)
    }

    // Check the response body is what we expect.
    expected := `{"alive": true}`
    if rr.Body.String() != expected {
        t.Errorf("handler returned unexpected body: got %v want %v",
            rr.Body.String(), expected)
    }
}
```

In the case that our routes have [variables](#examples), we can pass those in the request. We could write
[table-driven tests](https://dave.cheney.net/2013/06/09/writing-table-driven-tests-in-go) to test multiple
possible route variables as needed.

```go
// endpoints.go
func main() {
    r := mux.NewRouter()
    // A route with a route variable:
    r.HandleFunc("/metrics/{type}", MetricsHandler)

    log.Fatal(http.ListenAndServe("localhost:8080", r))
}
```

Our test file, with a table-driven test of `routeVariables`:

```go
// endpoints_test.go
func TestMetricsHandler(t *testing.T) {
    tt := []struct{
        routeVariable string
        shouldPass bool
    }{
        {"goroutines", true},
        {"heap", true},
        {"counters", true},
        {"queries", true},
        {"adhadaeqm3k", false},
    }

    for _, tc := range tt {
        path := fmt.Sprintf("/metrics/%s", tc.routeVariable)
        req, err := http.NewRequest("GET", path, nil)
        if err != nil {
            t.Fatal(err)
        }

        rr := httptest.NewRecorder()
	
	// Need to create a router that we can pass the request through so that the vars will be added to the context
	router := mux.NewRouter()
        router.HandleFunc("/metrics/{type}", MetricsHandler)
        router.ServeHTTP(rr, req)

        // In this case, our MetricsHandler returns a non-200 response
        // for a route variable it doesn't know about.
        if rr.Code == http.StatusOK && !tc.shouldPass {
            t.Errorf("handler should have failed on routeVariable %s: got %v want %v",
                tc.routeVariable, rr.Code, http.StatusOK)
        }
    }
}
```

## Full Example

Here's a complete, runnable example of a small `mux` based server:

```go
package main

import (
    "net/http"
    "log"
    "github.com/gorilla/mux"
)

func YourHandler(w http.ResponseWriter, r *http.Request) {
    w.Write([]byte("Gorilla!\n"))
}

func main() {
    r := mux.NewRouter()
    // Routes consist of a path and a handler function.
    r.HandleFunc("/", YourHandler)

    // Bind to a port and pass our router in
    log.Fatal(http.ListenAndServe(":8000", r))
}
```

## License

BSD licensed. See the LICENSE file for details.

1) comprehensive mode

In comprehensive mode incoming json treated as much versatile as possible to 
match the given jsonpath. Examples: strings that represent integers can be
 treated as such to access array elements by index, dot-notated numeric keys can select array elements as well as object keys, depending on the actual input. This mode tries to emulate JavaScript behaviour. 

path | input 1 | result 1<br>explanation | input 2 | result 2<br>explanation
--- | --- | --- | --- | ---
`$[2]` or `$['2']` or `$.2` or `$.'2'`  | `["a","b","c"]` | `"c"`<br>array element by index |  `{ "1":"a", "2":"b" }` | `"b"`<br>object key by name
`$.*`<br>`$[*]`<br>`$[:]`  | `["a","b","c"]` | `["a","b","c"]`<br>all elements of an array |  `{ "1":"a", "2":"b" }` | `["a", "b"]`<br>values of all the keys
`$[*].bar`<br>`$.*.bar` | `[{"foo":1},{"bar":2}]` | `[2]`<br>find a key in every element |  `{"a":{"foo":1},"b":{"bar":2}}` | `[2]`<br>find a key in every value
`$[1,2]`<br>`$['1','2']` | `["a","b","c"]` | `["b","c"]`<br>aggregate array elements by index |  `{"1":"a", "2":"b"}}` | `["a","b"]`<br>aggregate values of keys by name

2) strict mode

This mode is useful when you need to be thorough about the input data. In strict mode dot notation is used solely for accessing object fields by name thus allowing to distinguish between objects and arrays in ambiguous cases, while bracket notation is used for:  
a) indexing or aggregating array elements -- in this case all values inside brackets must be numeric  
b) aggregating key values in objects -- in this case all values inside brackets must be alphanumeric or quoted strings.  
c) selecting a single key value from object -- this does not imply aggregation and treated as a synonym for dot notation (for compatibility reasons).  
Similarly, dot-notated wildcard is only applicable to objects.

path | input 1 | result 1<br>explanation | input 2 | result 2<br>explanation
--- | --- | --- | --- | ---
`$[2]`  | `["a","b","c"]` | `"c"`<br>array element by index |  `{ "1":"a", "2":"b" }` | <span style="color:#DD4444">not applicable:<br>trying to index an object</span>
`$.2` or `$.'2'` | `["a","b","c"]` | <span style="color:#DD4444">dot notation is not applicable for array</span> |  `{ "1":"a", "2":"b" }` | `"b"`<br>object key by name
`$['2']` | `["a","b","c"]` | <span style="color:#DD4444">not applicable:<br>non-integer index</span> |  `{ "1":"a", "2":"b" }` | `"b"`<br>object key by name
`$['1','2']` | `["a","b","c"]` | <span style="color:#DD4444">not applicable:<br>non-integer indexes</span> |  `{"1":"a", "2":"b"}}` | `["a","b"]`<br>aggregate values of keys by name
`$[1,2]` | `["a","b","c"]` | `["b","c"]`<br>aggregate array elements by index |  `{"1":"a", "2":"b"}}` | `["a","b"]`<br>aggregate values of keys by name
`$[*]`  | `["a","b","c"]` | `["a","b","c"]`<br>all elements of an array |  `{ "1":"a", "2":"b" }` | `["a", "b"]`<br>values of all the keys (aggregation)
`$.*`  | `["a","b","c"]` | <span style="color:#DD4444">dot notation is not applicable for array</span> |  `{ "1":"a", "2":"b" }` | `["a", "b"]`<br>values of all the keys (aggregation)

### Schema

```jsonpath:
	${ref}{ref}...
	@{ref}{ref}...
ref:
	.|..{keyref}
	|.|..|{brackets}
keyref:
	key|brackets
brackets:
	[{someth}]
key:
	string
	word
	index
	*
someth:
	?({expr})
	{key}
	{key},{key}...
	{index},{index}...
	{start}:{end}
	{start}:{end}:{step}
index, start, end, step:
	integer
expr:
	{operand} {operator}{operand}
operand:
	jsonpath
	value
value:
	string
	integer
	float
	'null'
	regexp
operator:
	==,!=,>,<,>=,<=,&&,||
string:
	"..."
	'...'
	`...`
```

ref types:
example | type | applicable to | flags | notes (NF = not found)
--- | --- | --- | --- | ---
`$.key` or `$.'key'` | single word key | object | **common** | NF on arrays
`$.3` or `$[3]` or `$['3']` or `$.[3]` or `$.['3']` | single numeric key == index | object or array | **common**
`$[1,2]` | union | object or array | **aggregating** | 
`$[1,'a']` | union | object or array | **aggregating** | word keys NF on arrays
`$[1,'a']` | union | object or array | **aggregating** | word keys NF on arrays
`$['*']` | == `$.'*'` | object | **common** | syntax to get a value of a `*` key
`$.key.size()` | function | object or array | **function** | 
`$[xx:yy:zz]` | slice | array | **slice** | NF on objects
`$[:]` | slice | array | **slice** | == `$.*` in comprehensive mode (?)
`$..key` or `$..['key']` | sigle word key | object or array | **deepscan** | 
`$..[0]` or `$..['0']` | sigle word key == index | object or array | **deepscan** | 
`$..[0:2]` |  | array | **deepscan** | 
`$.*` or `$[*]` | wildcard | object or array | **wildcard**
`$..*` or `$..[*]` | deepscan wildcard | object or array | **deepscan** **wildcard**

- common ($.book or $[0]) -- cDot
	- array
		- word ref: NF
		- index: by index
	- object
		- word ref: by name
		- index: NF
	- base type
		- NF
- aggregating ($[1,2] or $['a','b']) -- cAgg
	- array
		- word ref: NF
		- index: AGG( by index )
	- object
		- word ref: AGG( by name )
		- index: NF
- slice ($[0:3]) -- cSlice
	- array 
		- AGG( slice elems )
	- object
		- NF
- function -- cFunction
	- array
		- length or size
	- string
		- string length
- wildcard (.* or [:] or [*]) -- cWild
	- array
		- AGG( all elements )
	- object
		- AGG( values of all keys )
- deepscan ($..book or $..[0] or $..['foo','bar']) -- cDeep
	- array
		- word ref: AGG( recurse on every elem )
		- index: AGG( get elem by index + recurse on every elem )
	- object
		- word ref: AGG( by name + recurse on every kvalue )
		- index: AGG( recurse on every kvalue )
- deepscan wildcard ($..* or $..[*] or $..[:])
	- array
		- AGG( all elems + recurse on every elem )
	- object
		- word ref: AGG( all kvalues + recurse on every kvalue )

ref flags:

	- common		- common node
	- terminal		- no more refs follow, return result
	- union
		- object: collect values of keys, return array
		- array: collect specified elems, return array
	- function		- apply function to the last value
	- slice			- slice elems (the subject must be array)
	- filter		- apply filter (the subject must be array)
	- wildcard		- wildcard for object or array. Result is array.
	- deepscan		- deepscan. Result is array

array modes:

	single
		[1]     - counting scan 
		[-1]    - fullscan & cut
	
	multi
		[1,2]   - counting scan 
		[-1,2]  - fullscan & cut
	
	slice
		[:]     - [cEmpty:cEmpty] fullscan & cut
		[2:]    - [INT:cEmpty]    fullscan & cut
		[:5]    - [cEmpty:INT]    fullscan & cut
		[-2:]   - [-INT:cEmpty]   fullscan & cut
		[:-5]   - [cEmpty:-INT]   fullscan & cut
		[2:5]   - [INT:INT]       counting scan
		[-5:2]  - [-INT:INT]      fullscan & cut
		[1:-3]  - [INT:-INT]      fullscan & cut
		[-5:-3] - [-INT:-INT]     fullscan & cut

	slice with step
		fullscan & cut
[![Build Status](https://travis-ci.org/bhmj/jsonslice.svg?branch=master)](https://travis-ci.org/bhmj/jsonslice)
[![Go Report Card](https://goreportcard.com/badge/github.com/bhmj/jsonslice)](https://goreportcard.com/report/github.com/bhmj/jsonslice)
[![GoDoc](https://godoc.org/github.com/bhmj/jsonslice?status.svg)](https://godoc.org/github.com/bhmj/jsonslice)

# JSON Slice

## What is it?

JSON Slice is a Go package which allows to execute fast jsonpath queries without unmarshalling the whole data.  

Sometimes you need to get a single value from incoming json using jsonpath, for example to route data accordingly or so. To do that you must unmarshall the whole data into interface{} struct and then apply some jsonpath library to it, only to get just a tiny little value. What a waste of resourses! Well, now there's `jsonslice`.

Simply call `jsonslice.Get` on your raw json data to slice out just the part you need. The `[]byte` received can then be unmarshalled into a struct or used as it is.

## Getting started

#### 1. install

```
$ go get github.com/bhmj/jsonslice
```

#### 2. use it

```golang
import "github.com/bhmj/jsonslice"
import "fmt"

func main() {
    var data = []byte(`
    { "sku": [ 
        { "id": 1, "name": "Bicycle", "price": 160, "extras": [ "flashlight", "pump" ] },
        { "id": 2, "name": "Scooter", "price": 280, "extras": [ "helmet", "gloves", "spare wheel" ] }
      ]
    } `)

    a, _ := jsonslice.Get(data, "$.sku[0].price")
    b, _ := jsonslice.Get(data, "$.sku[1].extras.count()")
    c, _ := jsonslice.Get(data, "$.sku[?(@.price > 200)].name")
    d, _ := jsonslice.Get(data, "$.sku[?(@.extras.count() < 3)].name")

    fmt.Println(string(a)) // 160
    fmt.Println(string(b)) // 3
    fmt.Println(string(c)) // ["Scooter"]
    fmt.Println(string(d)) // ["Bicycle"]
}
```
[Run in Go Playground](https://play.golang.org/p/fYv-Y12akvs)

## Package functions
  
`jsonslice.Get(data []byte, jsonpath string) ([]byte, error)`  
  - get a slice from raw json data specified by jsonpath

## Specs

See [Stefan Gössner's article](http://goessner.net/articles/JsonPath/index.html#e2) for original specs and examples.  

## Syntax features

1. Classic dot notation (`$.simple_key`) is limited to alphanumeric characters. For more complex cases use `$['complex key!']` or `$.'complex key!'`. 

2. A single index reference returns an element, not an array; a slice reference always returns array:  
```
> echo '[{"id":1}, {"id":2}]' | ./jsonslice '$[0].id' 
1
```
```
> echo '[{"id":1}, {"id":2}]' | ./jsonslice '$[0:1].id'
[1]
```

3. Indexing or slicing on root node is supported (assuming json is an array and not an object):  
```
./jsonslice '$[0].author' sample1.json
```

## Expressions

### Overview 
```
  $                   -- root node (can be either object or array)
  .node               -- dot-notated child
  .'some node'        -- dot-notated child (syntax extension)
  ['node']            -- bracket-notated child
  ['foo','bar']       -- bracket-notated children (aggregation)
  [5]                 -- array index
  [-5]                -- negative index means "from the end"
  [1:9]               -- array slice
  [1:9:2]             -- array slice (+step)
  .*  .[*]  .[:]      -- wildcard
  ..key               -- deepscan
```
### Functions
```
  $.obj.length()      -- number of elements in an array or string length, depending on the obj type
  $.obj.count()       -- same as above
  $.val.size()        -- value size in bytes (as is)
```
### Slices
```
  $.arr[start:end:step]
  $.arr[start:end]
```
Selects elements from `start` (inclusive) to `end` (exclusive), stepping by `step`. If `step` is omitted or zero, then 1 is implied. Out-of-bounds values are reduced to the nearest bound.

If `step` is positive:
  - empty `start` treated as the first element inclusive
  - empty `end` treated as the last element inclusive
  - `start` should be less then `end`, otherwise result will be empty

If `step` is negative:
  - empty `start` treated as last element inclusive
  - empty `end` treated as the first element inclusive
  - `start` should be greater then `end`, otherwise result will be empty

### Filters

```
  [?(<expression>)]  -- filter expression. Applicable to arrays only
  @                  -- the root of the current element of the array. Used only within a filter.
  @.val              -- a field of the current element of the array.
```

#### Filter operators

  Operator | Description
  --- | ---
  `==`  | Equal to<br>Use single or double quotes for string expressions.<br>`[?(@.color=='red')]` or `[?(@.color=="red")]`
  `!=`  | Not equal to<br>`[?(@.author != "Herman Melville")]`
  `>`   | Greater than<br>`[?(@.price > 10)]`
  `>=`  | Greater than or equal to
  `<`   | Less than
  `<=`  | Less than or equal to
  `=~`  | Match a regexp<br>`[?(@.name =~ /sword.*/i]`
  `!~` or `!=~`  | Don't match a regexp<br>`[?(@.name !~ /sword.*/i]`
  `&&`  | Logical AND<br>`[?(@.price < 10 && @isbn)]`
  `\|\|`  | Logical OR<br>`[?(@.price > 10 \|\| @.category == 'reference')]`

## Examples

Assuming `sample0.json` and `sample1.json` in the example directory:  

  `cat sample0.json | ./jsonslice '$.store.book[0]'`  
  `cat sample0.json | ./jsonslice '$.store.book[0].title'`  
  `cat sample0.json | ./jsonslice '$.store.book[0:-1]'`  
  `cat sample1.json | ./jsonslice '$[1].author'`  
  `cat sample0.json | ./jsonslice '$.store.book[?(@.price > 10)]'`  
  `cat sample0.json | ./jsonslice '$.store.book[?(@.price > $.expensive)]'`  

Much more examples can be found in `jsonslice_test.go`  

## Benchmarks (Core i5-7500)

```diff
$ go test -bench=. -benchmem -benchtime=4s
goos: linux
goarch: amd64
pkg: github.com/bhmj/jsonslice
++ usually you need to unmarshal the whole JSON to get an object by jsonpath (for reference):
Benchmark_Unmarshal-4                     500000             14712 ns/op            4368 B/op        130 allocs/op
++ and here's a jsonslice.Get:
Benchmark_Jsonslice_Get_Simple-4         2000000              3878 ns/op             128 B/op          4 allocs/op
++ Get() involves parsing a jsonpath, here it is:
Benchmark_JsonSlice_ParsePath-4         10000000               858 ns/op             160 B/op          5 allocs/op
++ in case you aggregate some non-contiguous elements, it may take a bit longer (extra mallocs involved):
Benchmark_Jsonslice_Get_Aggregated-4     1000000              5671 ns/op             417 B/op         13 allocs/op
++ usual unmarshalling a large json:
Benchmark_Unmarshal_10Mb-4                   100          50744817 ns/op             248 B/op          5 allocs/op
++ jsonslicing the same json, target element is near the start:
Benchmark_Jsonslice_Get_10Mb_First-4     3000000              1851 ns/op             128 B/op          4 allocs/op
++ jsonslicing the same json, target element is near the end: still beats Unmarshal
Benchmark_Jsonslice_Get_10Mb_Last-4          200          38286509 ns/op             133 B/op          4 allocs/op
PASS
ok      github.com/bhmj/jsonslice       83.152s

```

## Changelog

**1.0.4** (2020-05-07) -- bugfix: `$*` path generated panic.

**1.0.3** (2019-12-24) -- bugfix: `$[0].foo` `[{"foo":"\\"}]` generated "unexpected end of input".

**1.0.2** (2019-12-07) -- nested aggregation (`$[:].['a','b']`) now works as expected. TODO: add option to switch nested aggregation mode at runtime!

**1.0.1** (2019-12-01) -- "not equal" regexp operator added (`!=~` or `!~`).

**1.0.0** (2019-11-29) -- deepscan operator (`..`) added, slice with step `$[1:9:2]` is now supported, syntax extensions added. `GetArrayElements()` removed.

**0.7.6** (2019-09-11) -- bugfix: escaped backslash at the end of a string value.

**0.7.5** (2019-05-21) -- Functions `count()`, `size()`, `length()` work in filters.
> `$.store.bicycle.equipment[?(@.count() == 2)]` -> `[["light saber", "apparel"]]`  

**0.7.4** (2019-03-01) -- Mallocs reduced (see Benchmarks section).

**0.7.3** (2019-02-27) -- `GetArrayElements()` added.

**0.7.2** (2018-12-25) -- bugfix: closing square bracket inside a string value.

**0.7.1** (2018-10-16) -- bracket notation is now supported.
> `$.store.book[:]['price','title']` -> `[[8.95,"Sayings of the Century"],[12.99,"Sword of Honour"],[8.99,"Moby Dick"],[22.99,"The Lord of the Rings"]]`

**0.7.0** (2018-07-23) -- Wildcard key (`*`) added.
> `$.store.book[-1].*` -> `["fiction","J. R. R. Tolkien","The Lord of the Rings","0-395-19395-8",22.99]`  
> `$.store.*[:].price` -> `[8.95,12.99,8.99,22.99]`

**0.6.3** (2018-07-16) -- Boolean/null value error fixed.

**0.6.2** (2018-07-03) -- More tests added, error handling clarified.

**0.6.1** (2018-06-26) -- Nested array indexing is now supported.
> `$.store.bicycle.equipment[1][0]` -> `"peg leg"`

**0.6.0** (2018-06-25) -- Regular expressions added.
> `$.store.book[?(@.title =~ /(dick)|(lord)/i)].title` -> `["Moby Dick","The Lord of the Rings"]`

**0.5.1** (2018-06-15) -- Logical expressions added.
> `$.store.book[?(@.price > $.expensive && @.isbn)].title` -> `["The Lord of the Rings"]`

**0.5.0** (2018-06-14) -- Expressions (aka filters) added.
> `$.store.book[?(@.price > $.expensive)].title` -> `["Sword of Honour","The Lord of the Rings"]`

**0.4.0** (2018-05-16) -- Aggregating sub-queries added.
> `$.store.book[1:3].author` -> `["John","William"]`

**0.3.0** (2018-05-05) -- MVP.

## Roadmap

- [x] length(), count(), size() functions
- [x] filters: simple expressions
- [x] filters: complex expressions (with logical operators)
- [x] nested arrays support
- [x] wildcard operator (`*`)
- [x] bracket notation for multiple field queries (`$['a','b']`)
- [x] deepscan operator (`..`)
- [x] syntax extensions: `$.'keys with spaces'.price`
- [x] flexible syntax: `$[0]` works on both `[1,2,3]` and `{"0":"abc"}`
- [ ] Optionally unmarshal the result
- [ ] Option to select aggregation mode (nested or plain)

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :)

## Licence

[MIT](http://opensource.org/licenses/MIT)

## Author

Michael Gurov aka BHMJ
# errors [![Travis-CI](https://travis-ci.org/pkg/errors.svg)](https://travis-ci.org/pkg/errors) [![AppVeyor](https://ci.appveyor.com/api/projects/status/b98mptawhudj53ep/branch/master?svg=true)](https://ci.appveyor.com/project/davecheney/errors/branch/master) [![GoDoc](https://godoc.org/github.com/pkg/errors?status.svg)](http://godoc.org/github.com/pkg/errors) [![Report card](https://goreportcard.com/badge/github.com/pkg/errors)](https://goreportcard.com/report/github.com/pkg/errors) [![Sourcegraph](https://sourcegraph.com/github.com/pkg/errors/-/badge.svg)](https://sourcegraph.com/github.com/pkg/errors?badge)

Package errors provides simple error handling primitives.

`go get github.com/pkg/errors`

The traditional error handling idiom in Go is roughly akin to
```go
if err != nil {
        return err
}
```
which applied recursively up the call stack results in error reports without context or debugging information. The errors package allows programmers to add context to the failure path in their code in a way that does not destroy the original value of the error.

## Adding context to an error

The errors.Wrap function returns a new error that adds context to the original error. For example
```go
_, err := ioutil.ReadAll(r)
if err != nil {
        return errors.Wrap(err, "read failed")
}
```
## Retrieving the cause of an error

Using `errors.Wrap` constructs a stack of errors, adding context to the preceding error. Depending on the nature of the error it may be necessary to reverse the operation of errors.Wrap to retrieve the original error for inspection. Any error value which implements this interface can be inspected by `errors.Cause`.
```go
type causer interface {
        Cause() error
}
```
`errors.Cause` will recursively retrieve the topmost error which does not implement `causer`, which is assumed to be the original cause. For example:
```go
switch err := errors.Cause(err).(type) {
case *MyError:
        // handle specifically
default:
        // unknown error
}
```

[Read the package documentation for more information](https://godoc.org/github.com/pkg/errors).

## Roadmap

With the upcoming [Go2 error proposals](https://go.googlesource.com/proposal/+/master/design/go2draft.md) this package is moving into maintenance mode. The roadmap for a 1.0 release is as follows:

- 0.9. Remove pre Go 1.9 and Go 1.10 support, address outstanding pull requests (if possible)
- 1.0. Final release.

## Contributing

Because of the Go2 errors changes, this package is not accepting proposals for new functionality. With that said, we welcome pull requests, bug fixes and issue reports. 

Before sending a PR, please discuss your change by raising an issue.

## License

BSD-2-Clause
XPath
====
[![GoDoc](https://godoc.org/github.com/antchfx/xpath?status.svg)](https://godoc.org/github.com/antchfx/xpath)
[![Coverage Status](https://coveralls.io/repos/github/antchfx/xpath/badge.svg?branch=master)](https://coveralls.io/github/antchfx/xpath?branch=master)
[![Build Status](https://travis-ci.org/antchfx/xpath.svg?branch=master)](https://travis-ci.org/antchfx/xpath)
[![Go Report Card](https://goreportcard.com/badge/github.com/antchfx/xpath)](https://goreportcard.com/report/github.com/antchfx/xpath)

XPath is Go package provides selecting nodes from XML, HTML or other documents using XPath expression.

Implementation
===

- [htmlquery](https://github.com/antchfx/htmlquery) - an XPath query package for HTML document

- [xmlquery](https://github.com/antchfx/xmlquery) - an XPath query package for XML document.

- [jsonquery](https://github.com/antchfx/jsonquery) - an XPath query package for JSON document

Supported Features
===

#### The basic XPath patterns.

> The basic XPath patterns cover 90% of the cases that most stylesheets will need.

- `node` : Selects all child elements with nodeName of node.

- `*` : Selects all child elements.

- `@attr` : Selects the attribute attr.

- `@*` : Selects all attributes.

- `node()` : Matches an org.w3c.dom.Node.

- `text()` : Matches a org.w3c.dom.Text node.

- `comment()` : Matches a comment.

- `.` : Selects the current node.

- `..` : Selects the parent of current node.

- `/` : Selects the document node.

- `a[expr]` : Select only those nodes matching a which also satisfy the expression expr.

- `a[n]` : Selects the nth matching node matching a When a filter's expression is a number, XPath selects based on position.

- `a/b` : For each node matching a, add the nodes matching b to the result.

- `a//b` : For each node matching a, add the descendant nodes matching b to the result. 

- `//b` : Returns elements in the entire document matching b.

- `a|b` : All nodes matching a or b, union operation(not boolean or).

- `(a, b, c)` : Evaluates each of its operands and concatenates the resulting sequences, in order, into a single result sequence


#### Node Axes 

- `child::*` : The child axis selects children of the current node.

- `descendant::*` : The descendant axis selects descendants of the current node. It is equivalent to '//'.

- `descendant-or-self::*` : Selects descendants including the current node.

- `attribute::*` : Selects attributes of the current element. It is equivalent to @*

- `following-sibling::*` : Selects nodes after the current node.

- `preceding-sibling::*` : Selects nodes before the current node.

- `following::*` : Selects the first matching node following in document order, excluding descendants. 

- `preceding::*` : Selects the first matching node preceding in document order, excluding ancestors. 

- `parent::*` : Selects the parent if it matches. The '..' pattern from the core is equivalent to 'parent::node()'.

- `ancestor::*` : Selects matching ancestors.

- `ancestor-or-self::*` : Selects ancestors including the current node.

- `self::*` : Selects the current node. '.' is equivalent to 'self::node()'.

#### Expressions

 The gxpath supported three types: number, boolean, string.

- `path` : Selects nodes based on the path.

- `a = b` : Standard comparisons.

    * a = b	    True if a equals b.
    * a != b	True if a is not equal to b.
    * a < b	    True if a is less than b.
    * a <= b	True if a is less than or equal to b.
    * a > b	    True if a is greater than b.
    * a >= b	True if a is greater than or equal to b.

- `a + b` : Arithmetic expressions.

    * `- a`	Unary minus
    * a + b	Add
    * a - b	Substract
    * a * b	Multiply
    * a div b	Divide
    * a mod b	Floating point mod, like Java.

- `a or b` : Boolean `or` operation.

- `a and b` : Boolean `and` operation.

- `(expr)` : Parenthesized expressions.

- `fun(arg1, ..., argn)` : Function calls:

| Function | Supported |
| --- | --- |
`boolean()`| ✓ |
`ceiling()`| ✓ |
`choose()`| ✗ |
`concat()`| ✓ |
`contains()`| ✓ |
`count()`| ✓ |
`current()`| ✗ |
`document()`| ✗ |
`element-available()`| ✗ |
`ends-with()`| ✓ |
`false()`| ✓ |
`floor()`| ✓ |
`format-number()`| ✗ |
`function-available()`| ✗ |
`generate-id()`| ✗ |
`id()`| ✗ |
`key()`| ✗ |
`lang()`| ✗ |
`last()`| ✓ |
`local-name()`| ✓ |
`name()`| ✓ |
`namespace-uri()`| ✓ |
`normalize-space()`| ✓ |
`not()`| ✓ |
`number()`| ✓ |
`position()`| ✓ |
`replace()`| ✓ |
`reverse()`| ✓ |
`round()`| ✓ |
`starts-with()`| ✓ |
`string()`| ✓ |
`string-length()`| ✓ |
`substring()`| ✓ |
`substring-after()`| ✓ |
`substring-before()`| ✓ |
`sum()`| ✓ |
`system-property()`| ✗ |
`translate()`| ✓ |
`true()`| ✓ |
`unparsed-entity-url()` | ✗ |

Changelogs
===

2019-03-19 
- optimize XPath `|` operation performance. [#33](https://github.com/antchfx/xpath/issues/33). Tips: suggest split into multiple subquery if you have a lot of `|` operations.

2019-01-29
-  improvement `normalize-space` function. [#32](https://github.com/antchfx/xpath/issues/32)

2018-12-07
-  supports XPath 2.0 Sequence expressions. [#30](https://github.com/antchfx/xpath/pull/30) by [@minherz](https://github.com/minherz).jsonquery
====
[![Build Status](https://travis-ci.org/antchfx/jsonquery.svg?branch=master)](https://travis-ci.org/antchfx/jsonquery)
[![Coverage Status](https://coveralls.io/repos/github/antchfx/jsonquery/badge.svg?branch=master)](https://coveralls.io/github/antchfx/jsonquery?branch=master)
[![GoDoc](https://godoc.org/github.com/antchfx/jsonquery?status.svg)](https://godoc.org/github.com/antchfx/jsonquery)
[![Go Report Card](https://goreportcard.com/badge/github.com/antchfx/jsonquery)](https://goreportcard.com/report/github.com/antchfx/jsonquery)

Overview
===

jsonquery is an XPath query package for JSON document, lets you extract data from JSON documents through an XPath expression. Built-in XPath expression cache avoid re-compile XPath expression each query.

Getting Started
===

### Install Package
```
go get github.com/antchfx/jsonquery
```

#### Load JSON document from URL.

```go
doc, err := jsonquery.LoadURL("http://www.example.com/feed?json")
```

#### Load JSON document from string.

```go
s :=`{
    "name":"John",
    "age":31, 
    "city":"New York" 
    }`
doc, err := jsonquery.Parse(strings.NewReader(s))
```

#### Load JSON document from io.Reader.

```go
f, err := os.Open("./books.json")
doc, err := jsonquery.Parse(f)
```

#### Find authors of all books in the store.
```go
list := jsonquery.Find(doc, "store/book/*/author")
// or equal to
list := jsonquery.Find(doc, "//author")
// or by QueryAll()
nodes, err := jsonquery.QueryAll(doc, "//a")
```

#### Find the third book.

```go
book := jsonquery.Find(doc, "//book/*[3]")
```

#### Find the last book.

```go
book := jsonquery.Find(doc, "//book/*[last()]")
```

#### Find all books that have an isbn number.

```go
list := jsonquery.Find(doc, "//book/*[isbn]")
```

#### Find all books priced less than 10.

```go
list := jsonquery.Find(doc, "//book/*[price<10]")
```

Examples
===

```go
func main() {
	s := `{
		"name": "John",
		"age"      : 26,
		"address"  : {
		  "streetAddress": "naist street",
		  "city"         : "Nara",
		  "postalCode"   : "630-0192"
		},
		"phoneNumbers": [
		  {
			"type"  : "iPhone",
			"number": "0123-4567-8888"
		  },
		  {
			"type"  : "home",
			"number": "0123-4567-8910"
		  }
		]
	}`
	doc, err := jsonquery.Parse(strings.NewReader(s))
	if err != nil {
		panic(err)
	}
	name := jsonquery.FindOne(doc, "name")
	fmt.Printf("name: %s\n", name.InnerText())
	var a []string
	for _, n := range jsonquery.Find(doc, "phoneNumbers/*/number") {
		a = append(a, n.InnerText())
	}
	fmt.Printf("phone number: %s\n", strings.Join(a, ","))
	if n := jsonquery.FindOne(doc, "address/streetAddress"); n != nil {
		fmt.Printf("address: %s\n", n.InnerText())
	}
}
```

Implement Principle
===
If you are familiar with XPath and XML, you can easily figure out how to
write your XPath expression.

```json
{
"name":"John",
"age":30,
"cars": [
	{ "name":"Ford", "models":[ "Fiesta", "Focus", "Mustang" ] },
	{ "name":"BMW", "models":[ "320", "X3", "X5" ] },
	{ "name":"Fiat", "models":[ "500", "Panda" ] }
]
}
```
The above JSON document will be convert to similar to XML document by the *JSONQuery*, like below:

```XML
<name>John</name>
<age>30</age>
<cars>
	<element>
		<name>Ford</name>
		<models>
			<element>Fiesta</element>
			<element>Focus</element>
			<element>Mustang</element>
		</models>		
	</element>
	<element>
		<name>BMW</name>
		<models>
			<element>320</element>
			<element>X3</element>
			<element>X5</element>
		</models>		
	</element>
	<element>
		<name>Fiat</name>
		<models>
			<element>500</element>
			<element>Panda</element>
		</models>		
	</element>
</cars>
```

Notes: `element` is empty element that have no any name.

List of XPath query packages
===
|Name |Description |
|--------------------------|----------------|
|[htmlquery](https://github.com/antchfx/htmlquery) | XPath query package for the HTML document|
|[xmlquery](https://github.com/antchfx/xmlquery) | XPath query package for the XML document|
|[jsonquery](https://github.com/antchfx/jsonquery) | XPath query package for the JSON document|
# Alice 

[![GoDoc](https://godoc.org/github.com/golang/gddo?status.svg)](http://godoc.org/github.com/justinas/alice)
[![Build Status](https://travis-ci.org/justinas/alice.svg?branch=master)](https://travis-ci.org/justinas/alice)
[![Coverage](http://gocover.io/_badge/github.com/justinas/alice)](http://gocover.io/github.com/justinas/alice)

Alice provides a convenient way to chain 
your HTTP middleware functions and the app handler.

In short, it transforms

    Middleware1(Middleware2(Middleware3(App)))

to

    alice.New(Middleware1, Middleware2, Middleware3).Then(App)

### Why?

None of the other middleware chaining solutions
behaves exactly like Alice.
Alice is as minimal as it gets:
in essence, it's just a for loop that does the wrapping for you.

Check out [this blog post](http://justinas.org/alice-painless-middleware-chaining-for-go/)
for explanation how Alice is different from other chaining solutions.

### Usage

Your middleware constructors should have the form of

    func (http.Handler) http.Handler

Some middleware provide this out of the box.
For ones that don't, it's trivial to write one yourself.

```go
func myStripPrefix(h http.Handler) http.Handler {
    return http.StripPrefix("/old", h)
}
```

This complete example shows the full power of Alice.

```go
package main

import (
    "net/http"
    "time"

    "github.com/throttled/throttled"
    "github.com/justinas/alice"
    "github.com/justinas/nosurf"
)

func timeoutHandler(h http.Handler) http.Handler {
    return http.TimeoutHandler(h, 1*time.Second, "timed out")
}

func myApp(w http.ResponseWriter, r *http.Request) {
    w.Write([]byte("Hello world!"))
}

func main() {
    th := throttled.Interval(throttled.PerSec(10), 1, &throttled.VaryBy{Path: true}, 50)
    myHandler := http.HandlerFunc(myApp)

    chain := alice.New(th.Throttle, timeoutHandler, nosurf.NewPure).Then(myHandler)
    http.ListenAndServe(":8000", chain)
}
```

Here, the request will pass [throttled](https://github.com/PuerkitoBio/throttled) first,
then an http.TimeoutHandler we've set up,
then [nosurf](https://github.com/justinas/nosurf)
and will finally reach our handler.

Note that Alice makes **no guarantees** for
how one or another piece of  middleware will behave.
It executes all middleware sequentially so that if a
piece of middleware were to stop the chain,
the request will not reach the inner handlers.
This is intentional behavior.

Alice works with Go 1.0 and higher,
but running tests requires at least Go 1.1.

### Contributing

0. Find an issue that bugs you / open a new one.
1. Discuss.
2. Branch off, commit, test.
3. Make a pull request / attach the commits to the issue.
## `jwt-go` Version History

#### 3.0.0

* **Compatibility Breaking Changes**: See MIGRATION_GUIDE.md for tips on updating your code
	* Dropped support for `[]byte` keys when using RSA signing methods.  This convenience feature could contribute to security vulnerabilities involving mismatched key types with signing methods.
	* `ParseFromRequest` has been moved to `request` subpackage and usage has changed
	* The `Claims` property on `Token` is now type `Claims` instead of `map[string]interface{}`.  The default value is type `MapClaims`, which is an alias to `map[string]interface{}`.  This makes it possible to use a custom type when decoding claims.
* Other Additions and Changes
	* Added `Claims` interface type to allow users to decode the claims into a custom type
	* Added `ParseWithClaims`, which takes a third argument of type `Claims`.  Use this function instead of `Parse` if you have a custom type you'd like to decode into.
	* Dramatically improved the functionality and flexibility of `ParseFromRequest`, which is now in the `request` subpackage
	* Added `ParseFromRequestWithClaims` which is the `FromRequest` equivalent of `ParseWithClaims`
	* Added new interface type `Extractor`, which is used for extracting JWT strings from http requests.  Used with `ParseFromRequest` and `ParseFromRequestWithClaims`.
	* Added several new, more specific, validation errors to error type bitmask
	* Moved examples from README to executable example files
	* Signing method registry is now thread safe
	* Added new property to `ValidationError`, which contains the raw error returned by calls made by parse/verify (such as those returned by keyfunc or json parser)

#### 2.7.0

This will likely be the last backwards compatible release before 3.0.0, excluding essential bug fixes.

* Added new option `-show` to the `jwt` command that will just output the decoded token without verifying
* Error text for expired tokens includes how long it's been expired
* Fixed incorrect error returned from `ParseRSAPublicKeyFromPEM`
* Documentation updates

#### 2.6.0

* Exposed inner error within ValidationError
* Fixed validation errors when using UseJSONNumber flag
* Added several unit tests

#### 2.5.0

* Added support for signing method none.  You shouldn't use this.  The API tries to make this clear.
* Updated/fixed some documentation
* Added more helpful error message when trying to parse tokens that begin with `BEARER `

#### 2.4.0

* Added new type, Parser, to allow for configuration of various parsing parameters
	* You can now specify a list of valid signing methods.  Anything outside this set will be rejected.
	* You can now opt to use the `json.Number` type instead of `float64` when parsing token JSON
* Added support for [Travis CI](https://travis-ci.org/dgrijalva/jwt-go)
* Fixed some bugs with ECDSA parsing

#### 2.3.0

* Added support for ECDSA signing methods
* Added support for RSA PSS signing methods (requires go v1.4)

#### 2.2.0

* Gracefully handle a `nil` `Keyfunc` being passed to `Parse`.  Result will now be the parsed token and an error, instead of a panic.

#### 2.1.0

Backwards compatible API change that was missed in 2.0.0.

* The `SignedString` method on `Token` now takes `interface{}` instead of `[]byte`

#### 2.0.0

There were two major reasons for breaking backwards compatibility with this update.  The first was a refactor required to expand the width of the RSA and HMAC-SHA signing implementations.  There will likely be no required code changes to support this change.

The second update, while unfortunately requiring a small change in integration, is required to open up this library to other signing methods.  Not all keys used for all signing methods have a single standard on-disk representation.  Requiring `[]byte` as the type for all keys proved too limiting.  Additionally, this implementation allows for pre-parsed tokens to be reused, which might matter in an application that parses a high volume of tokens with a small set of keys.  Backwards compatibilty has been maintained for passing `[]byte` to the RSA signing methods, but they will also accept `*rsa.PublicKey` and `*rsa.PrivateKey`.

It is likely the only integration change required here will be to change `func(t *jwt.Token) ([]byte, error)` to `func(t *jwt.Token) (interface{}, error)` when calling `Parse`.

* **Compatibility Breaking Changes**
	* `SigningMethodHS256` is now `*SigningMethodHMAC` instead of `type struct`
	* `SigningMethodRS256` is now `*SigningMethodRSA` instead of `type struct`
	* `KeyFunc` now returns `interface{}` instead of `[]byte`
	* `SigningMethod.Sign` now takes `interface{}` instead of `[]byte` for the key
	* `SigningMethod.Verify` now takes `interface{}` instead of `[]byte` for the key
* Renamed type `SigningMethodHS256` to `SigningMethodHMAC`.  Specific sizes are now just instances of this type.
    * Added public package global `SigningMethodHS256`
    * Added public package global `SigningMethodHS384`
    * Added public package global `SigningMethodHS512`
* Renamed type `SigningMethodRS256` to `SigningMethodRSA`.  Specific sizes are now just instances of this type.
    * Added public package global `SigningMethodRS256`
    * Added public package global `SigningMethodRS384`
    * Added public package global `SigningMethodRS512`
* Moved sample private key for HMAC tests from an inline value to a file on disk.  Value is unchanged.
* Refactored the RSA implementation to be easier to read
* Exposed helper methods `ParseRSAPrivateKeyFromPEM` and `ParseRSAPublicKeyFromPEM`

#### 1.0.2

* Fixed bug in parsing public keys from certificates
* Added more tests around the parsing of keys for RS256
* Code refactoring in RS256 implementation.  No functional changes

#### 1.0.1

* Fixed panic if RS256 signing method was passed an invalid key

#### 1.0.0

* First versioned release
* API stabilized
* Supports creating, signing, parsing, and validating JWT tokens
* Supports RS256 and HS256 signing methodsA [go](http://www.golang.org) (or 'golang' for search engine friendliness) implementation of [JSON Web Tokens](http://self-issued.info/docs/draft-ietf-oauth-json-web-token.html)

[![Build Status](https://travis-ci.org/dgrijalva/jwt-go.svg?branch=master)](https://travis-ci.org/dgrijalva/jwt-go)

**BREAKING CHANGES:*** Version 3.0.0 is here. It includes _a lot_ of changes including a few that break the API.  We've tried to break as few things as possible, so there should just be a few type signature changes.  A full list of breaking changes is available in `VERSION_HISTORY.md`.  See `MIGRATION_GUIDE.md` for more information on updating your code.

**NOTICE:** A vulnerability in JWT was [recently published](https://auth0.com/blog/2015/03/31/critical-vulnerabilities-in-json-web-token-libraries/).  As this library doesn't force users to validate the `alg` is what they expected, it's possible your usage is effected.  There will be an update soon to remedy this, and it will likey require backwards-incompatible changes to the API.  In the short term, please make sure your implementation verifies the `alg` is what you expect.


## What the heck is a JWT?

JWT.io has [a great introduction](https://jwt.io/introduction) to JSON Web Tokens.

In short, it's a signed JSON object that does something useful (for example, authentication).  It's commonly used for `Bearer` tokens in Oauth 2.  A token is made of three parts, separated by `.`'s.  The first two parts are JSON objects, that have been [base64url](http://tools.ietf.org/html/rfc4648) encoded.  The last part is the signature, encoded the same way.

The first part is called the header.  It contains the necessary information for verifying the last part, the signature.  For example, which encryption method was used for signing and what key was used.

The part in the middle is the interesting bit.  It's called the Claims and contains the actual stuff you care about.  Refer to [the RFC](http://self-issued.info/docs/draft-jones-json-web-token.html) for information about reserved keys and the proper way to add your own.

## What's in the box?

This library supports the parsing and verification as well as the generation and signing of JWTs.  Current supported signing algorithms are HMAC SHA, RSA, RSA-PSS, and ECDSA, though hooks are present for adding your own.

## Examples

See [the project documentation](https://godoc.org/github.com/dgrijalva/jwt-go) for examples of usage:

* [Simple example of parsing and validating a token](https://godoc.org/github.com/dgrijalva/jwt-go#example_Parse_hmac)
* [Simple example of building and signing a token](https://godoc.org/github.com/dgrijalva/jwt-go#example_New_hmac)
* [Directory of Examples](https://godoc.org/github.com/dgrijalva/jwt-go#pkg-examples)

## Extensions

This library publishes all the necessary components for adding your own signing methods.  Simply implement the `SigningMethod` interface and register a factory method using `RegisterSigningMethod`.  

Here's an example of an extension that integrates with the Google App Engine signing tools: https://github.com/someone1/gcp-jwt-go

## Compliance

This library was last reviewed to comply with [RTF 7519](http://www.rfc-editor.org/info/rfc7519) dated May 2015 with a few notable differences: 

* In order to protect against accidental use of [Unsecured JWTs](http://self-issued.info/docs/draft-ietf-oauth-json-web-token.html#UnsecuredJWT), tokens using `alg=none` will only be accepted if the constant `jwt.UnsafeAllowNoneSignatureType` is provided as the key.

## Project Status & Versioning

This library is considered production ready.  Feedback and feature requests are appreciated.  The API should be considered stable.  There should be very few backwards-incompatible changes outside of major version updates (and only with good reason).

This project uses [Semantic Versioning 2.0.0](http://semver.org).  Accepted pull requests will land on `master`.  Periodically, versions will be tagged from `master`.  You can find all the releases on [the project releases page](https://github.com/dgrijalva/jwt-go/releases).

While we try to make it obvious when we make breaking changes, there isn't a great mechanism for pushing announcements out to users.  You may want to use this alternative package include: `gopkg.in/dgrijalva/jwt-go.v2`.  It will do the right thing WRT semantic versioning.

## Usage Tips

### Signing vs Encryption

A token is simply a JSON object that is signed by its author. this tells you exactly two things about the data:

* The author of the token was in the possession of the signing secret
* The data has not been modified since it was signed

It's important to know that JWT does not provide encryption, which means anyone who has access to the token can read its contents. If you need to protect (encrypt) the data, there is a companion spec, `JWE`, that provides this functionality. JWE is currently outside the scope of this library.

### Choosing a Signing Method

There are several signing methods available, and you should probably take the time to learn about the various options before choosing one.  The principal design decision is most likely going to be symmetric vs asymmetric.

Symmetric signing methods, such as HSA, use only a single secret. This is probably the simplest signing method to use since any `[]byte` can be used as a valid secret. They are also slightly computationally faster to use, though this rarely is enough to matter. Symmetric signing methods work the best when both producers and consumers of tokens are trusted, or even the same system. Since the same secret is used to both sign and validate tokens, you can't easily distribute the key for validation.

Asymmetric signing methods, such as RSA, use different keys for signing and verifying tokens. This makes it possible to produce tokens with a private key, and allow any consumer to access the public key for verification.

### JWT and OAuth

It's worth mentioning that OAuth and JWT are not the same thing. A JWT token is simply a signed JSON object. It can be used anywhere such a thing is useful. There is some confusion, though, as JWT is the most common type of bearer token used in OAuth2 authentication.

Without going too far down the rabbit hole, here's a description of the interaction of these technologies:

* OAuth is a protocol for allowing an identity provider to be separate from the service a user is logging in to.  For example, whenever you use Facebook to log into a different service (Yelp, Spotify, etc), you are using OAuth.
* OAuth defines several options for passing around authentication data. One popular method is called a "bearer token". A bearer token is simply a string that _should_ only be held by an authenticated user. Thus, simply presenting this token proves your identity. You can probably derive from here why a JWT might make a good bearer token.
* Because bearer tokens are used for authentication, it's important they're kept secret. This is why transactions that use bearer tokens typically happen over SSL.
 
## More

Documentation can be found [on godoc.org](http://godoc.org/github.com/dgrijalva/jwt-go).

The command line utility included in this project (cmd/jwt) provides a straightforward example of token creation and parsing as well as a useful tool for debugging your own integration.  You'll also find several implementation examples in to documentation.
## Migration Guide from v2 -> v3

Version 3 adds several new, frequently requested features.  To do so, it introduces a few breaking changes.  We've worked to keep these as minimal as possible.  This guide explains the breaking changes and how you can quickly update your code.

### `Token.Claims` is now an interface type

The most requested feature from the 2.0 verison of this library was the ability to provide a custom type to the JSON parser for claims. This was implemented by introducing a new interface, `Claims`, to replace `map[string]interface{}`.  We also included two concrete implementations of `Claims`: `MapClaims` and `StandardClaims`.

`MapClaims` is an alias for `map[string]interface{}` with built in validation behavior.  It is the default claims type when using `Parse`.  The usage is unchanged except you must type cast the claims property.

The old example for parsing a token looked like this..

```go
	if token, err := jwt.Parse(tokenString, keyLookupFunc); err == nil {
		fmt.Printf("Token for user %v expires %v", token.Claims["user"], token.Claims["exp"])
	}
```

is now directly mapped to...

```go
	if token, err := jwt.Parse(tokenString, keyLookupFunc); err == nil {
		claims := token.Claims.(jwt.MapClaims)
		fmt.Printf("Token for user %v expires %v", claims["user"], claims["exp"])
	}
```

`StandardClaims` is designed to be embedded in your custom type.  You can supply a custom claims type with the new `ParseWithClaims` function.  Here's an example of using a custom claims type.

```go
	type MyCustomClaims struct {
		User string
		*StandardClaims
	}
	
	if token, err := jwt.ParseWithClaims(tokenString, &MyCustomClaims{}, keyLookupFunc); err == nil {
		claims := token.Claims.(*MyCustomClaims)
		fmt.Printf("Token for user %v expires %v", claims.User, claims.StandardClaims.ExpiresAt)
	}
```

### `ParseFromRequest` has been moved

To keep this library focused on the tokens without becoming overburdened with complex request processing logic, `ParseFromRequest` and its new companion `ParseFromRequestWithClaims` have been moved to a subpackage, `request`.  The method signatues have also been augmented to receive a new argument: `Extractor`.

`Extractors` do the work of picking the token string out of a request.  The interface is simple and composable.

This simple parsing example:

```go
	if token, err := jwt.ParseFromRequest(tokenString, req, keyLookupFunc); err == nil {
		fmt.Printf("Token for user %v expires %v", token.Claims["user"], token.Claims["exp"])
	}
```

is directly mapped to:

```go
	if token, err := request.ParseFromRequest(tokenString, request.OAuth2Extractor, req, keyLookupFunc); err == nil {
		fmt.Printf("Token for user %v expires %v", token.Claims["user"], token.Claims["exp"])
	}
```

There are several concrete `Extractor` types provided for your convenience:

* `HeaderExtractor` will search a list of headers until one contains content.
* `ArgumentExtractor` will search a list of keys in request query and form arguments until one contains content.
* `MultiExtractor` will try a list of `Extractors` in order until one returns content.
* `AuthorizationHeaderExtractor` will look in the `Authorization` header for a `Bearer` token.
* `OAuth2Extractor` searches the places an OAuth2 token would be specified (per the spec): `Authorization` header and `access_token` argument
* `PostExtractionFilter` wraps an `Extractor`, allowing you to process the content before it's parsed.  A simple example is stripping the `Bearer ` text from a header


### RSA signing methods no longer accept `[]byte` keys

Due to a [critical vulnerability](https://auth0.com/blog/2015/03/31/critical-vulnerabilities-in-json-web-token-libraries/), we've decided the convenience of accepting `[]byte` instead of `rsa.PublicKey` or `rsa.PrivateKey` isn't worth the risk of misuse.

To replace this behavior, we've added two helper methods: `ParseRSAPrivateKeyFromPEM(key []byte) (*rsa.PrivateKey, error)` and `ParseRSAPublicKeyFromPEM(key []byte) (*rsa.PublicKey, error)`.  These are just simple helpers for unpacking PEM encoded PKCS1 and PKCS8 keys. If your keys are encoded any other way, all you need to do is convert them to the `crypto/rsa` package's types.

```go 
	func keyLookupFunc(*Token) (interface{}, error) {
		// Don't forget to validate the alg is what you expect:
		if _, ok := token.Method.(*jwt.SigningMethodRSA); !ok {
			return nil, fmt.Errorf("Unexpected signing method: %v", token.Header["alg"])
		}
		
		// Look up key 
		key, err := lookupPublicKey(token.Header["kid"])
		if err != nil {
			return nil, err
		}
		
		// Unpack key from PEM encoded PKCS8
		return jwt.ParseRSAPublicKeyFromPEM(key)
	}
```
# Change Log

**ATTN**: This project uses [semantic versioning](http://semver.org/).

## [Unreleased] -

## [1.0.0] - 2018-09-01

### Fixed
- `Logger` middleware now correctly handles paths containing a `%` instead of trying to treat it as a format specifier

## [0.3.0] - 2017-11-11
### Added
- `With()` helper for building a new `Negroni` struct chaining handlers from
  existing `Negroni` structs
- Format log output in `Logger` middleware via a configurable `text/template`
  string injectable via `.SetFormat`. Added `LoggerDefaultFormat` and
  `LoggerDefaultDateFormat` to configure the default template and date format
  used by the `Logger` middleware.
- Support for HTTP/2 pusher support via `http.Pusher` interface for Go 1.8+.
- `WrapFunc` to convert `http.HandlerFunc` into a `negroni.Handler`
- `Formatter` field added to `Recovery` middleware to allow configuring how
  `panic`s are output. Default of `TextFormatter` (how it was output in
  `0.2.0`) used. `HTMLPanicFormatter` also added to allow easy outputing of
  `panic`s as HTML.

### Fixed
- `Written()` correct returns `false` if no response header has been written
- Only implement `http.CloseNotifier` with the `negroni.ResponseWriter` if the
  underlying `http.ResponseWriter` implements it (previously would always
  implement it and panic if the underlying `http.ResponseWriter` did not.

### Changed
- Set default status to `0` in the case that no handler writes status -- was
  previously `200` (in 0.2.0, before that it was `0` so this reestablishes that
  behavior)
- Catch `panic`s thrown by callbacks provided to the `Recovery` handler
- Recovery middleware will set `text/plain` content-type if none is set
- `ALogger` interface to allow custom logger outputs to be used with the
  `Logger` middleware. Changes embeded field in `negroni.Logger` from `Logger`
  to `ALogger`.
- Default `Logger` middleware output changed to be more structure and verbose
  (also now configurable, see `Added`)
- Automatically bind to port specified in `$PORT` in `.Run()` if an address is
  not passed in. Fall back to binding to `:8080` if no address specified
  (configuable via `DefaultAddress`).
- `PanicHandlerFunc` added to `Recovery` middleware to enhance custom handling
  of `panic`s by providing additional information to the handler including the
  stack and the `http.Request`. `Recovery.ErrorHandlerFunc` was also added, but
  deprecated in favor of the new `PanicHandlerFunc`.

## [0.2.0] - 2016-05-10
### Added
- Support for variadic handlers in `New()`
- Added `Negroni.Handlers()` to fetch all of the handlers for a given chain
- Allowed size in `Recovery` handler was bumped to 8k
- `Negroni.UseFunc` to push another handler onto the chain

### Changed
- Set the status before calling `beforeFuncs` so the information is available to them
- Set default status to `200` in the case that no handler writes status -- was previously `0`
- Panic if `nil` handler is given to `negroni.Use`

## 0.1.0 - 2013-07-22
### Added
- Initial implementation.

[Unreleased]: https://github.com/urfave/negroni/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/urfave/negroni/compare/v0.1.0...v0.2.0
# Negroni
[![GoDoc](https://godoc.org/github.com/urfave/negroni?status.svg)](http://godoc.org/github.com/urfave/negroni)
[![Build Status](https://travis-ci.org/urfave/negroni.svg?branch=master)](https://travis-ci.org/urfave/negroni)
[![codebeat](https://codebeat.co/badges/47d320b1-209e-45e8-bd99-9094bc5111e2)](https://codebeat.co/projects/github-com-urfave-negroni)
[![codecov](https://codecov.io/gh/urfave/negroni/branch/master/graph/badge.svg)](https://codecov.io/gh/urfave/negroni)

**Notice:** This is the library formerly known as
`github.com/codegangsta/negroni` -- Github will automatically redirect requests
to this repository, but we recommend updating your references for clarity.

Negroni is an idiomatic approach to web middleware in Go. It is tiny,
non-intrusive, and encourages use of `net/http` Handlers.

If you like the idea of [Martini](https://github.com/go-martini/martini), but
you think it contains too much magic, then Negroni is a great fit.

Language Translations:
* [Deutsch (de_DE)](translations/README_de_de.md)
* [Português Brasileiro (pt_BR)](translations/README_pt_br.md)
* [简体中文 (zh_CN)](translations/README_zh_CN.md)
* [繁體中文 (zh_TW)](translations/README_zh_tw.md)
* [日本語 (ja_JP)](translations/README_ja_JP.md)
* [Français (fr_FR)](translations/README_fr_FR.md)

## Getting Started

After installing Go and setting up your
[GOPATH](http://golang.org/doc/code.html#GOPATH), create your first `.go` file.
We'll call it `server.go`.

<!-- { "interrupt": true } -->
``` go
package main

import (
  "fmt"
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    fmt.Fprintf(w, "Welcome to the home page!")
  })

  n := negroni.Classic() // Includes some default middlewares
  n.UseHandler(mux)

  http.ListenAndServe(":3000", n)
}
```

Then install the Negroni package (**NOTE**: &gt;= **go 1.1** is required):

```
go get github.com/urfave/negroni
```

Then run your server:

```
go run server.go
```

You will now have a Go `net/http` webserver running on `localhost:3000`.

### Packaging

If you are on Debian, `negroni` is also available as [a
package](https://packages.debian.org/sid/golang-github-urfave-negroni-dev) that
you can install via `apt install golang-github-urfave-negroni-dev` (at the time
of writing, it is in the `sid` repositories).

## Is Negroni a Framework?

Negroni is **not** a framework. It is a middleware-focused library that is
designed to work directly with `net/http`.

## Routing?

Negroni is BYOR (Bring your own Router). The Go community already has a number
of great http routers available, and Negroni tries to play well with all of them
by fully supporting `net/http`. For instance, integrating with [Gorilla Mux]
looks like so:

``` go
router := mux.NewRouter()
router.HandleFunc("/", HomeHandler)

n := negroni.New(Middleware1, Middleware2)
// Or use a middleware with the Use() function
n.Use(Middleware3)
// router goes last
n.UseHandler(router)

http.ListenAndServe(":3001", n)
```

## `negroni.Classic()`

`negroni.Classic()` provides some default middleware that is useful for most
applications:

* [`negroni.Recovery`](#recovery) - Panic Recovery Middleware.
* [`negroni.Logger`](#logger) - Request/Response Logger Middleware.
* [`negroni.Static`](#static) - Static File serving under the "public"
  directory.

This makes it really easy to get started with some useful features from Negroni.

## Handlers

Negroni provides a bidirectional middleware flow. This is done through the
`negroni.Handler` interface:

``` go
type Handler interface {
  ServeHTTP(rw http.ResponseWriter, r *http.Request, next http.HandlerFunc)
}
```

If a middleware hasn't already written to the `ResponseWriter`, it should call
the next `http.HandlerFunc` in the chain to yield to the next middleware
handler.  This can be used for great good:

``` go
func MyMiddleware(rw http.ResponseWriter, r *http.Request, next http.HandlerFunc) {
  // do some stuff before
  next(rw, r)
  // do some stuff after
}
```

And you can map it to the handler chain with the `Use` function:

``` go
n := negroni.New()
n.Use(negroni.HandlerFunc(MyMiddleware))
```

You can also map plain old `http.Handler`s:

``` go
n := negroni.New()

mux := http.NewServeMux()
// map your routes

n.UseHandler(mux)

http.ListenAndServe(":3000", n)
```

## `With()`

Negroni has a convenience function called `With`. `With` takes one or more
`Handler` instances and returns a new `Negroni` with the combination of the
receiver's handlers and the new handlers.

```go
// middleware we want to reuse
common := negroni.New()
common.Use(MyMiddleware1)
common.Use(MyMiddleware2)

// `specific` is a new negroni with the handlers from `common` combined with the
// the handlers passed in
specific := common.With(
	SpecificMiddleware1,
	SpecificMiddleware2
)
```

## `Run()`

Negroni has a convenience function called `Run`. `Run` takes an addr string
identical to [`http.ListenAndServe`](https://godoc.org/net/http#ListenAndServe).

<!-- { "interrupt": true } -->
``` go
package main

import (
  "github.com/urfave/negroni"
)

func main() {
  n := negroni.Classic()
  n.Run(":8080")
}
```
If no address is provided, the `PORT` environment variable is used instead.
If the `PORT` environment variable is not defined, the default address will be used. 
See [Run](https://godoc.org/github.com/urfave/negroni#Negroni.Run) for a complete description.

In general, you will want to use `net/http` methods and pass `negroni` as a
`Handler`, as this is more flexible, e.g.:

<!-- { "interrupt": true } -->
``` go
package main

import (
  "fmt"
  "log"
  "net/http"
  "time"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    fmt.Fprintf(w, "Welcome to the home page!")
  })

  n := negroni.Classic() // Includes some default middlewares
  n.UseHandler(mux)

  s := &http.Server{
    Addr:           ":8080",
    Handler:        n,
    ReadTimeout:    10 * time.Second,
    WriteTimeout:   10 * time.Second,
    MaxHeaderBytes: 1 << 20,
  }
  log.Fatal(s.ListenAndServe())
}
```

## Route Specific Middleware

If you have a route group of routes that need specific middleware to be
executed, you can simply create a new Negroni instance and use it as your route
handler.

``` go
router := mux.NewRouter()
adminRoutes := mux.NewRouter()
// add admin routes here

// Create a new negroni for the admin middleware
router.PathPrefix("/admin").Handler(negroni.New(
  Middleware1,
  Middleware2,
  negroni.Wrap(adminRoutes),
))
```

If you are using [Gorilla Mux], here is an example using a subrouter:

``` go
router := mux.NewRouter()
subRouter := mux.NewRouter().PathPrefix("/subpath").Subrouter().StrictSlash(true)
subRouter.HandleFunc("/", someSubpathHandler) // "/subpath/"
subRouter.HandleFunc("/:id", someSubpathHandler) // "/subpath/:id"

// "/subpath" is necessary to ensure the subRouter and main router linkup
router.PathPrefix("/subpath").Handler(negroni.New(
  Middleware1,
  Middleware2,
  negroni.Wrap(subRouter),
))
```

`With()` can be used to eliminate redundancy for middlewares shared across
routes.

``` go
router := mux.NewRouter()
apiRoutes := mux.NewRouter()
// add api routes here
webRoutes := mux.NewRouter()
// add web routes here

// create common middleware to be shared across routes
common := negroni.New(
	Middleware1,
	Middleware2,
)

// create a new negroni for the api middleware
// using the common middleware as a base
router.PathPrefix("/api").Handler(common.With(
  APIMiddleware1,
  negroni.Wrap(apiRoutes),
))
// create a new negroni for the web middleware
// using the common middleware as a base
router.PathPrefix("/web").Handler(common.With(
  WebMiddleware1,
  negroni.Wrap(webRoutes),
))
```

## Bundled Middleware

### Static

This middleware will serve files on the filesystem. If the files do not exist,
it proxies the request to the next middleware. If you want the requests for
non-existent files to return a `404 File Not Found` to the user you should look
at using [http.FileServer](https://golang.org/pkg/net/http/#FileServer) as
a handler.

Example:

<!-- { "interrupt": true } -->
``` go
package main

import (
  "fmt"
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    fmt.Fprintf(w, "Welcome to the home page!")
  })

  // Example of using a http.FileServer if you want "server-like" rather than "middleware" behavior
  // mux.Handle("/public", http.FileServer(http.Dir("/home/public")))

  n := negroni.New()
  n.Use(negroni.NewStatic(http.Dir("/tmp")))
  n.UseHandler(mux)

  http.ListenAndServe(":3002", n)
}
```

Will serve files from the `/tmp` directory first, but proxy calls to the next
handler if the request does not match a file on the filesystem.

### Recovery

This middleware catches `panic`s and responds with a `500` response code. If
any other middleware has written a response code or body, this middleware will
fail to properly send a 500 to the client, as the client has already received
the HTTP response code. Additionally, an `PanicHandlerFunc` can be attached
to report 500's to an error reporting service such as Sentry or Airbrake.

Example:

<!-- { "interrupt": true } -->
``` go
package main

import (
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    panic("oh no")
  })

  n := negroni.New()
  n.Use(negroni.NewRecovery())
  n.UseHandler(mux)

  http.ListenAndServe(":3003", n)
}
```

Will return a `500 Internal Server Error` to each request. It will also log the
stack traces as well as print the stack trace to the requester if `PrintStack`
is set to `true` (the default).

Example with error handler:

``` go
package main

import (
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    panic("oh no")
  })

  n := negroni.New()
  recovery := negroni.NewRecovery()
  recovery.PanicHandlerFunc = reportToSentry
  n.Use(recovery)
  n.UseHandler(mux)

  http.ListenAndServe(":3003", n)
}

func reportToSentry(info *negroni.PanicInformation) {
    // write code here to report error to Sentry
}
```

The middleware simply output the informations on STDOUT by default.
You can customize the output process by using the `SetFormatter()` function.

You can use also the `HTMLPanicFormatter` to display a pretty HTML when a crash occurs.

<!-- { "interrupt": true } -->
``` go
package main

import (
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    panic("oh no")
  })

  n := negroni.New()
  recovery := negroni.NewRecovery()
  recovery.Formatter = &negroni.HTMLPanicFormatter{}
  n.Use(recovery)
  n.UseHandler(mux)

  http.ListenAndServe(":3003", n)
}
```

## Logger

This middleware logs each incoming request and response.

Example:

<!-- { "interrupt": true } -->
``` go
package main

import (
  "fmt"
  "net/http"

  "github.com/urfave/negroni"
)

func main() {
  mux := http.NewServeMux()
  mux.HandleFunc("/", func(w http.ResponseWriter, req *http.Request) {
    fmt.Fprintf(w, "Welcome to the home page!")
  })

  n := negroni.New()
  n.Use(negroni.NewLogger())
  n.UseHandler(mux)

  http.ListenAndServe(":3004", n)
}
```

Will print a log similar to:

```
[negroni] 2017-10-04T14:56:25+02:00 | 200 |      378µs | localhost:3004 | GET /
```

on each request.

You can also set your own log format by calling the `SetFormat` function. The format is a template string with fields as mentioned in the `LoggerEntry` struct. So, as an example -

```go
l.SetFormat("[{{.Status}} {{.Duration}}] - {{.Request.UserAgent}}")
```

will show something like - `[200 18.263µs] - Go-User-Agent/1.1 `

## Third Party Middleware

Here is a current list of Negroni compatible middlware. Feel free to put up a PR
linking your middleware if you have built one:

| Middleware | Author | Description |
| -----------|--------|-------------|
| [authz](https://github.com/casbin/negroni-authz) | [Yang Luo](https://github.com/hsluoyz) | ACL, RBAC, ABAC Authorization middlware based on [Casbin](https://github.com/casbin/casbin) |
| [binding](https://github.com/mholt/binding) | [Matt Holt](https://github.com/mholt) | Data binding from HTTP requests into structs |
| [cloudwatch](https://github.com/cvillecsteele/negroni-cloudwatch) | [Colin Steele](https://github.com/cvillecsteele) | AWS cloudwatch metrics middleware |
| [cors](https://github.com/rs/cors) | [Olivier Poitrey](https://github.com/rs) | [Cross Origin Resource Sharing](http://www.w3.org/TR/cors/) (CORS) support |
| [csp](https://github.com/awakenetworks/csp) | [Awake Networks](https://github.com/awakenetworks) | [Content Security Policy](https://www.w3.org/TR/CSP2/) (CSP) support |
| [delay](https://github.com/jeffbmartinez/delay) | [Jeff Martinez](https://github.com/jeffbmartinez) | Add delays/latency to endpoints. Useful when testing effects of high latency |
| [New Relic Go Agent](https://github.com/yadvendar/negroni-newrelic-go-agent) | [Yadvendar Champawat](https://github.com/yadvendar) | Official [New Relic Go Agent](https://github.com/newrelic/go-agent) (currently in beta)  |
| [gorelic](https://github.com/jingweno/negroni-gorelic) | [Jingwen Owen Ou](https://github.com/jingweno) | New Relic agent for Go runtime |
| [Graceful](https://github.com/tylerb/graceful) | [Tyler Bunnell](https://github.com/tylerb) | Graceful HTTP Shutdown |
| [gzip](https://github.com/phyber/negroni-gzip) | [phyber](https://github.com/phyber) | GZIP response compression |
| [JWT Middleware](https://github.com/auth0/go-jwt-middleware) | [Auth0](https://github.com/auth0) | Middleware checks for a JWT on the `Authorization` header on incoming requests and decodes it|
| [JWT Middleware](https://github.com/mfuentesg/go-jwtmiddleware) | [Marcelo Fuentes](https://github.com/mfuentesg) | JWT middleware for golang |
| [logrus](https://github.com/meatballhat/negroni-logrus) | [Dan Buch](https://github.com/meatballhat) | Logrus-based logger |
| [oauth2](https://github.com/goincremental/negroni-oauth2) | [David Bochenski](https://github.com/bochenski) | oAuth2 middleware |
| [onthefly](https://github.com/xyproto/onthefly) | [Alexander Rødseth](https://github.com/xyproto) | Generate TinySVG, HTML and CSS on the fly |
| [permissions2](https://github.com/xyproto/permissions2) | [Alexander Rødseth](https://github.com/xyproto) | Cookies, users and permissions |
| [prometheus](https://github.com/zbindenren/negroni-prometheus) | [Rene Zbinden](https://github.com/zbindenren) | Easily create metrics endpoint for the [prometheus](http://prometheus.io) instrumentation tool |
| [render](https://github.com/unrolled/render) | [Cory Jacobsen](https://github.com/unrolled) | Render JSON, XML and HTML templates |
| [RestGate](https://github.com/pjebs/restgate) | [Prasanga Siripala](https://github.com/pjebs) | Secure authentication for REST API endpoints |
| [secure](https://github.com/unrolled/secure) | [Cory Jacobsen](https://github.com/unrolled) | Middleware that implements a few quick security wins |
| [sessions](https://github.com/goincremental/negroni-sessions) | [David Bochenski](https://github.com/bochenski) | Session Management |
| [stats](https://github.com/thoas/stats) | [Florent Messa](https://github.com/thoas) | Store information about your web application (response time, etc.) |
| [VanGoH](https://github.com/auroratechnologies/vangoh) | [Taylor Wrobel](https://github.com/twrobel3) | Configurable [AWS-Style](http://docs.aws.amazon.com/AmazonS3/latest/dev/RESTAuthentication.html) HMAC authentication middleware |
| [xrequestid](https://github.com/pilu/xrequestid) | [Andrea Franz](https://github.com/pilu) | Middleware that assigns a random X-Request-Id header to each request |
| [mgo session](https://github.com/joeljames/nigroni-mgo-session) | [Joel James](https://github.com/joeljames) | Middleware that handles creating and closing mgo sessions per request |
| [digits](https://github.com/bamarni/digits) | [Bilal Amarni](https://github.com/bamarni) | Middleware that handles [Twitter Digits](https://get.digits.com/) authentication |
| [stats](https://github.com/guptachirag/stats) | [Chirag Gupta](https://github.com/guptachirag/stats) | Middleware that manages qps and latency stats for your endpoints and asynchronously flushes them to influx db |
| [Chaos](https://github.com/falzm/chaos) | [Marc Falzon](https://github.com/falzm) | Middleware for injecting chaotic behavior into application in a programmatic way |

## Examples

[Alexander Rødseth](https://github.com/xyproto) created
[mooseware](https://github.com/xyproto/mooseware), a skeleton for writing a
Negroni middleware handler.

[Prasanga Siripala](https://github.com/pjebs) created an effective skeleton structure for web-based Go/Negroni projects: [Go-Skeleton](https://github.com/pjebs/go-skeleton) 

## Live code reload?

[gin](https://github.com/codegangsta/gin) and
[fresh](https://github.com/pilu/fresh) both live reload negroni apps.

## Essential Reading for Beginners of Go & Negroni

* [Using a Context to pass information from middleware to end handler](http://elithrar.github.io/article/map-string-interface/)
* [Understanding middleware](https://mattstauffer.co/blog/laravel-5.0-middleware-filter-style)

## About

Negroni is obsessively designed by none other than the [Code
Gangsta](https://codegangsta.io/)

[Gorilla Mux]: https://github.com/gorilla/mux
[`http.FileSystem`]: https://godoc.org/net/http#FileSystem
# envconfig

[![Build Status](https://travis-ci.org/kelseyhightower/envconfig.svg)](https://travis-ci.org/kelseyhightower/envconfig)

```Go
import "github.com/kelseyhightower/envconfig"
```

## Documentation

See [godoc](http://godoc.org/github.com/kelseyhightower/envconfig)

## Usage

Set some environment variables:

```Bash
export MYAPP_DEBUG=false
export MYAPP_PORT=8080
export MYAPP_USER=Kelsey
export MYAPP_RATE="0.5"
export MYAPP_TIMEOUT="3m"
export MYAPP_USERS="rob,ken,robert"
export MYAPP_COLORCODES="red:1,green:2,blue:3"
```

Write some code:

```Go
package main

import (
    "fmt"
    "log"
    "time"

    "github.com/kelseyhightower/envconfig"
)

type Specification struct {
    Debug       bool
    Port        int
    User        string
    Users       []string
    Rate        float32
    Timeout     time.Duration
    ColorCodes  map[string]int
}

func main() {
    var s Specification
    err := envconfig.Process("myapp", &s)
    if err != nil {
        log.Fatal(err.Error())
    }
    format := "Debug: %v\nPort: %d\nUser: %s\nRate: %f\nTimeout: %s\n"
    _, err = fmt.Printf(format, s.Debug, s.Port, s.User, s.Rate, s.Timeout)
    if err != nil {
        log.Fatal(err.Error())
    }

    fmt.Println("Users:")
    for _, u := range s.Users {
        fmt.Printf("  %s\n", u)
    }

    fmt.Println("Color codes:")
    for k, v := range s.ColorCodes {
        fmt.Printf("  %s: %d\n", k, v)
    }
}
```

Results:

```Bash
Debug: false
Port: 8080
User: Kelsey
Rate: 0.500000
Timeout: 3m0s
Users:
  rob
  ken
  robert
Color codes:
  red: 1
  green: 2
  blue: 3
```

## Struct Tag Support

Envconfig supports the use of struct tags to specify alternate, default, and required
environment variables.

For example, consider the following struct:

```Go
type Specification struct {
    ManualOverride1 string `envconfig:"manual_override_1"`
    DefaultVar      string `default:"foobar"`
    RequiredVar     string `required:"true"`
    IgnoredVar      string `ignored:"true"`
    AutoSplitVar    string `split_words:"true"`
    RequiredAndAutoSplitVar    string `required:"true" split_words:"true"`
}
```

Envconfig has automatic support for CamelCased struct elements when the
`split_words:"true"` tag is supplied. Without this tag, `AutoSplitVar` above
would look for an environment variable called `MYAPP_AUTOSPLITVAR`. With the
setting applied it will look for `MYAPP_AUTO_SPLIT_VAR`. Note that numbers
will get globbed into the previous word. If the setting does not do the
right thing, you may use a manual override.

Envconfig will process value for `ManualOverride1` by populating it with the
value for `MYAPP_MANUAL_OVERRIDE_1`. Without this struct tag, it would have
instead looked up `MYAPP_MANUALOVERRIDE1`. With the `split_words:"true"` tag
it would have looked up `MYAPP_MANUAL_OVERRIDE1`.

```Bash
export MYAPP_MANUAL_OVERRIDE_1="this will be the value"

# export MYAPP_MANUALOVERRIDE1="and this will not"
```

If envconfig can't find an environment variable value for `MYAPP_DEFAULTVAR`,
it will populate it with "foobar" as a default value.

If envconfig can't find an environment variable value for `MYAPP_REQUIREDVAR`,
it will return an error when asked to process the struct.  If
`MYAPP_REQUIREDVAR` is present but empty, envconfig will not return an error.

If envconfig can't find an environment variable in the form `PREFIX_MYVAR`, and there
is a struct tag defined, it will try to populate your variable with an environment
variable that directly matches the envconfig tag in your struct definition:

```shell
export SERVICE_HOST=127.0.0.1
export MYAPP_DEBUG=true
```
```Go
type Specification struct {
    ServiceHost string `envconfig:"SERVICE_HOST"`
    Debug       bool
}
```

Envconfig won't process a field with the "ignored" tag set to "true", even if a corresponding
environment variable is set.

## Supported Struct Field Types

envconfig supports these struct field types:

  * string
  * int8, int16, int32, int64
  * bool
  * float32, float64
  * slices of any supported type
  * maps (keys and values of any supported type)
  * [encoding.TextUnmarshaler](https://golang.org/pkg/encoding/#TextUnmarshaler)
  * [encoding.BinaryUnmarshaler](https://golang.org/pkg/encoding/#BinaryUnmarshaler)
  * [time.Duration](https://golang.org/pkg/time/#Duration)

Embedded structs using these fields are also supported.

## Custom Decoders

Any field whose type (or pointer-to-type) implements `envconfig.Decoder` can
control its own deserialization:

```Bash
export DNS_SERVER=8.8.8.8
```

```Go
type IPDecoder net.IP

func (ipd *IPDecoder) Decode(value string) error {
    *ipd = IPDecoder(net.ParseIP(value))
    return nil
}

type DNSConfig struct {
    Address IPDecoder `envconfig:"DNS_SERVER"`
}
```

Also, envconfig will use a `Set(string) error` method like from the
[flag.Value](https://godoc.org/flag#Value) interface if implemented.
# Go CORS handler [![godoc](http://img.shields.io/badge/godoc-reference-blue.svg?style=flat)](https://godoc.org/github.com/rs/cors) [![license](http://img.shields.io/badge/license-MIT-red.svg?style=flat)](https://raw.githubusercontent.com/rs/cors/master/LICENSE) [![build](https://img.shields.io/travis/rs/cors.svg?style=flat)](https://travis-ci.org/rs/cors) [![Coverage](http://gocover.io/_badge/github.com/rs/cors)](http://gocover.io/github.com/rs/cors)

CORS is a `net/http` handler implementing [Cross Origin Resource Sharing W3 specification](http://www.w3.org/TR/cors/) in Golang.

## Getting Started

After installing Go and setting up your [GOPATH](http://golang.org/doc/code.html#GOPATH), create your first `.go` file. We'll call it `server.go`.

```go
package main

import (
    "net/http"

    "github.com/rs/cors"
)

func main() {
    mux := http.NewServeMux()
    mux.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
        w.Header().Set("Content-Type", "application/json")
        w.Write([]byte("{\"hello\": \"world\"}"))
    })

    // cors.Default() setup the middleware with default options being
    // all origins accepted with simple methods (GET, POST). See
    // documentation below for more options.
    handler := cors.Default().Handler(mux)
    http.ListenAndServe(":8080", handler)
}
```

Install `cors`:

    go get github.com/rs/cors

Then run your server:

    go run server.go

The server now runs on `localhost:8080`:

    $ curl -D - -H 'Origin: http://foo.com' http://localhost:8080/
    HTTP/1.1 200 OK
    Access-Control-Allow-Origin: foo.com
    Content-Type: application/json
    Date: Sat, 25 Oct 2014 03:43:57 GMT
    Content-Length: 18

    {"hello": "world"}

### Allow * With Credentials Security Protection

This library has been modified to avoid a well known security issue when configured with `AllowedOrigins` to `*` and `AllowCredentials` to `true`. Such setup used to make the library reflects the request `Origin` header value, working around a security protection embedded into the standard that makes clients to refuse such configuration. This behavior has been removed with [#55](https://github.com/rs/cors/issues/55) and [#57](https://github.com/rs/cors/issues/57).

If you depend on this behavior and understand the implications, you can restore it using the `AllowOriginFunc` with `func(origin string) {return true}`.

Please refer to [#55](https://github.com/rs/cors/issues/55) for more information about the security implications.

### More Examples

* `net/http`: [examples/nethttp/server.go](https://github.com/rs/cors/blob/master/examples/nethttp/server.go)
* [Goji](https://goji.io): [examples/goji/server.go](https://github.com/rs/cors/blob/master/examples/goji/server.go)
* [Martini](http://martini.codegangsta.io): [examples/martini/server.go](https://github.com/rs/cors/blob/master/examples/martini/server.go)
* [Negroni](https://github.com/codegangsta/negroni): [examples/negroni/server.go](https://github.com/rs/cors/blob/master/examples/negroni/server.go)
* [Alice](https://github.com/justinas/alice): [examples/alice/server.go](https://github.com/rs/cors/blob/master/examples/alice/server.go)
* [HttpRouter](https://github.com/julienschmidt/httprouter): [examples/httprouter/server.go](https://github.com/rs/cors/blob/master/examples/httprouter/server.go)
* [Gorilla](http://www.gorillatoolkit.org/pkg/mux): [examples/gorilla/server.go](https://github.com/rs/cors/blob/master/examples/gorilla/server.go)
* [Buffalo](https://gobuffalo.io): [examples/buffalo/server.go](https://github.com/rs/cors/blob/master/examples/buffalo/server.go)
* [Gin](https://gin-gonic.github.io/gin): [examples/gin/server.go](https://github.com/rs/cors/blob/master/examples/gin/server.go)
* [Chi](https://github.com/go-chi/chi): [examples/chi/server.go](https://github.com/rs/cors/blob/master/examples/chi/server.go)

## Parameters

Parameters are passed to the middleware thru the `cors.New` method as follow:

```go
c := cors.New(cors.Options{
    AllowedOrigins: []string{"http://foo.com", "http://foo.com:8080"},
    AllowCredentials: true,
    // Enable Debugging for testing, consider disabling in production
    Debug: true,
})

// Insert the middleware
handler = c.Handler(handler)
```

* **AllowedOrigins** `[]string`: A list of origins a cross-domain request can be executed from. If the special `*` value is present in the list, all origins will be allowed. An origin may contain a wildcard (`*`) to replace 0 or more characters (i.e.: `http://*.domain.com`). Usage of wildcards implies a small performance penality. Only one wildcard can be used per origin. The default value is `*`.
* **AllowOriginFunc** `func (origin string) bool`: A custom function to validate the origin. It takes the origin as an argument and returns true if allowed, or false otherwise. If this option is set, the content of `AllowedOrigins` is ignored.
* **AllowOriginRequestFunc** `func (r *http.Request origin string) bool`: A custom function to validate the origin. It takes the HTTP Request object and the origin as argument and returns true if allowed or false otherwise. If this option is set, the content of `AllowedOrigins` and `AllowOriginFunc` is ignored
* **AllowedMethods** `[]string`: A list of methods the client is allowed to use with cross-domain requests. Default value is simple methods (`GET` and `POST`).
* **AllowedHeaders** `[]string`: A list of non simple headers the client is allowed to use with cross-domain requests.
* **ExposedHeaders** `[]string`: Indicates which headers are safe to expose to the API of a CORS API specification
* **AllowCredentials** `bool`: Indicates whether the request can include user credentials like cookies, HTTP authentication or client side SSL certificates. The default is `false`.
* **MaxAge** `int`: Indicates how long (in seconds) the results of a preflight request can be cached. The default is `0` which stands for no max age.
* **OptionsPassthrough** `bool`: Instructs preflight to let other potential next handlers to process the `OPTIONS` method. Turn this on if your application handles `OPTIONS`.
* **Debug** `bool`: Debugging flag adds additional output to debug server side CORS issues.

See [API documentation](http://godoc.org/github.com/rs/cors) for more info.

## Benchmarks

    BenchmarkWithout          20000000    64.6 ns/op      8 B/op    1 allocs/op
    BenchmarkDefault          3000000      469 ns/op    114 B/op    2 allocs/op
    BenchmarkAllowedOrigin    3000000      608 ns/op    114 B/op    2 allocs/op
    BenchmarkPreflight        20000000    73.2 ns/op      0 B/op    0 allocs/op
    BenchmarkPreflightHeader  20000000    73.6 ns/op      0 B/op    0 allocs/op
    BenchmarkParseHeaderList  2000000      847 ns/op    184 B/op    6 allocs/op
    BenchmarkParse…Single     5000000      290 ns/op     32 B/op    3 allocs/op
    BenchmarkParse…Normalized 2000000      776 ns/op    160 B/op    6 allocs/op

## Licenses

All source code is licensed under the [MIT License](https://raw.github.com/rs/cors/master/LICENSE).
