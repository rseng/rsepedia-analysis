## datapack: A Flexible Container to Transport and Manipulate Data and Associated Resources
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/datapack)](https://cran.r-project.org/package=datapack)
[![Build Status](https://travis-ci.org/ropensci/datapack.svg?branch=master)](https://travis-ci.org/ropensci/datapack)

- **Author**: Matthew B. Jones and Peter Slaughter ([NCEAS](https://www.nceas.ucsb.edu))
- [doi:10.5063/F1QV3JGM](https://doi.org/10.5063/F1QV3JGM)
- **License**: [Apache 2](https://opensource.org/licenses/Apache-2.0)
- [Package source code on Github](https://github.com/ropensci/datapack)
- [**Submit Bugs and feature requests**](https://github.com/ropensci/datapack/issues)

The *datapack* R package provides an abstraction for collating 
heterogeneous collections of data objects and metadata into a bundle that can 
be transported and loaded into a single composite file.  The methods in 
this package provide a convenient way to load data from common repositories 
such as DataONE into the R environment, and to document, serialize, and save 
data from R to data repositories worldwide.

> Note that this package ('datapack') is not related to the similarly named rOpenSci package 'DataPackageR'.
> Documentation from the DataPackageR github repository states that "DataPackageR is used to reproducibly
> process raw data into packaged, analysis-ready data sets."

## Installation Notes 

The *datapack* R package requires the R package *redland*. If you are installing on Ubuntu then the Redland C libraries
must be installed before the *redland* and *datapack* package can be installed. If you are installing on Mac OS X or Windows then installing these libraries is not required.

The following instructions illustrate how to install *datapack* and its requirements.

### Installing on Mac OS X

On Mac OS X datapack can be installed with the following commands:

```
install.packages("datapack")
library(datapack)
```

The *datapack* R package should be available for use at this point.

Note: if you wish to build the required *redland* package from source before installing *datapack*, please see the redland [installation instructions](https://github.com/ropensci/redland-bindings/tree/master/R/redland).

## Installing on Ubuntu

For Ubuntu, install the required Redland C libraries by entering the following commands 
in a terminal window:

```
sudo apt-get update
sudo apt-get install librdf0 librdf0-dev
```

Then install the R packages from the R console:

```
install.packages("datapack")
library(datapack)
```

The *datapack* R package should be available for use at this point

## Installing on Windows

For windows, the required redland R package is distributed as a binary release, so it is not
necessary to install any additional system libraries.

To install the R packages from the R console:

```
install.packages("datapack")
library(datapack)
```

## Quick Start

See the full manual for documentation, but once installed, the package can be run in R using:

```
library(datapack)
help("datapack")
```

Create a DataPackage and add metadata and data DataObjects to it:

```
library(datapack)
library(uuid)
dp <- new("DataPackage")
mdFile <- system.file("extdata/sample-eml.xml", package="datapack")
mdId <- paste("urn:uuid:", UUIDgenerate(), sep="")
md <- new("DataObject", id=mdId, format="eml://ecoinformatics.org/eml-2.1.0", file=mdFile)
addData(dp, md)

csvfile <- system.file("extdata/sample-data.csv", package="datapack")
sciId <- paste("urn:uuid:", UUIDgenerate(), sep="")
sciObj <- new("DataObject", id=sciId, format="text/csv", filename=csvfile)
dp <- addData(dp, sciObj)
ids <- getIdentifiers(dp)
```

Add a relationship to the DataPackage that shows that the metadata describes, or "documents", the science data:

```
dp <- insertRelationship(dp, subjectID=mdId, objectIDs=sciId)
relations <- getRelationships(dp)
```  

Create an Resource Description Framework representation of the relationships in the package:

```
serializationId <- paste("resourceMap", UUIDgenerate(), sep="")
filePath <- file.path(sprintf("%s/%s.rdf", tempdir(), serializationId))
status <- serializePackage(dp, filePath, id=serializationId, resolveURI="")
```  
Save the DataPackage to a file, using the BagIt packaging format:

```
bagitFile <- serializeToBagIt(dp) 
```

Note that the *dataone* R package can be used to upload a DataPackage to a DataONE Member Node
using the *uploadDataPackage* method. Please see the documentation for the *dataone* R package,
for example:

```
vignette("upload-data", package="dataone")
```

## Acknowledgements
Work on this package was supported by:

- NSF-ABI grant #1262458 to C. Gries, M. B. Jones, and S. Collins.
- NSF-DATANET grants #0830944 and #1430508 to W. Michener, M. B. Jones, D. Vieglais, S. Allard and P. Cruse
- NSF DIBBS grant #1443062 to T. Habermann and M. B. Jones
- NSF-PLR grant #1546024 to M. B. Jones, S. Baker-Yeboah, J. Dozier, M. Schildhauer, and A. Budden

Additional support was provided for working group collaboration by the National Center for Ecological Analysis and Synthesis, a Center funded by the University of California, Santa Barbara, and the State of California.

[![nceas_footer ](https://www.nceas.ucsb.edu/files/newLogo_0.png)](https://www.nceas.ucsb.edu)

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
## Test environments

* tested via rhub::check_for_cran()
  * Debian, R-devel, clang
  * Debian, R-devel, GCC
  * Debian, R-patched, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Fedora Linux, R-devel, GCC
  * macOS 10.13.6 High Sierra, R-release
  * macOS 10.13.6 High Sierra, R-release, CRAN's setup
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-release, 32/64 bit
* tested via devtools::check()
  * macOS 10.14.6, R 4.0.2
  * Ubuntu 17.10
* tested via win-builder
  * Windows: x86_64-w64-mingw32 (64-bit) R 4.0.2
  * Windows: x86_64-w64-mingw32 (64-bit) R 3.6.3
  * Windows: x86_64-w64-mingw32 (64-bit) R devel

## Changes since last release

* Use SHA-256 as the default hash algorithm (#117)
* Added 'checksumAlgorithm' argument to DataObject initialization method  (#117)
* Handle dc:creator in resource map properly #116
* Update tests for compatibility with testthat 3e (#87)
* Added 'targetPath' argument to DataObject to set 'prov:atLocation' for an object (#109)

## R CMD check results

* There were no ERRORs or WARNINGs or NOTES.
    
## Downstream dependencies

* revdepcheck:revdep_check() found no problems with the only downstream dependency: `dataone`.
---
title: "datapack R Package Overview"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
<!-- output: pdf_document -->
vignette: >
  %\VignetteIndexEntry{datapack R Package Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

## Overview

The *datapack* R package provides an abstraction for collating multiple data 
objects of different types and metadata describing those objects into a bundle that can be transported and loaded using a single composite file.  It is primarily meant as
a container to bundle together files for transport to or from DataONE data repositories.

The methods in this package provide a convenient way to load data from common repositories 
such as DataONE into the R environment, and to document, serialize, and save data from R to 
data repositories. 

## Create a Single Object 
The *datapack* DataObject class is a wrapper that contains both data and system metadata that describes the data. 
The data can be either R `raw` data or a data file, for example a CSV file. The system metadata includes attributes such as the 
object identifier, type, size, checksum, owner, version relationship to other objects, access rules, and other critical metadata. The DataObject class also holds additional metadata about the data file. For example, where the file should go when the package is downloaded. This is often times the same as the filepath however, care should be taken to not include any drive letters or unnecessary folders.

The following example shows how to create a DataObject locally from a CSV file:

```{r}
library(datapack)
library(uuid)
csvfile <- system.file("extdata/sample-data.csv", package="datapack")
myId <- paste("urn:uuid:", UUIDgenerate(), sep="")
myObj <- new("DataObject", id=myId, format="text/csv", filename=csvfile, targetPath="extdata/sample-data.csv")
```

The DataObject `myObj` now contains the CSV data as well as the system-level 
information describing the file, such as its identifier, type, and checksum. 
The *getData* method can be used to extract the data content of a DataObject. 
Using the example DataObject:

```{r}
rawData <- getData(myObj)
```

This raw data can be converted back to CSV format using the R commands:

```{r,eval=FALSE}
tf <- tempfile(fileext=".csv")
write.csv(rawToChar(rawData), tf, quote=F, row.names=F)
```

Alternatively, the CSV data could be converted into a data frame using standard R
functions:

```{r}
df <- read.csv(textConnection(rawToChar(rawData)))
head(df)
```

If the data were another format than CSV, such as PNG, JPEG, or NetCDF, the 
corresponding R packages could be used to handle the object.

Each DataObject has an identifier which can be used to refer to that object, and 
is meant to be globally unique so that it can be used in data repositories such 
as those from the DataONE federation. To retrieve the identifier associated with
a DataObject:

```{r}
getIdentifier(myObj)
```

In this case, the identifier was created in the UUID format, but other identifiers
such as DOIs (Digital Object Identifers) can also be used.
Each object also is associated with a specific format. To retrieve the format type:

```{r}
getFormatId(myObj)
```

### System Metadata

All system metadata information for a DataObject can be accessed directly from the SystemMetadata object contained in the DataObject. To access the `fileName` field, for example:

```{r}
myObj@sysmeta@fileName
```

The system metadata contains access policy information for the DataObject that could be used
by a data repository that the object is uploaded to. For example, when a DataObject is
uploaded to a [DataONE Member Node](https://www.dataone.org/current-member-nodes), 
the access policy is applied to the uploaded data
and controls access to the data on the Member Node by DataONE users.

#### Access Policy

Before the DataObject is uploaded, access can be set so that anyone can read the uploaded data:
```{r}
myObj <- setPublicAccess(myObj)
myObj@sysmeta@accessPolicy
```

Individual access rules can also be added one at a time.  The access rules are expressed
using the unique identifier for an individual, such as their ORCID identity, or whatever
form the repository supports.
```{r}
myObj <- addAccessRule(myObj, "https://orcid.org/0000-0003-0077-4738", "write")
myObj@sysmeta@accessPolicy
```

The permissions that can be set include:

- read: permission to read the data from the repository
- write: permission to perform update operations on the data
- changePermission: permission to control access to the data on the repository

Alternatively, multiple access rules can be added:
```{r}
accessRules <- data.frame(subject=c("uid=jsmith,o=Account,dc=example,dc=com",  
                                    "uid=jadams,o=Account,dc=example,dc=org"), 
                          permission=c("write", "changePermission"))
myObj <- addAccessRule(myObj, accessRules)
myObj@sysmeta@accessPolicy
```

The *dataone* R package can be used to upload or download DataObjects to a DataONE Member Node.
Please see the web page for the [*dataone*](https://github.com/DataONEorg/rdataone) R package and the
vignettes for more information:

```{r, eval=FALSE}
library(dataone)
vignette("download-data", package="dataone")
vignette("upload-data", package="dataone")
```

## Create a Collection of Objects

A DataPackage is a container for a set of DataObjects. DataObject is a class that is a proxy for data of any type, including traditional data like CSV, tabular data, and spatial rasters, but also for non-traditional objects like derived data, figures, and scripts in R and Python. A collection of related DataObjects can be placed in a DataPackage and actions can be performed on it, such as serializing the entire collection of objects into a package file, or uploading all package member objects to a data repository.

Figure 1. is a diagram of a typical DataPackage showing a metadata file that
describes, or `documents` the data granules that the package contains.

![](package-diagram.png)

This example creates a DataPackage with one DataObject containing metadata and two others containing science data. First the individual objects are created:

```{r}
metadataFile <- system.file("extdata/sample-eml.xml", package="datapack")
metadataId <- "metadataId"
metadataObj <- new("DataObject", id=metadataId, format="eml://ecoinformatics.org/eml-2.1.0", file=metadataFile)

csvfile <- system.file("extdata/sample-data.csv", package="datapack")
sciId <- "sciId1"
sciObj <- new("DataObject", id=sciId, format="text/csv", filename=csvfile)

outFile <- system.file("extdata/sample-data-filtered.csv", package="datapack")
sciId2 <- "sciId2"
sciObj2 <- new("DataObject", id=sciId2, filename=outFile, format="text/csv")
```

The identifier values used in this example are simple and easily recognizable for demonstration purposes. A more standard unique identifier can be created with the `uuid::UUIDgenerate()` function:

```{r}
myid <- paste("urn:uuid:", UUIDgenerate(), sep="")
myid
```

Next a DataPackage is created and the DataObjects are added to it. Note the `mo` argument in the `addMember` function when adding the data file. Including this argument specifies that the metadata object documents the data object. More information on relationships between DataObjects is included in the next section.
```{r}
dp <- new("DataPackage")
dp <- addMember(dp, do = metadataObj)
dp <- addMember(dp, do = sciObj, mo = metadataObj)
# The second object will be added in the next section 
```

Information can also be extracted from the DataPackage. To show the identifiers of the DataObjects that are in the package:

```{r}
getIdentifiers(dp)
```

To show the number of DataObjects in the package:
```{r}
getSize(dp)
```

To extract the data in a DataObject as raw data, ask for the data using the identifier of the DataObject:
```{r}
sciObjRaw <- getData(dp, sciId)
```

To get access to the full instance of the DataObject class representing a data object, 
use the `datapack::getMember` function and pass in the identifier of the desired object, 
which will return an instance of the DataObject class:
```{r}
mySciObj <- getMember(dp, sciId)
```

## Relationships Between DataObjects

The relationships between DataObjects in a DataPackage can be recorded in the DataPackage. 
For example, a typical relationship is that a DataObject containing a metadata 
document in a domain specific format such as Ecological Metadata Language (EML) 
or ISO19139 geospatial metadata can describe, or document, DataObjects containing 
associated science data. Adding relationship information about data package members 
may assist a consumer of the package in better understanding the contents of the 
package and how to make use of the package.

While the DataPackage can record any type of relationships that are important to
a community, we have provided functions to establish common relationships that 
are needed to understand scientific data in the DataONE federation.  These include 
the following typical provenance properties:

- `cito:documents`: for establishing that a metadata document provides descriptive
  information about one or more associated data objects
- `prov:wasDerivedFrom`: for asserting that a derived data object was created using 
  data from one or more source data objects
- `prov:used`: for asserting that when a program (such as an R script) was executed
  that it used one or more source data objects as inputs
- `prov:wasGeneratedBy`: for asserting that when a program (such as an R script) was executed
  that it generated one or more derived data objects as outputs
  
Figure 2. A DataPackage with provenance relationships.

![](package-diagram-with-prov.png)

### Linking a metadata file with one or more data files using cito:documents

As mentioned above, the fastest way to add the `cito:documents` relationship is to include the metadata object when a science data object is added to the package:

```{r}
dp <- addMember(dp, do = sciObj2, mo = metadataObj)
getRelationships(dp, condense=TRUE)
```

In that example, the `sciObj2` DataObject is added to the package using the `addMember` call, 
and the metadata object `metadataObj` is passed in to the function as well.  This
tells the DataPackage that `metadataId cito:documents sciId2`. The `cito:documents` relationship is defined by the [Citation Typing Ontology](https://sparontologies.github.io/cito/current/cito.html) (CITO)).

## Asserting data provenance relationships between objects

Relationships that describe the processing history of package members can be added. For example,
a program that performs a modeling calculation might read one or more source data files as inputs,
perform  a calculation based on the data read, and then write a data or graphics file 
characterizing the results of the model run. 

The `datapack` package uses the [ProvONE](https://purl.dataone.org/provone-v1-dev) data model to represent provenance relationships.

The following example demonstrates how to insert provenance relationships into a DataPackage
for the R program `logit-regression.R` that reads the source data file `binary.csv` and 
generates the derived image file `gre-predicted.png`. Using the example DataPackage for which 
DataObjects for the program input and output have already been added, we create a 
DataObject for the program, and call `describeWorkflow` to add the necessary provenance
relationships:

```{r}
dp <- new("DataPackage")

metadataFile <- system.file("extdata/sample-eml.xml", package="datapack")
metadataId <- "metadataId"
metadataObj <- new("DataObject", id=metadataId, format="eml://ecoinformatics.org/eml-2.1.0", file=metadataFile)

# This DataObject contains the program script that was executed
progObj <- new("DataObject", format="application/R", 
           filename=system.file("extdata/pkg-example/logit-regression-example.R", package="datapack"))
dp <- addMember(dp, progObj, mo = metadataObj)

doIn <- new("DataObject", format="text/csv", 
             filename=system.file("./extdata/pkg-example/binary.csv", package="datapack"))
dp <- addMember(dp, doIn, mo = metadataObj)

doOut <- new("DataObject", format="image/png", 
             filename=system.file("./extdata/pkg-example/gre-predicted.png", package="datapack"))
dp <- addMember(dp, doOut, mo = metadataObj)

# The arguments "sources" and "derivations" can also contain lists of "DataObjects"
dp <- describeWorkflow(dp, sources=doIn, program=progObj, derivations=doOut)


rels <- getRelationships(dp, condense=TRUE)
rels[grepl("prov:", rels$predicate),]
```



```{r, message = FALSE, warning = FALSE, fig.width=8, fig.height=8}
library(igraph)
plotRelationships(dp)
```


Note that in this example, the R script had previously been run and generated the image
file before `describeWorkflow()` was called.  The `sources` and `derivations` arguments 
for `describeWorkflow()` can be lists of either DataObjects or the identifiers of DataObjects.



### Inserting other (arbitrary) relationships

Other types of relationships between DataPackage member DataObjects can be recorded with the `insertRelationship` method. The main requirement is that each relationship to be described
needs to have a unique URI that is drawn from a controlled vocabulary like the Citation Typing Ontology described above.  The `cito:documents` relationship is the default used by `insertRelationship`, so the relationship type doesn't need to be specified in this case. For example, with the example DataPackage created above, we can add the `cito:documents` relationship:

```{r}
dp <- insertRelationship(dp, subjectID=metadataId, objectIDs=sciId)
relations <- getRelationships(dp, condense=TRUE)
relations[grepl("cito:documents", relations$predicate),]
```

Relationships can be fully specified using the URI of the concept, as shown in the following statement that adds a provenance relationship between two objects in the example package:

```{r, eval=F}
dp <- insertRelationship(dp, subjectID=sciId2, objectIDs=sciId,
                   predicate="https://www.w3.org/ns/prov#wasDerivedFrom")
relations <- getRelationships(dp, condense=TRUE)
relations[grepl("prov:wasDerivedFrom", relations$predicate),]
``` 

The relationships contained in a DataPackage conform to the [Resource Description Framework](https://www.w3.org/RDF/) (RDF), which is a [World Wide Web Consortium](https://www.w3.org/)
standard for describing web accessible resources.

## Describing The Contents of a DataPackage 

In order to transport a DataPackage, for example to a data repository, a 
description of the contents of the DataPackage is created so that the consumer 
of the DataPackage can determine how to extract and process the contents.

A DataPackage can produce a standard description of its members and relationships 
which conforms to the Open Archives Initiative [Object Reuse and Exchange](https://www.openarchives.org/ore/) (OAI-ORE) specification, 
which is a widely used standard to describe aggregations of web accessible 
resources. This OAI-ORE description is referred to as a *resource map*.

The `serializePackage` method will create Resource Description Framework 
serialization of a resource map, written to a file in this case, that conforms 
to the OAI-ORE specification.

To create a resource map for the example DataPackage:

```{r, eval=FALSE}
tf <- tempfile()
packageId <- paste("urn:uuid:", UUIDgenerate(), sep="")
serializePackage(dp, file=tf, id=packageId)
```
This example writes to a tempfile using the default serialization format of 
"rdfxml". Also the URLs for each package member are prepended with the default 
value of the DataONE resolve service, which would be the URL that could be used 
to access this data object if the package is uploaded to a DataONE member node. 

A different value to be prepended to each identifier can be specified with the 
`resoveURI` argument. To specify that no value be prepended to the identifier 
URLs, specify a zero-length character:

```{r, eval=FALSE}
tf <- tempfile()
packageId <- paste("urn:uuid:", UUIDgenerate(), sep="")
serializePackage(dp, file=tf, id=packageId, resolveURI="")
```

It is also possible to create a JSON serialization, if desired:

```{r, eval=FALSE}
tf <- tempfile()
packageId <- paste("urn:uuid:", UUIDgenerate(), sep="")
serializePackage(dp, file=tf, id=packageId, syntaxName="json", mimeType="application/json", resolveURI="")
```

## Saving DataPackage Contents to a File

The contents of a DataPackage can be saved to a file using the `serializeToBagIt` 
method. This creates a *BagIt* file, which is a hierarchical file packaging format. 

The created BagIt file contains the data from the DataPackage members as well 
as an OAI-ORE resource map that is automatically created by `serializeToBagIt`.
 
The following R command shows how to create the BagIt file for the example DataPackage:

```{r, eval=F}
bagitFilename <- serializeToBagIt(dp)
```

The variable `bagitFilename` contains the file path to the temporary BagIt file.
This file should be copied to another location before quitting or restarting R:

```{r, eval=F}
file.copy(bagitFilename, "~/myPackageFile.zip")
```

This serialized BagIt version of the file is an excellent way to transport all 
of the files and metadata of a DataPackage to a data repository or a collaborator.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\name{parseSystemMetadata}
\alias{parseSystemMetadata}
\alias{parseSystemMetadata,SystemMetadata-method}
\title{Parse an external XML document and populate a SystemMetadata object with the parsed data}
\usage{
parseSystemMetadata(x, ...)

\S4method{parseSystemMetadata}{SystemMetadata}(x, xml, ...)
}
\arguments{
\item{x}{The \code{SystemMetadata} object}

\item{...}{Additional arguments passed to other functions or methods}

\item{xml}{The XML representation of the capabilities, as an XMLInternalElementNode}
}
\value{
the SystemMetadata object representing an object
}
\description{
Parse an XML representation of system metadata, and set the object slots of a SystemMetadata object 
the with obtained values.
}
\examples{
library(XML)
doc <- xmlParseDoc(system.file("testfiles/sysmeta.xml", package="datapack"), asText=FALSE)
sysmeta <- new("SystemMetadata")
sysmeta <- parseSystemMetadata(sysmeta, xmlRoot(doc))
}
\seealso{
\code{\link{SystemMetadata-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R, R/DataPackage.R
\name{getData}
\alias{getData}
\alias{getData,DataObject-method}
\alias{getData,DataPackage-method}
\title{Get the data content of a specified data object}
\usage{
getData(x, ...)

\S4method{getData}{DataObject}(x)

\S4method{getData}{DataPackage}(x, id)
}
\arguments{
\item{x}{DataObject or DataPackage: the data structure from where to get the data}

\item{...}{Additional arguments}

\item{id}{Missing or character: if \code{'x'} is DataPackage, the identifier of the
package member to get data from}
}
\value{
raw representation of the data
}
\description{
Get the data content of a specified data object
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB")
bytes <- getData(do)
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do1 <- new("DataObject", id="id1", data, format="text/csv", user="smith", mnNodeId="urn:node:KNB")
dp <- addMember(dp, do1)
bytes <- getData(dp, "id1")
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{serializeRDF}
\alias{serializeRDF}
\alias{serializeRDF,ResourceMap-method}
\title{Serialize a ResouceMap.}
\usage{
serializeRDF(x, ...)

\S4method{serializeRDF}{ResourceMap}(
  x,
  file,
  syntaxName = "rdfxml",
  mimeType = "application/rdf+xml",
  namespaces = data.frame(namespace = character(), prefix = character(),
    stringsAsFactors = FALSE),
  syntaxURI = NA_character_
)
}
\arguments{
\item{x}{a ResourceMap}

\item{...}{Additional parameters}

\item{file}{the file to which the ResourceMap will be serialized}

\item{syntaxName}{name of the syntax to use for serialization - default is "rdfxml"}

\item{mimeType}{the mimetype of the serialized output - the default is "application/rdf+xml"}

\item{namespaces}{a data frame containing one or more namespaces and their associated prefix}

\item{syntaxURI}{A URI of the serialized syntax}
}
\value{
status of the serialization (non)
}
\description{
The Redland RDF library is used to serialize the ResourceMap RDF model
to a file as RDF/XML.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do1 <- new("DataObject", id="id1", data, format="text/csv")
do2 <- new("DataObject", id="id2", data, format="text/csv")
dp <- addMember(dp, do1)
dp <- addMember(dp, do2)
dp <- insertRelationship(dp, subjectID="id1", objectIDs="id2", 
  predicate="http://www.w3.org/ns/prov#wasDerivedFrom")
relations <- getRelationships(dp)
resmap <- new("ResourceMap")
resmap <- createFromTriples(resmap, relations, id="myuniqueid")
\dontrun{
tf <- tempfile(fileext=".xml")
serializeRDF(resmap, tf)
}
}
\seealso{
\code{\link{ResourceMap-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{recordDerivation}
\alias{recordDerivation}
\alias{recordDerivation,DataPackage-method}
\title{Record derivation relationships between objects in a DataPackage}
\usage{
recordDerivation(x, ...)

\S4method{recordDerivation}{DataPackage}(x, sourceID, derivedIDs, ...)
}
\arguments{
\item{x}{a DataPackage object}

\item{...}{Additional parameters}

\item{sourceID}{the identifier of the source object in the relationship}

\item{derivedIDs}{an identifier or list of identifiers of objects that were derived from the source}
}
\description{
Record a derivation relationship that expresses that a target object has been derived from a source object.
For use with DataONE, a best practice is to specify the subject and predicate as DataONE persistent identifiers 
(https://mule1.dataone.org/ArchitectureDocs-current/design/PIDs.html). If the objects are not known to DataONE, then local identifiers can be
used, and these local identifiers may be promoted to DataONE PIDs when the package is uploaded to a DataONE member node.
}
\details{
A derived relationship is created for each value in the list "objectIDs".  For each derivedId, one statement will be
added expressing that it was derived from the sourceId.  The predicate is will be an RDF property (as a IRI) from the W3C PROV
specification, namely, "http://www.w3.org/ns/prov#wasDerivedFrom"
}
\examples{
\dontrun{
dp <- new("DataPackage")
recordDerivation(dp, "doi:1234/_030MXTI009R00_20030812.40.1", 
                 "doi:1234/_030MXTI009R00_20030812.45.1")
                     }
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{freeResourceMap}
\alias{freeResourceMap}
\alias{freeResourceMap,ResourceMap-method}
\title{Free memory used by a ResouceMap.}
\usage{
freeResourceMap(x)

\S4method{freeResourceMap}{ResourceMap}(x)
}
\arguments{
\item{x}{a ResourceMap}
}
\description{
The resources allocated by the redland RDF package are freed. The ResourceMap
object should be deleted immediately following this call.
}
\seealso{
\code{\link{ResourceMap-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{setValue}
\alias{setValue}
\alias{setValue,DataPackage-method}
\title{Set values for selected DataPackage members.}
\usage{
setValue(x, ...)

\S4method{setValue}{DataPackage}(x, name, value, identifiers = NA_character_, ...)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(Not yet used)}

\item{name}{A DataObject slot name.}

\item{value}{A new value to assign to the slot for selected DataPackage members.}

\item{identifiers}{A list of identifiers of DataPackage members to update.}
}
\value{
A DataPackage with possibly updated DataObjects.
}
\description{
The \code{'setValue'} method is used to modify values stored in DataPackage members.
Each member in a DataPackage is a DataObject which is an R S4 object that contains a set of values (slots).
The available slots are described at \code{help("DataObject-class")}.
}
\details{
If the parameter \code{identifiers} is provided, then DataPackage members that
have identifiers specified in the list will be updated. If this parameter is not provided
then no members will be updated. To update all members in a package, specify the
value of \code{identifiers=getIdentifiers(pkg)} where \code{pkg} is the variable name
of the DataPackage to update. Note that this method can be used to update the
\code{data} or \code{filenane} slots, but it is instead recommended to us the
\code{replaceMember} method to achieve this, as the \code{replaceMember} method assists 
in properly setting the related SystemMetadata values.
}
\examples{
# First create a package that we can modify. 
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
# The next statment sets the format type incorrectly as an example, so we can correct it later
do <- new("DataObject", id="myNewId", dataobj=data, format="image/jpg", user="jsmith")
dp <- addMember(dp, do)
data <- charToRaw("7,8.9\n4,10,11")
# This next statement also sets the format type incorrectly
do <- new("DataObject", id="myNewId2", dataobj=data, format="image/jpg", user="jsmith")
dp <- addMember(dp, do)
# Change format types to correct value for both package members
# Careful! Specifying 'identifiers=getIdentifiers(dp) will update all package members!
dp <- setValue(dp, name="sysmeta@formatId", value="text/csv", identifiers=getIdentifiers(dp))
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{updateXML}
\alias{updateXML}
\alias{updateXML,DataObject-method}
\title{Update selected elements of the XML content of a DataObject}
\usage{
updateXML(x, ...)

\S4method{updateXML}{DataObject}(x, xpath = NA_character_, replacement = NA_character_, ...)
}
\arguments{
\item{x}{A DataObject instance}

\item{...}{Additional parameters (not yet used)}

\item{xpath}{A \code{character} value specifying the location in the XML to update.}

\item{replacement}{A \code{character} value that will replace the elements found with the \code{xpath}.}
}
\value{
The modified DataObject
}
\description{
The data content of the DataObject is updated by using the \code{xpath} 
argument to locate the elements to update with the character value specified in the 
\code{replacement} argument.
}
\examples{
\dontrun{
library(datapack)
dataObj <- new("DataObject", format="text/csv", file=sampleData)
sampleEML <- system.file("extdata/sample-eml.xml", package="datapack")
dataObj <- updateMetadata(dataObj, xpath="", replacement=)
}
library(datapack)
# Create the metadata object with a sample EML file
sampleMeta <- system.file("./extdata/sample-eml.xml", package="datapack")
metaObj <- new("DataObject", format="eml://ecoinformatics.org/eml-2.1.1", file=sampleMeta)
# In the metadata object, replace "sample-data.csv" with 'sample-data.csv.zip'
xp <- sprintf("//dataTable/physical/objectName[text()=\"\%s\"]", "sample-data.csv")
metaObj <- updateXML(metaObj, xpath=xp, replacement="sample-data.csv.zip")
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{initialize,DataPackage-method}
\alias{initialize,DataPackage-method}
\alias{DataPackage-initialize}
\title{Initialize a DataPackage object.}
\usage{
\S4method{initialize}{DataPackage}(.Object, packageId)
}
\arguments{
\item{.Object}{The object being initialized}

\item{packageId}{The package id to assign to the package}
}
\description{
Initialize a DataPackage object.
}
\examples{
# Create a DataPackage with undefined package id (to be set manually later)
pkg <- new("DataPackage")
# Alternatively, manually assign the package id when the DataPackage object is created
pkg <- new("DataPackage", "urn:uuid:4f953288-f593-49a1-adc2-5881f815e946")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datapack-package.r
\docType{package}
\name{datapack}
\alias{datapack}
\title{datapack, a container for packages of data and associated metadata}
\description{
The datapack R package provides an abstraction for collating 
heterogeneous collections of data objects and metadata into a bundle that can 
be transported and loaded into a single composite file.  The methods in
this package provide a convenient way
to load data from common repositories such as DataONE into the R environment, 
and to document, serialize, and save data from R to data repositories worldwide. 
A data package is represented as an instance of the S4 class \code{\link[=DataPackage-class]{DataPackage}}, which 
consists of one or more instances of the S4 DataObject class, which in turn contains
an instance of the S4 SystemMetadata class.  The SystemMetadata
class provides critical metadata about a data object that is needed to transport
it to an external repository, including the identifier for the object, its
format, its checksum and size, and information about which repositories the
data is associated with.  DataPackages can be loaded from and saved to the 
DataONE federated network of repositories using the dataone package, but they 
can also be used as standalone transport containers for other systems.

A DataPackage includes a manifest based on the OAI-ORE 
specification for describing aggregations of files as a ResourceMap. 
Resource maps are RDF documents that conform to the Open Archives Initiative
Object Reuse and Exchange (OAI-ORE) specification. Resource maps are generated 
by data providers to define data packages, and have a namespace of 
http://www.openarchives.org/ore/terms/.

A DataPackage is serialized as a zip file following the BagIt RFC specification,
which provides a consistent mechanism for a serialized representation of a 
group of opaque objects in a predictable structure. BagIt includes a 
specification for including metadata about each of the objects, the bag itself, 
and fixity attributes so that any BagIt implementation can validate the 
components contained within a package.  When expanded, a BagIt zipfile will
expand to a common directory structure with a predictable set of metadata that
describes the structure and content of the bag.  Conformance with the BagIt
specification is handled by the DataPackage class.
}
\section{Classes}{

\itemize{
 \item{\code{\link{DataPackage-class}}}{: A class representing a data package, which can contain data objects}
 \item{\code{\link{DataObject-class}}}{: DataObject wraps raw data with system-level metadata}
 \item{\code{\link{SystemMetadata-class}{SystemMetadata}}}{: A DataONE SystemMetadata object containing basic identification, ownership, access policy, replication policy, and related metadata.}
 \item{\code{\link{ResourceMap-class}{ResourceMap}}}{: ResourceMap provides methods to create, serialize and deserialize an OAI ORE resource map.}
}
}

\author{
Matthew B. Jones (NCEAS), Peter Slaughter (NCEAS)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\name{SystemMetadata}
\alias{SystemMetadata}
\alias{SystemMetadata,XMLInternalElementNode-method}
\title{Create DataONE SystemMetadata object}
\usage{
SystemMetadata(...)

\S4method{SystemMetadata}{XMLInternalElementNode}(x, ...)
}
\arguments{
\item{...}{Additional arguments}

\item{x}{A value of type \code{"XMLInternalElementNode"}, containing the parsed XML element with SystemMetadata fields.}
}
\description{
A class representing DataONE SystemMetadata, which is core information about objects stored in a repository
and needed to manage those objects across systems.  SystemMetadata contains basic identification, ownership,
access policy, replication policy, and related metadata.

If the *sysmeta* parameter is specified, then construct a new SystemMetadata instance by using the fields from 
an XML representation of the SystemMetadata.
}
\seealso{
\code{\link{SystemMetadata-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{createFromTriples}
\alias{createFromTriples}
\alias{createFromTriples,ResourceMap-method}
\title{Populate a ResourceMap with RDF relationships from data.frame.}
\usage{
createFromTriples(x, ...)

\S4method{createFromTriples}{ResourceMap}(
  x,
  relations,
  identifiers,
  resolveURI = NA_character_,
  externalIdentifiers = list(),
  creator = NA_character_,
  ...
)
}
\arguments{
\item{x}{a ResourceMap}

\item{...}{(Additional parameters)}

\item{relations}{A data.frame to read relationships from}

\item{identifiers}{A list of the identifiers of data objects contained in the associated data package}

\item{resolveURI}{A character string containing a URI to prepend to datapackage identifiers.}

\item{externalIdentifiers}{A list of identifiers that are referenced from the package, but are not package members.}

\item{creator}{A \code{character} string containing the creator of the package.}
}
\description{
RDF relationships are added to a ResourceMap object from a data.frame that
contains RDF triples. For example, relationships can be exported from a DataPackage via
\code{\link{getRelationships}}. The resulting data.frame is then read by \code{createFromTriples}
to create the ResourceMap.
}
\details{
The \code{identifiers} parameter contains the identifiers of all data objects in the DataPackage.
For each data objects, additional relationships will be added that are required by the OAI-ORE specification,
for example a Dublin Core identifier statement is added. The resolveURI string value is prepended to 
DataPackage member identifiers in the resulting resource map. If no resolveURI value
is specified, then 'https://cn.dataone.org/cn/v1/resolve' is used.
}
\examples{
library(datapack)
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do1 <- new("DataObject", id="id1", data, format="text/csv")
do2 <- new("DataObject", id="id2", data, format="text/csv")
dp <- addMember(dp, do1)
dp <- addMember(dp, do2)
dp <- insertRelationship(dp, subjectID="id1", objectIDs="id2", 
  predicate="http://www.w3.org/ns/prov#wasDerivedFrom")
relations <- getRelationships(dp)
resMapId <- sprintf("\%s\%s", "resourceMap_", uuid::UUIDgenerate())  
resMap <- new("ResourceMap", id=resMapId)
resMap <- createFromTriples(resMap, relations, getIdentifiers(dp)) 
}
\seealso{
\code{\link{ResourceMap-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{initialize,DataObject-method}
\alias{initialize,DataObject-method}
\alias{DataObject-initialize}
\title{Initialize a DataObject}
\usage{
\S4method{initialize}{DataObject}(
  .Object,
  id = NA_character_,
  dataobj = NA,
  format = NA_character_,
  user = NA_character_,
  mnNodeId = NA_character_,
  filename = NA_character_,
  seriesId = NA_character_,
  mediaType = NA_character_,
  mediaTypeProperty = list(),
  dataURL = NA_character_,
  targetPath = NA_character_,
  checksumAlgorithm = "SHA-256"
)
}
\arguments{
\item{.Object}{the DataObject instance to be initialized}

\item{id}{the identifier for the DataObject, unique within its repository. Optionally this can be an existing SystemMetadata object}

\item{dataobj}{the bytes of the data for this object in \code{'raw'} format, optional if \code{'filename'} is provided}

\item{format}{the format identifier for the object, e.g."text/csv", "eml://ecoinformatics.org/eml-2.1.1"}

\item{user}{the identity of the user owning the package, typically in X.509 format}

\item{mnNodeId}{the node identifier for the repository to which this object belongs.}

\item{filename}{the filename for the fully qualified path to the data on disk, optional if \code{'data'} is provided}

\item{seriesId}{A unique string to identifier the latest of multiple revisions of the object.}

\item{mediaType}{The When specified, indicates the IANA Media Type (aka MIME-Type) of the object. The value should include the media type and subtype (e.g. text/csv).}

\item{mediaTypeProperty}{A list, indicates IANA Media Type properties to be associated with the parameter \code{"mediaType"}}

\item{dataURL}{A character string containing a URL to remote data (a repository) that this DataObject represents.}

\item{targetPath}{An optional string that denotes where the file should go in a downloaded package}

\item{checksumAlgorithm}{A character string specifying the checksum algorithm to use}
}
\description{
When initializing a DataObject using passed in data, one can either pass 
in the \code{'id'} param as a \code{'SystemMetadata'} object, or as a \code{'character'} string 
representing the identifier for an object along with parameters for format, user,and associated member node.
If \code{'data'} is not missing, the \code{'data'} param holds the \code{'raw'} data.  Otherwise, the
\code{'filename'} parameter must be provided, and points at a file containing the bytes of the data.
}
\details{
If filesystem storage is used for the data associated with a DataObject, care must be
taken to not modify or remove that file in R or via other facilities while the DataObject exists in the R session.
Changes to the object are not detected and will result in unexpected results. Also, if the \code{'dataobj'} parameter
is used to specify the data source, then \code{'filename'} argument may also be specified, but in this case 
the value \code{'filename'} parameter is used to tell DataONE the filename to create when this file is
downloaded from a repository.
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB", targetPath="data/rasters/data.tiff")
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R, R/DataObject.R,
%   R/DataPackage.R
\name{removeAccessRule}
\alias{removeAccessRule}
\alias{removeAccessRule,SystemMetadata-method}
\alias{removeAccessRule,DataObject-method}
\alias{removeAccessRule,DataPackage-method}
\title{Remove an access rule from the specified object.}
\usage{
removeAccessRule(x, ...)

\S4method{removeAccessRule}{SystemMetadata}(x, y, ...)

\S4method{removeAccessRule}{DataObject}(x, y, ...)

\S4method{removeAccessRule}{DataPackage}(x, y, permission = NA_character_, identifiers = list(), ...)
}
\arguments{
\item{x}{The object instance to which to remove the rule}

\item{...}{Additional arguments
\itemize{
  \item{permission The permission to be remove to subject if x is character (read, write, changePermission)}
}}

\item{y}{The subject of the rule to be removed, or a data.frame containing access rules.}

\item{permission}{The permission to remove, if parameter \code{x} is a character string containing a \code{subject}.}

\item{identifiers}{A list of \code{character} values containing package member identifiers that the access rule will be 
applied to (default is all package members).}
}
\value{
The SystemMetadata object with the updated access policy.

The DataObject object with the updated access policy.

The Datapackage with members having updated access policies.
}
\description{
Remove access rules from the access policy of the specified object.
}
\examples{
#
# Remove access rules from a SystemMetadata object.
# Parameter "y" can be character string containing the subject of the access rule:
sysmeta <- new("SystemMetadata")
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "write")
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")
sysmeta <- removeAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")

# Alternatively, parameter "y" can be a data.frame containing one or more access rules:
# Add write, changePermission for uid=jones,...
sysmeta <- addAccessRule(sysmeta, "uid=jones,ou=Account,dc=example,dc=com", "write")
sysmeta <- addAccessRule(sysmeta, "uid=jones,ou=Account,dc=example,dc=com", "changePermission")
# Now take privs for uid=jones,... away
accessRules <- data.frame(subject=c("uid=jones,ou=Account,dc=example,dc=com", 
                                     "uid=jones,ou=Account,dc=example,dc=com"), 
                                     permission=c("write", "changePermission"))
sysmeta <- removeAccessRule(sysmeta, accessRules)
#
# Remove access rules form a DataObject.
library(datapack)
do <- new("DataObject", file=system.file("./extdata/sample-data.csv", package="datapack"), 
                        format="text/csv")
do <- setPublicAccess(do)
isPublic <- hasAccessRule(do, "public", "read")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
                          "uid=wiggens,o=unaffiliated,dc=example,dc=org"), 
                          permission=c("write", "changePermission"), 
                          stringsAsFactors=FALSE)
do <- addAccessRule(do, accessRules)
do <- removeAccessRule(do, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")
# hasAccessRule should return FALSE
hasWrite <- hasAccessRule(do, "smith", "write")

# Alternatively, parameter "y" can be a data.frame containing one or more access rules:
do <- addAccessRule(do, "uid=smith,ou=Account,dc=example,dc=com", "write")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
  "uid=slaughter,o=unaffiliated,dc=example,dc=org"), 
  permission=c("write", "changePermission"))
sysmeta <- removeAccessRule(do, accessRules)
# 
# Remove access rules from a DataPackage.
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="id1", dataobj=data, format="text/csv")
dp <- addMember(dp, obj)
data2 <- charToRaw("7,8,9\n4,10,11\n")
obj2 <- new("DataObject", id="id2", dataobj=data2, format="text/csv")
dp <- addMember(dp, obj2)
# Add access rule to all package members
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "write")
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "changePermission" )
# Now take 'changePermission' away for user 'uid=smith...', specifying parameter 'y' 
# as a character string containing a 'subject'.
dp <- removeAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "write")
dp <- removeAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")

# Alternatively, parameter "y" can be a data.frame containing one or more access rules:
# Add write, changePermission for uid=jones,...
dp <- addAccessRule(dp, "uid=jones,ou=Account,dc=example,dc=com", "write")
dp <- addAccessRule(dp, "uid=jones,ou=Account,dc=example,dc=com", "changePermission")
# Now take privs for uid=jones,... away
accessRules <- data.frame(subject=c("uid=jones,ou=Account,dc=example,dc=com", 
                                     "uid=jones,ou=Account,dc=example,dc=com"), 
                                     permission=c("write", "changePermission"))
dp <- removeAccessRule(dp, accessRules)
}
\seealso{
\code{\link{SystemMetadata-class}}

\code{\link{DataObject-class}}

\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{updateRelationships}
\alias{updateRelationships}
\alias{updateRelationships,DataPackage-method}
\title{Update package relationships by replacing an old identifier with a new one.}
\usage{
updateRelationships(x, ...)

\S4method{updateRelationships}{DataPackage}(x, id, newId, ...)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{(Not yet used)}

\item{id}{A character value containing the identifier to be replaced.}

\item{newId}{A character value containing the identifier that will replace the old identifier.}
}
\description{
When package members are updated, they receive a new identifier (replaceMember). It is therefor
necessary to update the package relationships to update occurrences of the old identifier
with the new one when the old identifier appears in the "subject" or "object" of a 
relationship.
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{removeMember}
\alias{removeMember}
\alias{removeMember,DataPackage-method}
\title{Remove the Specified Member from the Package}
\usage{
removeMember(x, ...)

\S4method{removeMember}{DataPackage}(x, do, removeRelationships = FALSE)
}
\arguments{
\item{x}{a DataPackage object}

\item{...}{(Not yet used)}

\item{do}{The package member to remove, either as a \code{"DataObject"} or \code{"character"} (for the object identifier)}

\item{removeRelationships}{A \code{logical} value. If TRUE, package relationships for this package member are removed. Default is FALSE.}
}
\description{
Given the identifier of a DataObject in a DataPackage, delete the DataObject
from the DataPackage.
}
\details{
The \code{removeMember} method removes the specified DataObject from the DataPackage. In 
addition, any package relationships that included the DataObject are removed.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", id="myNewId", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
# Remove the package member and any provenance relationships that reference it.
removeMember(dp, "myNewId", removeRelationships=TRUE)
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{describeWorkflow}
\alias{describeWorkflow}
\alias{describeWorkflow,DataPackage-method}
\title{Add data derivation information to a DataPackage}
\usage{
describeWorkflow(x, ...)

\S4method{describeWorkflow}{DataPackage}(
  x,
  sources = list(),
  program = NA_character_,
  derivations = list(),
  insertDerivations = TRUE,
  ...
)
}
\arguments{
\item{x}{The \code{DataPackage} to add provenance relationships to.}

\item{...}{Additional parameters}

\item{sources}{A list of DataObjects for files that were read by the program. Alternatively, a list 
of DataObject identifiers can be specified as a list of character strings.}

\item{program}{The DataObject created for the program such as an R script. Alternatively the DataObject identifier can
be specified.}

\item{derivations}{A list of DataObjects for files that were generated by the program. Alternatively, a list 
of DataObject identifiers can be specified as a list of character strings.}

\item{insertDerivations}{A \code{logical} value. If TRUE then the provenance relationship 
\code{prov:wasDerivedFrom} will be used to connect every source and derivation. The default value 
is TRUE.}
}
\description{
Add information about the relationships among DataObject members 
in a DataPackage, retrospectively describing the way in which derived data were 
created from source data using a processing program such as an R script.  These provenance
relationships allow the derived data to be understood sufficiently for users
to be able to reproduce the computations that created the derived data, and to
trace lineage of the derived data objects. The method \code{describeWorkflow} 
will add provenance relationships between a script that was executed, the files 
that it used as sources, and the derived files that it generated.
}
\details{
This method operates on a DataPackage that has had DataObjects for 
the script, data sources (inputs), and data derivations (outputs) previously 
added to it, or can reference identifiers for objects that exist in other DataPackage
instances. This allows a user to create a standalone package that contains all of
its source, script, and derived data, or a set of data packages that are chained
together via a set of derivation relationships between the members of those packages.

Provenance relationships are described following the the ProvONE data model, which
can be viewed at \url{https://purl.dataone.org/provone-v1-dev}.  In particular, 
the following relationships are inserted (among others):
\itemize{
 \item{\code{prov:used}} {indicates which source data was used by a program execution}
 \item{\code{prov:generatedBy}} {indicates which derived data was created by a program execution}
 \item{\code{prov:wasDerivedFrom}} {indicates the source data from which derived data were created using the program}
}
}
\examples{
library(datapack)
dp <- new("DataPackage")
# Add the script to the DataPackage
progFile <- system.file("./extdata/pkg-example/logit-regression-example.R", package="datapack")
progObj <- new("DataObject", format="application/R", filename=progFile)
dp <- addMember(dp, progObj)

# Add a script input to the DataPackage
inFile <- system.file("./extdata/pkg-example/binary.csv", package="datapack") 
inObj <- new("DataObject", format="text/csv", filename=inFile)
dp <- addMember(dp, inObj)

# Add a script output to the DataPackage
outFile <- system.file("./extdata/pkg-example/gre-predicted.png", package="datapack")
outObj <- new("DataObject", format="image/png", file=outFile)
dp <- addMember(dp, outObj)

# Add the provenenace relationshps, linking the input and output to the script execution
# Note: 'sources' and 'derivations' can also be lists of "DataObjects" or "DataObject' identifiers
dp <- describeWorkflow(dp, sources = inObj, program = progObj, derivations = outObj) 
# View the results
head(getRelationships(dp))
}
\seealso{
The R 'recordr' package for run-time recording of provenance relationships.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{removeRelationships}
\alias{removeRelationships}
\alias{removeRelationships,DataPackage-method}
\title{Remove relationships of objects in a DataPackage}
\usage{
removeRelationships(x, ...)

\S4method{removeRelationships}{DataPackage}(x, subjectID = NA_character_, predicate = NA_character_)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{(Additional parameters)}

\item{subjectID}{The identifier of the subject of the relationships to be removed}

\item{predicate}{The identifier of the predicate of the relationships to be removed}
}
\value{
the updated DataPackage object
}
\description{
Use this function to remove all or a subset of the relationships that have previously been added in a data package.
}
\details{
Remove a relationship of the form "subject -> predicate -> object", as defined by the Resource Description Framework (RDF), i.e.
an RDF triple. If neither subjectID nor predicate are provided, then all relationships are removed.  If one or both
are provided, they are used to select matching triples to be removed.
Note: This method updates the passed-in DataPackage object.
}
\examples{
dp <- new("DataPackage")
# Create a relationship
dp <- insertRelationship(dp, "/Users/smith/scripts/genFields.R",
    "https://knb.org/data_20030812.40.1",
    "http://www.w3.org/ns/prov#used")
# Create a relationshp with the subject as a blank node with an automatically assigned blank 
# node id
dp <- insertRelationship(dp, subjectID=NA_character_, objectIDs="thing6", 
    predicate="http://myns.org/wasThing")
# Create a relationshp with the subject as a blank node with a user assigned blank node id
dp <- insertRelationship(dp, subjectID="urn:uuid:bc9e160e-ca21-47d5-871b-4a4820fe4451", 
      objectIDs="thing7", predicate="http://myns.org/hadThing")
# Create multiple relationships with the same subject, predicate, but different objects
dp <- insertRelationship(dp, subjectID="https://myns.org/subject1", 
      objectIDs=c("thing4", "thing5"), predicate="http://myns.org/hadThing")
# Create multiple relationships with subject and object types specified
dp <- insertRelationship(dp, subjectID="orcid.org/0000-0002-2192-403X", 
    objectIDs="http://www.example.com/home", predicate="http://myns.org/hadHome",
                   subjectType="uri", objectType="literal")
nrow(getRelationships(dp)) 
dp <- removeRelationships(dp, predicate='http://myns.org/wasThing')
nrow(getRelationships(dp)) 
dp <- removeRelationships(dp, subjectID='orcid.org/0000-0002-2192-403X')
nrow(getRelationships(dp)) 
dp <- removeRelationships(dp, subjectID='https://myns.org/subject1', 
    predicate='http://myns.org/hadThing')
nrow(getRelationships(dp)) 
dp <- removeRelationships(dp)
nrow(getRelationships(dp)) 
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\name{validate}
\alias{validate}
\alias{validate,SystemMetadata-method}
\title{Validate a SystemMetadata object.}
\usage{
validate(x, ...)

\S4method{validate}{SystemMetadata}(x, ...)
}
\arguments{
\item{x}{the instance to be validated}

\item{...}{(Additional parameters)}
}
\value{
logical, \code{TRUE} if the SystemMetadata object is valid, else a list of strings detailing errors
}
\description{
Validate a system metadata object, ensuring that required fields are present and of the right type.
}
\examples{
library(XML)
doc <- xmlParseDoc(system.file("testfiles/sysmeta.xml", package="datapack"), asText=FALSE)
sysmeta <- new("SystemMetadata")
sysmeta <- parseSystemMetadata(sysmeta, xmlRoot(doc))
valid <- validate(sysmeta)
}
\seealso{
\code{\link{SystemMetadata-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{getFormatId}
\alias{getFormatId}
\alias{getFormatId,DataObject-method}
\title{Get the FormatId of the DataObject}
\usage{
getFormatId(x, ...)

\S4method{getFormatId}{DataObject}(x)
}
\arguments{
\item{x}{DataObject}

\item{...}{(not yet used)}
}
\value{
the formatId
}
\description{
Get the FormatId of the DataObject
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB")
fmtId <- getFormatId(do)
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{containsId}
\alias{containsId}
\alias{containsId,DataPackage-method}
\title{Returns true if the specified object is a member of the package}
\usage{
containsId(x, ...)

\S4method{containsId}{DataPackage}(x, identifier)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{(Not yet used)}

\item{identifier}{The DataObject identifier to check for inclusion in the DataPackage}
}
\value{
A logical - a value of TRUE indicates that the DataObject is in the DataPackage
}
\description{
Returns true if the specified object is a member of the package
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
id <- "myNewId"
do <- new("DataObject", id=id, dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
isInPackage <- containsId(dp, identifier="myNewId")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{getIdentifiers}
\alias{getIdentifiers}
\alias{getIdentifiers,DataPackage-method}
\title{Get the Identifiers of Package Members}
\usage{
getIdentifiers(x, ...)

\S4method{getIdentifiers}{DataPackage}(x)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(not yet used)}
}
\value{
A list of identifiers
}
\description{
The identifiers of the objects in the package are retrieved and returned as a list.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
getIdentifiers(dp)
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R, R/DataObject.R,
%   R/DataPackage.R
\name{clearAccessPolicy}
\alias{clearAccessPolicy}
\alias{clearAccessPolicy,SystemMetadata-method}
\alias{clearAccessPolicy,DataObject-method}
\alias{clearAccessPolicy,DataPackage-method}
\title{Clear the accessPolicy from the specified object.}
\usage{
clearAccessPolicy(x, ...)

\S4method{clearAccessPolicy}{SystemMetadata}(x, ...)

\S4method{clearAccessPolicy}{DataObject}(x, ...)

\S4method{clearAccessPolicy}{DataPackage}(x, identifiers = list(), ...)
}
\arguments{
\item{x}{the instance to clear access rules from.}

\item{...}{(Additional parameters)}

\item{identifiers}{A list of \code{character} values containing package member identifiers that the access rule will be applied to.}
}
\value{
The SystemMetadata object with the cleared access policy.

The DataObject with the cleared access policy.

The SystemMetadata object with the cleared access policy.
}
\description{
Clears the accessPolicy from the specified object by overwriting
all existing access rules set on the object with an empty set.
}
\examples{
# Clear access policy for a SystemMetadata object.
sysmeta <- new("SystemMetadata")
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "write")
sysmeta <- clearAccessPolicy(sysmeta)
# Clear access policy for a DataObject
do <- new("DataObject", format="text/csv", filename=system.file("extdata/sample-data.csv", 
          package="datapack"))
do <- addAccessRule(do, "uid=smith,ou=Account,dc=example,dc=com", "write")
do <- clearAccessPolicy(do)
# Clear access policy for a DataPackage
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", dataobj=data, format="text/csv")
dp <- addMember(dp, obj)
data2 <- charToRaw("7,8,9\n4,10,11\n")
obj2 <- new("DataObject", dataobj=data2, format="text/csv")
dp <- addMember(dp, obj2)

# Add the access rule to all package members
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", 
    permission="write")
# Now clear the access policy for just the second object 
dp <- clearAccessPolicy(dp, getIdentifier(obj2))

}
\seealso{
\code{\link{SystemMetadata-class}}

\code{\link{DataObject-class}}

\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{initialize,ResourceMap-method}
\alias{initialize,ResourceMap-method}
\alias{ResourceMap-initialize}
\title{Initialize a ResourceMap object.}
\usage{
\S4method{initialize}{ResourceMap}(.Object, id = NA_character_)
}
\arguments{
\item{.Object}{a ResourceMap object}

\item{id}{a unique identifier to identify this ResourceMap. This id will be used internally in the ResourceMap.}
}
\value{
the ResourceMap object
}
\description{
Create a ResourceMap object that contains relationships (RDF triples) of objects in the DataPackage.
}
\seealso{
\code{\link{ResourceMap-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{parseRDF}
\alias{parseRDF}
\alias{parseRDF,ResourceMap-method}
\title{Parse an RDF/XML resource map from a file.}
\usage{
parseRDF(x, rdf, ...)

\S4method{parseRDF}{ResourceMap}(
  x,
  rdf,
  asText = FALSE,
  name = "rdfxml",
  mimeType = "application/rdf+xml",
  ...
)
}
\arguments{
\item{x}{ResourceMap}

\item{rdf}{A file or character value containing a resource map that will be parsed into the ResourceMap object}

\item{...}{Additional parameters (not yet used).}

\item{asText}{A logical value. If TRUE, then the 'rdf' parameter is a character vector, if FALSE then it is the name of a file to read.}

\item{name}{The name of the RDF xml parser, the default is "rdfxml".}

\item{mimeType}{A character value containing the RDF format type. The default is "application/rdf+xml".}
}
\value{
x the ResourceMap containing the parsed RDF/XML content
}
\description{
parseRDF reads a file containing an RDF model in RDF/XML format and initializes
a ResourceMap based on this content.
}
\details{
This method resets the slot ResourceMap@world so any previously stored triples are discarded, allowing
for a clean model object in which to parse the new RDF content into. It is assumed that the content is a
valid ORE resource map therefor no validation checks specific to the OAI-ORE content model are performed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{getIdentifier}
\alias{getIdentifier}
\alias{getIdentifier,DataObject-method}
\title{Get the Identifier of the DataObject}
\usage{
getIdentifier(x, ...)

\S4method{getIdentifier}{DataObject}(x)
}
\arguments{
\item{x}{DataObject}

\item{...}{(not yet used)}
}
\value{
the identifier
}
\description{
Get the Identifier of the DataObject
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB")
id <- getIdentifier(do)
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{serializePackage}
\alias{serializePackage}
\alias{serializePackage,DataPackage-method}
\title{Create an OAI-ORE resource map from the package}
\usage{
serializePackage(x, ...)

\S4method{serializePackage}{DataPackage}(
  x,
  file,
  id = NA_character_,
  syntaxName = "rdfxml",
  mimeType = "application/rdf+xml",
  namespaces = data.frame(namespace = character(), prefix = character(),
    stringsAsFactors = FALSE),
  syntaxURI = NA_character_,
  resolveURI = NA_character_,
  creator = NA_character_
)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{Additional arguments}

\item{file}{The file to which the ResourceMap will be serialized}

\item{id}{A unique identifier for the serialization. The default value is the id assigned
to the DataPackage when it was created.}

\item{syntaxName}{The name of the syntax to use for serialization - default is "rdfxml"}

\item{mimeType}{The mimetype of the serialized output - the default is "application/rdf+xml"}

\item{namespaces}{A data frame containing one or more namespaces and their associated prefix}

\item{syntaxURI}{URI of the serialization syntax}

\item{resolveURI}{A character string containing a URI to prepend to datapackage identifiers}

\item{creator}{A \code{character} string containing the creator of the package.}
}
\description{
The DataPackage is serialized as a OAI-ORE resource map to the specified file.
}
\details{
The resource map that is created is serialized by default as RDF/XML. Other serialization formats
can be specified using the \code{syntaxName} and \code{mimeType} parameters. Other available formats
include: 
\tabular{ll}{
\strong{syntaxName}\tab \strong{mimeType}\cr
json\tab application/json\cr
ntriples\tab application/n-triples\cr
turtle\tab text/turtle\cr
dot\tab text/x-graphviz\cr
} 
Note that the \code{syntaxName} and \code{mimeType} arguments together specify o serialization format.

Also, for packages that will be uploaded to the DataONE network, "rdfxml" is the only 
accepted format.  

The resolveURI string value is prepended to DataPackage member identifiers in the resulting resource map. 
If no resolveURI value is specified, then 'https://cn.dataone.org/cn/v1/resolve' is used.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", id="do1", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
data2 <- charToRaw("7,8,9\n10,11,12")
do2 <- new("DataObject", id="do2", dataobj=data2, format="text/csv", user="jsmith")
dp <- addMember(dp, do2)
dp <- describeWorkflow(dp, sources=do, derivations=do2)
\dontrun{
td <- tempdir()
status <- serializePackage(dp, file=paste(td, "resmap.json", sep="/"), syntaxName="json",  
    mimeType="application/json")
status <- serializePackage(dp, file=paste(td, "resmap.xml", sep="/"), syntaxName="rdfxml", 
    mimeType="application/rdf+xml")
status <- serializePackage(dp, file=paste(td, "resmap.ttl", sep="/"), syntaxName="turtle", 
    mimeType="text/turtle")
}
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R, R/DataObject.R,
%   R/DataPackage.R
\name{addAccessRule}
\alias{addAccessRule}
\alias{addAccessRule,SystemMetadata-method}
\alias{addAccessRule,DataObject-method}
\alias{addAccessRule,DataPackage-method}
\title{Add access rules to the specified object.}
\usage{
addAccessRule(x, ...)

\S4method{addAccessRule}{SystemMetadata}(x, y, ...)

\S4method{addAccessRule}{DataObject}(x, y, ...)

\S4method{addAccessRule}{DataPackage}(x, y, ...)
}
\arguments{
\item{x}{The object instance to which to add the rules}

\item{...}{Additional arguments
\itemize{
  \item{permission The permission to be applied to subject if x is character (read, write, changePermission)}
}}

\item{y}{The subject of the rule to be added, or a data frame of subject/permission tuples}
}
\value{
The SystemMetadata object with the updated access policy.

The DataObject with the updated access policy

The DataPackage with updated DataObject access policies
}
\description{
Add one or more access rules to the access policy of the specified object.
}
\details{
If the \code{y} argument is specified as a character string containing a \code{subject},
then an optional \code{permission} parameter must be specified, that contains a character list
specifying the permissions to add for each \code{subject}.

Note that when \code{addAccessRule} is called with a `DataPackage` argument, the 
additional parameter \code{identifiers} can be used:
\itemize{
  \item{identifiers A list of \code{character} values containing package member identifiers that the access rule will be applied to (all members is the default)}.
}
}
\examples{
# Add an access rule to a SystemMetadata access policy.
# Parameter "y" can be character string containing the subject of the access rule:
sysmeta <- new("SystemMetadata")
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "write")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
  "uid=slaughter,o=unaffiliated,dc=example,dc=org"), permission=c("write", "changePermission"))
sysmeta <- addAccessRule(sysmeta, accessRules)
# Alternatively, parameter "y" can be a data.frame containing one or more access rules:
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "write")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
  "uid=slaughter,o=unaffiliated,dc=example,dc=org"), permission=c("write", "changePermission"))
sysmeta <- addAccessRule(sysmeta, accessRules)
# Add an access rule to a DataObject
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="1234", dataobj=data, format="text/csv")
obj <- addAccessRule(obj, "uid=smith,ou=Account,dc=example,dc=com", "write")
# Add an access rule to members of a DataPackage
# First create a sample DataPackage
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="id1", dataobj=data, format="text/csv")
dp <- addMember(dp, obj)
data2 <- charToRaw("7,8,9\n4,10,11\n")
obj2 <- new("DataObject", id="id2", dataobj=data2, format="text/csv")
dp <- addMember(dp, obj2)
# Add access rule to all package members
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "write", getIdentifiers(dp))
}
\seealso{
\code{\link{SystemMetadata-class}}

\code{\link{DataObject-class}}

\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\docType{class}
\name{DataPackage-class}
\alias{DataPackage-class}
\title{A class representing a data package}
\description{
The DataPackage class provides methods for adding and extracting
data objects from a data package. The contents of a data package
can include arbitrary types of objects, including data files, program code,
visualizations and images, animations, and any other type of file. The DataPackage class
stores the individual members of the data package along with key system-level metadata
about each object, including its size, checksum, identifier, and other key information
needed to effectively archive the members of the package.  In addition, the
DataPackage class can include key provenance metadata about the relationships among
the objects in the data package.  For example, the data package can document that one object
provides documentation for another (\code{cito:documents}), and that one object was
derived from another (\code{prov:wasDerivedFrom}) by executing a program that 
used source data (\code{prov:used}) to create a derived data object 
{\code{prov:wasGeneratedBy}}.  These relationships are integral to the data package,
and can be visualized by programs that understand the ProvONE provenance 
model (see \url{https://purl.dataone.org/provone-v1-dev}). 

The DataPackage class is an R representation of an underlying Open Archives 
Initiative ORE model (Object Reuse and Exchange; 
see \url{https://www.openarchives.org/ore/}), and follows the DataONE Data
Packaging model
(see \url{https://releases.dataone.org/online/api-documentation-v2.0.1/design/DataPackage.html}).
}
\section{Slots}{

\describe{
\item{\code{relations}}{A hash containing provenance relationships of package objects}

\item{\code{objects}}{A hash containing identifiers for objects in the DataPackage}

\item{\code{sysmeta}}{A SystemMetadata class instance describing the package}

\item{\code{externalIds}}{A list containing identifiers for objects associated with the DataPackage}

\item{\code{resmapId}}{A character string specifying the identifier for the package resource map. 
This is assigned after a package is uploaded or downloaded from a repository.}
}}

\section{Methods}{

\itemize{
 \item{\code{\link[=DataPackage-initialize]{initialize}}}{: Initialize a DataPackage object.}
 \item{\code{\link{addAccessRule}}}{: Add access rules to DataObjects in a DataPackage.}
 \item{\code{\link{addMember}}}{: Add a DataObject to a DataPackage.}
 \item{\code{\link{clearAccessPolicy}}}{: Clear access policies for DataObjects in a DataPackage.}
 \item{\code{\link{containsId}}}{: Returns true if the specified object is a member of the data package.}
 \item{\code{\link{describeWorkflow}}}{: Add data derivation information to a DataPackage.}
 \item{\code{\link{getData}}}{: Get the data content of a specified data object.}
 \item{\code{\link{getSize}}}{: Get the Count of Objects in the DataPackage.}
 \item{\code{\link{getIdentifiers}}}{: Get the Identifiers of DataPackage members.}
 \item{\code{\link{getMember}}}{: Return the DataPackage Member by Identifier.}
 \item{\code{\link{getRelationships}}}{: Retrieve relationships of data package objects.}
 \item{\code{\link{getValue}}}{: Get values for selected DataPackage members.}
 \item{\code{\link{hasAccessRule}}}{: Determine if access rules exists for DataObjects in a DataPackage.}
 \item{\code{\link{insertRelationship}}}{: Insert relationships between objects in a DataPackage.}
 \item{\code{\link{removeAccessRule}}}{: Remove an access rule from DataObject in a DataPackage.}
 \item{\code{\link{removeMember}}}{: Remove the specified DataObject from a DataPackage.}
 \item{\code{\link{removeRelationships}}}{: Remove relationships of objects in a DataPackage.}
 \item{\code{\link{replaceMember}}}{: Replace the raw data or file associated with a DataObject.}
 \item{\code{\link{selectMember}}}{: Select package members based on slot values.}
 \item{\code{\link{serializePackage}}}{: Create an OAI-ORE resource map from the DataPackage.}
 \item{\code{\link{serializeToBagIt}}}{: Serialize A DataPackage into a BagIt Archive File.}
 \item{\code{\link{setPublicAccess}}}{: Set the access policy to readable by anyone for DataObject in a DataPackage.}
 \item{\code{\link{setValue}}}{: Set values for selected DataPackage members}
 \item{\code{\link{show}}}{: Print DataPackage information in a formatted view.}
 \item{\code{\link{updateMetadata}}}{: Update selected elements of the XML content of a DataObject in a DataPackage}
 \item{\code{\link{updateRelationships}}}{: Update package relationships by replacing an old identifier with a new one.}
}
}

\seealso{
\code{\link{datapack}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{calculateChecksum}
\alias{calculateChecksum}
\alias{calculateChecksum,DataObject-method}
\title{Calculate a checksum for the DataObject using the specified checksum algorithm}
\usage{
calculateChecksum(x, ...)

\S4method{calculateChecksum}{DataObject}(x, checksumAlgorithm = "SHA256", ...)
}
\arguments{
\item{x}{A DataObject instance}

\item{...}{Additional parameters (not yet used)}

\item{checksumAlgorithm}{a \code{character} value specifying the checksum algorithm to use (i.e "MD5" or "SHA1" or "SHA256")}
}
\value{
The calculated checksum
}
\description{
calculates a checksum
}
\note{
this method is intended for internal package use only.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{plotRelationships}
\alias{plotRelationships}
\alias{plotRelationships,DataPackage-method}
\title{Plot derivation relationships obtained from getRelationships}
\usage{
plotRelationships(x, ...)

\S4method{plotRelationships}{DataPackage}(x, col = NULL, ...)
}
\arguments{
\item{x}{a DataPackage object}

\item{...}{other options passed to the igraph plot function}

\item{col}{vector of colors used for plotting}
}
\description{
Creates graph of dataPackage object generated from getRelationships
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\name{initialize,SystemMetadata-method}
\alias{initialize,SystemMetadata-method}
\alias{SystemMetadata-initialize}
\title{Initialize a DataONE SystemMetadata object with default values or values passed in to the constructor.}
\usage{
\S4method{initialize}{SystemMetadata}(
  .Object,
  identifier = NA_character_,
  formatId = NA_character_,
  size = NA_real_,
  checksum = NA_character_,
  checksumAlgorithm = "SHA-256",
  submitter = NA_character_,
  rightsHolder = NA_character_,
  accessPolicy = data.frame(subject = character(), permission = character()),
  replicationAllowed = TRUE,
  numberReplicas = 3,
  obsoletes = NA_character_,
  obsoletedBy = NA_character_,
  archived = FALSE,
  dateUploaded = NA_character_,
  dateSysMetadataModified = NA_character_,
  originMemberNode = NA_character_,
  authoritativeMemberNode = NA_character_,
  preferredNodes = list(),
  blockedNodes = list(),
  seriesId = NA_character_,
  mediaType = NA_character_,
  fileName = NA_character_,
  mediaTypeProperty = list()
)
}
\arguments{
\item{.Object}{The object being initialized}

\item{identifier}{value of type \code{"character"}, the identifier of the object that this system metadata describes.}

\item{formatId}{value of type \code{"character"}, the DataONE object format for the object.}

\item{size}{value of type \code{"numeric"}, the size of the object in bytes.}

\item{checksum}{value of type \code{"character"}, the checksum for the object using the designated checksum algorithm.}

\item{checksumAlgorithm}{value of type \code{"character"}, the name of the hash function used to generate a checksum, from the DataONE controlled list.}

\item{submitter}{value of type \code{"character"}, the Distinguished Name or identifier of the person submitting the object.}

\item{rightsHolder}{value of type \code{"character"}, the Distinguished Name or identifier of the person who holds access rights to the object.}

\item{accessPolicy}{value of type \code{"data.frame"} containing (subject, permission) tuples to constitute the access authorization rules.}

\item{replicationAllowed}{value of type \code{"logical"}, for replication policy allows replicas.}

\item{numberReplicas}{value of type \code{"numeric"}, for number of supported replicas.}

\item{obsoletes}{value of type \code{"character"}, the identifier of an object which this object replaces.}

\item{obsoletedBy}{value of type \code{"character"}, the identifier of an object that replaces this object.}

\item{archived}{value of type \code{"logical"}, a boolean flag indicating whether the object has been archived and thus hidden.}

\item{dateUploaded}{value of type \code{"character"}, the date on which the object was uploaded to a member node.}

\item{dateSysMetadataModified}{value of type \code{"character"}, the last date on which this system metadata was modified.}

\item{originMemberNode}{value of type \code{"character"}, the node identifier of the node on which the object was originally registered.}

\item{authoritativeMemberNode}{value of type \code{"character"}, the node identifier of the node which currently is authoritative for the object.}

\item{preferredNodes}{list of \code{"character"}, each of which is the node identifier for a node to which a replica should be sent.}

\item{blockedNodes}{list of \code{"character"}, each of which is the node identifier for a node blocked from housing replicas.}

\item{seriesId}{value of type \code{"character"}, a unique Unicode string that identifies an object revision chain. A seriesId will resolve to the latest version of an object.}

\item{mediaType}{value of type \code{"character"}, the IANA Media Type (aka MIME-Type) of the object, e.g. "text/csv".}

\item{fileName}{value of type \code{"character"}, the name of the file to create when this object is downloaded from DataONE.}

\item{mediaTypeProperty}{value of type a \code{"list"} of \code{"character"}, IANA Media Type properties for the \code{"mediaType"} argument}
}
\value{
the SystemMetadata instance representing an object
}
\description{
Initialize a SystemMetadata object by providing default values for core information 
needed to manage objects across repository systems. SystemMetadata contains basic identification, ownership,
access policy, replication policy, and related metadata.
}
\seealso{
\url{https://releases.dataone.org/online/api-documentation-v2.0/apis/Types.html}

\code{\link{SystemMetadata-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{getSize}
\alias{getSize}
\alias{getSize,DataPackage-method}
\title{Get the Count of Objects in the Package}
\usage{
getSize(x, ...)

\S4method{getSize}{DataPackage}(x)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(not yet used)}
}
\value{
The number of object in the Package
}
\description{
Get the Count of Objects in the Package
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
getSize(dp)
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{updateMetadata}
\alias{updateMetadata}
\alias{updateMetadata,DataPackage-method}
\title{Update selected elements of the XML content of a DataObject in a DataPackage (aka package member).}
\usage{
updateMetadata(x, do, ...)

\S4method{updateMetadata}{DataPackage}(x, do, xpath, replacement, newId = NA_character_, ...)
}
\arguments{
\item{x}{a DataPackage instance}

\item{do}{A DataObject instance object, or DataObject identifier}

\item{...}{(Not yet used)}

\item{xpath}{A \code{character} value specifying the location in the XML to update.}

\item{replacement}{A \code{character} value that will replace the elements found with the \code{xpath}.}

\item{newId}{A value of type \code{"character"} which will replace the identifier for this DataObject.}
}
\description{
A DataObject that contains an XML document can be edited by specifying a path
to the elements to edit (an XPath expression) and a value to replace the text node.
}
\details{
This method requires some knowledge of the structure of the metadata document as well
as facility with the XPath language. If the \code{newId} argument is used, the specified new 
identifier will be assigned to the object, and the previous identifier will be stored in the \code{oldId} slot, 
for possible use when updating the DataObject to a repository. If \code{newId} is not used, a new
identifier will be generated for the DataObject only the first time that updateMetadata is called for
a particular object in a DataPackage.
}
\examples{
# Create a DataObject and add it to the DataPackage
dp <- new("DataPackage")
sampleMeta <- system.file("./extdata/sample-eml.xml", package="datapack")
id <- "1234"
metaObj <- new("DataObject", id="1234", format="eml://ecoinformatics.org/eml-2.1.1", 
                file=sampleMeta)
dp <- addMember(dp, metaObj)

# In the metadata object, insert the newly assigned data 
xp <- sprintf("//dataTable/physical/distribution[../objectName/text()=\"\%s\"]/online/url", 
              "sample-data.csv") 
newURL <- sprintf("https://cn.dataone.org/cn/v2/resolve/\%s", "1234")
dp <- updateMetadata(dp, id, xpath=xp, replacement=newURL)
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\docType{class}
\name{DataObject-class}
\alias{DataObject-class}
\title{DataObject wraps raw data with system-level metadata}
\description{
DataObject is a wrapper class that associates raw data or a data file with system-level metadata 
describing the data.  The system metadata includes attributes such as the object's identifier, 
type, size, checksum, owner, version relationship to other objects, access rules, and other critical metadata.
The SystemMetadata is compliant with the DataONE federated repository network's definition of SystemMetadata, and
is encapsulated as a separate object of type \code{\link{SystemMetadata}} that can be manipulated as needed. Additional science-level and
domain-specific metadata is out-of-scope for SystemMetadata, which is intended only for critical metadata for
managing objects in a repository system.
}
\details{
A DataObject can be constructed by passing the data and SystemMetadata to the new() method, or by passing
an identifier, data, format, user, and DataONE node identifier, in which case a SystemMetadata instance will
be generated with these fields and others that are calculated (such as size and checksum).

Data are associated with the DataObject either by passing it as a \code{'raw'} value to the \code{'dataobj'}
parameter in the constructor, which is then stored in memory, or by passing a fully qualified file path to the 
data in the \code{'filename'} parameter, which is then stored on disk.  One of dataobj or filename is required.
Use the \code{'filename'} approach when data are too large to be managed effectively in memory.  Callers can
access the \code{'filename'} slot to get direct access to the file, or can call \code{'getData()'} to retrieve the
contents of the data or file as a raw value (but this will read all of the data into memory).
}
\section{Slots}{

\describe{
\item{\code{sysmeta}}{A value of type \code{"SystemMetadata"}, containing the metadata about the object}

\item{\code{data}}{A value of type \code{"raw"}, containing the data represented in this object}

\item{\code{filename}}{A character value that contains the fully-qualified path to the object data on disk}

\item{\code{dataURL}}{A character value for the URL used to load data into this DataObject}

\item{\code{updated}}{A hash containing logical values which indicate if system metadata or the data object have been updated since object creation.}

\item{\code{oldId}}{A character string containing the previous identifier used, before a \code{"replaceMember"} call.}

\item{\code{targetPath}}{A character string holding the path of where the file should be in an exported package}
}}

\section{Methods}{

\itemize{
  \item{\code{\link[=DataObject-initialize]{initialize}}}{: Initialize a DataObject}
  \item{\code{\link{addAccessRule}}}{: Add a Rule to the AccessPolicy}
  \item{\code{\link{canRead}}}{: Test whether the provided subject can read an object.}
  \item{\code{\link{getData}}}{: Get the data content of a specified data object}
  \item{\code{\link{getFormatId}}}{: Get the FormatId of the DataObject}
  \item{\code{\link{getIdentifier}}}{: Get the Identifier of the DataObject}
  \item{\code{\link{hasAccessRule}}}{: Determine if an access rules exists for a DataObject.}
  \item{\code{\link{setPublicAccess}}}{: Add a Rule to the AccessPolicy to make the object publicly readable.}
  \item{\code{\link{updateXML}}}{: Update selected elements of the xml content of a DataObject}
}
}

\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
targetPath <- "myData/time-trials/trial_data.csv"
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB", targetPath=targetPath)
getIdentifier(do)
getFormatId(do)
getData(do)
canRead(do, "uid=anybody,DC=example,DC=com")
do <- setPublicAccess(do)
canRead(do, "public")
canRead(do, "uid=anybody,DC=example,DC=com")
# Also can create using a file for storage, rather than memory
\dontrun{
tf <- tempfile()
con <- file(tf, "wb")
writeBin(data, con)
close(con)
targetPath <- "myData/time-trials/trial_data.csv"
do <- new("DataObject", "id1", format="text/csv", user="uid=jones,DC=example,DC=com", 
  mnNodeId="urn:node:KNB", filename=tf, targetPath=targetPath)
}
}
\seealso{
\code{\link{datapack}}
}
\keyword{classes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{serializeToBagIt}
\alias{serializeToBagIt}
\alias{serializeToBagIt,DataPackage-method}
\title{Serialize A DataPackage into a BagIt Archive File}
\usage{
serializeToBagIt(x, ...)

\S4method{serializeToBagIt}{DataPackage}(
  x,
  mapId = NA_character_,
  syntaxName = NA_character_,
  namespaces = data.frame(),
  mimeType = NA_character_,
  syntaxURI = NA_character_,
  resolveURI = NA_character_,
  creator = NA_character_,
  ...
)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{Additional arguments}

\item{mapId}{A unique identifier for the package resource map. If not specified, one will be 
automatically generated.}

\item{syntaxName}{The name of the syntax to use for the resource map serialization, defaults to "rdfxml"}

\item{namespaces}{An optional data frame containing one or more namespaces and their associated prefix for 
the resource map serialization.}

\item{mimeType}{The mimetype for the resource map serialization, defaults to "application/rdf+xml".}

\item{syntaxURI}{An optional string specifying the URI for the resource map serialization.}

\item{resolveURI}{A character string containing a URI to prepend to datapackage identifiers for the resource map.}

\item{creator}{A \code{character} string containing the creator of the package.}
}
\value{
The file name that contains the BagIt zip archive.
}
\description{
The BagIt packaging format \url{https://tools.ietf.org/html/draft-kunze-bagit-08}
is used to prepare an archive file that contains the contents of a DataPackage.
}
\details{
A BagIt Archive File is created by copying each member of a DataPackage, and preparing
files that describe the files in the archive, including information about the size of the files
and a checksum for each file. An OAI-ORE resource map is automatically created and added to the
archive. These metadata files and the data files are then packaged into
a single zip file.
}
\examples{
# Create the first data object
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", id="do1", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
# Create a second data object
data2 <- charToRaw("7,8,9\n4,10,11")
do2 <- new("DataObject", id="do2", dataobj=data2, format="text/csv", user="jsmith")
dp <- addMember(dp, do2)
# Create a relationship between the two data objects
dp <- describeWorkflow(dp, sources="do2", derivations="do2")
# Write out the data package to a BagIt file
\dontrun{
bagitFile <- serializeToBagIt(dp, syntaxName="json", mimeType="application/json")
}
}
\seealso{
\code{\link{DataPackage-class}}

For more information and examples regarding the parameters specifying the creation of the resource map,
see \link{serializePackage}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\docType{class}
\name{SystemMetadata-class}
\alias{SystemMetadata-class}
\title{A DataONE SystemMetadata object containing basic identification, ownership, access policy, replication policy, and related metadata.}
\description{
A class representing DataONE SystemMetadata, which is core information about objects stored in a repository
and needed to manage those objects across systems.  SystemMetadata contains basic identification, ownership,
access policy, replication policy, and related metadata.
}
\section{Slots}{

\describe{
\item{\code{serialVersion}}{value of type \code{"numeric"}, the current version of this system metadata; only update the current version}

\item{\code{identifier}}{value of type \code{"character"}, the identifier of the object that this system metadata describes.}

\item{\code{replicationAllowed}}{value of type \code{"logical"}, replication policy allows replicas.}

\item{\code{numberReplicas}}{value of type \code{"numeric"}, for number of supported replicas.}

\item{\code{preferredNodes}}{value of type \code{"list"}, of preferred member nodes.}

\item{\code{blockedNodes}}{value of type \code{"list"}, of blocked member nodes.}

\item{\code{formatId}}{value of type \code{"character"}, the DataONE object format for the object.}

\item{\code{size}}{value of type \code{"numeric"}, the size of the object in bytes.}

\item{\code{checksum}}{value of type \code{"character"}, the checksum for the object using the designated checksum algorithm.}

\item{\code{checksumAlgorithm}}{value of type \code{"character"}, the name of the hash function used to generate a checksum, from the DataONE controlled list.}

\item{\code{submitter}}{value of type \code{"character"}, the Distinguished Name or identifier of the person submitting the object.}

\item{\code{rightsHolder}}{value of type \code{"character"}, the Distinguished Name or identifier of the person who holds access rights to the object.}

\item{\code{accessPolicy}}{value of type \code{"data.frame"}, a list of access rules as (subject, permission) tuples to be applied to the object.}

\item{\code{obsoletes}}{value of type \code{"character"}, the identifier of an object which this object replaces.}

\item{\code{obsoletedBy}}{value of type \code{"character"}, the identifier of an object that replaces this object.}

\item{\code{archived}}{value of type \code{"logical"}, a boolean flag indicating whether the object has been archived and thus hidden.}

\item{\code{dateUploaded}}{value of type \code{"character"}, the date on which the object was uploaded to a member node.}

\item{\code{dateSysMetadataModified}}{value of type \code{"character"}, the last date on which this system metadata was modified.}

\item{\code{originMemberNode}}{value of type \code{"character"}, the node identifier of the node on which the object was originally registered.}

\item{\code{authoritativeMemberNode}}{value of type \code{"character"}, the node identifier of the node which currently is authoritative for the object.}

\item{\code{seriesId}}{value of type \code{"character"}, a unique Unicode string that identifies an object revision chain. A seriesId will resolve to the latest version of an object.}

\item{\code{mediaType}}{value of type \code{"character"}, the IANA Media Type (aka MIME-Type) of the object, e.g. "text/csv".}

\item{\code{fileName}}{value of type \code{"character"}, the name of the file to create when this object is downloaded from DataONE.}

\item{\code{mediaTypeProperty}}{value of type a \code{"list"} of \code{"character"}, IANA Media Type properties for the \code{"mediaType"} argument}
}}

\section{Methods}{

\itemize{
 \item{\code{\link[=SystemMetadata-initialize]{initialize}}}{: Initialize a DataONE SystemMetadata object with default values or values passed in to the constructor object}
 \item{\code{\link{SystemMetadata}}}{: Create a SystemMetadata object, with all fields set to the value found in an XML document}
 \item{\code{\link{parseSystemMetadata}}}{: Parse an external XML document and populate a SystemMetadata object with the parsed data}
 \item{\code{\link{serializeSystemMetadata}}}{: Get the Count of Objects in the Package}
 \item{\code{\link{validate}}}{: Validate a SystemMetadata object}
 \item{\code{\link{addAccessRule}}}{: Add access rules to an object such as system metadata}
 \item{\code{\link{hasAccessRule}}}{: Determine if a particular access rules exists within SystemMetadata.}
 \item{\code{\link{clearAccessPolicy}}}{: Clear the accessPolicy from the specified object.}
}
}

\seealso{
\code{\link{datapack}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R, R/DataObject.R,
%   R/DataPackage.R
\name{hasAccessRule}
\alias{hasAccessRule}
\alias{hasAccessRule,SystemMetadata-method}
\alias{hasAccessRule,DataObject-method}
\alias{hasAccessRule,DataPackage-method}
\title{Determine if an access rules exists}
\usage{
hasAccessRule(x, ...)

\S4method{hasAccessRule}{SystemMetadata}(x, subject, permission)

\S4method{hasAccessRule}{DataObject}(x, subject, permission)

\S4method{hasAccessRule}{DataPackage}(x, subject, permission, identifiers = list(), ...)
}
\arguments{
\item{x}{the object to check for presence of the access rule.}

\item{...}{Additional arguments}

\item{subject}{of the rule to be checked}

\item{permission}{the permission to be checked}

\item{identifiers}{A list of \code{character} values containing package member identifiers for which the access rule will be checked.}
}
\value{
A logical value: if TRUE the access rule was found, if FALSE it was not found.

When called for SystemMetadata, boolean TRUE if the access rule exists already, FALSE otherwise

When called for a DataObject, boolean TRUE if the access rule exists already, FALSE otherwise

When called for a DataPackage, boolean TRUE if the access rule exists in all specified package members already, FALSE otherwise
}
\description{
Each SystemMetadata document may contain a set of (subject, permission) tuples
that represent the access rules for its associated object. This method determines
whether a particular access rule already exists within the set.

If called for a DataObject, then the SystemMetadata for the DataObject is checked.

If called for a DataPackage, then the SystemMetadata for DataObjects in the DataPackage are checked.
}
\examples{
#
# Check access rules for a SystemMetadata object.
sysmeta <- new("SystemMetadata")
sysmeta <- addAccessRule(sysmeta, "uid=smith,ou=Account,dc=example,dc=com", "write")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
  "uid=slaughter,o=unaffiliated,dc=example,dc=org"), permission=c("write", "changePermission"))
sysmeta <- addAccessRule(sysmeta, accessRules)
ruleExists <- hasAccessRule(sysmeta, subject="uid=smith,ou=Account,dc=example,dc=com", 
  permission="write")
#
# Check access rules for a DataObject
data <- system.file("extdata/sample-data.csv", package="datapack")
do <- new("DataObject", file=system.file("./extdata/sample-data.csv", package="datapack"), 
                                         format="text/csv")
do <- setPublicAccess(do)
isPublic <- hasAccessRule(do, "public", "read")
accessRules <- data.frame(subject=c("uid=smith,ou=Account,dc=example,dc=com", 
                          "uid=wiggens,o=unaffiliated,dc=example,dc=org"), 
                          permission=c("write", "changePermission"), 
                          stringsAsFactors=FALSE)
do <- addAccessRule(do, accessRules)
SmithHasWrite <- hasAccessRule(do, "uid=smith,ou=Account,dc=example,dc=com", "write")
#
# Check access rules for member DataObjects of a DataPackage.
# First create an example DataPackage
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="id1", dataobj=data, format="text/csv")
dp <- addMember(dp, obj)
data2 <- charToRaw("7,8,9\n4,10,11\n")
obj2 <- new("DataObject", id="id2", dataobj=data2, format="text/csv")
dp <- addMember(dp, obj2)
# Add access rules to all package members
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "write")
dp <- addAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")
hasWrite <- hasAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "write")
hasChange <- hasAccessRule(dp, "uid=smith,ou=Account,dc=example,dc=com", "changePermission")
}
\seealso{
\code{\link{SystemMetadata-class}}

\code{\link{DataObject-class}}

\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{getMember}
\alias{getMember}
\alias{getMember,DataPackage-method}
\title{Return the Package Member by Identifier}
\usage{
getMember(x, ...)

\S4method{getMember}{DataPackage}(x, identifier)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(Not yet used)}

\item{identifier}{A DataObject identifier}
}
\value{
A DataObject if the member is found, or NULL if not
}
\description{
Given the identifier of a member of the data package, return the DataObject
representation of the member.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", id="myNewId", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
do2 <- getMember(dp, "myNewId")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R
\name{canRead}
\alias{canRead}
\alias{canRead,DataObject-method}
\title{Test whether the provided subject can read an object.}
\usage{
canRead(x, ...)

\S4method{canRead}{DataObject}(x, subject)
}
\arguments{
\item{x}{DataObject}

\item{...}{Additional arguments}

\item{subject}{: the subject name of the person/system to check for read permissions}
}
\value{
boolean TRUE if the subject has read permission, or FALSE otherwise
}
\description{
Using the AccessPolicy, tests whether the subject has read permission
for the object.  This method is meant work prior to submission to a repository, 
and will show the permissions that would be enforced by the repository on submission.
Currently it only uses the AccessPolicy to determine who can read (and not the rightsHolder field,
which always can read an object).  If an object has been granted read access by the
special "public" subject, then all subjects have read access.
}
\details{
The subject name used in both the AccessPolicy and in the \code{'subject'}
argument to this method is a string value, but is generally formatted as an X.509
name formatted according to RFC 2253.
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="1234", dataobj=data, format="text/csv")
obj <- addAccessRule(obj, "smith", "read")
access <- canRead(obj, "smith")
}
\seealso{
\code{\link{DataObject-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataObject.R, R/DataPackage.R
\name{setPublicAccess}
\alias{setPublicAccess}
\alias{setPublicAccess,DataObject-method}
\alias{setPublicAccess,DataPackage-method}
\title{Add a Rule to the AccessPolicy to make the object publicly readable.}
\usage{
setPublicAccess(x, ...)

\S4method{setPublicAccess}{DataObject}(x)

\S4method{setPublicAccess}{DataPackage}(x, identifiers = list())
}
\arguments{
\item{x}{DataObject}

\item{...}{(not yet used)}

\item{identifiers}{A list of \code{character} values containing package member identifiers that will be updated (default is all package members).}
}
\value{
A DataObject with modified access rules.

A DataPackage with modified access rules.
}
\description{
To be called prior to creating the object in DataONE.  When called before 
creating the object, adds a rule to the access policy that makes this object
publicly readable.  If called after creation, it will only change the system
metadata locally, and will not have any effect on remotely uploaded copies of
the DataObject.
}
\examples{
data <- charToRaw("1,2,3\n4,5,6\n")
do <- new("DataObject", "id1", dataobj=data, "text/csv", 
  "uid=jones,DC=example,DC=com", "urn:node:KNB")
do <- setPublicAccess(do)
# First create a sample package with two DataObjects
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6\n")
obj <- new("DataObject", id="id1", dataobj=data, format="text/csv")
dp <- addMember(dp, obj)
data2 <- charToRaw("7,8,9\n4,10,11\n")
obj2 <- new("DataObject", id="id2", dataobj=data2, format="text/csv")
dp <- addMember(dp, obj2)
# Now add public read to all package members ("id1", "id2")
dp <- setPublicAccess(dp)
}
\seealso{
\code{\link{DataObject-class}}

\code{\link{DataObject-class}}

\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\docType{class}
\name{ResourceMap-class}
\alias{ResourceMap-class}
\title{ResourceMap provides methods to create, serialize and deserialize an OAI ORE resource map.}
\description{
The Open Archives Initiative Object Reuse and Exchange (OAI-ORE) defines standards for the description
and exchange of aggregations of web resources, such as a DataPackage. A Resource Map describes the objects
in a DataPackage and the relationships between these objects.
}
\section{Slots}{

\describe{
\item{\code{relations}}{value of type \code{"data.frame"}, containing RDF triples representing the relationship between package objects}

\item{\code{world}}{a Redland RDF World object}

\item{\code{storage}}{a Redland RDF Storage object}

\item{\code{model}}{a Redland RDF Model object}

\item{\code{id}}{a unique identifier for a ResourceMap instance}
}}

\section{Methods}{

\itemize{
 \item{\code{\link[=ResourceMap-initialize]{initialize}}}{: Initialize a ResourceMap object.}
 \item{\code{\link{createFromTriples}}}{: Populate a ResourceMap with RDF relationships from data.frame.}
 \item{\code{\link{getTriples}}}{: Get the RDF relationships stored in the ResourceMap.}
 \item{\code{\link{parseRDF}}}{: Parse an RDF/XML resource map from a file.}
 \item{\code{\link{serializeRDF}}}{: Write the ResourceMap relationships to a file.}
}
}

\examples{
dp <- new("DataPackage")
dp <- insertRelationship(dp, "/Users/smith/scripts/genFields.R",
    "http://www.w3.org/ns/prov#used",
    "https://knb.ecoinformatics.org/knb/d1/mn/v1/object/doi:1234/_030MXTI009R00_20030812.40.1")
relations <- getRelationships(dp)
resMap <- new("ResourceMap")
resMap <- createFromTriples(resMap, relations, getIdentifiers(dp))
\dontrun{
tf <- tempfile(fileext=".rdf")
serializeRDF(resMap, file=tf)
}
}
\seealso{
\code{\link{datapack}}
}
\keyword{resourceMap}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResourceMap.R
\name{getTriples}
\alias{getTriples}
\alias{getTriples,ResourceMap-method}
\title{Get the RDF relationships stored in the ResourceMap.}
\usage{
getTriples(x, ...)

\S4method{getTriples}{ResourceMap}(x, filter = TRUE, identifiers = list(), ...)
}
\arguments{
\item{x}{ResourceMap}

\item{...}{Additional parameters (not yet implemented).}

\item{filter}{A \code{logical} value. If TRUE, then DataONE packaging relationships are omitted.}

\item{identifiers}{A list of \code{character} values of the identifiers of DataPackage members.}
}
\value{
x A data.frame containing the relationships from the ResourceMap
}
\description{
The \code{getTriples} method extracts the RDF relationships from a ResourceMap.
}
\details{
The \code{filter} argument causes DataONE packaging relationships to be removed. 
A description of these can be viewed at https://purl.dataone.org/architecture/design/DataPackage.html. 
The \code{identifiers} parameter can contain a list of DataPackage members for which the 
identifiers will be 'demoted', that is any relationship that has these identifiers as a 
URL as the subject or object will be changed to the 'bare' identifier. The intent of these two parameter is to
transform the DataPackage to a 'local' state, so that it can be more easily updated locally.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{selectMember}
\alias{selectMember}
\alias{selectMember,DataPackage-method}
\title{Return identifiers for objects that match search criteria}
\usage{
selectMember(x, ...)

\S4method{selectMember}{DataPackage}(x, name, value, as = "character")
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(Not yet used)}

\item{name}{The name of the DataObject slot to inspect, for example "sysmeta@formatId".}

\item{value}{A character or logical value to match. If specified as a character value, PERL style regular expressions can be used (see ?grepl).}

\item{as}{A character value to specify the return type, either "DataObject" or "character" (the default)}
}
\value{
A list of matching DataObjects or DataObject identifiers. The default is to return a list of 
DataObject identifiers.
}
\description{
Return DataObjects or DataObject identifiers that match search terms.
}
\details{
The \code{"selectMember"} method inspects the DataObject slot \code{"name"} for a match with \code{"value"}
for each DataObject in a DataPackage. Matching DataObjects are returned as a list containing either package member
identifiers (character) or the DataObjects themselves, depending on the value of the \code{as} parameter.
}
\examples{
#' library(datapack)
dp <- new("DataPackage")
# Add the script to the DataPackage
progFile <- system.file("./extdata/pkg-example/logit-regression-example.R", package="datapack")
# An 'id' parameter is not specified, so one will be generated automatically.
progObj <- new("DataObject", format="application/R", filename=progFile)
dp <- addMember(dp, progObj)

# Add a script input to the DataPackage
inFile <- system.file("./extdata/pkg-example/binary.csv", package="datapack") 
inObj <- new("DataObject", format="text/csv", filename=inFile)
dp <- addMember(dp, inObj)

# Add a script output to the DataPackage
outFile <- system.file("./extdata/pkg-example/gre-predicted.png", package="datapack")
outObj <- new("DataObject", format="image/png", file=outFile)
dp <- addMember(dp, outObj)

# Now determine the package member identifier for the R script
progIds  <- selectMember(dp, name="sysmeta@formatId", value="application/R", as="character")
inputId <- selectMember(dp, name="sysmeta@fileName", value="binary.csv")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{addData}
\alias{addData}
\alias{addData,DataPackage,DataObject-method}
\title{Add a DataObject to the DataPackage}
\usage{
addData(x, do, ...)

\S4method{addData}{DataPackage,DataObject}(x, do, mo = NA_character_)
}
\arguments{
\item{x}{A DataPackage instance}

\item{do}{A DataObject instance}

\item{...}{(Additional parameters)}

\item{mo}{A DataObject (containing metadata describing \code{"do"} ) to associate with the science object.}
}
\value{
the updated DataPackage object
}
\description{
The DataObject is added to the DataPackage.
}
\details{
The DataObject \code{"do"} is added to the DataPackage. If the optional \code{"mo"} parameter is specified, then it is 
assumed that the DataObject \code{"mo"} is a metadata
object that describes the science object \code{"do"} that is being added. The \code{addData} function will add a relationship
to the DataPackage resource map that indicates that the metadata object describes the science object using the 
Citation Typing Ontology (CITO).
Note: this method updates the passed-in DataPackage object.
\code{documents} and \code{isDocumentedBy} relationship.
}
\examples{
dpkg <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
metadata <- charToRaw("EML or other metadata document text goes here\n")
md <- new("DataObject", id="md1", dataobj=metadata, format="text/xml", user="smith", 
  mnNodeId="urn:node:KNB")
do <- new("DataObject", id="id1", dataobj=data, format="text/csv", user="smith", 
  mnNodeId="urn:node:KNB")
# Associate the metadata object with the science object. The 'mo' object will be added 
# to the package  automatically, since it hasn't been added yet.
# This method is now deprecated, so suppress warnings if desired. 
suppressWarnings(dpkg <- addData(dpkg, do, md))
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{insertRelationship}
\alias{insertRelationship}
\alias{insertRelationship,DataPackage-method}
\title{Record relationships of objects in a DataPackage}
\usage{
insertRelationship(x, ...)

\S4method{insertRelationship}{DataPackage}(
  x,
  subjectID,
  objectIDs,
  predicate = NA_character_,
  subjectType = NA_character_,
  objectTypes = NA_character_,
  dataTypeURIs = NA_character_
)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{(Additional parameters)}

\item{subjectID}{The identifier of the subject of the relationship}

\item{objectIDs}{A list of identifiers of the object of the relationships (a relationship is recorded for each objectID)}

\item{predicate}{The IRI of the predicate of the relationship}

\item{subjectType}{the type to assign the subject, values can be 'uri', 'blank'}

\item{objectTypes}{the types to assign the objects (cal be single value or list), each value can be 'uri', 'blank', or 'literal'}

\item{dataTypeURIs}{An RDF data type that specifies the type of the object}
}
\value{
the updated DataPackage object
}
\description{
Record a relationship of the form "subject -> predicate -> object", as defined by the Resource Description Framework (RDF), i.e.
an RDF triple.
}
\details{
For use with DataONE, a best practice is to specify the subject and predicate as DataONE persistent identifiers 
(https://mule1.dataone.org/ArchitectureDocs-current/design/PIDs.html). If the objects are not known to DataONE, then local identifiers can be
used, and these local identifiers may be promoted to DataONE PIDs when the package is uploaded to a DataONE member node.
The predicate is typically an RDF property (as a IRI) from a schema supported by DataONE, i.e. "http://www.w3.org/ns/prov#wasGeneratedBy"
If multiple values are specified for argument objectIDS, a relationship is created for each value in the list "objectIDs". IF a value
is not specified for subjectType or objectType, then NA is assigned. Note that if these relationships are fetched via the getRelationships()
function, and passed to the createFromTriples() function to initialize a ResourceMap object, the underlying redland package will assign
appropriate values for subjects and objects.
Note: This method updates the passed-in DataPackage object.
}
\examples{
dp <- new("DataPackage")
# Create a relationship
dp <- insertRelationship(dp, "/Users/smith/scripts/genFields.R",
    "https://knb.ecoinformatics.org/knb/d1/mn/v1/object/doi:1234/_030MXTI009R00_20030812.40.1",
    "http://www.w3.org/ns/prov#used")
# Create a relationshp with the subject as a blank node with an automatically assigned blank 
# node id
dp <- insertRelationship(dp, subjectID=NA_character_, objectIDs="thing6", 
    predicate="http://www.myns.org/wasThing")
# Create a relationshp with the subject as a blank node with a user assigned blank node id
dp <- insertRelationship(dp, subjectID="urn:uuid:bc9e160e-ca21-47d5-871b-4a4820fe4451", 
      objectIDs="thing7", predicate="http://www.myns.org/hadThing")
# Create multiple relationships with the same subject, predicate, but different objects
dp <- insertRelationship(dp, subjectID="urn:uuid:95055dc1-b2a0-4a00-bdc2-05c16d048ca2", 
      objectIDs=c("thing4", "thing5"), predicate="http://www.myns.org/hadThing")
# Create multiple relationships with subject and object types specified
dp <- insertRelationship(dp, subjectID="orcid.org/0000-0002-2192-403X", 
    objectIDs="http://www.example.com/home", predicate="http://www.example.com/hadHome",
                   subjectType="uri", objectType="literal")                
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{datapack-deprecated}
\alias{datapack-deprecated}
\title{Deprecated Methods}
\description{
The following items are deprecated in this release of datapack and will be
marked as Defunct and removed in a future version.
}
\section{These methods are deprecated}{

\itemize{
 \item{\code{\link{recordDerivation}}}{: Record derivation relationships between objects in a DataPackage.}
 \item{\code{\link{addData}}}{: Add a DataObject to the DataPackage.}
}
}

\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{getRelationships}
\alias{getRelationships}
\alias{getRelationships,DataPackage-method}
\title{Retrieve relationships of package objects}
\usage{
getRelationships(x, ...)

\S4method{getRelationships}{DataPackage}(x, condense = F, ...)
}
\arguments{
\item{x}{A DataPackage object}

\item{...}{(Not yet used)}

\item{condense}{A logical value, if TRUE then a more easily viewed version of relationships are returned.}
}
\description{
Relationships of objects in a package are defined using the \code{'insertRelationship'} call and retrieved
using \code{getRetaionships}. These relationships are returned in a data frame with \code{'subject'}, \code{'predicate'}, \code{'objects'}
as the columns, ordered by "subject"
}
\examples{
dp <- new("DataPackage")
insertRelationship(dp, "/Users/smith/scripts/genFields.R",
    "http://www.w3.org/ns/prov#used",
    "https://knb.ecoinformatics.org/knb/d1/mn/v1/object/doi:1234/_030MXTI009R00_20030812.40.1")
rels <- getRelationships(dp)
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{getValue}
\alias{getValue}
\alias{getValue,DataPackage-method}
\title{Get values for selected DataPackage members.}
\usage{
getValue(x, ...)

\S4method{getValue}{DataPackage}(x, name, identifiers = NA_character_)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(Not yet used)}

\item{name}{A name of a DataObject slot.}

\item{identifiers}{A list of DataPackage member identifiers}
}
\value{
A list of values for matching slot names and included identifiers.
}
\description{
Given a slot name and set of package member identifiers, return slot values.
}
\details{
If the parameter \code{identifiers} is provided, then only the DataPackage
members that have identifiers in the provided list will have there values fetched.
If this parameter is not provided, then the values for all DataPackage members are returned.
}
\examples{
dp <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
do <- new("DataObject", id="myNewId", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
data <- charToRaw("7,8.9\n4,10,11")
do <- new("DataObject", id="myNewId2", dataobj=data, format="text/csv", user="jsmith")
dp <- addMember(dp, do)
formats <- getValue(dp, name="sysmeta@formatId")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmsg.R
\name{dmsg}
\alias{dmsg}
\title{Print a debugging message to stderr.}
\usage{
dmsg(msg)
}
\arguments{
\item{msg}{the message to be printed}
}
\description{
Print a debugging message to stderr.
}
\details{
Only print the message if the option "datapack.debugging_mode" is TRUE.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SystemMetadata.R
\name{serializeSystemMetadata}
\alias{serializeSystemMetadata}
\alias{serializeSystemMetadata,SystemMetadata-method}
\title{Serialize a SystemMetadata object to an XML representation}
\usage{
serializeSystemMetadata(x, ...)

\S4method{serializeSystemMetadata}{SystemMetadata}(x, version = "v1", ...)
}
\arguments{
\item{x}{The SystemMetadata instance to be serialized.}

\item{...}{(Not currently used)}

\item{version}{A character string representing the DataONE API version that this system will be used with (e.g. "v1", "v2").}
}
\value{
A character value of the filename that the XML representation of the SystemMetadata object was written to.

the character string representing a SystemMetadata object
}
\description{
The SystemMetadata object is converted to XML and 
written to a file.
}
\details{
If the \code{'version'} parameter is specified as *v2* then the SystemMetadata
object is serialized according to the DataONE version 2.0 system metadata format.
}
\examples{
library(XML)
doc <- xmlParseDoc(system.file("testfiles/sysmeta.xml", package="datapack"), asText=FALSE)
sysmeta <- new("SystemMetadata")
sysmeta <- parseSystemMetadata(sysmeta, xmlRoot(doc))
sysmetaXML <- serializeSystemMetadata(sysmeta, version="v2")
}
\seealso{
\code{\link{SystemMetadata-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{replaceMember}
\alias{replaceMember}
\alias{replaceMember,DataPackage-method}
\title{Replace the raw data or file associated with a DataObject}
\usage{
replaceMember(x, do, ...)

\S4method{replaceMember}{DataPackage}(
  x,
  do,
  replacement,
  formatId = NA_character_,
  mediaType = NA_character_,
  mediaTypeProperty = NA_character_,
  newId = NA_character_,
  ...
)
}
\arguments{
\item{x}{A DataPackage instance}

\item{do}{A DataObject instance}

\item{...}{(Not yet used)}

\item{replacement}{A \code{raw} object or \code{character} (for filename) that will replace the current value in the DataObject \code{do}.}

\item{formatId}{A value of type \code{"character"}, the DataONE object format for the object.}

\item{mediaType}{A value of type \code{"character"}, the IANA Media Type (aka MIME-Type) of the object, e.g. "text/csv".}

\item{mediaTypeProperty}{A value of type \code{"list"} of \code{"character"}, IANA Media Type properties for the \code{"mediaType"} argument.}

\item{newId}{A value of type \code{"character"} which will replace the identifier for this DataObject.}
}
\description{
A DataObject is a container for data that can be either an R raw object or
a file on local disk. The \code{replaceMember} method can be used to update the
date that a DataObject contains, for a DataObject that is a member of a DataPackage, 
substituting a new file or raw object in the specified DataObject.
}
\details{
The data that is replacing the existing DataObject data may be of a different
format or type than the existing data. Because the data type and format may change, the
system metadata that describes the data can be updated as well. The \code{replaceMember}
method will update the SystemMetadata \code{size}, \code{checksum} values automatically, 
but does not update the \code{formatId}, \code{mediaType}, \code{mediaTypeProperty}
unless requested, so these should be specified in the call to \code{replaceMember} if necessary. 
If the \code{newId} argument is used, the specified new identifier will be assigned to the 
object, otherwise one will be generated if necessary. This new identifier will be used
if the DataPackage is uploaded to DataONE, and this object is updating an existing object in DataONE.
}
\examples{
# Create a DataObject and add it to the DataPackage
dp <- new("DataPackage")
doIn <- new("DataObject", format="text/csv", 
            filename=system.file("./extdata/pkg-example/binary.csv", package="datapack"))
dp <- addMember(dp, doIn)

# Use the zipped version of the file instead by updating the DataObject
dp <- replaceMember(dp, doIn, 
          replacement=system.file("./extdata/pkg-example/binary.csv.zip", 
          package="datapack"),
                    formatId="application/zip")
}
\seealso{
\code{\link{DataPackage-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataPackage.R
\name{addMember}
\alias{addMember}
\alias{addMember,DataPackage-method}
\title{Add a DataObject to the DataPackage}
\usage{
addMember(x, ...)

\S4method{addMember}{DataPackage}(x, do, mo = NA_character_)
}
\arguments{
\item{x}{A DataPackage instance}

\item{...}{(Additional parameters)}

\item{do}{The DataObject to add.}

\item{mo}{A DataObject (containing metadata describing \code{"do"} ) to associate with the science object. If this DataObject 
has already been added to the package, the argument can be a \code{"character"} containing the DataObject identifier.}
}
\value{
the updated DataPackage object
}
\description{
The DataObject is added to the DataPackage.
}
\details{
The DataObject \code{"do"} is added to the DataPackage. If the optional \code{"mo"} parameter is specified, then it is 
assumed that the DataObject \code{"mo"} is a metadata
object that describes the science object \code{"do"} that is being added. The \code{addMember} function will add a relationship
to the DataPackage resource map that indicates that the metadata object describes the science object using the 
Citation Typing Ontology (CITO).
Note: this method updates the passed-in DataPackage object.
\code{documents} and \code{isDocumentedBy} relationship.
}
\examples{
dpkg <- new("DataPackage")
data <- charToRaw("1,2,3\n4,5,6")
metadata <- charToRaw("EML or other metadata document text goes here\n")
md <- new("DataObject", id="md1", dataobj=metadata, format="text/xml", user="smith", 
  mnNodeId="urn:node:KNB")
do <- new("DataObject", id="id1", dataobj=data, format="text/csv", user="smith", 
  mnNodeId="urn:node:KNB")
# Associate the metadata object with the science object. The 'mo' object will be added 
# to the package  automatically, since it hasn't been added yet.
dpkg <- addMember(dpkg, do, md)
}
\seealso{
\code{\link{DataPackage-class}}
}
