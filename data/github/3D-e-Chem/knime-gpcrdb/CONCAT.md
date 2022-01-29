# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [1.4.2] - 2020-04-06

Nodes in existing workflows giving not found errors should be replaced.

### Fixed

- Not found error [#32](https://github.com/3D-e-Chem/knime-gpcrdb/issues/320

## [1.4.1] - 2019-06-27

### Changed

- Compatible with KNIME 4 [#30](https://github.com/3D-e-Chem/knime-gpcrdb/issues/30)

## [1.4.0] - 2018-04-04

### Fixed

- Add all fields that are returned by mutant ws endpoint to KNIME node (#29)

## [1.3.0] - 2017-09-01

### Fixed

- Mutant api changed (#29)

## [1.2.0] - 2017-03-01

### Added

- Node to fetch all GPCRDB structures (#27)

### Fixed

- Use rowkey generator instead of custom rowkey (#28)

## [1.1.0] - 2017-01-12

### Added

- Test coverage (#25)
- Timeout option (#21)
- Workflow test for all nodes (#24)

### Changed

- Nicer message for Timeout and not found exceptions (#23)

### Fixed

- Change default port names to specific names (#22)

## [1.0.15] - 2016-10-18

### Fixed

- Protein with single structure produce error (#8)

## [1.0.14] - 2016-10-18

### Added

- Test workflow using WireMock server (#17)

### Fixed

- Protein with single structure produce error (#8)

## [1.0.13] - 2016-09-13

### Changed

- Client no longer needs to be modified after generation, now the Swagger spec is adjusted before generation (#16)

### Fixed

- Based client on okhttp+gson library, a non jaxrs based implementation (#16)

## [1.0.12] - 2016-07-18

### Changed

- Moved nodes under /community/3d-e-chem (#15)

## [1.0.11] - 2016-06-07

### Added

- Additional input and node configuration checks in line with KNIME specification

## [1.0.10] - 2016-05-20

### Changed

- Shortened maven modules names

### Fixed

- Optional publication when requesting structures of protein (#12)

## [1.0.9] - 2016-05-13

### Fixed

- Protein with single structure produce error (#8)

### Changed

- Input is case-insensitive (#10)

## [1.0.8] - 2016-04-18

### Added

- Node for retrieving protein similarity (#5)

## [1.0.7] - 2016-03-29

### Added

- Node for retrieving mutations of protein (#4)
- Node for retrieving interactions of structure with its ligands (#6)

## [1.0.6] - 2016-01-27

### Changed

- Return alternative numbers as map (#7)
# GPCRDB node for KNIME

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-gpcrdb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-gpcrdb)
[![Build status](https://ci.appveyor.com/api/projects/status/4n4bjgaq04dbem0u?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/knime-gpcrdb)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/116701411bee4b92a9f265f1a0a9efaf)](https://www.codacy.com/app/3D-e-Chem/knime-gpcrdb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/knime-gpcrdb&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/116701411bee4b92a9f265f1a0a9efaf)](https://www.codacy.com/app/3D-e-Chem/knime-gpcrdb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/knime-gpcrdb&amp;utm_campaign=Badge_Coverage)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3257985.svg)](https://doi.org/10.5281/zenodo.3257985)

KNIME plugin for retrieving data from [https://gpcrdb.org](https://gpcrdb.org), GPCRdb website contains data, web tools and diagrams for G protein-coupled receptors (GPCRs).

## Installation

Requirements:

* KNIME, https://www.knime.org, version 4.0 or higher

Steps to get GPCRDB nodes inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/updates`
4. Select --all sites-- in work with pulldown
5. Open KNIME 3D-e-Chem Contributions folder
6. Select GPCRDB
7. Install software & restart

## Usage

See example workflow in `examples` folder.

## Build

Requires Maven and JDK 8.

```shell
mvn verify
```

Jar has been made in `plugin/target` folder.
An Eclipse update site will be made in `p2/target/repository` repository.

## Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.molviewer.targetplatform/KNIME-AP-4.0.target` target definition.

During import the Tycho Eclipse providers must be installed.

## New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Test node by installing it from local update site
5. Append new release to 3D-e-Chem update site
    1. Make clone of [https://github.com/3D-e-Chem/3D-e-Chem.github.io](https://github.com/3D-e-Chem/3D-e-Chem.github.io) repo
    2. Append release to 3D-e-Chem update site with `mvn install -Dtarget.update.site=<3D-e-Chem repo/updates>`
6. Commit and push changes in this repo and 3D-e-Chem.github.io repo
7. Create a Github release
8. Update Zenodo entry
    1. Fix authors
    2. Fix license
9. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release

## Create GPCRDB client

1. Download swagger code generator

  ```shell
  wget http://repo1.maven.org/maven2/io/swagger/swagger-codegen-cli/2.3.0/swagger-codegen-cli-2.3.0.jar
  ```

2. Download and unpack the swagger rewriter

3. Generate a Swagger spec for the client

  Install the swagger rewriter from https://github.com/3D-e-Chem/swagger-rewriter

  ```shell
  swagger-rewriter/bin/swagger-rewriter \
  https://gpcrdb.org/services/reference/api-docs/ \
  client-config/swagger-rewriter.config.yml \
  client-config/gpcrdb.swagger-spec.json
  ```

3.1 Optionally, make manual changes to client-config/gpcrdb.swagger-spec.json

4. Generate a client for GPCRDB web service using the rewritten spec

  ```shell
  java -jar swagger-codegen-cli-2.3.0.jar generate \
  --input-spec client-config/gpcrdb.swagger-spec.json \
  --output client \
  --lang java \
  --config client-config/swagger-codegen.config.json
  ```

5. Compile client

  ```shell
  cd client
  mvn package
  ```

6. Make client jar and it's dependencies available in plugin

  ```shell
  cp -r target/lib/* target/*jar ../plugin/lib/
  ```

7. Remove test dependencies

  ```shell
  rm plugin/lib/*-tests.jar plugin/lib/junit* plugin/lib/hamcrest*
  ```

8. Update `plugin/META-INF/MANIFEST.MF`, `plugin/build.properties` files to reflect contents of lib/

## Create stub recordings for integration tests

The test workflow are tested against a mocked web server and not the actual https://gpcrdb.org site.
The mock server is called [WireMock](http://WireMock.org/) and normally gives empty responses.
To have WireMock server return filled responses, stubs stored in `tests/src/test/resources/` directory must be provided.
The stubs can be recorded by starting a WireMock server in recording mode by:

```shell
java -jar tests/lib/wiremock-standalone-2.5.0.jar --proxy-all="https://gpcrdb.org/" \
--port=8089 --record-mappings --verbose --root-dir=tests/src/test/resources/
```

Then in a KNIME workflow in the GPCRDB nodes set the base path to http://localhost:8089.
Executing the workflow will fetch data from https://gpcrdb.org/ via the WireMock server and cause new stubs to be recorded in the `tests/src/test/resources/` directory.

To run the test workflows from inside KNIME desktop environment start the WireMock server in mock mode by:

```shell
java -jar tests/lib/wiremock-standalone-2.5.0.jar --port=8089 --verbose --root-dir=tests/src/test/resources/
```

Then import the test workflows in `tests/src/knime/` directory, select the workflow in the KNIME explorer and in the context menu (right-click) select `Run as workflow test`.

## References

* V Isberg, S Mordalski, C Munk, K Rataj, K Harpsøe, AS Hauser, B Vroling, AJ Bojarski, G Vriend, DE Gloriam. “GPCRdb: an information system for G protein-coupled receptors”, 2016, Nucleic Acids Res., 44, D356-D364. [10.1093/nar/gkv1178](http://dx.doi.org/10.1093/nar/gkv1178)
* V Isberg, B Vroling, R van der Kant, K Li, G Vriend* and DE Gloriam*, “GPCRDB: an information system for G protein-coupled receptors”, 2014, Nucleic Acids Res., 42 (D1), D422-D425. [10.1093/nar/gkv1178](http://dx.doi.org/10.1093/nar/gkv1178)
# nl.esciencecenter.e3dchem.gpcrdb.client

## Requirements

Building the API client library requires [Maven](https://maven.apache.org/) to be installed.

## Installation

To install the API client library to your local Maven repository, simply execute:

```shell
mvn install
```

To deploy it to a remote Maven repository instead, configure the settings of the repository and execute:

```shell
mvn deploy
```

Refer to the [official documentation](https://maven.apache.org/plugins/maven-deploy-plugin/usage.html) for more information.

### Maven users

Add this dependency to your project's POM:

```xml
<dependency>
    <groupId>nl.esciencecenter.e3dchem.gpcrdb</groupId>
    <artifactId>nl.esciencecenter.e3dchem.gpcrdb.client</artifactId>
    <version>2.0.3</version>
    <scope>compile</scope>
</dependency>
```

### Gradle users

Add this dependency to your project's build file:

```groovy
compile "nl.esciencecenter.e3dchem.gpcrdb:nl.esciencecenter.e3dchem.gpcrdb.client:2.0.3"
```

### Others

At first generate the JAR by executing:

    mvn package

Then manually install the following JARs:

* target/nl.esciencecenter.e3dchem.gpcrdb.client-2.0.3.jar
* target/lib/*.jar

## Getting Started

Please follow the [installation](#installation) instruction and execute the following Java code:

```java

import nl.esciencecenter.e3dchem.gpcrdb.client.*;
import nl.esciencecenter.e3dchem.gpcrdb.client.auth.*;
import nl.esciencecenter.e3dchem.gpcrdb.client.model.*;
import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;

import java.io.File;
import java.util.*;

public class ServicesalignmentApiExample {

    public static void main(String[] args) {
        
        ServicesalignmentApi apiInstance = new ServicesalignmentApi();
        String slug = "slug_example"; // String | 
        try {
            Object result = apiInstance.familyAlignmentGET(slug);
            System.out.println(result);
        } catch (ApiException e) {
            System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentGET");
            e.printStackTrace();
        }
    }
}

```

## Documentation for API Endpoints

All URIs are relative to *https://gpcrdb.org*

Class | Method | HTTP request | Description
------------ | ------------- | ------------- | -------------
*ServicesalignmentApi* | [**familyAlignmentGET**](docs/ServicesalignmentApi.md#familyAlignmentGET) | **GET** /services/alignment/family/{slug}/ | Get a full sequence alignment of a protein family including a consensus sequence
*ServicesalignmentApi* | [**familyAlignmentGET_0**](docs/ServicesalignmentApi.md#familyAlignmentGET_0) | **GET** /services/alignment/family/{slug}/statistics/ | Get a full sequence alignment of a protein family including a consensus sequence
*ServicesalignmentApi* | [**familyAlignmentPartialGET**](docs/ServicesalignmentApi.md#familyAlignmentPartialGET) | **GET** /services/alignment/family/{slug}/{segments}/ | Get a partial sequence alignment of a protein family
*ServicesalignmentApi* | [**familyAlignmentPartialGET_0**](docs/ServicesalignmentApi.md#familyAlignmentPartialGET_0) | **GET** /services/alignment/family/{slug}/{segments}/statistics/ | Get a partial sequence alignment of a protein family
*ServicesalignmentApi* | [**familyAlignmentPartialSpeciesGET**](docs/ServicesalignmentApi.md#familyAlignmentPartialSpeciesGET) | **GET** /services/alignment/family/{slug}/{segments}/{latin_name}/ | Get a partial sequence alignment of a protein family
*ServicesalignmentApi* | [**familyAlignmentPartialSpeciesGET_0**](docs/ServicesalignmentApi.md#familyAlignmentPartialSpeciesGET_0) | **GET** /services/alignment/family/{slug}/{segments}/{latin_name}/statistics/ | Get a partial sequence alignment of a protein family
*ServicesalignmentApi* | [**proteinAlignmentGET**](docs/ServicesalignmentApi.md#proteinAlignmentGET) | **GET** /services/alignment/protein/{proteins}/ | Get a full sequence alignment of two or more proteins
*ServicesalignmentApi* | [**proteinAlignmentPartialGET**](docs/ServicesalignmentApi.md#proteinAlignmentPartialGET) | **GET** /services/alignment/protein/{proteins}/{segments}/ | Get a partial sequence alignment of two or more proteins
*ServicesalignmentApi* | [**proteinAlignmentStatisticsGET**](docs/ServicesalignmentApi.md#proteinAlignmentStatisticsGET) | **GET** /services/alignment/protein/{proteins}/statistics/ | Add a /statics at the end of an alignment in order to
*ServicesalignmentApi* | [**proteinAlignmentStatisticsGET_0**](docs/ServicesalignmentApi.md#proteinAlignmentStatisticsGET_0) | **GET** /services/alignment/protein/{proteins}/{segments}/statistics/ | Add a /statics at the end of an alignment in order to
*ServicesalignmentApi* | [**proteinSimilaritySearchAlignmentGET**](docs/ServicesalignmentApi.md#proteinSimilaritySearchAlignmentGET) | **GET** /services/alignment/similarity/{proteins}/{segments}/ | Get a segment sequence alignment of two or more proteins ranked by similarity
*ServicesmutantsApi* | [**mutantListGET**](docs/ServicesmutantsApi.md#mutantListGET) | **GET** /services/mutants/{entry_name}/ | Get a list of mutants of single protein instance by entry name
*ServicesplotApi* | [**helixBoxGET**](docs/ServicesplotApi.md#helixBoxGET) | **GET** /services/plot/helixbox/{entry_name}/ | Get SVG source code for a protein&#39;s helix box plot
*ServicesplotApi* | [**snakePlotGET**](docs/ServicesplotApi.md#snakePlotGET) | **GET** /services/plot/snake/{entry_name}/ | Get SVG source code for a protein&#39;s snake plot
*ServicesproteinApi* | [**proteinByAccessionDetailGET**](docs/ServicesproteinApi.md#proteinByAccessionDetailGET) | **GET** /services/protein/accession/{accession}/ | Get a single protein instance by accession
*ServicesproteinApi* | [**proteinDetailGET**](docs/ServicesproteinApi.md#proteinDetailGET) | **GET** /services/protein/{entry_name}/ | Get a single protein instance by entry name
*ServicesproteinfamilyApi* | [**proteinFamilyChildrenListGET**](docs/ServicesproteinfamilyApi.md#proteinFamilyChildrenListGET) | **GET** /services/proteinfamily/children/{slug}/ | Get a list of child families of a protein family
*ServicesproteinfamilyApi* | [**proteinFamilyDescendantListGET**](docs/ServicesproteinfamilyApi.md#proteinFamilyDescendantListGET) | **GET** /services/proteinfamily/descendants/{slug}/ | Get a list of descendant families of a protein family
*ServicesproteinfamilyApi* | [**proteinFamilyDetailGET**](docs/ServicesproteinfamilyApi.md#proteinFamilyDetailGET) | **GET** /services/proteinfamily/{slug}/ | Get a single protein family instance
*ServicesproteinfamilyApi* | [**proteinFamilyListGET**](docs/ServicesproteinfamilyApi.md#proteinFamilyListGET) | **GET** /services/proteinfamily/ | Get a list of protein families
*ServicesproteinfamilyApi* | [**proteinsInFamilyListGET**](docs/ServicesproteinfamilyApi.md#proteinsInFamilyListGET) | **GET** /services/proteinfamily/proteins/{slug}/ | Get a list of proteins in a protein family
*ServicesproteinfamilyApi* | [**proteinsInFamilySpeciesListGET**](docs/ServicesproteinfamilyApi.md#proteinsInFamilySpeciesListGET) | **GET** /services/proteinfamily/proteins/{slug}/{latin_name}/ | Get a list of proteins in a protein family
*ServicesresiduesApi* | [**residuesExtendedListGET**](docs/ServicesresiduesApi.md#residuesExtendedListGET) | **GET** /services/residues/extended/{entry_name}/ | Get a list of residues of a protein, including alternative generic numbers
*ServicesresiduesApi* | [**residuesListGET**](docs/ServicesresiduesApi.md#residuesListGET) | **GET** /services/residues/{entry_name}/ | Get a list of residues of a protein
*ServicesspeciesApi* | [**speciesDetailGET**](docs/ServicesspeciesApi.md#speciesDetailGET) | **GET** /services/species/{latin_name}/ | Get a single species instance
*ServicesspeciesApi* | [**speciesListGET**](docs/ServicesspeciesApi.md#speciesListGET) | **GET** /services/species/ | Get a list of species
*ServicesstructureApi* | [**representativeStructureListGET**](docs/ServicesstructureApi.md#representativeStructureListGET) | **GET** /services/structure/representative/ | Get a list of representative structures (one for each protein and activation state)
*ServicesstructureApi* | [**representativeStructureListProteinGET**](docs/ServicesstructureApi.md#representativeStructureListProteinGET) | **GET** /services/structure/protein/{entry_name}/representative/ | Get a list of representative structures of a protein (one for each activation state)
*ServicesstructureApi* | [**structureAssignGenericNumbersPOST**](docs/ServicesstructureApi.md#structureAssignGenericNumbersPOST) | **POST** /services/structure/assign_generic_numbers | Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file
*ServicesstructureApi* | [**structureDetailGET**](docs/ServicesstructureApi.md#structureDetailGET) | **GET** /services/structure/{pdb_code}/ | Get a single structure instance
*ServicesstructureApi* | [**structureLigandInteractionsGET**](docs/ServicesstructureApi.md#structureLigandInteractionsGET) | **GET** /services/structure/{pdb_code}/interaction/ | Get a list of interactions between structure and ligand
*ServicesstructureApi* | [**structureListGET**](docs/ServicesstructureApi.md#structureListGET) | **GET** /services/structure/ | Get a list of structures
*ServicesstructureApi* | [**structureListProteinGET**](docs/ServicesstructureApi.md#structureListProteinGET) | **GET** /services/structure/protein/{entry_name}/ | Get a list of structures of a protein
*ServicesstructureApi* | [**structureSingleProteinGET**](docs/ServicesstructureApi.md#structureSingleProteinGET) | **GET** /services/structure/protein/{entry_name}/?single | Get the structure of a protein
*ServicesstructureApi* | [**structureTemplateGET**](docs/ServicesstructureApi.md#structureTemplateGET) | **GET** /services/structure/template/{entry_name}/ | Get the most similar structure template for a protein using a 7TM alignment
*ServicesstructureApi* | [**structureTemplatePartialGET**](docs/ServicesstructureApi.md#structureTemplatePartialGET) | **GET** /services/structure/template/{entry_name}/{segments}/ | Get the most similar structure template for a protein using a partial alignment


## Documentation for Models

 - [Ligand](docs/Ligand.md)
 - [MutationSerializer](docs/MutationSerializer.md)
 - [ParentProteinFamilySerializer](docs/ParentProteinFamilySerializer.md)
 - [ProteinFamilySerializer](docs/ProteinFamilySerializer.md)
 - [ProteinSerializer](docs/ProteinSerializer.md)
 - [ProteinSimilarities](docs/ProteinSimilarities.md)
 - [ProteinSimilarity](docs/ProteinSimilarity.md)
 - [ResidueExtendedSerializer](docs/ResidueExtendedSerializer.md)
 - [ResidueGenericNumberSerializer](docs/ResidueGenericNumberSerializer.md)
 - [ResidueSerializer](docs/ResidueSerializer.md)
 - [SpeciesSerializer](docs/SpeciesSerializer.md)
 - [Structure](docs/Structure.md)
 - [StructureLigandInteractionSerializer](docs/StructureLigandInteractionSerializer.md)
 - [WriteMutationSerializer](docs/WriteMutationSerializer.md)
 - [WriteParentProteinFamilySerializer](docs/WriteParentProteinFamilySerializer.md)
 - [WriteProteinFamilySerializer](docs/WriteProteinFamilySerializer.md)
 - [WriteProteinSerializer](docs/WriteProteinSerializer.md)
 - [WriteResidueExtendedSerializer](docs/WriteResidueExtendedSerializer.md)
 - [WriteResidueGenericNumberSerializer](docs/WriteResidueGenericNumberSerializer.md)
 - [WriteResidueSerializer](docs/WriteResidueSerializer.md)
 - [WriteSpeciesSerializer](docs/WriteSpeciesSerializer.md)
 - [WriteStructureLigandInteractionSerializer](docs/WriteStructureLigandInteractionSerializer.md)


## Documentation for Authorization

All endpoints do not require authorization.
Authentication schemes defined for the API:

## Recommendation

It's recommended to create an instance of `ApiClient` per thread in a multithreaded environment to avoid any potential issues.

## Author




# ResidueSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequenceNumber** | **Integer** |  | 
**aminoAcid** | **String** |  | 
**proteinSegment** | **String** |  | 
**displayGenericNumber** | **String** |  | 




# WriteMutationSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**reference** | **String** |  | 
**protein** | **String** |  | 
**mutationPos** | **Integer** |  | 
**mutationFrom** | **String** |  | 
**mutationTo** | **String** |  | 
**ligandName** | **String** |  | 
**ligandIdtype** | **String** |  | 
**ligandId** | **String** |  | 
**ligandClass** | **String** |  | 
**expType** | **String** |  | 
**expFunc** | **String** |  | 
**expWtValue** | **Float** |  | 
**expWtUnit** | **String** |  | 
**expMuEffectSign** | **String** |  | 
**expMuEffectType** | **String** |  | 
**expMuEffectValue** | **Float** |  | 
**expFoldChange** | **Float** |  | 
**expMuEffectQual** | **String** |  | 
**expMuEffectLigandProp** | **String** |  | 
**expMuLigandRef** | **String** |  | 
**optType** | **String** |  | 
**optWt** | **Float** |  | 
**optMu** | **Float** |  | 
**optSign** | **String** |  | 
**optPercentage** | **Float** |  | 
**optQual** | **String** |  | 
**optAgonist** | **String** |  | 




# ProteinSimilarity

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**similarity** | **Integer** |  | 
**identity** | **Integer** |  | 
**AA** | **String** |  | 




# ResidueExtendedSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequenceNumber** | **Integer** |  | 
**aminoAcid** | **String** |  | 
**proteinSegment** | **String** |  | 
**displayGenericNumber** | **String** |  | 
**alternativeGenericNumbers** | [**List&lt;ResidueGenericNumberSerializer&gt;**](ResidueGenericNumberSerializer.md) |  | 



# ServicesproteinfamilyApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**proteinFamilyChildrenListGET**](ServicesproteinfamilyApi.md#proteinFamilyChildrenListGET) | **GET** /services/proteinfamily/children/{slug}/ | Get a list of child families of a protein family
[**proteinFamilyDescendantListGET**](ServicesproteinfamilyApi.md#proteinFamilyDescendantListGET) | **GET** /services/proteinfamily/descendants/{slug}/ | Get a list of descendant families of a protein family
[**proteinFamilyDetailGET**](ServicesproteinfamilyApi.md#proteinFamilyDetailGET) | **GET** /services/proteinfamily/{slug}/ | Get a single protein family instance
[**proteinFamilyListGET**](ServicesproteinfamilyApi.md#proteinFamilyListGET) | **GET** /services/proteinfamily/ | Get a list of protein families
[**proteinsInFamilyListGET**](ServicesproteinfamilyApi.md#proteinsInFamilyListGET) | **GET** /services/proteinfamily/proteins/{slug}/ | Get a list of proteins in a protein family
[**proteinsInFamilySpeciesListGET**](ServicesproteinfamilyApi.md#proteinsInFamilySpeciesListGET) | **GET** /services/proteinfamily/proteins/{slug}/{latin_name}/ | Get a list of proteins in a protein family


<a name="proteinFamilyChildrenListGET"></a>
# **proteinFamilyChildrenListGET**
> ProteinFamilySerializer proteinFamilyChildrenListGET(slug)

Get a list of child families of a protein family

Get a list of child families of a protein family&lt;br/&gt;/proteinfamily/children/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
String slug = "slug_example"; // String | 
try {
    ProteinFamilySerializer result = apiInstance.proteinFamilyChildrenListGET(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinFamilyChildrenListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

[**ProteinFamilySerializer**](ProteinFamilySerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinFamilyDescendantListGET"></a>
# **proteinFamilyDescendantListGET**
> ProteinFamilySerializer proteinFamilyDescendantListGET(slug)

Get a list of descendant families of a protein family

Get a list of descendant families of a protein family&lt;br/&gt;/proteinfamily/descendants/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
String slug = "slug_example"; // String | 
try {
    ProteinFamilySerializer result = apiInstance.proteinFamilyDescendantListGET(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinFamilyDescendantListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

[**ProteinFamilySerializer**](ProteinFamilySerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinFamilyDetailGET"></a>
# **proteinFamilyDetailGET**
> ProteinFamilySerializer proteinFamilyDetailGET(slug)

Get a single protein family instance

Get a single protein family instance&lt;br/&gt;/proteinfamily/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
String slug = "slug_example"; // String | 
try {
    ProteinFamilySerializer result = apiInstance.proteinFamilyDetailGET(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinFamilyDetailGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

[**ProteinFamilySerializer**](ProteinFamilySerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinFamilyListGET"></a>
# **proteinFamilyListGET**
> List&lt;ProteinFamilySerializer&gt; proteinFamilyListGET()

Get a list of protein families

Get a list of protein families&lt;br/&gt;/proteinfamily/

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
try {
    List<ProteinFamilySerializer> result = apiInstance.proteinFamilyListGET();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinFamilyListGET");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

[**List&lt;ProteinFamilySerializer&gt;**](ProteinFamilySerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinsInFamilyListGET"></a>
# **proteinsInFamilyListGET**
> List&lt;ProteinSerializer&gt; proteinsInFamilyListGET(slug)

Get a list of proteins in a protein family

Get a list of proteins in a protein family&lt;br/&gt;/proteinfamily/proteins/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
String slug = "slug_example"; // String | 
try {
    List<ProteinSerializer> result = apiInstance.proteinsInFamilyListGET(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinsInFamilyListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

[**List&lt;ProteinSerializer&gt;**](ProteinSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinsInFamilySpeciesListGET"></a>
# **proteinsInFamilySpeciesListGET**
> List&lt;ProteinSerializer&gt; proteinsInFamilySpeciesListGET(slug, latinName)

Get a list of proteins in a protein family

Get a list of proteins in a protein family&lt;br/&gt;/proteinfamily/proteins/{slug}/{species}&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001&lt;br/&gt;{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinfamilyApi;


ServicesproteinfamilyApi apiInstance = new ServicesproteinfamilyApi();
String slug = "slug_example"; // String | 
String latinName = "latinName_example"; // String | 
try {
    List<ProteinSerializer> result = apiInstance.proteinsInFamilySpeciesListGET(slug, latinName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinfamilyApi#proteinsInFamilySpeciesListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |
 **latinName** | **String**|  |

### Return type

[**List&lt;ProteinSerializer&gt;**](ProteinSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined


# WriteResidueGenericNumberSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**label** | **String** |  | 




# StructureLigandInteractionSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdbCode** | **String** |  | 
**ligandName** | **String** |  | 
**aminoAcid** | **String** |  | 
**sequenceNumber** | **Integer** |  | 
**displayGenericNumber** | **String** |  | 
**interactionType** | **String** |  | 




# WriteResidueSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequenceNumber** | **Integer** |  | 
**aminoAcid** | **String** |  | 



# ServicesproteinApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**proteinByAccessionDetailGET**](ServicesproteinApi.md#proteinByAccessionDetailGET) | **GET** /services/protein/accession/{accession}/ | Get a single protein instance by accession
[**proteinDetailGET**](ServicesproteinApi.md#proteinDetailGET) | **GET** /services/protein/{entry_name}/ | Get a single protein instance by entry name


<a name="proteinByAccessionDetailGET"></a>
# **proteinByAccessionDetailGET**
> ProteinSerializer proteinByAccessionDetailGET(accession)

Get a single protein instance by accession

Get a single protein instance by accession&lt;br/&gt;/protein/accession/{accession}/&lt;br/&gt;{accession} is a protein identifier from Uniprot, e.g. P07550

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinApi;


ServicesproteinApi apiInstance = new ServicesproteinApi();
String accession = "accession_example"; // String | 
try {
    ProteinSerializer result = apiInstance.proteinByAccessionDetailGET(accession);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinApi#proteinByAccessionDetailGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **accession** | **String**|  |

### Return type

[**ProteinSerializer**](ProteinSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinDetailGET"></a>
# **proteinDetailGET**
> ProteinSerializer proteinDetailGET(entryName)

Get a single protein instance by entry name

Get a single protein instance by entry name&lt;br/&gt;/protein/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesproteinApi;


ServicesproteinApi apiInstance = new ServicesproteinApi();
String entryName = "entryName_example"; // String | 
try {
    ProteinSerializer result = apiInstance.proteinDetailGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesproteinApi#proteinDetailGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**ProteinSerializer**](ProteinSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

# ServicesplotApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**helixBoxGET**](ServicesplotApi.md#helixBoxGET) | **GET** /services/plot/helixbox/{entry_name}/ | Get SVG source code for a protein&#39;s helix box plot
[**snakePlotGET**](ServicesplotApi.md#snakePlotGET) | **GET** /services/plot/snake/{entry_name}/ | Get SVG source code for a protein&#39;s snake plot


<a name="helixBoxGET"></a>
# **helixBoxGET**
> Object helixBoxGET(entryName)

Get SVG source code for a protein&#39;s helix box plot

Get SVG source code for a protein&#39;s helix box plot&lt;br/&gt;/plot/helixbox/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesplotApi;


ServicesplotApi apiInstance = new ServicesplotApi();
String entryName = "entryName_example"; // String | 
try {
    Object result = apiInstance.helixBoxGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesplotApi#helixBoxGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="snakePlotGET"></a>
# **snakePlotGET**
> Object snakePlotGET(entryName)

Get SVG source code for a protein&#39;s snake plot

Get SVG source code for a protein&#39;s snake plot&lt;br/&gt;/plot/snake/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesplotApi;


ServicesplotApi apiInstance = new ServicesplotApi();
String entryName = "entryName_example"; // String | 
try {
    Object result = apiInstance.snakePlotGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesplotApi#snakePlotGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined


# ProteinSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**entryName** | **String** |  | 
**name** | **String** |  | 
**accession** | **String** |  | 
**family** | **String** |  | 
**species** | **String** |  | 
**source** | **String** |  | 
**residueNumberingScheme** | **String** |  | 
**sequence** | **String** |  | 




# Structure

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**publicationDate** | **String** |  | 
**preferredChain** | **String** |  | 
**type** | **String** |  | 
**species** | **String** |  | 
**protein** | **String** |  | 
**resolution** | **Float** |  | 
**ligands** | [**List&lt;Ligand&gt;**](Ligand.md) |  | 
**publication** | **String** |  |  [optional]
**pdbCode** | **String** |  | 
**family** | **String** |  | 




# WriteStructureLigandInteractionSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------




# SpeciesSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**latinName** | **String** |  | 
**commonName** | **String** |  | 




# ProteinSimilarities

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------




# MutationSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**reference** | **String** |  | 
**protein** | **String** |  | 
**mutationPos** | **Integer** |  | 
**mutationFrom** | **String** |  | 
**mutationTo** | **String** |  | 
**ligandName** | **String** |  | 
**ligandIdtype** | **String** |  | 
**ligandId** | **String** |  | 
**ligandClass** | **String** |  | 
**expType** | **String** |  | 
**expFunc** | **String** |  | 
**expWtValue** | **Float** |  | 
**expWtUnit** | **String** |  | 
**expMuEffectSign** | **String** |  | 
**expMuEffectType** | **String** |  | 
**expMuEffectValue** | **Float** |  | 
**expFoldChange** | **Float** |  | 
**expMuEffectQual** | **String** |  | 
**expMuEffectLigandProp** | **String** |  | 
**expMuLigandRef** | **String** |  | 
**optReceptorExpression** | **String** |  | 
**optBasalActivity** | **String** |  | 
**optGainOfActivity** | **String** |  | 
**optLigandEmax** | **String** |  | 
**optAgonist** | **String** |  | 



# ServicesalignmentApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**familyAlignmentGET**](ServicesalignmentApi.md#familyAlignmentGET) | **GET** /services/alignment/family/{slug}/ | Get a full sequence alignment of a protein family including a consensus sequence
[**familyAlignmentGET_0**](ServicesalignmentApi.md#familyAlignmentGET_0) | **GET** /services/alignment/family/{slug}/statistics/ | Get a full sequence alignment of a protein family including a consensus sequence
[**familyAlignmentPartialGET**](ServicesalignmentApi.md#familyAlignmentPartialGET) | **GET** /services/alignment/family/{slug}/{segments}/ | Get a partial sequence alignment of a protein family
[**familyAlignmentPartialGET_0**](ServicesalignmentApi.md#familyAlignmentPartialGET_0) | **GET** /services/alignment/family/{slug}/{segments}/statistics/ | Get a partial sequence alignment of a protein family
[**familyAlignmentPartialSpeciesGET**](ServicesalignmentApi.md#familyAlignmentPartialSpeciesGET) | **GET** /services/alignment/family/{slug}/{segments}/{latin_name}/ | Get a partial sequence alignment of a protein family
[**familyAlignmentPartialSpeciesGET_0**](ServicesalignmentApi.md#familyAlignmentPartialSpeciesGET_0) | **GET** /services/alignment/family/{slug}/{segments}/{latin_name}/statistics/ | Get a partial sequence alignment of a protein family
[**proteinAlignmentGET**](ServicesalignmentApi.md#proteinAlignmentGET) | **GET** /services/alignment/protein/{proteins}/ | Get a full sequence alignment of two or more proteins
[**proteinAlignmentPartialGET**](ServicesalignmentApi.md#proteinAlignmentPartialGET) | **GET** /services/alignment/protein/{proteins}/{segments}/ | Get a partial sequence alignment of two or more proteins
[**proteinAlignmentStatisticsGET**](ServicesalignmentApi.md#proteinAlignmentStatisticsGET) | **GET** /services/alignment/protein/{proteins}/statistics/ | Add a /statics at the end of an alignment in order to
[**proteinAlignmentStatisticsGET_0**](ServicesalignmentApi.md#proteinAlignmentStatisticsGET_0) | **GET** /services/alignment/protein/{proteins}/{segments}/statistics/ | Add a /statics at the end of an alignment in order to
[**proteinSimilaritySearchAlignmentGET**](ServicesalignmentApi.md#proteinSimilaritySearchAlignmentGET) | **GET** /services/alignment/similarity/{proteins}/{segments}/ | Get a segment sequence alignment of two or more proteins ranked by similarity


<a name="familyAlignmentGET"></a>
# **familyAlignmentGET**
> Object familyAlignmentGET(slug)

Get a full sequence alignment of a protein family including a consensus sequence

Get a full sequence alignment of a protein family including a consensus sequence&lt;br/&gt;/alignment/family/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentGET(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="familyAlignmentGET_0"></a>
# **familyAlignmentGET_0**
> Object familyAlignmentGET_0(slug)

Get a full sequence alignment of a protein family including a consensus sequence

Get a full sequence alignment of a protein family including a consensus sequence&lt;br/&gt;/alignment/family/{slug}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentGET_0(slug);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentGET_0");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="familyAlignmentPartialGET"></a>
# **familyAlignmentPartialGET**
> Object familyAlignmentPartialGET(slug, segments)

Get a partial sequence alignment of a protein family

Get a partial sequence alignment of a protein family&lt;br/&gt;/alignment/family/{slug}/{segments}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
String segments = "segments_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentPartialGET(slug, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentPartialGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |
 **segments** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="familyAlignmentPartialGET_0"></a>
# **familyAlignmentPartialGET_0**
> Object familyAlignmentPartialGET_0(slug, segments)

Get a partial sequence alignment of a protein family

Get a partial sequence alignment of a protein family&lt;br/&gt;/alignment/family/{slug}/{segments}/&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
String segments = "segments_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentPartialGET_0(slug, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentPartialGET_0");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |
 **segments** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="familyAlignmentPartialSpeciesGET"></a>
# **familyAlignmentPartialSpeciesGET**
> Object familyAlignmentPartialSpeciesGET(slug, segments, latinName)

Get a partial sequence alignment of a protein family

Get a partial sequence alignment of a protein family&lt;br/&gt;/alignment/family/{slug}/{segments}/{species}&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50&lt;br/&gt;{species} is a species identifier from Uniprot, e.g. Homo sapiens

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
String segments = "segments_example"; // String | 
String latinName = "latinName_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentPartialSpeciesGET(slug, segments, latinName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentPartialSpeciesGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |
 **segments** | **String**|  |
 **latinName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="familyAlignmentPartialSpeciesGET_0"></a>
# **familyAlignmentPartialSpeciesGET_0**
> Object familyAlignmentPartialSpeciesGET_0(slug, segments, latinName)

Get a partial sequence alignment of a protein family

Get a partial sequence alignment of a protein family&lt;br/&gt;/alignment/family/{slug}/{segments}/{species}&lt;br/&gt;{slug} is a protein family identifier, e.g. 001_001_001&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50&lt;br/&gt;{species} is a species identifier from Uniprot, e.g. Homo sapiens

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String slug = "slug_example"; // String | 
String segments = "segments_example"; // String | 
String latinName = "latinName_example"; // String | 
try {
    Object result = apiInstance.familyAlignmentPartialSpeciesGET_0(slug, segments, latinName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#familyAlignmentPartialSpeciesGET_0");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **slug** | **String**|  |
 **segments** | **String**|  |
 **latinName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinAlignmentGET"></a>
# **proteinAlignmentGET**
> Object proteinAlignmentGET(proteins)

Get a full sequence alignment of two or more proteins

Get a full sequence alignment of two or more proteins&lt;br/&gt;/alignment/protein/{proteins}/&lt;br/&gt;{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String proteins = "proteins_example"; // String | 
try {
    Object result = apiInstance.proteinAlignmentGET(proteins);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#proteinAlignmentGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **proteins** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinAlignmentPartialGET"></a>
# **proteinAlignmentPartialGET**
> Object proteinAlignmentPartialGET(proteins, segments)

Get a partial sequence alignment of two or more proteins

Get a partial sequence alignment of two or more proteins&lt;br/&gt;/alignment/protein/{proteins}/{segments}/&lt;br/&gt;{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String proteins = "proteins_example"; // String | 
String segments = "segments_example"; // String | 
try {
    Object result = apiInstance.proteinAlignmentPartialGET(proteins, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#proteinAlignmentPartialGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **proteins** | **String**|  |
 **segments** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinAlignmentStatisticsGET"></a>
# **proteinAlignmentStatisticsGET**
> Object proteinAlignmentStatisticsGET(proteins)

Add a /statics at the end of an alignment in order to

Add a /statics at the end of an alignment in order to      receive an additional residue property statistics output e.g.:&lt;br/&gt;/alignment/protein/{proteins}/{segments}/statistics&lt;br/&gt;{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String proteins = "proteins_example"; // String | 
try {
    Object result = apiInstance.proteinAlignmentStatisticsGET(proteins);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#proteinAlignmentStatisticsGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **proteins** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinAlignmentStatisticsGET_0"></a>
# **proteinAlignmentStatisticsGET_0**
> Object proteinAlignmentStatisticsGET_0(proteins, segments)

Add a /statics at the end of an alignment in order to

Add a /statics at the end of an alignment in order to      receive an additional residue property statistics output e.g.:&lt;br/&gt;/alignment/protein/{proteins}/{segments}/statistics&lt;br/&gt;{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String proteins = "proteins_example"; // String | 
String segments = "segments_example"; // String | 
try {
    Object result = apiInstance.proteinAlignmentStatisticsGET_0(proteins, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#proteinAlignmentStatisticsGET_0");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **proteins** | **String**|  |
 **segments** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="proteinSimilaritySearchAlignmentGET"></a>
# **proteinSimilaritySearchAlignmentGET**
> ProteinSimilarities proteinSimilaritySearchAlignmentGET(proteins, segments)

Get a segment sequence alignment of two or more proteins ranked by similarity

Get a segment sequence alignment of two or more proteins ranked by similarity&lt;br/&gt;/alignment/similarity/{proteins}/&lt;br/&gt;{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human,cxcr4_human,     where the first protein is the query protein and the following the proteins to compare it to&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers and/ or     generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesalignmentApi;


ServicesalignmentApi apiInstance = new ServicesalignmentApi();
String proteins = "proteins_example"; // String | 
String segments = "segments_example"; // String | 
try {
    ProteinSimilarities result = apiInstance.proteinSimilaritySearchAlignmentGET(proteins, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesalignmentApi#proteinSimilaritySearchAlignmentGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **proteins** | **String**|  |
 **segments** | **String**|  |

### Return type

[**ProteinSimilarities**](ProteinSimilarities.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined


# WriteProteinSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**entryName** | **String** |  | 
**name** | **String** |  | 
**accession** | **String** |  |  [optional]
**sequence** | **String** |  | 




# WriteParentProteinFamilySerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**slug** | **String** |  | 
**name** | **String** |  | 




# Ligand

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**name** | **String** |  | 
**type** | **String** |  | 
**function** | **String** |  | 




# ParentProteinFamilySerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**slug** | **String** |  | 
**name** | **String** |  | 



# ServicesspeciesApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**speciesDetailGET**](ServicesspeciesApi.md#speciesDetailGET) | **GET** /services/species/{latin_name}/ | Get a single species instance
[**speciesListGET**](ServicesspeciesApi.md#speciesListGET) | **GET** /services/species/ | Get a list of species


<a name="speciesDetailGET"></a>
# **speciesDetailGET**
> SpeciesSerializer speciesDetailGET(latinName)

Get a single species instance

Get a single species instance&lt;br/&gt;/species/{latin_name}/&lt;br/&gt;{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesspeciesApi;


ServicesspeciesApi apiInstance = new ServicesspeciesApi();
String latinName = "latinName_example"; // String | 
try {
    SpeciesSerializer result = apiInstance.speciesDetailGET(latinName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesspeciesApi#speciesDetailGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **latinName** | **String**|  |

### Return type

[**SpeciesSerializer**](SpeciesSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="speciesListGET"></a>
# **speciesListGET**
> SpeciesSerializer speciesListGET()

Get a list of species

Get a list of species&lt;br/&gt;/species/

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesspeciesApi;


ServicesspeciesApi apiInstance = new ServicesspeciesApi();
try {
    SpeciesSerializer result = apiInstance.speciesListGET();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesspeciesApi#speciesListGET");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

[**SpeciesSerializer**](SpeciesSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined


# WriteProteinFamilySerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**slug** | **String** |  | 
**name** | **String** |  | 




# WriteSpeciesSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**latinName** | **String** |  | 
**commonName** | **String** |  |  [optional]




# WriteResidueExtendedSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequenceNumber** | **Integer** |  | 
**aminoAcid** | **String** |  | 




# ProteinFamilySerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**slug** | **String** |  | 
**name** | **String** |  | 
**parent** | [**ParentProteinFamilySerializer**](ParentProteinFamilySerializer.md) |  | 



# ServicesresiduesApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**residuesExtendedListGET**](ServicesresiduesApi.md#residuesExtendedListGET) | **GET** /services/residues/extended/{entry_name}/ | Get a list of residues of a protein, including alternative generic numbers
[**residuesListGET**](ServicesresiduesApi.md#residuesListGET) | **GET** /services/residues/{entry_name}/ | Get a list of residues of a protein


<a name="residuesExtendedListGET"></a>
# **residuesExtendedListGET**
> List&lt;ResidueExtendedSerializer&gt; residuesExtendedListGET(entryName)

Get a list of residues of a protein, including alternative generic numbers

Get a list of residues of a protein, including alternative generic numbers&lt;br/&gt;/residues/extended/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesresiduesApi;


ServicesresiduesApi apiInstance = new ServicesresiduesApi();
String entryName = "entryName_example"; // String | 
try {
    List<ResidueExtendedSerializer> result = apiInstance.residuesExtendedListGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesresiduesApi#residuesExtendedListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**List&lt;ResidueExtendedSerializer&gt;**](ResidueExtendedSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="residuesListGET"></a>
# **residuesListGET**
> List&lt;ResidueSerializer&gt; residuesListGET(entryName)

Get a list of residues of a protein

Get a list of residues of a protein&lt;br/&gt;/residues/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesresiduesApi;


ServicesresiduesApi apiInstance = new ServicesresiduesApi();
String entryName = "entryName_example"; // String | 
try {
    List<ResidueSerializer> result = apiInstance.residuesListGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesresiduesApi#residuesListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**List&lt;ResidueSerializer&gt;**](ResidueSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined


# ResidueGenericNumberSerializer

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**scheme** | **String** |  | 
**label** | **String** |  | 



# ServicesstructureApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**representativeStructureListGET**](ServicesstructureApi.md#representativeStructureListGET) | **GET** /services/structure/representative/ | Get a list of representative structures (one for each protein and activation state)
[**representativeStructureListProteinGET**](ServicesstructureApi.md#representativeStructureListProteinGET) | **GET** /services/structure/protein/{entry_name}/representative/ | Get a list of representative structures of a protein (one for each activation state)
[**structureAssignGenericNumbersPOST**](ServicesstructureApi.md#structureAssignGenericNumbersPOST) | **POST** /services/structure/assign_generic_numbers | Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file
[**structureDetailGET**](ServicesstructureApi.md#structureDetailGET) | **GET** /services/structure/{pdb_code}/ | Get a single structure instance
[**structureLigandInteractionsGET**](ServicesstructureApi.md#structureLigandInteractionsGET) | **GET** /services/structure/{pdb_code}/interaction/ | Get a list of interactions between structure and ligand
[**structureListGET**](ServicesstructureApi.md#structureListGET) | **GET** /services/structure/ | Get a list of structures
[**structureListProteinGET**](ServicesstructureApi.md#structureListProteinGET) | **GET** /services/structure/protein/{entry_name}/ | Get a list of structures of a protein
[**structureSingleProteinGET**](ServicesstructureApi.md#structureSingleProteinGET) | **GET** /services/structure/protein/{entry_name}/?single | Get the structure of a protein
[**structureTemplateGET**](ServicesstructureApi.md#structureTemplateGET) | **GET** /services/structure/template/{entry_name}/ | Get the most similar structure template for a protein using a 7TM alignment
[**structureTemplatePartialGET**](ServicesstructureApi.md#structureTemplatePartialGET) | **GET** /services/structure/template/{entry_name}/{segments}/ | Get the most similar structure template for a protein using a partial alignment


<a name="representativeStructureListGET"></a>
# **representativeStructureListGET**
> Object representativeStructureListGET()

Get a list of representative structures (one for each protein and activation state)

Get a list of representative structures (one for each protein and activation state)&lt;br/&gt;/structure/representative/

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
try {
    Object result = apiInstance.representativeStructureListGET();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#representativeStructureListGET");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="representativeStructureListProteinGET"></a>
# **representativeStructureListProteinGET**
> Object representativeStructureListProteinGET(entryName)

Get a list of representative structures of a protein (one for each activation state)

Get a list of representative structures of a protein (one for each activation state)&lt;br/&gt;/structure/protein/{entry_name}/representative/

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String entryName = "entryName_example"; // String | 
try {
    Object result = apiInstance.representativeStructureListProteinGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#representativeStructureListProteinGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureAssignGenericNumbersPOST"></a>
# **structureAssignGenericNumbersPOST**
> Object structureAssignGenericNumbersPOST()

Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file

Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file.&lt;br/&gt;/structure/assign_generic_numbers&lt;br/&gt;    e.g.     curl -X POST -F \&quot;pdb_file&#x3D;@myfile.pdb\&quot; https://gpcrdb.org/services/structure/assign_generic_numbers

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
try {
    Object result = apiInstance.structureAssignGenericNumbersPOST();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureAssignGenericNumbersPOST");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureDetailGET"></a>
# **structureDetailGET**
> Object structureDetailGET(pdbCode)

Get a single structure instance

Get a single structure instance&lt;br/&gt;/structure/{pdb_code}/&lt;br/&gt;{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String pdbCode = "pdbCode_example"; // String | 
try {
    Object result = apiInstance.structureDetailGET(pdbCode);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureDetailGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **pdbCode** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureLigandInteractionsGET"></a>
# **structureLigandInteractionsGET**
> List&lt;StructureLigandInteractionSerializer&gt; structureLigandInteractionsGET(pdbCode)

Get a list of interactions between structure and ligand

Get a list of interactions between structure and ligand&lt;br/&gt;/structure/{pdb_code}/interaction/&lt;br/&gt;{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String pdbCode = "pdbCode_example"; // String | 
try {
    List<StructureLigandInteractionSerializer> result = apiInstance.structureLigandInteractionsGET(pdbCode);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureLigandInteractionsGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **pdbCode** | **String**|  |

### Return type

[**List&lt;StructureLigandInteractionSerializer&gt;**](StructureLigandInteractionSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureListGET"></a>
# **structureListGET**
> List&lt;Structure&gt; structureListGET()

Get a list of structures

Get a list of structures&lt;br/&gt;/structure/

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
try {
    List<Structure> result = apiInstance.structureListGET();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureListGET");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

[**List&lt;Structure&gt;**](Structure.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureListProteinGET"></a>
# **structureListProteinGET**
> List&lt;Structure&gt; structureListProteinGET(entryName)

Get a list of structures of a protein

Get a list of structures of a protein&lt;br/&gt;/structure/protein/{entry_name}

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String entryName = "entryName_example"; // String | 
try {
    List<Structure> result = apiInstance.structureListProteinGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureListProteinGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**List&lt;Structure&gt;**](Structure.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureSingleProteinGET"></a>
# **structureSingleProteinGET**
> Structure structureSingleProteinGET(entryName)

Get the structure of a protein

Get the structure of a protein&lt;br/&gt;/structure/protein/{entry_name}

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String entryName = "entryName_example"; // String | 
try {
    Structure result = apiInstance.structureSingleProteinGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureSingleProteinGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**Structure**](Structure.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureTemplateGET"></a>
# **structureTemplateGET**
> Object structureTemplateGET(entryName)

Get the most similar structure template for a protein using a 7TM alignment

Get the most similar structure template for a protein using a 7TM alignment&lt;br/&gt;/structure/template/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String entryName = "entryName_example"; // String | 
try {
    Object result = apiInstance.structureTemplateGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureTemplateGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

<a name="structureTemplatePartialGET"></a>
# **structureTemplatePartialGET**
> Object structureTemplatePartialGET(entryName, segments)

Get the most similar structure template for a protein using a partial alignment

Get the most similar structure template for a protein using a partial alignment&lt;br/&gt;/structure/template/{entry_name}/{segments}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human&lt;br/&gt;{segments} is a comma separated list of protein segment identifiers, e.g. TM3,TM5,TM6

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesstructureApi;


ServicesstructureApi apiInstance = new ServicesstructureApi();
String entryName = "entryName_example"; // String | 
String segments = "segments_example"; // String | 
try {
    Object result = apiInstance.structureTemplatePartialGET(entryName, segments);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesstructureApi#structureTemplatePartialGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |
 **segments** | **String**|  |

### Return type

**Object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

# ServicesmutantsApi

All URIs are relative to *https://gpcrdb.org*

Method | HTTP request | Description
------------- | ------------- | -------------
[**mutantListGET**](ServicesmutantsApi.md#mutantListGET) | **GET** /services/mutants/{entry_name}/ | Get a list of mutants of single protein instance by entry name


<a name="mutantListGET"></a>
# **mutantListGET**
> List&lt;MutationSerializer&gt; mutantListGET(entryName)

Get a list of mutants of single protein instance by entry name

Get a list of mutants of single protein instance by entry name&lt;br/&gt;/mutant/{entry_name}/&lt;br/&gt;{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human

### Example
```java
// Import classes:
//import nl.esciencecenter.e3dchem.gpcrdb.client.ApiException;
//import nl.esciencecenter.e3dchem.gpcrdb.client.ServicesmutantsApi;


ServicesmutantsApi apiInstance = new ServicesmutantsApi();
String entryName = "entryName_example"; // String | 
try {
    List<MutationSerializer> result = apiInstance.mutantListGET(entryName);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling ServicesmutantsApi#mutantListGET");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **entryName** | **String**|  |

### Return type

[**List&lt;MutationSerializer&gt;**](MutationSerializer.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

