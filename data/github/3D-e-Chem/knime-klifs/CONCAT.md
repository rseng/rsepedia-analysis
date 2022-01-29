# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/.

## [Unreleased]

## [1.1.2] - 2019-07-02

### Changes

- Requires KNIME 4.0 [#6](https://github.com/3D-e-Chem/knime-klifs/issues/6)
- Updated testflow

## [1.1.1] - 2017-03-05

### Changed
- Added PDB download option (test)
- Updated test settings for Mac OS X
- Updated the repository/update site
- Handling structures without IFPs by the Interaction Fingerprint Retriever node

## [1.1.0] - 2017-01-20
- Time-out option (#1)
- Usual suspects for settings configuration (#2)

## [1.0.16] - 2016-12-21

### Changed
- Error messages handling KLIFS server
- Cleaned code
- Updated testflow

## [1.0.15] - 2016-12-20

### Added
- Base URL option to specify the location the KLIFS web services (e.g. for proxy/forwarding settings)

### Changed
- Updated testflows
- Type of SMILES columns from the "Ligands Overview Retriever" from StringCell to SmilesCell

## [1.0.14] - 2016-10-21

### Changed

- Added KNIME testflow

### Fixed

- Fixed handling of empty tables by the Interaction Fingerprint Retriever node
- Cleaned code based on Codacy flags

## [1.0.13] - 2016-10-05

### Changed

- Client is now based on okhttp/gson libraries to prevent conflicts within KNIME
- Project is linked to Codacy and Travis-CI

### Fixed

- When and empty input table is provided the nodes will now return an emtpy list instead of polling the KLIFS webservice.
# KLIFS KNIME nodes

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-klifs.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-klifs)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e48d6b99a76d479ea37f13f6f3189b5b)](https://www.codacy.com/app/AJK-dev/knime-klifs?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/knime-klifs&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/20180/3D-e-Chem/knime-klifs.svg)](https://zenodo.org/badge/latestdoi/20180/3D-e-Chem/knime-klifs)

KNIME nodes for retrieving data from KLIFS (http://klifs.vu-compmedchem.nl). KLIFS is a structural kinase-ligand interaction database. For more information regarding KLIFS see [the website](http://klifs.vu-compmedchem.nl) and the references at the bottom.

# Installation

Requirements:

* KNIME, https://www.knime.org

Steps to get KLIFS nodes inside KNIME:

1. Goto Help > Install new software ... menu
2. Press add button
3. Fill text fields with `https://3d-e-chem.github.io/knime-node-collection`
4. Select --all sites-- in work with pulldown
6. Select KLIFS knime nodes
7. Install software & restart (for now an "Unsigned Content" warning can popup during the installation, you can safely ignore this)

# Usage

See the example workflow in `examples` folder.

# Build

```
mvn verify
```

Jar has been made in `plugin/target` folder.
An Eclipse update site will be made in `p2/target/repository` repository.

# Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.sstea.targetplatform/KNIME-AP-4.0.target` target definition.

During import the Tycho Eclipse providers must be installed.

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Append new release to 3D-e-Chem update site
  1. Make clone of https://github.com/3D-e-Chem/3D-e-Chem.github.io repo
  2. Append release to 3D-e-Chem update site with `mvn install -Dtarget.update.site=<3D-e-Chem repo/updates>`
5. Commit and push changes in this repo and 3D-e-Chem.github.io repo
6. Make nodes available to 3D-e-Chem KNIME feature by following steps at https://github.com/3D-e-Chem/knime-node-collection#new-release
7. Create a GitHub release
8. Update DOI in CITATION.cff

# Create KLIFS client

1. Compile client
```
cd nl.vu_compmedchem.klifs.client
mvn package
```

2. Make client jar and it's dependencies available in plugin
```
cp -r target/lib/* target/*jar ../nl.vu_compmedchem.klifs/lib/
```

3. Update `plugin/META-INF/MANIFEST.MF`, `plugin/build.properties` files to reflect contents of lib/

# References

* Kooistra, A. J.; Kanev, G. K.; van Linden, O. P.; Leurs, R.; de Esch, I. J.; de Graaf, C. Klifs: A Structural Kinase-Ligand Interaction Database. *Nucleic Acids Res.* **2016**, *44*, D365-371. [10.1093/nar/gkv1082](http://dx.doi.org/10.1093/nar/gkv1082)
* van Linden, O. P.; Kooistra, A. J.; Leurs, R.; de Esch, I. J.; de Graaf, C. Klifs: A Knowledge-Based Structural Database to Navigate Kinase–Ligand Interaction Space. *J. Med. Chem.* **2013**, *57*, 249-277. [10.1021/jm400378w](http://dx.doi.org/10.1021/jm400378w)
# swagger-java-client

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
    <groupId>io.swagger</groupId>
    <artifactId>swagger-java-client</artifactId>
    <version>1.0.0</version>
    <scope>compile</scope>
</dependency>
```

### Gradle users

Add this dependency to your project's build file:

```groovy
compile "io.swagger:swagger-java-client:1.0.0"
```

### Others

At first generate the JAR by executing:

    mvn package

Then manually install the following JARs:

* target/swagger-java-client-1.0.0.jar
* target/lib/*.jar

## Getting Started

Please follow the [installation](#installation) instruction and execute the following Java code:

```java

import io.swagger.client.*;
import io.swagger.client.auth.*;
import io.swagger.client.model.*;
import io.swagger.client.api.InformationApi;

import java.io.File;
import java.util.*;

public class InformationApiExample {

    public static void main(String[] args) {
        
        InformationApi apiInstance = new InformationApi();
        String kinaseGroup = "kinaseGroup_example"; // String | Optional: Name (or multiple names separated by a comma) of the kinase group for which the kinase families are requested (e.g. TKL,STE).
        try {
            List<String> result = apiInstance.kinaseFamiliesGet(kinaseGroup);
            System.out.println(result);
        } catch (ApiException e) {
            System.err.println("Exception when calling InformationApi#kinaseFamiliesGet");
            e.printStackTrace();
        }
    }
}

```

## Documentation for API Endpoints

All URIs are relative to *http://klifs.vu-compmedchem.nl/api*

Class | Method | HTTP request | Description
------------ | ------------- | ------------- | -------------
*InformationApi* | [**kinaseFamiliesGet**](docs/InformationApi.md#kinaseFamiliesGet) | **GET** /kinase_families | Kinase families
*InformationApi* | [**kinaseGroupsGet**](docs/InformationApi.md#kinaseGroupsGet) | **GET** /kinase_groups | Kinase groups
*InformationApi* | [**kinaseIDGet**](docs/InformationApi.md#kinaseIDGet) | **GET** /kinase_ID | Kinase ID
*InformationApi* | [**kinaseInformationGet**](docs/InformationApi.md#kinaseInformationGet) | **GET** /kinase_information | Kinase information
*InformationApi* | [**kinaseNamesGet**](docs/InformationApi.md#kinaseNamesGet) | **GET** /kinase_names | Kinase names
*InteractionsApi* | [**interactionsGetIFPGet**](docs/InteractionsApi.md#interactionsGetIFPGet) | **GET** /interactions_get_IFP | Get structure IFP
*InteractionsApi* | [**interactionsGetTypesGet**](docs/InteractionsApi.md#interactionsGetTypesGet) | **GET** /interactions_get_types | Get interaction types
*InteractionsApi* | [**interactionsMatchResiduesGet**](docs/InteractionsApi.md#interactionsMatchResiduesGet) | **GET** /interactions_match_residues | Match IFP residues
*LigandsApi* | [**ligandsListGet**](docs/LigandsApi.md#ligandsListGet) | **GET** /ligands_list | Get all co-crystallized ligands optionally restricted to a set of kinase IDs
*LigandsApi* | [**ligandsListStructuresGet**](docs/LigandsApi.md#ligandsListStructuresGet) | **GET** /ligands_list_structures | Get all structures in complex with one of the provided ligand IDs
*StructuresApi* | [**structureGetComplexGet**](docs/StructuresApi.md#structureGetComplexGet) | **GET** /structure_get_complex | Get full complex
*StructuresApi* | [**structureGetLigandGet**](docs/StructuresApi.md#structureGetLigandGet) | **GET** /structure_get_ligand | Get ligand from structure
*StructuresApi* | [**structureGetPocketGet**](docs/StructuresApi.md#structureGetPocketGet) | **GET** /structure_get_pocket | Get pocket from structure
*StructuresApi* | [**structureGetProteinGet**](docs/StructuresApi.md#structureGetProteinGet) | **GET** /structure_get_protein | Get protein from structure
*StructuresApi* | [**structuresListGet**](docs/StructuresApi.md#structuresListGet) | **GET** /structures_list | Get all structures based on a kinase ID
*StructuresApi* | [**structuresPdbListGet**](docs/StructuresApi.md#structuresPdbListGet) | **GET** /structures_pdb_list | Get all structures based on a set of PDB-codes


## Documentation for Models

 - [Error](docs/Error.md)
 - [IDlist](docs/IDlist.md)
 - [IFPList](docs/IFPList.md)
 - [InteractionList](docs/InteractionList.md)
 - [KinaseInformation](docs/KinaseInformation.md)
 - [LigandDetails](docs/LigandDetails.md)
 - [MatchList](docs/MatchList.md)
 - [StructureDetails](docs/StructureDetails.md)


## Documentation for Authorization

All endpoints do not require authorization.
Authentication schemes defined for the API:

## Recommendation

It's recommended to create an instance of `ApiClient` per thread in a multithreaded environment to avoid any potential issues.

## Author



# LigandsApi

All URIs are relative to *http://klifs.vu-compmedchem.nl/api*

Method | HTTP request | Description
------------- | ------------- | -------------
[**ligandsListGet**](LigandsApi.md#ligandsListGet) | **GET** /ligands_list | Get all co-crystallized ligands optionally restricted to a set of kinase IDs
[**ligandsListStructuresGet**](LigandsApi.md#ligandsListStructuresGet) | **GET** /ligands_list_structures | Get all structures in complex with one of the provided ligand IDs


<a name="ligandsListGet"></a>
# **ligandsListGet**
> List&lt;LigandDetails&gt; ligandsListGet(kinaseID)

Get all co-crystallized ligands optionally restricted to a set of kinase IDs

The Ligands List endpoint returns a list of co-crystallized ligands in KLIFS optiontally restricted to a set of specific kinase IDs. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.LigandsApi;


LigandsApi apiInstance = new LigandsApi();
List<Integer> kinaseID = Arrays.asList(56); // List<Integer> | ID(s) of the kinase(s) for which all structures are requested.
try {
    List<LigandDetails> result = apiInstance.ligandsListGet(kinaseID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling LigandsApi#ligandsListGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseID** | [**List&lt;Integer&gt;**](Integer.md)| ID(s) of the kinase(s) for which all structures are requested. | [optional]

### Return type

[**List&lt;LigandDetails&gt;**](LigandDetails.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="ligandsListStructuresGet"></a>
# **ligandsListStructuresGet**
> List&lt;StructureDetails&gt; ligandsListStructuresGet(ligandID)

Get all structures in complex with one of the provided ligand IDs

The Ligands Get Structures endpoint returns a list of monomers of crystal structures in KLIFS in complex with provided ligand IDs. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.LigandsApi;


LigandsApi apiInstance = new LigandsApi();
List<Integer> ligandID = Arrays.asList(56); // List<Integer> | ID(s) of the ligand(s) for which all structures are requested.
try {
    List<StructureDetails> result = apiInstance.ligandsListStructuresGet(ligandID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling LigandsApi#ligandsListStructuresGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **ligandID** | [**List&lt;Integer&gt;**](Integer.md)| ID(s) of the ligand(s) for which all structures are requested. |

### Return type

[**List&lt;StructureDetails&gt;**](StructureDetails.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json


# MatchList

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**index** | **Integer** |  |  [optional]
**xrayPosition** | **String** |  |  [optional]
**kLIFSPosition** | **String** |  |  [optional]




# LigandDetails

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligandID** | **Integer** | KLIFS ID for a specific kinase structure | 
**pDBCode** | **String** | 3-letter PDB-code for a specific ligand | 
**name** | **String** | Name of the ligand |  [optional]
**SMILES** | **String** | SMILES code for a specific ligand |  [optional]
**inChIKey** | **String** | InChiKey of a specific ligand |  [optional]




# InteractionList

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**position** | **Integer** |  |  [optional]
**name** | **String** |  |  [optional]




# IDlist

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**kinaseID** | **Integer** |  |  [optional]
**name** | **String** |  |  [optional]
**fullName** | **String** |  |  [optional]
**species** | **String** |  |  [optional]




# StructureDetails

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**structureID** | **Integer** | KLIFS ID for a specific kinase structure | 
**kinase** | **String** | HGNC gene symbol for a specific kinase structure |  [optional]
**species** | **String** | Species of a specific kinase structure (e.g. human) | 
**kinaseID** | **Integer** | KLIFS ID for a specific kinase | 
**pdb** | **String** | 4-letter PDB-code for a given structure | 
**alt** | **String** | Indicates the (alternate) model of the given structure |  [optional]
**chain** | **String** | Indicates the chain of the given structure | 
**rmsd1** | **Float** | RMSD compared to the reference template for the superpose pocket | 
**rmsd2** | **Float** | RMSD compared to the reference template for the full pocket | 
**pocket** | **String** | Alignment of the 85 pocket residues | 
**resolution** | **Float** | Resolution of the crystal struture |  [optional]
**qualityScore** | **Float** | The quality score estimates the quality of the structure with respect to the binding pocket as well as the quality of the processing by KLIFS. | 
**missingResidues** | **Integer** | Number of residues missing in the pocket | 
**missingAtoms** | **Integer** | Number of atoms missing from pocket residues (not including the missing residues) | 
**ligand** | **String** | 3-letter PDB-code of the ligand within the main pocket |  [optional]
**allostericLigand** | **String** | 3-letter PDB-code of the ligand outside the main pocket |  [optional]
**DFG** | **String** | Conformation of the DFG motif |  [optional]
**aCHelix** | **String** | Conformation of the alpha-C helix |  [optional]
**grichDistance** | **Float** | Conformation of the G-rich loop - distance |  [optional]
**grichAngle** | **Float** | Conformation of the G-rich loop - angle |  [optional]
**grichRotation** | **Float** | Conformation of the G-rich loop - rotation |  [optional]
**front** | **Boolean** | The ligand binds to the front pocket |  [optional]
**gate** | **Boolean** | The ligand binds to the gate area |  [optional]
**back** | **Boolean** | The ligand binds to the back pocket |  [optional]
**fpI** | **Boolean** | The ligand binds to FP-I |  [optional]
**fpII** | **Boolean** | The ligand binds to FP-II |  [optional]
**bpIA** | **Boolean** | The ligand binds to BP-I-A |  [optional]
**bpIB** | **Boolean** | The ligand binds to BP-I-B |  [optional]
**bpIIIn** | **Boolean** | The ligand binds to BP-II-in |  [optional]
**bpIIAIn** | **Boolean** | The ligand binds to BP-II-A-in |  [optional]
**bpIIBIn** | **Boolean** | The ligand binds to BP-II-B-in |  [optional]
**bpIIOut** | **Boolean** | The ligand binds to BP-II-out |  [optional]
**bpIIB** | **Boolean** | The ligand binds to BP-II-B |  [optional]
**bpIII** | **Boolean** | The ligand binds to BP-III |  [optional]
**bpIV** | **Boolean** | The ligand binds to BP-IV |  [optional]
**bpV** | **Boolean** | The ligand binds to BP-V |  [optional]




# Error

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**code** | **Integer** |  |  [optional]
**message** | **String** |  |  [optional]




# IFPList

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**structureID** | **Integer** |  |  [optional]
**IFP** | **String** |  |  [optional]



# InformationApi

All URIs are relative to *http://klifs.vu-compmedchem.nl/api*

Method | HTTP request | Description
------------- | ------------- | -------------
[**kinaseFamiliesGet**](InformationApi.md#kinaseFamiliesGet) | **GET** /kinase_families | Kinase families
[**kinaseGroupsGet**](InformationApi.md#kinaseGroupsGet) | **GET** /kinase_groups | Kinase groups
[**kinaseIDGet**](InformationApi.md#kinaseIDGet) | **GET** /kinase_ID | Kinase ID
[**kinaseInformationGet**](InformationApi.md#kinaseInformationGet) | **GET** /kinase_information | Kinase information
[**kinaseNamesGet**](InformationApi.md#kinaseNamesGet) | **GET** /kinase_names | Kinase names


<a name="kinaseFamiliesGet"></a>
# **kinaseFamiliesGet**
> List&lt;String&gt; kinaseFamiliesGet(kinaseGroup)

Kinase families

The Kinase families endpoint returns a list of all available kinase families in KLIFS. When a kinase group is specified only those kinase families that are within that kinase group are returned. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InformationApi;


InformationApi apiInstance = new InformationApi();
String kinaseGroup = "kinaseGroup_example"; // String | Optional: Name (or multiple names separated by a comma) of the kinase group for which the kinase families are requested (e.g. TKL,STE).
try {
    List<String> result = apiInstance.kinaseFamiliesGet(kinaseGroup);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InformationApi#kinaseFamiliesGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseGroup** | **String**| Optional: Name (or multiple names separated by a comma) of the kinase group for which the kinase families are requested (e.g. TKL,STE). | [optional]

### Return type

**List&lt;String&gt;**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="kinaseGroupsGet"></a>
# **kinaseGroupsGet**
> List&lt;String&gt; kinaseGroupsGet()

Kinase groups

The Kinase groups endpoint returns a list of all available kinase groups in KLIFS. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InformationApi;


InformationApi apiInstance = new InformationApi();
try {
    List<String> result = apiInstance.kinaseGroupsGet();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InformationApi#kinaseGroupsGet");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

**List&lt;String&gt;**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="kinaseIDGet"></a>
# **kinaseIDGet**
> List&lt;KinaseInformation&gt; kinaseIDGet(kinaseName, species)

Kinase ID

The Kinase ID endpoint returns the KLIFS ID of a specific kinase based on the abbreviated kinase name (e.g. ABL1) or UniProt ID (e.g. P00519). The kinase names from Manning et al., the HUGO Gene Nomenclature Committee, and the MGI are recognized. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InformationApi;


InformationApi apiInstance = new InformationApi();
String kinaseName = "kinaseName_example"; // String | The name (or multiple names separated by a comma) of (a) specfic kinase(s) (e.g. ABL1).
String species = "species_example"; // String | Optional: Species for which the kinase ID are requested (e.g. HUMAN OR MOUSE).
try {
    List<KinaseInformation> result = apiInstance.kinaseIDGet(kinaseName, species);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InformationApi#kinaseIDGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseName** | **String**| The name (or multiple names separated by a comma) of (a) specfic kinase(s) (e.g. ABL1). |
 **species** | **String**| Optional: Species for which the kinase ID are requested (e.g. HUMAN OR MOUSE). | [optional]

### Return type

[**List&lt;KinaseInformation&gt;**](KinaseInformation.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="kinaseInformationGet"></a>
# **kinaseInformationGet**
> List&lt;KinaseInformation&gt; kinaseInformationGet(kinaseID, species)

Kinase information

The Kinase information endpoint returns a list of information related to the requested kinase (Uniprot, pocket sequence, etc.). 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InformationApi;


InformationApi apiInstance = new InformationApi();
List<Integer> kinaseID = Arrays.asList(56); // List<Integer> | ID (or IDs separated by a comma) of the kinase for which the kinase information is requested.
String species = "species_example"; // String | Optional: Species for which the kinase names are requested (e.g. HUMAN OR MOUSE).
try {
    List<KinaseInformation> result = apiInstance.kinaseInformationGet(kinaseID, species);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InformationApi#kinaseInformationGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseID** | [**List&lt;Integer&gt;**](Integer.md)| ID (or IDs separated by a comma) of the kinase for which the kinase information is requested. | [optional]
 **species** | **String**| Optional: Species for which the kinase names are requested (e.g. HUMAN OR MOUSE). | [optional]

### Return type

[**List&lt;KinaseInformation&gt;**](KinaseInformation.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="kinaseNamesGet"></a>
# **kinaseNamesGet**
> List&lt;IDlist&gt; kinaseNamesGet(kinaseGroup, kinaseFamily, species)

Kinase names

The Kinase names endpoint returns a list of all available kinases in KLIFS according using the HGNC gene symbols. When a kinase group or kinase family is specified only those kinase names that are within that kinase group or kinase family are returned. When both a group and a family are specified, only the family is used to process the request. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InformationApi;


InformationApi apiInstance = new InformationApi();
String kinaseGroup = "kinaseGroup_example"; // String | Optional: Name (or multiple names separated by a comma) of the kinase group for which the kinase families are requested (e.g. TKL,STE).
String kinaseFamily = "kinaseFamily_example"; // String | Optional: Name (or multiple names separated by a comma) of the kinase family for which the kinase names are requested (e.g. AUR,WEE).
String species = "species_example"; // String | Optional: Species for which the kinase names are requested (e.g. HUMAN OR MOUSE).
try {
    List<IDlist> result = apiInstance.kinaseNamesGet(kinaseGroup, kinaseFamily, species);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InformationApi#kinaseNamesGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseGroup** | **String**| Optional: Name (or multiple names separated by a comma) of the kinase group for which the kinase families are requested (e.g. TKL,STE). | [optional]
 **kinaseFamily** | **String**| Optional: Name (or multiple names separated by a comma) of the kinase family for which the kinase names are requested (e.g. AUR,WEE). | [optional]
 **species** | **String**| Optional: Species for which the kinase names are requested (e.g. HUMAN OR MOUSE). | [optional]

### Return type

[**List&lt;IDlist&gt;**](IDlist.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

# StructuresApi

All URIs are relative to *http://klifs.vu-compmedchem.nl/api*

Method | HTTP request | Description
------------- | ------------- | -------------
[**structureGetComplexGet**](StructuresApi.md#structureGetComplexGet) | **GET** /structure_get_complex | Get full complex
[**structureGetLigandGet**](StructuresApi.md#structureGetLigandGet) | **GET** /structure_get_ligand | Get ligand from structure
[**structureGetPocketGet**](StructuresApi.md#structureGetPocketGet) | **GET** /structure_get_pocket | Get pocket from structure
[**structureGetProteinGet**](StructuresApi.md#structureGetProteinGet) | **GET** /structure_get_protein | Get protein from structure
[**structuresListGet**](StructuresApi.md#structuresListGet) | **GET** /structures_list | Get all structures based on a kinase ID
[**structuresPdbListGet**](StructuresApi.md#structuresPdbListGet) | **GET** /structures_pdb_list | Get all structures based on a set of PDB-codes


<a name="structureGetComplexGet"></a>
# **structureGetComplexGet**
> File structureGetComplexGet(structureID)

Get full complex

The Get Kinase Complex endpoint returns the full structure (including solvent, cofactors, ligands, etc.) in MOL2 format 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
Integer structureID = 56; // Integer | ID(s) of the structure(s) that is/are requested.
try {
    File result = apiInstance.structureGetComplexGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structureGetComplexGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | **Integer**| ID(s) of the structure(s) that is/are requested. |

### Return type

[**File**](File.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: chemical/x-mol2

<a name="structureGetLigandGet"></a>
# **structureGetLigandGet**
> File structureGetLigandGet(structureID)

Get ligand from structure

The Get kinase ligand endpoint returns the orthosteric ligand of a specific structure in MOL2 format 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
Integer structureID = 56; // Integer | ID(s) of the structure(s) that is/are requested.
try {
    File result = apiInstance.structureGetLigandGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structureGetLigandGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | **Integer**| ID(s) of the structure(s) that is/are requested. |

### Return type

[**File**](File.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: chemical/x-mol2

<a name="structureGetPocketGet"></a>
# **structureGetPocketGet**
> File structureGetPocketGet(structureID)

Get pocket from structure

The Get kinase pocket endpoint returns only the KLIFS pocket of a specific structure (excluding solvent, cofactors, ligands, etc.) in MOL2 format 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
Integer structureID = 56; // Integer | ID(s) of the structure(s) that is/are requested.
try {
    File result = apiInstance.structureGetPocketGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structureGetPocketGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | **Integer**| ID(s) of the structure(s) that is/are requested. |

### Return type

[**File**](File.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: chemical/x-mol2

<a name="structureGetProteinGet"></a>
# **structureGetProteinGet**
> File structureGetProteinGet(structureID)

Get protein from structure

The Get kinase protein endpoint returns the full protein of a specific structure (excluding solvent, cofactors, ligands, etc.) in MOL2 format 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
Integer structureID = 56; // Integer | ID(s) of the structure(s) that is/are requested.
try {
    File result = apiInstance.structureGetProteinGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structureGetProteinGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | **Integer**| ID(s) of the structure(s) that is/are requested. |

### Return type

[**File**](File.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: chemical/x-mol2

<a name="structuresListGet"></a>
# **structuresListGet**
> List&lt;StructureDetails&gt; structuresListGet(kinaseID)

Get all structures based on a kinase ID

The Structures list endpoint returns a list of available kinase structures in KLIFS based on a specific kinase ID (e.g. kinase ID 392, which is the ID for ABL1). 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
List<Integer> kinaseID = Arrays.asList(56); // List<Integer> | ID(s) of the kinase(s) for which all structures are requested.
try {
    List<StructureDetails> result = apiInstance.structuresListGet(kinaseID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structuresListGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **kinaseID** | [**List&lt;Integer&gt;**](Integer.md)| ID(s) of the kinase(s) for which all structures are requested. |

### Return type

[**List&lt;StructureDetails&gt;**](StructureDetails.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="structuresPdbListGet"></a>
# **structuresPdbListGet**
> List&lt;StructureDetails&gt; structuresPdbListGet(pdbCodes)

Get all structures based on a set of PDB-codes

The Structures PDB list endpoint returns a list of available kinase structures in KLIFS based on a set of 4-letter PDB-codes. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.StructuresApi;


StructuresApi apiInstance = new StructuresApi();
List<String> pdbCodes = Arrays.asList("pdbCodes_example"); // List<String> | PDB-codes for which all structures are requested.
try {
    List<StructureDetails> result = apiInstance.structuresPdbListGet(pdbCodes);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling StructuresApi#structuresPdbListGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **pdbCodes** | [**List&lt;String&gt;**](String.md)| PDB-codes for which all structures are requested. |

### Return type

[**List&lt;StructureDetails&gt;**](StructureDetails.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

# InteractionsApi

All URIs are relative to *http://klifs.vu-compmedchem.nl/api*

Method | HTTP request | Description
------------- | ------------- | -------------
[**interactionsGetIFPGet**](InteractionsApi.md#interactionsGetIFPGet) | **GET** /interactions_get_IFP | Get structure IFP
[**interactionsGetTypesGet**](InteractionsApi.md#interactionsGetTypesGet) | **GET** /interactions_get_types | Get interaction types
[**interactionsMatchResiduesGet**](InteractionsApi.md#interactionsMatchResiduesGet) | **GET** /interactions_match_residues | Match IFP residues


<a name="interactionsGetIFPGet"></a>
# **interactionsGetIFPGet**
> List&lt;IFPList&gt; interactionsGetIFPGet(structureID)

Get structure IFP

The Get structure IFP endpoint returns the full IFP a specific kinase structure. This IFP consists of 7 interactions for each of the 85 pocket residues. The presence or absence of an interaction is annotation with a 1 and 0 respectively. This method returns a string of 85 residues x 7 interaction types&#x3D;595 bits. 

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InteractionsApi;


InteractionsApi apiInstance = new InteractionsApi();
List<Integer> structureID = Arrays.asList(56); // List<Integer> | ID(s) of the structure(s) that is/are requested.
try {
    List<IFPList> result = apiInstance.interactionsGetIFPGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InteractionsApi#interactionsGetIFPGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | [**List&lt;Integer&gt;**](Integer.md)| ID(s) of the structure(s) that is/are requested. |

### Return type

[**List&lt;IFPList&gt;**](IFPList.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="interactionsGetTypesGet"></a>
# **interactionsGetTypesGet**
> List&lt;InteractionList&gt; interactionsGetTypesGet()

Get interaction types

The Get interaction types endpoint returns a list of the interactions (1-7) and the type of interaction. IFP encodes 7 types of interactions between the ligand and each pocket residue. E.g. 1000010 represents a hydrophobic interaction (position 1) and an ionic interaction: protein cation - ligand anion (position 6).

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InteractionsApi;


InteractionsApi apiInstance = new InteractionsApi();
try {
    List<InteractionList> result = apiInstance.interactionsGetTypesGet();
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InteractionsApi#interactionsGetTypesGet");
    e.printStackTrace();
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

[**List&lt;InteractionList&gt;**](InteractionList.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

<a name="interactionsMatchResiduesGet"></a>
# **interactionsMatchResiduesGet**
> List&lt;MatchList&gt; interactionsMatchResiduesGet(structureID)

Match IFP residues

The Match IFP residues endpoint returns a list of the pocket definition of IFP including the X-ray numbering for a specific structure. This list can be used to decompose an IFP into specific residue interactions and to match Xray numbering to KLIFS numbering and vice versa.

### Example
```java
// Import classes:
//import io.swagger.client.ApiException;
//import io.swagger.client.api.InteractionsApi;


InteractionsApi apiInstance = new InteractionsApi();
Integer structureID = 56; // Integer | ID of the structure(s) that is requested.
try {
    List<MatchList> result = apiInstance.interactionsMatchResiduesGet(structureID);
    System.out.println(result);
} catch (ApiException e) {
    System.err.println("Exception when calling InteractionsApi#interactionsMatchResiduesGet");
    e.printStackTrace();
}
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **structureID** | **Integer**| ID of the structure(s) that is requested. |

### Return type

[**List&lt;MatchList&gt;**](MatchList.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json


# KinaseInformation

## Properties
Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**kinaseID** | **Integer** |  |  [optional]
**name** | **String** |  |  [optional]
**HGNC** | **String** |  |  [optional]
**family** | **String** |  |  [optional]
**group** | **String** |  |  [optional]
**kinaseClass** | **String** |  |  [optional]
**species** | **String** |  |  [optional]
**fullName** | **String** |  |  [optional]
**uniprot** | **String** |  |  [optional]
**iuphar** | **String** |  |  [optional]
**pocket** | **String** |  |  [optional]



