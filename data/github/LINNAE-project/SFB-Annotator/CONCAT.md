# Changelog

## [0.1.1] - 2021-04-23
- one servlet to handle `GET`, `POST` and `DELETE` requests (#87)
- fixed showing & deleting annotations via GUI (#88, #89, #91)
- consolidate I/O JSON data (#90)

## [0.1.0] - 2021-03-12
- initial release

# Semantic Field Book Annotator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4602263.svg)](https://doi.org/10.5281/zenodo.4602263)
![Build](https://github.com/LINNAE-project/SFB-Annotator/workflows/continuous_integration/badge.svg)
[![Published in J. Web Semant.](https://img.shields.io/badge/published%20in-JWebSemant-blue.svg)](https://doi.org/10.1016/j.websem.2018.06.002)

The Semantic Field Book Annotator is a web application developed for domain experts to harvest structured annotations from field books, drawings and specimen labels of natural history collections. Users can draw bounding boxes over (zoomable) image scans of historical field notes, to which annotations can be attached. All metadata regarding an annotation event, annotation provenance, transcription and semantic interpretation of the text are stored in a knowledge base, which is accessible via a SPARQL endpoint and Web API.

## Prerequisites
- [Docker CE](https://docs.docker.com/install/)
- [Docker Compose](https://docs.docker.com/compose/install/)

## Software used
- [Eclipse RDF4J Server and Workbench](https://github.com/eclipse/rdf4j)
- [grlc](https://github.com/CLARIAH/grlc) Web API
- [Cantaloupe](https://github.com/cantaloupe-project/cantaloupe) IIIF image server
- [Mirador](https://github.com/projectmirador/mirador) IIIF image viewer
- [Annotorius](https://github.com/recogito/annotorious)
- [OpenSeadragon](https://github.com/openseadragon/openseadragon)

## Vocabularies, ontologies & specs used
- [Web Annotation Vocabulary](https://www.w3.org/TR/annotation-vocab/)
- [International Image Interoperability Framework](https://iiif.io/) (IIIF)
- [Darwin Core](https://dwc.tdwg.org/) (DwC)
- [Darwin Semantic Web](http://purl.org/dsw/) (Darwin-SW)
- [Dublin Core](http://purl.org/dc/terms/) (DC)

## Install & deploy

**1. Clone this repository.**

```bash
git clone https://github.com/LINNAE-project/SFB-Annotator.git
```
**2. Install Docker Compose.**

```bash
pip install docker-compose
```

**3. Start Docker service(s).**

| Service | Port | Docker Image | Description |
| ------- | ---- | -------------| ----------- |
| `sea` | `8080` | [`linnae/sfb-annotator`](https://hub.docker.com/r/linnae/sfb-annotator) | Semantic Field Book Annotator |
| `melon` | `8182` | [`linnae/cantaloupe`](https://hub.docker.com/r/linnae/cantaloupe) | Cantaloupe image server |
| `mirador` | `8000` | [`linnae/mirador`](https://hub.docker.com/r/linnae/mirador) | Mirador image viewer |
| `grlc` | `8088` | [`clariah/grlc`](https://hub.docker.com/r/clariah/grlc) | Web API (optional)|

```bash
cd SFB-Annotator
# list available services
docker-compose config --services

# start all services or one-by-one
docker-compose up -d # or append [SERVICE]
```

**4. Configure service(s).**

```bash
# populate RDF store (repository) with example data
docker-compose exec sea ./init.sh

# configure sea to use a remote image archive (optional)
# default: data-local.json
BASE_DIR=/usr/local/tomcat/webapps/semanticAnnotator/data/
IMG_SRC=remote  # default: local
docker exec -it sea bash -c "cp $BASE_DIR/data-$IMG_SRC.json $BASE_DIR/data.json"

# configure grlc to use local path (optional)
git clone https://github.com/LINNAE-project/queries
docker cp ./queries grlc:/home/grlc/
```

**5. Build Docker image and deploy container locally (development)**

```bash
docker build -t linnae/sfb-annotator:local .
docker run --name sea -d -p 8080:8080 linnae/sfb-annotator:local
docker exec sea ./init.sh

# generate RDF triples for example inputs (annotation events)
for json in $(ls data/json/$IMG_SRC/*.json | sort)
do
  prefix="$(basename "$json" .json)"
  suffix=ttl  # or jsonld
  rdf="$prefix.$suffix"
  ./run.sh "$json" "$rdf"
done
```

## Web apps & API endpoints
- http://localhost:8080/semanticAnnotator/
  - requires user/password: `tomcat/tomcat`
  - _Register_->_Save_->_Collections_
- http://localhost:8080/rdf4j-workbench/
  - includes an empty repository: `mem-rdf`
- http://localhost:8080/rdf4j-server/
- http://localhost:8088/ followed by
  - remote path [`/api-git/LINNAE-project/queries/`](http://localhost:8088/api-git/LINNAE-project/queries/) or
    - requires `GRLC_GITHUB_ACCESS_TOKEN` to be set in [`docker-compose.yml`](https://github.com/LINNAE-project/SFB-Annotator/blob/master/docker-compose.yml#L27)
  - local path [`/api-local/`](http://localhost:8088/api-local/)
- http://localhost:8182/iiif/2
  - sample image [`info.json`](http://localhost:8182/iiif/2/900c341c1c10fff7:MMNAT01_PM_NNM001001033_001/info.json)
  - get a [JPG](http://localhost:8182/iiif/2/900c341c1c10fff7:MMNAT01_PM_NNM001001033_001/full/max/0/default.jpg) version of this [TIF](https://trng-repository.surfsara.nl/deposit/900c341c1c10fff7/files/MMNAT01_PM_NNM001001033_001.tif) image in the archive
- http://localhost:8000 includes an example [`manifest.json`](data/manifest.json)
  - _Add Item_->_Sample Field Book_->_Change view type_ to _Book View_
# Jupyter Notebooks

**1. Install dependencies.**

```
conda env create -n nb -f environment.yml
conda activate nb
```

**2. Run examples.**

```
jupyter notebook
```

### **Annotation notes:**
Field notes are collections of observation records ([`dwc:HumanObservation`](http://rs.tdwg.org/dwc/terms/HumanObservation)) that describe the occurrence ([`dwc:Occurrence`](http://rs.tdwg.org/dwc/terms/Occurrence)) of an organism ([`dwc:Organism`](http://rs.tdwg.org/dwc/terms/Organism)) at a certain location ([`dwc:Location`](http://rs.tdwg.org/dwc/terms/Location)) and date ([`nhc:Date`](http://makingsense.liacs.nl/rdf/nhc-content/2018-04-04.html#Date)). An observation record most often includes  an identification ([`dwc:Identification`](https://dwc.tdwg.org/terms/#dwc:Identification)) to a taxon ([`dwc:Taxon`](https://dwc.tdwg.org/terms/#dwc:Taxon)), and accompanying measurements and facts (`dwc:MeasurementOrFact`, `UBERON_0001062`, `ncit:C20189`) that were used for identification. These observation records usually span multiple field notes, but can also be very short utterances such as on page 3: _Sci. Diard Buitenzorg_, most likely referring to the occurrence of another organism with the same identification.

The different variations for annotation are:

|class| variation | example |
|---|---|---
|`dwc:Taxon` | `rdf:type` | [`example 1_1`](#Example-1_1)
|       |`dsw:hasIdentification` | [`example 1_2`](#Example-1_2)
|       |`nhc:additionalIdentification`| [`example 1_3`](#Example-1_3)
|`foaf:Person` | `rdf:type` | [`example 2_1`](#Example-2_1)
| | `nhc:scientificNameAuthorship` |  [`example 2_2`](#Example-2_2)
|       | `dwciri:identifiedBy` | [`example 2_3`](#Example-2_3)
|       | `dwciri:recordedBy`| [`example 2_4`](#Example-2_4)
|`dwc:Location` | `rdf:type` | [`example 3_1`](#Example-3_1)
| | `dsw:locatedAt` | [`example 3_2`](#Example-3_2)
|`dwc:MeasurementOrFact` | `rdf:type` | [`example 4_1`](#Example-4_1)
| | `dsw:derivedFrom` | [`example 4_2`](#Example-4_2) |
|`ncit:C20189` (propertyOrAttribute)| `rdf:type` | [`example 5_1`](#Example-5_1)
| |  `nhc:measuresOrDescribes` | [`example 5_2`](#Example-5_2)|
|`uberon:0001062` (anatomicalEntity)| `rdf:type` | [`example 6_1`](#Example-6_1)
| | `nhc:measuresOrDescribes` | [`example 6_2`](#Example-6_2) |
|`nhc:Date`| `rdf:type` | [`example 7_1`](#Example-7_1)
| | `nhc:verbatimDate` | [`example 7_2`](#Example-7_2)|


Figure 1: Example occurrence record annotated with the NHC-Ontology.
![alt text](https://github.com/LINNAE-project/SFB-Annotator/blob/iss22/data/rdf/rdf_examples_low.png)

### Example 1_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Pteropus minimus_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Taxon` |  The handwritten text contains a taxon name. | Should be a drop-down menu with all possible classes
| verbatim text | _Pteropus minimus_ | The verbatim text as written in the bounding box |
| language | _la_ | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  add autocomplete/drop down menu
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a taxon. |
| type | `taxon` |  auto-fill from entity type |
| Belongs to taxon | `https://identifiers.org/taxonomy:9397` | refers to _Chiroptera_ (order) in NCBI Taxonomy. Choice should be made whether this can be an IRI from an external database, or whether it should point to an already annotated verbatim entity. (such as _Chiroptera_ on page 2) | **question:** who do we allow to annotate, and who to interpret the data
| Taxon rank | `scientificName` | from list (`kingdom`,`class`,`order`,`family`,`genus`, `specificEpithet`, `scientificName`) | should be drop-down menu

\* **_nc_** refers to our namespace (natural committee for research in the Netherlands' Indies)

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|0|

!FIXME!

### Example 1_2

Variation: `dsw:hasIdentification`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Pteropus minimus_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Taxon` | The handwritten text contains a taxon name. | Should be a drop-down menu with all possible classes
| verbatim text | _Pteropus minimus_ | The verbatim text as written in the bounding box |
| language | _la_ | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ | add autocomplete/drop down menu
| Select property | `Organism identification to` | The annotation of an entity that refers to the identification of an organism (a taxon name as the initialisation of a field observation record, usually a new genus or scientificName). |
| type | `taxon` | auto-fill from entity type |
| Belongs to taxon | `https://identifiers.org/taxonomy:9397` | refers to _Chiroptera_ (order) in NCBI Taxonomy. Choice should be made whether this can be an IRI from an external database, or whether it should point to an already annotated verbatim entity. (such as _Chiroptera_ on page 2) | **question:** who do we allow to annotate, and who to interpret the data
| Taxon rank | `scientificName` | from list (`kingdom`,`class`,`order`,`family`,`genus`, `specificEpithet`, `scientificName`) | should be drop-down menu
| Identified by | `https://viaf.org/viaf/45106482/` | Corresponds to meaning of [`dwc:identifiedBy`](`https://dwc.tdwg.org/terms/#dwc:identifiedBy`) (people, groups, or organizations who assigned the Taxon to the subject). Preferably enter an instance of the class `foaf:Person` (or organization). Here, the writer of the field book is known (_Heinrich Kuhl_), but leave empty if unknown. | should be retrievable from pre-populated triple store
| Organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically, used to link all information belonging to an observation record of a single organism occurrence together. `1` is a placeholder; an unique ID should be generated for the organism record.\* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?

\* In the absence of a persistent global unique identifier, construct one from a combination of identifiers in the record that will most closely make the ID globally unique. Examples: http://arctos.database.museum/guid/MSB:Mamm:233627, 000866d2-c177-4648-a200-ead4007051b9, urn:catalog:UWBM:Bird:89776

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|55|

See [`example_1_2.ttl`](/data/rdf/example_1_2.ttl) file.

**Remarks:**

**Deleting annotations:**
* Delete the annotated box (page 2) via GUI.

* Check if the annotations (triples) were deleted from the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o
}
```
|?n|
|--|
|35|

!FIXME!

See [`example_1_del.ttl`](/data/rdf/example_1_del.ttl) file.

* Delete triples manually from the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/update).

```
delete where {
  ?s ?p ?o .
}
```

### Example 1_3

Variation: `nhc:additionalIdentification`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 8

* Annotate -> Click and Drag -> draw a bounding box around the written text _Pteropus_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Taxon` | The handwritten text contains a taxon name. |
| verbatim text | _Pteropus_ | The verbatim text as written in the bounding box |  
| language | _la_ | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ | add autocomplete/drop down menu
| Select property | `Additional identification to` | The annotation of an entity that refers to an _additional_ (secondary) identification of an organism (a scientific name as the initialisation of a field observation record). Usually, this is a second identification performed at a later stage (with availability of more knowledge for identification). Within this field note, the organism is first identified `in the field' as a new taxon _Gymnonotus_. |
| type | `taxon` | auto-fill from entity type |
| Belongs to taxon | `https://identifiers.org/taxonomy:9397` | refers to _Chiroptera_ (order) in NCBI Taxonomy. Choice should be made whether this can be an IRI from an external database, or whether it should point to an already annotated verbatim entity. (such as _Chiroptera_ on page 2) (question: who do we allow to annotate, and who to interpret the data) |
| Taxon rank | `genus` | from list (`kingdom`,`class`,`order`,`family`,`genus`, `specificEpithet`, `scientificName`)|
| Identified by | | Corresponds to meaning of [`dwc:identifiedBy`](`https://dwc.tdwg.org/terms/#dwc:identifiedBy`) (people, groups, or organizations who assigned the Taxon to the subject). Preferably enter an instance of the class `foaf:Person` (or organization). _Pteropus_ is most likely added at a later stage; therefore, we do not know (for sure) who performed the identification.  |
| Organism ID | `2` | Matches to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `2` is a placeholder; an unique ID should be generated for the organism record.\*| This additional identification should refer to the same organism ID as the first identification (`dsw:hasIdentification`). **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?
| Occurrence ID | `1` | Maps to [`dwc:occurrenceID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; some unique ID should be generated for the organism record.\* | Should potentially be removed, or should be the same as the occurrence ID generated for the first identification, as it is based on information from the same organism occurrence. Right now, the pop-up table for example_1_2 does not prompt for occurrenceID.


* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|0|

!FIXME!

**Remarks:** Triples not generated. Output should be similar to output of [`example_1_2.ttl`](/data/rdf/example_1_2.ttl), but `dsw:hasIdentification` should be replaced with `nhc:additionalIdentification`, and the additional identification should refer to the same organism ID as the first identification (`dsw:hasIdentification`).

### Example 2_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Geoff_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Person` |  The handwritten text contains a person name. | Should be a drop-down menu with all possible classes
| verbatim text | _Geoff_ | The verbatim text as written in the bounding box |
| language |  | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  not relevant
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a person. |
| type | `person` |  auto-fill from entity type |
| instance | `http://viaf.org/viaf/39377694` | Link the bounding box to the IRI if known. Here the person Étienne Geoffroy-Saint-Hilaire. | Preferably these can be retrieved with semantic autocomplete

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|20|

* See [`example_2_1.ttl`](/data/rdf/example_2_1.ttl) file.


### Example 2_2

Variation `nhc:scientificNameAuthorship`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Geoff_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Person` | The handwritten text contains a person name. | Should be a drop-down menu with all possible classes
| verbatim text | _Geoff_ | The verbatim text as written in the bounding box |
| language | | [ISO code](https://www.iso.org/iso-639-language-codes.html), here not relevant |can be removed
| Select property | `Author of scientific name` | The annotation of an entity that refers to a person that was the (published) author of a scientific name. Maps to [`dwc:scientificNameAuthorship`](https://dwc.tdwg.org/terms/#dwc:scientificNameAuthorship) The author is a crucial part of the taxon. Preferably, here, it refers to the persistent IRI of a person, such as `http://viaf.org/viaf/39377694` |
| type | `taxon` | auto-fill from entity type |
| Belongs to taxon | `nc:taxon1` | This field should contain the IRI that was generated when annotating the annotated scientific name, the subject of [`dwc:scientificNameAuthorship`](https://dwc.tdwg.org/terms/#dwc:scientificNameAuthorship) (here the annotated text: _Pteropus minimus_ on page 3) | **important note** The way it is setup now requires prior annotation of the taxon to which this name belongs. We should consider whether this makes sense.
| viaf IRI | `http://viaf.org/viaf/39377694` | the persistent identifier for the author of the scientific name. |  `nc:taxon1` [`dwc:scientificNameAuthorship`](https://dwc.tdwg.org/terms/#dwc:scientificNameAuthorship) `http://viaf.org/viaf/39377694`. Prompt should be changed to, e.g., person IRI instead of viaf IRI, as it could also be an orcid or other persistent identifier.

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|0|

!FIXME!

### Example 2_3

Variation: `dwciri:identifiedBy`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Diard_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Person` | The handwritten text contains a person name. | Should be a drop-down menu with all possible classes
| verbatim text | _Diard_ | The verbatim text as written in the bounding box |
| language | | [ISO code](https://www.iso.org/iso-639-language-codes.html) | can be removed
| Select property | `Organism identified by` | Maps to `dwciri:identifiedBy`
| type | `taxon` | auto-fill from entity type |
| Organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically, used to link all information belonging to an observation record of a single organism occurrence together. `1` is a placeholder; an unique ID should be generated for the organism record.\* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?
| viaf IRI | `http://viaf.org/viaf/39377694` | the persistent identifier for the person. | Prompt should be changed to, e.g., person IRI instead of viaf IRI, as it could also be an orcid or other persistent identifier.

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|21\*|

* See [`example_2_3.ttl`](/data/rdf/example_2_3.ttl) file.

### Example 2_4

Variation: `dwciri:recordedBy`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Diard_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Person` | The handwritten text contains a person name. |  Should be a drop-down menu with all possible classes
| verbatim text | _Diard_ | The verbatim text as written in the bounding box |
| language | | [ISO code](https://www.iso.org/iso-639-language-codes.html) | can be removed |
| Select property | `Occurrence recorded by` | The annotation of an entity that refers to a person that recorded the occurrence of the organism. Maps to `dwciri:recordedBy`|
| type | `taxon` | auto-fill from entity type |
| viaf IRI | `http://viaf.org/viaf/39377694` | the persistent identifier for the person. |  Prompt should be changed to, e.g., person IRI instead of viaf IRI, as it could also be an orcid or other persistent identifier.

* See [`example_2_4.ttl`](/data/rdf/example_2_4.ttl) file.

### Example 3_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Buitenzorg_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Location` |  The handwritten text contains a location name. | Should be a drop-down menu with all possible classes
| verbatim text | _Buitenzorg_ | The verbatim text as written in the bounding box |
| language |  nl | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a location. |
| type | `location` |  auto-fill from entity type |
| instance | `http://sws.geonames.org/1648473/` | Link the bounding box to the IRI if known. Here the location _Buitenzorg_, currently called _Bogor_. | Preferably these instances can be retrieved with semantic autocomplete


* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|20|

* See [`example_3_1.ttl`](/data/rdf/example_3_1.ttl) file.

### Example 3_2

Variation: `dsw:locatedAt`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Bagalonga_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `Location` | The handwritten text contains a location name. | Should be a drop-down menu with all possible classes|  
| verbatim text | _Bagalonga_ | The verbatim text as written in the bounding box |
| language | | [ISO code](https://www.iso.org/iso-639-language-codes.html)  |add autocomplete/drop down menu
| Select property | `Occurrence located at` | The annotation of an entity that refers to the location where the observation of the organism occurrence took place. Maps to `dsw:locatedAt`. The prompt `Additional occurrence located at` should be merged with this one. |
| type | `taxon` | auto-fill from entity type |
| gn:geonamesfeature IRI |  | if known, the persistent identifier for the location from the geonames ontology, e.g., `http://sws.geonames.org/1648473/` (stands for _Buitenzorg_, also _Bogor_) |
| organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; an unique ID should be generated for the organism record.* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?

* See [`example_3_2.ttl`](/data/rdf/example_3_2.ttl) file.

### Example 4_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the table

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `MeasurementOrFact` |  The handwritten text contains a measurement or fact. | Should be a drop-down menu with all possible classes
| verbatim text |  | leave empty |
| language |   | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ | should be removed here
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a measurement or fact. |
| type | `measurementorfact` |  auto-fill from entity type |
| instance |  | Should be left empty | This field should be removed here.


* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|22|

* See [`example_4_1.ttl`](/data/rdf/example_4_1.ttl) file.

### Example 4_2

Variation: `dsw:derivedFrom`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the measurement or fact (such as a table or a statement about the animal)

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|----
| Entity type | `MeasurementOrFact` | The handwritten text contains a table or other measurement or fact. | Should be a drop-down menu with all possible classes|  
| verbatim text | leave empty |  |  |
| language | leave empty | [ISO code](https://www.iso.org/iso-639-language-codes.html) | add autocomplete/drop down menu |
| Select property | `Identification based on (table)` |  should be table _or_ paragraph| |
| type | `measurementorfact` | auto-fill from entity type |
| organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; an unique ID should be generated for the organism record.* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?  |

See [`example_4_2.ttl`](/data/rdf/example_4_2.ttl) file.

### Example 5_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Color_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `PropertyOrAttribute` |  The handwritten text contains a property or attribute name. | Should be a drop-down menu with all possible classes
| verbatim text | _Color_ | The verbatim text as written in the bounding box |
| language |  nl | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a property or attribute name. |
| type | `propertyOrAttribute` |  auto-fill from entity type |
| instance | `http://identifiers.org/ncit/C37927` | Link the bounding box to the IRI if known. Here the property _Color_. | Preferably these instances can be retrieved with semantic autocomplete

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|26|

* See [`example_5_1.ttl`](/data/rdf/example_5_1.ttl) file.


### Example 5_2

Variation: `nhc:measuresOrDescribes` a `propertyOrAttribute`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Color_

* Fill in the pop-up form/table:

|Key|Value |Notes  | Diff Expected
|---|-----|-----|-----
| Entity type | `propertyOrAttribute` | The handwritten text contains an indication of a property or attribute of the described organism. | Should be a drop-down menu with all possible classes|  
| verbatim text | _Color_ |  |  |
| language | la | [ISO code](https://www.iso.org/iso-639-language-codes.html)  |add autocomplete/drop down menu
| Select property | `Table/paragraph measures or describes` | Indicating that a table or a paragraph describes the certain property or attribute, here the _Color_ of the animal | |
| type | `propertyOrAttribute` | auto-fill from entity type |
| `ncit:propertyorattribute subclass IRI` | `http://identifiers.org/ncit/C37927` | A subclass of `http://identifiers.org/ncit/C20189` or `propertyOrAttribute` from the NCIT ontology.
| organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; an unique ID should be generated for the organism record.* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?  |

See [`example_5_2.ttl`](/data/rdf/example_5_2.ttl) file.

### Example 6_1

Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Dentibus_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `AnatomicalEntity` |  The handwritten text contains the name of an anatomical entity. | Should be a drop-down menu with all possible classes
| verbatim text | _Dentibus_ | The verbatim text as written in the bounding box |
| language |  nl | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains the name of an anatomical entity. |
| type | `anatomicalentity` |  auto-fill from entity type |
| instance | `http://purl.obolibrary.org/obo/UBERON_0003672` | Link the bounding box to the IRI if known. Here the anatomical entity _Dentibus_. | Preferably these instances can be retrieved with semantic autocomplete

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|26|

* See [`example_6_1.ttl`](/data/rdf/example_6_1.ttl) file.

### Example 6_2

Variation: `nhc:measuresOrDescribes` an `anatomicalEntity`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Dentibus_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `anatomicalEntity` | The handwritten text contains the name of an anatomical entity of the described organism. | Should be a drop-down menu with all possible classes|  
| verbatim text | _Dentibus_ |  | |
| language | la | [ISO code](https://www.iso.org/iso-639-language-codes.html)   | add autocomplete/drop down menu
| Select property | `Table/paragraph measures or describes` | Indicating that a table or a paragraph describes the certain anatomical entity, here the _Dentibus_ (teeth) of the animal | |
| type | `anatomicalEntity` | auto-fill from entity type |
| `uberon:anatomicalentity subclass IRI` | `http://purl.obolibrary.org/obo/UBERON_0003672` | A subclass of `http://purl.obolibrary.org/obo/UBERON_0001062` or `anatomical entity` from the uberon ontology.
| organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; an unique ID should be generated for the organism record.* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend? |

See [`example_6_2.ttl`](/data/rdf/example_6_2.ttl) file.

### Example 7_1


Variation: `rdf:type`

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 4

* Annotate -> Click and Drag -> draw a bounding box around the handwritten words _10 April 1821_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `date` |  The handwritten text contains a date. | Should be a drop-down menu with all possible classes
| verbatim text | _10 April 1821_ | The verbatim text as written in the bounding box |
| language |  nl | [ISO code](https://www.iso.org/iso-639-language-codes.html) for _latin_ |  
| Select property | `Type`  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a date. |
| type | `date` |  auto-fill from entity type |
| instance | | this field should be removed | instead of the instance field, we would like to see the same fields as when the property `Organism described on` is selected, minus the organism ID field. See example [`example 7_2`](#Example-7_2) and corresponding ttl file [`example_7_2.ttl`](/data/rdf/example_7_2.ttl) file.

* Check generated triples in the [RDF store](http://localhost:8080/rdf4j-workbench/repositories/mem-rdf/query).

```
select (count(*) as ?n)
where {
  ?s ?p ?o .
}
```
|?n|
|--|
|17|

See [`example_7_1.ttl`](/data/rdf/example_7_1.ttl) file.


### Example 7_2

Variation: `nhc:verbatimDate` (and `rdf:type`\*)

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 4

* Annotate -> Click and Drag -> draw a bounding box around the handwritten date _10 april 1821_

* Fill in the pop-up form/table:

|Key|Value |Notes | Diff Expected
|---|-----|-----|-----
| Entity type | `date` | The handwritten text contains a date, either a day, month year, or a combination thereof. | Should be a drop-down menu with all possible classes|  
| verbatim text | _10 april 1821_ |  Text as is, no formatting| |
| language | nl | [ISO code](https://www.iso.org/iso-639-language-codes.html)   | add autocomplete/drop down menu
| Select property | `Organism described on` | Indicates that the observation of the animal was on a certain date |  |
| type | `date` | auto-fill from entity type |
| Year (yyyy) | 1821 || should not be auto-filled with a 0
| Month (mm) | 04 | |should not be auto-filled with a 0
| Day (dd) | 10 | |should not be auto-filled with a 0
| organism ID | `1` | Maps to [`dwc:organismID`](https://dwc.tdwg.org/terms/#dwc:occurrenceID). Should be generated automatically. `1` is a placeholder; an unique ID should be generated for the organism record.* | **important note**: when is this organism and possibly also occurrence ID generated? It allows all information belonging to the occurrence or organism observation to be linked together. Hence, how will we enforce this in the backend?  |

See [`example_7_2.ttl`](/data/rdf/example_7_2.ttl) file.

\* For the `rdf:type` variant of `nhc:Date`, the fields should be the same, only omitting the field `organism ID`.

**General notes:**

* (back-end) Most annotations (except for variety 1, a class annotation) link to an annotation record, and atm, some annotations therefore depend on the prior annotation of other named entities. However, these annotation events should be able to occur un any order and the code should reflect this.  

* (relates to front-end) The entity type field should produce a drop-down menu with the possible classes  (`taxon`, `person`, `location`, `date`, `anatomical entity`, `measurementOrFact`) rather than a free-text field. Similarly, the entry to this field should limit the number of possibilities for the `select property` field. (E.g., `person` -> `Type`, `Author of scientific name`, `Organism identified by`, `Occurrence recorded by`/`Additional Occurrence recorded by`)
| dc:identifier	|	dcmitype:Image	| foaf:Person	| foaf:firstName	| foaf:lastName	| rdfs:comment |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | 
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_001.jpg	 |	https://viaf.org/viaf/45106482/	|Heinrich	|Kuhl	|Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata|
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_002.jpg |	https://viaf.org/viaf/45106482/	| Heinrich | 	Kuhl	| Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_003.jpg |	https://viaf.org/viaf/45106482/ |	Heinrich | Kuhl |	Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_004.jpg	 |	https://viaf.org/viaf/45106482/ |	Heinrich |	Kuhl |	Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_005.jpg |	https://viaf.org/viaf/45106482/ |	Heinrich	|Kuhl	| Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_006.jpg |	https://viaf.org/viaf/45106482/	| Heinrich |	Kuhl	| Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_007.jpg|	https://viaf.org/viaf/45106482/ |	Heinrich	| Kuhl |	Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_008.jpg |	https://viaf.org/viaf/45106482/ |	Heinrich |	Kuhl |	Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |
|NC_a_Ku_012_002|	MMNAT01_PM_NNM001001033_009.jpg  |	https://viaf.org/viaf/45106482/	 |Heinrich |	Kuhl | Bundel aantekeningen: Mammalien. Bevat mapjes met de volgende titels: Chiropterae/Lemur volitans sijn. Galeopithecus variegatus Aud./Quadrumanae/Falculatae/Glires/Rumantia/Nantes/Edentuli/Anatomie des Elephas indicus fem. Buitenzorg 1821 d.16 februar/Pachydermata |



![alt text](https://i.creativecommons.org/p/mark/1.0/88x31.png "Creative Commons License")

This work is licensed under a [Creative Commons Public Domain Mark 1.0 License](http://creativecommons.org/publicdomain/mark/1.0/).
### Use cases
Field notes are collections of observation records ([`dwc:HumanObservation`](https://dwc.tdwg.org/terms/#humanobservation)) that describe the occurrence ([`dwc:Occurrence`](https://dwc.tdwg.org/terms/#occurrence)) of an organism ([`dwc:Organism`](http://rs.tdwg.org/dwc/terms/Organism)) at some location ([`dwc:Location`](https://dwc.tdwg.org/terms/#location)) during some time ([`dwc:Event`](https://dwc.tdwg.org/terms/#event)). An observation record most often includes an identification ([`dwc:Identification`](https://dwc.tdwg.org/terms/#identification)) to a taxon ([`dwc:Taxon`](https://dwc.tdwg.org/terms/#taxon)), and accompanying measurements and facts ([`dwc:MeasurementOrFact`](https://dwc.tdwg.org/terms/#measurementorfact)) that were used for identification. These observation records usually span multiple field notes, but can also be very short utterances such as _Sci. Diard Buitenzorg_ (on [page 3](data/jpg/MMNAT01_AF_NNM001001033_003.jpg)), most likely referring to the occurrence of another organism with the same identification.

**Table 1.** List of examples with annotation classes and properties.
| class                                                                                                                                                                             | property   | example                                                         | revised            |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | --------------------------------------------------------------- | ------------------ |
| [`dwc:Taxon`](https://dwc.tdwg.org/terms/#taxon)                                                                                                                                  | `rdf:type` | [`1_1`](#Example-1_1)                                           | :heavy_check_mark: |
| [`foaf:Person`](http://xmlns.com/foaf/spec/#term_Person)                                                                                                                          | `rdf:type` | [`2_1`](#Example-2_1)                                           | :heavy_check_mark: |
| [`dwc:Location`](https://dwc.tdwg.org/terms/#location), [`dcterms:Location`](https://www.dublincore.org/specifications/dublin-core/dcmi-terms/#http://purl.org/dc/terms/Location) | `rdf:type` | [`3_1`](#Example-3_1)                                           | :heavy_check_mark: |
| [`dwc:MeasurementOrFact`](https://dwc.tdwg.org/terms/#measurementorfact)                                                                                                          | `rdf:type` | [`4_1`](#Example-4_1)[`4_2`](#Example-4_2)[`4_3`](#Example-4_3) | :heavy_check_mark: |
| [`dwc:Event`](https://dwc.tdwg.org/terms/#event)                                                                                                                                  | `rdf:type` | [`5_1`](#Example-5_1)                                           | :heavy_check_mark: |


### Example 1_1

Property: `rdf:type`

Input: [JSON](/data/json/local/example_1_1.json)

Output: [RDF/Turtle](/data/rdf/local/example_1_1.ttl) | [SVG](/doc/models/example_1_1.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Pteropus minimus_

* Fill in the pop-up form/table:

| Key                               | Value                                                      | Notes                                                                                                                                                                                                                                                                                                                                                                           | Diff Expected                                                          |
| --------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| Entity type/`text`                | `Taxon`                                                    | The handwritten text contains a taxon name.                                                                                                                                                                                                                                                                                                                                     | Should be a drop-down menu with all possible classes                   |
| verbatim text/`verbatim`          | `Pteropus minimus`                                         | The verbatim text as written in the bounding box                                                                                                                                                                                                                                                                                                                                |
| language/`language`               | `lat`                                                      | [ISO 639-3 code](https://iso639-3.sil.org/code/lat) for _Latin_                                                                                                                                                                                                                                                                                                                 | add autocomplete/drop down menu                                        |
| Select property/`property`        | `Type`                                                     | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a taxon.                                                                                                                                                                                                                                                   |
| type/`type`                       | `taxon`                                                    | auto-fill according to the Entity type                                                                                                                                                                                                                                                                                                                                          |
| Belongs to taxon/`belongstotaxon` | `http://www.gbif.org/species/9180402`                      | Point to an IRI, if known (e.g., using GBIF or NCBI Taxonomy), or to an exisiting (verbatim) annotation such as _Chiroptera_ (on [page 2](data/jpg/MMNAT01_AF_NNM001001033_002.jpg)).                                                                                                                                                                                           | A distinction should be made between transcription and interpretation. |
| Taxon rank/`rank`                 | `http://purl.obolibrary.org/obo/TAXRANK_0000006` (species) | from list ([`kingdom`](http://purl.obolibrary.org/obo/TAXRANK_0000017),[`class`](http://purl.obolibrary.org/obo/TAXRANK_0000002),[`order`](http://purl.obolibrary.org/obo/TAXRANK_0000003),[`family`](http://purl.obolibrary.org/obo/TAXRANK_0000004),[`genus`](http://purl.obolibrary.org/obo/TAXRANK_0000005) or [`species`](http://purl.obolibrary.org/obo/TAXRANK_0000006)) | should be drop-down menu according to the taxon IRI above              |


### Example 2_1

Property: `rdf:type`

Input: [JSON](/data/json/local/example_2_1.json)

Output: [RDF/Turtle](/data/rdf/local/example_2_1.ttl) | [SVG](/doc/models/example_2_1.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Geoff_

* Fill in the pop-up form/table:

| Key                        | Value                           | Notes                                                                                                                          | Diff Expected                                        |
| -------------------------- | ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------- |
| Entity type/`text`         | `Person`                        | The handwritten text contains a person name.                                                                                   | Should be a drop-down menu with all possible classes |
| verbatim text/`verbatim`   | `Geoff`                         | The verbatim text as written in the bounding box                                                                               |
| language                   |                                 | NA (default: und)                                                                                                              |
| Select property/`property` | `Type`                          | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a person. |
| type/`type`                | `person`                        | auto-fill from entity type                                                                                                     |
| Instance/`instance`        | `http://viaf.org/viaf/39377694` | Add IRI if known (optional). For example, it refers to the person _Étienne Geoffroy-Saint-Hilaire_.                            |                                                      |


### Example 3_1

Property: `rdf:type`

Input: [JSON](/data/json/local/example_3_1.json)

Output: [RDF/Turtle](/data/rdf/local/example_3_1.ttl) | [SVG](/doc/models/example_3_1.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the written text _Buitenzorg_

* Fill in the pop-up form/table:

| Key             | Value                             | Notes                                                                                                                            | Diff Expected                                        |
| --------------- | --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Entity type     | `Location`                        | The handwritten text contains a location name.                                                                                   | Should be a drop-down menu with all possible classes |
| verbatim text   | `Buitenzorg`                      | The verbatim text as written in the bounding box                                                                                 |
| language        | `dut`                             | [ISO 639-3 code](https://iso639-3.sil.org/code/latl) for _Dutch_                                                                 |
| Select property | `Type`                            | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a location. |
| type            | `location`                        | auto-fill from entity type                                                                                                       |
| Instance        | `http://sws.geonames.org/1648473` | Add IRI if known (optional). For example, the location _Buitenzorg_ currently refers to _Bogor_.                                 |


### Example 4_1

Property: `rdf:type`

Input: [JSON](/data/json/local/example_4_1.json)

Output: [RDF/Turtle](/data/rdf/local/example_4_1.ttl) | [SVG](/doc/models/example_4_1.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the table

* Fill in the pop-up form/table:

| Key             | Value                                                                                                                                   | Notes                                                                                                                                       | Diff Expected                                        |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Entity type     | `MeasurementOrFact`                                                                                                                     | The handwritten text contains a measurement or fact.                                                                                        | Should be a drop-down menu with all possible classes |
| verbatim text   | `- Longitud. tota --- 1,0,3\n- corporis --- 0,2,9\n- tota --- 0,2,9\n- capitis --- 0,1,1\n- caudae --- 0,0,1\n- digit. medii --- 0,3,0` | [Markdown](https://www.markdownguide.org/) list or table                                                                                    |
| language        | `lat`                                                                                                                                   | [ISO 639-3 code](https://iso639-3.sil.org/code/lat) for _Latin_                                                                             |
| Select property | `Type`                                                                                                                                  | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a measurement or fact. |
| type            | `measurementorfact`                                                                                                                     | auto-fill from entity type                                                                                                                  |
| Instance        |                                                                                                                                         | Add IRI if known (optional).                                                                                                                |


### Example 4_2

Property: `rdf:type`

Input: [JSON](/data/json/local/example_4_2.json)

Output: [RDF/Turtle](/data/rdf/local/example_4_2.ttl) | [SVG](/doc/models/example_4_2.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Color_

* Fill in the pop-up form/table:

| Key             | Value                                                 | Notes                                                                                                                                              | Diff Expected                                        |
| --------------- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Entity type     | `MeasurementOrFact`                                   | The handwritten text contains a property or attribute name.                                                                                        | Should be a drop-down menu with all possible classes |
| verbatim text   | `Color`                                               | The verbatim text as written in the bounding box                                                                                                   |
| language        | `lat`                                                 | [ISO 639-3 code](https://iso639-3.sil.org/code/lat) for _Latin_                                                                                    |
| Select property | `Type`                                                | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a property or attribute name. |
| type            | `measurementorfact`                                   | auto-fill from entity type                                                                                                                         |
| Instance        | `http://purl.obolibrary.org/obo/PATO_0000014` (color) | Add IRI if known (optional).                                                                                                                       |


### Example 4_3

Property: `rdf:type`

Input: [JSON](/data/json/local/example_4_3.json)

Output: [RDF/Turtle](/data/rdf/local/example_4_3.ttl) | [SVG](/doc/models/example_4_3.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 3

* Annotate -> Click and Drag -> draw a bounding box around the handwritten word _Dentibus_

* Fill in the pop-up form/table:

| Key             | Value                                                       | Notes                                                                                                                                                  | Diff Expected                                        |
| --------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------- |
| Entity type     | `MeasurementOrFact`                                         | The handwritten text contains the name of an anatomical entity.                                                                                        | Should be a drop-down menu with all possible classes |
| verbatim text   | `Dentibus`                                                  | The verbatim text as written in the bounding box                                                                                                       |
| language        | `lat`                                                       | [ISO 639-3 code](https://iso639-3.sil.org/code/lat) for _Latin_                                                                                        |
| Select property | `Type`                                                      | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains the name of an anatomical entity. |
| type            | `measurementorfact`                                         | auto-fill from entity type                                                                                                                             |
| Instance        | `http://purl.obolibrary.org/obo/UBERON_0001062` (dentition) | Add IRI if known (optional).                                                                                                                           |


### Example 5_1

Property: `rdf:type`

Input: [JSON](/data/json/local/example_5_1.json)

Output: [RDF/Turtle](/data/rdf/local/example_5_1.ttl) | [SVG](/doc/models/example_5_1.svg)

Steps:

* Go to http://localhost:8080/semanticAnnotator/

* Register -> Save

* Collections -> manuscript MMNAT01_AF -> page 4

* Annotate -> Click and Drag -> draw a bounding box around the handwritten words _10 April 1821_

* Fill in the pop-up form/table:

| Key             | Value           | Notes                                                                                                                        | Diff Expected                                        |
| --------------- | --------------- | ---------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Entity type     | `date`          | The handwritten text contains a date.                                                                                        | Should be a drop-down menu with all possible classes |
| verbatim text   | `10 April 1821` | The verbatim text as written in the bounding box                                                                             |
| language        | `ger`           | [ISO 639-3 code](https://iso639-3.sil.org/code/ger) for _German_                                                             |
| Select property | `Type`          | This refers to an entity annotation without further interpretation, merely specifying that the bounding box contains a date. |
| type            | `date`          | auto-fill from entity type                                                                                                   |
| Instance        | `1821-04-10`    | [ISO 8601](https://www.w3.org/TR/NOTE-datetime) date/time format (`YYYY-MM-DD`)                                              |
