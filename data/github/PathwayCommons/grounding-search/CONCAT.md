---
title: 'A flexible search system for high-accuracy identification of biological entities and molecules'
tags:
  - biological entities
  - small molecules
  - grounding
  - identifiers
  - identification
  - search
  - JavaScript
  - Docker
authors:
  - name: Max Franz
    orcid: 0000-0003-0169-0480
    affiliation: 1
  - name: Jeffrey V. Wong
    affiliation: 1
    orcid: 0000-0002-8912-5699
  - name: Metin Can Siper
    orcid: 0000-0002-7556-093X
    affiliation: 2
  - name: Christian Dallago
    orcid: 0000-0003-4650-6181
    affiliation: "3, 4, 5"
  - name: John Giorgi
    orcid: 0000-0001-9621-5046
    affiliation: 1
  - name: Emek Demir
    orcid: 0000-0002-3663-7113
    affiliation: 2
  - name: Chris Sander
    orcid: 0000-0001-6059-6270
    affiliation: "3, 6, 7"
  - name: Gary D. Bader^[Corresponding author]
    orcid: 0000-0003-0185-8861
    affiliation: "1, 8, 9, 10, 11"
affiliations:
 - name: The Donnelly Centre, University of Toronto, Toronto, Ontario, M5S 3E1, Canada
   index: 1
 - name: Computational Biology Program, Oregon Health and Science University, Portland, OR 97239, USA
   index: 2
 - name: Department of Cell Biology, Harvard Medical School, Boston, MA, 02215, USA
   index: 3
 - name: Department of Systems Biology, Harvard Medical School, Boston, MA, 02215, USA
   index: 4
 - name: Department of Informatics, Technische Universität München, 85748 Garching, Germany
   index: 5
 - name: Department of Data Sciences, Dana-Farber Cancer Institute, Boston, MA, 02215, USA
   index: 6
 - name: Broad Institute of MIT and Harvard, Boston, MA, 02142, USA
   index: 7
 - name: Department of Computer Science, University of Toronto, Ontario, M5S 2E4, Canada
   index: 8
 - name: Department of Molecular Genetics, University of Toronto, Ontario, M5S 1A8, Canada
   index: 9
 - name: The Lunenfeld-Tanenbaum Research Institute, Mount Sinai Hospital, Toronto, Ontario, M5G 1X5, Canada
   index: 10
 - name: Princess Margaret Cancer Centre, University Health Network, Toronto, Ontario, M5G 2C1, Canada
   index: 11
date: 9 September 2021
bibliography: paper.bib
---

# Summary

Identifying subcellular biological entities (genes, gene products, and small molecules) is essential in using and creating bioinformatics analysis tools, text mining, and accessible biological research apps. When research information is uniquely and unambiguously identified, it enables data to be accurately retrieved, cross-referenced, and integrated. In practice, biological entities are identified when they are associated with a matching record from a knowledge base that specialises in collecting and organising information of that type (e.g. genes in NCBI Gene). Our search service increases the efficiency and ease of use for identifying biological entities compared to prior approaches [@reimand2007g; @cote2007protein; @berriz2008synergizer; @juty2012identifiers]. A comparison of similar tools is available in the software documentation[@githubRepoReadme]. This identification service may be used to power research apps and tools, such as Biofactoid, GeneMANIA, and STRING, where colloquial entity names may be provided as input by users [@wong2021capturing; @mostafavi2008genemania; @mering2003string].

# Statement of need

Most biologists are unaware of the concept of grounding data to database identifiers (i.e., a value, usually a piece of text, that uniquely identifies the entity) [@abeysooriya2021gene]. When an author labels a biological entity as ‘IL6’, for instance, they may not consider that this label could be ambiguous. Is this IL6 for *Homo sapiens* or for *Mus musculus*? If for *Mus musculus*, there is more than one gene that is called by that name. By mapping the user’s entity to a database identifier, in this case, 3569 in NCBI Gene, the data becomes disambiguated.

Because many biologists are unaware of the utility of database identifiers, they often have the perception that the common or canonical name of an entity is sufficient for use with analysis tools. These users can be confused by traditional database grounding interfaces, where it is the user’s responsibility to select a particular identifier from a long list manually. This sort of grounding interface is incongruent with users’ mental model: What purpose could this list of entities have when the entity has already been identified by name?

Our search service can be used to power interfaces that allow the user to identify an entity by name, as per their mental model. This makes grounding accessible to a wider set of researchers, with high ease of use. The service returns its results quickly, i.e., in less than 100 milliseconds on a 2.4 GHz dual-core processor with 8 GB of RAM, so that the result may be shown to a user interactively and without impeding the user’s actions. As output, the service returns a ranked list of possible groundings in descending order of relevance. The first entry in this list is the predicted identifier with the highest confidence given to the user’s input. The remainder of the list exists only to allow the user to recover from an incorrect first identifier. The service is customisable to accommodate various use cases, e.g., interactive grounding interfaces that can dynamically build up heuristics during a user’s session with the system.

# Performance evaluation

In the rare case where the service does not return the correct result as the first result, the app which uses the service may present one or more of the following entries in the list. For instance, the app may allow the user to indicate a correction manually — e.g., “I meant this IL6, not that one”. To facilitate this, each identifier in the results list has included a corresponding set of descriptions and metadata.

To verify the accuracy of the service’s results, a test suite was created. The tests include entity names used in PubMed Central, Pathway Commons, and Biofactoid projects. Primarily, the test cases were prioritised based on popularity to measure how well commonly-researched entities are correctly assigned by the service. In all, there are currently over 750 test cases. Of those test cases, 91% returned the expected result as the first entry in the returned list in our tests. In nearly 98% of cases, the expected result was within the top ten entries in the returned list.

# Mechanism

The database is built using a dynamic indexing approach with Elasticsearch. The system includes facilities to download the latest set of identifiers and associated metadata from NCBI Gene (including DNA, genes, and proteins in *H. sapiens*, *M. musculus*, *S. cerevisiae*, *D. melanogaster*, *E. coli*, *C. elegans*, *A. thaliana*, *R. norvegicus*, *D. rerio*, *SARS-CoV-2*, and extensible to other organisms), ChEBI (including small molecules, e.g., drugs), and UniProt (including proteins for cross-referencing NCBI Gene identifiers) to build the database. This indexing process is exposed as a top-level command, integrated with the search server itself, so that data sources can automatically be kept fresh regularly. The indices of the grounding service can be exported to an external repository (e.g. Zenodo) and directly imported in order to provide a means of referencing a particular version of an index and to reduce the need to manually build an index for researchers that reuse the grounding service in their own projects.

The grounding service operates by making queries on the database, with the query string normalised -- i.e., the string has punctuation removed, the case is normalised (e.g., “TNF-a” is considered the same as “tnf alpha”). A fuzzy query and a precise query are made to the database to ensure that both exact matches and near matches are included in the initial results. These initial results are processed with a multithreaded ranking approach. The Sørensen–Dice coefficient is used to rank the initial results based on each entity’s official name and synonyms, considering each entity's best-case score.

The Sørensen–Dice coefficient, $s_i$, for the $i$th synonym, is given as follows, where:
* $b_t$ is the total number of bigrams in both the $i$th synonym and the query string,
* $b_i$ is the total number of bigrams in the $i$th synonym, and
* $b_q$ is the total number of bigrams in the query string:

$$s_i = (\frac{2b_t}{b_i + b_q})$$

For tie-breaking, a series of further measures are used to order the results. These measures include an organism ranking preference, molecular charge (with a preference for charged, aqueous molecules), and the Sørensen–Dice coefficient of the official name. Finally, the data source (namespace) filter may be applied, if specified, in order to include results only of a particular type (e.g., small molecules from ChEBI).


# Usage

With a complete index, the grounding service’s server can be started. The server exposes a REST-like API to enable client applications to query the index [@groundingswagger]. The main endpoint in the API is search. A search contains a number of input parameters, the user’s typed entity name chief among them. Other parameters are optional. They may be implicitly specified by the user, or they may be specified by the client application. The target namespace is an optional parameter that filters the search result by a particular data source (e.g., NCBI Gene). When unspecified, no namespace filtering is applied. A second optional parameter is an organism ranking, which can inform the service of the relative likelihood of a search’s pertinence to a particular organism. An organism ranking based on general popularity, measured by PubMed mentions, is used by default [@gene2pubmed].

# Discussion

An existing app, Biofactoid, provides an example of how this grounding service empowers novice users to ground entities to database identifiers automatically. Biofactoid users typically take one to five minutes to summarise the molecular interactions with a paper. The database identifiers in author’s Biofactoid documents have been accurate thus far, as evaluated by extensive manual spot-checking, and these novice users require no training in or understanding of database identifiers in order to create these documents.

Apps that previously required a user to explicitly specify an organism may instead use heuristics paired with the grounding service’s API to provide intelligent results without user intervention. On the other hand, organism-specific apps may leverage the organism ordering to provide results only for the relevant organism: The first instance of a non-specified organism in one of the returned entries indicates that there are no further grounding entries for the specified organism.

Further, more natural search systems may emerge as a result of the approaches of this grounding service. A biological data search engine may allow for natural language queries, powered by the grounding service. A user may type a search, such as “interactions of tnf in human”, in order to get relevant results. If the user types an ambiguous query, such as “tnf”, the search engine may use the ranked results to provide intelligent follow-up questions: “Did you mean tnf (mouse)?” These results and follow-up questions may go so far as to be user-personalised. A user that predominantly searches for mouse genes may have the search engine provide custom-tailored results, based on a search history used to inform the grounding service.

Modern research apps and tools motivate the need for robust, fast, reusable grounding tools that allow for easily identifying biological entities from their common names. Our grounding service can be used in apps to provide users with an easy-to-use experience in line with their mental model of biological entities.


# Acknowledgements

This project was funded by the US National Institutes of Health (NIH) [U41 HG006623, U41 HG003751, R01 HG009979 and P41 GM103504].

# References
# grounding-search

[![DOI](https://zenodo.org/badge/170363308.svg)](https://zenodo.org/badge/latestdoi/170363308)
[![status](https://joss.theoj.org/papers/dfad3b1bc874b0f813f6814723a646f5/status.svg)](https://joss.theoj.org/papers/dfad3b1bc874b0f813f6814723a646f5)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/PathwayCommons/grounding-search/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/PathwayCommons/grounding-search.svg?branch=master)](https://travis-ci.org/PathwayCommons/grounding-search)

## Description

The identification of sub-cellular biological entities is an important consideration in the use and creation of bioinformatics analysis tools and accessible biological research apps.  When research information is uniquely and unambiguously identified, it enables data to be accurately retrieved, cross-referenced, and integrated.  In practice, biological entities are “identified” when they are associated with a matching record from a knowledge base that specialises in collecting and organising information of that type (e.g. gene sequences).  Our search service increases the efficiency and ease of use for identifying biological entities.  This identification may be used to power research apps and tools where common entity synonyms may be provided as input.

For instance, [Biofactoid](https://biofactoid.org) uses this grounding service to allow users to simply specify their preferred synonyms to identify biological entities (e.g. proteins):

https://user-images.githubusercontent.com/989043/140164756-1fa22796-1c60-4f13-9a2a-d65393b89155.mp4

## Citation

To cite the Pathway Commons Grounding Search Service in a paper, please cite the Journal of Open Source Software paper:

Franz et al., (2021). A flexible search system for high-accuracy identification of biological entities and molecules. Journal of Open Source Software, 6(67), 3756, https://doi.org/10.21105/joss.03756

[View the paper at JOSS](https://joss.theoj.org/papers/10.21105/joss.03756) or [view the PDF directly](https://www.theoj.org/joss-papers/joss.03756/10.21105.joss.03756.pdf).

## Maintenance

The Pathway Commons Grounding Search Service is an academic project built and maintained by: 
<a href="https://baderlab.org" target="_blank">Bader Lab at the University of Toronto</a>
, 
<a href="http://sanderlab.org" target="_blank">Sander Lab at Harvard</a>
, and the 
<a href="https://www.ohsu.edu/people/emek-demir/AFE06DC89ED9AAF1634F77D11CCA24C3" target="_blank">Pathway and Omics Lab at the Oregon Health &amp; Science University</a>
.

## Funding

This project was funded by the US National Institutes of Health (NIH) [U41 HG006623, U41 HG003751, R01 HG009979 and P41 GM103504].

## Quick start

### Via Docker

Install [Docker](https://docs.docker.com/) (>=20.10.0) and [Docker Compose](https://docs.docker.com/compose/) (>=1.29.0).

Clone this remote or at least the `docker-compose.yml` file then run:

```
docker-compose up
```

Swagger documentation can be accessed at [`http://localhost:3000`](http://localhost:3000).

NB: Server start will take some time in order for Elasticsearch to initialize and for the grounding data to be retrieved and the index restored. If it takes more than 10 minutes consider increasing the allocated memory for Docker: `Preferences` > `Resources` > `Memory` and remove this line in docker-compose.yml: `ES_JAVA_OPTS=-Xms2g -Xmx2g` 

### Via source

With [Node.js](https://nodejs.org/en/) (>=8) and [Elasticsearch](https://www.elastic.co/products/elasticsearch) (>=6.6.0, <7) installed with default options, run the following in a cloned copy of the repository:

- `npm install`: Install npm dependencies
- `npm run update`: Download and index the data
- `npm start`: Start the server (by default on port 3000)

## Documentation

Swagger documentation is available on a publicly-hosted instance of the service at [https://grounding.baderlab.org](https://grounding.baderlab.org).  You can run queries to test the API on this instance.

Please do not use `https://grounding.baderlab.org` for your production apps or scripts.

## Example usage

Here, we provide usage examples in common languages for the main search API.  For more details, please refer to the Swagger documentation at [https://grounding.baderlab.org](https://grounding.baderlab.org), which is also accessible when running a local instance.

### Example search in JS

```js
const response = await fetch('http://hostname:port/search', {
  method: 'POST',
  headers: {
    'Content-Type': 'application/json'
  },
  body: JSON.stringify({ // search options here
    q: 'p53'
  })
});

const responseJSON = await response.json();
```

### Example search in Python

```python
import requests

url = 'http://hostname:port/search'
body = {'q': 'p53'}

response = requests.post(url, data = body)

responseJSON = response.json()
```

### Example in shell script via curl

```bash
curl -X POST "http://hostname:port/search" -H  "accept: application/json" -H  "Content-Type: application/json" -d "{  \"q\": \"p53\" }"
```

## Tool comparison

Here, we summarise a set of tools that overlap to some degree with the main use case of the Pathway Commons Grounding Search Service, where a user searches for a biological entity grounding by providing only a commonly-used synonym.  This table was last updated on 25 October 2021 (2021-10-25).

If you have developed a new tool in this space or your tool supports new features, let us know by making a pull request, and we'll add your revision to this table.

|                                                  | PC Grounding Search | GProfiler | GNormPlus (PubTator) | Gilda                                                 | BridgeDB  |
|--------------------------------------------------|---------------------|-----------|----------------------|-------------------------------------------------------|-----------|
| Allows for searching by synonym                  | ●                   |           | ●                    | ●                                                     |           |
| Supports multiple organisms                      | ●                   | ●         | ●                    | ●                                                     | ●         |
| Accepts organism ranking preference              | ●                   |           |                      |                                                       |           |
| Multiple organisms per query                     | ●                   |           | ●                    | Partial support (only one organism returned)          |           |
| Multiple results per query                       | ●                   |           |                      | One per type (e.g. protein)                           | ●         |
| Multiple results are ranked based on relevance   | ●                   |           |                      | ●                                                     |           |
| Speed/Throughput                                 | < 100 ms            | < 100 ms  | < 100ms              | < 100 ms                                              | < 1000 ms |
| Allows querying for a particular grounding by ID | ●                   | ●         | ●                    | ●                                                     | ●         |


## Grounding data

`grounding-search` uses data files provided by three public databases:

- [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene)
  - Information about genes
  - Alias: `ncbi`
  - Data file: [gene_info.gz](https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz)
- [ChEBI](https://www.ebi.ac.uk/chebi/) (`chebi`)
  - Information about small molecules of biological interest
  - Alias: `chebi`
  - Data file: [chebi.owl](https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl)
- [UniProt](https://www.uniprot.org/) (`uniprot`)
  - Information about proteins
  - Alias: `uniprot`
  - Data file: [uniprot_sprot.xml.gz](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz)

### Build index from source database files

If you have followed the Quick Start ("Run from source"), you can download and index the data provided by the  source databases `ncbi`, `chebi` and `uniprot` by running:

```
npm run update
```

### Restore index from Elasticsearch dump files

Downloading and building the index from source ensures that the latest information is indexed. Alternatively, to quickly retrieve and recreate the index a dump of a previously indexed Elasticsearch instance has been published on [Zenodo](https://zenodo.org/) under the following DOI:

[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4495013.svg)](https://doi.org/10.5281/zenodo.4495013)

This data is published under the [Creative Commons Zero v1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/legalcode) license.

To restore, create a running Elasticsearch instance and run:

```
npm run restore
```

To both restore and start the grounding-search server run:

```
npm run boot
```

NB: Index dump published on Zenodo is offered for demonstration purposes only. We do not guarantee that this data will be up-to-date or that releases of grounding-search software will be compatible with any previously published version of the dump data. To ensure you are using the latest data compatible with grounding-search, follow instructions in "Build the index database from source database files".

## Issues & feedback

To let us know about an issue in the software or to provide feedback, please [file an issue on GitHub](https://github.com/PathwayCommons/grounding-search/issues/new).

## Contributing

To make a contribution to this project, please start by please [filing an issue on GitHub that describes your proposal](https://github.com/PathwayCommons/grounding-search/issues/new).  Once your proposal is ready, you can make a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

## Configuration

The following environment variables can be used to configure the server:

- `NODE_ENV` : the environment mode, either `production` or `development` (default)
- `LOG_LEVEL` : the level for the log file (`info`, `warn`, `error`)
- `PORT` : the port on which the server runs (default 3000)
- `ELASTICSEARCH_HOST` : the `host:port` that points to elasticsearch
- `MAX_SEARCH_ES` : the maximum number of results to return from elasticsearch
- `MAX_SEARCH_WS` : the maximum number of results to return in json from the webservice
- `CHUNK_SIZE` : how many grounding entries make up a chunk that gets bulk inserted into elasticsearch
- `MAX_SIMULT_CHUNKS` : maximum number of chunks to insert simulteneously into elasticsearch
- `INPUT_PATH` : the path to the input folder where data files are located
- `INDEX` : the elasticsearch index name to store data from all data sources
- `UNIPROT_FILE_NAME` : name of the file where uniprot data will be read from
- `UNIPROT_URL` : url to download uniprot file from
- `CHEBI_FILE_NAME` : name of the file where chebi data will be read from
- `CHEBI_URL` : url to download chebi file from
- `NCBI_FILE_NAME` : name of the file where ncbi data will be read from
- `NCBI_URL` : url to download ncbi file from
- `NCBI_EUTILS_BASE_URL` : url for NCBI EUTILS
- `NCBI_EUTILS_API_KEY` : NCBI EUTILS API key
- `ESDUMP_LOCATION` : The location (URL, file path) of elasticdump files (note: terminate with '/')
- `ZENODO_API_URL`: base url for Zenodo
- `ZENODO_ACCESS_TOKEN`: access token for Zenodo REST API (Scope: `deposit:actions`, `deposit:write`)
- `ZENODO_BUCKET_ID`: id for Zenodo deposition 'bucket' (Files API)
- `ZENODO_DEPOSITION_ID`: id for Zenodo deposition (for a published dataset)

## Run targets

- `npm start` : start the server
- `npm stop` : stop the server
- `npm run watch` : watch mode (debug mode enabled, autoreload)
- `npm run refresh` : run clear, update, then start
- `npm test` : run tests for read only methods (e.g. search and get) assuming that data is already existing
- `npm test:sample` : run tests with sample data
- `npm run test:quality` : run the search quality tests (expects full db)
- `npm run test:quality:csv` : run the search quality tests and output a csv file
- `npm run lint` : lint the project
- `npm run benchmark` : run all benchmarking
- `npm run benchmark:source` : run benchmarking for `source` (i.e. `ncbi`, `chebi`)
- `npm run clear` : clear all data
- `npm run clear:source` : clear data for `source` (i.e. `ncbi`, `chebi`)
- `npm run update` : update all data (download then index)
- `npm run update:source` : update data for `source` (i.e. `ncbi`, `chebi`) in elasticsearch
- `npm run download` : download all data
- `npm run download:source` download data for `source` (i.e. `ncbi`, `chebi`)
- `npm run index` : index all data
- `npm run index:source` : index data for `source` (i.e. `ncbi`, `chebi`) in elasticsearch
- `npm run test:inputgen` : generate input test file for each `source` (i.e. `uniprot`, ...)
- `npm run test:inputgen` : generate input test file for `source` (i.e. `uniprot`, ...)
- `npm run dump` : dump the information for `INDEX` to `ESDUMP_LOCATION`
- `npm run restore` : restore the information for `INDEX` from `ESDUMP_LOCATION`
- `npm run boot` : run `clear`, `restore` then `start`; exit on errors

## Using Zenodo to store index dumps

[Zenodo](https://zenodo.org/) lets you you to store and retrieve digital artefacts related to a scientific project or publication. Here, we use Zenodo to store Elasticsearch index dump data used to quickly recreate the index used by grounding-search.

### Create and publish a new record deposition

Briefly, using their [RESTful web service API](https://developers.zenodo.org/), you can create a 'Deposition' for a record that has a 'bucket' referenced by a `ZENODO_BUCKET_ID` to which you can upload and download 'files' (i.e. `<ZENODO_API_URL>api/files/<ZENODO_BUCKET_ID>/<filename>`; list them with `https://zenodo.org/api/deposit/depositions/<deposition id>/files`). In particular, there are three files required to recreate an index, corresponding to the elasticsearch types: `data`; `mapping` and `analyzer`.

To setup follow these steps:

1. Get a `ZENODO_ACCESS_TOKEN` by creating a 'Personal access token' ([see docs for details](https://sandbox.zenodo.org/account/settings/applications/)). Be sure to add the `deposit:actions` and `deposit:write` scopes.
2. Create a recrod 'Deposition' by POSTing to `https://zenodo.org/api/deposit/depositions` with at least the following information, keeping in mind to set the header `Authorization = Bearer <ZENODO_ACCESS_TOKEN>`:
```json
{
	"metadata": {
		"title": "Elasticsearch data for biofactoid.org grounding-search service",
		"upload_type": "dataset",
		"description": "This deposition contains files with data describing an Elasticsearch index (https://github.com/PathwayCommons/grounding-search). The files were generated from the elasticdump npm package (https://www.npmjs.com/package/elasticdump). The data are the neccessary and sufficient information to populate an Elasticsearch index.",
		"creators": [
			{
				"name": "Biofactoid",
				"affiliation": "biofactoid.org"
			}
		],
		"access_right": "open",
		"license": "cc-zero"
	}
}
```
3. The POST response should have a 'bucket' (e.g. `"bucket": "https://zenodo.org/api/files/<uuid>"` ) within the `links` object. The variable `ZENODO_BUCKET_ID` is the value `<uuid>` in the example URL.
4. Publish. You'll want to dump the index and upload to Zenodo (`npm run dump`). You can publish this from the API by POSTing to `https://zenodo.org/api/deposit/depositions/<deposition id>/actions/publish`. Alternatively, log in to the Zenodo [web page](https://zenodo.org/deposit) and click 'Publish' to make the deposition public.

Once published, a deposition cannot be updated or altered. However, you can create a new version of a record (below).

### Create and publish a new version of a record

In this case, you already have a record which points to a published deposition (i.e. elasticsearch index files) and wish to create a new version for that record. Here, you'll create a new deposition under the same record:

1. Make a POST request to `https://zenodo.org/api/deposit/depositions/<deposition id>/actions/newversion` to create a new version. Alternatively, visit `https://zenodo.org/record/<deposition id>` where `deposition id` is that of the latest published version (default).
2. Fetch `https://zenodo.org/api/deposit/depositions?all_versions` to list all your depositions and identify the new deposition bucket id.
3. Proceed to upload (i.e. dump) your new files as described in "Create a new deposition", Step 3.

- Notes:
  - New version's files must differ from all previous versions
  - See https://help.zenodo.org/#versioning and https://developers.zenodo.org/#new-version for more info


## Testing

All files `/test` will be run by [Mocha](https://mochajs.org/).  You can `npm test` to run all tests, or you can run `npm test -- -g specific-test-name` to run specific tests.

[Chai](http://chaijs.com/) is included to make the tests easier to read and write.



## Publishing a release

1. Make sure the tests are passing: `npm test`
1. Make sure the linting is passing: `npm run lint`
1. Bump the version number with `npm version`, in accordance with [semver](http://semver.org/).  The `version` command in `npm` updates both `package.json` and git tags, but note that it uses a `v` prefix on the tags (e.g. `v1.2.3`).
  1. For a bug fix / patch release, run `npm version patch`.
  1. For a new feature release, run `npm version minor`.
  1. For a breaking API change, run `npm version major.`
  1. For a specific version number (e.g. 1.2.3), run `npm version 1.2.3`.
1. Push the release: `git push && git push --tags`
1. [Publish a GitHub release](https://github.com/PathwayCommons/grounding-search/releases/new) so that Zenodo creates a DOI for this version.


