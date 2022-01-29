# Linked Data Platform for Animal Breeding & Genomics

[![Build Status](https://travis-ci.org/candYgene/abg-ld.svg?branch=master)](https://travis-ci.org/candYgene/abg-ld)

This software provides semantically integrated genotypic/phenotypic data on animals to enable ranking of candidate genes associated with traits of interest (e.g. nipple quantity in pig).

**1. Clone this git repo.**

`git clone https://github.com/candYgene/abg-ld.git`

**2. Start a [Docker container](https://hub.docker.com/r/candygene/docker-virtuoso/) with [Virtuoso Universal Server](http://virtuoso.openlinksw.com/) & ingest data in [RDF](https://www.w3.org/RDF/).**

```
cd abg-ld/src
make all # with defaults: CONTAINER_NAME=virtuoso and CONTAINER_PORT=8890 (in development)
make -e all CONTAINER_NAME=abg-ld CONTAINER_PORT=80 # override defaults (in production)
```

Note: other `make` rules: `pull-image`, `build-image`, `start-srv`, `stop-srv`, `restart-srv`, `install-pkgs`, `get-rdf`, `import-rdf`, `update-rdf`, `post-install` and `clean`.

**3. [Login](http://localhost:8890/conductor) to running Virtuoso instance for admin tasks.**

Use `dba` for both account name and password.

**4. Run queries via Virtuoso [SPARQL endpoint](http://localhost:8890/sparql) or browse data via [Faceted Browser](http://localhost:8890/fct/) (no login required).**

RDF graphs:IRIs (_A-Box_)
  * Pig QTLdb: `http://www.animalgenome.org/QTLdb/pig`
  * Ensembl: `http://www.ensembl.org/pig`, `http://www.ensembl.org/human`
  * UniProt: `http://www.uniprot.org/proteomes/pig`
  * OMIM: `http://bio2rdf.org/omim_resource:bio2rdf.dataset.omim.R4`

RDF graphs:IRIs (_T-Box_)
  * FALDO: `http://biohackathon.org/resource/faldo.rdf`
  * SO[FA]: `http://purl.obolibrary.org/obo/so.owl`
  * SIO: `http://semanticscience.org/ontology/sio.owl`
  * RO: `http://purl.obolibrary.org/obo/ro.owl`
  * VT: `http://purl.obolibrary.org/obo/vt.owl`
  * CMO: `http://purl.obolibrary.org/obo/cmo.owl`
  * LPT: `http://purl.bioontology.org/ontology/LPT`
  * LBO: `http://purl.bioontology.org/ontology/LBO`
  * Uniprot Core: `http://purl.uniprot.org/core/`

Further details can be found on the [wiki](https://github.com/candYgene/abg-ld/wiki/Home).
