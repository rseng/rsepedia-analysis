# SIGA.py

[![Build Status](https://travis-ci.org/candYgene/siga.svg?branch=master)](https://travis-ci.org/candYgene/siga)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1076437.svg)](https://doi.org/10.5281/zenodo.1076437)

*SIGA.py* is a command-line tool written in Python to generate *Semantically Interoperable Genome Annotations* from
text files in the Generic Feature Format ([GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)) according to the Resource Description Framework ([RDF](https://www.w3.org/TR/rdf11-concepts/)) specification.

<div align="center">
  <figure>
    <p>
      <img src ="doc/SIGA.png" alt="SIGA software architecture." />
      <figcaption>Fig. SIGA software architecture.</figcaption>
    </p>
  </figure>
</div>

## Key features ##
- _Input_:
  - one or more files in the GFF format (version 2 or 3)
  - `config.ini` file with ontology mappings and feature type amendments (if applicable)
- _Output_: genomic features stored in a [SQLite](https://sqlite.org/) database or serialized in one of the RDF formats:
  - [XML](https://www.w3.org/TR/rdf-syntax-grammar/)
  - [N-Triples](https://www.w3.org/TR/n-triples/)
  - [Turtle](https://www.w3.org/TeamSubmission/turtle/)
  - [Notation3](https://www.w3.org/DesignIssues/Notation3.html) (N3)
- check referential integrity for parent-child feature relationships in SQLite
- controlled vocabularies and ontologies used:
  - [DCMI](http://dublincore.org/documents/dcmi-terms/) terms (e.g. _creator_, _hasVersion_, _license_)
  - Sequence Ontology ([SO](http://www.sequenceontology.org/)) to describe feature types (e.g. _genome_, _chromosome_, _gene_, _transcript_) and their relationships (e.g. _has part_/_part of_, _genome of_, _transcribed to_, _translated_to_)
  - Feature Annotation Location Description Ontology ([FALDO](https://github.com/JervenBolleman/FALDO))

## Requirements ##

    python (>=2.7)
    docopt (0.6.2
    RDFLib (4.2.2)
    gffutils (https://github.com/arnikz/gffutils)
    optional: RDF store to query ingested data using SPARQL (e.g. using Virtuoso or Berkeley DB)


## Installation ##

```
git clone https://github.com/candYgene/siga.git
cd siga
virtualenv .sigaenv
source .sigaenv/bin/activate
pip install -r requirements.txt
```

## How to use ##

**Command-line interface**

```
Usage:
  SIGA.py -h|--help
  SIGA.py -v|--version
  SIGA.py db [-ruV] [-d DB_FILE | -e DB_FILEXT] GFF_FILE...
  SIGA.py rdf [-V] [-o FORMAT] [-c CFG_FILE] DB_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF version 2 or 3.
  DB_FILE...       Input database file(s) in SQLite.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Show verbose output in debug mode.
  -c FILE          Set the path of config file [default: config.ini]
  -d DB_FILE       Create a database from GFF file(s).
  -e DB_FILEXT     Set the database file extension [default: .db].
  -r               Check the referential integrity of the database(s).
  -u               Generate unique IDs for duplicated features.
  -o FORMAT        Output RDF graph in one of the following formats:
                     turtle (.ttl) [default: turtle]
                     nt (.nt),
                     n3 (.n3),
                     xml (.rdf)
 ```

**Input files**

Small test set in `examples/features.gff3` including `config.ini`. Alternatively, download tomato or potato genome annotations.

```
wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/ITAG2.4_gene_models.gff3
wget http://solanaceae.plantbiology.msu.edu/data/PGSC_DM_V403_genes.gff.zip
```

**Generate RDF graph**

1. GFF->DB

    ```
    python SIGA.py db -rV ../examples/features.gff3 # output *.db
    ```

2. DB->RDF (default: `turtle`)

    ```
    python SIGA.py rdf -c config.ini ../examples/features.db # output *.ttl
    ```

Summary of I/O files:

- config file: `config.ini`
- GFF file: `features.gff3`
- SQLite DB file: `features.db`
- RDF Turtle file: `features.ttl`

**Import RDF graph into Virtuoso RDF Quad Store**

See the [documentation](http://virtuoso.openlinksw.com/dataspace/doc/dav/wiki/Main/VirtBulkRDFLoader) on bulk data loading.

Edit `virtuoso.ini` config file by adding _/mydir/_ to _DirsAllowed_.

Connect to db server as `dba` user:

`isql 1111 dba dba`

Delete (existing) RDF graph if necessary:

`SPARQL CLEAR GRAPH <http://solgenomics.net/genome/Solanum_lycopersicum> ;`

Delete any previously registered data files:

`DELETE FROM DB.DBA.load_list ;`

Register data file(s):

`ld_dir('/mydir/', 'features.ttl', 'http://solgenomics.net/genome/Solanum_lycopersicum') ;`

List registered data file(s):

`SELECT * FROM DB.DBA.load_list ;`

Bulk data loading:

`rdf_loader_run() ;`

Re-index triples for full-text search (via Faceted Browser):

`DB.DBA.VT_INC_INDEX_DB_DBA_RDF_OBJ() ;`

Note: A single data file can be uploaded using the following command:

`SPARQL LOAD "file:///mydir/features.ttl" INTO "http://solgenomics.net/genome/Solanum_lycopersicum" ;`

Count imported RDF triples:

```
SPARQL
SELECT COUNT(*)
FROM <http://solgenomics.net/genome/Solanum_lycopersicum>
WHERE { ?s ?p ?o } ;
```


**Alternatively, import RDF graph into Berkeley DB (requires [Redland](http://librdf.org/) RDF processor)**

```
rdfproc features parse features.ttl turtle
rdfproc features serialize turtle
```

## License ##
The software is released under Apache License, Version 2.0.
