[![DOI](https://zenodo.org/badge/138730107.svg)](https://zenodo.org/badge/latestdoi/138730107)
![Website](https://img.shields.io/website-up-down-green-red/http/lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080.svg?label=fdp)
![Website](https://img.shields.io/website-up-down-green-red/http/lofar-ld.lofarpipelines-astron.surf-hosted.nl:8890.svg?label=virtuoso)

# LOFAR Linked Data Platform

This repository contains services to improve the FAIRness of the [LOFAR](http://www.lofar.org/) [Long Term Archive](https://lta.lofar.eu/) (meta)data using semantic web technologies:

* [Virtuoso Universal Server](https://virtuoso.openlinksw.com/) (including some code from [pbg-ld](https://www.research-software.nl/software/linked-data-platform-for-plant-breeding-genomics))
* [FAIR Data Point](https://www.research-software.nl/software/fairdatapoint) (FDP)

**1. Clone this git repo.**

`git clone https://github.com/EOSC-LOFAR/lofar-ld.git`

**2. Install [Ansible](https://www.ansible.com/).**

`pip install ansible`

**3. Deploy services.**

`ansible-playbook -i inventory playbook.yml` # on a remote host specified in the 'inventory' file

`make -C src` # or locally

**4. [Login](http://localhost:8890/conductor) to running Virtuoso instance for admin tasks.**

Use `dba` for both account name and password.

**5. Run queries via Virtuoso [SPARQL endpoint](http://localhost:8890/sparql) or browse data via [Faceted Browser](http://localhost:8890/fct/) (no login required).**

RDF graphs:IRIs (_A-Box_)
 * (subset of) LOFAR DBView: `https://lta-dbview.lofar.eu/`
 * FDP metadata:
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/fdp`
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/catalog/astron-01`
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/dataset/lofar-lta-dbview`
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/distribution/lofar-lta-dbview-sparql`
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/distribution/lofar-lta-dbview-rdfdump` 	
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/distribution/lofar-lta-dbview-sqldump` 	
   * `http://lofar-ld.lofarpipelines-astron.surf-hosted.nl:8080/distribution/lofar-lta-dbview-csvdump`
  
RDF graphs:IRIs (_T-Box_)
 * Astronomy: `http://www.astro.umd.edu/~eshaya/astro-onto/owl/astronomy.owl#`
 * Physics: `http://www.astro.umd.edu/~eshaya/astro-onto/owl/physics.owl#`
 * Instruments: `http://www.astro.umd.edu/~eshaya/astro-onto/owl/instruments.owl#`
 * IVOAT: `http://www.ivoa.net/rdf/Vocabularies/IVOAT`
 * IVOAO: `http://www.astro.umd.edu/~eshaya/astro-onto/owl/IVOAO.owl#`
 * UCD: `http://www.ivoa.net/rdf/Vocabularies/UCD`
 * IAUT93: `http://www.ivoa.net/rdf/Vocabularies/vocabularies-20091007/IAUT93/IAUT93.rdf#`
 * SEM: `https://semanticweb.cs.vu.nl/2009/11/sem/sem.rdf#`
 
See [wiki](../../wiki) for further details.
