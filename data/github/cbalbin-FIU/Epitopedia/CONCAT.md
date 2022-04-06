# Epitopedia



[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2021.08.26.457577%20-blue)](https://doi.org/10.1101/2021.08.26.457577)

![Epitopedia Screencap Example](example_output/epitopedia_screencap.gif)


## Getting started

The quickest way to start using Epitopedia is by downloading the docker container which contains all the dependencies preinstalled:

```shell
git clone https://github.com/cbalbin-bio/Epitopedia.git

docker pull cbalbin/epitopedia
```

Epitopedia requires the [PDB in mmCIF](https://www.wwpdb.org/ftp/pdb-ftp-sites) format, EpitopediaDB and EPI-SEQ DB. EpitopediaDB and EPI-SEQ DB can be downloaded [here](https://fiudit-my.sharepoint.com/:u:/g/personal/cbalbin_fiu_edu/EZWcy2ki66dGrUGKPhwXYwUBUcENzab55CA0lZ4_laOuBQ?e=zJqGrC).


To download the entirety of PDB in mmCIF format:
```shell
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
```

**OR**
<br>

To download the only the PDB files present in EpitopediaDB (EPI-PDB) you can supply the pdb_id_list.txt to rsync:
```shell
rsync -rlpt -v -z --delete --port=33444 --include-from=/path/to/pdb_id_list.txt \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
```

To run Epitopedia provide the paths to the various directories discussed below.

The data directory should contain Epitopedia DB (epitopedia.sqlite3) and EPI-SEQ (EPI-SEQ.fasta*) which can be downloaded [here](https://fiudit-my.sharepoint.com/:u:/g/personal/cbalbin_fiu_edu/EZWcy2ki66dGrUGKPhwXYwUBUcENzab55CA0lZ4_laOuBQ?e=zJqGrC).

The mmcif directory should point to the sharded PDB directory in mmCIF format as downloaded above.

NOTE: you may need to unzip the mmCIF directory:
```shell
gunzip -r mmCIF
```

The output directory is where the output files will be written.

Replace the the paths on the left side of the colon with the actual **absolute** path on your local system. The paths on the right side of the colon are internal and should not be altered.

```shell
python3 Epitopedia/docker/run_epitopedia.py \
/Path/to/Output/Dir/ \
/Path/to/PDB/Dir/ \
/Path/to/Data/Dir/ \
--afdb-dir /Path/to/AlphaFold/Dir/ \
--taxid-filter 11118 --PDB-IDS 6VXX_A
```

NOTE: on some systems you may need to run docker with sudo.

It is recommended to use the flag taxid_filter to prevent the input protein from finding itself or other versions of itself. For example, if we wnted to find mimics of the SARS-CoV-2 spike protien (6VXX) is a SARS-CoV-2 protein
we could use a taxid_filter of 11118 to prevent finding mimics in other Coronaviridae. The [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) will be helpful in determining what taxid to specify.




Epitopedia can run on multiple input structures to represent a conformational ensemble. To do so, simply provide a list of structures in the format PDBID_CHAINID as shown below.
```shell
run_epitopedia.py --PDB-IDS 6VXX_A 6VXX_B 6XR8_A 6XR8_B
```

Epitopedia defaults to a span length of 5, surface accesbility cutoff of 20% surface accesbility span legnth of 3, and no taxa filter, but these parameters can be set using the follow flags:

Flag | Description
------------ | -------------
--span | Minimum span length for a hit to progress
--rasa | Cutoff for relative accessible surface area
--rasa_span | Minimum consecutive accesssible residues to consider a hit a SeqBMM
--taxid_filter | taxa filter; example to filter out all Coronaviridae --taxid_filter 11118
--rmsd | Max RMSD to still be considered a structural mimic
--view | View results from a previous run
--port | Port to be used by webserver
--use-afdb | Include AFDB in search
--pplddt | Minimum protein pLDDT score a structure predicted by alphafold must have to be considered
--mplddt | Minimum average local pLDDT score a region predicted by alphafold must have to be considered
## Output

Example output files 6VXX_A with a taxid_filt of 11118 as an input can be found [here](example_output).
<br>
Definitions for the output file headers can be found [here](headers.md).
### Intermediate Output

Epitopedia will output the following files at various stages of execution:

File Name | Description
------------ | -------------
EPI_SEQ_hits_{pdb_id(s)}.tsv  | Contains the raw results from the BLAST search of the input structure against EPI-SEQ
EPI_SEQ_span_filt_hits_{pdb_id(s)}.tsv | Contains hits  with consecutive spans that meet the set minimum span length
EPI_SEQ_span_filt_acc_hits_{pdb_id(s)}.tsv | Contains the above spans that contain the minimum span of accessible residues
EPI_PDB_hits_{pdb_id(s)}.tsv" | Contains epitope source sequences against EPI_PDB hits
EPI_PDB_fragment_pairs_{pdb_id(s)}.tsv | Contains structurally aligned fragment pairs consisting of spans of the input structure aligned against the structural representatives
EPI_PDB_fragment_pairs_{pdb_id(s)}_ranked.tsv | Contains the above but ranked from best to worst RMSD

### Final Output


Epitopedia will show the best hit per epitope motif if there are redundant source sequences at the final stage of the execution. There results can be viewed in a tsv file ([**Example**](example_output/EPI_PDB_fragment_pairs_6VXX_A_best.tsv)) or a more legible HTML file ([**Example**](https://cbalbin-fiu.github.io/)).

<!-- [**Click here**](https://cbalbin-fiu.github.io/) for an example of the HTML output using input structure 6VXX_A with a taxid_filt of 11118. -->



## Epitopedia database generation

Epitopedia uses IEDB and PDB to generate EpitopediaDB, which is used in the molecular mimicry search.

Generation of the database takes some time (~12 hours). Thus, the EpitopediaDB is provided above.

To create the EpitopediaDB, download [IEDB](https://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz) and a mmCIF version of PDB.

Point the container to the approriate paths for the IEDB, PDB (mmCIF format) and a data directory where the databases will be written.


```shell
docker run --rm -it \
-v /Path/To/iedb_public.sql:/app/iedb \
-v /Path/to/mmCIF/Dir/:/app/mmcif \
-v /Path/to/Data/Dir/:/app/data \
cbalbin/epitopedia generate_database.py
```

## License

This software is released under the [MIT License](LICENSE).

Software and databases used in Epitopedia may be released under various licenses:

Software:

* [DSSP](https://github.com/cmbi/dssp/blob/master/COPYING)
* [TM-align](https://zhanggroup.org/TM-align/TMalign.f)
* [NCBI BLAST](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE)
* [Entrez](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE)
* [MMseqs2](https://github.com/soedinglab/MMseqs2/blob/master/LICENSE.md)
* [mysql2sqlite](https://github.com/dumblob/mysql2sqlite/blob/master/LICENSE)
* [sqlite](https://www.sqlite.org/copyright.html)
* [Docker](https://github.com/docker/cli/blob/master/LICENSE)
* [Flask](https://github.com/pallets/flask/blob/main/LICENSE.rst)
* [gemmi](https://github.com/project-gemmi/gemmi/blob/master/LICENSE.txt)
* [rich](https://github.com/willmcgugan/rich/blob/master/LICENSE)
* [biopython](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
* [dataclasses-json](https://github.com/lidatong/dataclasses-json/blob/master/LICENSE)
* [python](https://github.com/python/cpython/blob/main/LICENSE)

Databases:
* [IEDB](http://www.iedb.org/terms_of_use_v3.php)
* [PDB](https://www.rcsb.org/pages/usage-policy)


## Reference

If you use Epitopedia in your work, please cite:
```
Epitopedia: identifying molecular mimicry of known immune epitopes
Christian Andrew Balbin, Janelle Nunez-Castilla, Jessica Siltberg-Liberles
bioRxiv 2021.08.26.457577; doi: https://doi.org/10.1101/2021.08.26.457577
```## Description of Headers

This document contains the description for the headers used across various output files

### EPI-SEQ Headers

Header | Description
------------ | -------------
EPI_SEQ Input Structure | The structure provided to Epitopedia that this row relates to. Its sequence is extracted and used as a query against EPI_SEQ
EPI_SEQ Epitope ID | The Epitope ID for the subject in the resulting BLAST alignment
EPI_SEQ Input Structure Seq Start Pos | The starting position of the query sequence, extracted from the input structure, in the resulting BLAST Alignment (NOTE: due to inconsistencies in the PDB/mmCIF format, this is not necessarily the starting position in the Structure itself)
EPI_SEQ Input Structure Seq End Pos | Ending position for what is described above
EPI_SEQ Epitope Start Pos | The starting position for the subject (Epitope) in the resulting BLAST alignment
EPI_SEQ Epitope Stop Pos | Ending position for what is described above
EPI_SEQ Aln Input Struc Seq | The query portion of the aligned sequence resulting from the BLAST alignment
EPI_SEQ Aln Epitope Seq | The subject portion of the aligned sequence resulting from the BLAST alignment
EPI_SEQ Evalue | The Evalue for the resulting BLAST alignment
EPI_SEQ Qcov | The query coverage for the resulting BLAST alignment
EPI_SEQ Pident | The percent identity for the resulting BLAST alignment
EPI_SEQ Epitope Taxid | The taxid for the organism of which the epitope is sourced from
EPI_SEQ Aln Cigar | the cigar sequence for the BLAST alignment
EPI_SEQ Span Ranges | List of lists containing the starting and ending position for spans of consecutive matches in the BLAST alignment
EPI_SEQ Span Lengths | List of the lengths of the EPI_SEQ Span Ranges
EPI_SEQ Span Seqs | List of amino acid sequences resulting from these spans
PDB_DSSP Input Struc ASA | List of the ASA values as calculaed from DSSP for the aligned portion of the input structure
mmCIF_SEQ Input Struc Solv Seq | String where "?" describes unsolved residues in the aligned input structure (query) sequence
mmCIF_SEQ Input Struc Res Nums | List mapping the aligned input structure (query) sequence to the associated residue number, if present, in the PDB structure


### EPI-PDB Headers

Header | Description
------------ | -------------
SeqBMM Motif | motif of the SeqBMM found in the EPI_SEQ search step
SeqBMM Input Struc Res Nums | List containing the residue numbers in the input structure that correspond to the reidues in the SeqBMM motif
SeqBMM Acc | List contianing the surface accessibility classification for each residue in the SeqBMM motif; A = accessible, B = Buried
EPI_PDB Epi Source Acc | The epitope soruce sequence accession which is used as a query against EPI_PDB
EPI_PDB Rep PDB | The PDB_ID for the resulting representative strucure (target) from the query against EPI_PDB
EPI_PDB Qcov | The query coverage for the mmSEQ result
EPI_PDB Pident | The percent identity for the mmSEQ result
EPI_PDB Evalue | The Evalue for the mmSEQ result
mmCIF_SEQ Rep Res | The residue sequence of the structural representative
mmCIF_SEQ Rep Solv | String where "?" describes unsolved residues in the structural representative
mmCIF_SEQ Rep Num | List mapping the sequence of the structural representative to the associated residue number, if present, in the PDB structure
EPI_PDB Rep Res Nums | List containing the residue numbers in the structural representative that correspond to the reidues in the EPI_PDB motif
EPI_PDB Input Dice Path | Path to the SeqBMM motif dice from the input structure 
EPI_PDB Rep Dice Path | Path to the SeqBMM motif dice from the representative srructure
EPI_PDB TMalign RMSD | resulting RMSD value for the alignment of the dices described above
EPI_PDB TMalign TMscore | TMscore for above
EPI_PDB TMalign PDB | Path to the PDB file containing the structural superposition of the two dices described above
EPI_PDB Rep Acc | List contianing the surface accessibility classification for each residue in the representative motif; A = accessible, B = Buried
EPI_PDB Input Motif Perc Acc | Percentage of accessible residue in the motif span of the input structure
EPI_PDB Rep Motif Perc Acc | Percentage of accessible residue in the motif span of the representative structure
EPI_PDB Perc Acc Agree | Percentage of surface accessibility agreement between the motif span of the input structure and the motif span of the representative structure
IEDB_FILT Epitope ID | Epitope ID of the epitope which contains the SeqBMM
IEDB_FILT Epitope Seq | Full sequence of the Epitope
IEDB_FILT Source Seq Acc | Acession for the source sequence of which the epitope is a subsequence of
IEDB_FILT Start Pos | Starting position of the epitope in the epitope source sequence
IEDB_FILT Stop Pos | Stoping position for above
IEDB_FILT Source DB | Database from which the source sequence was obtained
IEDB_FILT Source Title | Title for the source sequence protein
IEDB_FILT Source Taxid | Taxid for the organism associated with the source sequence
IEDB_FILT Source Org | Organism name the source sequence is associated with
IEDB_FILT Source Seq | Epitope Source sequence
IEDB_FILT Iacc | Internal acession number

[For remaining headers in this file see EPI-SEQ headers](#EPI-SEQ-Headers)


