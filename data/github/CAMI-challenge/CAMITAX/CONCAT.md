# CAMITAX: Taxon labels for microbial genomes

![The CAMITAX taxonomic assignment workflow](workflow.png "The CAMITAX taxonomic assignment workflow")

**The CAMITAX taxonomic assignment workflow.**
CAMITAX assigns one NCBI Taxonomy ID (taxID) to an input genome *G* by combining genome distance-, 16S rRNA gene-, and gene homology-based taxonomic assignments with phylogenetic placement.
**(A) Genome distance-based assignment.**
CAMITAX uses Mash to estimate the average nucleotide identity (ANI) between *G* and more than a hundred thousand microbial genomes in RefSeq, and assigns the lowest common ancestor (LCA) of genomes showing >95% ANI, which was found to be a clear species boundary.
**(B) 16S rRNA gene-based assignment.**
CAMITAX uses Dada2 to label *G*'s 16S rRNA gene sequences using the naïve Bayesian classifier method to assign taxonomy across multiple ranks (down to genus level), and exact sequence matching for species-level assignments, against the SILVA or RDP database.
**(C) Gene homology-based assignments.**
CAMITAX uses Centrifuge and Kaiju to perform gene homology searches against nucleotide and amino acid sequences in NCBI's nr and nt (or proGenomes' genes and proteins datasets), respectively. CAMITAX determines the interval-union LCA (iuLCA) of gene-level assignments and places *G* on the lowest taxonomic node with at least 50% coverage.
**(D) Phylogenetic placement.**
CAMITAX uses Pplacer to place *G* onto a fixed reference tree, as implemented in CheckM, and estimates genome completeness and contamination using lineage-specific marker genes.
**(E) Classification algorithm.**
CAMITAX considers the lowest consistent assignment as the longest unambiguous root-to-node path in the taxonomic tree spanned by the five taxIDs derived in (A)–(D), i.e. it retains the most specific, yet consistent taxonomic label among all tools.

## Requirements

All you need is [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) (or [Singularity](https://singularity.lbl.gov/)). This is the recommended way to run CAMITAX and (by default) CAMITAX requires 8 CPU cores and 24 GB of memory.

**Plan B:** You may run CAMITAX without software containers. However, this is not recommended and you have to install [all software dependencies](requirements.txt) by yourself.

*If you need any help or further guidance: Please [get in touch](https://github.com/CAMI-challenge/CAMITAX/issues)!*

## User Guide

### Installation

CAMITAX relies on multiple reference databases (which we do not bundle by default, due to their sheer size). You can either [build them from scratch](https://github.com/CAMI-challenge/CAMITAX/blob/master/db/README.md) or simply use the latest of our "official" [releases](https://doi.org/10.5281/zenodo.1250043). To do so, please run:
```
nextflow pull CAMI-challenge/CAMITAX
nextflow run CAMI-challenge/CAMITAX/init.nf --db /path/to/db/folder
```
**Warning:** This will download ~30 GB of data, expect this to run a while! `/path/to/db/folder` should have >100 GB of available disk space. Note that you have to do this only once; specify the location in all future CAMITAX runs.

**Warning:** To foster reproducibility, we strongly recommend that you use our "official" releases and we will continue to provide stable and versioned updates in the future.

### Input

CAMITAX expects all input genomes in (genomic/nucleotide multi-)FASTA format.
If your input genomes are in the folder `input/` with file extension `.fasta`, please run:
```
nextflow run CAMI-challenge/CAMITAX -profile docker --db /path/to/db/folder --i input --x fasta
```
If you want to use Singularity instead of Docker (without sudo), please replace `docker` with `singularity` as profile.

### Output

CAMITAX outputs a tab-seperated file `camitax.tsv` containing the individual taxon assignments in the `data` folder.

**Warning:** While CAMITAX is built around computational reproducibility, results might sometimes be slightly different from run to run (but of comparable quality) because software used within (e.g. CheckM) are non-deterministic.

*Again, if you need any help or further guidance: Please [get in touch](https://github.com/CAMI-challenge/CAMITAX/issues)!*

## Citation

* Bremges, Fritz & McHardy (2020). **CAMITAX: Taxon labels for microbial genomes.** *GigaScience*, 9, 1:1–7. doi:[10.1093/gigascience/giz154](https://doi.org/10.1093/gigascience/giz154)
* Sczyrba, Hofmann, Belmann, *et al.* (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
##  Mash sketch database

```
# Follow the Genomes Download FAQ (minus the 'complete' part): https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete
wget -O assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F "\t" '$11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
xargs -n 1 -P 32 wget -q < ftpfilepaths

# Double-check that all genomes were downloaded, please manually download missing file(s) – if any:
diff <(sed 's/.*\///' < ftpfilepaths | sort) <(ls -1 | grep "genomic.fna.gz" | sort)

# Prefix the filename with the NCBI Taxonomy ID (field 6 in assembly_summary.txt):
mkdir mash_genomes && cut -f6,20 assembly_summary.txt | sed 's/ftp.*\///' | while read i; do cp $(echo $i | cut -d' ' -f2)\_genomic.fna.gz mash_genomes/$(echo $i | tr ' ' '_'); done

# Create the Mash index
cd mash_genomes && mash sketch -p 32 -o ../RefSeq * && cd ..
mash info -t RefSeq.msh | tail -n +2 | cut -f3 | cut -f1 -d'_' | sort | uniq > RefSeq.ids
```

## DADA2 reference databases

```
# Download RDP trainset 16: https://doi.org/10.5281/zenodo.801827
wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz
wget https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz

# Download SILVA version 132: https://doi.org/10.5281/zenodo.1172782
wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz
wget https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz

# Source: https://benjjneb.github.io/dada2/training.html
```

## Kaiju index (and NCBI taxonomy)

```
# Kaiju's makeDB.sh downloads the proGenomes and NCBI Taxonomy databases, and builds the index:
makeDB.sh -p
```

## Centrifuge index

```
# Download the proGenomes database (genes):
wget http://progenomes.embl.de/data/repGenomes/representatives.genes.fasta.gz

# Replace taxon IDs found in merged.dmp by their updated IDs:
gunzip -c representatives.genes.fasta.gz | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/>(\d+)\.(\S+)/){print ">",defined($h{$1})?$h{$1}:$1,".",$2;}else{print}' -- -m=merged.dmp > centrifuge_db.fna

# Generate mapping file (tab-separated file mapping sequence IDs to taxon IDs):
perl -lne 'if(/>(\d+)\.(\S+)/){print $1,".",$2,"\t",$1}' < centrifuge_db.fna > centrifuge_db.conv

# Build the index (creates proGenomes.[1-4].cf):
centrifuge-build --conversion-table centrifuge_db.conv --taxonomy-tree nodes.dmp --name-table names.dmp centrifuge_db.fna proGenomes
```

## CheckM database

```
# Download and unpack CheckM database
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xfv checkm_data_2015_01_16.tar.gz

# Source: https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm
```
