[![DOI](https://zenodo.org/badge/327435617.svg)](https://zenodo.org/badge/latestdoi/327435617)

# Meta-Fish-Lib: A generalised, dynamic DNA reference library for metabarcoding of fishes

Collins, R.A., Trauzzi, G., Maltby, K.M., Gibson, T.I., Ratcliffe, F.C., Hallam, J., Rainbird, S., Maclaine, J., Henderson, P.A., Sims, D.W., Mariani, S. & Genner, M.J. (2021). Meta-Fish-Lib: A generalised, dynamic DNA reference library pipeline for metabarcoding of fishes. _Journal of Fish Biology_. [https://doi.org/10.1111/jfb.14852](https://doi.org/10.1111/jfb.14852).

![SeaDNA Logo](assets/logo.svg)

This repository hosts a comprehensive multi-locus mitochondrial DNA reference library dataset for fish species of the United Kingdom (UK), derived from the [NCBI GenBank](https://www.ncbi.nlm.nih.gov/nucleotide) and [Barcode of Life BOLD](http://www.boldsystems.org/index.php) databases. The dataset includes both freshwater and marine species, and can be employed in a variety of applications. e.g. DNA barcoding of human food products using full COI barcodes, to metabarcoding of gut or environmental samples using fragments of 12S. The library will be updated with each new GenBank release. Both common and rare UK fish species are included. A species coverage report for all primer sets can be found at [assets/reports-tables.md](assets/reports-tables.md). This UK reference library is curated and ready-to-use, but the code provided here can easily generate a new reference library for a different location (see [FAQ](#FAQ)).

In addition to providing quality controlled and curated fish references, this reference library has several unique features that make it useful to the wider DNA barcoding and DNA metabarcoding communities:

* Flexible - the library is not limited to any particular metabarcode locus or primer set. I have included the most popular ones (Table 1), but new ones can be added as required.
* Comprehensive - seaching by single gene names can often miss critical results due to poorly annotated records, but using hidden Markov models it is simple to extract homologous DNA fragments from large dumps of sequence data.
* Exhaustive - searching by species names can exclude potential hits because of changes in taxonomy, but here we search for all species synonyms, and then subsequently validate those names to provide a taxonomically up-to-date reference library. 
* Reliable - sequence data on GenBank are frequently misannotated with incorrect species names, but we have created a list of dubious quality accessions that are automatically excluded when the reference library is loaded each time. We use phylogenetic quality control methods to assist in screening each new GenBank version and update this list accordingly.
* Dynamic - it's easy to update to each new GenBank release (see code [below](#bash-terminal)), and the versioning of this repository reflects the GenBank release on which it was made.
* Quick - the final reference library can be downloaded onto your computer in just a few seconds with only two packages loaded and seven lines of R code ([below](#retrieve-latest-reference-library)). Generating this reference library from scratch ([below](#bash-terminal)) takes a couple of hours, with the phylogenetic quality control steps completing overnight.
* Customisable - by forking or cloning the repository, custom modifications can be made, e.g. excluding particular species, making taxonomic changes, or using a completely different list of species.
* Self contained - to recreate the reference libraries, all code and R package versions are found within in this self contained project, courtesy of [renv](https://rstudio.github.io/renv/articles/renv.html). This means less risk of clashing installations or broken code when packages and R versions upgrade.
* Citable - DOIs are issued with each new GenBank release.

This README outlines the contents of the repository and a brief description of the workflow involved in creating/updating a metabarcoding reference library, as well instructions to simply access the current data immediately. If an error is apparent, raise a ticket in [Issues](https://github.com/genner-lab/meta-fish-lib/issues) or submit a pull request.

The work is part of the NERC funded [SeaDNA Project](https://twitter.com/SeaDNAproject), and should be cited using version appropriate DOIs that are found in [Releases](https://github.com/genner-lab/meta-fish-lib/releases), or the Collins et al. (2021) _Journal of Fish Biology_ article ([https://doi.org/10.1111/jfb.14852](https://doi.org/10.1111/jfb.14852)) describing the software.

### TL;DR (give me the data)

If you require simply the final reference library file for immediate use, it can be downloaded directly using the R code below, and converted into FASTA and CSV formats for any of the available primer sets in Table 1.

##### Retrieve latest reference library:

```r
### START A FRESH R SESSION ###

# load packages (install if required)
library("tidyverse")
library("ape")

# load remote references and scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")

# choose a metabarcode fragment (primer set) from the following:
print(c("coi.lerayxt","coi.ward","12s.miya","12s.riaz","12s.valentini","12s.taberlet","16s.berry","cytb.minamoto"))
# change 'frag' argument as appropriate:
reflib.sub <- subset_references(df=reflib.orig, frag="12s.miya")

# convert to fasta file
# uses the standard database id field ('dbid') as a label
# 'dbid' is the GenBank GI number, or the BOLD processID number
# custom labels can be created with 'mutate()' and 'paste()' using other fields and changing 'namecol' argument - see FAQ
reflib.fas <- tab2fas(df=reflib.sub, seqcol="nucleotides", namecol="dbid")

# write out fasta file and corresponding csv table
ape::write.FASTA(reflib.fas, file="references.fasta")

# write out corresponding csv table
# this requires readr v1.4 or higher; change to 'path="references.csv"' for older versions, or use write.csv() instead
readr::write_csv(reflib.sub, file="references.csv")
```

Particular attention should be paid to cleaning steps in `scripts/references-clean.R`; sequences flagged as unreliable (using phylogenetic quality control) are listed in `assets/exclusions.csv` and excluded, while sequences flagged by NCBI as "unverified" are also removed. Taxonomic changes are also made, automatically via validating names against FishBase, and also custom changes (*Cottus perifretum* relabelled as *Cottus cottus*, *Atherina presbyter* relabelled as *Atherina boyeri*, and *Pungitius laevis* relabelled as *Pungitius pungitius*. Where are changes are made, both the original GenBank names and the validated FishBase names are provided (see Table 2).

**Table 1: Available primer sets**

Study | Official name | Nickname | Locus
----- | ----- | ----- | -----
[Miya et al. (2015)](http://dx.doi.org/10.1098/rsos.150088) | MiFish U/E | 12s.miya | 12S
[Taberlet et al. (2018)](http://dx.doi.org/10.1093/oso/9780198767220.001.0001) | Tele02 | 12s.taberlet | 12S
[Valentini et al. (2016)](http://dx.doi.org/10.1111/mec.13428) | L1848/H1913 | 12s.valentini | 12S
[Riaz et al. (2011)](http://dx.doi.org/10.1093/nar/gkr732) | 12S-V5 | 12s.riaz | 12S
[Wangensteen et al. (2018)](http://dx.doi.org/10.7717/peerj.4705) | Leray-XT | coi.leray | COI
[Ward et al. (2005)](http://dx.doi.org/10.1098/rstb.2005.1716) | FishF1/R1 | coi.ward | COI
[Berry et al. (2017)](http://dx.doi.org/10.1002/ece3.3123) | Fish16sF/D | 16s.berry | 16S
[Minamoto et al. (2012)](http://dx.doi.org/10.1007/s10201-011-0362-4) | L14912-CYB | cytb.minamoto | cytb

### Create/update the reference library manually

You don't need to run this code below if you just want a copy of the reference library (run code above). This code below is if you want to update it yourself or want to modify and make a new library. I will endeavour to keep this repository up-to-date with GenBank, but if hasn't been updated, email me. 

System requirements: [R](https://cran.r-project.org/), [git](https://git-scm.com/), [hmmer](http://hmmer.org/), [mafft](https://mafft.cbrc.jp/alignment/software/) and [raxml-ng](https://github.com/amkozlov/raxml-ng) need to be installed on your system, and available on your [$PATH](https://www.howtogeek.com/658904/how-to-add-a-directory-to-your-path-in-linux/). With the exception of raxml-ng, the programs can be installed from Ubuntu repositories using `sudo apt install`. In case of difficulties, check the developer's website and update to newer versions if required. Unfortunately, these scripts are optimised for a Unix system, and I'm unable to offer any Windows support here ([Windows is now able to run Ubuntu Linux ](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0)).

You will also require an API key from NCBI in order to access GenBank data at a decent rate. See info here for how to get a key: [ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)

##### Bash terminal:

```bash
### admin - clone the repository and create temporary directories ###
git clone https://github.com/genner-lab/meta-fish-lib.git meta-fish-lib
cd meta-fish-lib
mkdir -p reports temp/fasta-temp

### create NCBI API key (this is ignored by git) ###
# substitute the 'my-ncbi-key' part for your actual key obtained from NCBI
echo 'ncbi.key <- "my-ncbi-key"' > assets/ncbi-key.R

### install required R packages using renv (only need to run this once) ###
Rscript -e "renv::restore()"

### check the genbank versions remotely (github), locally (your machine), and on genbank itself ###
# you can update if the remote or local versions are behind genbank
scripts/check-genbank.R

### replace the species table with your custom list ###
# change to file location on your machine
# this will overwrite the current table 
cp ~/path/to/my-species-table.csv assets/species-table.csv

### search GenBank ###
# argument "-q" [2000] is the max search batch query length (in characters)
#     the maximum search query length allowed by NCBI is 2500
#     bigger queries will be faster, but more prone to server errors
# argument "-t" [4] is the number of processing threads to run in parallel
#     more threads are faster, but more prone to errors by overloading the server with requests
#     most laptops and desktops are hyperthreaded with two virtual CPUs (threads) for each physical CPU (cores)
#     run "lscpu | grep -E '^Thread|^Core|^Socket|^CPU\('" to obtain info on available cores/threads
#     make sure not to request more threads than are present on your system
# argument "-e" [true] is to run an exhaustive ("true") or simple search ("false")
#     the simple search just uses the terms "mitochondrion,mitochondrial"  
#     the exhaustive search in addition uses "COI,CO1,cox1,cytb,cytochrome,subunit,COB,CYB,12S,16S,rRNA,ribosomal"
#     the simple search will pick up 99% of mtDNA sequences, but may miss older sequences that were not well annotated
#     the simple search is faster, less prone to error, and may help if you need to search for a lot of species 
scripts/sequences-download.R -q 2000 -t 4 -e true

### assemble the reference library with hidden Markov models and obtain metadata ###
# argument "-t" [4] is the number of processing threads to run in parallel
#     do not request more threads than metabarcode markers
#     do not request more threads than are present on your system
# argument "-m" [all] is the metabarcode marker, here choosing "all" eight available metabarcodes
#     if you don't require all metabarcodes, it's strongly recommended to specify only the one(s) you want
#     to choose a specific metabarcode marker(s), use the codes in Table 1 and separate with a comma and no space
#     e.g. "-m 12s.miya,coi.ward"
# note this script overwrites the local master reference library 'assets/reference-library-master.csv.gz'
scripts/references-assemble.R -t 4 -m all

### phylogenetic quality control (QC) ###
# argument "-t" [4] is the number of processing threads to run in parallel
#    note that the threads argument doesn't influence the performance operation of raxml-ng or mafft
#    these applications automatically determine the optimum number of threads for your system
#    the argument only applies to preparing/plotting
#    do not request more threads than metabarcode markers
#    do not request more threads than are present on your system
# argument "-v" [false] is verbosity (information printed to screen) from the alignment and phylogenetic steps
#    a value of "true" prints program info, a value of "false" prevents printing of info
#    seeing the output on screen is useful if you run into problems and need to debug, otherwise it's not required
scripts/qc.R -t 4 -v false

# now manually review the phylogenetic tree PDFs output into 'reports/qc_GBVERSION_MONTH-YEAR' 
# if found, add suspect accessions manually to 'assets/exclusions.csv'

### compile species coverage reports ###
make -f scripts/Makefile

### update GitHub repository (updates remote version - only works if you have forked the repository rather than cloned) ###
# change x to actual version
cp reports/reports-tables.md assets/reports-tables.md
git add assets/reference-library-master.csv.gz assets/exclusions.csv assets/reports-tables.md
git commit -m "updated master to genbank version x"
git push
git tag -a vX -m "GenBank vX"
git push --tags
# now make a release using this tag on GitHub
# if Zenodo has been set up correctly, a DOI should become available immediately

### confirm changes are made ###
scripts/check-genbank.R

### if all worked, you can clean up to save disk space ###
# be sure that you want to do this!
rm -r temp
```

### FAQ

* **How do I cite the reference library?** - Zenodo DOIs for each version are in see [Releases](https://github.com/genner-lab/meta-fish-lib/releases). An important note: the reference library and code presented here supercedes a previous iteration hosted at [github.com/boopsboops/reference-libraries](https://github.com/boopsboops/reference-libraries). The new version here starts at v241, but I have archived only the final reference library file (`assets/reference-library-master.csv.gz`) for the previous versions here also. Therefore, while the library files are here, the old code used to generate these libraries prior to v241 are not archived together with that library version here. You can also cite the Collins et al. (2021) _Journal of Fish Biology_ article ([https://doi.org/10.1111/jfb.14852](https://doi.org/10.1111/jfb.14852)) describing the software.
* **Can I make a reference library for fishes of my country/region?** - Yes, very easily. Just change the list of species in `assets/species-table.csv`. You can provide this list yourself, but make sure the format of the table is the same. If not interested in synonyms, you use the same species name for 'speciesName' and 'validName' and set 'status' set to "accepted name". The 'commonSpecies' field can be all set to TRUE if that is not of interest either, and the other information can be obtained from FishBase ('fbSpecCode' is the FishBase species code). Alternatively, follow the [tutorial here](assets/species-list-synonyms.md) to generate an annotated species/synonyms list using the [rfishbase](https://docs.ropensci.org/rfishbase/index.html) package from scratch.
* **What if I don't know which species I need?** - This is a common problem in diverse tropical regions where there is often poorly resolved taxonomy and lots of undescribed species. Here you will want to search for genera instead of species. Copy this format below for the `assets/species-table.csv` table:

speciesName | status | fbSpecCode | validName | class | order | family | genus | commonName | commonSpecies
----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | -----
Anguilla | accepted name | NA | Anguilla | Actinopterygii | Anguilliformes | Anguillidae | Anguilla | NA | NA

* **Can I make a reference library for a taxonomic group that isn't fishes?** - At the moment I  use the [rfishbase](https://docs.ropensci.org/rfishbase/index.html) package to generate higher taxonomic ranks and validate scientific names etc, because this is the best source of data for fishes. However, more general solutions could easily be employed using the [taxize](https://docs.ropensci.org/taxize/) package or [taxadb](https://docs.ropensci.org/taxadb/index.html) package, with minimal changes to the code. One thing to bear in mind is that the hidden Markov models were designed on fishes, and while these would probably work well for other vertebrates, they would need to be recreated for other groups (see below).
* **How do I get synonyms?** - I used the `synonyms()` function in the [rfishbase](https://docs.ropensci.org/rfishbase/index.html) package, but the [taxize](https://docs.ropensci.org/taxize/) package or [taxadb](https://docs.ropensci.org/taxadb/index.html) package would achieve similar results for other groups.
* **How do I generate the hidden Markov models?** - First I downloaded the fish mitochondrial genomes and annotations from Prof. Masaki Miya's MitoFish website at [mitofish.aori.u-tokyo.ac.jp/](http://mitofish.aori.u-tokyo.ac.jp/), and extracted the genes of interest and aligned them with mafft. Then I searched for primer sequences and cut out the fragments of interest using [Geneious](https://www.geneious.com/prime/) and exported as fasta. Then I ran the hmmer function `hmmbuild` to create the hidden Markov models. Unfortunately, I did not include the code to perform these steps as it is not really general enough to be useful (requires manual actions and checking). Please contact me if you need specific help with these.
* **What if I want more than fishes?** - Indeed, for many metabarcoding applications you would want to identify 'off-target' reads, so a wider reference library is required as a supplement to the one presented here. I use the NCBI RefSeq mitochondrial DNA database ([ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion)), which should have a sufficiently broad coverage to roughly classify most eukaryote mtDNA.
* **How does the exclusions blacklist file work?** -  The file `assets/exclusions.csv` is permanently available on this GitHub repository, and contains all the accessions that have been flagged by me as potentially erroneous, as part of the work on the UK fish reference library. New records are added manually each time a new GenBank version becomes available and the quality control steps are performed. When the `scripts/references-clean.R` script is run, the exclusion file is called and these blacklisted accessions are removed from the library. The user does not need to regenerate or interact with this exclusions file if they are simply wanting to use the UK reference library as provided. If the user wishes to create their own custom reference library then they have the option of tailoring the contents of this exclusions file to their own requirements by keeping, deleting, or adding accessions to it.
* **Is the reference library guaranteed error free?** - LOL, no! I have tried to curate a reliable reference library as best as I can. However, the phylogenetic quality control step is tedious and subjective and takes a lot of effort. Here, phylogenetic QC trees for each primer set need to be manually checked. To help with this tips are coloured by monophyly and haplotype sharing to visually assist identifying dubious accessions. This is a much easier task for loci where taxa are well differentiated and large numbers of sequences exist (such as for COI). It is not easy for ribosomal fragments with fewer informative nucleotides and fewer sequences. The choice of which accessions to blacklist in `assets/exclusions.csv` has been entirely at my discretion thus far. However, I hope I have caught the majority of the most egregious examples. As a rule of thumb, an accession is blacklisted if: (a) individual(s) of species x are identical to or nested within a cluster of sequences of species y, but with other individuals of species x forming an independent cluster; and (b) the putatively spurious sequences coming from a single study, while the putatively correct sequences of species x and y coming from multiple studies. It is important to note that this is far from foolproof, and many species will naturally be non-monophyletic and/or share haplotypes with other species. I tried to be conservative, and not remove too many sequences if there was doubt, and especially so for taxa that I am not familar with. Mistakes certainly remain, so I recommend running the QC step to check yourself (or email me for the trees). [NCBI blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is also useful for checking for matches against species not in the reference library. 
* **The search step takes too long, hangs, or errors!** - The GenBank search step relies on NCBI servers, and if they are overloaded then the searches can fail. I suggest: (i) reducing the number of threads ("-t" option in the `scripts/sequences-download.R` step) to lower the frequency of requests; (ii) running searches at USA off-peak times; (iii) disabling the exhaustive search option (use "-e false" instead of "-e true" in the `scripts/sequences-download.R` step), which searches for fewer search terms, will take less time and consume less RAM and disk space, but will give you 99% of the sequences; (iv) making each concatenated search string length shorter with the "-q" option; and (v) requesting only the metabarcodes that you are interested in (use the "-m" option in the `scripts/references-assemble.R` step).
* **The phylogenetic quality control steps takes too long!** - Making ML trees for many taxa can take a very long time. Here, the largest one (Ward COI) is over 9,000 haplotypes, and runs overnight. If your dataset is too big, I suggest: (i) skipping this step if you aren't sure you need it; (ii) requesting only the metabarcodes that you are interested in (use the "-m" option in the `scripts/references-assemble.R` step); or (iii) maybe break up the species input list into smaller chunks and merge the tables later.
* **Why not use [sativa](https://github.com/amkozlov/sativa) for automated quality control?** - Good question.  Software such as [sativa](https://github.com/amkozlov/sativa) is available to automate the process, and while I may investigate this option in the future, for the meantime I think it is always a good idea to eyeball and become familiar with your data and develop an informed judgement.
* **Why are the sequence labels in the `references.fasta` file just numbers?** - When you download the reference library as shown above, the `references.fasta` file will use the 'dbid' column which is the database identification numbers. For NCBI these are 'GI' numbers (GenInfo Identifiers); these are equivalent to NCBI accession numbers, and will resolve accordingly on NCBI services; for BOLD, these are the 'processid' numbers. As there are many possible formats required for various taxonomy assignment software, I am unable to know which ones you will require, and have therefore chosen a sensible default. To make your own custom labels, just use the dplyr `mutate()` and `paste()` functions to join columns in the table to make a new label column. Below I make a label of format 'dbid_family_genus_species'. Table 2 explains the fields in the reference library table. 

```r
reflib.label <- reflib.sub %>% 
    mutate(label=paste(dbid,family,str_replace_all(sciNameValid," ","_"),sep="_"))
reflib.fas <- tab2fas(df=reflib.label,seqcol="nucleotides", namecol="label")
```

**Table 2: Key to reference library table fields**

Field (column name) | Description
----- | -----
source | source of record (GenBank or BOLD)
dbid | GenBank or BOLD database ID (GI, processid)
gbAccession | GenBank accession
sciNameValid | FishBase validated scientific name 
subphylum | taxonomic subphylum
class | FishBase taxonomic class
order | FishBase taxonomic order
family | FishBase taxonomic family
genus | FishBase taxonomic genus
sciNameOrig | original scientific name from GenBank/BOLD
fbSpecCode | FishBase species code
country | sample voucher collection country
catalogNumber | sample voucher catalogue number
institutionCode | sample voucher institution code
decimalLatitude | sample voucher collection latitude (decimal degrees)
decimalLongitude | sample voucher collection longitude (decimal degrees)
publishedAs | publication title
publishedIn | publication journal
publishedBy | publication lead author
date | date of sequence publication
notesGenBank | title of GenBank record
genbankVersion | version of GenBank used to generate this reference library
searchDate | date of reference library search
length | number of nucleotides in full record
nucleotides | nucleotides for full record
nucleotidesFrag.GENE.FRAGMENT.noprimers | nucleotides for gene fragment primer subset
lengthFrag.GENE.FRAGMENT.noprimers | number nucleotides in gene fragment primer subset


### Contents (A-Z)

* **`assets/`** - Required file and reference library.
    - **`hmms/`** - Hidden Markov models (HMMs) of gene markers of interest.
    - `exclusions.csv` - unreliable accessions to be excluded 
    - `logo.svg` - project logo
    - `reference-library-master.csv.gz` - master copy of the reference library
    - `reports-tables.md` - species coverage reports
    - `species-list-synonyms.md` - tutorial and R code to generate species lists and synonyms
    - `species-table.csv` - list of species to search for
    - `species-table-testing.csv` - a test list of goby species
* **`renv/`** - Settings for the R environment.
* **`reports/`** - Location of QC reports. Temporary directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.
* **`scripts/`** - R and shell scripts.
    - `check-genbank.R` - script to get genbank versions
    - `load-libs.R` - script to load all required packages and functions
    - `Makefile` - makefile to generate the species coverage reports
    - `qc.R` - quality control a reference library
    - `references-assemble.R` - extract and annotate reference libraries from ncbi/bold dumps 
    - `references-clean.R` - quality filter reference library
    - `references-load-local.R` - load reference library locally
    - `references-load-remote.R` - load reference library remotely
    - `reports-tables.Rmd` - knitr file to prepare species coverage reports
    - `sequences-download.R` - pulls all mitochondrial DNA from NCBI and BOLD for a list of species
* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.
Generating species and synonym lists
================
Rupert A. Collins
18 January 2021


Here is a tutorial for obtaining, annotating, cleaning and formatting a fish species list for a given country. We use the [rfishbase](https://docs.ropensci.org/rfishbase/index.html) package to generate both the list and then augment with data from its numerous tables. The [rgbif](https://docs.ropensci.org/rgbif/index.html) package is also very useful for generating species lists for particular regions.


### Species lists

Use ISO country codes to define countries. These can be found at [en.wikipedia.org/wiki/List_of_ISO_3166_country_codes](https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes). 

```r
# load libs
library("tidyverse")
library("rfishbase")

# load countries
# filter on country of interest - use ISO country codes
# ISO country code "826" is "United Kingdom"
# also remove subspecific names
species.list <- rfishbase::country(server="fishbase") %>% 
    filter(C_Code=="826") %>% 
    mutate(Species=paste(str_split_fixed(Species," ",3)[,1],str_split_fixed(Species," ",3)[,2])) %>%
    distinct(SpecCode,Species)
```


### Synonyms

Now we get synonyms for all species in our list. We also clean the table to remove those with unusual, non-alphabetic characters, indicating the name is not a standard format bionomial/trinomial species name. We also remove other types of names (e.g. misapplied names), and also remove any duplicates with >1 accepted scientific name.

```r
# load synonyms
# clean up to get only synonyms and accepted names
# remove records with non alphabetic characters
# remove species with >1 accepted names
fishbase.synonyms.clean <- rfishbase::synonyms(server="fishbase") %>% 
    select(synonym,Status,SpecCode,TaxonLevel) %>% 
    mutate(Status=str_replace_all(Status,"Synonym","synonym")) %>% 
    filter(Status=="synonym" | Status=="accepted name") %>% 
    filter(!str_detect(synonym,"[^a-zA-Z\\d\\s:]")) %>%
    group_by(SpecCode) %>% 
    mutate(nacc=length(unique(synonym[Status=="accepted name"]))) %>%
    ungroup() %>%
    mutate(dup=if_else(nacc>1 & TaxonLevel=="Nominotypical" & Status=="accepted name",TRUE,FALSE)) %>%
    filter(dup==FALSE) %>%
    select(!dup)

# join the countries and synonyms tables
species.list.syn <- species.list %>% left_join(distinct(fishbase.synonyms.clean,synonym,Status,SpecCode))
```


### Taxonomy and common names

Now we get higher taxonomy and common names and add these fields to the species list.

```r
# load taxonomy and common name tables
fishbase.taxonomy <- rfishbase::load_taxa(server="fishbase")
fishbase.species <- rfishbase::species(server="fishbase")

# add the taxonomy
species.list.tax <- species.list.syn %>% left_join(distinct(fishbase.taxonomy,SpecCode,Genus,Family,Order,Class))

# add the common names
species.list.com <- species.list.tax %>% left_join(distinct(fishbase.species,SpecCode,Species,FBname))
```


### Format and write out

Now we format the species list in the same way as the species list file (`assets/species-table.csv`) in order for it to work in the pipeline. Then we write it out to a CSV formatted file (move it to `assets/species-table.csv` when you are happy with it). The 'commonSpecies' field was set as all 'TRUE', so it you wish to use this feature, then manually change these values in the table at your own discretion.

```r
# format
species.list.form <- species.list.com %>% 
    rename(speciesName=synonym,status=Status,fbSpecCode=SpecCode,validName=Species,class=Class,order=Order,family=Family,genus=Genus,commonName=FBname) %>% 
    mutate(commonSpecies=TRUE) %>%
    relocate(speciesName,status,fbSpecCode,validName,class,order,family,genus,commonName,commonSpecies) %>% 
    arrange(class,order,family,genus,validName,status,speciesName) 

# write out
species.list.form %>% write_csv(file="species-table.csv")
```
Reference library coverage report
================
Rupert A. Collins
20 January 2022

##### Methods and description

This document describes the contents of the UK fish reference library, generated from public databases. The document is a dynamic knitr document and can be updated quickly using the Makefile in `scripts/`. A list of species from the UK was generated from three sources: GBIF, FishBase, and the Water Framework Directive list of transitional species. This list was filtered to identify synonyms and duplicates, and annotated with FishBase taxonomic classification and FishBase common names. Next a sub-list of "common" species was generated. These were species that we believe are likely to be encountered in eDNA surveys of inshore and transitional waters of the UK, and comprise most of the species in Henderson (2015). Most of the remaining are either introduced species, rarely encountered migrants, oceanic pelagics, or deep sea organisms.

The search was performed on the NCBI nucleotide and BOLD sequences databases. Because of inconsistencies in how researchers annotate their GenBank submissions and the differing internal coverage of primer pairs for particular gene fragments, we performed a search requesting mitochondrial DNA using multiple search relevant search terms (COI, 12S, 16S, rRNA, ribosomal, cytb, CO1, cox1, cytochrome, subunit, COB, CYB, mitochondrial, mitochondrion). Then we pulled out fragments of interest using a hidden Markov model. This enabled us to have greater confidence that useful sequences had not been missed. For the resulting sequences we then tabulate all their metadata from GenBank in order to allow us the capability to later tailor a custom reference library according to any criteria required (e.g. must have reference specimen or locality data etc).

##### Results

The total number of accepted UK species is estimated to be around 531, with 176 common species, and 4280 total names including synonyms. The NCBI GenBank and BOLD databases were searched on 13 Jan 2022 (GenBank version 247), and the search retrieved 52970 accessions from 496 unique species. Below is presented a summary table of reference library coverage (Table 1), numbers of individuals per common species (Table 2), and the sequences added to the reference library in the most recent update (Table 3).

**Table 1. Summary of coverage. Locus = mitochondrial gene; Fragment = metabarcode primer set; Total = total number of sequences; Cov. (all) = proportion of all species with at least one sequence; Cov. (common) = proportion of common species with at least one sequence; Cov. (rare) = proportion of rare species with at least one sequence; Singletons = proportion of species represented by only one sequence, only including those with &gt;0 sequences; Haps (mean) = mean number unique haplotypes per species; Haps (median) = median number unique haplotypes per species.**

<table>
<colgroup>
<col width="6%" />
<col width="10%" />
<col width="6%" />
<col width="11%" />
<col width="14%" />
<col width="12%" />
<col width="11%" />
<col width="12%" />
<col width="14%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Locus</th>
<th align="left">Fragment</th>
<th align="right">Total</th>
<th align="right">Cov. (all)</th>
<th align="right">Cov. (common)</th>
<th align="right">Cov. (rare)</th>
<th align="right">Singletons</th>
<th align="right">Haps (mean)</th>
<th align="right">Haps (median)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">12S</td>
<td align="left">Miya</td>
<td align="right">2813</td>
<td align="right">0.77</td>
<td align="right">0.97</td>
<td align="right">0.68</td>
<td align="right">0.17</td>
<td align="right">1.8</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">12S</td>
<td align="left">Riaz</td>
<td align="right">3112</td>
<td align="right">0.74</td>
<td align="right">0.93</td>
<td align="right">0.65</td>
<td align="right">0.22</td>
<td align="right">1.8</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">12S</td>
<td align="left">Taberlet</td>
<td align="right">2813</td>
<td align="right">0.77</td>
<td align="right">0.97</td>
<td align="right">0.68</td>
<td align="right">0.17</td>
<td align="right">1.8</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">12S</td>
<td align="left">Valentini</td>
<td align="right">2020</td>
<td align="right">0.61</td>
<td align="right">0.73</td>
<td align="right">0.54</td>
<td align="right">0.28</td>
<td align="right">1.3</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">16S</td>
<td align="left">Berry</td>
<td align="right">4948</td>
<td align="right">0.79</td>
<td align="right">0.97</td>
<td align="right">0.70</td>
<td align="right">0.14</td>
<td align="right">3.5</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td align="left">COI</td>
<td align="left">Lerayxt</td>
<td align="right">30853</td>
<td align="right">0.92</td>
<td align="right">0.98</td>
<td align="right">0.89</td>
<td align="right">0.02</td>
<td align="right">11.2</td>
<td align="right">7</td>
</tr>
<tr class="odd">
<td align="left">COI</td>
<td align="left">Ward</td>
<td align="right">31029</td>
<td align="right">0.92</td>
<td align="right">0.98</td>
<td align="right">0.89</td>
<td align="right">0.02</td>
<td align="right">19.1</td>
<td align="right">10</td>
</tr>
<tr class="even">
<td align="left">CYTB</td>
<td align="left">Minamoto</td>
<td align="right">17773</td>
<td align="right">0.69</td>
<td align="right">0.89</td>
<td align="right">0.59</td>
<td align="right">0.15</td>
<td align="right">8.7</td>
<td align="right">2</td>
</tr>
</tbody>
</table>

**Table 2. Numbers of sequences represented per species for each primer set metabarcode fragment. Only common species are shown.**

<table>
<colgroup>
<col width="8%" />
<col width="16%" />
<col width="13%" />
<col width="6%" />
<col width="6%" />
<col width="8%" />
<col width="8%" />
<col width="6%" />
<col width="7%" />
<col width="6%" />
<col width="8%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Family</th>
<th align="left">Scientific Name</th>
<th align="left">Common Name</th>
<th align="right">12S (Miya)</th>
<th align="right">12S (Riaz)</th>
<th align="right">12S (Taberlet)</th>
<th align="right">12S (Valentini)</th>
<th align="right">16S (Berry)</th>
<th align="right">COI (Lerayxt)</th>
<th align="right">COI (Ward)</th>
<th align="right">CYTB (Minamoto)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Anguillidae</td>
<td align="left"><em>Anguilla anguilla</em></td>
<td align="left">European eel</td>
<td align="right">63</td>
<td align="right">66</td>
<td align="right">63</td>
<td align="right">59</td>
<td align="right">231</td>
<td align="right">234</td>
<td align="right">234</td>
<td align="right">141</td>
</tr>
<tr class="even">
<td align="left">Congridae</td>
<td align="left"><em>Conger conger</em></td>
<td align="left">European conger</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">34</td>
<td align="right">34</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Atherinidae</td>
<td align="left"><em>Atherina boyeri</em></td>
<td align="left">Big-scale sand smelt</td>
<td align="right">10</td>
<td align="right">169</td>
<td align="right">10</td>
<td align="right">1</td>
<td align="right">68</td>
<td align="right">29</td>
<td align="right">29</td>
<td align="right">103</td>
</tr>
<tr class="even">
<td align="left">Belonidae</td>
<td align="left"><em>Belone belone</em></td>
<td align="left">Garfish</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">6</td>
<td align="right">43</td>
<td align="right">43</td>
<td align="right">9</td>
</tr>
<tr class="odd">
<td align="left">Clupeidae</td>
<td align="left"><em>Alosa alosa</em></td>
<td align="left">Allis shad</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">32</td>
<td align="right">32</td>
<td align="right">197</td>
</tr>
<tr class="even">
<td align="left">Clupeidae</td>
<td align="left"><em>Alosa fallax</em></td>
<td align="left">Twaite shad</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">10</td>
<td align="right">18</td>
<td align="right">18</td>
<td align="right">520</td>
</tr>
<tr class="odd">
<td align="left">Clupeidae</td>
<td align="left"><em>Clupea harengus</em></td>
<td align="left">Atlantic herring</td>
<td align="right">107</td>
<td align="right">115</td>
<td align="right">107</td>
<td align="right">109</td>
<td align="right">120</td>
<td align="right">191</td>
<td align="right">191</td>
<td align="right">196</td>
</tr>
<tr class="even">
<td align="left">Clupeidae</td>
<td align="left"><em>Sardina pilchardus</em></td>
<td align="left">European pilchard</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">11</td>
<td align="right">4</td>
<td align="right">34</td>
<td align="right">233</td>
<td align="right">238</td>
<td align="right">99</td>
</tr>
<tr class="odd">
<td align="left">Clupeidae</td>
<td align="left"><em>Sprattus sprattus</em></td>
<td align="left">European sprat</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">8</td>
<td align="right">3</td>
<td align="right">15</td>
<td align="right">75</td>
<td align="right">75</td>
<td align="right">20</td>
</tr>
<tr class="even">
<td align="left">Engraulidae</td>
<td align="left"><em>Engraulis encrasicolus</em></td>
<td align="left">European anchovy</td>
<td align="right">6</td>
<td align="right">20</td>
<td align="right">6</td>
<td align="right">8</td>
<td align="right">75</td>
<td align="right">265</td>
<td align="right">265</td>
<td align="right">1016</td>
</tr>
<tr class="odd">
<td align="left">Cobitidae</td>
<td align="left"><em>Cobitis taenia</em></td>
<td align="left">Spined loach</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">21</td>
<td align="right">21</td>
<td align="right">189</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Abramis brama</em></td>
<td align="left">Freshwater bream</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">91</td>
<td align="right">91</td>
<td align="right">27</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Alburnus alburnus</em></td>
<td align="left">Bleak</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">10</td>
<td align="right">205</td>
<td align="right">206</td>
<td align="right">37</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Barbus barbus</em></td>
<td align="left">Barbel</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">6</td>
<td align="right">132</td>
<td align="right">132</td>
<td align="right">80</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Blicca bjoerkna</em></td>
<td align="left">White bream</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">10</td>
<td align="right">53</td>
<td align="right">53</td>
<td align="right">8</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Carassius auratus</em></td>
<td align="left">Goldfish</td>
<td align="right">55</td>
<td align="right">47</td>
<td align="right">55</td>
<td align="right">44</td>
<td align="right">78</td>
<td align="right">305</td>
<td align="right">305</td>
<td align="right">721</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Carassius carassius</em></td>
<td align="left">Crucian carp</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">43</td>
<td align="right">43</td>
<td align="right">128</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Cyprinus carpio</em></td>
<td align="left">Common carp</td>
<td align="right">71</td>
<td align="right">54</td>
<td align="right">71</td>
<td align="right">46</td>
<td align="right">122</td>
<td align="right">708</td>
<td align="right">709</td>
<td align="right">158</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Gobio gobio</em></td>
<td align="left">Gudgeon</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">12</td>
<td align="right">128</td>
<td align="right">128</td>
<td align="right">45</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Leuciscus idus</em></td>
<td align="left">Ide</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">6</td>
<td align="right">57</td>
<td align="right">57</td>
<td align="right">9</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Leuciscus leuciscus</em></td>
<td align="left">Common dace</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">47</td>
<td align="right">151</td>
<td align="right">151</td>
<td align="right">154</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Phoxinus phoxinus</em></td>
<td align="left">Eurasian minnow</td>
<td align="right">19</td>
<td align="right">17</td>
<td align="right">19</td>
<td align="right">17</td>
<td align="right">21</td>
<td align="right">753</td>
<td align="right">753</td>
<td align="right">359</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Pseudorasbora parva</em></td>
<td align="left">Stone moroko</td>
<td align="right">14</td>
<td align="right">9</td>
<td align="right">14</td>
<td align="right">8</td>
<td align="right">47</td>
<td align="right">264</td>
<td align="right">264</td>
<td align="right">1129</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Rutilus rutilus</em></td>
<td align="left">Roach</td>
<td align="right">3</td>
<td align="right">7</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">15</td>
<td align="right">111</td>
<td align="right">111</td>
<td align="right">274</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Scardinius erythrophthalmus</em></td>
<td align="left">Rudd</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">7</td>
<td align="right">92</td>
<td align="right">92</td>
<td align="right">27</td>
</tr>
<tr class="even">
<td align="left">Cyprinidae</td>
<td align="left"><em>Squalius cephalus</em></td>
<td align="left">Chub</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">8</td>
<td align="right">265</td>
<td align="right">265</td>
<td align="right">165</td>
</tr>
<tr class="odd">
<td align="left">Cyprinidae</td>
<td align="left"><em>Tinca tinca</em></td>
<td align="left">Tench</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">12</td>
<td align="right">172</td>
<td align="right">172</td>
<td align="right">35</td>
</tr>
<tr class="even">
<td align="left">Nemacheilidae</td>
<td align="left"><em>Barbatula barbatula</em></td>
<td align="left">Stone loach</td>
<td align="right">9</td>
<td align="right">8</td>
<td align="right">9</td>
<td align="right">8</td>
<td align="right">11</td>
<td align="right">135</td>
<td align="right">135</td>
<td align="right">344</td>
</tr>
<tr class="odd">
<td align="left">Esocidae</td>
<td align="left"><em>Esox lucius</em></td>
<td align="left">Northern pike</td>
<td align="right">9</td>
<td align="right">18</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">31</td>
<td align="right">194</td>
<td align="right">195</td>
<td align="right">325</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Gadiculus argenteus</em></td>
<td align="left">Silvery pout</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">7</td>
<td align="right">38</td>
<td align="right">38</td>
<td align="right">11</td>
</tr>
<tr class="odd">
<td align="left">Gadidae</td>
<td align="left"><em>Gadus morhua</em></td>
<td align="left">Atlantic cod</td>
<td align="right">172</td>
<td align="right">178</td>
<td align="right">172</td>
<td align="right">170</td>
<td align="right">182</td>
<td align="right">520</td>
<td align="right">521</td>
<td align="right">1099</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Melanogrammus aeglefinus</em></td>
<td align="left">Haddock</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">12</td>
<td align="right">246</td>
<td align="right">247</td>
<td align="right">44</td>
</tr>
<tr class="odd">
<td align="left">Gadidae</td>
<td align="left"><em>Merlangius merlangus</em></td>
<td align="left">Whiting</td>
<td align="right">10</td>
<td align="right">12</td>
<td align="right">10</td>
<td align="right">5</td>
<td align="right">21</td>
<td align="right">109</td>
<td align="right">110</td>
<td align="right">41</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Micromesistius poutassou</em></td>
<td align="left">Blue whiting</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">8</td>
<td align="right">111</td>
<td align="right">111</td>
<td align="right">22</td>
</tr>
<tr class="odd">
<td align="left">Gadidae</td>
<td align="left"><em>Pollachius pollachius</em></td>
<td align="left">Pollack</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">17</td>
<td align="right">17</td>
<td align="right">12</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Pollachius virens</em></td>
<td align="left">Saithe</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">85</td>
<td align="right">85</td>
<td align="right">27</td>
</tr>
<tr class="odd">
<td align="left">Gadidae</td>
<td align="left"><em>Raniceps raninus</em></td>
<td align="left">Tadpole fish</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Trisopterus esmarkii</em></td>
<td align="left">Norway pout</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">36</td>
<td align="right">36</td>
<td align="right">10</td>
</tr>
<tr class="odd">
<td align="left">Gadidae</td>
<td align="left"><em>Trisopterus luscus</em></td>
<td align="left">Pouting</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">25</td>
<td align="right">25</td>
<td align="right">20</td>
</tr>
<tr class="even">
<td align="left">Gadidae</td>
<td align="left"><em>Trisopterus minutus</em></td>
<td align="left">Poor cod</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">2</td>
<td align="right">9</td>
<td align="right">31</td>
<td align="right">31</td>
<td align="right">42</td>
</tr>
<tr class="odd">
<td align="left">Lotidae</td>
<td align="left"><em>Ciliata mustela</em></td>
<td align="left">Fivebeard rockling</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">22</td>
<td align="right">22</td>
<td align="right">92</td>
</tr>
<tr class="even">
<td align="left">Lotidae</td>
<td align="left"><em>Ciliata septentrionalis</em></td>
<td align="left">Northern rockling</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Lotidae</td>
<td align="left"><em>Enchelyopus cimbrius</em></td>
<td align="left">Fourbeard rockling</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">39</td>
<td align="right">39</td>
<td align="right">12</td>
</tr>
<tr class="even">
<td align="left">Lotidae</td>
<td align="left"><em>Gaidropsarus mediterraneus</em></td>
<td align="left">Shore rockling</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">23</td>
<td align="right">23</td>
<td align="right">7</td>
</tr>
<tr class="odd">
<td align="left">Lotidae</td>
<td align="left"><em>Gaidropsarus vulgaris</em></td>
<td align="left">Three-bearded rockling</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">25</td>
<td align="right">25</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left">Lotidae</td>
<td align="left"><em>Molva molva</em></td>
<td align="left">Ling</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">22</td>
<td align="right">22</td>
<td align="right">21</td>
</tr>
<tr class="odd">
<td align="left">Merlucciidae</td>
<td align="left"><em>Merluccius merluccius</em></td>
<td align="left">European hake</td>
<td align="right">4</td>
<td align="right">10</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">33</td>
<td align="right">341</td>
<td align="right">341</td>
<td align="right">99</td>
</tr>
<tr class="even">
<td align="left">Gasterosteidae</td>
<td align="left"><em>Gasterosteus aculeatus</em></td>
<td align="left">Three-spined stickleback</td>
<td align="right">18</td>
<td align="right">18</td>
<td align="right">18</td>
<td align="right">13</td>
<td align="right">24</td>
<td align="right">294</td>
<td align="right">294</td>
<td align="right">489</td>
</tr>
<tr class="odd">
<td align="left">Gasterosteidae</td>
<td align="left"><em>Pungitius pungitius</em></td>
<td align="left">Ninespine stickleback</td>
<td align="right">11</td>
<td align="right">29</td>
<td align="right">11</td>
<td align="right">11</td>
<td align="right">48</td>
<td align="right">164</td>
<td align="right">164</td>
<td align="right">321</td>
</tr>
<tr class="even">
<td align="left">Gasterosteidae</td>
<td align="left"><em>Spinachia spinachia</em></td>
<td align="left">Sea stickleback</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left">Gobiesocidae</td>
<td align="left"><em>Apletodon dentatus</em></td>
<td align="left">Small-headed clingfish</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Gobiesocidae</td>
<td align="left"><em>Diplecogaster bimaculata</em></td>
<td align="left">Two-spotted clingfish</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Gobiesocidae</td>
<td align="left"><em>Lepadogaster candolii</em></td>
<td align="left">Connemarra clingfish</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Gobiesocidae</td>
<td align="left"><em>Lepadogaster purpurea</em></td>
<td align="left">Cornish sucker</td>
<td align="right">2</td>
<td align="right">12</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Lophiidae</td>
<td align="left"><em>Lophius piscatorius</em></td>
<td align="left">Angler</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">11</td>
<td align="right">90</td>
<td align="right">90</td>
<td align="right">24</td>
</tr>
<tr class="even">
<td align="left">Mugilidae</td>
<td align="left"><em>Chelon auratus</em></td>
<td align="left">Golden grey mullet</td>
<td align="right">4</td>
<td align="right">8</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">19</td>
<td align="right">64</td>
<td align="right">64</td>
<td align="right">18</td>
</tr>
<tr class="odd">
<td align="left">Mugilidae</td>
<td align="left"><em>Chelon labrosus</em></td>
<td align="left">Thicklip grey mullet</td>
<td align="right">7</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">15</td>
<td align="right">28</td>
<td align="right">28</td>
<td align="right">13</td>
</tr>
<tr class="even">
<td align="left">Mugilidae</td>
<td align="left"><em>Chelon ramada</em></td>
<td align="left">Thinlip grey mullet</td>
<td align="right">7</td>
<td align="right">10</td>
<td align="right">7</td>
<td align="right">1</td>
<td align="right">12</td>
<td align="right">48</td>
<td align="right">48</td>
<td align="right">10</td>
</tr>
<tr class="odd">
<td align="left">Osmeridae</td>
<td align="left"><em>Osmerus eperlanus</em></td>
<td align="left">European smelt</td>
<td align="right">4</td>
<td align="right">8</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">10</td>
<td align="right">36</td>
<td align="right">36</td>
<td align="right">21</td>
</tr>
<tr class="even">
<td align="left">Ammodytidae</td>
<td align="left"><em>Ammodytes marinus</em></td>
<td align="left">Lesser sand-eel</td>
<td align="right">5</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">52</td>
<td align="right">52</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">Ammodytidae</td>
<td align="left"><em>Ammodytes tobianus</em></td>
<td align="left">Small sandeel</td>
<td align="right">7</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">13</td>
<td align="right">13</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Ammodytidae</td>
<td align="left"><em>Gymnammodytes semisquamatus</em></td>
<td align="left">Smooth sandeel</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Ammodytidae</td>
<td align="left"><em>Hyperoplus immaculatus</em></td>
<td align="left">Greater sand-eel</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Ammodytidae</td>
<td align="left"><em>Hyperoplus lanceolatus</em></td>
<td align="left">Great sandeel</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">39</td>
<td align="right">39</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Anarhichadidae</td>
<td align="left"><em>Anarhichas lupus</em></td>
<td align="left">Atlantic wolffish</td>
<td align="right">90</td>
<td align="right">89</td>
<td align="right">90</td>
<td align="right">89</td>
<td align="right">93</td>
<td align="right">196</td>
<td align="right">196</td>
<td align="right">95</td>
</tr>
<tr class="even">
<td align="left">Blenniidae</td>
<td align="left"><em>Blennius ocellaris</em></td>
<td align="left">Butterfly blenny</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">24</td>
<td align="right">24</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Blenniidae</td>
<td align="left"><em>Coryphoblennius galerita</em></td>
<td align="left">Montagu's blenny</td>
<td align="right">2</td>
<td align="right">62</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">69</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Blenniidae</td>
<td align="left"><em>Lipophrys pholis</em></td>
<td align="left">Shanny</td>
<td align="right">5</td>
<td align="right">14</td>
<td align="right">5</td>
<td align="right">0</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Blenniidae</td>
<td align="left"><em>Parablennius gattorugine</em></td>
<td align="left">Tompot blenny</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Callionymidae</td>
<td align="left"><em>Callionymus lyra</em></td>
<td align="left">Dragonet</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">31</td>
<td align="right">31</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Callionymidae</td>
<td align="left"><em>Callionymus maculatus</em></td>
<td align="left"></td>
<td align="right">5</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">12</td>
<td align="right">12</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Callionymidae</td>
<td align="left"><em>Callionymus reticulatus</em></td>
<td align="left">Reticulated dragonet</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Caproidae</td>
<td align="left"><em>Capros aper</em></td>
<td align="left">Boarfish</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">33</td>
<td align="right">33</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td align="left">Carangidae</td>
<td align="left"><em>Trachurus trachurus</em></td>
<td align="left">Atlantic horse mackerel</td>
<td align="right">8</td>
<td align="right">6</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">23</td>
<td align="right">172</td>
<td align="right">173</td>
<td align="right">35</td>
</tr>
<tr class="odd">
<td align="left">Cepolidae</td>
<td align="left"><em>Cepola macrophthalma</em></td>
<td align="left">Red bandfish</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">27</td>
<td align="right">27</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Aphia minuta</em></td>
<td align="left">Transparent goby</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">15</td>
<td align="right">18</td>
<td align="right">18</td>
<td align="right">12</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Crystallogobius linearis</em></td>
<td align="left">Crystal goby</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">12</td>
<td align="right">12</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Gobius cobitis</em></td>
<td align="left">Giant goby</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">11</td>
<td align="right">11</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Gobiusculus flavescens</em></td>
<td align="left">Two-spotted goby</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">6</td>
<td align="right">14</td>
<td align="right">14</td>
<td align="right">5</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Gobius niger</em></td>
<td align="left">Black goby</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">11</td>
<td align="right">68</td>
<td align="right">68</td>
<td align="right">5</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Gobius paganellus</em></td>
<td align="left">Rock goby</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">17</td>
<td align="right">17</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Lesueurigobius friesii</em></td>
<td align="left">Fries's goby</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Pomatoschistus lozanoi</em></td>
<td align="left">Lozano's goby</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">3</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Pomatoschistus microps</em></td>
<td align="left">Common goby</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">9</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">17</td>
<td align="right">17</td>
<td align="right">39</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Pomatoschistus minutus</em></td>
<td align="left">Sand goby</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">1</td>
<td align="right">7</td>
<td align="right">14</td>
<td align="right">14</td>
<td align="right">134</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Pomatoschistus norvegicus</em></td>
<td align="left">Norway goby</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Gobiidae</td>
<td align="left"><em>Pomatoschistus pictus</em></td>
<td align="left">Painted goby</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Gobiidae</td>
<td align="left"><em>Thorogobius ephippiatus</em></td>
<td align="left">Leopard-spotted goby</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">Labridae</td>
<td align="left"><em>Centrolabrus exoletus</em></td>
<td align="left">Rock cook</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">3</td>
</tr>
<tr class="even">
<td align="left">Labridae</td>
<td align="left"><em>Ctenolabrus rupestris</em></td>
<td align="left">Goldsinny-wrasse</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">11</td>
<td align="right">11</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left">Labridae</td>
<td align="left"><em>Labrus bergylta</em></td>
<td align="left">Ballan wrasse</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">163</td>
<td align="right">165</td>
<td align="right">5</td>
</tr>
<tr class="even">
<td align="left">Labridae</td>
<td align="left"><em>Labrus mixtus</em></td>
<td align="left">Cuckoo wrasse</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">28</td>
<td align="right">28</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Labridae</td>
<td align="left"><em>Symphodus bailloni</em></td>
<td align="left">Baillon's wrasse</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Labridae</td>
<td align="left"><em>Symphodus melops</em></td>
<td align="left">Corkwing wrasse</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Moronidae</td>
<td align="left"><em>Dicentrarchus labrax</em></td>
<td align="left">European seabass</td>
<td align="right">5</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">14</td>
<td align="right">65</td>
<td align="right">65</td>
<td align="right">27</td>
</tr>
<tr class="even">
<td align="left">Mullidae</td>
<td align="left"><em>Mullus surmuletus</em></td>
<td align="left">Surmullet</td>
<td align="right">6</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">21</td>
<td align="right">97</td>
<td align="right">97</td>
<td align="right">19</td>
</tr>
<tr class="odd">
<td align="left">Percidae</td>
<td align="left"><em>Gymnocephalus cernua</em></td>
<td align="left">Ruffe</td>
<td align="right">5</td>
<td align="right">8</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">12</td>
<td align="right">88</td>
<td align="right">88</td>
<td align="right">9</td>
</tr>
<tr class="even">
<td align="left">Percidae</td>
<td align="left"><em>Perca fluviatilis</em></td>
<td align="left">European perch</td>
<td align="right">21</td>
<td align="right">15</td>
<td align="right">21</td>
<td align="right">10</td>
<td align="right">41</td>
<td align="right">156</td>
<td align="right">156</td>
<td align="right">54</td>
</tr>
<tr class="odd">
<td align="left">Percidae</td>
<td align="left"><em>Sander lucioperca</em></td>
<td align="left">Pike-perch</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">7</td>
<td align="right">97</td>
<td align="right">97</td>
<td align="right">30</td>
</tr>
<tr class="even">
<td align="left">Pholidae</td>
<td align="left"><em>Pholis gunnellus</em></td>
<td align="left">Rock gunnel</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">27</td>
<td align="right">27</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td align="left">Scombridae</td>
<td align="left"><em>Scomber scombrus</em></td>
<td align="left">Atlantic mackerel</td>
<td align="right">6</td>
<td align="right">28</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">16</td>
<td align="right">430</td>
<td align="right">430</td>
<td align="right">78</td>
</tr>
<tr class="even">
<td align="left">Sparidae</td>
<td align="left"><em>Pagellus bogaraveo</em></td>
<td align="left">Blackspot seabream</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">20</td>
<td align="right">28</td>
<td align="right">16</td>
</tr>
<tr class="odd">
<td align="left">Sparidae</td>
<td align="left"><em>Pagrus pagrus</em></td>
<td align="left">Red porgy</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">68</td>
<td align="right">72</td>
<td align="right">13</td>
</tr>
<tr class="even">
<td align="left">Sparidae</td>
<td align="left"><em>Sparus aurata</em></td>
<td align="left">Gilthead seabream</td>
<td align="right">8</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">3</td>
<td align="right">15</td>
<td align="right">186</td>
<td align="right">186</td>
<td align="right">70</td>
</tr>
<tr class="odd">
<td align="left">Sparidae</td>
<td align="left"><em>Spondyliosoma cantharus</em></td>
<td align="left">Black seabream</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">47</td>
<td align="right">47</td>
<td align="right">77</td>
</tr>
<tr class="even">
<td align="left">Stichaeidae</td>
<td align="left"><em>Chirolophis ascanii</em></td>
<td align="left">Yarrell's blenny</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Stichaeidae</td>
<td align="left"><em>Lumpenus lampretaeformis</em></td>
<td align="left">Snakeblenny</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">33</td>
<td align="right">33</td>
<td align="right">65</td>
</tr>
<tr class="even">
<td align="left">Trachinidae</td>
<td align="left"><em>Echiichthys vipera</em></td>
<td align="left">Lesser weever</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">20</td>
<td align="right">20</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left">Trachinidae</td>
<td align="left"><em>Trachinus draco</em></td>
<td align="left">Greater weever</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">7</td>
<td align="right">34</td>
<td align="right">34</td>
<td align="right">7</td>
</tr>
<tr class="even">
<td align="left">Zoarcidae</td>
<td align="left"><em>Zoarces viviparus</em></td>
<td align="left">Eelpout</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">82</td>
</tr>
<tr class="odd">
<td align="left">Bothidae</td>
<td align="left"><em>Arnoglossus laterna</em></td>
<td align="left">Mediterranean scaldfish</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">0</td>
<td align="right">10</td>
<td align="right">60</td>
<td align="right">60</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Glyptocephalus cynoglossus</em></td>
<td align="left">Witch flounder</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">60</td>
<td align="right">60</td>
<td align="right">17</td>
</tr>
<tr class="odd">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Hippoglossoides platessoides</em></td>
<td align="left">American plaice</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">7</td>
<td align="right">62</td>
<td align="right">62</td>
<td align="right">16</td>
</tr>
<tr class="even">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Hippoglossus hippoglossus</em></td>
<td align="left">Atlantic halibut</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">27</td>
<td align="right">27</td>
<td align="right">15</td>
</tr>
<tr class="odd">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Limanda limanda</em></td>
<td align="left">Common dab</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">16</td>
<td align="right">45</td>
<td align="right">45</td>
<td align="right">78</td>
</tr>
<tr class="even">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Microstomus kitt</em></td>
<td align="left">Lemon sole</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">38</td>
<td align="right">38</td>
<td align="right">13</td>
</tr>
<tr class="odd">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Platichthys flesus</em></td>
<td align="left">European flounder</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">1</td>
<td align="right">20</td>
<td align="right">105</td>
<td align="right">105</td>
<td align="right">78</td>
</tr>
<tr class="even">
<td align="left">Pleuronectidae</td>
<td align="left"><em>Pleuronectes platessa</em></td>
<td align="left">European plaice</td>
<td align="right">5</td>
<td align="right">7</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">17</td>
<td align="right">99</td>
<td align="right">99</td>
<td align="right">108</td>
</tr>
<tr class="odd">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Lepidorhombus whiffiagonis</em></td>
<td align="left">Megrim</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">34</td>
<td align="right">34</td>
<td align="right">11</td>
</tr>
<tr class="even">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Phrynorhombus norvegicus</em></td>
<td align="left">Norwegian topknot</td>
<td align="right">3</td>
<td align="right">1</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">5</td>
</tr>
<tr class="odd">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Scophthalmus maximus</em></td>
<td align="left">Turbot</td>
<td align="right">3</td>
<td align="right">6</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">18</td>
<td align="right">108</td>
<td align="right">110</td>
<td align="right">79</td>
</tr>
<tr class="even">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Scophthalmus rhombus</em></td>
<td align="left">Brill</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">14</td>
<td align="right">30</td>
<td align="right">30</td>
<td align="right">18</td>
</tr>
<tr class="odd">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Zeugopterus punctatus</em></td>
<td align="left">Topknot</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">6</td>
</tr>
<tr class="even">
<td align="left">Scophthalmidae</td>
<td align="left"><em>Zeugopterus regius</em></td>
<td align="left">Eckstr<f6>m's topknot</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Soleidae</td>
<td align="left"><em>Buglossidium luteum</em></td>
<td align="left">Solenette</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">37</td>
<td align="right">37</td>
<td align="right">8</td>
</tr>
<tr class="even">
<td align="left">Soleidae</td>
<td align="left"><em>Microchirus variegatus</em></td>
<td align="left">Thickback sole</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">35</td>
<td align="right">35</td>
<td align="right">11</td>
</tr>
<tr class="odd">
<td align="left">Soleidae</td>
<td align="left"><em>Pegusa lascaris</em></td>
<td align="left">Sand sole</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">59</td>
</tr>
<tr class="even">
<td align="left">Soleidae</td>
<td align="left"><em>Solea solea</em></td>
<td align="left">Common sole</td>
<td align="right">7</td>
<td align="right">11</td>
<td align="right">7</td>
<td align="right">2</td>
<td align="right">24</td>
<td align="right">151</td>
<td align="right">151</td>
<td align="right">252</td>
</tr>
<tr class="odd">
<td align="left">Salmonidae</td>
<td align="left"><em>Oncorhynchus mykiss</em></td>
<td align="left">Rainbow trout</td>
<td align="right">21</td>
<td align="right">27</td>
<td align="right">21</td>
<td align="right">22</td>
<td align="right">43</td>
<td align="right">454</td>
<td align="right">454</td>
<td align="right">76</td>
</tr>
<tr class="even">
<td align="left">Salmonidae</td>
<td align="left"><em>Salmo salar</em></td>
<td align="left">Atlantic salmon</td>
<td align="right">13</td>
<td align="right">31</td>
<td align="right">13</td>
<td align="right">14</td>
<td align="right">25</td>
<td align="right">371</td>
<td align="right">373</td>
<td align="right">34</td>
</tr>
<tr class="odd">
<td align="left">Salmonidae</td>
<td align="left"><em>Salmo trutta</em></td>
<td align="left">Sea trout</td>
<td align="right">24</td>
<td align="right">36</td>
<td align="right">24</td>
<td align="right">19</td>
<td align="right">30</td>
<td align="right">289</td>
<td align="right">289</td>
<td align="right">363</td>
</tr>
<tr class="even">
<td align="left">Salmonidae</td>
<td align="left"><em>Thymallus thymallus</em></td>
<td align="left">Grayling</td>
<td align="right">27</td>
<td align="right">28</td>
<td align="right">27</td>
<td align="right">26</td>
<td align="right">31</td>
<td align="right">96</td>
<td align="right">96</td>
<td align="right">34</td>
</tr>
<tr class="odd">
<td align="left">Agonidae</td>
<td align="left"><em>Agonus cataphractus</em></td>
<td align="left">Hooknose</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">24</td>
<td align="right">24</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td align="left">Cottidae</td>
<td align="left"><em>Cottus gobio</em></td>
<td align="left">Bullhead</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">86</td>
<td align="right">86</td>
<td align="right">7</td>
</tr>
<tr class="odd">
<td align="left">Cottidae</td>
<td align="left"><em>Micrenophrys lilljeborgii</em></td>
<td align="left">Norway bullhead</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Cottidae</td>
<td align="left"><em>Myoxocephalus scorpius</em></td>
<td align="left">Shorthorn sculpin</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">9</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">103</td>
<td align="right">103</td>
<td align="right">111</td>
</tr>
<tr class="odd">
<td align="left">Cottidae</td>
<td align="left"><em>Taurulus bubalis</em></td>
<td align="left">Longspined bullhead</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left">Cyclopteridae</td>
<td align="left"><em>Cyclopterus lumpus</em></td>
<td align="left">Lumpfish</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">71</td>
<td align="right">71</td>
<td align="right">7</td>
</tr>
<tr class="odd">
<td align="left">Liparidae</td>
<td align="left"><em>Liparis liparis</em></td>
<td align="left">Striped seasnail</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">15</td>
<td align="right">15</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left">Liparidae</td>
<td align="left"><em>Liparis montagui</em></td>
<td align="left">Montagus seasnail</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td align="left">Triglidae</td>
<td align="left"><em>Chelidonichthys cuculus</em></td>
<td align="left">Red gurnard</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">51</td>
<td align="right">51</td>
<td align="right">8</td>
</tr>
<tr class="even">
<td align="left">Triglidae</td>
<td align="left"><em>Chelidonichthys lastoviza</em></td>
<td align="left">Streaked gurnard</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">8</td>
<td align="right">35</td>
<td align="right">35</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left">Triglidae</td>
<td align="left"><em>Chelidonichthys lucerna</em></td>
<td align="left">Tub gurnard</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">17</td>
<td align="right">71</td>
<td align="right">71</td>
<td align="right">11</td>
</tr>
<tr class="even">
<td align="left">Triglidae</td>
<td align="left"><em>Eutrigla gurnardus</em></td>
<td align="left">Grey gurnard</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">2</td>
<td align="right">7</td>
<td align="right">30</td>
<td align="right">30</td>
<td align="right">9</td>
</tr>
<tr class="odd">
<td align="left">Siluridae</td>
<td align="left"><em>Silurus glanis</em></td>
<td align="left">Wels catfish</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">9</td>
<td align="right">76</td>
<td align="right">76</td>
<td align="right">6</td>
</tr>
<tr class="even">
<td align="left">Syngnathidae</td>
<td align="left"><em>Entelurus aequoreus</em></td>
<td align="left">Snake pipefish</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">180</td>
</tr>
<tr class="odd">
<td align="left">Syngnathidae</td>
<td align="left"><em>Hippocampus guttulatus</em></td>
<td align="left">Long-snouted seahorse</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left">Syngnathidae</td>
<td align="left"><em>Hippocampus hippocampus</em></td>
<td align="left">Short snouted seahorse</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">6</td>
<td align="right">2</td>
<td align="right">23</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">17</td>
</tr>
<tr class="odd">
<td align="left">Syngnathidae</td>
<td align="left"><em>Nerophis lumbriciformis</em></td>
<td align="left">Worm pipefish</td>
<td align="right">124</td>
<td align="right">119</td>
<td align="right">124</td>
<td align="right">1</td>
<td align="right">120</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">119</td>
</tr>
<tr class="even">
<td align="left">Syngnathidae</td>
<td align="left"><em>Nerophis ophidion</em></td>
<td align="left">Straightnose pipefish</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Syngnathidae</td>
<td align="left"><em>Syngnathus acus</em></td>
<td align="left">Greater pipefish</td>
<td align="right">5</td>
<td align="right">15</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">8</td>
<td align="right">43</td>
<td align="right">43</td>
<td align="right">17</td>
</tr>
<tr class="even">
<td align="left">Syngnathidae</td>
<td align="left"><em>Syngnathus rostellatus</em></td>
<td align="left">Nilsson's pipefish</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">14</td>
<td align="right">14</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td align="left">Syngnathidae</td>
<td align="left"><em>Syngnathus typhle</em></td>
<td align="left">Broadnosed pipefish</td>
<td align="right">4</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">6</td>
<td align="right">50</td>
<td align="right">50</td>
<td align="right">38</td>
</tr>
<tr class="even">
<td align="left">Balistidae</td>
<td align="left"><em>Balistes capriscus</em></td>
<td align="left">Grey triggerfish</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">59</td>
<td align="right">59</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left">Molidae</td>
<td align="left"><em>Mola mola</em></td>
<td align="left">Ocean sunfish</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">4</td>
<td align="right">10</td>
<td align="right">19</td>
<td align="right">19</td>
<td align="right">12</td>
</tr>
<tr class="even">
<td align="left">Zeidae</td>
<td align="left"><em>Zeus faber</em></td>
<td align="left">John dory</td>
<td align="right">10</td>
<td align="right">8</td>
<td align="right">10</td>
<td align="right">6</td>
<td align="right">25</td>
<td align="right">117</td>
<td align="right">117</td>
<td align="right">19</td>
</tr>
<tr class="odd">
<td align="left">Petromyzontidae</td>
<td align="left"><em>Lampetra fluviatilis</em></td>
<td align="left">River lamprey</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">8</td>
<td align="right">50</td>
<td align="right">50</td>
<td align="right">45</td>
</tr>
<tr class="even">
<td align="left">Petromyzontidae</td>
<td align="left"><em>Lampetra planeri</em></td>
<td align="left">European brook lamprey</td>
<td align="right">2</td>
<td align="right">3</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">19</td>
<td align="right">72</td>
<td align="right">72</td>
<td align="right">123</td>
</tr>
<tr class="odd">
<td align="left">Petromyzontidae</td>
<td align="left"><em>Petromyzon marinus</em></td>
<td align="left">Sea lamprey</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">58</td>
<td align="right">58</td>
<td align="right">10</td>
</tr>
<tr class="even">
<td align="left">Carcharhinidae</td>
<td align="left"><em>Prionace glauca</em></td>
<td align="left">Blue shark</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">6</td>
<td align="right">7</td>
<td align="right">9</td>
<td align="right">862</td>
<td align="right">866</td>
<td align="right">310</td>
</tr>
<tr class="odd">
<td align="left">Scyliorhinidae</td>
<td align="left"><em>Scyliorhinus canicula</em></td>
<td align="left">Lesser spotted dogfish</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">585</td>
<td align="right">585</td>
<td align="right">84</td>
</tr>
<tr class="even">
<td align="left">Scyliorhinidae</td>
<td align="left"><em>Scyliorhinus stellaris</em></td>
<td align="left">Nursehound</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">22</td>
<td align="right">22</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left">Triakidae</td>
<td align="left"><em>Galeorhinus galeus</em></td>
<td align="left">Tope shark</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">31</td>
<td align="right">31</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Triakidae</td>
<td align="left"><em>Mustelus asterias</em></td>
<td align="left">Starry smooth-hound</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">45</td>
<td align="right">45</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td align="left">Alopiidae</td>
<td align="left"><em>Alopias vulpinus</em></td>
<td align="left">Thresher</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">6</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">53</td>
<td align="right">53</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left">Cetorhinidae</td>
<td align="left"><em>Cetorhinus maximus</em></td>
<td align="left">Basking shark</td>
<td align="right">25</td>
<td align="right">24</td>
<td align="right">25</td>
<td align="right">24</td>
<td align="right">25</td>
<td align="right">82</td>
<td align="right">82</td>
<td align="right">25</td>
</tr>
<tr class="odd">
<td align="left">Lamnidae</td>
<td align="left"><em>Lamna nasus</em></td>
<td align="left">Porbeagle</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">118</td>
<td align="right">118</td>
<td align="right">3</td>
</tr>
<tr class="even">
<td align="left">Dasyatidae</td>
<td align="left"><em>Dasyatis pastinaca</em></td>
<td align="left">Common stingray</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">5</td>
<td align="right">41</td>
<td align="right">41</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left">Rajidae</td>
<td align="left"><em>Amblyraja radiata</em></td>
<td align="left">Starry ray</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">206</td>
<td align="right">206</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td align="left">Rajidae</td>
<td align="left"><em>Leucoraja naevus</em></td>
<td align="left">Cuckoo ray</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">43</td>
<td align="right">43</td>
<td align="right">6</td>
</tr>
<tr class="odd">
<td align="left">Rajidae</td>
<td align="left"><em>Raja brachyura</em></td>
<td align="left">Blonde ray</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">8</td>
<td align="right">45</td>
<td align="right">45</td>
<td align="right">9</td>
</tr>
<tr class="even">
<td align="left">Rajidae</td>
<td align="left"><em>Raja clavata</em></td>
<td align="left">Thornback ray</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">3</td>
<td align="right">12</td>
<td align="right">194</td>
<td align="right">194</td>
<td align="right">16</td>
</tr>
<tr class="odd">
<td align="left">Rajidae</td>
<td align="left"><em>Raja microocellata</em></td>
<td align="left">Small-eyed ray</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">3</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left">Rajidae</td>
<td align="left"><em>Raja montagui</em></td>
<td align="left">Spotted ray</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">7</td>
<td align="right">95</td>
<td align="right">95</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left">Rajidae</td>
<td align="left"><em>Raja undulata</em></td>
<td align="left">Undulate ray</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">3</td>
<td align="right">5</td>
<td align="right">12</td>
<td align="right">12</td>
<td align="right">3</td>
</tr>
<tr class="even">
<td align="left">Squalidae</td>
<td align="left"><em>Squalus acanthias</em></td>
<td align="left">Picked dogfish</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">4</td>
<td align="right">5</td>
<td align="right">16</td>
<td align="right">279</td>
<td align="right">279</td>
<td align="right">16</td>
</tr>
</tbody>
</table>

**Table 3. Numbers of new sequences for latest reference library version compared to previous. Current version is GenBank v247 (13 Jan 2022); previous version is GenBank v245 (21 Aug 2021).**

<table>
<colgroup>
<col width="21%" />
<col width="8%" />
<col width="8%" />
<col width="10%" />
<col width="11%" />
<col width="8%" />
<col width="10%" />
<col width="8%" />
<col width="11%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Scientific Name</th>
<th align="right">12S (Miya)</th>
<th align="right">12S (Riaz)</th>
<th align="right">12S (Taberlet)</th>
<th align="right">12S (Valentini)</th>
<th align="right">16S (Berry)</th>
<th align="right">COI (Lerayxt)</th>
<th align="right">COI (Ward)</th>
<th align="right">CYTB (Minamoto)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><em>Abramis brama</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Acipenser ruthenus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Agonus cataphractus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Alburnus alburnus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Alepisaurus ferox</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Ameiurus nebulosus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Ammodytes marinus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Ammodytes tobianus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Anguilla anguilla</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Anoplogaster cornuta</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Antimora rostrata</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Apristurus profundorum</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">16</td>
<td align="right">16</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Argyropelecus hemigymnus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Atherina boyeri</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Auxis rochei</em></td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td align="left"><em>Auxis thazard</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">11</td>
<td align="right">11</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Balistes capriscus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Bathypterois dubius</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Belone belone</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Blennius ocellaris</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Borostomias antarcticus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Brosme brosme</em></td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left"><em>Buglossidium luteum</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Callionymus reticulatus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Carassius auratus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">13</td>
<td align="right">13</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Centrolophus niger</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Centrophorus squamosus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Ceratias holboelli</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Ceratoscopelus warmingii</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">24</td>
<td align="right">24</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Chauliodus sloani</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Chelidonichthys cuculus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Chelidonichthys lucerna</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Chelon auratus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Chelon ramada</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Cobitis taenia</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Conger conger</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Coptodon zillii</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">10</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Coregonus albula</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Coris julis</em></td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td align="left"><em>Cyclothone microdon</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Cyprinus carpio</em></td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">29</td>
<td align="right">29</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Dactylopterus volitans</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Dasyatis pastinaca</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Dipturus batis</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Dipturus nidarosiensis</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Electrona risso</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Enchelyopus cimbrius</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Engraulis encrasicolus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">13</td>
<td align="right">13</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Entelurus aequoreus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Esox lucius</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td align="left"><em>Euthynnus alletteratus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">55</td>
<td align="right">55</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Eutrigla gurnardus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Evermannella balbo</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Gadus morhua</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Gaidropsarus mediterraneus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Gaidropsarus vulgaris</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Galeorhinus galeus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Gobio gobio</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">16</td>
<td align="right">16</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Gobius cobitis</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Gobiusculus flavescens</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left"><em>Gobius niger</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Gobius paganellus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Grammicolepis brachiusculus</em></td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Gymnammodytes semisquamatus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Helicolenus dactylopterus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Hippocampus guttulatus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Hyperoplus lanceolatus</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Ictalurus punctatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Isurus oxyrinchus</em></td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">3</td>
</tr>
<tr class="even">
<td align="left"><em>Katsuwonus pelamis</em></td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">44</td>
<td align="right">44</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td align="left"><em>Lamna nasus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lampetra planeri</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Lampris guttatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lepidocybium flavobrunneum</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Lepidorhombus whiffiagonis</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lepomis gibbosus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Leucaspius delineatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Leuciscus idus</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Limanda limanda</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Liparis liparis</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Liparis montagui</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lophius budegassa</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">19</td>
<td align="right">19</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Lota lota</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">13</td>
<td align="right">13</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lycenchelys sarsii</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Lycodes terraenovae</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Lycodes vahlii</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Malacosteus niger</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Melanonus zugmayeri</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Merlangius merlangus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">12</td>
<td align="right">12</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Merluccius merluccius</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">41</td>
<td align="right">41</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Micromesistius poutassou</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Microstomus kitt</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Misgurnus fossilis</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Mobula mobular</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Mola mola</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Molva dypterygia</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Molva molva</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Mugil cephalus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">13</td>
<td align="right">13</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Mullus barbatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">30</td>
<td align="right">30</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Mullus surmuletus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">33</td>
<td align="right">33</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Mustelus mustelus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Myctophum nitidulum</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Myliobatis aquila</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Myoxocephalus scorpius</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Naucrates ductor</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Neocyttus helgae</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Neogobius melanostomus</em></td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">6</td>
<td align="right">6</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Nerophis lumbriciformis</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Nerophis ophidion</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Nessorhamphus ingolfianus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Notolychnus valdiviae</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Oncorhynchus gorbuscha</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Oncorhynchus mykiss</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Oreochromis niloticus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">130</td>
<td align="right">130</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Paraliparis hystrix</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">4</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Perca fluviatilis</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">11</td>
</tr>
<tr class="odd">
<td align="left"><em>Pholis gunnellus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Phoxinus phoxinus</em></td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Phrynorhombus norvegicus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Phycis blennoides</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Platichthys flesus</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Pleuronectes platessa</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Poecilia reticulata</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Polyprion americanus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Pomatoschistus microps</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Pomatoschistus minutus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td align="left"><em>Pomatoschistus pictus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Prionace glauca</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">29</td>
<td align="right">29</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Psenes maculatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Pseudorasbora parva</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">28</td>
<td align="right">28</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Pseudoscopelus altipinnis</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Pterycombus brama</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Raja microocellata</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Raja miraletus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Raja undulata</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Raniceps raninus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Regalecus glesne</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Reinhardtius hippoglossoides</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Rhodeus amarus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Salmo salar</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Salvelinus alpinus</em></td>
<td align="right">96</td>
<td align="right">94</td>
<td align="right">96</td>
<td align="right">94</td>
<td align="right">94</td>
<td align="right">94</td>
<td align="right">94</td>
<td align="right">94</td>
</tr>
<tr class="even">
<td align="left"><em>Sander lucioperca</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">10</td>
<td align="right">10</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Sarda sarda</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">5</td>
<td align="right">5</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Sardina pilchardus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Schedophilus medusophagus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Scomber scombrus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Scopelogadus beanii</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Scopelosaurus lepidus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Scophthalmus maximus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Scophthalmus rhombus</em></td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">2</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Scorpaena porcus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Sebastes norvegicus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Seriola dumerili</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">2</td>
<td align="right">2</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Silurus glanis</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Solea solea</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Sphoeroides pachygaster</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Sphyrna zygaena</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">11</td>
<td align="right">11</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Sprattus sprattus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">24</td>
<td align="right">24</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Squatina squatina</em></td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
<td align="right">7</td>
</tr>
<tr class="even">
<td align="left"><em>Stomias boa</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Symphodus melops</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Thunnus albacares</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">15</td>
<td align="right">15</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td align="left"><em>Torpedo marmorata</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td align="left"><em>Trachinotus ovatus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">8</td>
<td align="right">8</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Trachinus draco</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Trachipterus arcticus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Trachurus trachurus</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">9</td>
<td align="right">9</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Trichiurus lepturus</em></td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">106</td>
<td align="right">72</td>
<td align="right">72</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Umbrina cirrosa</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Xiphias gladius</em></td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">58</td>
<td align="right">58</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td align="left"><em>Zeus faber</em></td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">3</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left"><em>Zoarces viviparus</em></td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">0</td>
</tr>
</tbody>
</table>
