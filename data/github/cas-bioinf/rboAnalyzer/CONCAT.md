# rboAnalyzer
A tool for analyzing BLAST search output for RNA sequences.

## Short description
rboAnalyzer is a tool complementing the BLAST algorithm when searching for a query
 sequence that is RNA with a secondary structure (which does not have to be known).

The high-scoring pairs (HSPs) in BLAST output are often incomplete
 (ie. the alignment in HSP does not cover the whole query sequence).
This is a major drawback when trying to characterize the potential ncRNA
 indicated by the HSP.

Therefore, rboAnalyzer tries to find full-length RNA sequences
 from the incomplete HSPs from the BLAST output and predict their secondary structures
 with one or more methods.
Score for similarity (as proxy to homology) between the estimated full-length sequence
 and query sequence is also computed.
The BLAST output is combined with computed data and presented in form of an interactive HTML page.

To achieve this, rboAnalyzer takes as input:
- the query sequence
- the BLAST output
- the BLAST database containing sequences within the output

The paper for this tool is available [here](https://doi.org/10.3389/fgene.2020.00675).

You can also try out the [rboanalyzer webserver](http://rboanalyzer.elixir-czech.cz) implementing interactive HSP analysis with selected secondary structure prediction methods.

## Installation
This tool was tested on 64-bit linux (Ubuntu 14, 18 and Centos 7). 

### Install via Conda
 The most convenient way to install this pipeline is to use Conda package manager. 
 
 If you don't have Conda for `Python3`, install it. We recommend the [miniconda3](https://conda.io/en/latest/miniconda.html).

The `rboanalyzer` package is available from `schwarz.marek` channel. You also need to include the `bioconda` and `conda-forge` channels when installing the `rboanalyzer`.

 
__Installation to new virtual environment__

The `rboAnalyzer` and its dependencies will be available only in the shell session for which the virtual environment was activated. You will need to `activate` the virtual environment for each session.

Following commands will install the rboAnalyzer into new conda environment named "rbo". If you don't wish to install rboAnalyzer to new environment, omit commands `1)` and `2)`.

```shell
# 1) update conda
conda update conda

# 2) create virtual environment
#conda create -n YOUR_VIRTUAL_ENV_NAME
conda create -n rbo

# 3) activate it
# conda activate YOUR_VIRTUAL_ENV_NAME
conda activate rbo

# 4) run installation
# ommit "-c conda-forge" and/or "-c bioconda" if those channels are in your .condarc file
conda install -c conda-forge -c bioconda -c schwarz.marek rboanalyzer
```

### Install from source

 __Prerequisites__
- python >= 3.4, <3.8 [link](https://www.python.org/downloads/)
  Verify that you have latest compatible `pip3`. If you try to use EOL python, please see the `setup.py` script on compatible `pip3` and `setuptools` versions and use them. 

- ncbi-blast+ >= 2.8.1, <2.10 [link](http://ftp.ncbi.nih.gov/blast/executables/blast+/2.9.0/)
  (The pipeline can use blast from version 2.6.0, however this version is not compatible with blast dbv5)
- locarna >= 1.9.2, <2 [link](https://github.com/s-will/LocARNA/releases/tag/v1.9.2.2)
- infernal >= 1.1, <1.2 [link](http://eddylab.org/infernal/)
- clustalo >= 1.2.4, <2 [link](http://www.clustal.org/omega/)
- muscle >= 3.8.31, <4 [link](https://www.drive5.com/muscle/downloads.htm)

For prediction:
- viennarna >=2.3.5, <3 (with refold.pl in PATH) [link](https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_3_x/ViennaRNA-2.3.5.tar.gz)
  Don't forget to add `refold.pl` to your `PATH`. The `refold.pl` script is located in the `ViennaRNA-[version]/src/Utils/`.
- centroid_homfold >= 0.0.15, <0.1 [link](https://github.com/satoken/centroid-rna-package/releases/tag/v0.0.15)
- RNAstructure >= 6.0, <7 (TurboFold - Text (Command Line) Interfaces ) [link](https://rna.urmc.rochester.edu/RNAstructure.html)
  Don't forget to set the `DATAPATH` environment variable [link](http://rna.urmc.rochester.edu/Text/Thermodynamics.html).

Optional (some prediction methods are not available without):
- UNAFold >= 3.8, <4 [link](http://www.unafold.org/)

Download this repository (or release), unpack it if needed. Go to directory with the source code for the rboAnalyzer and run

```shell
python3 setup.py install --user
```
Note that the `--user` switch puts the executables in `$HOME/.local/bin` and you may need to add it to `PATH`. 

The rboAnalyzer executable should be created.
To test it, restart terminal (close and open new) and run

```shell
rboAnalyzer --version
```
which should return the version number.

## Preparation

<a name="rfamdownload" id="rfamdownload"></a>

### Obtain Rfam database
For correct function the rboAnalyzer needs a copy of Rfam database.

There are 2 ways:
1. Run rboAnalyzer with `--download_rfam` flag.
    ```shell
    rboAnalyzer --download_rfam
    ```
This will download Rfam covariance models to default directory
(`[INSTALL_LOCATION]/rna_blast_analyze/3rd_party_source/rfam`).

2. Alternatively download the `Rfam.cm.gz` file from
[Rfam CURRENT](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT),
unpack it and add the path to directory containing `Rfam.cm` file to your `config.txt` file.
Note that running rboAnalyzer with `--download_rfam` will overwrite this manually installed file.

### Installing UNAFold (optional)
Prediction methods using suboptimal structures need UNAFold software to work.
 It is available here [http://www.unafold.org](http://www.unafold.org).
 Follow installation instructions. The pipeline uses the `hybrid-ss-min` program.
 Either add it to PATH or add path to directory containing the executable to `config.txt` file.

### Shell autocomplete (optional)
The rboAnalyzer is equipped with argument completion for bash shell.
To enable this feature the `argcomplete` package needs to be installed (version >1.6, <2), and script registered with your shell.
It is installed by default in conda package. To install it manually run
```shell
pip3 install --user "argcomplete>1.6, <2a0"
```
Note that the `--user` switch puts the executables in `$HOME/.local/bin` and you may need to add it to `PATH`.

To register the script to your shell (more info [here](https://pypi.org/project/argcomplete/)).
```shell
register-python-argcomplete rboAnalyzer >> ~/.bashrc
```


### BLAST database
The rboAnalyzer needs to get relevant 5' and 3' regions of subject sequence of HSPs, for this we use the BLAST database used in the BLAST search.
For each analysis you need to provide the nucleotide BLAST database containing whole sequences (complete genomes, etc.) for the sequence ids present in the BLAST output.

The procedure on how to get the BLAST databases for the examples is described in the Example section. For the general information, please see the BLAST databases [section](#blastdatabase).

## Basic Usage
### Help
```shell
rboAnalyzer -h
```
<details>
<summary>help output</summary>

```
usage: rboAnalyzer [-h] -in PATH -q PATH (-db path | --entrez ENTREZ)
                   [--db_type {blastdb,fasta,gb,server}]
                   [--b_type {guess,xml,plain}] [--blast_regexp BLAST_REGEXP]
                   [--mode {simple,locarna,meta}] [--turbo_fast_preset]
                   [--centroid_fast_preset] [--html PATH] [--threads N]
                   [--csv CSV] [--json PATH] [--cm_file CM_file | --use_rfam]
                   [--download_rfam] [--version] [--config_file PATH]
                   [-pm [prediction_method_name [prediction_method_name ...]]]
                   [--pm_param_file PATH] [--logfile logfile]
                   [--subseq_window_locarna SUBSEQ_WINDOW_LOCARNA]
                   [--locarna_anchor_length LOCARNA_ANCHOR_LENGTH]
                   [--filter_by_eval FILTER_BY_EVAL | --filter_by_bitscore FILTER_BY_BITSCORE]
                   [-v] [--skip_missing] [--show_HSP]

rboAnalyzer - tool for analyzing BLAST search output for RNA sequences

optional arguments:
  -h, --help            show this help message and exit
  -db path, --blast_db path
                        Provide path to blast database, that is the complete
                        path with blast db name without any extension (*.nin,
                        nsd, nog, nsi, nhr, nsq, nal).
  --entrez ENTREZ       EMAIL - Indicate that you want to use NCBI Entrez
                        service to download required regions of sequences at
                        runtime. To comply with NCBI service rules you are
                        required to provide valid email address at which the
                        NCBI staff could contact you if they need to.

INPUT:
  -in PATH, --blast_in PATH
                        BLAST output file with hits to analyze.
  -q PATH, --blast_query PATH
                        The Blast query fasta file.

OUTPUT:
  --html PATH           Output html file with secondary structure pictures and
                        other useful stuff.
  --csv CSV             Output in csv table, infered sequence and structure
                        present.
  --json PATH           Dump all stored data to JSON (developer only - it is
                        possible to convert to all other output formats).

PARAMETERS:
  --mode {simple,locarna,meta}
                        Choose mode of hit elongation: simple (extend by
                        unaligned parts of query) locarna (run locarna
                        algorithm - uses secondary structure for better
                        alignment) meta (uses both methods and chooses the
                        alignment which has better RSEARCH score).
  --turbo_fast_preset   Act's as parameter preset for Turbo-fast setting the
                        max_seqs_in_prediction to 2. This means that only the
                        query sequence is used as homologous sequence. It is
                        useful if analyzing very distant BLAST HITs.
  --centroid_fast_preset
                        Parameter preset for centroid-fast. Set's the
                        max_seqs_in_prediction to 1. This means that only the
                        query sequence is used as homologous sequence for
                        prediction. It is useful if analyzing very distant
                        BLAST HITs.
  --config_file PATH    Provide config file if tools and data are in non-
                        default paths.
  -pm [prediction_method_name [prediction_method_name ...]], --prediction_method [prediction_method_name [prediction_method_name ...]]
                        Prediction method to use. Multiple prediction methods
                        are allowed. Possible values: C-A-U-r-Rc centroid
                        rnafold M-A-sub M-A-r-Rc fq-sub Turbo-fast TurboFold
                        centroid-fast C-A-sub C-A-r-Rc M-A-U-r-Rc rfam-sub
                        rfam-Rc rfam-centroid
  --pm_param_file PATH  Path to file with parameters for prediction methods in
                        JSON. Prediction methods not declared within provided
                        file are used with default values. File is in json
                        format. Default values (also example how to provide
                        parameters) are stored in '[install location]/rna_blas
                        t_analyze/BR_core/prediction_parameters.json'
  --subseq_window_locarna SUBSEQ_WINDOW_LOCARNA
                        N of nucleotides to add to expected start/end of
                        sequence before realignement. The unaligned
                        nucleotides are not included in reported sequence.
  --locarna_anchor_length LOCARNA_ANCHOR_LENGTH
                        Minimal number of adjacent matching bases in BLAST hit
                        to create an anchor for Locarna.

MISC:
  --db_type {blastdb,fasta,gb,server}
                        Type of a database provided. If 'fasta' or 'gb' then
                        --blast_db must be directory containing files with
                        names in accession.version format. Example '/home/my-
                        best-db/' with files like 'CP000001.1'.
  --b_type {guess,xml,plain}
  --blast_regexp BLAST_REGEXP
                        Provide python valid regular expression which capture
                        the index key to blastdb (usualy the accession.version
                        number).
  --threads N           Number of threads to use (default = N of logical cores
                        detected).
  --cm_file CM_file     Provided covariance model will be used for homology
                        inference instead of RSEARCH model.
  --use_rfam            Search in rfam database for covariance model to infer
                        homology with instead of RSEARCH model.
  --download_rfam       Retrieve RFAM covariance models database. Will
                        download only if new version avalible.
  --version             show program's version number and exit
  --logfile logfile     Path to where logfile should be written.
  --filter_by_eval FILTER_BY_EVAL
                        Filter the input blast by E-value. Only hits following
                        the rule will be kept. Example ">10e-10" will keep
                        only hits with eval greater then 10e-10. Interval can
                        be specified with "," e.g. ">10e-100, <10e-1". The
                        homologous sequences used with certain prediction
                        methods are taken from all hits (regardless of the
                        filtering).
  --filter_by_bitscore FILTER_BY_BITSCORE
                        Filter the input blast by bit score. Only hits
                        following the rule will be kept. Example "<20" will
                        keep only hits with bit score less then 20. Interval
                        can be specified with "," e.g. ">30, <45". The
                        homologous sequences used with certain prediction
                        methods are taken from all hits (regardless of the
                        filtering).
  -v, --verbose         output verbosity -> most detailed -vv (lot of output)
  --skip_missing        If given, the missing records in given blast database
                        will be skipped. This may alter the results of
                        bit_score computation (for homology prediction) and
                        secondary structure prediction for several methods.
  --show_HSP            Show HSP marker in NCBI sequence viewer.
```
</details>

### Usage
```shell
rboAnalyzer -in BLAST_OUTPUT.xml -db USED_DATABASE_PATH -q BLAST_QUERY.fasta
```
Note that only __one__ query sequence in the query FASTA file is expected.

## Example
Examples are provided in example directory.

To try examples you will need to:

  1. Install the rboAnalyzer and download Rfam. [How to?](#rfamdownload).

  2. Obtain a copy of `example` directory. If you've cloned or downloaded the program you should already have it.
   Otherwise it is [here](https://github.com/cas-bioinf/rna_blast_analyze/tree/master/example).

  3. cd to `example` directory.

#### Example 1:

Analyzing subset of NCBI blast HITs for [6S RNA](https://doi.org/10.1038%2F229147a0).

1) Now you need to obtain a copy of the BLAST database with all
 accessions which are in the BLAST output.
 As above, you can either get the NCBI nt database or download only the sequences (genomes) from the BLAST output.
 Here we describe the variant with downloading only the necessary sequences.
 This is done by using the `genomes_from_blast` by calling:
    ```shell
    genomes_from_blast -e YOUR_EMAIL_ADDRESS -in 6S_super_short.xml -o genomes.bdb
    ```
    The parameter `-e` `YOUR_EMAIL_ADDRESS` should be your valid email address on which NCBI staff could contact you
    if they need to. It is not logged by the tool.

    The parameter `-in` is path to file (in this example `6S_super_short.xml`) containing the BLAST output.

    The parameter `-o` is output file path. In this command the BLAST database with name `genomes.bdb` was created for you if everything was successful.
    You will need it in the next step.

    The intermediate file `genomes.bdb.fasta` was also created and contains all sequences added to the BLAST database.
    When another BLAST output is analyzed (and sequences are needed) then only those sequences not present in the intermediate FASTA file are downloaded (assuming same BLAST db name is used).

2) Now you can run the pipeline itself:
    ```shell
    rboAnalyzer -in 6S_super_short.xml -q 6S_query.fasta -db genomes.bdb --html 6S_out.html
    ```
3) The output is single html file. You can scroll through analyzed HSPs, show the genomic loci of the HSP and select data to export.

#### Example 2:
Analyzing possible remote homologs for [MS1 RNA](https://doi.org/10.1093%2Fnar%2Fgku793). This take about 10 minutes on average pc.

```shell
# update the genome database with new sequences
genomes_from_blast -e YOUR_EMAIL_ADDRESS -in MS1_BLAST_output.txt -o genomes.bdb

# run the rboAnalyzer
rboAnalyzer -in MS1_BLAST_output.txt -q MS1_query.fasta -db genomes.bdb --html MS1_out.html --prediction_method rnafold rfam-Rc Turbo-fast --turbo_fast_preset
```

The BLAST was run with database where Streptomycetaceae and Mycobacteriaceae families where excluded.
As the MS1 RNA is primarily known from Mycobacteriaceae family we can expect
  incomplete HITs and many false positives.

Also you can notice, that the BLAST output was in text format, the rboAnalyzer accepts BLAST output in plain text and xml.

The BLAST database is obtained with similar command as in previous example.
Since the output BLAST database file is the same as before the `genomes_from_blast`
 will check the intermediate fasta file (`genomes.bdb.fasta`) and will
 download only sequences which are not present.

We can expect that the BLAST output contain many false positive HSPs,
so we selected prediction methods from those, which do not rely on information in the BLAST output itself.
These are:
- rnafold
- rfam-Rc
- rfam-centroid
- Turbo-fast with "max_seqs_in_prediction" parameter set to 2 (`--turbo_fast_preset` flag)
- centroid-fast with "max_seqs_in_prediction" parameter set to 1 (`--centroid_fast_preset` flag)
- rfam-sub (If `UnaFold` is installed.)


### Solving issues:
- __BLAST txt output from WEB not recognized__    
  Reason: the NCBI changed the `txt` output for the WEB BLAST when they switched to new design. Our txt parser is compatible with commandline `txt` output (`-outfmt 0`) which was also the `txt` output for the WEB BLAST.   
  Solution: Download the `xml` output or switch to the "Traditional result page" and download `txt` there.
- __One or more records not found__   
  Reason: the blastdbcmd was not able to find sequence(s) with respective id(s) in provided database.
  This is due to inconsistency between the sequence accessions and the BLAST database.
  The inconsistency may rise from:
    1. __sequence is not in the database__   
      Solution: Provide correct BLAST database (update current or create new with `genomes_from_blast`).
    2. __capturing regexp does not capture the accession number__   
      Solution: Provide capturing regular expression (python 3 syntax) for capturing the sequence id from the fasta header (it must match the id to the BLAST database used)
    3. __the BLAST database was created without the `-parse_seqids` flag__   
      Solution: Create new database from the sequences used to create new one, this time with `-parse_seqids` flag.
    4. __inconsistent accession mapping__   
      Detailed cause: In certain NCBI's BLAST databases there are some sequences with accession numbers for which ENTREZ 
       return sequences with different accession numbers. (We observed this for several sequences with accession numbers starting with `GPS_`.) 
       The `genomes_from_blast` script depends on matching accession numbers and can't handle the inconsistency.   
      Solution: We recommend obtaining the database from NCBI which was used to generate the BLAST output. 
       Other option is to manually add the sequence in question to the blast database.  

  Another option is to call pipeline with `--skip_missing` flag.
  This will skip the missing sequences.

  Note that no HSP for the missing sequence will be included in pipeline output
  and some prediction methods may be influenced by the missing sequence.

- __The `genomes_from_blast` failed__   
  The `genomes_from_blast` script has build in handling of failed downloads,
  but by default it tries only 10 times. If you are on unstable connection
  you might get better results by setting the `--retry` to some larger number.
  Also check if NCBI ENTREZ services are functional.

- __ValueError: could not convert string to float__    
  (raised with traceback containing `infer_homology.py ... dtypes/cast.py)    
  Please check, that there is only one query sequence in the input FASTA file. 

<a name="blastdatabase" id="blastdatabase"></a>

### BLAST databases

#### BLAST on NCBI web

If you used the BLAST using the NCBI web service against one of preformatted databased, you can download the whole database or use a `genomes_from_blast` script to download only the sequences in your blast output.

1. downloading whole database (~50GB)   
The latest databases are provided here [NCBI LATEST](ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST).
Note that databases included in the BLAST database releases are not the latest ones.
This code snippet can be used to obtain and update the database:
    ```shell
    # cd to directory to which you want to download the database

    # for the "nt" database
    wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/nt*

    # for other databases provided by NCBI (insert database name without square brackets)
    wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/[database name]*
    ```

2. downloading only relevant sequences   
  If you do not wish to download whole blastdb you may use prepared script
 `genomes_from_blast`, which downloads only the needed sequences
  (those in the BLAST output) and build the blastdb from them.
  This command will download all needed genomes and create BLAST database for you.
    ```shell
    # The `YOUR_EMAIL_ADDRESS` is needed so the NCBI would contact you in case of misuse of their resources.

    genomes_from_blast -e YOUR_EMAIL_ADDRESS -in BLAST_OUT_FILE -o BLAST_DATABASE_OUTFILE_NAME
    ```

#### Custom BLAST database
If custom database was used for the BLAST search you need to ensure multiple things for the rboAnalyzer to find the sequences correctly.
1. custom database is nucleotide and it was created with `-parse_seqids` (this makes sequences retrievable by their ids).
2. provide regular expression capturing the sequence ids. By default the rboAnalyzer captures the Accession.Version as documented [here](https://www.ncbi.nlm.nih.gov/Sequin/acc.html).

## Learn more
Description of the pipeline logic is avialable in the [publication](https://doi.org/10.3389/fgene.2020.00675).

Details can also be found [here](docs/help.md).

This readme and documentation in `docs` are valid for latest release 0.1.4 and dev 0.1.5a1.

## References
- ViennaRNA: Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
 ViennaRNA Package 2.0. <https://doi.org/10.1186/1748-7188-6-26>.
 We include copy of refold.pl script with Viennarna license here for convinience. [website](https://www.tbi.univie.ac.at/RNA/)
- Locarna: Sebastian Will, Tejal Joshi, Ivo L. Hofacker, Peter F. Stadler, and Rolf Backofen.
LocARNA-P: Accurate boundary prediction and improved detection of structural RNAs
RNA, 18 no. 5, pp. 900-14, 2012. <https://doi.org/10.1261/rna.029041.111>, [website](http://rna.informatik.uni-freiburg.de/LocARNA/Input.jsp)
- NCBI BLAST: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008)
 "BLAST+: architecture and applications." BMC Bioinformatics 10:421. <https://doi.org/10.1186/1471-2105-10-421>, [website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- RSEARCH: Finding Homologs of Single Structured RNA Sequences.
 R. J. Klein, S. R. Eddy. BMC Bioinformatics, 4:44, 2003. <https://doi.org/10.1186/1471-2105-4-44>, [website](http://eddylab.org/software.html#rsearch)
- Infernal: Infernal 1.1: 100-fold Faster RNA Homology Searches.
 E. P. Nawrocki, S. R. Eddy. Bioinformatics, 29:2933-2935, 2013. <https://doi.org/10.1093/bioinformatics/btt509>, [website](http://eddylab.org/infernal/)
- Clustal omega: Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011).
 Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.
  Molecular Systems Biology 7:539. <https://doi.org/10.1038/msb.2011.75>, [website](http://www.clustal.org/omega/)
- MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput
 Nucleic Acids Res. 32(5):1792-1797. <https://doi.org/10.1093/nar/gkh340>, [website](https://www.drive5.com/muscle/)
- Centroid homfold: Michiaki Hamada, Koichiro Yamada, Kengo Sato, Martin C. Frith, Kiyoshi Asai;
 CentroidHomfold-LAST: accurate prediction of RNA secondary structure using automatically collected homologous sequences,
 Nucleic Acids Research, Volume 39, Issue suppl_2, 1 July 2011, Pages W100–W106.
  <https://doi.org/10.1093/nar/gkr290>, [github](https://github.com/satoken/centroid-rna-package/)
- TurboFold (RNAstructure): Tan, Z., Fu, Y., Sharma, G., & Mathews, D. H. (2017).
 TurboFold II: RNA structural alignment and secondary structure prediction informed by multiple homologs.
  Nucleic Acids Research. 45: 11570-11581. <https://doi.org/10.1093/nar/gkx815>, [website](http://rna.urmc.rochester.edu/RNAstructure.html)
- UNAFold: Markham N.R., Zuker M. (2008) UNAFold. In: Keith J.M. (eds) Bioinformatics.
 Methods in Molecular Biology™, vol 453. Humana Press. <https://doi.org/10.1007/978-1-60327-429-6_1>, [website](http://www.unafold.org)
- Biopython: Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics.
 Bioinformatics 2009 Jun 1; 25(11) 1422-3. <http://dx.doi.org/10.1093/bioinformatics/btp163>, [website](https://biopython.org/)
- NumPy [website](https://numpy.org)
- Pandas [website](https://pandas.pydata.org)
- Jinja2 [website](https://jinja.palletsprojects.com)
- Matplotlib [website](https://matplotlib.org)

## Citation
If you find this software useful please cite this [paper](https://doi.org/10.3389/fgene.2020.00675).

Schwarz, M., Vohradský, J., Modrák, M. and Pánek, J., 2020. rboAnalyzer: A Software to Improve Characterization of Non-coding RNAs From Sequence Database Search Output. Frontiers in genetics, 11, p.675.

The rboAnalyzer version used for the paper is `0.1.4`.

## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](docs/ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund - Project "ELIXIR-CZ: Budování kapacit" (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](docs/logolink_OP_VVV_hor_barva_eng.jpg)

This work was also supported by Czech Science Foundation GA15-00885S.k
# Run rboAnalyzer using Docker

Docker image containing the rboAnalyzer and essential data.

## Install Docker

Install compatible docker version from https://www.docker.com/.


## rboAnalyzer container structure
- The container is at https://hub.docker.com/repository/docker/schwarzmarek/rboanalyzer
- The working directory of the rboAnalyzer in the container is `/data`. You are expected to mount your data there.
- The Rfam database is in the container and it's version is 14.1.
- The container is ready to use - no other preparation should be necessary

## Running the examples inside the docker container
### Prepare for running the examples
- Pull rboAnalyzer image
    ```
    docker pull schwarzmarek/rboanalyzer:0.1.4
    ```

- Obtain the `examples` directory

    Obtain a copy of `examples` directory. (Download and unpack https://github.com/cas-bioinf/rboAnalyzer/releases/download/v0.1.0/example.zip)

- In the terminal, go into the unpacked `examples` directory and run docker command(s) from there.

### Run interactive docker session 
Inside the session you can run rboAnalyzer commands. Take care to point the output files to the mounted directory (otherwise they will not be saved).

The commnand below expects that you are inside the directory with the data you want to analyze (e.g. the `examples` directory) and does following:
- Runs the `/bin/bash` in the specified container (`schwarzmarek/rboanalyzer:0.1.4`) in the interactive mode (`-it`) thus providing interactive session
- Mounts the current directory (`source="$(pwd)"`) to `/data` directory (`target=/data`) inside the container. This will cause that any changes in the mounted directory inside the container will be propagated to the directory on host. (for other options see https://docs.docker.com/storage/)

```
docker container run -it --mount type=bind,source="$(pwd)",target=/data schwarzmarek/rboanalyzer:0.1.4 /bin/bash
```

Now you should see the content of the `examples` dictionary. If so, you can execute the commands from the [Example](../readme.md#Example) section.# Prediction methods in more detail

Structure of this document:

First parameters and their meaning is described.

Second, prediction methods are listed. Each prediction method is described
  (what and how the structures are predicted).
For each method, there is section called "parameters",
  which serves as an example for parameter definition in json format.
Parameters are specified for each method independently on other prediction methods.

## Parameters and their meaning

### Third party tools commandline argument specification
The pipeline offer a way to interact with most parts of tools used for
  secondary structure prediction.

These tools are used for specific tasks and the parameters controlling
  such behaviour are always set and are not exposed for user control.
  (io flags, mode flags, etc.)


#### RNAfold: "RNAFOLD PARAMETERS"
is string with commandline arguments for `RNAfold`
  (see RNAFold [documentation](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)).
  It must be specified with double quotes.

Default: No parameters specified.

#### alifold: "ALIFOLD PARAMETERS"
is string with commandline arguments for `RNAalifold`
  (see RNAalifold [documentation](https://www.tbi.univie.ac.at/RNA/RNAalifold.1.html)).
It must be specified with double quotes.

Default: No parameters specified.

#### cmscan: "CMSCAN PARAMETERS"
is string with commandline arguments for `cmscan`
  (see cmscan [documentation](http://eddylab.org/infernal/)).
  It must be specified with double quotes.

Default: No parameters specified.

#### cmalign: "CMALIGN PARAMETERS"
is string with commandline arguments for `cmalign`
  (see cmalign [documentation](http://eddylab.org/infernal/)).
  It must be specified with double quotes.

Default: No parameters specified.

#### mfold: \[P, W, X\]
specifies parameters for `hybrid-ss-min` program where P, W, and X are integers
  see UNAFold [documenation](http://www.unafold.org/).
  It must be specified without double quotes.

The hybrid-ss-min is additionally used with `--suffix=DAT --NA=RNA --noisolate` parameters.

#### clustalo: "CLUSTALO PARAMETERS"
is string of commandline arguments for `clustalo`
  (see clustal omega [documentation](http://www.clustal.org/omega/))

Default: No parameters specified.

#### clustalo_profile: "CLUSTALO_PROFILE"
is string of commandline argument for `clustalo`. Allows specification
  of different `clustalo` parameters for the profile alignment stage.
  (see clustal omega [documentation](http://www.clustal.org/omega/))

Default: `clustalo` called with `--profile1` only profile align compatible option can be specified.

#### centroid_homfold: "CENTROID_HOMFOLD PARAMETERS"
is string of commandline argument for `centroid_homfold`
  (see centroid rna package [documentation](https://github.com/satoken/centroid-rna-package/)).

By default `-g -1` is used. If `-g` or `-t` is specified by user then the `-g -1` is not used.

#### muscle: "MUSCLE PARAMETERS"
is a string argument for `muscle` aligner. For more details, see muscle
  [documentation](https://www.drive5.com/muscle/manual/).

By default `-seqtype rna -quiet -clwstrict` are used and cannot be changed.

### Selection of estimated full-length sequences to provide reference for secondary structure prediction (reference sequences)
Parameters used for selection of sequences similar to the query sequence.
The selection parameters influence how much each sequence can differ from covariance model,
  how much the sequences can be similar to each other and what is maximum
  accepted length difference between the estimated full-length sequence and the query sequence.

Selection of estimated full-length sequences:
All estimated full-length sequences are filtered based on cm bit-score,
  then the remaining sequences are filtered according to their length to
  be within specified length difference from query sequence.
  After that, the remaining sequences are filtered based on sequence similarity with each other
  (from pairwise sequence identity), so only those sequences with similarity
  less then defined similarity threshold are retained.

#### cmscore_percent: CMSCORE_PERCENT
Defines inclusion threshold in percent of bit-score value obtained by
  aligning query sequence to covariance model.
It serves as filter for not sufficiently related sequences.
The inclusion threshold is computed by `CMSCORE_PERCENT`&nbsp;*&nbsp;`QUERY_BITSCORE`&nbsp;/&nbsp;`100`,
  only estimated full-length sequences with CM alignment bit-score higher then inclusion
  threshold are considered for further computation.
`CMSCORE_PERCENT` ranges from 0 to 100, allowed values are integers.

The higher the threshold, the more conservative setting for trusted sequences
  (i.e. more similarity to cm is required).

#### query_max_len_diff: MAX_LEN_DIFF_to_QUERY
Defines a maximum length difference between the query sequence and the estimated full-length sequence.
This serves complementary to cmscore_percent to allow setting low cmscore_percent,
while preventing very short or very long estimated full-length sequences (for any reason)
 to be part of selected sequences set.

`MAX_LEN_DIFF_to_QUERY` ranges from 0 to 1, allowed values are floating point (with decimal dot).

The higher the threshold, the more difference is allowed.

#### pred_sim_threshold: PRED_SIM_THRESHOLD
Defines exclusion threshold in percent of sequence similarity.
This serves as a protection from populating set of selected sequences with
  too many too similar sequences (as it may skew alignment and other prediction methods).
`PRED_SIM_THRESHOLD` ranges from 0 to 100, allowed values are integers.

The query is included as first sequence and all sequences which are similar to it
  down to defined threshold are excluded. Then next sequence from the remaining sequences is
  analyzed in same manner.

The higher the similarity the more similar sequences are accepted.
  (i.e. `PRED_SIM_THRESHOLD = 100` means that only exact duplicate sequences are removed)

#### repred_unpaired_tr: REPRED_UNPAIRED_TR
Serves for prediction methods where the conserved unpaired positions in alignment are taken as constraints.
Defines how much MSA column must be conserved at to denote the position
  as single-strand constraint for `RNAfold`.
This parameter is always used together with conseq_conserved.

#### conseq_conserved: CONSEQ_CONSERVED
Serves for prediction method where the conserved unpaired positions in alignment are taken as constraints.
Defines how many bases in a row must be conserved to denote the position
  as single-strand constraint for `RNAfold`.
This parameter is always used together with repred_unpaired_tr.

#### max_seqs_in_prediction: MAX_SEQS_IN_PREDICTION
Used with TurboFold to set required (and maximum) number of sequences in
  prediction. Allowed values are integers >= 2. Beware of setting this too high.
  Then the prediction is __very__ memory and time expensive.


## rnafold
Predict secondary structure for each estimated full-length sequence with RNAFold.

parameters:
```
    {"rnafold" : {
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## rfam-Rc
Obtain the related covariance model (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), align estimated full-length sequences to it (`cmalign`),
  then take conserved bp as constrains for for folding with `RNAfold -C`.

parameters:
```
  {"rfam-Rc": {
      "RNAfold": "RNAFOLD PARAMETERS",
      "cmscan": "CMSCAN PARAMETERS",
      "cmalign": "CMALIGN PARAMETERS"
    }
  }
```

Non optional parameters:
- `RNAFold`: `-C` (constraints)
- `cmscan`: `-g` (global aligment of models to query sequence)
- `cmalign`: `--notrunc` (disable truncated hits detection - speeds up computation).

## rfam-sub
Obtain the related CM (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), extract reference secondary structure,
  run `hybrid-ss-min` (UNAFold) for suboptimal structures and select the most similar
  one to the one from CM.

parameters:
```
    {"rfam-sub" : {
        "mfold": [P, W, X]
        }
    }

```

## fq-sub
Predict structure of query with RNAFold and take it as reference structure,
  then predict suboptimal structures with `hybrid-ss-min` (UNAFold) and select
  structure most similar to the predicted structure of query (by `RNAdistance` score).

parameters:
```
    {"fq-sub" : {
        "RNAfold": "RNAFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## C-A-sub
clustalo - RNAalifold - hybrid-ss-min

Select reference estimated full-length sequences, align them with `clustalo` (Clustal Omega),
  predict consensus structure with `RNAalifold`.
For each estimated full-length sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure (by `RNAdistance` score).

parameters:
```
    {"C-A-sub" : {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## M-A-sub
clustalo - RNAalifold - hybrid-ss-min

Select reference estimated full-length sequences, align them with `muscle`,
  predict consensus structure with `RNAalifold`.
For each estimated full-length sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure (by `RNAdistance` score).

parameters:
```
    {"M-A-sub" : {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## C-A-r-Rc
clustalo - RNAalifold - refold.pl - rnafold-C

Select reference estimated full-length sequences, compute alignment with `clustalo` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences
  and add the consensus structure to all estimated full-length sequences from profile alignment,
  finally run `refold.pl` and `RNAFold` for the result.

parameters:
```
    {"C-A-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## C-A-U-r-Rc
clustalo - RNAalifold - unpaired conserved - refold.pl - rnafold-C

Select reference estimated full-length sequences, compute alignment with `clustalo` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences,
  and add the consensus structure to all estimated full-length sequences from profile alignment.
Then select conserved parts of alignment where the consensus secondary
  structure annotation is single-strand (i.e. no base-pairs) and use them as constrains for `RNAFold -c`.

parameters:
```
    {"C-A-U-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```



## M-A-r-Rc
muscle - RNAalifold - refold.pl - rnafold -C

Select reference estimated full-length sequences, compute alignment with `muscle` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences
  and add the consensus structure to all estimated full-length sequences from profile alignment,
  finally run `refold.pl` and `RNAFold` for the result.

parameters:
```
    {"M-A-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## M-A-U-r-Rc
muscle - RNAalifold - unpaired conserved - refold.pl - rnafold-C

Select reference estimated full-length sequences, compute alignment with `muscle` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences,
  and add the consensus structure to all estimated full-length sequences from profile alignment.
Then select conserved parts of alignment where the consensus secondary
  structure annotation is single-strand (i.e. no base-pairs) and use them as constrains for `RNAFold -c`.

parameters:
```
    {"M-A-U-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## centroid
Select reference estimated full-length sequences and pass them to `centroid_homfold` as
  homologous sequences, predict secondary structure for all estimated full-length sequences.
The `centroid_homfold` is by default called with `-g -1` parameter.
With this setting the `centroid_homfold` predicts multiple structures
 (see [cetroid_homfold documentation](https://github.com/satoken/centroid-rna-package/)),
 the structure with best score is chosen.

parameters:
```
    {"centroid": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS"
        }
    }
```

## centroid-fast
With centroid-fast the reference estimated full-length sequences are selected as follows:
 we take non-redundant estimated full-length sequences without ambiguous
 bases with cm_score > 0 within allowed length-to-query-length difference.
Then we take up to `N` reference sequences as homologous sequences for `centroid_homfold`.
 The sequences are added in order of the original BLAST HSPs, starting with query.
For this method there is also parameter preset avalible which sets `N` to 1
 meaning that only query will be used as homologous sequence in prediction.

parameters:
```
    {"centroid-fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS",
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        }
    }
```

## rfam-centroid
Use covariance model (CM) to generate set of reference sequences for `centroid_homfold`.
By default the highest scoring CM in Rfam is found by `cmscan`
 (or provided one is used `--cm_file` option).
Then requested number of sequences (`N_SEQS`) is generated from the covariance model with `cmemit`.
The sequences are used as homologous sequences for `centroid_homfold`.
The generated sequences can differ between runs. If repeatable behaviour
 is desired the `cmemit` can be seeded (see it's options).

parameters:
```
    {"rfam-centroid": {
        "n_seqs": N_SEQS,
        "cmscan": "CMSCAN PARAMETERS",
        "cmemit": "CMEMIT PARAMETERS",
        }
    }
```

## Turbo-fast
With Turbo-fast the reference estimated full-length sequences are selected as follows:
we take non-redundant estimated full-length sequences without ambiguous
 bases with cm_score > 0 within allowed length-to-query-length difference.
For each estimated full-length sequence, we make non-redundant group of sequences
 consisting of the estimated full-length sequence for which we want to predict secondary structure
 and up to `N`-1 reference sequences.
The sequences are added in order of the original BLAST HSPs, starting with query.
That means that if `N` is 2 and query does not contain ambiguous bases,
 then each estimated full-length sequence secondary structure is computed with query sequence as a reference.
This setting is also available as commandline argument with `--turbo_fast_preset` flag and
 will override the `"max_seqs_in_prediction"` value from prediction_parameters file.

The `N` can be defined in parameters as `"max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION`.

parameters:
```
    {"Turbo-fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        "TurboFold": "TurboFold PARAMETERS"
        }
    }
```

## TurboFold
Select reference estimated full-length sequences without ambiguous bases.
Then continue as with Turbo-fast.

parameters:
```
    {"TurboFold": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        "TurboFold": "TurboFold PARAMETERS"
        }
    }
```
# rboAnalyzer
## Introduction
rboAnalyzer is a pipeline meant as a complementary step when BLAST algorithm
 was used to search for query that is non-coding RNA (ncRNA) with secondary structure (does not have to be known).
As the BLAST is general sequence alignment algorithm, it's results (output)
 is missing some very useful features in context of ncRNAs.

Also, the output apart from well scoring hits, usually contains sequence
 fragments (subject sequence only partially aligned to query).
Such hit may bring valuable information by capturing distant homology, or
it may be nothing.

With our rboAnalyzer we add information to such BLAST search to help researcher
 decide which hits are real ncRNAs and what their secondary structure might be.

## Getting help
Commandline reference can be viewed [here](../readme.md#help) (click on the dropdown) or accessed by
```shell
rboAnalyzer -h
```

In case of any questions, comments, bugreports or suggestions, create an issue or write to `marek.schwarz AT biomed.cas.cz`.


## Functionality overview
<img src="RBA_pipeline_overview.svg" width="700px" />

The rboAnalyzer has 3 stages:
1) Estimation of full-length RNA sequence from HSP (extension).
2) Estimation of homology of estimated full-length sequences to query sequence.
3) Prediction of secondary structures.

Each of these stages has dedicated section.

### Methods for estimation of full-length sequences from HSPs
The pipeline has 3 methods for estimating the full-length sequences from BLAST HSPs.

1) __simple__

    This means that location of estimated full-length sequence is computed from
     unaligned parts of query sequence on 5' and 3' ends of the HSP.

2) __locarna__

    In this method, the loci containing hit with flanking regions at the subject sequence is realigned to
     the query sequence with Locarna algorithm. The sequence
     aligned to the query is considered to be the estimated full-length sequence.

3) __meta__

    Here the two aforementioned methods are combined and estimated full-length sequences
     are scored with covariance model. The better scoring sequence is chosen.

#### ad simple)
<img src="figure_blast_simple.svg" width="700px" />

In the __simple__ mode we compute the location of the extended subject sequence according 
to the unaligned parts of the query sequence, i.e. those that were not aligned in HSP 
and flank the HSP alignment at both 5’ and 3’ ends. In this toy example we have the Plus/Plus 
BLAST HSP with a section of the query sequence between nucleotides 10 and 21 aligned 
to a section of the subject sequence between nucleotides 1000 and 1009. 
The task is to extend the partial HSP subject sequence between nucleotides 1000 and 1009 
to the length of the query sequence.

Suppose that the query sequence (the red bar in the figure of the example) is 50 bases long. 
Then the length of the unaligned part of the query at 5’ end is 9 nucleotides 
(subtract nucleotide positions 10 - 1) and the length of the unaligned part of 
the query at 3’ end is 29 nucleotides (subtract positions 50 - 21). 
The positions of the extended subject sequence at the whole subject sequence is computed by 
adding/subtracting the lengths of the unaligned parts of the query sequence to/from 3’/5’ 
ends of the partial HSP subject sequence, respectively. 
Then, 5’ and 3’ ends of the extended HSP sequence lie at nucleotides 991 (1000 – 9) 
and 1038 (1009 + 29), respectively. Because there are 2 gaps in HPS subject sequence, 
the resulting extended subject sequence will be 2 nucleotides shorter than the query sequence. 

Theoretically, the 2 extra nucleotides could be added at 3’ end of the extended subject 
sequence to make it as long as the query sequence. 
But this might not be biologically relevant as the gap could occur naturally instead 
of being caused by the alignment.


#### ad locarna)
<img src="figure_blast_locarna.svg" width="700px" />

With __locarna__ mode we first extract so called _supersequence_, which is
 region on subject sequence as with __simple__,
 additionally padded on 5' and 3' ends by extra sequence from the subject sequence.
This _supersequence_ is then realigned with Locarna algorithm to obtain the estimated full-length sequence.

The Locarna algorithm utilises possible pairings in it's computations,
 thus it is better suited to align RNAs then BLAST algorithm.
The Locarna is by default called with `struct-local=0`,
 `sequ-local=0` and `free-endgaps=++++` parameters.
Additionally, the information about matching nucleotides from BLAST HSPs
is used to construct so called anchor for the Locarna algorithm.
The anchor defines columns of alignment which are considered aligned.
As the anchor we consider consecutive series of matches of length at
 least `L` in BLAST alignment.
The default value of `L` is 7.
This way the alignment is anchored and the Locarna algorithm can align
 query to the _supersequence_. With the `free-endgaps=++++` option,
 the algorithm does not put penalty to unaligned ends of _supersequence_.
The estimated full-length sequence is the continuous part of _supersequence_ aligned to the query sequence
 (i.e. the subject sequence between the bases, inclusive, on subject sequence matching to the 5' terminal and 3' terminal bases).

#### add meta)
This approach combines the __simple__ and __locarna__. It computes both and 
 for each HSPs it chooses the estimated full-length sequence with higher score to covariance model.

### Estimation of homology
Here we compute score for relation between the estimated full-length sequence and query sequence.
The computation is based on aligning covariance model (CM) to each estimated full-length
 sequence with `cmalign` program from the Infernal package.

We've implemented 3 options on how to provide covariance model:

1) build with RSEARCH (default)

    By default, we build the covariance model from the query sequence (secondary structure predicted by RNAfold) and RIBOSUM matrix.
    The RIBOSUM is RIBOSUM65 by default and it can be changed in [alternative](config_how_to.md) `config.txt` file.

2) supply your own model (the `--cm_file` option)

    If the covariance model is known, it can be provided with `--cm_file` option.
    Only one model per file is allowed.
    
    Note that if you provide the covariance model, it will also be used in all methods for prediction of secondary structures using covariance models (those starting with `rfam`).

3) infer from Rfam (the `--use_rfam` option)

    The Rfam database is searched with query sequence for the best matching
    model (`cmscan`).

### Prediction of secondary structures
The rboAnalyzer can use multiple approaches (prediction methods) to predict secondary structures.
The prediction methods can be (roughly) divided to following groups:

- Predict structure independently of other estimated full-length sequences
    The advantage for these methods is robustness to possible improper parameter choice.
    - [rnafold](prediction_methods.md#rnafold)
    - [fq-sub](prediction_methods.md#fq-sub)
    - [rfam-Rc](prediction_methods.md#rfam-Rc)
    - [rfam-centroid](prediction_methods.md#rfam-centroid)
    - [rfam-sub](prediction_methods.md#rfam-sub)

- Use of selected estimated full-length sequences as reference
    - [centroid](prediction_methods.md#centroid)
    - [TurboFold](prediction_methods.md#TurboFold)
    - [Turbo-fast](prediction_methods.md#Turbo-fast)
    - [centroid-fast](prediction_methods.md#centroid-fast)

- Use of selected estimated full-length sequences to build consensus secondary structure
    - [C-A-r-Rc](prediction_methods.md#C-A-r-Rc)
    - [M-A-r-Rc](prediction_methods.md#M-A-r-Rc)
    - [C-A-U-r-Rc](prediction_methods.md#C-A-U-r-Rc)
    - [M-A-U-r-Rc](prediction_methods.md#M-A-U-r-Rc)
    - [C-A-sub](prediction_methods.md#C-A-sub)
    - [M-A-sub](prediction_methods.md#M-A-sub)

## Output

### Output formats
The rboAnalyzer is able to produce several output formats, most handy
  being the `.html`.
- html
    Stand-alone web page containing estimated full-length sequences and predicted secondary structures.
    If internet connection is available, it can be used to view respective
    genome loci for each BLAST HSP using NCBI SeqViewer.
- json
    Json-readable rboAnalyzer output (contains all data).
- csv
    Output table in comma separated values. Contains all important information
    including original HSP data, estimated full-length sequence location,
    sequence and predicted secondary structure(s).

## HTML output
In the head section there is report on basic input data and name of Rfam
 covariance model with best score to the provided query sequence.

The html output is organized around BLAST output.
Each BLAST HSP gets it's separate section with five parts:

  1) the text representation of BLAST HSP

  2) rboAnalyzer report with estimated full-length sequence indices and RSEARCH bit score

  3) the estimated full-length sequence itself

  4) one or multiple predicted secondary structures

  5) NCBI Sequence viewer (optional - by default only load button is shown)

The `html` outputs offers sorting, selecting sequences and structures and
  their export to fasta format or fasta-like format with predicted secondary structures in dot-bracket notation.
If internet connection is available, the NCBI genome browser can be used
  to explore synteny and known features of current genome.

The header for each BLAST HSP contains Accession.Version number (based on provided regular expression).
The header is also color-coded on color scale from green to red based on RSEARCH score.
This allows rapid identification of interesting or suspicious HSPs differing from others.

### Example output ideal case
<img src="html_ba_1.png" width="695px" />

The black arrows points to the HSP header (color indicating homology) and control buttons respectively.

### Example output with  loaded NCBI sequence viewer
<img src="html_ba_2.png" width="824px" />

### Notes
#### Control buttons
At the bottom of the view there are general control buttons which allow
  selecting and deselecting of sequences and structures, sorting by E-value
  and bulk initialization of NCBI Sequence Viewer.
- Select/Unselect all Seqs. (will select/unselect (check checkbox) all estimated full-length sequences)
- Select/Unselect all Structs (will select/unselect (check checkbox) all predicted structures)
- Export sel. Structs (will trigger download of selected predicted secondary structures in fasta-like format)
- Export sel. Seqs (will trigger download of selected estimated full-length sequences in fasta-like format)
- Sort Eval desc/asc (will sort BLAST HSPs according to E-value Ascending or Descending)
- View all Regions (will trigger loading of NCBI viewer for all (not yet loaded) HSPs)

#### Report structure

1) Inputs: query input file and BLAST input file
  Best matching mode from Rfam

2) Estimated full-length sequences, predicted secondary structures and other data

3) Command and parameters
  - executed commandline string
  - date and time of run
  - parameters

#### The fasta-like format containing secondary structures
```
>uid:N|ACCESSION.VERSIONdirection-method_name START-END (genome location)
SEQUENCE
SECONDARY_STRUCTURE

# - where the N is serial number of BLAST HSP
# - direcion can be "fw" for plus strand and "rc" for minus strand
# - genome-location is then location of found sequence on original genome
#   in START-END format where START is always lower index then END (direction is defined by "direction")
# the sequence is always 5' to 3' direction

>uid:104|CP006976.1fw-rfam-Rc 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((.(((..(((((((..((((.((((((.(....((((......(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))...))))....).)))))).))))...))))))).))))))))))))))......
>uid:104|CP006976.1fw-rnafold 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((((((.((((((......)))......................(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))(((((((((....)))))))))......))).)))..))))))))))))......
```
#### The NCBI sequence Viewer
The NCBI sequence viewer works only if internet connection is available.
It may take some time to load (especially with large genomes) and when the report
 contains many BLAST hits it may require more substantial amount of RAM.
The data for the sequence viewer are not saved across browser sessions.

## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund - Project "ELIXIR-CZ: Budování kapacit" (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](logolink_OP_VVV_hor_barva_eng.jpg)
# The configuration file
 The configuration file contains pipeline wide settings such are paths to tools
 and databases used.

 It should be used if non-default installation path is used.

 ## Configuration
 Configuration is done using the file `config.txt` placed in
 `INSTALL_LOCATION/rna_blast_analyze/BR_CORE/config.txt` or in custom location
 and running the pipeline with the `--config_file PATH_to_custom_LOCATION` argument.

 ## Specification

There are 3 sections

- TOOL_PATHS

  defines path to executable __parent__ directory

- DATA

  defines paths to data and databases

Each section is specified by its name in square brackets.

## Example
 ```
[TOOL_PATHS]
clustal = /usr/bin/
[DATA]
rfam_dir = /home/user/rfamdir/
 ```

## Available settings
 Setting name with executable(s) which should be accessible from provided location (in brackets).

### TOOL_PATHS
 - refold (`refold.pl`)
 - infernal (`cmalign`, `cmbuild`, `cmfetch`, `cmscan`)
 - muscle (`muscle`)
 - clustal (`clustalo`)
 - locarna (`locarna`)
 - viennarna_bin (`RNAfold`, `RNAalifold`, `RNAplot`)
 - mfold (`hybrid-ss-min`)
 - centroid (`centroid_homfold`)
 - turbofold (`TurboFold`)

### DATA
 - rsearch_ribosum

  specify full path to RSEARCH matrix file (default: `RIBOSUM65.mat`).

 - rfam_dir

  custom path where Rfam database should be stored. (default: `INSTALL_LOCATION/rna_blast_analyze/3rd_party_source/rfamdb/`)

 - rfam_url

  custom url from where the `rfam.cm.gz` (database dump) should be downloaded.
  This url is also checked when update is requested (wget timestamp is checked).
  (default: `ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/rfam.cm.gz`)
  

 - html_template_dir

  custom location for html jinja2 template (default: `INSTALL_LOCATION/rna_blast_analyze/BR_core/output/`)

 - rnastructure_datapath

  datapath for RNAstructure (see installation notes for RNAstructure https://rna.urmc.rochester.edu/Text/Thermodynamics.html)
