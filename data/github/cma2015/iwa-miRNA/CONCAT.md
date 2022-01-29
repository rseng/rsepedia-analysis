
# iwa-miRNA:  A Web-based Platform for Interactive Annotation of Plant MiRNAs

[![webserver](https://img.shields.io/badge/Web_server-ready-red.svg)](http://iwa-miRNA.omicstudio.cloud/) [![docker](https://img.shields.io/badge/Docker_image-ready-red.svg)](https://hub.docker.com/r/malab/iwa-mirna/) [![source](https://img.shields.io/badge/Source_code-support-blue.svg)](https://github.com/cma2015/iwa-miRNA/tree/master/Source_code) [![testdata](https://img.shields.io/badge/Test_data-support-blue.svg)](https://github.com/cma2015/iwa-miRNA/tree/master/Test_data)

- iwa-miRNA allows users to generate a comprehensive collection of miRNA candidates, and to interrogate miRNA annotation in a straightforward way, without the need for computational skills.

- iwa-miRNA Docker image is available at https://hub.docker.com/r/malab/iwa-mirna. Source codes and user manual are available at https://github.com/cma2015/iwa-miRNA. The web server of iwa-miRNA is accessible at http://iwa-miRNA.omicstudio.cloud/.

- The iwa-miRNA is composed with three functional modules: MiRNA Compilation, MiRNA Selection, and Manual Curation. It also provides users with some useful tools for downstream exploratory analysis. More details regarding these functional modules can be found [here](https://github.com/cma2015/iwa-miRNA/blob/master/Tutorials/Modules.md).

<img src="Tutorials/_images/Graphical_summary.png" alt="Graphical summary of iwa-miRNA" style="zoom:18%">

## How to use iwa-miRNA

- Tutorial for iwa-miRNA: http://iwa-miRNA.omicstudio.cloud/static/assets/html/index.html
- Test datasets referred in the tutorials for iwa-miRNA: https://github.com/cma2015/iwa-miRNA/tree/master/Test_data

## Abbreviations in iwa-miRNA

- SVM: support vector machine;
- HT criteria: high-throughput criteria;
- PsRNA: [plant small RNA genes](http://plantsmallrnagenes.science.psu.edu/);
- MFE: minimal free energy;
- AMFE: adjusted minimal free energy;
- sRNA-Seq: small RNA sequencing.

## Changelog

- 2020-11-30: Update the description of output results in the tutorial.
- 2020-06-24: A demo server of iwa-miRNA was released for users running small RNA sequencing datasets.
- 2020-06-20: Source codes and Docker image of iwa-miRNA were released for the first time.

## How to access help

- For any feedback and tool suggestions, please feel free to leave a message at Github [issues](https://github.com/cma2015/iwa-miRNA/issues). We will try our best to deal with all issues as soon as possible.
- In addition, if any suggestions are available, feel free to contact: **Ting Zhang** [zting135@gmail.com](mailto:zting135@gmail.com) or **Chuang Ma** [chuangma2006@gmail.com](mailto:chuangma2006@gmail.com)

## Citation

T. Zhang, J. Zhai, X. Zhang, L. Ling, M. Li, S. Xie, M. Song, C. Ma, Interactive Web-based Annotation of Plant MicroRNAs with iwa-miRNA, *Genomics, Proteomics & Bioinformatics* (2021), doi: https://doi.org/10.1016/j.gpb.2021.02.010

### Module I 

This module generates a comprehensive collection of miRNA candidates by aggregating already annotated miRNAs from four plant miRNA databases (i.e., miRBase, PmiREN, sRNAanno, and PsRNA) and predicted miRNAs from user-submitted sRNA-Seq data.

|       Tool       |                     **Input**                      |                            Output                            |                         Applications                         |
| :--------------: | :------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| *miRNARetrival*  |        Name of species and miRNA databases         | Already annotated miRNAs<br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/Overview.html">Overview</a> | Aggregate annotated miRNAs provided by different miRNA databases |
|  *miRNAPredict*  |    SRA accession numbers or uploaded fastq file    | Predicted miRNAs<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNAPredict_output.txt">Prediction</a> |              Predict miRNAs from sRNA-Seq data               |
| *genomeRetrival* | Name of species or genome sequences and annotation |      Path of formatted genome sequences and annotation       |             Get genome sequences and annotation              |
| *miRNATranslate* |     Output from miRNARetrival and miRNAPredict     | miRNA and miRNA precursors with a uniform format<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNATranslate_output.txt">Aggregation</a><br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/Aggregation.html">Report</a> | Translate annotated and predicted miRNAs into the genomic coordinate system |

### Module II

This module selects a subset of miRNA candidates according to the high-throughput criteria and/or using an ML-based approach.

|    **Tools**     |         **Input**          |                          **Output**                          |           Applications            |
| :--------------: | :------------------------: | :----------------------------------------------------------: | :-------------------------------: |
| *miRNASelection* | Output from miRNATranslate | Selected miRNAs<br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNASelection_output.txt">Output</a> | Select promising miRNA candidates |

### Module III

This module provides the information for all miRNA candidates  for rapid curation of the quality of selected miRNAs.

|    **Tools**     |          **Input**          |                          **Output**                          |          Applications           |
| :--------------: | :-------------------------: | :----------------------------------------------------------: | :-----------------------------: |
| *manualCuration* | Output from MiRNA Selection | Summary and report pages<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/manualCuration_output.html">Output</a> | Determine the quality of miRNAs |



> The table above briefly describes each function in this module. The sample data are displayed through a hyperlink. Users can preview the results by clicking it. Due to the limitation of file size, some results may not be displayed properly. For more detailed descriptions about inputs and outputs, please refer to [tutorial for iwa-miRNA](http://iwa-miRNA.omicstudio.cloud/static/assets/html/index.html) or the interpretation of each functional analysis interface in [web server](http://iwa-miRNA.omicstudio.cloud/).

# enaBrowserTools

enaBrowserTools is a set of scripts that interface with the ENA web services to download data from ENA easily, without any knowledge of scripting required.

*Important: Python 2 based scripts in the python/ folder are deprecated. Use the scripts in python3/

Please note that v1.5.5 is the last update that will support python 2.7. Due to the 2019 scheduled retirement of python 2.7, the next release (1.6.0) will be restructured to hold only python 3 scripts.*

# License

Copyright 2017 EMBL - European Bioinformatics Institute Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Installation and Setup

## Python installation

Both Python 2 and Python 3 scripts are available.  The Python 2 scripts can be found in the "python" folder, and the Python 3 scripts in the "python3" folder.   

To run these scripts you will need to have Python installed.  You can download either Python 2 or Python 3 from [here](https://www.python.org/downloads/). If you already have Python installed, you can find out which version when you start the python interpreter.  If using Python 2, we suggest you upgrade to the latest version if you don't already have it: 2.7.

Note that EBI now uses HTTPS servers. This can create a problem when using Python 3 on a Mac due to an oft-missed
installation step. Please run the "Install Certificates.command" command to ensure your Python 3 installation on
the Mac can correctly authenticate against the servers. To do this, run the following from a terminal window, updating
the Python version with the correct version of Python 3 that you have installed:
open "/Applications/Python 3.6/Install Certificates.command"

We have had a report from a user than when Python 3 was installed using homebrew, the following code needed to be run instead:
```
# install_certifi.py
#
# sample script to install or update a set of default Root Certificates
# for the ssl module.  Uses the certificates provided by the certifi package:
#       https://pypi.python.org/pypi/certifi

import os
import os.path
import ssl
import stat
import subprocess
import sys

STAT_0o775 = ( stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
             | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP
             | stat.S_IROTH |                stat.S_IXOTH )

openssl_dir, openssl_cafile = os.path.split(
    ssl.get_default_verify_paths().openssl_cafile)

print(" -- pip install --upgrade certifi")
subprocess.check_call([sys.executable,
    "-E", "-s", "-m", "pip", "install", "--upgrade", "certifi"])

import certifi

# change working directory to the default SSL directory
os.chdir(openssl_dir)
relpath_to_certifi_cafile = os.path.relpath(certifi.where())
print(" -- removing any existing file or link")
try:
    os.remove(openssl_cafile)
except FileNotFoundError:
    pass
print(" -- creating symlink to certifi certificate bundle")
os.symlink(relpath_to_certifi_cafile, openssl_cafile)
print(" -- setting permissions")
os.chmod(openssl_cafile, STAT_0o775)
print(" -- update complete")
```

## Installing and running the scripts

Download the [latest release](https://github.com/enasequence/enaBrowserTools/releases/latest) and extract it to the preferred location on your computer. You will now have the enaBrowserTools folder, containing both the python 2 and 3 option scripts.  If you are using a Unix/Linux or Mac computer, we suggest you add the following aliases to your .bashrc or .bash_profile file. Where INSTALLATION_DIR is the location where you have saved the enaBrowserTools to and PYTHON_CHOICE will depend on whether you are using the Python 2 or Python 3 scripts.

```
alias enaDataGet=INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaDataGet
alias enaGroupGet=INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaGroupGet
```

This will allow you to run the tools from any location on your computer.

You can run install and run these scripts on Windows as you would in Unix/Linux using [Cygwin](https://cygwin.com). If you have python installed on your Windows machine, you can run the python scripts directly without Cygwin. However the call is a bit more complicated.

For example, instead of calling ```enaDataGet```

you would need to call ```python INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaDataGet.py```

We will look more into the Windows equivalent of aliases to run batch files from the command line and hopefully be able to provide a better solution to Windows users shortly.

## Setting up for Aspera

*Important: there has been a change to using aspera in version 1.4.*

If you wish to use Aspera to download read or analysis files, you will need to use the aspera_settings.ini file.  Please save it to a static location on your local computer, and edit the file to include the location of your aspera binary (ASPERA_BIN) and the private key file (ASPERA_PRIVATE_KEY):

```
[aspera]
ASPERA_BIN = /path/to/ascp
ASPERA_PRIVATE_KEY = /path/to/aspera_dsa.openssh
ASPERA_OPTIONS =
ASPERA_SPEED = 100M
```

There are two command line flags/options available if you wish to use aspera for download. These are:
1. -a (or --aspera): use this flag if you'd like to download with aspera.
2. -as (or --aspera-settings): use this option if you'd like to specify the location of your aspera settings file. If you use this option, you don't need to use the --aspera flag.

If you use the --aspera-settings option, you don't need to also use the --aspera flag, e.g:

```
enaDataGet -f fastq -as /path/to/aspera_settings.ini ACCESSION
```

If you don't wish to specify the location of the aspera settings file each time you use the scripts, you have the option to either set an ENA_ASPERA_INIFILE environment variable to save the location:

```
export ENA_ASPERA_INIFILE="/path/to/aspera_settings.ini"
```

or you can use the default location for the file, this is the enaBrowserTools directory.  Note that if you use this option, you will have to be careful about how you update your scripts so that you don't overwrite your aspera settings file.

Using just the --aspera flag will result in the scripts looking first for the ENA_ASPERA_INIFILE environment variable, and second for the default file location.

```
enaDataGet -f fastq -a ACCESSION
```

Regardless of which option you have selected, if the aspera settings file cannot be found or the licence key file declared within your settings file does not exist, the scripts will default to using FTP for the download.

# Command line

There are two main tools for downloading data from ENA:  enaDataGet and enaGroupGet.  

## enaDataGet

This tool will download all data for a given sequence, assembly, read or analysis accession or WGS set.  Usage of this tool is described below.  Note that unless a destination directory is provided, the data will be downloaded to the directory from which you run the command. When using an assembly, run, experiment or analysis accession, a subdirectory will be created using that accession as its name.

Accepted WGS set accession formats are:
- AAAK03
- AAAK03000000
- AAAK
- AAAK00000000

```
usage: enaDataGet [-h] [-f {embl,fasta,submitted,fastq,sra}] [-d DEST] [-w]
                  [-m] [-i] [-a] [-as ASPERA_SETTINGS] [-v]
                  accession

Download data for a given accession

positional arguments:
  accession             Sequence, coding, assembly, run, experiment or
                        analysis accession or WGS prefix (LLLLVV) to download

optional arguments:
  -h, --help            show this help message and exit
  -f {embl,fasta,submitted,fastq,sra}, --format {embl,fasta,submitted,fastq,sra}
                        File format required. Format requested must be
                        permitted for data type selected. sequence, assembly
                        and wgs accessions: embl(default) and fasta formats.
                        read group: submitted, fastq and sra formats. analysis
                        group: submitted only.
  -d DEST, --dest DEST  Destination directory (default is current running
                        directory)
  -w, --wgs             Download WGS set for each assembly if available
                        (default is false)
  -e, --extract-wgs     Extract WGS scaffolds for each assembly if available
                        (default is false)
  -exp, --expanded      Expand CON scaffolds when downloading embl format
                        (default is false)
  -m, --meta            Download read or analysis XML in addition to data
                        files (default is false)
  -i, --index           Download CRAM index files with submitted CRAM files,
                        if any (default is false). This flag is ignored for
                        fastq and sra format options.
  -a, --aspera          Use the aspera command line client to download,
                        instead of FTP.
  -as ASPERA_SETTINGS, --aspera-settings ASPERA_SETTINGS
                        Use the provided settings file, will otherwise check
                        for environment variable or default settings file
                        location.
  -v, --version         show program's version number and exit
```

## enaGroupGet

This tool will allow you to download all data of a particular group (sequence, WGS, assembly, read or analysis) for a given sample or study accession. You can also download all sequence, WGS or assembly data for a given NCBI tax ID. When fetching data for a tax ID, the default is to only search for the specific tax ID, however you can use the subtree option to download the data associated with either the requested taxon or any of its subordinate taxa in the NCBI taxonomy tree.

Downloading read and analysis data by tax ID is currently disabled. We will be adding a data volume sanity check in place before we enable this functionality.

Usage of this tool is described below.  A new directory will be created using the provided accession as the name, and all data will be downloaded here. There will also be a separate subdirectory created for each assembly, run and analysis being fetched.  Note that unless a destination directory is provided, this group directory will be created in the directory from which you run the command.  

```
usage: enaGroupGet [-h] [-g {sequence,wgs,assembly,read,analysis}]
                   [-f {embl,fasta,submitted,fastq,sra}] [-d DEST] [-w] [-m]
                   [-i] [-a] [-as ASPERA_SETTINGS] [-t] [-v]
                   accession

Download data for a given study or sample, or (for sequence and assembly) taxon

positional arguments:
  accession             Study or sample accession or NCBI tax ID to fetch data
                        for

optional arguments:
  -h, --help            show this help message and exit
  -g {sequence,wgs,assembly,read,analysis}, --group {sequence,wgs,assembly,read,analysis}
                        Data group to be downloaded for this
                        study/sample/taxon (default is read)
  -f {embl,fasta,submitted,fastq,sra}, --format {embl,fasta,submitted,fastq,sra}
                        File format required. Format requested must be
                        permitted for data group selected. sequence, assembly
                        and wgs groups: embl and fasta formats. read group:
                        submitted, fastq and sra formats. analysis group:
                        submitted only.
  -d DEST, --dest DEST  Destination directory (default is current running
                        directory)
  -w, --wgs             Download WGS set for each assembly if available
                        (default is false)
  -e, --extract-wgs     Extract WGS scaffolds for each assembly if available
                        (default is false)
  -exp, --expanded      Expand CON scaffolds when downloading embl format
                        (default is false)
  -m, --meta            Download read or analysis XML in addition to data
                        files (default is false)
  -i, --index           Download CRAM index files with submitted CRAM files,
                        if any (default is false). This flag is ignored for
                        fastq and sra format options.
  -a, --aspera          Use the aspera command line client to download,
                        instead of FTP.
  -as ASPERA_SETTINGS, --aspera-settings ASPERA_SETTINGS
                        Use the provided settings file, will otherwise check
                        for environment variable or default settings file
                        location.
  -t, --subtree         Include subordinate taxa (taxon subtree) when querying
                        with NCBI tax ID (default is false)
  -v, --version         show program's version number and exit
```

# Tips

From version 1.4, when downloading read data if you use the default format (that is, don't use the format option), the scripts will look for available files in the following priority: submitted, sra, fastq.

A word of advice for read formats:
- submitted: only read data submitted to ENA have files available as submitted by the user.
- sra:  this is the NCBI SRA format, and is the format in which all NCBI/DDBJ data is mirrored to ENA.
- fastq:  not all submitted format files can be converted to FASTQ

# Problems

For any problems, please contact the ENA helpdesk (https://www.ebi.ac.uk/ena/browser/support) with 'enaBrowserTools' in your subject line.

We have had a couple of reports that the R2 FASTQ files are not always successfully downloading for paired runs. We have been unable to replicate this problem and have therefore exposed the error message to you should there be any download failure of read/analysis files via FTP or Aspera. If you get one of these errors, please copy the error into the support form.
**sRNA-Seq_test_data**

This file contains raw sequencing reads in FASTQ format, which can be compressed into **sRNA-Seq_test_data.zip** and used as input of **miRNAPredict** module. In detail, two Arabidopsis sRNA-seq data (SRR11347201 and SRR11829907) were selected, and the first 10 million rows were packaged and named as test1.fastq.gz and test2.fastq.gz.

**sample_information.txt**

This file contains the tissue name of corresponding samples involved in **sRNA-Seq_test_data**. It has two columns with tab-delimited format. The first column is sample names, and the second column is tissue names.

**gene_description.txt**

This file contains genes descriptions and is used as input of **manualCuration** function. It has two columns with tab-delimited format. The first column is the gene names and the second column is the gene description. Gene descriptions were extracted from [Arabidopsis GFF3 annotation file](ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.gff3.gz).

**genome_annotation.txt**

This file contains gene function annotation information, which can be uploaded to iwa-miRNA server and used as inputs of the **manualCuration** and **associationAnalysis** modules. It has seven columns with tab-delimited format, which represent chromosomes, start sites, end sites, ID, strand, gene types, and detailed types. These gene attributes were extracted from [Arabidopsis GFF3 annotation file](ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.gff3.gz). TEs (Transposed elements) were downloaded from [TAIR10_Transposable_Elements](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt).

**SNPs_in_Arabidopsis.txt**

This file contains part of SNPs from The 1001 Genomes Project and is used as input of the **sequenceVariation** function. It has five columns with tab-delimited format, representing chromosomes, start sites, ID, reference, and alternative alleles, respectively. This information can be downloaded from [Ensembl Plants Variations](https://plants.ensembl.org/biomart/martview).