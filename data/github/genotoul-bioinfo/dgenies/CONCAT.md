[![GitHub release](https://img.shields.io/github/release/genotoul-bioinfo/dgenies.svg)](https://GitHub.com/genotoul-bioinfo/dgenies/releases/) [![PyPI version](https://badge.fury.io/py/dgenies.svg)](https://badge.fury.io/py/dgenies) [![GitHub license](https://img.shields.io/github/license/genotoul-bioinfo/dgenies.svg)](https://github.com/genotoul-bioinfo/dgenies/blob/master/LICENSE.txt)
 [![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/) [![Documentation Status](https://readthedocs.org/projects/dgenies/badge/?version=latest)](http://dgenies.readthedocs.io/?badge=latest) [![Github all releases](https://img.shields.io/github/downloads/genotoul-bioinfo/dgenies/total.svg)](https://GitHub.com/genotoul-bioinfo/dgenies/releases/) [![Downloads](http://pepy.tech/badge/dgenies)](http://pepy.tech/project/dgenies)



D-Genies
--------

Dot plots are widely used to quickly compare sequence sets. They provide a synthetic similarity overview, highlighting repetitions, breaks and inversions. Different tools have been developed to easily generated genomic alignment dot plots, but they are often limited in the input sequence size. D-GENIES is a standalone and WEB application performing large genome alignments using [minimap2](https://github.com/lh3/minimap2) software package and generating interactive dot plots. It enables users to sort query sequences along the reference , zoom in the plot and download several image, alignment or sequence files. D-GENIES is an easy to install open source software package (GPL) developed in Python and JavaScript.

How to use?
-----------

You can use the demo instance [here](http://dgenies.toulouse.inra.fr).

Or you can install your own instance. The install documentation is available [here](http://dgenies.toulouse.inra.fr/install).

Documentation
-------------

Full documentation, including code API, is available [here](https://dgenies.readthedocs.io/en/latest/index.html).

How to cite?
------------

Cabanettes F, Klopp C. (2018) D-GENIES: dot plot large genomes in an interactive, efficient and simple way. PeerJ 6:e4958 https://doi.org/10.7717/peerj.4958

Credits
-------

D-Genies is developed by the [Bioinfo team](http://bioinfo.genotoul.fr/index.php/about-us/) of the [Genotoul platform](http://www.genotoul.fr/) ([INRA](http://www.inra.fr/)).Install your own instance
=========================

{% if version != "" %}
Latest available version: **{{version}}**
{% endif %}

Linux
-----

### Install

Install in 1 step:

As root:

    pip3 install dgenies

Or as simple user:

    pip3 install dgenies --user

Alternatively, you can install it manually:

    git clone https://github.com/genotoul-bioinfo/dgenies
    cd dgenies
    pip3 install -r requirements.txt
    python3 setup.py install

### Upgrade

#### Standalone mode

    pip3 install dgenies --upgrade
    
Add `--user` flag if you have not root access.

#### Webserver mode

    dgenies clear -c
    pip3 install dgenies --upgrade

Then, you need to restart your webserver.



### Requirements

D-Genies requires python >= 3.5, < 3.10 to run.

We use minimap2 (and mashmap2 on linux) for mapping. A binary of each of them is shipped with the program. Alternatively, you can
install your own from [minimap2](https://github.com/lh3/minimap2) and [mashmap2](https://github.com/marbl/MashMap/) repositories.

Some python modules are required (will be automatically installed by the commands above):

    Flask==1.0.*
    Jinja2~=2.11.3
    numpy
    requests~=2.20.1
    biopython>=1.70
    psutil~=5.6.6
    tendo==0.2.*
    matplotlib>=2.1.*
    intervaltree==2.1.*
    Markdown==2.6.*
    pyyaml~=5.4.1

Additional modules for webserver mode:

    Flask-Mail==0.9.*
    peewee==2.10.2
    python-crontab>=2.2.*

And if you use a cluster (webserver mode):

    drmaa==0.7.*

In webserver mode, you must install `mysqlclient` python module (will not be installed automatically) if you use mysql as RDBM.


Windows
-------

We provide an installer to install D-Genies. You can download it 
[here]({%if win32 %}{{win32}}{% else %}https://github.com/genotoul-bioinfo/dgenies/releases{% endif %}).

All requirements are present inside the package, so you don't have to do anything else.

### System requirements

You need Windows 7 or newer, 64 bits architecture.


How to start
-------------

You can launch DGenies in `standalone` mode or in `webserver` mode. You can use standalone mode if
you launch it locally and only one user launch jobs at once. If you are several users to use it
simultaneously or if you run it on a server, you must run it in webserver mode.

### Standalone mode

Unix: start with the command below:

    dgenies run

Optional arguments:

`-p <port>` run in a specified port (default: 5000)
`--no-browser` don't start the browser automatically

Windows: just click on the launcher in the desktop or into the install folder.

### Webserver mode

*Note: this mode is only available for Unix systems and will NOT work on MS Windows.*

#### Recommended method

Flask webserver (which is used in standalone mode) is not recommended in production servers.
So, we recommend using the WSGY module of Apache (or µWSGI + nginx, not documented here).

Once dgenies is installed, you just need to use the `/var/www/dgenies.wsgi` file into your apache
virtualhost file.

Here is an example of configuration file for apache:

    <VirtualHost *>
        ServerName <url>

        WSGIDaemonProcess dgenies user=<user> group=<group> threads=8
        WSGIScriptAlias / /var/www/dgenies/dgenies.wsgi

        <Directory /var/www/dgenies>
            WSGIProcessGroup dgenies
            WSGIApplicationGroup %{GLOBAL}
            Order deny,allow
            Allow from all
        </Directory>
    </VirtualHost>

With:
`<url>`: the URL of your instance
`<user>`: the user who launch the server
`<group>`: the group who launch the server

#### Debug method

For debug or for development only, you can launch dgenies through flask in webserver mode:

    dgenies run -m webserver

Optional parameters:

`-d` run in debug mode
`-o <IP>` specify the host into run the application (default: 127.0.0.1, set 0.0.0.0 for distant access)
`-p <port>` run in a specified port (default: 5000)
`--no-crons` don't run the crons automatically
`--no-browser` don't start the browser automatically (always true if *-d* option is given)



Running with a cluster
----------------------

If you want to run jobs on a cluster, some configuration is required. We only support SLURM and SGE schedulers. But please note that only SLURM scheduler has been fully tested.

Note: submitting jobs on a cluster is only available for webserver mode.

Jobs are submitted throw the DRMAA library. So you need it for your scheduler. Please see [configuration below](#cluster) to define path to this library.

Also, scripts for preparing data must be moved in a location accessible by all nodes of your cluster. You must move them in a same folder and set the full path to `all_prepare.py` script (full path on the nodes of the cluster) in the configuration file ([see below](#cluster)).

To get these scripts, follow the commands below:

    curl https://raw.githubusercontent.com/genotoul-bioinfo/dgenies/v{{version}}/get_cluster_scripts.py > get_cluster_scripts.py
    python get_cluster_scripts.py -d <dir>

With `<dir>`: the folder into save the scripts (must be accessible by cluster nodes).



Configuration
-------------

Changing the default configuration is not required for standalone mode, but you can want to custom some parts of the program.

Configuration file location:  
* Linux:  
    * `/etc/dgenies/application.properties` if installed with root access  
    * `~/.dgenies/application.properties` else  
* Windows:  
    * `application.properties` file of the install folder

The file is divided in 9 parts described below.

To change this file, please copy it into `application.properties.local` (at the same location) to avoid erase of the file on upgrades.

### Global

Main parameters are stored into this section:

* `config_dir`: where configuration file will be stored.
* `upload_folder`: where uploaded files will be stored.
* `data_folder`: where data files will be stored (PAF files and other files used for the dotplot).
* `threads_local`: number of threads to use for local jobs.
* `web_url`: public URL of your website.
* `max_upload_size`: max size allowed for query and target file (-1 to avoid the limit) - size uncompressed.
* `max_upload_size_ava`: max size allowed for target file for all-vs-all mode (only target given, -1 to avoid the limit) - size uncompressed.
* `max_upload_file_size`: max size of the uploaded size (real size of the file, compressed or not, -1 to avoid the limit).

For webserver mode only (ignored in standalone mode):

* `batch_system_type`: local for run all jobs locally, sge or slurm to use a cluster scheduler.

### Debug

Some parameters for debug:

* `enable`: True to enable debug
* `log_dir`: folder into store debug logs

### Cluster

This section concerns only the webserver mode with *batch_system_type* not set to *local*.

* `drmaa_lib_path`: absolute path to the drmaa library. Required as we use the DRMAA library to submit jobs to the cluster.
* `native_specs`: how to set memory, time and number of CPU on the cluster (should be kept as default).

By default, small jobs are still launched locally even if *batch_system_type* is not set to *local*, and if not too much of these jobs are running or waiting. This limit can be customized:

* `max_run_local`: max number of jobs running locally (if this number is reached, future jobs will be submitted to the cluster regardless of their size). Set to 0 to run all jobs on the cluster.
* `max_wait_local`; max number of jobs waiting for a local run (if this number is reached, future jobs will be submitted to the cluster regardless of their size). Set to 0 to run all jobs on the cluster.

You can also customize the size from which jobs are submitted on the cluster. If only one of these limit is reached, the job is submitted on the cluster.

* `min_query_size`: minimum size for the query (uncompressed).
* `min_size_target`: minimum size for the target (uncompressed).

Other parameters:

* `prepare_script`: absolute path to the all_prepare.py script downloaded in the section [above](#running-with-a-cluster).
* `python3_exec`: path to python3 executable on the cluster.
* `memory`: max memory to reserve on the cluster.
* `memory_ava`: max memory to reserve on the cluster un all-vs-all mode (should be higher than memory).
* `threads`: number of threads for launching jobs on the cluster (must be a divider of the memory).

### Database

This section concerns only the webserver mode.

In webserver mode, we use a database to store jobs.

* `type`: sqlite or mysql. We recommend mysql for better performances.
* `url`: path to the sqlite file, or url to the mysql server (localhost if the mysql server is on the same machine).

If type is mysql, some other parameters must be filled:

* `port`: port to connect to mysql (3306 as default).
* `db`: name of the database.
* `user`: username to connect to the database.
* `password`: the associated password.

### Mail

This section concerns only the webserver mode.

At the end of the job, a mail is send to the user to advise him of the end of the job.

* `status`: mail to use for status mail.
* `reply`: mail to use as reply to.
* `org`: name of the organisation who send the mail.
* `send_mail_status`: True to send mail status, False else. Should be True in production.

### Cron

This section concerns only the webserver mode.

We use crons to launch the local scheduler who start jobs, and some script that clean old jobs.

* `clean_time`: time at which we launch the clean job (example: 1h00)
* `clean_freq`: frequency of the clean job execution

### Jobs

This section concerns only the webserver mode.

Several parameters for jobs:

* `run_local`: max number of concurrent jobs launched locally.
* `data_prepare`: max number of data prepare jobs launched locally.
* `max_concurrent_dl`: max number of concurrent upload of files allowed.

### Example

Here, you can fill example data. At least target is required to enable example data.

Fill for target and query the absolute local path of the file. This path will not be shown to the client. Only the file name will be shown.

If at least target is filled, a button "Load example" will be shown in the run form. Click on it will load example data in the form.

### Analytics

Set `enable_logging_runs` to True will enable storage of analytics data. It stores for each job creation date, user mail, size of query and target, and batch type.


Customize your installation
---------------------------

You can make easily make some changes to the application, described below. To do so, you must first clone the D-Genies repository (except if changes can be done in the current installation, see below):

    git clone https://github.com/genotoul-bioinfo/dgenies
    
The created folder is named the `D-Genies repository` in the text below.

Then, you make changes described below. When done, you can easily install the cusomized version with this command:

    pip install .
    
Add the `--upgrade` flag if you already have D-Genies installed.


### Add or change alignment tools

D-Genies uses minimap2 as default aligner, but you can add other tools, or replace minimap2 by another tool. You can also change minimap2 executable path.

If your tool needs a parser (see below), you must customize the installation (see above). For other changes, you can make them on your current installation.

Tools definition YAML file:  
* Linux:  
    * `/etc/dgenies/tools.yaml` if installed with root access  
    * `~/.dgenies/tools.yaml` else  
* Windows:  
    * `tools.yaml` file of the install folder
    
To change this file, please copy it into `tools.yaml.local` (at the same location) to avoid to erase the file on upgrades.

Note: you can also edit the `tools.yaml` file of the D-Genies repository, which will be installed then (if you customize the installation). Edit it without renaming it.

For each tool, you can or must define several properties described bellow.

#### Executable path

Required.

Always set tool binary path in the `exec` property. If it differs for cluster execution, set it in the `exec_cluster` property.

#### Command line skeleton

Required.

The command line skeleton defines how to launch the program, in the `command_line` property. The skeleton must include the following tags:
  
* `{exe}`: will be replaced by the executable path
* `{target}`: will be replaced by the target fasta file
* `{query}`: will be replaced by the query fasta file
* `{out}`: will be replaced by the output file

If the tool can be multithreaded, add the `{threads}` tag which will be replaced by the number of threads to use.

If your program is able to use all-vs-all mode (target versus itself), define the same skeleton in the `all_vs_all` property. All tag described hereover must be set except the `{query}` tag.

#### Threads

Required.

Defines how much threads to use for this tool. Set it in the `threads` property for local executions, and the `threads_cluster` property for cluster execution (if different from the local one).

#### Max memory

Optional.

If the tool requires less memory than defined in the configuration file, you can add the `max_memory` property. The unit is Gb.

#### Parser

Optional.

If the tool does not output PAF format, you must define a function in the `src/dgenies/lib/parsers.py` file of the D-Genies repository. This function get the input file as first parameter and outputs a valid PAF file (output file, second parameter). You must reference the function name in the `parser` property.

#### Split before

Optional. Default: False

For some tools (like minimap2), splitting the query on 10 Mb blocks improves performances (blocks are merged after mapping). To enable this for the tool, set the `split_before` property to True.

#### Help

Optional.

Define a message to show aside the tool name in the run form. Set it in the `help` property.

#### Order

Optional. Default: random.

Define in which order we show tools in the run form. Set it in the `order` property.


### Add new formats

In `Plot alignment` mode in run form ([see here](/documentation/run#plot-alignment-mode)), we propose by default only PAF and MAF formats. It's easy to add new formats.

Just define 2 functions:

* Add the first one in the `src/dgenies/lib/validators.py` file of the D-Genies repository. It takes only one argument: the input file. It checks if the file format is correct and returns True in this case, else it returns False. The function name must be the same as the expected input file extension.  
* Add the second one in the `src/dgenies/lib/parsers.py` file of the D-Genies repository. It takes two arguments: the input and the output file. It convert the input file into a valid PAF file. The function name must be same as the previous one.


Maintenance
-----------

The `dgenies` command can be used to do some maintenance staff.

**Clear all jobs:**

    dgenies clear -j [--max-age <age>]

`--max-age` (opt): set the max age of jobs to delete (default: 0, for all)

**Clear all log:**

    dgenies clear -l

**Clear crons (webserver mode):**

    dgenies clear -c
    
**Display message on run form:**

You can display a message at the top of the run form. It can be used to add extra informations for user, or for prevent him for problem or for a maintenance on your instance.

    dgenies inforun -m "message to display" -t [warn, critical, info, success]
    
`-m <message>`: message to display. Html allowed.  
`-t <type>`: type of message: warn (orange background), critical (red background), info (blue background) or success (green background).

Remove the message by:

    dgenies inforun -c

Gallery
-------

Note: gallery is only available in webserver mode.

To add a job to the gallery, copy illustrating picture file into the *gallery* folder inside the data folder (*~/.dgenies/data/gallery* as default, create it if not exists). Then use the *dgenies* command:

    dgenies gallery add -i <id_job> -n <name> -q <query_name> -t <target_name> -p <pict_filename>

With:
`id_job`: the name of the job
`name`: name of the job to show in the gallery
`query_name`: name of the query
`target_name`: name of the target
`pict_filename`: filename added in the gallery folder (without path)

You can also delete an item from the gallery:

    dgenies gallery del -i <id_job>

or:

    dgenies gallery del -n <name>

With `id_job` and `name` as described above. You can add the `--remove-pict` option to remove the picture file from the gallery folder.

Note: first item of the gallery will be shown on home page.How to run a job?
-----------------

### New alignment mode

Launch a new mapping between two fasta files and dot plot it.

{% if mode == "webserver" %}
![illustrating](/static/images/D-GENIES-run_na.png)
{% else %}
![illustrating](/static/images/D-GENIES-run-standalone_na.png)
{% endif %}

{% set puce=1 %}

#### ({{puce}}) Main menu

You just need to click on the main menu *run* tab, and fill the fields. All results will be stored in the result menu.

{% set puce=puce+1 %}

#### ({{puce}}) Updatable job name

Required field

A unique job name is set automatically. You can change it. Note that if a job already exists with the same name, it will be automatically renamed.

{% set puce=puce+1 %}

{% if mode == "webserver" %}

#### ({{puce}}) User email

Required field

Please enter your email. When the job is finished, you will receive a mai. Some features of the result page will also send you a mail to this address (see [manual](/documentation/result)).

{% set puce=puce+1 %}

{% endif %}

#### ({{puce}}) Target fasta

Required field

With the selector at the left, you can choose to select a local file or enter an URL. For a local file, click on the button at the right to select it.

Files must be in fasta format. We recommend using gzipped files to preserve bandwidth and speed job submission.

Allowed extensions: fa, fasta, fna, fa.gz, fasta.gz, fna.gz

Max file size: {{size}} ({{size_unc}} once uncompressed, {{size_ava}} in all-vs-all mode)

{% set puce=puce+1 %}

#### ({{puce}}) Query fasta

Optional field

Works like the target fasta. If not given, target file will be mapped to itself, in all-vs-all mode.

Max file size: {{size}} ({{size_unc}} once uncompressed)

{% set puce=puce+1 %}

#### ({{puce}}) Aligner

You can choose aligner to use for mapping. By default, it's minimap2.

If your job fails due to memory limit reached, you can try mashmap. It uses less resources. But is only suitable for highly similar genomes as it only detect matches with more than 75% of identity.

{% set puce=puce+1 %}

### Plot alignment mode

Dot plot an existing alignment file.

{% if mode == "webserver" %}
![illustrating](/static/images/D-GENIES-run_pa.png)
{% else %}
![illustrating](/static/images/D-GENIES-run-standalone_pa.png)
{% endif %}

{% if mode == "webserver" %}
{% set puce=4 %}
{% else %}
{% set puce=3 %}
{% endif %}

For numbers from 1 to {% if mode == "webserver" %}3{% else %}2{% endif %}, see previous section.

#### ({{puce}}) Alignment file

Required field (except if backup file is filled, see bellow)

An alignment file in PAF or MAF format.

Allowed extensions: paf, maf

With the selector at the left, you can choose to select a local file or enter an URL. For a local file, click on the button at the right to select it.

{% set puce=puce+1 %}

#### ({{puce}}) Target file

Required field (except if backup file is filled, see bellow)

Can be a fasta file or the corresponding index file.

To improve bandwidth and computation time, we recommend to use the index file. This file format is described [here](/documentation/formats#index-file). You can use [our tool](https://raw.githubusercontent.com/genotoul-bioinfo/dgenies/v{{version}}/src/dgenies/bin/index.py) to build it.

Allowed extensions:  
Fasta: fa, fasta, fna, fa.gz, fasta.gz, fna.gz  
Index: idx

With the selector at the left, you can choose to select a local file or enter an URL. For a local file, click on the button at the right to select it.

{% set puce=puce+1 %}

#### ({{puce}}) Query file

Optional field

Can be the fasta file or the corresponding index file.

Works like the target file.

{% set puce=puce+1 %}

#### ({{puce}}) Backup file

Optional field

If you downloaded the backup file from a previous job, you can enter it here to restore the dot plot. In this case, don't fill previous fields, only this one is required.

With the selector at the left, you can choose to select a local file or enter an URL. For a local file, click on the button at the right to select it.# Dot plot

In bioinformatics a dot plot is a graphical method that allows the comparison of two biological sequences and identify regions of close similarity between them. It is a type of recurrence plot.

More details of dot plot [here](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)). Below, some examples of events which can be detected by dot plots.

## Match

When two samples sequence are identical, it's a match.

![match](/static/images/dotplot/match.png)

## Gap

Dot plots can be used to detect a gap between two samples: small sequence which exists only in one sample, between two matching regions.

![gap](/static/images/dotplot/gap.png)

## Inversion

Sequence which exists in the two samples but not in the same order.

![inversion](/static/images/dotplot/inversion.png)

## Repeats

Dot plot can be used to detect repeated regions: a sequence which is repeated several times in a sample.

![repeats](/static/images/dotplot/repeat2.png)# Description of file formats used in D-Genies

##PAF (Pairwise mApping Format)

Default output format of minimap2.

It's a tabulated file. Description of columns below.

Col | Type   | Description                               
--- | ------ | -----------------------------------------
1   | string | Query sequence name                       
2   | int    | Query sequence length                     
3   | int    | Query start coordinate (0-based)                     
4   | int    | Query end coordinate (0-based)                       
5   | char   | ‘+’ if query/target on the same strand; ‘-’ if opposite               
6   | string | Target sequence name                      
7   | int    | Target sequence length                    
8   | int    | Target start on original strand (0-based) 
9   | int    | Target end on original strand (0-based)   
10  | int    | Number of matching bases in the mapping                 
11  | int    | Number bases, including gaps, in the mapping                    
12  | int    | Mapping quality (0-255; 255 for missing)  

Column 11 gives the total number of sequence matches, mismatches and gaps in the alignment; column 10 divided by column 11 gives the BLAST-like alignment identity.

PAF may optionally have additional fields in the SAM-like typed key-value format. Minimap2 may output the following tags:
 
Tag | Type | Description
--- | ---- | -----------------------------------------
tp  | A	   | Type of aln: P/primary, S/secondary and I,i/inversion
cm	| i    | Number of minimizers on the chain
s1	| i	   | Chaining score
s2	| i	   | Chaining score of the best secondary chain
NM	| i	   | Total number of mismatches and gaps in the alignment
AS	| i	   | DP alignment score
ms	| i	   | DP score of the max scoring segment in the alignment
nn	| i	   | Number of ambiguous bases in the alignment

Source: [minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html).


## Maf (Multiple Alignment File)

Description of the format is available [here](http://www.bx.psu.edu/~dcking/man/maf.xhtml).

## Index file

Index files used in D-Genies are built as follow.

First line contains the name of the sample. Next lines describes contigs of the sample. They are composed of two columns, tab separated. First it the name of the contig, second it's size in bases.

Example:

    Homo sapiens  
    chr1    248956422  
    chr2    242193529  
    chr3    198295559

## Backup file

Backup file is a TAR archive. It contains three files:

* The alignment file, in paf format, named `map.paf`.
* The target index, named `target.idx`.
* The query index, named `query.idx`.

Names of files must be kept. Otherwise, the backup file will not be accepted by the run form.Install your own instance
=========================

{% if version != "" %}
Latest available version: **{{version}}**
{% endif %}

Linux
-----

### Install

Install in 1 step:

As root:

    pip3 install dgenies

Or as simple user:

    pip3 install dgenies --user

Alternatively, you can install it manually:

    git clone https://github.com/genotoul-bioinfo/dgenies
    cd dgenies
    pip3 install -r requirements.txt
    python3 setup.py install

### Upgrade

#### Standalone mode

    pip3 install dgenies --upgrade
    
Add `--user` flag if you have not root access.

#### Webserver mode

    dgenies clear -c
    pip3 install dgenies --upgrade

Then, you need to restart your webserver.



### Requirements

D-Genies requires python >= 3.5, < 3.10 to run.

We use minimap2 (and mashmap2 on linux) for mapping. A binary of each of them is shipped with the program. Alternatively, you can
install your own from [minimap2](https://github.com/lh3/minimap2) and [mashmap2](https://github.com/marbl/MashMap/) repositories.

Some python modules are required (will be automatically installed by the commands above):

    Flask==1.0.*
    Jinja2~=2.11.3
    numpy
    requests~=2.20.1
    biopython>=1.70
    psutil~=5.6.6
    tendo==0.2.*
    matplotlib>=2.1.*
    intervaltree==2.1.*
    Markdown==2.6.*
    pyyaml~=5.4.1

Additional modules for webserver mode:

    Flask-Mail==0.9.*
    peewee==2.10.2
    python-crontab>=2.2.*

And if you use a cluster (webserver mode):

    drmaa==0.7.*

In webserver mode, you must install `mysqlclient` python module (will not be installed automatically) if you use mysql as RDBM.


Windows
-------

We provide an installer to install D-Genies. You can download it 
[here]({%if win32 %}{{win32}}{% else %}https://github.com/genotoul-bioinfo/dgenies/releases{% endif %}).

All requirements are present inside the package, so you don't have to do anything else.

### System requirements

You need Windows 7 or newer, 64 bits architecture.


How to start
-------------

You can launch DGenies in `standalone` mode or in `webserver` mode. You can use standalone mode if
you launch it locally and only one user launch jobs at once. If you are several users to use it
simultaneously or if you run it on a server, you must run it in webserver mode.

### Standalone mode

Unix: start with the command below:

    dgenies run

Optional arguments:

`-p <port>` run in a specified port (default: 5000)
`--no-browser` don't start the browser automatically

Windows: just click on the launcher in the desktop or into the install folder.

### Webserver mode

*Note: this mode is only available for Unix systems and will NOT work on MS Windows.*

#### Recommended method

Flask webserver (which is used in standalone mode) is not recommended in production servers.
So, we recommend using the WSGY module of Apache (or µWSGI + nginx, not documented here).

Once dgenies is installed, you just need to use the `/var/www/dgenies.wsgi` file into your apache
virtualhost file.

Here is an example of configuration file for apache:

    <VirtualHost *>
        ServerName <url>

        WSGIDaemonProcess dgenies user=<user> group=<group> threads=8
        WSGIScriptAlias / /var/www/dgenies/dgenies.wsgi

        <Directory /var/www/dgenies>
            WSGIProcessGroup dgenies
            WSGIApplicationGroup %{GLOBAL}
            Order deny,allow
            Allow from all
        </Directory>
    </VirtualHost>

With:
`<url>`: the URL of your instance
`<user>`: the user who launch the server
`<group>`: the group who launch the server

#### Debug method

For debug or for development only, you can launch dgenies through flask in webserver mode:

    dgenies run -m webserver

Optional parameters:

`-d` run in debug mode
`-o <IP>` specify the host into run the application (default: 127.0.0.1, set 0.0.0.0 for distant access)
`-p <port>` run in a specified port (default: 5000)
`--no-crons` don't run the crons automatically
`--no-browser` don't start the browser automatically (always true if *-d* option is given)



Running with a cluster
----------------------

If you want to run jobs on a cluster, some configuration is required. We only support SLURM and SGE schedulers. But please note that only SLURM scheduler has been fully tested.

Note: submitting jobs on a cluster is only available for webserver mode.

Jobs are submitted throw the DRMAA library. So you need it for your scheduler. Please see [configuration below](#cluster) to define path to this library.

Also, scripts for preparing data must be moved in a location accessible by all nodes of your cluster. You must move them in a same folder and set the full path to `all_prepare.py` script (full path on the nodes of the cluster) in the configuration file ([see below](#cluster)).

To get these scripts, follow the commands below:

    curl https://raw.githubusercontent.com/genotoul-bioinfo/dgenies/v{{version}}/get_cluster_scripts.py > get_cluster_scripts.py
    python get_cluster_scripts.py -d <dir>

With `<dir>`: the folder into save the scripts (must be accessible by cluster nodes).



Configuration
-------------

Changing the default configuration is not required for standalone mode, but you can want to custom some parts of the program.

Configuration file location:  
* Linux:  
    * `/etc/dgenies/application.properties` if installed with root access  
    * `~/.dgenies/application.properties` else  
* Windows:  
    * `application.properties` file of the install folder

The file is divided in 9 parts described below.

To change this file, please copy it into `application.properties.local` (at the same location) to avoid erase of the file on upgrades.

### Global

Main parameters are stored into this section:

* `config_dir`: where configuration file will be stored.
* `upload_folder`: where uploaded files will be stored.
* `data_folder`: where data files will be stored (PAF files and other files used for the dotplot).
* `threads_local`: number of threads to use for local jobs.
* `web_url`: public URL of your website.
* `max_upload_size`: max size allowed for query and target file (-1 to avoid the limit) - size uncompressed.
* `max_upload_size_ava`: max size allowed for target file for all-vs-all mode (only target given, -1 to avoid the limit) - size uncompressed.
* `max_upload_file_size`: max size of the uploaded size (real size of the file, compressed or not, -1 to avoid the limit).

For webserver mode only (ignored in standalone mode):

* `batch_system_type`: local for run all jobs locally, sge or slurm to use a cluster scheduler.

### Debug

Some parameters for debug:

* `enable`: True to enable debug
* `log_dir`: folder into store debug logs

### Cluster

This section concerns only the webserver mode with *batch_system_type* not set to *local*.

* `drmaa_lib_path`: absolute path to the drmaa library. Required as we use the DRMAA library to submit jobs to the cluster.
* `native_specs`: how to set memory, time and number of CPU on the cluster (should be kept as default).

By default, small jobs are still launched locally even if *batch_system_type* is not set to *local*, and if not too much of these jobs are running or waiting. This limit can be customized:

* `max_run_local`: max number of jobs running locally (if this number is reached, future jobs will be submitted to the cluster regardless of their size). Set to 0 to run all jobs on the cluster.
* `max_wait_local`; max number of jobs waiting for a local run (if this number is reached, future jobs will be submitted to the cluster regardless of their size). Set to 0 to run all jobs on the cluster.

You can also customize the size from which jobs are submitted on the cluster. If only one of these limit is reached, the job is submitted on the cluster.

* `min_query_size`: minimum size for the query (uncompressed).
* `min_size_target`: minimum size for the target (uncompressed).

Other parameters:

* `prepare_script`: absolute path to the all_prepare.py script downloaded in the section [above](#running-with-a-cluster).
* `python3_exec`: path to python3 executable on the cluster.
* `memory`: max memory to reserve on the cluster.
* `memory_ava`: max memory to reserve on the cluster un all-vs-all mode (should be higher than memory).
* `threads`: number of threads for launching jobs on the cluster (must be a divider of the memory).

### Database

This section concerns only the webserver mode.

In webserver mode, we use a database to store jobs.

* `type`: sqlite or mysql. We recommend mysql for better performances.
* `url`: path to the sqlite file, or url to the mysql server (localhost if the mysql server is on the same machine).

If type is mysql, some other parameters must be filled:

* `port`: port to connect to mysql (3306 as default).
* `db`: name of the database.
* `user`: username to connect to the database.
* `password`: the associated password.

### Mail

This section concerns only the webserver mode.

At the end of the job, a mail is send to the user to advise him of the end of the job.

* `status`: mail to use for status mail.
* `reply`: mail to use as reply to.
* `org`: name of the organisation who send the mail.
* `send_mail_status`: True to send mail status, False else. Should be True in production.

### Cron

This section concerns only the webserver mode.

We use crons to launch the local scheduler who start jobs, and some script that clean old jobs.

* `clean_time`: time at which we launch the clean job (example: 1h00)
* `clean_freq`: frequency of the clean job execution

### Jobs

This section concerns only the webserver mode.

Several parameters for jobs:

* `run_local`: max number of concurrent jobs launched locally.
* `data_prepare`: max number of data prepare jobs launched locally.
* `max_concurrent_dl`: max number of concurrent upload of files allowed.

### Example

Here, you can fill example data. At least target is required to enable example data.

Fill for target and query the absolute local path of the file. This path will not be shown to the client. Only the file name will be shown.

If at least target is filled, a button "Load example" will be shown in the run form. Click on it will load example data in the form.

### Analytics

Set `enable_logging_runs` to True will enable storage of analytics data. It stores for each job creation date, user mail, size of query and target, and batch type.


Customize your installation
---------------------------

You can make easily make some changes to the application, described below. To do so, you must first clone the D-Genies repository (except if changes can be done in the current installation, see below):

    git clone https://github.com/genotoul-bioinfo/dgenies
    
The created folder is named the `D-Genies repository` in the text below.

Then, you make changes described below. When done, you can easily install the cusomized version with this command:

    pip install .
    
Add the `--upgrade` flag if you already have D-Genies installed.


### Add or change alignment tools

D-Genies uses minimap2 as default aligner, but you can add other tools, or replace minimap2 by another tool. You can also change minimap2 executable path.

If your tool needs a parser (see below), you must customize the installation (see above). For other changes, you can make them on your current installation.

Tools definition YAML file:  
* Linux:  
    * `/etc/dgenies/tools.yaml` if installed with root access  
    * `~/.dgenies/tools.yaml` else  
* Windows:  
    * `tools.yaml` file of the install folder
    
To change this file, please copy it into `tools.yaml.local` (at the same location) to avoid to erase the file on upgrades.

Note: you can also edit the `tools.yaml` file of the D-Genies repository, which will be installed then (if you customize the installation). Edit it without renaming it.

For each tool, you can or must define several properties described bellow.

#### Executable path

Required.

Always set tool binary path in the `exec` property. If it differs for cluster execution, set it in the `exec_cluster` property.

#### Command line skeleton

Required.

The command line skeleton defines how to launch the program, in the `command_line` property. The skeleton must include the following tags:
  
* `{exe}`: will be replaced by the executable path
* `{target}`: will be replaced by the target fasta file
* `{query}`: will be replaced by the query fasta file
* `{out}`: will be replaced by the output file

If the tool can be multithreaded, add the `{threads}` tag which will be replaced by the number of threads to use.

If your program is able to use all-vs-all mode (target versus itself), define the same skeleton in the `all_vs_all` property. All tag described hereover must be set except the `{query}` tag.

#### Threads

Required.

Defines how much threads to use for this tool. Set it in the `threads` property for local executions, and the `threads_cluster` property for cluster execution (if different from the local one).

#### Max memory

Optional.

If the tool requires less memory than defined in the configuration file, you can add the `max_memory` property. The unit is Gb.

#### Parser

Optional.

If the tool does not output PAF format, you must define a function in the `src/dgenies/lib/parsers.py` file of the D-Genies repository. This function get the input file as first parameter and outputs a valid PAF file (output file, second parameter). You must reference the function name in the `parser` property.

#### Split before

Optional. Default: False

For some tools (like minimap2), splitting the query on 10 Mb blocks improves performances (blocks are merged after mapping). To enable this for the tool, set the `split_before` property to True.

#### Help

Optional.

Define a message to show aside the tool name in the run form. Set it in the `help` property.

#### Order

Optional. Default: random.

Define in which order we show tools in the run form. Set it in the `order` property.


### Add new formats

In `Plot alignment` mode in run form ([see here](/documentation/run#plot-alignment-mode)), we propose by default only PAF and MAF formats. It's easy to add new formats.

Just define 2 functions:

* Add the first one in the `src/dgenies/lib/validators.py` file of the D-Genies repository. It takes only one argument: the input file. It checks if the file format is correct and returns True in this case, else it returns False. The function name must be the same as the expected input file extension.  
* Add the second one in the `src/dgenies/lib/parsers.py` file of the D-Genies repository. It takes two arguments: the input and the output file. It convert the input file into a valid PAF file. The function name must be same as the previous one.


Maintenance
-----------

The `dgenies` command can be used to do some maintenance staff.

**Clear all jobs:**

    dgenies clear -j [--max-age <age>]

`--max-age` (opt): set the max age of jobs to delete (default: 0, for all)

**Clear all log:**

    dgenies clear -l

**Clear crons (webserver mode):**

    dgenies clear -c
    
**Display message on run form:**

You can display a message at the top of the run form. It can be used to add extra informations for user, or for prevent him for problem or for a maintenance on your instance.

    dgenies inforun -m "message to display" -t [warn, critical, info, success]
    
`-m <message>`: message to display. Html allowed.  
`-t <type>`: type of message: warn (orange background), critical (red background), info (blue background) or success (green background).

Remove the message by:

    dgenies inforun -c

Gallery
-------

Note: gallery is only available in webserver mode.

To add a job to the gallery, copy illustrating picture file into the *gallery* folder inside the data folder (*~/.dgenies/data/gallery* as default, create it if not exists). Then use the *dgenies* command:

    dgenies gallery add -i <id_job> -n <name> -q <query_name> -t <target_name> -p <pict_filename>

With:
`id_job`: the name of the job
`name`: name of the job to show in the gallery
`query_name`: name of the query
`target_name`: name of the target
`pict_filename`: filename added in the gallery folder (without path)

You can also delete an item from the gallery:

    dgenies gallery del -i <id_job>

or:

    dgenies gallery del -n <name>

With `id_job` and `name` as described above. You can add the `--remove-pict` option to remove the picture file from the gallery folder.

Note: first item of the gallery will be shown on home page.How to use?
-----------

![illustrating](/static/images/D-GENIES-result.png)

The result page (see above) presents the dot plot following the fasta files sequence order*. The alignment matches are presented as colored lines on the graphical panel. The colors correspond to similarity values which have been binned in four groups (less than 25%, between 25 and 50%, between 50 and 75% and over 75% similarity).

The top and right margin of the graphical panel show the sequence names. Depending on the sequence and name length, the names will be fully or partially presented. In order to ease visualization, all sequences smaller than 0.2 percent of the total length are merged in a unique super-sequence for which the margin is grayed. The left and bottom margins show the sequence size scales.

\* Except if more than 75% of the query length if composed of contigs with size less than 1% of the query length. In this case, contigs of the query are already sorted according to the reference.


### (1) Main menu

You can access to your previous results by clicking on the main menu in the `Results` item. A list will appear with all jobs you have launched on the server.


### (2) Select query and target

You can zoom on the graph by push the CRTL key while turning the mouse wheel forward to zoom in and backward to zoom out. Drag and drop to move the graph (with CTRL key still pushed).

You can zoom to a specific zone on the graph by clicking on the dotplot. This will zoom into the associated rectangle formed by the query-target contigs association.

Or, you can select the zone by selecting a query and a target in the dropdown menus at the top, and click Apply.

To come back to the initial view, click on the icon at the top right of the dotplot or press ESC.

### (3) Matches details

For each match, you can view positions on query and on target and the precise identity value by placing mouse cursor over it.

Due to technical limits, it doesn't work for too small matches. You can zoom to make it work for them.

### (4) Export

Several export options are available:

* Export as image (SVG or PNG - suitable for publication). All changes applied to the dot plot will be kept on export, including zoom.
* Download the PAF file generated by minimap2.
* Download the association table: TSV file with, for each contig of the query, the associated chromosome of the target with position of the match.
* List of contigs of the query which have no match with any chromosome of the target.
* List of chromosomes of the target which have no match with any contig of the query.
* Export the job as a tar file which can be re-uploaded in the run form to restore the job ([see here](/documentation/run#plot-alignment-mode)).

And, if you sorted the dotplot:

* Download the Fasta query file with contigs in the same order as in the dotplot.
* Download all contigs of the query assembled like the chromosomes of the target. We take the diagonal match line, and for all contigs that match the same chromosome, we stick them together, separated by a 100-N block.

### (5) Color scheme

You can change the default color scheme. Five other color schemes are available:

* Colorblind colors: colors more distinguishable for colorblind people.
* Black & White: for black & white printing.
* Reverse default: default colors but colors with the lower mapping identity have the lower brightness.
* Reverse Coloblind: same as reverse default, but with colorblind colors.
* All black: all in black.

To change color scheme, click on the legend.

### (6) Match size filtering

You can remove too small matches by moving the slider. By increments, it remove matches with size of 0.001 to 0.2% of the dotplot width. Too small matches are also removed by the `Remove noise` button (see below).

### (7) Match identity filtering

Set the minimal identity to show. All matches with a lower identity value will be hidden.

### (8) Strong precision

When checked, the strong precision check-box reduces match borders removing small matches from the graphical panel, often showing gaps between non contiguous matches.

### (9) Line breadth

Change the match lines thickness with the slider.

### (10) Chrom. border breadth

Change the visibility of the chromosomes borders with the slider.

### (11) Sort

You can sort (or unsort) contigs by clicking on the button. Contigs of the query will be sorted according to the reference. It will take few seconds.

How it works? For each contig of the query we search the region which have the biggest matches with the target and store these coordinates. Then, we sort contigs by their associated coordinates.

### (12) Hide noise

To remove noise. A match is considered noise if its size is small and its size frequency is quite high. Therefore we group matches by size bins, the number of bins corresponds to one tenth of the number of alignments, the bins are scanned in increasing size order to find the most represented one and from this one the one corresponding to one percent of its count is searched. All the alignments in bins smaller in size than this one are considered noise. It will take few seconds.

### (13) Similarity summary

To ease dot plot comparison, clicking the summary button generates a bar graph presenting the reference similarity profile, meaning the sums of the projections of the matches on the reference per similarity category divided by the total reference length. This graph is produced after sorting the query along the reference, removing included matches and noise filtering; result not shown on the graphical panel. It gives a realistic view of the overall reference and query similarity which is often not very precisely measured through visual inspection.

### (14) Delete job

By clicking on the button, your job will be definitively removed on the server. Be careful, this operation can not be undone!---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Computer (please complete the following information):**
 - OS: [e.g. linux, windows]
 - Browser [e.g. chrome, firefox]
 - Version [e.g. 50]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
