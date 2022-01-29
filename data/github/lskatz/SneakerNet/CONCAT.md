#DISCLAIMER 

Use of this service is limited only to **non-sensitive and publicly available data**. Users must not use, share, or store any kind
of sensitive data like health status, provision or payment of healthcare, Personally Identifiable Information (PII) and/or 
Protected Health Information (PHI), etc. under **ANY** circumstance. 

Administrators for this service reserve the right to moderate all information used, shared, or stored with this service at any time.
Any user that cannot abide by this disclaimer and Code of Conduct  may be subject to action, up to and including revoking access to 
services. 
# SneakerNet

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02334/status.svg)](https://doi.org/10.21105/joss.02334)

## Synopsis

A pipeline for processing reads from a sequencing run. Currently supports Illumina or Ion Torrent,
but it can be expanded to other platforms.

    # Run SneakerNet on the example data
    SneakerNetPlugins.pl --numcpus 4 t/M00123-18-001-test

<p align='center'>
  <img src='./docs/images/overview.png' alt='SneakerNet workflow' width='400' />
</p>

### Main steps

This is the default workflow in v0.14
but there are other workflows available as described in
[PLUGINS.md](/docs/PLUGINS.md#workflows).
 
* [Parse sample entries](/docs/plugins/sn_detectContamination-kraken.pl.md) - create an input file `samples.tsv`
* [Read metrics](/docs/plugins/addReadMetrics.pl.md) - get raw read yields and raw read QC summary (CG-Pipeline)
* [Assembly](/docs/plugins/assembleAll.pl.md) - assemble each genome (Shovill/skesa)
* [MLST](/docs/plugins/sn_mlst.pl.md) - 7-gene MLST (_mlst_)
* [Run Kraken](/docs/plugins/sn_kraken.pl.md)
* [Contamination detection](/docs/plugins/sn_detectContamination-kraken.pl.md) - check that all reads come from one taxon for each genome (Kraken)
* [Contamination detection](/docs/plugins/sn_detectContamination-mlst.pl.md) - check that all seven MLST genes have only one instance in the genome as expected (ColorID)
* [Base balance](/docs/plugins/baseBalance.pl.md) - check that the ratio of A/T is approximately 1 and same with C/T
* [Antimicrobial resistance gene prediction](/docs/plugins/sn_staramr.pl.md) - detect genotype and predict phenotype (staramr)
* [Pass/fail](/docs/plugins/sn_passfail.pl.md) - list all genomes that have failed Q/C
* [Transfer Files](/docs/plugins/transferFilesToRemoteComputers.pl.md) - files are copied to a remote folder
* [HTML summary report](/docs/plugins/sn_report.pl.md)
* [Email](/docs/plugins/emailWhoever.pl.md) the report

### Quick start

1. Install and configure SneakerNet - [from source](docs/INSTALL.md) or [with a container](docs/CONTAINERS.md)
2. Make an input folder from your MiSeq run [docs/SneakerNetInput.md](docs/SneakerNetInput.md)
3. Run [`SneakerNetPlugins.pl`](docs/SneakerNetPlugins.pl.md) on the input folder.

## Installation

See [docs/INSTALL.md](docs/INSTALL.md)

_NOTE_: to ensure all dependencies are met, please follow
the dependencies section under the [installation document](docs/INSTALL.md).

### Container installation

SneakerNet has been containerized and is at [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).
For more information, please see our [containers documentation](docs/CONTAINERS.md).

Here is a summary of Docker commands, from the [containers documentation](docs/CONTAINERS.md).

    # Pull image
    docker pull lskatz/sneakernet:latest
    # Import data directly from the MiSeq machine, where $MISEQ is a raw run folder exported by the MiSeq machine
    # and $INDIR is the newly created SneakerNet input folder
    docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
    # Run SneakerNet on the $INDIR (SneakerNet formatted folder)
    docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR

## Workflow

### Creating a SneakerNet project directory

_For more information on a SneakerNet-style folder, see [docs/SneakerNetInput.md](docs/SneakerNetInput.md)_

SneakerNet requires a project directory that is in a certain format already.
To create the project, you can use `SneakerNet.roRun.pl`.  For example,

    SneakerNet.roRun.pl --createsamplesheet -o M1234-18-001-test miseq/working/directory

M01234-19-01-test is a project folder name, where it is dash-delimited and contains
machine name, year, ordinal, and optionally a name.
Fastq files must be in the format of `_R1_` instead of `_1` and `_R2_` instead of `_2` for this particular script to parse the files properly.

### Running SneakerNet

It is generally a good idea to edit a file `snok.txt` to configure the run further.
For more information on the workflow, see the configuration section in `INSTALL.md`.
For example,

    echo "emails = example@example.com, blah@example.com" > t/data/M00123-18-001/snok.txt
    echo "workflow = default" >> t/data/M00123-18-001/snok.txt

And then run SneakerNet like so (optionally following the log with `tail -f`):

    SneakerNetPlugins.pl --numcpus 8 t/data/M00123-18-001 > t/data/M00123-18-001/SneakerNet.log 2>&1 &
    tail -f t/data/M00123-18-001/SneakerNet.log

#### Containers

SneakerNet has been containerized and is at [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).
For more information, please see our [containers documentation](docs/CONTAINERS.md).

## Output

_For more information, please see [docs/SneakerNetOutput.md](docs/SneakerNetOutput.md)_

SneakerNet produces a subfolder `SneakerNet/` in your run directory.
It also emails a report. To view a sample report, please go to 
t/report.html
in this repository.

## Plugins

SneakerNet is based on plugins.  In this context, a plugin is an independent script
that can run an analysis on a run directory, accept standard inputs (e.g., `--help`),
and create standard output files.

For more details, see the [plugins readme](docs/PLUGINS.md).

### Plugins for developers

You too can develop for SneakerNet!  For more information, 
please look at the [readme for plugins](docs/PLUGINSDEV.md)
and the [contributing](CONTRIBUTING.md) doc.

## Further reading

Please see the docs subfolder for more specific documentation.

# Contributing

Thank you for contributing! This is a collaborative environment.

## Requesting Changes: 

To request a change, please file an issue and/or a pull request.

### Open an issue in the repository 
                
If you do not have specific language to submit but would like to suggest a change or have something 
addressed, you can open an issue in this repository.

### Submit a pull request 
               
If you would like to contribute, please start off a discussion in issues.

# Creating a Culture of Innovation 

We aspire to create a culture where people work joyfully, communicate openly about things that matter, and 
provide great services globally. We would like our team and communities (both government and private sector)
to reflect on diveristy of all kinds, not just the classes protected in law. Diversity fosters innovation. 
Diverse teams are creative teams. We need a diversity of perpective to create solutions for the challenges
we face. 
    
This is our code of conduct (adapted from [18F's Code of Conduct](https://github.com/18F/code-of-conduct)). 
We follow all Equal Employment Opportunity laws and we expect everyone we work with to adhere to the [GSA Anti-harrasment Policy](http://www.gsa.gov/portal/directive/d0/content/512516), even if they do not work for the Centers for Disease Control 
and Prevention or GSA. We expect every user to follow this code of conduct and the laws and policies mentioned above.
    
## Be Empowering 
  
Consider what you can do to encourage and support others. Make room for quieter voices to contribute. Offer support 
and enthusiasm for great ideas. Leverage the low cost of experimentation to support your colleagues' ideas, and take 
care to acknowledge the original source. Look for ways to contribute and collaborate, even in situations where you 
normally wouldn't. Share your knowledge and skills. Prioritize access for and input from those who are traditionally 
excluded from the civic process. 


## Rules of Behavior

 * I understand that I must complete security awareness and records management training annually in order to comply with the 
   latest security and records management policies. 
 * I understand that I must also follow the [Rules of Behavior for use of HHS Information Resources](http://www.hhs.gov/ocio/policy/hhs-rob.html) 
 * I understand that I must not use, share, or store any kind of sensitive data (health status, provision or payment 
   of healthcare, PII, etc.) under ANY circumstance. 
 * I will not knowingly conceal, falsify, or remove information. 
 * I understand that I can only use non-sensitive and/or publicly available data. 
 * I understand that all passwords I create to set up accounts need to comply with CDC's password policy. 
 * I understand that the stewards reserves the right to moderate all data at any time. 

## Boundaries 

Create boundaries to your own behavior and consider how you can create a safe space that helps prevent unacceptable
behavior by others. We can't list all instances of unacceptable behavior, but we can provide examples to help guide 
our community in thinking through how to respond when we experience these types of behavior, whether directed at 
ourselves or others 

If you are unsure if something is appropriate behavior, it probably is not. Each person we interact with can define
where the line is for them. Impact matters more than intent. Ensuring that your behavior does not have a negative impact
is your responsibility. Problems usually arise when we assume that our way of thinking or behavior is the norm for everyone. 
    

### Here are some examples of unacceptable behavior: 
 * Negative or offensive remarks based on the protected classes as listed in the GSA Anti-harrasment Policy of
   race, religion, color, sex, national origin, age, disability, genetric information, sexual orientation, gender
   identity, parental status, maritual status, and political affiliation as well as gender expression, mental 
   illness, socioeconomic status or backgrounds, neuro(a)typicality, physical appearance, body size, or clothing. 
   Consider that calling attention to differences can feel alienating. 
 * Sustained disruption of meetings, talks, or discussions, including chatrooms. 
 * Patronizing language or behavior.
 * Aggresive behavior, such as unconstructive criticism, providing correction that do not improve the conversation 
   (sometimes referred to as "well actually's), repeatedly interrupting or talking over someone else, feigning 
    surprise at someone's lack of knowledge or awareness about a topic, or subtle prejudice. 
 * Referring to people in a way that misidentifies their gender and/or rejects the validity of their gender 
   identity; for instance by using incorrect pronouns or forms of address (misgendering) 
 * Retaliating against anyone who files a formal complaint that someone has violated these codes or laws. 
 
## Background 
CDC Scientific Clearance is the process of obtaining approvals by appropriate CDC officials before a CDC information product
is released to the public or CDC's external public health partners. Information products that require formal clearance include 
print, electronic, or oral materials, that CDC employees author or co-author, whether published by CDC or outside CDC. 
CDC contractors developing content on behalf of CDC for the public or CDC's external public health partners are also required
to put their content through the formal clearance process. The collaborative functions related to the projects in the R&D lab include blogs, wikis, forums, bug tracking sites, source control and others deemed as necessary. 

For those individuals within the CDC, adherence to the following policies are required: 
* CDC ["Clearance of Information Products Disseminated Outside CDC for Public Use"](http://www.cdc.gov/maso/Policy/PublicUse.pdf) 
* [HHS "Ensuring the Quality of Information Disseminated by HHS agencies"](http://aspe.hhs.gov/infoquality) 

All collaborativer materials will be controlled by the rules contained within this document. This will allow for the real-time collaboration opportunities among CDC employees, CDC contractors and CDC public health partners. 


## Credit: 
This code of conduct was mainly adapted from [18F's Code of Conduct](https://github.com/18F/code-of-conduct) and the [CDC's Informatics Innovation Unit R&D Lab's code of conduct.](http://www.phiresearchlab.org/?page_id=1715)

## Relevant Legal Considerations: 

[Laws enforced by the Equal Employment Opportunity Commission](http://www.eeoc.gov/laws/statutes/index.cfm) 

[Types of discrimination prohibited by law](http://www.eeoc.gov/laws/types) 

[New and proposed regulations](http://www.eeoc.gov/laws/regulations/index.cfm) 


  
          
          
---
title: 'SneakerNet: A modular quality assurance and quality check workflow for primary genomic and metagenomic read data'
authors:
- affiliation: 1
  name: Taylor Griswold^[co-first author]
- affiliation: "1, 2"
  name: Curtis Kapsak^[co-first author]
  orcid: 0000-0002-8735-1190
- affiliation: 1
  name: Jessica C. Chen
- affiliation: 3
  name: Henk C. den Bakker
  orcid: 0000-0002-4086-1580
- affiliation: 1
  name: Grant Williams
- affiliation: "2, 4"
  name: Alyssa Kelley
- affiliation: "1, 5"
  name: Eshaw Vidyaprakash
- affiliation: "1, 3"
  name: Lee S. Katz^[corresponding author]
  orcid: 0000-0002-2533-9161
date: "29 April, 2020"
bibliography: paper.bib
tags:
- QA/QC
affiliations:
- index: 1
  name: Enteric Diseases Laboratory Branch (EDLB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 2
  name: Weems Design Studio, Inc., Suwanee, GA, USA
- index: 3
  name: Center for Food Safety, University of Georgia, Griffin, GA, USA
- index: 4
  name: Waterborne Disease Prevention Branch (WDPB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 5
  name: IHRC, Atlanta, GA, USA
  
---

# Summary

Laboratories that run Whole Genome Sequencing (WGS) produce a tremendous amount of data,
up to 10 gigabytes for some common instruments.
There is a need to standardize the quality assurance and quality control process (QA/QC).
Therefore we have created SneakerNet to automate the QA/QC for WGS.

# Statement of need

Receiving a set of primary data from whole genome sequencing or metagenomics sequencing has become commonplace and perhaps ubiquitous in bioinformatics.
For a laboratory that performs WGS several times a week, some automation is necessary for both consistency and high throughput.

There are very few published workflows for performing an analysis on primary sequence data that span the breadth of initial standardized QA/QC (e.g., sequence yields, contamination checks, and subtyping).
For example, the Pandoo pipeline can be given a set of genomes to run analyses: species inference, 7-gene multilocus sequence typing (MLST), resistance gene profile, plasmid profile, virulence profile, and raw read QC [@Pandoo].
The Nullarbor pipeline is similar to Pandoo, but focused on public health datasets [@Nullarbor].
Another example is the ASA3P pipeline that runs raw read trimming, assembly, annotation, taxonomic classification, MLST, antibiotic resistance detection, virulence factor detection, reference mapping, and single nucleotide polymorphism (SNP) detection [@Schwengers654319].
However, no existing "broad stroke" QA/QC pipelines seem to be focused on a plugins-based architecture for batches of unrelated bacterial sequences or for batches of bacteria from different species.
To that end, we have created SneakerNet.
The major design principles are centered around the ability to collaboratively design plugins.
With the plugins architecture, SneakerNet can dynamically change for current and future needs
with input from the bioinformatics and public health community.

# Implementation

## Plugin design

SneakerNet has a modular plugin design, where the main program calls each plugin in an ordered succession.
Each plugin, in turn, reads a set of genomes as input.
Each plugin accepts specific flagged parameters, such that
the main program can call each plugin in a standardized way.
Workflows are thus defined as a specified order of plugins. 
An example workflow order might be genome assembly, followed by MLST, and finalized with report generation.
At the time of this writing, 25 plugins are available.
These plugins are listed in the documentation in a summary table,
and each plugin has its own documentation page.

## Plugin development

The plugin system has drastically lowered the activation energy needed to develop a new step in a
SneakerNet workflow. Documentation has been provided on how to develop a new plugin,
and 'Hello World' plugins have been published in three different languages: Perl, Python, and Bash.
Because plugins are not tied to any specific language, SneakerNet collaborators do not have to be bound by any specific language.

## Configuration

SneakerNet is highly configurable as described in the installation documentation.
There are many configurations.
We would like to highlight some ways that SneakerNet can be configured.

For some genera, SneakerNet comes packaged with some recommended configurations (e.g., _Salmonella_ or _Legionella_),
and an example genus with all options commented.
These options include the minimum coverage needed for a sample to pass QC
and even some detailed options to help customize a taxon for a particular plugin such as the antimicrobial resistance plugin.
Therefore, a user could easily add a taxon to customize the workflow for his or her instance of SneakerNet.
In fact, SneakerNet has been recently configured to accommodate the protist _Cryptosporidium_ successfully with input from the CDC WDPB Molecular Epidemiology Laboratory.

Users can also customize the workflow.
SneakerNet comes packaged with a default workflow which specifies the order of plugins that are run.
However, if a certain analysis is not needed, e.g., 7-gene MLST, then it can be removed from the configuration.
Likewise, if a new plugin is needed, it can be added into the workflow.

# Acknowledgements

The authors would like to thank the following SneakerNet users in the Enteric Diseases Laboratory Branch for useful feedback: Heather Carleton, Katie Dillon, Blake Dinsmore, Yang Gao, Jessica Halpin, Jasmine Hensley, Monica Im, Justin Kim, Charlotte Lane, Ana Lauer, Rebecca Lindsey, Tori McIntosh, Angela Poates, Zachary Rigney, Katie Roache, Ashley Sabol, Peyton Smith, Cheryl Tarr, Jenny Truong, Maryann Turnsek; Additionally, SneakerNet users from other branches: Dhwani Batra, Shatavia Morrison. The authors would also like to thank Aaron Petkau for useful feedback for staramr.

# References

# Example and tutorial for SneakerNet

This directory contains an example for how you would use SneakerNet. There is an 'inbox'
and a Rover Spreadsheet. First, this document will describe Rover; second, this
document will describe how to use SneakerNet.

## The Rover spreadsheet

Rover is a spreadsheet that can explore your reads and return basic information
for you. This is a double meaning between the NASA planatary vehicles and a dog
that might fetch for you! This example will not show you how to perform a 
sequencing run because the wet lab is outside of the scope of this project.
However, this example will show you how to fill in the appropriate sample
spreadsheet for the MiSeq, which SneakerNet will also read.

### How to use the Rover spreadsheet

There are two rover spreadsheet examples in this directory, and they differ
only on the number of samples.

* Rover_NexteraXTv2_16samples.xlsx
* Rover_NexteraXTv2_48samples.xlsx

Thank you to the National Enteric Reference Laboratory (NERL) for the basis of
this documentation and spreadsheets.

1. Ensure that the full header (i.e., "Run Name", "Sample Plate Name", "Sample 
Sheet Name", "Library Prep Date", and "Technician" fields) and 
sample-related fields for each sample (i.e., "Sample Name", "Organism", "Index
1", "Index 2", and "Genome Size Estimate") are filled in on the "Initial
Dilution" tab.

2. Optional: If raw read files are to be automatically transferred from the storage space to another location (e.g., Calculation Engine of BioNumerics), then enter "Yes" in the appropriate column for each sample to transfer on the "Raw Read Routing" tab.

3. Navigate to the "SampleSheet" tab of the Library Prep Workbook, and "Save
As..." a "CSV (Comma delimited) (*.csv)" file on a USB portable hard drive. 
  * Note: Because the workbook contains several tabs, a message might pop up along
the lines of "The selected file type does not support workbooks that contain 
multiple sheets..." Click OK to save only the SampleSheet tab.  
![Save as...](images/saveAs.jpg)
  * Excel will continue to ask you the OK/Cancel question about there being multiple sheets 
until you click Cancel. The file is saved where you specify when you click OK 
the first time, so hit Cancel when prompted the second time.  
![Warning](images/warning1.jpg)
  * Then, because the SampleSheet tab contains formulas, a second warning will pop up
  "Some features in your workbook might be lost... Do you want to keep using
  that format?" Click No, and the file will be auto-saved to your selected
  drive. (See warning above about second OK/Cancel prompt)  
![Warning](images/warning2.jpg)
4. On the MiSeq, open the Illumina Experiment Manager, and choose "Edit Sample Sheet"  
![EditSheet](images/MiSeqSampleEditSheet.jpg)
5. Navigate to your sample sheet csv file, and make sure all relevant fields are filled in, including the "Use Adapter Trimming" checkbox on the right side of the screen, and all indices used are present and compatible.

## Using SneakerNet

If you have installed SneakerNet correctly including editing the config files, 
then all you need to do is run it on the simulated example in the inbox folder.

Use `--help` to get SneakerNet help. You can also use `--help` on any SneakerNet plugin
to get help.

    $ SneakerNet.pl --help
    $ addReadMetrics.pl --help

To run SneakerNet on the example, use this syntax:

    $ SneakerNet.pl --inbox example/inbox --now

Using the `--now` flag assumes that no one can modify the run directory before you use it. Without 
this flag, SneakerNet monitors the inbox for two minutes to be sure that no one is 
currently adding a run, so that the run does not get corrupted. Therefore you should not 
use the `--now` flag in everyday use when others have access to the inbox.

### Requirements for a run

Currently, SneakerNet can only be used on Illumina-based run directories. The
directory name must be in the form of `machineName-year-runNumber-comment`
where machine name is a custom name, e.g., M0347 for a MiSeq. The year is a 
2 or 4-digit year. The run number is an integer. The comment is optional
but can be used to help indicate the date of a run or whatever else.
Collectively, `machineName-year-runNumber-comment` is the runId that might be
referred elsewhere in the documentation.

The following files must be present in a run directory before SneakerNet 
will accept it.

* QC information
  * `QC` - this is a folder that you create manually
  * `QC/runInfo.xml`
  * `QC/runParameters.xml`
  * `QC/CompletedJobInfo.xml`
  * `QC/GenerateFASTQRunStatistics.xml`
  * `QC/InterOp/` - this is a folder from an Illumina run with QC information. Keep all files in it intact.
* Data and metadata
  * `*.fastq.gz`
  * `SampleSheet.csv`, usually derived from Rover
  * `config.xml`

### Optional files for a run

At least one file at this time is optional

* `snok.txt`
  * Its presence helps SneakerNet understand that a run is ready to be read by SneakerNet. It can be present in the root directory of a full run instead of a folder with extracted contents.
  * It can also have key value pairs.  This file is still being developed and can accept the following key/value pairs, without ini-style headers.
    * `emails = email1@example.com, email2@example.com, ...`
    * `workflow = default`
      * can also be metagenomics or iontorrent

PLUGINS
=======

Please see the docs folder for more information

* [PLUGINS](../docs/PLUGINS.md)
* [DEVELOPMENT](../docs/PLUGINSDEV.md)

All files in this directory will be read as a configuration file. Although you can
try to use advanced syntax, the basic syntax is `key\tvalue`. These varables are
loaded into the scripts using the perl syntax `my $settings=readConfig();`.

**What type of PR is this?**

* Plugin addition
* Bug fix
* Other feature
* Other

**Describe the change**

**Does it change the usage of SneakerNetPlugins.pl?**

**Did you add any unit testing? Please describe**

---
name: Plugin contribution
about: Contribute a new plugin
title: "[PLUGIN CONTRIBUTION]"
labels: ''
assignees: ''

---

DOCUMENTATION
============

How to develop a plugin - (PLUGINSDEV.md)[https://github.com/lskatz/SneakerNet/blob/master/docs/PLUGINSDEV.md]

What is a plugin - (PLUGINS.md)[https://github.com/lskatz/SneakerNet/blob/master/docs/PLUGINSDEV.md]

CHECKLIST
=======

* [ ] First positional argument is the SneakerNet-formatted directory
* [ ] All flags are accepted in PLUGINS.md
* Soft-coded variables
  * [ ] Added to a conf file in `settings/`
  * [ ] Documented in the plugin usage
* [ ] Sample information is read from `samples.tsv`
* Outputs
  * [ ] File outputs are in `$rundirectory/SneakerNet/pluginName`
  * [ ] Any files to be attached to the SneakerNet report are copied into `$rundirectory/SneakerNet/forEmail`
* properties recorded after the plugin is run each time, in `$rundirectory/SneakerNet/properties.tsv`.
  * [ ] version is added to properties.tsv
  * [ ] any table paths are added to properties.tsv
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

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

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# SneakerNet Input folder

## Scripts

Usually you can run this script to import a MiSeq run:

    SneakerNet.roRun.pl path/to/dir -o SneakerNetDir --createsamplesheet

where `SneakerNetDir` is a properly formatted directory for SneakerNet.
For this particular script, fastq filenames must be in the format of
`_R1_`/`_R2_` instead of `_1`/`_2`.

## Building from scratch

Or you can just build the directory without an automated script.

### Required

#### \*.fastq.gz

Raw read files, both R1 and R2 if available. Single end raw reads are okay too.

#### SampleSheet.csv and/or samples.tsv

See below for samples.tsv, if SampleSheet.csv is not found.

### Recommended

#### samples.tsv

This file can be created from SampleSheet.csv by the plugin [sn_parseSampleSheet.pl](../docs/plugins/sn_parseSampleSheet.pl.md), from `SneakerNet.roRun.pl --createsamplesheet`, from `sn_createSampleSheet.pl`, or manually.
It has three columns, separated by tab: 

* sample - Name for a sample
* info - semicolon-delimited information. Keys are separated by an equals sign. For example: `taxon=Vibrio;route=calcEngine` means that the taxon is Vibrio and the ultimate destination for the files are in calcEngine. Vibrio is the key for the taxon in `config/taxonProperties.conf` and will not have meaning if there is no entry.
* filenames - semicolon-delimited fastq files. For example, `./2018AW-0585_S318_L001_R1_001.fastq.gz;./2018AW-0585_S318_L001_R2_001.fastq.gz`

#### snok.txt

This file stands for "SneakerNet Okay" which is an indicator that the run is ready to be analyzed.
The file is only necessary if a custom script needs to know whether the run is ready to be analyzed in an automated fashion.
At minimum, `snok.txt` is a zero-byte file.

However, this file can also have content in a key/value format separated by an equals sign.
The values can be comma-separated. Whitespace around each value is ignored.  For example:

    emails   = nobody@gatech.edu, nobody@cdc.gov, noreply@example.com
    workflow = default

The current keys are `emails` and `workflow`.
Emails describe where to send reports to, in addition to those in `config/emails.conf`.
Workflow describes the set of plugins, and their order.
Workflows are further described in [config/plugins.conf](PLUGINS.md).

### Optional

These files are found in the Illumina run directory although they are optional.

* CompletedJobInfo.xml
* config.xml
* QC/CompletedJobInfo.xml
* QC/GenerateFASTQRunStatistics.xml
* QC/RunInfo.xml
* QC/runParameters.xml
* QC/InterOp/ - this directory contains other files from the default Illumina run directory

# Advanced installation

This is an installation method to create a whole system centered
around cronjobs, logs, etc.
Only follow this method if you really know what you are doing.

## As root

For these steps, log in as root.

    $ sudo su

### Create a dedicated user

There should be a special user that has permissions for sequencing runs.  Ours is 
'sequencermaster'. All files made by sequencermaster should have privileges for 
both group and user. In this way, if you want to give special privileges to other
users for sequencing runs, you will add them to the sequencermaster group.
    
    $ useradd -m -s bash sequencermaster
    $ passwd sequencermaster

### Create a log file

    $ logfile=/var/log/SneakerNet.log
    $ touch $logfile
    $ chown sequencermaster.sequencermaster $logfile
    $ chmod 664 $logfile

### Take care of the log file with logrotate
    
Cat the following contents into a logrotate config file and use ctrl-D to signify the end of the file.  This config file will key into the logrotate command so that you do not have to worry about maintaining the sneakernet log.

**Tips:** The name of the file does not matter except it must be in the `/etc/logrotate.d` directory.  A good tutorial is at http://www.thegeekstuff.com/2010/07/logrotate-examples/.

    $ cat > /etc/logrotate.d/sneakernet
      /var/log/SneakerNet.log {
          missingok
          notifempty
          compress
          delaycompress
          copytruncate
          size 1M
          create 0644 sequencermaster sequencermaster
          dateext
          monthly
          maxage 9999999999999
          compressext .gz
      }
      ^D


### Create an 'inbox'

The inbox is some directory where sequencing runs will be deposited. In our lab, this
is a samba folder, but it could simply be a directory that user `cp` files to. The
inbox must have permissions like so.  In our example, we have a special group name
for this inbox so that other users can contribute to it. However, the group name
`sequencermaster` would also work.

    $ inbox=/path/to/inbox
    $ chmod g+s $inbox                             
    $ chown sequencermaster.edlb $inbox 
    $ ls -ld $inbox
      #  drwxrwsr-x. 4 sequencermaster edlb 4 Apr 18 12:03 /path/to/inbox


## As 'sequencermaster'

For these steps, log in as sequencermaster.

    $ ssh sequencermaster@localhost

### Download the software

    $ mkdir ~/bin
    $ cd bin
    $ git clone https://github.com/lskatz/SneakerNet.git
    $ cd SneakerNet

There are also a couple of prerequisites that the sequencermaster needs to install:

* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* Multithreaded Perl (already installed on most computers)
* Kraken: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/

When ready with the prerequisites, run the Makefile

    $ make

### Set up the cron job (optional)

Log in again as sequencermaster and set up a cronjob to run once an hour. SneakerNet
will check for run directories once an hour. It will do nothing if there is no 
run directory present.

    $ crontab -e

Add this line to the crontab file, and then save it.

      1  *  *  *   *     ((date; /usr/bin/perl /home/sequencermaster/bin/SneakerNet.pl; echo;) >> /var/log/SneakerNet.log 2>&1)

### Configuration

You will need to edit some files for configuration before using SneakerNet.

    $ cp -r config.bak config
    $ cd config

#### emails

List any emails here, comma-separated. These emails will be sent reports by default for each
SneakerNet run.

#### taxonProperties

Each taxon that you might have sequences for is defined here. If not defined here, nothing bad
will happen though.  For each taxon, you can show the minimum quality and coverage thresholds;
the genome size in bp; the regular expression to match filename to taxon; and the destination
folder on the remote computer.

#### settings

This file has certain key/values and should be left alone, unless you are a developer.

# Documentation folder

The following documentation is available

| File                         | Description |
|------------------------------|-------------|
|CONTAINERS.md                 | How to install / use the SneakerNet container |
|INSTALL.advanced.md           | Some older installation notes; might have some esoteric details |
|INSTALL.md                    | Installation instructions |
|INSTALL.taxonProperties.md    | How to add a new organism so that SneakerNet recognizes it and analyzes it appropriately |
|PLUGINSDEV.md                 | How to develop your own plugin |
|PLUGINS.md                    | The plugins catalog |
|SneakerNetInput.md            | SneakerNet input description |
|SneakerNetOutput.md           | SneakerNet output description |

# SYNOPSIS

SneakerNet is a set of plugins. Each plugin runs a distinct analysis
on a MiSeq run and accepts specific flags when it is run.
Therefore, any given plugin can run independently of the others (aside
from any prerequisite files, e.g., genome assemblies from a previous
plugin).

Quick table of contents:
- [Workflows](#workflows)
- [Which plugins are running?](#workflows-and-their-plugins)
- [Plugins on the command line](#command-line)
- [Plugins catalog](#catalog)

# Workflows

SneakerNet workflows define a particular order for the plugins to run.
They help resolve dependencies like ensuring that genome assemblies
are present before analyzed or enforcing that a report is generated only
after all plugins have created their outputs.

Workflows are defined in [plugins.conf](../config.bak/plugins.conf).

The exact order of plugins for all workflows can be found by running the command `SneakerNet.checkdeps.pl --list`.

To make your own custuom workflow, edit the file under `config/plugins.conf`.
Plugins are run in the order specified for any given workflow.  For example:

    default = pluginA.pl, pluginB.pl, pluginZ.pl

In this example, in the default workflow, 
if `pluginB` has a dependency on `pluginZ`, you might want to change the order
so that `pluginB` runs last.

    default = pluginA.pl, pluginZ.pl, pluginB.pl

## Default

This workflow runs most plugins and assumes that you have some flavor
of Illumina (MiSeq, HiSeq, MiniSeq).

## Ion Torrent

This workflow runs plugins designed for ion torrent.

## Metagenomics

For metagenomics runs. 

## Assembly

For assembly-only runs (ie, only assemblies and not raw reads in the folder).

## sarscov2

For running the SARS-CoV-2 workflow. Plugin(s) are prefixed with `sn_sars_`.

# Workflows and their plugins

Workflows in SneakerNet are defined by which plugins are run and in which order.
You can run `SneakerNet.checkdeps.pl --list` to see which plugins are run and in which order, for any workflow.
This script pulls from `config/plugins.conf`. 
Below is the output for SneakerNet version 0.11.2.
By default, only `default` will be run in a SneakerNet analysis, but other workflows are available.
The pseudo-workflow `all` is an alphabetical listing of all available plugins.

    $ SneakerNet.checkdeps.pl --list
        all
        addReadMetrics.pl, assembleAll.pl, baseBalance.pl, emailWhoever.pl, sn_assemblyWorkflow_init.pl, sn_crypto_assembleAll.pl, sn_crypto_gp60.pl, sn_detectContamination-kraken.pl, sn_detectContamination-mlst.pl, sn_detectContamination.pl, sn_helloWorld.pl, sn_helloWorld.py, sn_helloWorld.sh, sn_immediateStatus.pl, sn_iontorrent_assembleAll.pl, sn_kraken-metagenomics.pl, sn_kraken.pl, sn_mlst.pl, sn_mlst-wg.pl, sn_parseSampleSheet.pl, sn_passfail.pl, sn_report.pl, sn_SalmID.pl, sn_saveFailedGenomes.pl, sn_staramr.pl, transferFilesToRemoteComputers.pl

        assembly
        sn_assemblyWorkflow_init.pl, sn_mlst.pl, sn_staramr.pl, sn_passfail.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_report.pl, emailWhoever.pl

        cryptosporidium
        sn_parseSampleSheet.pl, addReadMetrics.pl, sn_crypto_assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_passfail.pl, transferFilesToRemoteComputers.pl, emailWhoever.pl

        default
        sn_parseSampleSheet.pl, addReadMetrics.pl, assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_detectContamination-mlst.pl, baseBalance.pl, sn_staramr.pl, sn_passfail.pl, transferFilesToRemoteComputers.pl, sn_report.pl, emailWhoever.pl

        iontorrent
        addReadMetrics.pl, sn_iontorrent_assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_passfail.pl, sn_staramr.pl, transferFilesToRemoteComputers.pl, emailWhoever.pl

        metagenomics
        sn_parseSampleSheet.pl, addReadMetrics.pl, sn_kraken.pl, sn_kraken-metagenomics.pl, sn_passfail.pl, sn_report.pl

![Default workflow](/docs/images/defaultworkflow.png)

# Command line

Each plugin can accept the following options. The first positional
parameter must be the SneakerNet run.

|Flag|Default value|description|
|:---|:------------|:-----------|
|`--help`|         |generate a help menu|
|`--numcpus`|     1|Parallelization|
|`--debug`|        |generate more messages or any other debugging|
|`--tempdir`|automatically generated, e.g., with `File::Temp` or `mktemp`|Where temporary files are located|
|`--force`|        |This is loosely defined but can be used for many things like overwriting output files|
|`--version`|      |Print a version in the format of X.Y or X.Y.Z|
|`--citation`|     | Print a citation statement. | 
|`--check-dependencies`|     | check all executable dependencies. Print executable dependencies to stdout and version information to stderr. Run `SneakerNet.checkdeps.pl` to check dependencies on all plugins. | 


# Catalog

Except for the legacy plugins, all plugins are prefixed with `sn_`.
The plugins are not specific to any one language, although the majority
are in Perl.

Contributions are welcome for the following plugin documents.

| Plugin                         | description |
|:-------------------------------|:------------|
|[sn_SalmID.pl](plugins/sn_SalmID.pl.md)          | Salmonella subspecies identification |
|[sn_staramr.pl](plugins/sn_staramr.pl.md)          | staramr antimicrobial resistance determinant analysis |
|[sn_passfail.pl](plugins/sn_passfail.pl.md)      | Table of pass/fail for each sample   |
|[sn_iontorrent_assembleAll.pl](plugins/sn_iontorrent_assembleAll.pl.md)    | Assembly for ion torrent data        |
|[addReadMetrics.pl](plugins/addReadMetrics.pl.md)| Raw read metrics                     |
|[sn_helloWorld.pl](plugins/sn_helloWorld.pl.md)               | Example plugin in Perl               |
|[sn_helloWorld.sh](plugins/sn_helloWorld.sh.md)               | Example plugin in Bash               |
|[sn_helloWorld.py](plugins/sn_helloWorld.py.md)  | Example plugin in Python             |
|[baseBalance.pl](plugins/baseBalance.pl.md)      | Dividing all As by Ts and all Cs by Gs to see if we get a ratio of 1 for each|
|[sn_mlst.pl](plugins/sn_mlst.pl.md)              | Runs 7-gene MLST on assemblies       |
|[sn_mlst-wg.pl](plugins/sn_mlst-wg.pl.md)              | Runs whole-genome MLST on assemblies       |
|[transferFilesToRemoteComputers.pl](plugins/transferFilesToRemoteComputers.pl.md)|Transfers files to a remote computer |
|[sn_detectContamination.pl](plugins/sn_detectContamination.pl.md)       | Detects potential contamination by kmer counting|
|[emailWhoever.pl](plugins/emailWhoever.pl.md)                 | Emails all results                   |
|[sn_detectContamination-mlst.pl](plugins/sn_detectContamination-mlst.pl.md)  | Runs 7-gene MLST on raw reads, checking for abnormal number of alleles |
|[sn_iontorrent_parseSampleSheet.pl](plugins/sn_iontorrent_parseSampleSheet.pl.md)|Turns the sample sheet for ion torrent into SneakerNet format |
|[sn_immediateStatus.pl](plugins/sn_immediateStatus.pl.md)           | Emails an immediate report           |
|[guessTaxon.pl](plugins/guessTaxon.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample |
|[sn_kraken.pl](plugins/sn_kraken.pl.md)                   | Runs metagenomics classifier on raw reads, or on assemblies if reads are not present. No secondary analysis is performed by this exact plugin. E.g., `sn_detectContamination-kraken.pl`. |
|[sn_detectContamination-kraken.pl](plugins/sn_detectContamination-kraken.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample and list at most a single major contaminant |
|[sn_kraken-metagenomics.pl](plugins/sn_kraken-metagenomics.pl.md)                   | Analyzes kraken results for a metagenomics sample |
|[assembleAll.pl](plugins/assembleAll.pl.md)                  | Assembles Illumina data              |
|[sn_assemblyWorkflow_init.pl](plugins/sn_assemblyWorkflow_init.pl.md)                  | For workflows that only have assembly data. Initializes the workflow so that other plugins can function properly. |
|[sn_crypto_assembleAll.pl](plugins/sn_crypto_assembleAll.pl.md)                  | Assembles Illumina data for Cryptosporidium        |
|[sn_crypto_gp60.pl](plugins/sn_crypto_gp60.pl.md)                  | Provides the gp60 profile for Cryptosporidium        |
|[sn_parseSampleSheet.pl](plugins/sn_parseSampleSheet.pl.md)          | Turns the sample sheet for Illumina into SneakerNet format |
|[sn_report.pl](plugins/sn_report.pl.md)                    | Creates an HTML report from all other plugins |
|[sn_sarscov2_assembleAll.pl](plugins/sn_sarscov2_assembleAll.pl.md)              | Runs assembly for SARS-CoV-2 amplicon-based genomes |
|[sn_assembleAll_reference.pl.md](plugins/sn_assembleAll_reference.pl.md)              | Runs reference assembly |
|[sn_saveFailedGenomes.pl](plugins/sn_saveFailedGenomes.pl.md)                    | Saves genomes into the destination folder, into a QC_Fails subfolder|

# SYNOPSIS

These are instructions on how to make a SneakerNet plugin.  It does not
matter which language the plugin is coded in.  All that matters is that
it is executable by the SneakerNet user and that it can accept certain
parameters.

# EXAMPLES

Please see below for language specific tips for Hello World examples.

# How to make a plugin

Two major steps are described below for making a plugin.

## Create the script.
1. The first positional argument must be the MiSeq run directory
2. The script must accept the following flags with the following example 
   values (or no values). The flags must adhere to original meaning behind
   the flags. For example, `--help` must provide a help menu.
   The script does not necessarily need to _use_ these flags however.
3. Add any desired soft-coded variables such as the location of a blast database
   into `config.bak/settings` and `config/settings`.
   `config.bak` is tracked with `git`, and `config` is created upon
   SneakerNet installation.
4. Inputs: the plugin can and should read sample names and properties
   from `samples.tsv`.
5. If the plugin generates any files, please organize them into 
   `runDirectory/SneakerNet/customdirectory` (where `customdirectory`
   is a name of your choice), and add any results for the
   resulting email to `runDirectory/SneakerNet/forEmail`. Any files under
   this directory will be emailed with the SneakerNet email.
6. List your results in `runDirectory/SneakerNet/properties.tsv`. This will 
   ensure that your results will be in the email report.
   1. `properties.tsv` is three columns: plugin name, key, value
   2. You should include `version`
     
|Flag|Default value|description|
|:---|:------------|:-----------|
|`--help`|         |generate a help menu|
|`--numcpus`|     1|Parallelization|
|`--debug`|        |generate more messages or any other debugging|
|`--tempdir`|automatically generated, e.g., with `File::Temp` or `mktemp`|Where temporary files are located|
|`--force`|        |This is loosely defined but can be used for many things like overwriting output files|
|`--version`|      |Print a version in the format of X.Y or X.Y.Z|
|`--citation`|     | Print a citation statement to give yourself credit. | 
|`--check-dependencies`|     | check all executable dependencies. Print executable dependencies to stdout and version information to stderr. | 

## Available inputs

How do you specify where the input assemblies are?  Or where the kraken inputs are?
After all, the only parameter is the input directory right?

Never fear!  Here are some basic inputs that you have access to as a developer.
The most important input is `samples.tsv`.
Here are the most common inputs you might want to have access to:

|What    | path from root input folder [^1] | Common source script or plugin [^2] | For more information |
|--------|----------------------------------|-------------------------------------|----------------------|
|Sample information | `samples.tsv`         | `SneakerNet.roRun.pl --createsamplesheet` or `sn_createSampleSheet.pl` | [SneakerNetInput.md](/docs/SneakerNetInput.md#samplestsv) under the samples.tsv section |
|Assemblies and genbank files               | `SneakerNet/assemblies/[samplename]/` | [`assembleAll.pl`](/docs/plugins/assembleAll.pl.md) | |
|Assembly and predicted gene metrics        | `SneakerNet/forEmail/assemblyMetrics.tsv` | [`assembleAll.pl`](/docs/plugins/assembleAll.pl.md) | |
|Kraken output                              | `SneakerNet/kraken/[samplename]`    | [`sn_kraken.pl`](/docs/plugins/sn_kraken.pl.md)| |
|read metrics | `SneakerNet/forEmail/readMetrics.tsv` | [`addReadMetrics.pl`](/docs/plugins/addReadMetrics.pl.md)|

[^1] `[samplename]` is the name of the sample in `samples.tsv`.  
[^2] It is possible that other sources would generate this information and so these are the most common sources.

## Activate the script as a plugin

1. Place it in the SneakerNet.plugins folder
2. chmod the script to be executable
3. Add the plugin to the list of plugins in `config.bak/plugins` and `config/plugins` 

## Unit tests

Create a file under `t/` and copy a previous example unit test from that directory.

To run all unit tests:

    perl Makefile.PL
    make
    make test

## language-specific tips

There are three Hello World example plugins. For perl, there is
an API. However, I would kindly take any help in other languages.

### perl

Please see the perl hello world plugin. Also, please run `perldoc lib/perl5/SneakerNet.pm`
for help on the library API.

### python

Please see the python hello world plugin

### shell

Please see the hello world shell plugin
# Containers

We have containerized SneakerNet into a docker image that is publicly available on [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).

Some test data can be found in this repository in
t/M00123-18-001-test/
and is generally what the variable `INDIR` would hold in this document.

## Requirements
Docker, Singularity,  or another Docker-compatible container software must be installed e.g. shifter (untested)

## Database(s)

As of SneakerNet version 0.11.0, the Kraken database must be installed separately.
You will need to choose a directory to keep the database in.
An example path is given in the instructions below, with `KRAKEN_DEFAULT_DB`.
```bash
KRAKEN_DEFAULT_DB=$HOME/db/kraken1/minikraken_20171013_4GB
mkdir -pv $KRAKEN_DEFAULT_DB
pushd $KRAKEN_DEFAULT_DB
cd ..
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz
tar -zxvf minikraken_20171019_4GB.tgz
rm -vf minikraken_20171019_4GB.tgz
popd
```

## Singularity

### Singularity image installation
Check to make sure Singularity is installed:
```bash
singularity --help
```
If it is not installed, please visit the following link for installing the latest version of Singularity (at time of writing this). You will need administrator priveleges and to install some additional dependencies before installing singularity itself. [https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps](https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps)


(OPTIONAL) Make a fast scratch directory if your current drive is not fast.
This step, including the `SINGULARITY_TMPDIR` variable is optional.

    export SINGULARITY_TMPDIR=/scratch/$USER/singularity-tmp
    mkdir -pv $SINGULARITY_TMPDIR

Navigate into a directory where you would like your Singularity image to be stored.
Build the image with `singularity build`.

    cd your/containers/directory/
    singularity build sneakernet.sif docker://lskatz/sneakernet:latest

### Running SneakerNet using Singularity

Create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
The variable `MISEQ` represents the original MiSeq directory,
which is then converted into a new format in a new directory, `INDIR`.
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below, where file transfer and email is disabled.

    export MISEQ=my/miseq/run/dir
    export INDIR=my/run/dir

    # This line is here to catch you if you missed the KRAKEN_DEFAULT_DB step
    [[ -e "$KRAKEN_DEFAULT_DB" ]] || echo "ERROR: please set KRAKEN_DEFAULT_DB"
    
    # this assumes $MISEQ and $INDIR are in your $PWD
    singularity exec -B $PWD:/data sneakernet.sif SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
    singularity exec -B $PWD:/data -B $KRAKEN_DEFAULT_DB:/kraken-database sneakernet.sif SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR

## Docker

### Docker image installation

Docker CE must first be installed onto your system. Check that it is installed by running:
```bash
docker info
```
If it is not installed, visit: [https://docs.docker.com/install/](https://docs.docker.com/install/) and follow the install instructions according to your operating system.

For Ubuntu users, please see concise install instructions here: [https://github.com/StaPH-B/scripts/blob/master/image-information.md#docker-ce](https://github.com/StaPH-B/scripts/blob/master/image-information.md#docker-ce)

Download the docker image using the `latest` docker image tag
```bash
docker pull lskatz/sneakernet:latest
```

### Running SneakerNet using Docker

As in the Singularity section,
create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
The variable `MISEQ` represents the original MiSeq directory,
which is then converted into a new format in a new directory, `INDIR`.
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below, where file transfer and email is disabled.

```bash
export MISEQ=my/miseq/run/dir
export INDIR=my/sneakernet/run/dir

# This line is here to catch you if you missed the KRAKEN_DEFAULT_DB step
[[ -e "$KRAKEN_DEFAULT_DB" ]] || echo "ERROR: please set KRAKEN_DEFAULT_DB"

# -v flag will mount your PWD into the /data directory inside the container
# -u flag preserves your user/group when executing commands in the container
# --rm flag will remove/delete the container after it exits 
# make sure output files are written to /data so you don't lose them after the container exits!
docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR
```
# SneakerNet Output

SneakerNet creates a subdirectory `SneakerNet/` in the input folder.
In this document, this subfolder will be referred to as `SneakerNet/`.
This document describes those files.

## Files under `SneakerNet/`

### properties.txt

This file records special values from each plugin in the format of
plugin name, key, value.
Duplicate plugin name and key combinations are allowed; however,
the last plugin name and key combination overrides the earlier instances.
Therefore for example, a warning from a plugin will be "erased"
if the plugin is run again and it records a value describing zero warnings.
There are some special keys that are treated differently:

* table - the relative path to a `.tsv` file from the main folder that describes a summary table from a plugin, e.g., `SneakerNet/forEmail/kraken.tsv`
* version - the version number of the plugin that was run in the format of simple integers separated by periods, e.g., `2.3.1`
* warnings - a summary of warnings from the plugin, e.g., `5 samples do not have fastq files`. Use "0" for no warnings.
* errors - a summary of errors from the plugin, e.g., `Error reading summary table`. Use "0" for no errors.
* date - date of analysis in YYYY-MM-DD format
* time - time of analysis in HH:MM:SS format
* other - other keys are allowed but not not specifically described here, e.g., `resfinder_gene_drug_version` for the staramr plugin. For other special keys, please see the individual plugin documentation pages.

### forEmail/passfail.tsv

This file records sample names and categories on which
a given sample can pass or fail. Each value for a category
is a simple integer.

* 1 - This sample failed on this category
* 0 - This sample passed on this category
* -1 - It is unknown whether this sample passed. For example, if the taxon is unknown and therefore the coverage is unknown.

### Other files under `forEmail/`

Please see individual plugin documentation pages on other files

## Subfolders

There are some special subfolders in `SneakerNet/`.

### forEmail

`SneakerNet/forEmail/` is the directory that all email attatchments live.
Plugins will create summary files, usually tables, and then they will get attached to the report email.

### assemblies

`SneakerNet/assemblies` is the directory that all genome assemblies and gene predictions live.
For each sample, e.g., `sample1`, there is a subfolder, e.g., `SneakerNet/assemblies/sample1`.
There is at least one `.fasta` and `.gbk` file in each sample subfolder.
Each file is named after the basic assembly method, e.g., `SneakerNet/assemblies/sample1/sample1.skesa.fasta`.

### kraken

`SneakerNet/kraken` contains a subfolder for each sample name,
e.g., `SneakerNet/kraken/sample1`.
In each sample folder, there are
kraken output files, or output files compatible with Kraken output files.
These files are:

* kraken.report - a tab-separated file whose values have optional space padding. Its fields are
  * percent of reads
  * number of reads
  * number of reads specific to this exact taxon
  * rank - one letter representation, e.g., `S` for species
  * taxid
  * Name for the taxon
* kraken.taxonomy - a tab-separated file whose values are
  * number of reads specific to this exact taxon
  * taxonomy lineage separated by tabs, starting with `root` and ending with the genus/species name
* report.html - the Krona html file for visualization

`kraken.out` is _not_ always present due to space limitations but
would describe the taxonomic classification for each read.

### Other

Other plugins can create subfolders. For example, the plugin `baseBalance.pl`
creates a subfolder `SneakerNet/baseBalance`.
Please see individual plugin documentation pages for more information.


# Taxon Properties

Each organism has an entry in the file `config/taxonProperties.conf`.
This is useful for configuring each taxon.
For example, if your taxon is _Klebsiella_, you can add add a new entry 
`Klebsiella` into `taxonProperties.conf` and describe exactly the
thresholds you expect and all the properties of the taxon, such that
SneakerNet understands exactly how to treat this taxon.

To invoke the correct taxon, simply use it in `samples.tsv`.
See the taxon key in the info field in [SneakerNetInput.md](SneakerNetInput.md) for more information.

## Example

There is an example taxon (`SAMPLE_TAXON`) in `taxonProperties.conf` that you can follow.
Here is the entry as of version 0.9.9

    [SAMPLE_TAXON]
    # What coverage level you are setting as a threshold.
    # Below this coverage, the sample fails.
    coverage=30
    # If a taxon is not specified in samples.tsv, then
    # create a regex to see if we can guess what the sample is
    # based on the filename
    regex='^\d+q\-\d+'
    # What is the genome size in bp
    genomesize=5000000
    # What is the minimum average Q score you accept
    # before failing the sample
    quality=28
    # When the reads are transferred to a remote location
    # as specified in settings.conf, what subfolder do they
    # go to?
    dest_subfolder=Example
    # What subfolder do you have a wgMLST scheme in?
    wgMLST=Example
    # Which option to use for staramr for pointfinder
    pointfinder=example

## Specification

### Format

Each entry has to start with a header in brackets. In the example, the header is `[SAMPLE_TAXON]`.
The header corresponds to the entry in `samples.tsv`.

Comments are allowed per line.

Keys and values are separated by equals sign. Whitespace is allowed around the equals sign.

Values can be comma separated and the order is retained.

Further documentation on the format can be found in [`Config::Simple`](https://metacpan.org/pod/Config::Simple).

### Specific properties

| Key            | Description |
|:---------------|:------------|
|coverage        |The minimum coverage required for a sample to pass. Without this coverage, a genome will be marked as "failed" which will have implications such as not transferring to the destination in the plugin `transferFilesToRemoteComputers.pl`.|
|regex           | If a taxon is not specified, SneakerNet will use this regex (Regular expression pattern) on the filename(s) to try to guess which taxon it is. Must be specified but you can give some nonsense value to ensure it does not match anything. |
|genomesize      | The estimated genome size in bp. Will be used to estimate genome coverage.|
|quality         | Minimum average phred score accepted before a sample will fail. Any individual file can cause the sample to fail, even if all others pass.|
|dest_subfolder  | Subfolder that raw reads will be transferred to in the plugin `transferFilesToRemoteComputers.pl`. See `transfer_destination_string` in `settings.conf` for more information. |
|wgMLST          | A subfolder in `SneakerNet/db/wgMLST` that contains a chewBBACA-formatted MLST scheme.|
|pointfinder     | A value to give to staramr for the pointfinder option.|


# SYNOPSIS

`SneakerNetPlugins.pl`

Runs all plugins for a SneakerNet-formatted run

# Software requirements

See: [dependencies](/docs/INSTALL.md#dependencies)

# Algorithm

Reads the workflow from `snok.txt` in an input folder.
If one is not provided, will run the default workflow.
Workflows are described in [PLUGINS.md](/docs/PLUGINS.md).

Each workflow's plugins and the order of plugins
are determined by `plugins.conf`, documented in [PLUGINS.md](/docs/PLUGINS.md).

## Passing of flags

Some options are passed onto the plugin scripts

* `--force`
* `--numcpus`
* `--tempdir`

# Outputs

Virtually all workflows create:

* [`report.html`](/docs/plugins/sn_report.pl.md)
* an [email](/docs/plugins/emailWhoever.pl.md) with the report

This script also produces a variety of outputs described
in their own documentation [SneakerNetOutput.md](/docs/SneakerNetOutput.md).

# Installation

## Quick installation

    mkdir ~/bin
    cd bin
    git clone https://github.com/lskatz/SneakerNet.git
    cd SneakerNet
    
    # the following two lines are for local installations
    cpanm --local-lib=$HOME local::lib && eval $(perl -I $HOME/lib/perl5/ -Mlocal::lib)
    export PERL5LIB=$PERL5LIB:$HOME/lib/perl5:$HOME/lib/perl5/x86_64-linux-gnu-thread-multi:$HOME/lib/perl5/x86_64-linux-gnu-thread-multi/auto
    
    # The following lines are regardless of local or global installation
    cpanm version::vpp
    cpanm --installdeps --notest . -l .
    perl Makefile.PL
    make

## Test the installation

    make test

## Dependencies

## Resources
The minimum requirements are about 4G RAM, but 8G RAM is recommended.
An average run with only 1 thread is about 200 minutes; however with 4 threads, the run time gets down to about 9 minutes.
For more information: https://github.com/lskatz/SneakerNet/issues/42#issuecomment-672900992

## Software
There are few core dependencies. However, SneakerNet is
based on plugins that have individual dependencies.
Here is a list of dependencies across the board, although
it is possible that all of the following are not needed
for a complete SneakerNet analysis.

To find the right dependencies for you, run the following

    SneakerNet.checkdeps.pl --list  # and find the workflow, e.g., 'default'
    SneakerNet.checkdeps.pl default # Check deps for all plugins for the workflow 'default'

#### Comprehensive list

This list was created using `SneakerNet.checkdeps.pl [iontorrent metagenomics cryptosporidium default]`
on version 0.8.14.
Dependencies may or may not have changed since then but they can still be checked using this script.

All workflows require perl v5.12 or higher, compiled with multithreading.
This version of perl is already installed in most modern Linux operating systems.

##### Default workflow

* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* Python3
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* sendmail
* Skesa
* staramr
* `zip`

##### metagenomics

* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* `zip`

##### cryptosporidium

* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* Shovill
* `countGP60repeats.pl`: currently in development in a private repo. To exclude, remove `sn_crypto_gp60.pl` from `config/plugins.conf` (already not included by default)
* `zip`

##### iontorrent

* Shovill
* SPAdes
* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* Python3
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* sendmail
* Skesa
* staramr
* `zip`

##### Non-workflow

###### SneakerNet.roRun.pl

* bcl2fastq

## Configuration

You will need to edit some files for configuration before using SneakerNet.

    $ cp -r config.bak config
    $ cd config

Some settings are necessary to change in the config/\*.conf files.
Please edit these files accordingly.

Special variables are allowed and have a __VAR__ format.
These are the variables currently allowed

* `__BASE__` evaluates to the installation directory of SneakerNet. Useful for file paths, e.g., location of where a bed file is.

### emails

List any emails here, comma-separated. These emails will be sent reports by default for each
SneakerNet run.

### taxonProperties

_see [INSTALL.taxonProperties.md](INSTALL.taxonProperties.md) for more information_

Each taxon that you might have sequences for is defined here. If not defined here, nothing bad
will happen though.  For each taxon, you can show the minimum quality and coverage thresholds;
the genome size in bp; the regular expression to match filename to taxon; and the destination
folder on the remote computer.

### settings

This file has certain key/values and should be carefully edited.
Some individual plugin settings are found here.

### plugins

This file defines workflows.
The first keyword is the workflow. After the equals sign,
a number of plugins are listed and delimited by commas.
They are shown in order of execution. For example, all
assembly-based plugins must appear after the genome
assembly plugin runs.

You are able to create your own workflow or use the 
pre-defined ones here.

Each workflow is listed in the SneakerNet `snok.txt` file
in the format of `workflow = default` where in this example
the default workflow is invoked. If no workflow is given
in `snok.txt` or if `snok.txt` is missing, the default
workflow will be used.

# SYNOPSIS

`sn_saveFailedGenomes.pl`

Sends semi-failed genome reads to a QC_Fails subfolder

# Software requirements

* rsync

# Algorithm

If a genome has >= 10x coverage but less than the required coverage,
sends the genome reads to a subfolder, e.g., /mnt/CalculationEngineReads.test/Campy/QC_Fails/

# Outputs

none

# SYNOPSIS

`sn_immediateStatus.pl`

Emails an immediate report on whether a run has warnings

# Software requirements

# Algorithm

This plugin is meant to be run as one of the first plugins in a workflow.
It reads `samples.tsv` and the fastq files in the directory and determines:

1. Does every sample have fastq file(s)?
2. Is every fastq file associated with a sample?

# Outputs

## Table

Table of the immediate reaction, with columns

* ErrType
  * fastq - A fastq file was found, but no sample "owns" it.
  * sample - A sample was listed in `samples.tsv` but no fastq files were found.
* Sample - either the sample or the fastq file in question
* ErrKeyword - A shorthand way of writing out the Error column
* Error - A specific error message. It might offer help such as suggesting the correct fastq files or the correct sample.

## email

An email with the table is sent to those listed in `snok.txt` and `emails.conf`.

# SYNOPSIS

`sn_iontorrent_parseSampleSheet.pl.md`

Parses an ion torrent run for the sample sheet.
Sets up sample sheet appropriately for downstream
plugins.

# Software requirements

# Algorithm

`TODO`

# Outputs

`samples.tsv`

# SYNOPSIS

`sn_passfail.pl`

Marks each sample as passing or failing

# Software requirements

# Algorithm

Reads upstream plugins for pass/fail criteria.
Currently these criteria are only coverage and
read quality, from the read metrics plugin.

# Outputs

Table with columns

* sample
* assembly - if all sub quality checks pass, then this passes. If all are unknown, then this is unknown. If any fail, then this fails.
  * number of Ns (threshold defined in `taxonProperties.conf`) or number of genes
  * longest contiguous contig - longest tract of genome between two Ns
* coverage - pass or fail
* quality  - pass or fail
* kraken - pass or fail 

1 indicates failure; 0 indicates pass; -1 indicates unknown (e.g., if Kraken was never run)

# SYNOPSIS

`sn_detectContamination-kraken.pl`

Attempts to guess a taxon for each sample and then lists
at most one major contaminant for the sample.

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

## species ID

Reads Kraken results from raw reads or from assemblies.
Starting with species, if at least 25% of the reads correspond
with one species taxon, then records the majority species for
the sample.  If 25% are not captured at species, then it moves
onto genus, etc.

## contamination ID

For the rank that is identified (species, genus, etc),
if at least 1% of reads conflict with the species identified,
lists the major conflicting taxon.

# Outputs

Table with columns

* sample
* assumed taxon
* best-fitting taxon
* percent of reads that match the best-fitting taxon
* rank (S is species, G is genus, F is family, ...)
* major conflicting taxon
* percent of reads supporting the major conflicting taxon

# SYNOPSIS

`sn_assemblyWorkflow_init.pl`

Initializes an assembly-only workflow so that other plugins can read the proper data.

# Software requirements

* CG-Pipeline
* Prodigal

# Algorithm

Copies each assembly in root SneakerNet-formatted folder into `$dir/SneakerNet/assembly/$sample`.
Runs assembly metrics.
Runs prodigal on all assemblies.
Runs gene prediction metrics.

# Outputs

Table with columns

* sample
* genomeLength
* CDS - number of coding sequences
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC


# SYNOPSIS

`sn_kraken-metagenomics.pl`

Looks at Kraken results in the context of a metagenome

# Software requirements

* Krona

# Algorithm

Merges individual Krona plots into a single html file.

TODO: create table of top X taxa for each sample

# Outputs

HTML file of Krona plots

# SYNOPSIS

`sn_SalmID.pl`

Attempts to guess a taxon for each sample

# Software requirements

* SalmID

# Algorithm

Runs SalmID on each sample to determine Salmonella
subspecies estimates. Ideally, one isolate would have
100% of one subspecies or another.

# Outputs

Table with columns

* `TODO`

# SYNOPSIS

`emailWhoever.pl`

Emails all SneakerNet results

# Software requirements

* sendmail (already installed on most systems)

# Algorithm

Attaches all files under `SneakerNet/forEmail` and emails
them to the list under `config/emails.conf`. If any emails
are found in `snok.txt` or `SampleSheet.csv`, then they 
will also receive the email.

# Outputs

email with attachments

# SYNOPSIS

`sn_mlst`

Runs 7-gene MLST on each sample

# Software requirements

* Torsten _mlst_

# Algorithm

For each genome assembly, runs _mlst_ and reports
an individual file in each sample subdirectory.

# Outputs

Table with columns

* sample
* mlst_scheme
* 7-gene_ST
* locus1
* locus2
* locus3
* locus4
* locus5
* locus6
* locus7

# SYNOPSIS

`sn_helloWorld.py`

Python example plugin for helping developers.
Because this is an example, please see the Algorithm section below for more details.

# Software requirements

None

# Algorithm

This script is an example plugin for the basics.
The steps are:

1. Argument parsing with `argparse` to conform to the required flags and positional parameter of every plugin
2. In `main()`, check on whether there are some options where we print something and exit: `--help`, `--version`, or `--check-dependencies`.  This is done in `setup(args)`.
3. If we have not exited already, move onto `main(args)`
4. Make a temporary directory to avoid collision with other instances of this script. However, the user could have also supplied `--tempdir` and so the script has to accept that parameter if supplied.
5. Make the required subfolders just in case they are not there: `dir/SneakerNet/forEmail`
6. The example analysis
    1. `samples = readWriteSamples(args.dir)`: rewrite `SampleSheet.csv` to `dir + "/SneakerNet/forEmail/helloworld.py.tsv"`
    2. `readWriteFlags(args.dir, args)`: record all parameters to the output table at `dir + "/SneakerNet/forEmail/helloworld.py.tsv"`
7. Wrap up: write all properties to the central properties table. Because other plugins write to this, the plugin _appends_ and does not truncate.
    1. `writeProperties(args.dir, samples)`
    2. The central table for all plugins is at `dir + "/SneakerNet/properties.txt"`
    3. Two entries are placed into `properties.txt`: the version and the table path.

## Epilogue

1. `properties.txt` will be read by `sn_report.pl` which will email a formatted table in the plugin `emailWhoever.pl`.
2. `helloworld.py.tsv` will be included in the email because it is in the folder `SneakerNet/forEmail`. All files in that folder will be attached to the email when `emailWhoever.pl` runs.

# Outputs

## Table

This is a combined table of keys/values and samples/sampleCounter.
An example table is shown below, when the script was run with `--debug`.

|sample             | sampleCounter |
|-------------------|---------------|
|2010EL-1786        |1|
|FA1090             |2|
|LT2                |3|
|Philadelphia_CDC   |4|
|contaminated       |5|
|version            |False|
|check_dependencies |False|
|citation           |False|
|debug              |True|
|force              |False|
|tempdir            ||
|numcpus            |1|
|dir                |../t/M00123-18-001-test|

# SYNOPSIS

`sn_iontorrent_assembleAll.pl`

Assembles ion torrent raw read data and
predicts genes

# Software requirements

* SPAdes
* CG-Pipeline
* Prodigal

# Algorithm

Assembles raw read data with SPAdes with optimal options,
predicts genes with Prodigal,
then reports assembly and gene prediction metrics.

# Outputs

Table with columns

* sample
* genomeLength
* CDS - number of coding sequences
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC

Please view the whole readme in the main directory in [README.md](../../README.md)
# SYNOPSIS

`sn_detectContamination-mlst.pl`

Detects contamination using 7-gene MLST.

# Software requirements

* ColorID
* _mlst_

# Algorithm

Uses Torsten _mlst_ for the database. Maps raw reads
against the MLST database using ColorID. Because the
seven genes are housekeeping genes, then a non-seven
allele answer for any isolate indicates possible 
contamination.

# Outputs

Table with columns `mlst-contamination-detection.tsv`

* sample
* Scheme - scheme that this script estimates to be the correct scheme
* NumLociFound
* questionableLoci

# SYNOPSIS

`addReadMetrics.pl`

Adds raw read metrics for each fastq.gz file

# Software requirements

* CG-Pipeline

# Algorithm

Quantifies statistics like average read length.  Only measures 1% of the reads
and extrapolates for the whole file, and so some inconsistencies might be found
such as a single sample of a set that has a minimum read length of 38 when
all the others have 35 -- in this case, a minimum read length of 35 might be 
very rare but present.

# Outputs

Table with columns describing metrics of each reads file.
Reads are not interleaved and are quantified independently of either R1 or R2.

* File
* avgReadLength
* totalBases - i.e., total number of nucleotides
* minReadLength
* maxReadLength
* avgQuality
* numReads
* PE? - usually will be `no` because each split read is analyzed separately
* coverage - genome coverage calculated by `totalBases` divided by expected genome size found in `taxonProperties.conf`
* readScore - TODO
* medianFragmentLength - TODO
# SYNOPSIS

`sn_mlst-wg.pl`

Runs a whole-genome MLST-style analysis. Scheme can be
a whole genome MLST scheme, core genome MLST scheme,
or anything smaller.

# Software requirements

* chewBBACA
* BLAST+
* Database

## The database

The database must be listed under taxonProperties.conf. At least
one example is provided at this time for Salmonella.

The raw database must be downloaded. One method is with https://github.com/Public-Health-Bioinformatics/pubmlst_client

For example,

    pubmlst_download --scheme_name salmonella --scheme_id 4 --outdir ./salmonella_enterobase.raw

The database must be formatted according to chewBBACA instructions
at https://github.com/B-UMMI/chewBBACA/wiki/1.-Schema-Creation#14-using-an-external-schema

    chewBBACA.py PrepExternalSchema -i salmonella_enterobase.raw -o wgMLST/salmonella.enterobase.chewBBACA --cpu 4
    rm -rf salmonella_enterobase.raw

Next, the database must be located under your SneakerNet installation
directory, under `SneakerNet/db/wgMLST/something` where _something_
is the name of the scheme you provided in taxonProperties.conf.

# Algorithm

Runs chewBBACA which is a very very smart BLAST. For more information,
please see https://github.com/B-UMMI/chewBBACA

# Outputs

Summary table

For more information, please see https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling#allele-call-statistics-output-results_statisticstxt

* Genome (input fasta file)
* db (the database listed under taxonProperties.conf)
* EXC alleles which have exact matches (100% DNA identity) with previously identified alleles
* INF inferred new alleles using Prodigal CDS predictions
* LNF loci not found.
* PLOT possible loci on the tip of the query genome contigs.
* NIPH non-informative paralogous hit
* NIPHEM similar to NIPH classification (NIPH with exact match), but specifically referring to exact matches
* ALM alleles 20% larger than length mode of the distribution of the matched loci
* ASM similar to ALM but for alleles 20% smaller than length mode distribution of the matched loci 
# SYNOPSIS

`sn_report.pl`

Creates an HTML report for a SneakerNet run. Usually,
this is the second-to-last plugin to run, right before
the email plugin.

# Software requirements

# Algorithm

## key value pairs

Reads all properties in `SneakerNet/properties.txt`, which
has three columns: plugin, key, value. If a duplicate
plugin/key combination is found, then the latest value
is used.

## tables

If the value has a suffix `.csv` or `.tsv`, then the corresponding table is 
converted to HTML and included in the report. All other
key/value combinations are included in the HTML report,
under their corresponding plugins.

## images

_For this plugin version >= 2.7_

If the value has a suffix `.png` or `.gif`, then the corresponding
file path will be converted to base64 and embedded in the HTML report.

# Outputs

HTML report

## emoji output

This script itself outputs a table and is therefore converted to a table in the HTML report.
The table has columns:

* sample - sample name
* emoji - reflective of the score.  Happiest emojis reflect 100%.
  * As of v0.15, emoticons range from &#128515; (best), &#129320;, &#128556;, and &#128561; (worst).
* score - a percentage, starting from 100.  Each item under the failure_code column subtracts an equal percentage from 100%.  These possible failures are shown as columns in the [passfail plugin](sn_passfail.pl.md).  If there are three possible items, then each penalty is 33%.  By default in SneakerNet version 0.10, there are three possible items: coverage, quality, and kraken.
* qual - quality scores of R1 and R2, separated by space.
* cov - genome coverages of R1 and R2, separated by space.
* taxon - the calculated taxon. The calculated taxon is, in priority order: the taxon listed on the sample spreadsheet, guessed from Kraken, pattern matching on the filename, and lastly "UNKNOWN". For example, if the taxon is specified on the sample spreadsheet, it will be used and not overwritten.

# SYNOPSIS

`sn_parseSampleSheet.pl`

Parses a run to create a sample sheet.
Sets up sample sheet appropriately for downstream
plugins.

# Software requirements

# Algorithm

transforms the SampleSheet.csv file into samples.tsv

* Skips samples with the keyword 'empty'

# Outputs

samples.tsv

# SYNOPSIS

`sn_helloWorld.sh`

An example plugin, for helping developers

# Software requirements

# Algorithm

Creates a table with all options supplied to
this script reported into the table.

# Outputs

Table with columns

* key
* value

# SYNOPSIS

`sn_assembleAll_reference.pl.md`

Runs a reference assembly against a reference genome assembly
and produces assembly metrics.

# Software requirements

* Bowtie2
* samtools
* bcftools
* trimmomatic
* seqtk

# Algorithm

## Database building

* Uses `reference_fasta_id` to download a fasta and genbank file from NCBI
  * Found in taxonProperties.conf
  * Can use comma separated list, e.g., `reference_fasta_id=MH185784,MH185777,MH185772,MN367319,MN367321,MN367323,MN367326,MH430075`
* Downloads to SneakerNet installation, into a subfolder `db/fasta`
* Indexes the fasta file with `samtools` and `bowtie2-build`
* Gene annotations are not indexed and can be read on the fly

## Assembly

* Downloads the reference genome into your installation directory
* Trims adapters with trimmomatic
  * primers are found in taxonProperties.conf
  * e.g., primers_bed_url=https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.primer.bed
  * If no primers, then will not trim and will produce a warning
* Maps reads against genome with bowtie2
* Calls SNPs with samtools
* Filters SNPs
* Consensus fasta from SNPs and the reference genome

## Annotation

* copies coordinates of reference genome genes to assembly
* Counts intact coding sequences designated by the CDS tag in the genbank file
* Counts intact coding sequences in the new assembly

# Outputs

Assemblies are found in `SneakerNet/assemblies/[samplename]/*.fasta`

Some intermediate files are in `SneakerNet/assemblies/[samplename]/consensus
* depth.tsv - output of `samtools depth`
* out.vcf.gz.masked.vcf.gz and related `csi` file
* sorted.bam and related `bai` file

A graph of all depths of coverage in `SneakerNet/forEmail/depth.png`

Table with columns

* File - the sample name with some file extensions
* genomeLength
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC
* effectiveCoverage - the number of mapped base pairs against the assembly. Will be lower than the coverage noted in readMetrics.
* percentNs - the percentage of Ns in the assembly, from 0 to 1
* altCdsCount - number of CDS in new assembly
* refCdsCount - number of CDS in reference assembly
* expectedCdsPercentage - altCdsCount / refCdsCount
* longestContiguous - the longest contiguous part of the genome between any two Ns



# SYNOPSIS

`assembleAll.pl`

Assembles raw reads into assembly and then
predicts genes.

# Software requirements

In version 2, Shovill was added and many requirements were added

* Shovill - has many dependencies itself
  * BWA
  * Flash
  * Java
  * Lighter
  * Mash (before shovill 1.1.0)
  * KMC (in shovill 1.1.0)
  * Megahit (not used but checked when shovill runs)
  * pigz
  * Pilon
  * Samclip
  * Samtools
  * Seqtk
  * Skesa
  * Trimmomatic
  * Velvet  (not used but checked when shovill runs)
* CG-Pipeline
* Prodigal

# Algorithm

Assembles genomes with Skesa, predicts genes with prodigal,
then produces assembly and gene prediction metrics.
Contigs with length < 500 are not considered in the assembly metrics.

# Outputs

Assemblies are found in `SneakerNet/assemblies/[samplename]/*.fasta`

Table with columns

* sample
* genomeLength
* CDS - number of coding sequences
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC
* effectiveCoverage - the number of mapped base pairs against the assembly. Will be lower than the coverage noted in readMetrics.

Starting with v2.5, four more output files are found in
`SneakerNet/assemblies/[samplename]/` but not included in the final report.
These files are used to make the aforementioned combined metrics table.

* predictionMetrics.tsv - gene prediction metrics such as CDS count
* assemblyMetrics.tsv - N50, genomeLength, etc
* depth.tsv.gz - `samtools depth` output: a three column file with contig, pos, depth of coverage
* effectiveCoverage.tsv - contains effective coverage

# SYNOPSIS

`sn_detectContamination.pl`

Runs a kmer histogram to detect potential contamination.

# Software requirements

* Perl module `Bio::Kmer`

# Algorithm

Runs a kmer histogram. The null hypothesis is that there will be
a peak around 1x coverage due to errors in sequencing, and a
second peak around the expected genome coverage. If something
else were sequenced at the same time at a different coverage,
then there will be another peak at that coverage level.
This other peak will be indicative of potential contamination.

# Outputs

Table with columns

* File
* numPeaks
* finalDelta - the amplitude difference between one coverage and the next to determine what a valley is
* hist       - ascii representation
* firstPeak  - coverage for the first peak
* firstValley - coverage for the first valley
* secondPeak
* secondValley...

More columns are possible with additional valleys

# SYNOPSIS

`sn_staramr.pl`

Runs antimicrobial resistance analysis

# Software requirements

* StarAMR >= 0.7.0
* StarAMR database

## StarAMR database

Build the database into the correct folder under SneakerNet
like so

    staramr db build --dir SneakerNet/db/staramr

Under taxonProperties.conf, add the following line under each taxon for pointfinder (SNP)
resistance detection. Most taxa are not supported right now. At the time of this documentation,
only these are supported: salmonella, campylobacter, enterococcus faecalis and enterococcus faecium.

    pointfinder=salmonella

# Algorithm

Runs StarAMR on genome assemblies.  Pointfinder is not currently enabled.

# Outputs

Table with columns

* Sample - sample name
* Assembly - which assembly was input
* Genotype 
* Predicted Phenotype

# SYNOPSIS

`sn_helloWorld.pl`

An example plugin, for helping developers

# Software requirements

# Algorithm

Creates a table with all options supplied to
this script reported into the table.

# Outputs

Table with columns

* key
* value

# SYNOPSIS

`guessTaxon.pl`

Attempts to guess a taxon for each sample

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

Runs sample raw reads through Kraken1. Quantifies percentage of
reads that match each taxon.
Compares the assigned taxon vs the assumed taxon.

# Outputs

Table with columns

* sample
* assumed taxon
* best-fitting taxon
* percent of reads that match the best-fitting taxon

# SYNOPSIS

`baseBalance.pl`

Compares the ratio of nucleotides vs expected value

# Software requirements

# Algorithm

For each read set, divides all As by Ts and all
Gs by Cs. The expected ratio for each set of raw reads

# Outputs

`TODO`

# SYNOPSIS

`transferFilesToRemoteComputers.pl`

Transfers files to a remote computer as specified in
config/settings.conf

# Software requirements

* rsync

# Algorithm

Will transfer each set of raw reads to the remote
computer.  Will not transfer if transfer was not
requested in samples.tsv. Also will not transfer
if did not pass in `sn_passfail.pl`.

# Outputs

none

# SYNOPSIS

`sn_crypto_assembleAll.pl`

Assembles raw reads into assembly and then
predicts genes.

# Software requirements

* Skesa
* CG-Pipeline
* Prodigal
* Shovill

# Algorithm

Assembles genomes with Shovill/Skesa, predicts genes with prodigal,
then produces assembly and gene prediction metrics.

# Outputs

Table with columns

* sample
* genomeLength
* CDS - number of coding sequences
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC


# SYNOPSIS

`sn_kraken.pl`

Runs Kraken1 on a sample.

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

Runs Kraken1 on a sample.
It looks for raw sequence reads (fastq) as input.
If not found, runs on genome assemblies under
runDir/SneakerNet/assemblies/sampleName/something.fasta.

# Outputs

No table output.
Creates files under runDir/SneakerNet/kraken/sampleName:

* kraken.report - a summary of counts per taxon with the following fields. More information here: https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports
  * Percentage of reads covered by the clade rooted at this taxon
  * Number of reads covered by the clade rooted at this taxon
  * Number of reads assigned directly to this taxon
  * A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
  * NCBI taxonomy ID
  * indented scientific name
* kraken.filtered.report - the same as above, but all hits < 1% are removed
* kraken.taxonomy - a tab-delimited file
  * first column is the number of reads
  * second column, and all subsequent columns are the taxonomy starting with root and going down as far as possible, e.g., scientific name.
* report.html - graphical output in Krona format

# SYNOPSIS

`sn_crypto_gp60`

Runs gp60 on each genome

# Software requirements

* blast+
* GP60_Counter (currently a private repo)

# Algorithm

For each genome assembly, provides the GP60
profile.

# Outputs

Table with columns

* sample
* GP60 type
# SYNOPSIS

`sn_sarscov2_assembleAll.pl`

Runs a reference assembly against Wuhan-1
and produces assembly metrics.

# Software requirements

* Bowtie2
* samtools
* bcftools
* trimmomatic
* seqtk

# Algorithm

## Database building

* Uses `reference_fasta_id` to download a fasta and genbank file from NCBI
* Downloads to SneakerNet installation, into a subfolder `db/fasta`
* Indexes the fasta file with `samtools` and `bowtie2-build`
* Gene annotations are not indexed and can be read on the fly

## Assembly

* Downloads the reference genome into your installation directory
* Trims adapters with trimmomatic
* Maps reads against genome with bowtie2
* Calls SNPs with samtools
* Filters SNPs
* Consensus fasta from SNPs and the reference genome

## Annotation

* copies coordinates of reference genome genes to assembly
* Counts intact coding sequences designated by the CDS tag in the genbank file
* Counts intact coding sequences in the new assembly

# Outputs

Assemblies are found in `SneakerNet/assemblies/[samplename]/*.fasta`

Some intermediate files are in `SneakerNet/assemblies/[samplename]/consensus
* depth.tsv - output of `samtools depth`
* out.vcf.gz.masked.vcf.gz and related `csi` file
* sorted.bam and related `bai` file

A graph of all depths of coverage in `SneakerNet/forEmail/depth.png`

Table with columns

* File - the sample name with some file extensions
* genomeLength
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC
* effectiveCoverage - the number of mapped base pairs against the assembly. Will be lower than the coverage noted in readMetrics.
* percentNs - the percentage of Ns in the assembly, from 0 to 1
* altCdsCount - number of CDS in new assembly
* refCdsCount - number of CDS in reference assembly
* expectedCdsPercentage - altCdsCount / refCdsCount
* longestContiguous - the longest contiguous part of the genome between any two Ns



