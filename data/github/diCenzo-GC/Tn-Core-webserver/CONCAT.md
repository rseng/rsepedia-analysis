# Tn-Core Webserver

Transposon-sequencing (TnSeq) and genome-scale constraint-based metabolic modelling are highly complementary experimental and *in silico* approaches. Workflows combining these two approaches have the potential to provide insights not obtainable by either method alone. The Tn-Core Toolbox was designed to provide tools for integrating these methods, facilitating the contextualization of TnSeq data and the formation of context-specific core metabolic reconstructions. The output models provide a framework to assist in the functional interpretation of TnSeq data and in identifying gaps in the data, in addition to being useful in downstream modelling pipelines.

The Tn-Core webserver is a web-based implementation of one of the main functionalities of Tn-Core: extraction of a context-specific core metabolic network on the basis of TnSeq data. The webserver returns the core metabolic network in two formats: i) as an easy-to-understand Excel file, and ii) as a COBRA-formatted model as a MATLAB file. The webserver was developed to (hopefully) be accessible even to those with little to no experience in metabolic modelling. Below, you will find a description of the input requirements and the outputs. Additionally, this repository contains a folder with five files for use as sample data to test the Tn-Core webserver.


## Contacts

If you have questions about the use of the Tn-Core webserver, please send inquiries to George diCenzo (george.dicenzo [at] queensu.ca) or Marco Fondi (marco.fondi [at] unifi.it).


## Principle

The main idea behind the Tn-Core webserver is to use metabolic reconstructions to functionally interpret TnSeq data. A genome-scale metabolic reconsruction is meant to capture the entire metabolism of the cell, with every reaction linked to the gene(s) whose gene product(s) catalyze the reaction. One of these reactions represents the formation of biomass by combining the relevant molecules (protein, DNA, RNA, lipids, etc.) in the appropriate ratio. Mathematical methods, such as FBA, can be used to run simulations with these reconstructions, with the aim of the simulations being to produce as much biomass as possible with a given set of input nutrients.

Essentially, Tn-Core takes a genome-scale metabolic reconstruction and attempts to produce a reduced 'core' reconstruction that is still able to produce biomass, using the tranTnSeq data (and optionally also RNA-sequencing data) as a guide to select which genes (and thus reactions) to remove from the reconstruction. The final output is therefore a context-specific core reconstruction, containing just the metabolism that is essential for growth in a given environment. These reconstructions, or essentially a list of metabolic reactions, can be used to gain a functional understanding of the TnSeq data; e.g., to uncover the metabolic functions that are essential in a given environemt. Importantly, they can identify metabolic reactions that are essential in a given environment that for some reason (e.g., functional redundancy) did not appear essential in the TnSeq data. In this way, Tn-Core can help provide a more comprehensive understanding of the metabolic processes underlying growth in an environment of interest.


## Required Inputs

### Metabolic Reconstruction
A metabolic reconstruction for your organism of interest must be obtained. This can be your own metabolic reconstruction, or can be downloaded from the web (for example, from the BiGG or VMH databases). The reconstruction should be in SBML format and saved as a .xml file. Additionally, it is important that the reconstruction have a biomass reaction, and that the gene names in the reconstruction match the gene names in the TnSeq data.

### Exchange Reactions
Exchange reactions provide a way of setting the growth medium for the simulations. Here, you should provide a tab-delineated text file providing the reaction identifiers of the exchange reactions in the 1st column and the maximal influx rate (a negative value) in the 2nd column. If unsure of what values to use in the 2nd column, using -10 for all reactions should provide a reasonable starting point). Tn-Core will still run and provide output if you simply provide an empty text file here. However, the closer the *in silico* medium resembles the experimental conditions, the better the output is likely to be. We recommend at minimum including the exchange reactions for the carbon source, the nitrogen source, the phosphorus source, and the sulphur source.

### Objective Reaction
This should indicate the biomass reaction. Here, you should provide a tab-delineated text file consisting of one line: the reaction identifier of the biomass reaction, followed by the objective coefficient (should always be 1).

### TnSeq Data
This should be a tab-delineated text file containing the TnSeq data. The 1st column should be a summary statistic of the TnSeq data, and the 2nd column be the gene names (ensure the gene names are the same as the gene names in the metabolic reconstruction). The TnSeq data should not be log-transformed, and should be in a format where the lower the number, the more important the gene. For example, it could consist of the number of reads mapping to the gene divided by the gene length. This file can contain data for genes that are not in the reconstruction, and it does not have to contain data for every gene in the model.


## Optional Inputs

### RNAseq Data
This should be a tab-delineated text file containing the RNAseq data. The 1st column should be the RNAseq data as either RPKM or TPM, and the 2nd column should be the gene names (ensure the gene names are the same as the gene names in the metabolic reconstruction). This file can contain data for genes that are not in the reconstruction, and it does not have to contain data for every gene in the model.

### Minimum Growth Fraction
This should be a value greater than 0 and lesser than or equal to 1 (i.e., 0 < n ≤ 1). This number sets the minimal rate that the output model can produce biomass relative to the input genome-scale metabolic model. For example, if the input model produces 0.3 g of biomass per hour per gram of biomass, and this value is set to 0.5, then the output model must produce at least 0.15 g of biomass per hour per gram of biomass. If no value is entered, the default is 0.5.

### Expression Threshold
This should be the minimum RPKM or TPM value for a gene to be considered 'highly-expressed'. This could be, for example, the average RPKM value of all genes. If no value is given, the default is to calculate the threshold as 0.02% the sum of all the expression values.


## Outputs

### Excel File
The Excel file contains three sheets summarizing the core metabolic reconstruction.

The first sheet include information on the reactions. The following columns are included: RxnID - the reaction identifier; RxnName - name of the reaction (i.e. brief description of the activity); Reaction_A - the reaction formula using the metabolite identifiers; Reaction_B - the reaction formula using actual metabolite names; Reversible - indicates if the reaction is reversible (true) or not reversible (false); KEGG_RID - the KEGG reaction ID of the reaction (only included if provided in the input model); EnzymeClass - the enzyme classification number (only included if provided in the input model); Genes - the genes associated with the reaction; Proteins - the proteins encoded by the genes associated with the reaction (only included if provided in the input model).

The second sheet includes information on the genes present in the core metabolic reconstruction. The following columns are included: Gene - the gene name; Protein - the name of the protein encoded by the gene (only included if provided in the input model); Tn_Seq_Data - the TnSeq data for the gene from the input TnSeq file; RNAseq_Data - the RNAseq data for the gene from the input RNAseq file (only included if RNAseq data are provided).

The third sheed includes information on the metabolites present in the core metabolic reconstruction. The following columns are included: MetIDs - the metabolite identifiers; MetNames - the names of the metabolties.

### MATLAB File
The MATLAB file contains the core metabolic reconstruction as a COBRA-formatted model.

## Citation

If you find the output of the Tn-Core webserver to be helpful, we ask you to consider including a statement such as the following in the methods section of your manuscript:

Transposon-sequeing data was integrated with the metabolic reconstruction using the Tn-Core webserver (1), which is dependent the Tn-Core Toolbox 2.1 (1), the COBRA Toolbox (2), the TIGER Toolbox 1.2.0-beta (3), FASTCORE 1.0 (4), GIMME (5), the SBMLToolbox 4.1.0 (6), libSBML 5.13.0 (7), MATLAB 2016b (mathworks.com), and the iLOG CPLEX Studio 12.7.1 solver (ibm.com).

(1) diCenzo GC, et al. (2019) Tn-Core: a toolbox for integrating Tn-seq gene essentiality data and constraint-based metabolic modelling. ACS Synth Biol. 8: 158-169.

(2) Schellenberger J, et al. (2011) Quantitative prediction of cellular metabolism with constraint-based models: the COBRA Toolbox v2.0. Nat Protoc. 6: 1290-1307.

(3) Jensen PA, et al. (2011) TIGER: Toolbox for integrating genome-scale metabolic models, expression data, and transcriptional regulatory networks. BMC Syst Biol. 5: 147.

(4) Vlassis N, et al. (2014). Fast reconstruction of compact context-specific metabolic network models. PLOS Comput Biol. 10: e1003424.

(5) Becker SA, Palsson BO. (2008) Context-specific metabolic networks are consistent with experiments. PLOS Comput Biol. 4: e1000082.

(6) Keating SM, et al. (2006) SBMLToolbox: an SBML toolbox for MATLAB users. Bioinformatics. 22: 1275-1277.

(7) Bornstein BJ, et al. (2008) LibSBML: an API library for SBML. Bioinformatics. 24: 880-881.
# Tn-Core Webserver

This directory contains the files used to run the Tn-Core web application. These files are based heavily on the files used to run the MeDuSa web application (https://github.com/combogenomics/medusa-webapp) prepared by Marco Galardini.

## Test the server

"Redis" and "atd deamon" need to be running prior to testing the server. It is also necessary to place a 3.x version of bootstrap in the "static" directory, as well as jquery.min.js in the "static/js" directory. MATLAB must also be installed, and the iLOG CPLEX solver must be on the MATLAB path.

Create a directory called "Software". Within this directory, place the following directories: i) "cobratoolbox" containing the COBRA Toolbox , ii) "FastCore" containing FASTCORE, iii) "libSBML" containing libSBML, iv) "SBML-Toolbox" containing the SBML Toolbox, v) "tiger" containing the TIGER Toolbox, and vi) "Tn-Core" containing the Tn-Core Toolbox. 

Set up the virtual environment with the following code:

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt
    python tncore.py

Open a browser on the same machine and go to 127.0.0.1:5000

## Production

First, install all the dependencies:

    sudo pip install -r requirements.txt

Then install in apache the tncore.conf file (changing the paths).

Create a production.py file which can then be used to override the settings.py debug options.

Restart apache and start redis.

If the tncore.py script is modified, please resart apache to initialize the changes.

Edit the mail_log.py file to setup the error logging through email.

You may also want to set up a cron job to wipe out the uploads directly every now and then
# Sample Data

Five files are provided for testing the Tn-Core webserver.

### metabolic_reconstruction.xml
A sample file to enter in the 'Metabolic reconstruction' field.

### exchange_reactions.txt
A sample file to enter in the 'Exchange reactions' field.

### objective_reaction.txt
A sample file to enter in the 'Objective reaction' field.

### tnseq_data.txt
A sample file to enter in the 'TnSeq Data' field.

### rnaseq_data.txt
A sample file to enter in the 'RNAseq Data' field.
