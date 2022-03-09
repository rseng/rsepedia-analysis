# NetCom 
This code is a local version of the tool described in the paper:
NetCom: a network-based tool for predicting metabolic activities of
microbial communities based on interpretation of metagenomics data
Authors: Ofir Tal, Rotem Bartuv, Maria Vetcos, Shlomit Medina, Jiandong
Jiang and Shiri Freilich

Online version
https://freilich-lab-tools.com/netcom/

# Instructions
Input file: EdgeR results file, tab separated consist the columns [enzyme	logFC	logCPM	PValue	FDR	association].
After a file is uploaded, a new section appears, enables the selection of the control and the experiment treatments. 
The selection of control and experiment treatments drives the calculation of basic statistics which are now shown.
Entities number in a pathway: select a range of minimum and maximum number of entities linked with a pathway in the enrichment analysis.  
The selection of the range above drives the enrichment analysis calculation and a new option is then enabled - 
Select pathways to dropout: select which pathways to drop out of the network.
Environmental resource node color: select the color of the compounds which are essential for the network to develop (seeds).
Unique node color: select the color of the compounds which are differentially abundant due to the treatment.
Limit node hubness: Select the maximum allowed number of connected edges to a node.
Set network layout iterations: NetworkX network layout iterations. Higher numbers would produce an arranged network, but consume more time to calculate (up to a few minutes).

# The tool is implemented in Python 3.8 (Linux OS, Anaconda) using the following packages:
pandas 1.0.5
plotly 4.14.1
dash 1.18.1
dash_core_components 1.14.1
dash_html_components 1.1.1
networkx 2.4
dash_bootstrap_components 0.11.1
matplotlib 3.2.2
dash_extensions 0.0.41
numpy 1.18.5
scipy 1.5.0
statsmodels 0.11.1

# Parameters
fisher exact test - scipy.stats.fisher_exact, alternative hypothesis 'greater'

Test results and p-value correction for multiple tests - statsmodels.stats.multitest.multipletests default parameters, alpha=0.05

# Output files:
-T1/T2_ECs: enzymatic reactions that are differentially abundant in one of the treatments (according to user's input file)
-T1/T2_resources: a list of compounds (KEGG accessions) that were predicted as treatment-specific environmental resources
-T1/T2_compounds: all compounds included in the treatment-specific network expanded for the treatment (KEGG accessions)
-T1/T2_Enymes_pathways: pathways that are enriched with environmental resources that are unique to the treatment
-T1/T2_resource_pathways: pathways that are enriched with treatment-specific environmental resources
-T1/T2_pathways: pathways that are enriched with compounds that are unique to the treatment
-3D_network_T1/T2: html files with 3D networks of the full treatment network or of specific pathways that are enriched with compounds that are unique to the treatment
-T1/T1_Network: high resolution figure file of the 2D network
-3D_network_<treatment_name>.html

# run locally - recommended approach

mkdir netcom

cd netcom

virtualenv netcom

source netcom/bin/activate 

pip install -r requirements.txt


execute by:

python app.py

open browser at http://127.0.0.1:8050/netcom/

upload file for analysis



