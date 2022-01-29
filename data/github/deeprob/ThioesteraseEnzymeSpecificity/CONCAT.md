

[![DOI](https://zenodo.org/badge/286865245.svg)](https://zenodo.org/badge/latestdoi/286865245)

# EnZymClass Tutorial


## Environment creation

### Setting up a conda virtual environment
1.  Install conda preferably through anaconda. For installation instructions, please visit [link](https://docs.anaconda.com/anaconda/install/)

2. Set the conda channel priority list by editing .condarc file to the following:

    channel_priority: flexible \
    channels:\
     \- conda-forge\
     \- bioconda\
     \- defaults

    This step is optional but is recommended as it can circumvent some installation errors. The .condarc file can be found in the home directory. In macOS/linux it can be edited from the terminal using the following command: 
    
    \$nano ~/.condarc
    
3. Create a conda environment with the following command:

   $conda create -n te_env python=3.9 scikit-learn pandas jupyter blast bioconductor-kebabs=1.24.0

4. Install ifeatpro [link](https://pypi.org/project/ifeatpro/), ngrampro [link](https://pypi.org/project/ngrampro/) and pssmpro [link](https://pypi.org/project/pssmpro/) using pip after activating the conda environment. Run the following commands in MacOS/Linux. If you have Windows, use WSL.
   
   \$ conda activate te_env
   
   \$ pip install ifeatpro
   
   \$ pip install ngrampro
   
   \$ pip install pssmpro


## Reproducing TE substrate specificity results 
The jupyter notebook, *TE_SubstrateSpecificityAnalysis.ipynb* present in the *notebooks/* directory provides step by step instructions on how we obtained the current results. Open a jupyter session and rerun the notebook. Some steps take several hours to run; please use multiple cores (I used 24) to attain results within a reasonable amount of time.  


## Applying EnZymClass on other protein sequence classification applications
1. Create train and test dataset of the same format as the csv files given in *data/raw/* directory. The format should be as follows:

    protein_name, protein_sequence, protein_label (for training dataset)
    
    protein_name, protein_sequence (for test dataset)
    
    
2. Duplicate and rename the *TE_SubstrateSpecificityAnalysis.ipynb* present in the *notebooks/* directory according to your application area. Rename *train_raw* and *test_raw* variables in the notebook

3. Run the notebook step by step. Please note that the pssm based features require creation of pssm profiles of protein sequences which in turn require the psiblast program path and a blast database. Blast database creation is described in the pssmpro tutorial [link](https://pypi.org/project/pssmpro/).  
