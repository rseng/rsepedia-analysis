# Opioid Receptors
*This paper describes an exciting big data analysis compiled in a freely available database, which can be applied to characterize the coupling of different G-Protein coupled receptors (GPCRs)  families with their intracellular partners. Opioid receptor (OR) family was used as case study in order to gain further insights into the physiological properties of these important drug targets, known to be associated with the opioid crisis, a huge socio-economic issue directly related to drug abuse. An extensive characterization of all members of the ORs family (μ (MOR), δ (DOR), κ (KOR), nociceptin (NOP)) and their corresponding binding partners (ARRs: Arr2, Arr3; G-protein: Gi1, Gi2, Gi3, Go, Gob, Gz, Gq, G11, G14, G15, G12, Gssh, Gslo) was carried out. A multi-step approach including models’ construction (multiple sequence alignment, homology modeling), complex assembling (protein complex refinement with HADDOCK and complex equilibration), and protein-protein interface (PPI) characterization (including both structural and dynamics analysis) were performed. Our database can be easily applied to several GPCR sub-families, to determine the key structural and dynamical determinants involved in GPCR coupling selectivity.*

If you use the code or tools available in this repository, please cite:
1. Carlos A.V. Barreto&ast;, Salete J. Baptista&ast;, A. J. Preto, Daniel Silvério, Rita Melo, Irina S. Moreira. Decoding partner specificity of opioid receptor family. Frontiers in Molecular Biosciences, 2021. &ast; joint first authors
2. A.J. Preto, Carlos A.V. Barreto, Salete J. Baptista, José Guilherme de Almeida, Agostinho Lemos, André Melo, M. Natália D. S. Cordeiro, Zeynep Kurkcuoglu, Rita Melo and Irina S. Moreira. Understanding the Binding Specificity of G-Protein Coupled Receptors toward G-Proteins and Arrestins: Application to the Dopamine Receptor Family. Journal of Chemical Information and Modeling, 2020.

Furthermore, if you which to access the website for the Opioid receptor analysis, access http://www.moreiralab.com/resources/oxr/. If you which to access the website with the Dopamine receptor analysis, go to http://www.moreiralab.com/resources/dxr/.

# Deployment instructions
This is the repository for Opiod Receptor Partners analysis code. By calling following this instructions, the user should be able to run the full pipeline for the opioid receptors profile. In order to make sure everything runs smoothly, the user should follow the upcoming steps:
1. Clone the repository
2. In the **gpcr_variables.py** script, change the **DEFAULT_FOLDER** variable into the users' location of the repository
3. Create a conda environment: `conda create --name opioid_receptors R=3.6 python=3.9`
4. Activate the conda environment: `conda activate opioid_receptors`
5. Add some Python packages: `pip install pandas selenium numpy toolz bs4 requests biopython`
6. Prepare for webscrapping:
	- Confirm your installed Chrome version
	- Download the corresponding chromedriver executable (check: https://chromedriver.chromium.org/downloads)
	- Move the executable to the folder
7. Add some R packages: `RScript  -e "install.packages(c('circlize','stringr','ggplot2','bio3d','tidyverse','ggrepel','ggsci','cowplot','svglite'), repos='https://cran.rstudio.com/')"`
8. Run the pipeline: `python call.py`

# Folder structure
The **pdb** files containing the dimeric structures should be in the root folder. Subdirectories should be:
- results
- processed_results
- images
- summary
- templates
- dymeric_complexes
- structural_complexes

After deploying the pipeline, all this folders should, however, be automatically created.
