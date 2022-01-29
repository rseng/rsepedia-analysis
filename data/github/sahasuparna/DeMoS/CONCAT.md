# DeMos
#CESCclinicalMatrix.csv---> The file contains the clinical information of cervical cancer patients.

#genomic_miRNA---> contains the miRNA expression profile of cervical cancer.

#SCC_SAMPLEID.csv-->contails the sample id for SCC group patients.

#cervical.ADENO.vs.SCC.DRmiRNA.csv--->contains the downregulated miRNA for ADENO vs SCC group.

#cervical.ADENO.vs.SCC.URmiRNA.csv--->contains the upregulated miRNA for ADENO vs SCC group.

DEGstat_NewmiRNA.R--> Function for finding the significantly expressed miRNA

DEGstat_New_sm.R--->Function for finding the significantly expressed mRNA

classification.R--->Function for classification of the class labels(ADENO,SCC) by our gene signature as well as comparison with others.

new_surv.R--->Function for prediction of survival of CESC patients.
# Pre-requisites
R version 3.5.3
