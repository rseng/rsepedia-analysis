# DeeReCT-APA: Prediction of Alternative Polyadenylation Site Usage through Deep Learning
## Overview
<div align="center">
  <img src="./resources/DeeReCT-APA.png" width="600" height="450">
</div>
This repository contains the implementation of DeeReCT-APA from 

Li,Zhongxiao, et al. "DeeReCT-APA: Prediction of Alternative Polyadenylation Site Usage Through Deep Learning" bioRxiv, 2020

If you use our work in your research, please cite our paper, which is now available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.26.009373v1)
```
@article {Li2020.03.26.009373,
	author = {Li, Zhongxiao and Li, Yisheng and Zhang, Bin and Li, Yu and Long, Yongkang and Zhou, Juexiao and Zou, Xudong and Zhang, Min and Hu, Yuhui and Chen, Wei and Gao, Xin},
	title = {DeeReCT-APA: Prediction of Alternative Polyadenylation Site Usage Through Deep Learning},
	elocation-id = {2020.03.26.009373},
	year = {2020},
	doi = {10.1101/2020.03.26.009373},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```

## Dataset Description
The datasets used in training is stored under `APA_ML/`. The parental datasets are in `APA_ML/Parental` and the F1 datasets are in `APA_ML/F1`.

### Parental Dataset
In folder ``Parental/``:
- ``bl.pAs.sequence.txt``: Contains genomic sequences of PAS sites in BL mouse. The rows are indexed by PAS ID. The PAS ID consists of the Ensembl ID of the gene and PAS number. As an example, the PAS ID ``ENSMUSG00000056763.13:3`` indicates the third PAS of gene ``ENSMUSG00000056763.13``. The PAS are numbered from **downstream** to **upstream**, therefore ``ENSMUSG00000056763.13:1`` is the last PAS of gene ``ENSMUSG00000056763.13``. The first column is ``pas_type``, there are five PAS types: ``altLastExon``, ``internalExon``, ``intronicPAS``, ``tandemUTR``, ``terminal``. The second column contains the coordinates of the PAS in BL and SP mouse. The +/- indicates which strand in the genome is the *template strand* of the transcript containing the PAS. *When only one strain's coordinate is present, it means that the PAS is only identified in that parental strain by 3'READS method. If the coordinates of both two parental strains are present, that PAS is identified in both parental strains.* The third column contains the sequences of the PAS, same as the sequences from their transcribed mRNA sequences, except that T should be U. There are 31888 PAS in total.
- ``sp.pAs.sequence.txt``: Contains genomic sequences of PAS sites in SP mouse. There are 31888 rows in total.
- ``parental.pAs.usage.txt``:  The percent usages computed from 3'-mRNA-seq. There are 30940 rows in total because some PAS, due to insufficient amount of supporting reads, are discarded.

### F1 Dataset
In folder ``F1/``:
-  ``control.f1.usage.txt``: Use data of 'controlled genes', whose usage do not change much between two strains. In many genes, some of their PAS are discarded due to insufficient amount of supporting reads. The supporting read counts are much lower than that in parental strains because many reads don't have allelic specific SNPs and therefore one cannot distinguish which strain they come from. 
- ``differential.f1.usage.txt``: Data of 'differentially expressed genes', whose usage change between two strains.

## Prerequisites
The code is tested with the following dependencies:
- Python 3.6
- Pytorch 1.2.0
- Numpy 1.17.2
- Scipy 1.3.1
- Pandas 0.25.1
- PyTables 3.5.2

The code is not guaranteed to work if different versions are used.

## Dataset Preprocessing 
**Note: All steps are mandatory and should be executed in the order presented**

To preprocess the parental dataset, one should run the following script in the repository folder
```
python parental_data_preprocess.py
```
To preprocess the F1 dataset, one should run
```
python F1_data_preprocess.py
```
To generate the 5 fold cross-validation dataset, one should run
```
python cross_validation_data_preprocess.py
```
The three steps will take some time and progress will be shown on the command line.
## Training
By default, we are using the first four folds to train and the remaining one to test. One can modify the parameter at the beginning of `main.py` to change the hyperparameters, output paths, train/test splits, etc. To start training, one simply needs to run
```
python main.py
```
The evaluation results will be displayed on the command line.
