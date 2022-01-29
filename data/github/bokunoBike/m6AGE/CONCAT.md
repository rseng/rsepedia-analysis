# m6AFinder
## Description
M6AFinder is a predictor for N6-methyladenosine (m6A) sites identification utilizing sequence characteristics and graph embedding-based geometrical information.

## Requirements
The code is based on python 3.7 It depends on several python packages, such as numpy, scikit-learn, networkx, pandas, catboost, karateclub.
* conda install numpy, pandas, scikit-learn, biopython, networkx
* pip install catboost, karateclub

## Usage

Command line usage:
```
$python main.py [-pos_fa pos_fa] [-neg_fa neg_fa] [-test_fa test_fa]
         [-dataset dataset] [-out out]
```
for example:
```
$python main.py -pos_fa ./datasets/A101/A101_Train_P.fasta -neg_fa ./datasets/A101/A101_Train_N.fasta -test_fa ./datasets/A101/A101_Test.fasta -dataset A101 -out ./results/A101_results.csv
```