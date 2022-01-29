# MHCVision
The tool for global and local false discovery rate (FDR) estimation for MHC-peptide binding prediction
### **Introduction**
MHCVision is the model for global and local FDR estimation for the predicted scores of MHC-peptide binding affinity. The model utlilises the approach of the Expectation Maximisation (EM) algorithm with the method of moments to estimate the parameters of data distribution for determining the relative true and false data for global and local FDR calculation. The current version is feasible for the predicted data using NetMHCpan (version >= 4.0) or MHCflurry. 

If you find MHCVision is useful in your research, please cite:

Pearngam, Phorutai, Sira Sriswasdi, Trairak Pisitkun, and Andrew R. Jones. "MHCVision: estimation of global and local false discovery rate for MHC class I peptide binding prediction." Bioinformatics (2021). https://doi.org/10.1093/bioinformatics/btab479

### **How to install?**
The model requires Python 3 ( >= 3.7) and the following python packages:
```
pandas (>= 1.1.2)
numpy (>= 1.19.1)
scipy (>=1.5.2)
```
For python installing packages, please see [here](https://packaging.python.org/tutorials/installing-packages/#use-pip-for-installing)

If your system has both Python 2 and Python 3, please ensure that Python 3 is being used when following these instructions.

To install the model:

1. Clone this repository
```
git clone https://github.com/PGB-LIV/MHCVision
```
For other methods for cloning a GitHub repository, please see  [here](https://help.github.com/articles/cloning-a-repository/)

2. Install the latest version of 'pip' and 'setuptools' packages for Python 3 if your system does not already have them
```
python -m ensurepip --default-pip
pip install setuptools
```
For more information, please see [here](https://packaging.python.org/tutorials/installing-packages/#install-pip-setuptools-and-wheel)

3.  Run Setup.py inside MHCVision directory to install the model
```
cd MHCVision
python Setup.py install
```

### **Usage**
```
usage: mhcvision.py [options] input_file.csv -o/--output output_file.csv
options:
-a, --allele   REQUIRED: type the allele name i.e. HLA-A0101, which are supported in the "supplied_alleles.txt"
-t, --tool     REQUIRED: Specify the MHC-peptide prediction tool you used, type NetMHCpan or MHCflurry
-i, --input    REQUIRED: specify the input filename, the input file must be in ".CSV" format (comma-separated values), the column headers must contain 'Peptide', 'IC50', '%Rank'
-o, --output   Optional: specify the output filename 
-h, --help     Print the usage information
```

### **Sample scripts**
You can use sample.csv as the input file
```
python mhcvision.py -a HLA-A1101 -t NetMHCpan -i sample.csv
```
