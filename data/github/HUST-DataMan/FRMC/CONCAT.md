# FRMC 

Python code for Wu and Wang et al. FRMC: A fast and robust method for the imputation of scRNA-seq data

FRMC can not only precisely distinguish "true zeros" from dropout events and correctly impute missing values 
attributed to technical noises, but also effectively enhance intracellular and intergenic connections and achieve 
accurate clustering of cells in biological applications.

* **FRMC.py**  
Last update: v1.0.0

This is a fast and robust imputation software based on matrix completion, 
called FRMC, for single cell RNA-Seq data. When dealing with single-cell datasets with more than 100k cells, 
please use FRMC_LowMemoryVersion.py.


* **FRMC_LowMemoryVersion.py**  
Version LowMemoryVersion 2.0.0

This LowMemoryVersion is an upgrade to the previous version v1.0.0, which significantly 
reduces the memory requirements and can be used to handle large-scale single-cell datasets 
with 100k or even more cells.


## Requirements
* Python 3
* numpy
* scipy
* matplotlib
* sklearn
* pandas
* math


## Usage
```
python FRMC.py  <Normalized Matrix csv file(cell-by-gene) >  <out_prefix> 

python FRMC_LowMemoryVersion.py  <Normalized Matrix csv file(cell-by-gene) >   <out_prefix> 

```

