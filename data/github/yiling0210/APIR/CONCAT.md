# APIR: flexible and powerful FDR-control framework for aggregating peptides identified by multiple database search algorithms from mass spectrometry data

Yiling Elaine Chen, Kyla Woyshner, MeiLu McDermott, Leo David Wang, and Jingyi Jessica Li

## Latest News!
> 2021-04-27
- APIR is now online!

## Introduction and tutorial
APIR is a statistical framework that aggregates peptide-spectrum-matches (PSMs) identified by multiple database search algorithms with FDR control. APIR takes as input a target FDR threshold (such as 0.05) and the PSM level output from the search algorithms of interest and outputs a PSM level, a peptide level, and a protein level results file. 

Users can refer to [our manuscript](XXX) for a detailed description of the method and applications. For a detailed user manual, please see [User Manual](XXX).

Any suggestions on the package are welcome! For technical problems, please report to Issues. For suggestions and comments on the method, please contact Yiling (yiling0210@ucla.edu) or Dr. Jingyi Jessica Li (jli@stat.ucla.edu).

## Installation
The package is not on CRAN yet. For installation please use the following codes in R
```
if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("yiling0210/APIR")

```

