# deMeta
Reverse the Meta-analysis mostly used by GWAS.

A typical scenario is for using GWAS consortium meta-analysis results for a contributing group.
For example, a group has performed GWAS on its own samples, and contributed the GWAS summary statistics to a Consortium. If the group wants to perform polygenic risk prediction on it sample, it will need the summary statistics from all contributing groups except its own. Two options may achieve this goal, 

 1. The analysts of the consortium can re-do the meta exlcuding the sample from requesting group; 

 2. the consortium share sub-study summary statistics with each contributing groups.

The first option will incur heavey workload to the analysts if more than 1 substudies wanted to be removed. The second option may sometimes not practical due to data regulation.

Other cases during daily research may also poped up. For example, you have two summary statistics from two GWAS, but one include another, you want to computed genetic correlation of the two. In this case, you need to computer the non-overlapping suammary statistics for the large GWAS.

# Getting Started
- Clone this repository using the following git command:

  `git clone https://github.com/Computational-NeuroGenetics/deMeta`

- Dependencies:
  
  pandas
  
  numpy
  
  matplotlib
  
  scipy
  
# Using deMeta

  Run `$python Subtract_meta.py -h or ./Subtract_meta.py --help` 
  
  This will give

 `usage: Subtract_meta.py [-h] --masf MASF [--masA1 MASA1] [--masA2 MASA2]
                     [--masEff MASEFF] [--masOR MASOR] [--masisORSE]
                     [--masP MASP] [--masSE MASSE] [--masV MASV]
                     [--masCHR MASCHR] [--masPOS MASPOS] [--masN MASN]
                     [--masZ MASZ] [--masSNP MASSNP] --ssf SSF [--ssA1 SSA1]
                     [--ssA2 SSA2] [--ssEff SSEFF] [--ssOR SSOR] [--ssP SSP]
                     [--ssSE SSSE] [--ssisORSE] [--ssV SSV] [--ssN SSN]
                     [--ssZ SSZ] [--ssSNP SSSNP] [--top1 TOP1] [--top2 TOP2]
                     [--flip] [--noIVW] --out OUT`        
```
  From meta-results remove one contributing study.

    Applicable to:
        1. Inverse variance weighted meta-analysis
        2. Sample size weighted meta-analysis
    Figures:
        1. Manhattan plots for before and after removing the sub-study
        2. QQ-plots for before and after removing the sub-study
    Notes:
        1. Only common SNPs in the two file will be analyzed. The SNPs that
        exist only in the original meta-analysis results should be added back.
```

Optional arguments:
 ```
   -h, --help       show this help message and exit
  
  --masf MASF      (required) Meta-analysis result file name
  
  --masSNP MASSNP  (required) Meta-analysis result SNP column name
  
  --masA1 MASA1    (required) Meta-analysis result effect allele column name, default='A1'
  
  --masA2 MASA2    (required) Meta-analysis result the other allele column name, default='A2'
  
  --ssf SSF        (required) Sub-study result file name
  
  --ssSNP SSSNP    (required) sub-study result SNP column name
  
  --ssA1 SSA1      (required) sub-study result effect allele column name
  
  --ssA2 SSA2      (required) sub-study result the other allele column name
  
  --flip           (optional) whether flip strand, using meta-analysis result as reference
  
  --out OUT        (required) Result file prefix (required)
  ```
  
- Sepcific arguments for sample size weighted meta-analysis: 
``` 
  Either effect or Odds ratio should be given for meta and sub-study, respectively

  --masEff MASEFF  (conditional) Meta-analysis result effect (of A1) column name, default='Beta'
  
  --masOR MASOR    (conditional) Meta-analysis result Odds ratio (of A1) column name, default='OR'
  
  --ssEff SSEFF    (conditional) sub-study result effect (of A1) column name, default='beta'
  
  --ssOR SSOR      (conditional) sub-study result Odds ratio (of A1) column name, default='or'
  
  Either SE of effect/Odds ratio or variance should be given for meta and sub-study, respectively
  
  --masisORSE      (conditional) Is the Meta-analysis result SE on OR scale, default on ln(OR) scale
  
  --masSE MASSE    (conditional) Meta-analysis result standard error (of Beta) column name, default='SE'
  
  --ssSE SSSE      (conditional) sub-study result standard error (of Beta) column name, default='se'
  
  --ssisORSE       (conditional) Is the sub-study result SE on OR scale, default on ln(OR) scale

  --masV MASV      (conditional) Meta-analysis result variance (of Beta) column name, default='VAR'
  
  --ssV SSV        (conditional) sub-study result variance (of Beta) column name, default='var'
```
 
- Sepcific arguments for sample size weighted meta-analysis
```
  --masN MASN      (conditional) Meta-analysis result Sample size column name, default='N'
  
  --masZ MASZ      (conditional) Meta-analysis result Z score column name, default='Zscore'
  
  --ssN SSN        (conditional) sub-study result Sample size column name, default='n'
  
  --ssZ SSZ        (conditional) sub-study result Z score column name, default='zscore'
  
  --noIVW          (required) whether meta-analysis result is inverse variance weighted? Otherwise using sample size weighted
```

-  Arguments for Manhattan and Q-Q plot
```
  --masCHR MASCHR  (required) Chromosome number column name in original Meta-analysis results, required for Manhattan plot, default='CHR'
  
  --masPOS MASPOS  (required) Genomic position column name in original Meta-analysis results, required for Manhattan plot, default='POS'
   
  --masP MASP      (required) Meta-analysis result p value (of A1) column name, required for Manhattan plot, default='P'

  --ssP SSP        (required) sub-study result p value (of A1) column name, default='p'
  
  --top1 TOP1      (optional) max -log10(P) for original meta-analysis to plot in the manhattan (default from the data)
  
  --top2 TOP2      (optional) max -log10(P) for subtracted results to plot in the manhattan (default from the data)
 ```
# Output

Output files:

1. Obtained summary statstics

2. Manhattan plot for before and after removing the sub-study

3. QQ-plot for before and after removing the sub-study

# Example

As a demonstration we applied deMeta to the summary statistics of the GWAS for body mass index (BMI) (Locke, et al., 2015; Yengo, et al., 2018). The Yengo et al data is a meta-analysis results of the Locke et al study and the UK biobank data. Both the Yengo et al and the Locke et al data were ware downloaded from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files. 

The inverse of the inverse variance weighted function of deMeta was applied to obtain summary statistics for UK biobank data.

- QQ-plot for before and after removing the UKB samples from the GIANT BMI GWAS studies 
![Image](../master/test/BMI_qq.png?raw=true)

- Manhattan plot for before and after removing the UKB samples from the GIANT BMI GWAS studies 
![Image](../master/test/BMI_manhattan.png?raw=true)

