# Cross-Validated Rank-based Marker Selection (CVRMS)

CVRMS is an R package designed to extract marker subset from ranked marker dataset induced by genome-wide association study (GWAS) or marker effect analysis for genome-wide prediction.

## Requirement
 * R (>= 3.5.1)
 * pre-installed library in R
   - argparse
   - ggplot2
   - rrBLUP
   - ggpmisc
   - caret
 * How to install libraries (in R console)
   > install.packages(c("argparse", "ggplot2", "rrBLUP", "ggpmisc", "caret"))
   
 ## Quick example (in linux, ios, windows command console)
  * Extract rice_geno.zip file
  
  $ Rscript CVRMS_v1.3.R -g rice_geno.txt -p rice_pheno.txt -pn 1 -gw gwas_result.txt -gw_snp_cn 1 -gw_chr_cn 2 -gw_pos_cn 3 -gw_pv_cn 4 -min 10 -max 300 -cv 5 -a 0.9 -d 0.001 -ss 1 -m rrblup -t 1
  
  $ cd phenotype1
  
  $ Rscript Create_final_model_v1.0.R snp_info.txt
  
## Detail options
 -g : genotypes input (row : markers, column : samples) - first row is the header, the first column includes marker names and the sample names should be included from the second column. !!Caution - genotypes are encoding to 0, 1, 2!!
 
 -p : phenotypes input (row : samples, column : phenotypes) - the first column includes the samples names and the traits of interests are from the second to end column)
 
 -pn : phenotype column number
 
 -gw : GWAS result file (e.g. plink, GAPIT, etc.)
 
 -gw_snp_cn : column number of SNPs in GWAS result file
 
 -gw_chr_cn : column number of chromosome in GWAS result file
 
 -gw_pos_cn : column number of position of SNP in GWAS result file
 
 -gw_pv_cn : column number of p-value of SNP in GWAS result file
 
 -min : Minimum number of markers user wants
 
 -max : Maximum number of markers user wants
 
 -cv : how many subsets for cross-validation in grid search
 
 -a : Goal of prediction accuracy (correlation coefficient - regression, Correct decision rate - classification)
 
 -d : minimum increasing rate of prediction accuracy for selecting markers
 
 -ss : how many selected markers for each iteration
 
 -m : prediction method (rrblup and rf are only available)
 
 -t : Limitation of time (hour) for iteration
 
 -SNP_info.txt : SNP information file similar with plink map file which have to contains SNP name, Chromosome, Position (BP).
 
 ## Contact
 
 lovemun@kribb.re.kr
 
 jyoon@kribb.re.kr
