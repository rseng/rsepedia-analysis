# Code for the SparsePro Paper

This repo contains simulation pipelines and scripts for generating figures used in the [SparsePro paper](https://doi.org/10.1101/2021.10.04.463133).

## Softwares

We have used the following softwares in our study.

* [GCTA](https://cnsgenomics.com/software/gcta/#Download) (version==1.93.2beta)

* [FINEMAP](http://www.christianbenner.com) (version==v1.4)

* [SuSiE](https://stephenslab.github.io/susieR/) (version==0.11.50)

* [PolyFun](https://github.com/omerwe/polyfun) (version==1.0.0)

* [SparsePro](https://github.com/zhwm/Sparse_Pro) (version==1.0.0)

* Python and Python modules:

  - [Python](https://www.python.org) (version==3.9.7)
  - [numpy](http://www.numpy.org/) (version==1.21.3)
  - [scipy](http://www.scipy.org/) (version==1.7.1)
  - [pandas](https://pandas.pydata.org/getpandas.html) (version==1.3.4)
  - [PyTorch](https://pytorch.org) (version==1.7.1)
  - [rpy2](https://rpy2.github.io/doc/latest/html/index.html) (version==3.4.5)

* R and R packages:

  - [R](https://www.r-project.org) (version==4.0.0)
  - [ggplot2]() (version==3.3.5)
  - [PRROC](https://cran.r-project.org/web/packages/PRROC/index.html) (version==1.3.1)

## Data

1. UK biobank LD information provided by PolyFun can be downloaded from [here](https://alkesgroup.broadinstitute.org/UKBB_LD/).

2. Individual level phenotype and genotype data from the UK Biobank are available upon successful application to its research committee.

3. Annotation files can be downloaded from the [baselineLF_v2.2 model](https://alkesgroup.broadinstitute.org/LDSCORE/).

## Simulation

We conducted simulations to showcase the efficiency and utility of our method. The pipeline and codes used for simulation are in the [sim](sim/) directory. Here we include a brief description of the simulation procedures:

1. Sample 50 causal SNPs per chromosome amongst genetic variants with MAF>0.001 and INFO>0.6;

```
python ../gen_Cidx.py --anno ../UKBB.UCSC.22.anno --save ./ --K 50 --h2 0.01 --aw 2.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 2.0 0.0 --seed GENO
```

2. Simulate trait with 1% heritability explained by these causal SNPs using GCTA;

```
~/utils/gcta_1.93.2beta/gcta64 --bfile ../UKB.22.snp --simu-qt  --simu-causal-loci causal --simu-hsq 0.01 --simu-rep 1 --out test --threads 1
```

3. Perform GWAS with GCTA-fastGWA to obtain summary statistics;

```
~/utils/gcta_1.93.2beta/gcta64 --bfile ../UKB.22.snp --pheno test.phen --fastGWA-lr --out test --threads 1
```

4. Perform GCTA-COJO to identify lead SNPs;

```
~/utils/gcta_1.93.2beta/gcta64 --cojo-file test.fastGWA.ma --cojo-slct --bfile ../UKB.22.snp --out test --threads 1
```

5. Run FINEMAP on loci defined by lead SNPs identified by GCTA-COJO;

```
~/utils/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --in-files master --sss > run_finemap_no.log
```

6. Run SuSiE with an in-house script using the "susie\_rss" function on sliding windows;

```
python ../susie.py > run_susie.log
```

7. Run SparsePro on the same sliding windows with var_Y estimated from summary statistics;

```
python ~/scratch/SparsePro/sparsepro_all.py --ss ss --N $(cat N) --LDdir /home/wmzh22/scratch/SparsePro/UKBBLD/ --save 5/ --CHR 22 --K 5 --var_Y $(cat var_Y) > run_sp_5.log
```

8. Repeat steps 1-7 22 times to generate whole-genome summary statistics; 

```
for i in $(seq 1 22); do sed "s/GENO/$i/g" gen_sim.sh > gen_sim_$i.sh; done
```

9. Run PolyFun to obtain per-SNP prior probability of being causal;

```
python ~/utils/polyfun/polyfun.py --compute-h2-L2 --output-prefix test --sumstats poly_ss.parquet --ref-ld-chr ~/scratch/SparsePro/POLY/polysim/fake. --w-ld-chr ~/scratch/SparsePro/POLY/polysim/weights.
```

10. Estimate relative enrichment of annotations and update prior probability of being causal with SparsePro;

```
python ~/scratch/SparsePro/sparsepro_all_plus.py --ss ss --var_Y $(cat var_Y) --N $(cat N) --anno ~/scratch/SparsePro/idpANNO/UKBB.UCSC.22.anno --LDdir ~/scratch/SparsePro/UKBBLD/ --save 5 --K 5 --CHR 22 --W ../W_sig_5.csv > run_asp_5.log
```

11. Run FINEMAP with PolyFun-informed prior;

```
~/utils/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --in-files master --sss --prior-snps > run_finemap_poly.log
```

12. Run SuSiE with PolyFun-informed prior;

```
python ../psusie.py --snpvar ../test_polyfun/test.GENO.snpvar_ridge_constrained.gz > run_psusie.log
```

## Plotting figures

Data and scripts for generating Figures 2-4 are deposited in the [plt](plt/) directory.

1. Figure 1 was created with PowerPiont;

2. Figure 2:

```
<in RStudio> plot_Fig2.R
```

3. Figure 3: 

```
<in RStudio> plot_Fig3.R
```

4. Figure 4 was generated with [PhenoGram](http://visualization.ritchielab.org/phenograms/plot) and [LocusZoom](http://locuszoom.org/).




