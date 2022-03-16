
####  T7-Adapter_trimming
Python version: Python 3.7.4

Packages required: Bio.Seq, fuzzywuzzy, Levenshtein, multiprocessing, pandas

Just keep all parameters as default

- For single-end data (for R1):
```shell
	Usage: python adapter_trim.T7.SE.py -i R1.fastq.gz
	Usage: python adapter_trim.T7.SE.multiprocessor.py -i R1.fastq.gz -p 5
```

- For pair-end data (R1 and R2)
```shell
	Usage: python adapter_trim.T7.PE.py -r1 R1.fastq.gz -r2 R2.fastq.gz
	Usage: python adapter_trim.T7.PE.multiprocessor.py -r1 R1.fastq.gz -r2 R2.fastq.gz -p 5
```

#### Pearson correlation between two different samples
```shell
cat $T7_peak1 $T7_peak2 $normal_peak1 $normal_peak2 |sort -k1,1 -k2,2n |bedtools merge -i stdin > ${out}_normal.merge.bed
bedtools multicov -bams $T7bam1 $T7bam2 $normal1 $normal2 -bed ${out}_normal.merge.bed > ${out}_normal.merge.bed.cov
Rscript Pearson_correlation.r ${out}_normal.merge.bed.cov ${out} Normal
```
####  mapping and peak calling
```shell
bash fastq2peak.sh 
```

#####  annotation for peaks

```shell
Rscript Genomic_annotation.r
Rscript relativetotss_annotation.r
```
#### superenhancer calling
```
bash Super_enhancer.sh
```

#### library complexity figure plot
```
python preseq_plot.py output
```

```shell
bedtools multicov -bams FA-HG-A5000-ac-1.q2.sort.bam FA-HG-A5000-ac-2.q2.sort.bam FA-HG-C5001-ac-1.q2.sort.bam FA-HG-C5001-ac-2.q2.sort.bam HCRC-1132.rep1.q2.sort.bam HCRC-1132.rep2.q2.sort.bam HGBM-A1-rep1.q2.sort.bam HGBM-A1-rep2.q2.sort.bam HGBM-C1-ac-1.q2.sort.bam HGBM-C1-ac-2.q2.sort.bam  -bed all.bed > all.merge.bed.cov
Rscript clusterbynmf.r all.merge.bed.cov
```
##### genes were divided into four grous by its expression 
```
Rscript Quantile_by_gene_expression.r
```



