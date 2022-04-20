# CircleBase
An integrated resource and analysis platform for human eccDNAs. Also see [homepage](http://circlebase.maolab.org/)


# Scoring system
### Dependencies
- bedtools 2.0 or higher [doc](http://bedtools.readthedocs.io/)
- python 3.7 [www.python.org](https://www.python.org/)
- numpy [www.numpy.org](http://www.numpy.org/)
- scipy [www.scipy.org](https://www.scipy.org/)
- other common packages: *multiprocessing* and *argparse*

tips: Anaconda is always a good choice to install the dependencies.

### Input files (bed file for six regulatory categories and eccDNAs)
1. Chromatin_access.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Chromatin_access.bed.gz)
2. Chromatin_interaction.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Chromatin_interaction.bed.gz)
3. Epigenetic_regulation.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Epigenetic_regulation.bed.gz)
4. Genetic_variant.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Genetic_variant.bed.gz)
5. Regulatory_elements.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Regulatory_elements.bed.gz)
6. Targeting_genes.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/Targeting_genes.bed.gz)
7. eccDNA_core.hg19.bed [download](http://159.226.67.237/sun/oncobase/assets/data/download/eccDNA_core.hg19.bed.gz)

### How to run
1. Go to the *scoring system* directory and set all the dependencies
2. Download all the input files listed above and decompress them
3. Run the *run.sh* shell script

### Output
- hits.stat.* files are annotated hits (records) count for each eccDNA in six regulatory categories. The last field is the count number.
- *.score files inlcude score for each eccDNA corresponding Gaussian mode in six regulatory categories. Here are the fields:
1. eccDNA id.
2. Chromosome to which the eccDNA belongs.
3. Hits number for the eccDNA.
4. Hits number after Box-Cox transformation for the eccDNA.
5. Mean of the hits number for all eccDNAs at chromosome list on the second field (i.e., 𝜇 of the Gaussian distribution).
6. Standard Deviation of the hits number for all eccDNAs at chromosome list on the second field (i.e., 𝜎 of the Gaussian distribution).
7. Probability greater than the hits number in the corresponding Gaussian  distribution.
8. The score for the eccDNA (i.e., negative of the base 10 logarithm of the Probability).
- *.nor files include normalized score of each category. The first 8 columns are same as *.score files, column 9 is the Z-score of the regulatory category and column 10 is the normalized score.
 
- final.score.txt file is the final result we want. Here are the fields:
1. eccDNA id.
2. Average of normalized scores for all six regulatory categories.

You can download all output [here](http://circlebase.maolab.org/assets/download/score_system_result.tar.gz)

