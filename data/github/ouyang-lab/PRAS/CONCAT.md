# Examples
All sample datasets can be downloaded from the [Downloads](zipped_code/Downloads.md) page.
## Sample input
PRAS requires a binding table, an annotation table, and an ID table as input files (described in the [Instructions](Instructions.md)).
Here, we present a set of toy datasets as an example for PRAS usage.
The binding table is a subset of the reliable cross-linking sites from CELF4 iCLIP-seq dataset used in our paper. The annotation table is the GTF file of Refseq mm9 annotation. The ID file is the transcript IDs and gene names of Refseq mm9 annotation.
The files looks like the screenshots as below:
### 1. Sample binding table: test_binding_tab.bed

![alt text](figures/sample_binding.png)

Note that the fifth column could have sign to indicate the strand.
### 2. Sample annotation table: mm9_refseq.gtf

![alt text](figures/sample_annotation.png)

### 3. Sample ID table: mm9_refseq_ID.txt

![alt text](figures/sample_ID.png)

In addition, PRAS support calculating the scores on a customized gene list via -l option. If it's not given, PRAS will perform the calculation on all the genes in the annotation files.
## Command example
### 1. run PRAS "check" mode to examine the binding profiles around the candidate reference sites.
python PRAS_1.0.py -g mm9_refseq.gtf -t mm9_refseq_id.txt -i test_binding_tab.bed -m check -s transcript -a test_binding_tab_assign.txt -w 10 -c 500 -p check_test
### 2. run PRAS with known reference site.
python PRAS_1.0.py -g mm9_refseq.gtf -t mm9_refseq_id.txt -m score -s 3UTR -i test_binding_tab.bed -a test_binding_tab_assign.txt -w 10 -r 3 -d 1000

## Sample output
### 1. Sample output table 1: test_binding_tab_utr3.assign.txt
In "check" mode, peaks in 5' UTR, CDS and 3' UTR are all annotated. In "score" mode, only peaks in selected genomic region are annotated.
The format of the annotated peak file is consistent for the two modes, which is shown as below.

![alt text](figures/sample_assign.png)

### 2. Sample output plot: check_test_500nt_around_TSS_TIS_TTS_TES.pdf
If PRAS "check" mode is enabled, the user will get the binding profile plot around the candidate reference sites.
Here is the output plot of the toy dataset when -s transcript is enabled.

![alt text](figures/check_profile.png)

We can see the 3'UTR enrichment from the binding profile plots.
In addition, PRAS also allows splice site binding check via -s splice option. Users can try it using RBPs with splicing-related function.
### 3. Sample output suggestion file: check_test_500nt_PRAS_option_suggestions.txt
If PRAS "check" mode is enabled, the user will get the parameter suggestion file based on the binding profile plot around the candidate reference sites.
These parameters can be used in the "score" mode if there is no pre-knowledge of the target RBP.

![alt text](figures/suggestion_file.png)

### 4. Sample output table: binding_PRAS.txt
Given known reference site, PRAS will generate a score table with four columns which are gene name, PRAS score, log2 of PRAS score, and the rank of the gene.

![alt text](figures/sample_score.png)

One will get something similar to the following lines if PRAS run successfully.

If PRAS "check" option is enabled, the process is as below:

![alt text](figures/check_time.png)

If PRAS is run with known reference site in "score" mode, the process is as below:

![alt text](figures/score_time.png)

# PRAS
PRAS: predicting functional targets of RNA binding proteins based on CLIP-seq peaks

We designed a tool, Protein-RNA Association Strength (PRAS), for RBP functional targets prediction.
Given the binding  (from CLIP-seq) and reference site (pre-known or interested site) information, it can score genes by combining the binding intensity and position information.
PRAS outputs the PRAS score of each interested gene.

![alt text](figures/flowchart1.png)

## [Instructions](Instructions.md)
1. Prerequisites
2. Run PRAS

## [Downloads](zipped_code/Downloads.md)
1. Source code
2. Sample data

## [Examples](Examples.md)
1. Sample input
2. Command example
3. Sample output

## Reference
[Lin, J., Zhang, Y., Frankel, W. N., & Ouyang, Z. (2019). PRAS: Predicting functional targets of RNA binding proteins based on CLIP-seq peaks. PLoS computational biology, 15(8), e1007227.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007227)

## License
Use of PRAS is free for academic users under the GNU General Public License (GPL). Commercial users please contact the authors.

## Contact
Zhengqing Ouyang: ouyang@schoolph.umass.edu
Jianan Lin: jianan.lin@jax.org

# Instructions
# Prerequisites
## 1. Prepare the data files.
The three required input files of PRAS are the binding file (from CLIP-seq), the annotation file and the ID file.
### a. The binding file.
The binding file has to include the intensity and position information of the binding sites in the 6-column BED (Browser Extensible Data) format. It should be a file with either the reliable cross-linking sites or the reliable binding peaks called by any peak calling algorithm. The fifth column of the BED file should be the intensity of the binding site, read counts for example, and the sixth column should be the strand information, "+" or "-". Consistent with the standard BED file, the columns of the binding file should be tab-seperated. The columns information are shown below:

chr   start   end   read_name    read_counts   strand
...

As for the reliable cross-linking sites or the peak regions, users can generate them from any peak-calling tool. For example, the reliable cross-linking sites in our study (iCLIP-seq) were generated from iCount (http://icount.fri.uni-lj.si) with FDR cutoff as 0.05.
For detailed example of the data files, please refer to Examples.
### b. The annotation file.
The annotation file should be formatted as GTF (Gene Transfer Format). We suggest to download the GTF file of certain genome from the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables).
### c. The ID file.
Considering that some of the downloaded GTF file doesn't include the gene name information, we require this additional file with transcript ID and gene name to be linked together. There should be two columns in this file, where the first column is the transcript ID, and the second column is the gene name or gene ID.
One easy way to generate such file is from the UCSC Table Browser. Here are the steps:
1st, select the genome build you want, mm9, for example;
2nd, select the annotation source you want, Refseq, for example;
3rd, select "selected fields from primary and related tables" in the output format option;

![alt text](figures/step3.png)

4th, check the boxes of "name" and "name2" from the table;
5th, click the "get output".

![alt text](figures/step4.png)

## 2. Download and install necessary tools.
### a. Bedtools
Bedtools can be downloaded from: http://bedtools.readthedocs.io/en/latest/
Bedtools installation instruction can be found from: http://bedtools.readthedocs.io/en/latest/content/installation.html
### b. Python
Python can be downloaded from: https://www.python.org/downloads/
### c. R
R can be downloaded from: https://www.r-project.org/
Bedtools commands, python and R should be executable.
### d. gtfToGenePred toolkit.
gtfToGenePred can be downloaded from UCSC utility library (http://hgdownload.soe.ucsc.edu/admin/exe/). This tool can transfer the annotation file in GTF format to that in GenePred format. PAS calls this tool in the annotation generating step, so please make it excutable.

# Run PRAS
## 1. [Download the source code of PRAS.](zipped_code/Downloads.md)
PRAS source code can be downloaded from the [Downloads](zipped_code/Downloads.md) page.
## 2. Check help page of PRAS for usage.
PRAS help page can be generated simply by typing:

PRAS_1.0.py -h

in your command window. And you should get the following lines:

[username]$ python PRAS_1.0.py -h
![alt text](figures/helppage.png)

If you get any errors in this command line, please go back to the prerequisites and check.
## 3. Run PRAS on your own datasets.
PRAS have four required arguments: GTF file name, Binding file name, ID file name and the filename for the output peak annotation file.
PRAS have seven optional arguments: running mode, genomic region, half window size, interested gene list, reference site direction (5', 3' or check), distance parameter, and PRAS score output filename.
One can follow the help page of PRAS to run on your own datasets.
If you have problems running the python code, please go to [Examples](Examples.md) to check the sample command line to run PRAS.
## 4. Collect PRAS output files.
If PRAS is run in 'check' mode instead of 'score' mode, it will output two files. One is a peak annotation file according to the -a option, and the other is the binding profile plot according to the -p option. 'check' will output parameter suggestions that can be used in 'score' mode based on the CLIP peak distribution around the selected reference sites. The suggested parameters include genomic region, reference site, and decay parameter, corresponding to -s, -r, and -d options, respectively. Users can use these suggestions if they don't have any pre-knowledge of the studied RBP.
Given known reference site, PRAS typically will output two files. One is a peak annotation file according to the -a option, and the other is the PRAS score table according to the -o option.
PRAS can work on a set of interested genes provided by a gene list customized by user via -l option. If there's no interested gene list as input, PRAS will take all the genes from the annotation file as the gene list and at the same time PRAS will generate a file called "_all_genels.txt", which can be deleted if it's not being used in the following analysis.
# Description of the files in this folder
They are images used in the PRAS/README.md for instructions of the tool.
# Downloads
## Source code
[PRAS source code.](script.zip)
## Sample data
[PRAS sample datasets.](sample.zip)
