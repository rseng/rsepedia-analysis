
# fMLC

fMLC is the official implementation of the MultiLevel Clustering (MLC) algorithm decribed in [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) , used to cluster massive DNA sequences. fMLC was initially implemented by Szaniszlo Szoke and further developed by Duong Vu. It is written in C++ and supports multi-threaded parallelism. fMLC is also integrated with an interactive web-based tool called [DIVE](https://github.com/NLeSC/DiVE) to visualize the resulting DNA sequences based embeddings in 2D or 3D. The work is financially supported by the Westerdijk Fungal Biodiversity Institute and the Netherlands eScience Center.

# Citation

Please cite the following paper if you are using fMLC:

D Vu, S Georgievska, S Szoke, A Kuzniar, V Robert. fMLC: Fast Multi-Level Clustering and Visualization of Large Molecular Datasets, Bioinformatics, btx810, https://doi.org/10.1093/bioinformatics/btx810 

[Pdf verion](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx810/4747887?guestAccessKey=da7a1811-354a-4445-8084-cae44ccafd6f)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.926820.svg)](https://doi.org/10.5281/zenodo.926820)

## Install

[Windows](https://github.com/FastMLC/fMLC/tree/master/Windows)

[Linux](https://github.com/FastMLC/fMLC/tree/master/Linux)

## Data
There are two datasets available as inputs for fMLC. The "small" dataset contains ~4000 ITS yeast sequences, checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute. This dataset were analyzed and released in [Vu D. et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5192050/). The "large" dataset contains ~350K ITS fungal sequences downloaded from GenBank (https://www.ncbi.nlm.nih.gov/) which was used in [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) to evaluate the speed of MLC.

[Download](http://www.westerdijkinstitute.nl/Download/LargeDatasetOf350KITSSequences.zip) the large demo dataset. 

## Results

After clustering the DNA sequences by fMLC, the groupings of the sequences can be saved as output of fMLC. A sparse (or complete) similarity matrix (in .sim format) can be saved in the folder where the dataset is given, to capture the similarity structure of the sequences. Based on this similarity matrix, the coordiates of the sequences can be computed and saved (in .outLargeVis format) using LargeVis. Finally, a json file containing the coordinates and metadata of the sequences is resided in the folder DiVE/data folder as an input of DiVE to visualize the data. This json file can be used for visualization by external applications as well.The clustering and visualization results of the two datasets can be found at https://github.com/FastMLC/fMLC/tree/master/data.

## Contact person 

Duong Vu (d.vu@westerdijkinstitute.nl)


## References

Bolten, E., Schliep, A., Schneckener, S., Schomburg D. & Schrader, R (2001). Clustering protein sequences- structure prediction by transitive homology. Bioinformatics 17, 935-941.

Edgar, R.C (2010). Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26, 2460-2461.
Paccanaro, P., Casbon, J.A. & Saqi, M.A (2006). Spectral clustering of proteins sequences.  Nucleic Acids Res 34, 1571.

Vu D. et al. (2014). Massive fungal biodiversity data re-annotation with multi-level clustering. Scientific Reports 4: 6837.


# fMLC (Windows version)


## Dependencies : 

- [Boost 1.60.0](http://www.boost.org/users/history/version_1_60_0.html)

- [Eigen 3.2](http://eigen.tuxfamily.org/dox-3.2/)

- [LargeVis](https://github.com/lferry007/LargeVis) 

- [DIVE](https://github.com/NLeSC/DiVE)


## Installation:

The folders [LargeVis](https://github.com/lferry007/LargeVis) and [DIVE](https://github.com/NLeSC/DiVE) and the file BioScience.x64.dll should be put in the same folder where the application file MfcCluster.exe is.	
The folder LargeVis should contain LargeVis.exe as compiled from [LargeVis](https://github.com/lferry007/LargeVis). 


## The main windows:

-Input file path: The path of the fasta file of sequences to be clustered. In this fasta file, a sequence is represented by two lines. The title line starting with character ">" containing multiple information fields separated by the pipe character "|". The first information field is the index of the sequence starting at 1. The second line contains the sequence. Another input file with the same name of the fasta file in the format .title describing the information fields can be given optionally. The two input files Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas and Yeast_CBS_GB_ITS_1seqfor1strain_Species.title in the Working folder can be seen as given examples.

-Reference field: When the input file is given, the drop down list on second position will display all possible fields found between the pipe characters. If the .title file is given, then the drop down list will display the information given in this file, otherwise, it displays the information given in the first sequence of the fasta file.
Output file path: The path to save the clustering result, or to save the prediction of the optimal threshold to cluster the given dataset.

-Algorithms: The algorithm used to cluster. There are three clustering algorithms to be selected for clustering: MLC (multilevel clustering, Vu et al. 2014), CCBC (connected components based cluster, Bolten et al. 2001) and GC (the greedy clustering, Edgar 2010) using multi threads or single threads.

-Thresholds: For CCBC and GC, this is a threshold to cluster the dataset. Its value is in between 0 and 1. For MLC, this is a list of increasing thresholds between 0 and 1. The final threshold of the list is the actual threshold that we want to cluster the dataset with. For example, the thresholds to cluster a dataset can be 0.95;0.98. MLC will first cluster the dataset with the threshold of 0.95, and then the obtained groups will be clustered with the threshold of 0.98.

-Cluster: Cluster the given dataset with the selected algorithm and thresholds. 

-Computing Fmeasure checkbox: Compute Fmeasure  immediately after clustering.

-Compute Fmeasure: Compute Fmeasure based on the clustering result.

-Group #: The group number obtained after clustering.

-Fmeasure: The Fmeasure (Paccanaro et al. 2006) obtained by comparing the clustering result with the classification of the given dataset based on the selected field.

-The bottom grid: The grid at the bottom of the window displays the obtained clusters after clustering.

-Save result in input file: The titles of the sequences in the input file will be extended with the centrality indexes of the groups that the sequences belong to at each level.

-Save result in output file: The clustering result are saved in the output file in tab delimited format.

-Save sparse SM: Save a sparse similarity matrix based on the clustering result.

-Save complete SM: Save a complete similarity matrix.

-From: The lower boundary threshold for the prediction. 

-To: The upper boundary threshold for the prediction.

-Step: The incremental step of the thresholds in the prediction.

-Predict OPT: To find an optimal threshold between the lower and upper boundaries that produces the best Fmeasure for clustering.

-Show final groups: Display only the grouping at the final level of MLC.

-Visualize: Visualize the dataset using DiVE. For this action, a sparse/full similarity matrix is required.

-More options: Open an option windows to modify parameters used in clustering, saving similarity matrix and visualizing.

The options windows:

-Minimum overlap: This parameter is used to recompute the similarity score between two sequences if the overlap obtained by BLAST when aligning them is shorter than this value (see Vu et al. 2014).

-Minimum sequence number for MLC: The minimum sequence number for that MLC can be applied.

-Minimum similarity: The minimum similarity score to be saved for a sparse or full similarity matrix. 

-K-neighbor number: This K-neighbor number parameter set up for LargeVis. It’s default value is 150. The remaining parameters of LargeVis are set as default.

-Visualization dimension: 2D or 3D.




# fMLC (Linux version)


Dependencies :

    ## Dependencies : 

- [Boost 1.60.0](http://www.boost.org/users/history/version_1_60_0.html)

- [Eigen 3.2](http://eigen.tuxfamily.org/dox-3.2/)

- [LargeVis](https://github.com/lferry007/LargeVis) 

- [DIVE](https://github.com/NLeSC/DiVE)

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Installation:

The folders LargeVis and DIVE folders should be put in the same folder where the application file ./cluster is.

Parameters:

-input inputfilename: The fasta input file of DNA sequences. Examples can be found at https://github.com/FastMLC/MLC/tree/master/Working

-title titlefilename: The input file describing what are the properties in the titles of the sequences in the input file (optional).

-output outputfilename: The path of the output file to save the clustering result in tab delimited format.

-algo algorithmname(MLC/CCBC/GC/MLC_ST/CCBC_ST/GC_ST): the algorithm used to cluster the dataset. It can be chosen between MLC (MultiLevel Clustering, Vu et al. 2014), CCBC (the connected components based clustering, Bolten et al. 2001), and GC (the greedy clustering, Edgar 2010) using multiple threads. For single threading, they are MLC_ST, CCBC_ST and GC_ST.

-thresholds t: The threshold used to cluster the dataset. For GC and CCBC, it is a value between 0 and 1. For MLC, it is a list of increasing thresholds between 0 and 1. For example, the threshold to cluster the given dataset using MLC can be given as -thresholds 0.95,0.9913.

-fmeasure p: To evaluate the clustering result compared with the classification based on a property given in the titles of the sequences at the position p.

-saveCSM simfilename: To save a complete similarity matrix for the sequences. This option can be run without clustering, i.e. the thresholds argument is not required.

-saveSSM simfilename: To save a sparse similarity matrix based on the clustering result by MLC.

-minsim s: The minimum similarity score to be saved in the simfilename, default value is 0.5.

-K neighs: The K-neighbor number specified for LargeVis, default value is 150. This parameter is provided to save a sparse similarity matrix, in order to improve the accuracy of LargeVis. 

-m minnumber: the minimum number of sequences for applying MLC, otherwise CCBC will be applied. Default value is 100.

-sim simfilename: To specify a similatiry file for visualization. If the arguments -saveSSM simfilename or -saveCSM simfilename are provided, then this argument is not required for visualization.

-visualize d(3D/2D): To visualize the sequences with 2D or 3D.

-predictOpt t1 t2 s n: To predict an optimal threshold that produces the best Fmeasure for clustering, where t1 is the starting threshold, t2 is the end threshold, s is an incremental step of the prediction, and p is the position of the property given in the sequence title to classify the sequences.

Examples:

To cluster:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -title Yeast_CBS_GB_ITS_1seqfor1strain_Species.title -algo MLC -thresholds 0.95,0.9913 -output result.txt

To compute Fmeasure:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -title Yeast_CBS_GB_ITS_1seqfor1strain_Species.title -algo MLC -thresholds 0.95,0.9913 -output result.txt -fmeasure 2

To save a sparse similarity matrix:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -algo MLC -thresholds 0.95,0.9913 -saveSSM yeastSSM.sim

To save a complete similarity matrix:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -saveCSM yeastCSM.sim

To visualize the dataset:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -sim yeastSSM.sim -visualize 3D

To cluster and then visualize:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -algo MLC -thresholds 0.95,0.9913 -saveSSM yeastSSM.sim -visualize 3D

To predict an optimal threshold to cluster the dataset using CCBC:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -algo CCBC -predictOpt 0.98 1 0.001 2

To predict an optimal threshold to cluster the dataset using MLC:

./cluster -input Yeast_CBS_GB_ITS_1seqfor1strain_Species.fas -algo MLC -predictOpt 0.95,0.98 1 0.001 2



There are two datasets available as inputs for fMLC. The "small" dataset contains ~4000 ITS yeast sequences, checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute. This dataset were analyzed and released in [Vu D. et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5192050/). The "large" dataset contains ~350K ITS fungal sequences downloaded from GenBank (https://www.ncbi.nlm.nih.gov/) which was used in  [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) to evaluate the speed of MLC. 

The clustering results of the two datasets can be found at https://github.com/FastMLC/fMLC/tree/master/data/ClusteringResults.

The visualization results (the coordinates computed by LargeVis and the embedded files in .json format) of the two datasets can be found at https://github.com/FastMLC/fMLC/tree/master/data/VisualizationResults. For the small dataset, a sparse similarity matrix is also given. For the large dataset, the similarity matrix cannot be uploaded due to the big size. The embedded files (in .json format) can be uploaded directly to https://nlesc.github.io/DiVE/ for a quick inspection. Manual for DiVE is at https://github.com/NLeSC/DiVE . 
