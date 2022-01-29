# DeepAc4C
DeepAc4C, which identifies ac4C using convolutional neural networks (CNNs) using hybrid features composed of physico-chemical patterns and a distributed representation of nucleic acids.
## Webserver and datasets
A webserver is available at: http://lab.malab.cn/~wangchao/softs/DeepAc4C/.

The source code and datasets(both training and testing datasets) can be freely download from the github and the webserver page.

## Brife tutorial

### 1. Environment requirements
Before running, please make sure the following packages are installed in Python environment:

gensim==3.4.0  
pandas==1.0.3  
tensorflow==2.3.0  
python==3.7.3  
biopython==1.7.8  
numpy==1.19.2  

For convenience, we strongly recommended users to install the Anaconda Python 3.7.3 (or above) in your local computer.

### 2. Running
Changing working dir to DeepAc4C, and then running the following command:  
python DeepAc4C.py -i ./sequence/input_query.fasta -o prediction_results.csv  
-i: input file in fasta format  
-o output file name 

### 3. Output
The output file (in ".csv" format) can be found in results folder, which including sequence name, predicted probability and redicted result.  
Sequence with predicted probability > 0.5 was regared as ac4C site.  

### 4. References
Chao Wang et al. 2021. DeepAc4C: A convolutional neural network model with hybrid features com-posed of physicochemical patterns and distributed representation information for identification of N4-acetylcytidine in mRNA. Bioinformatics (Accepted).  

