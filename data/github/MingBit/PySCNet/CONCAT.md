# PySCNet
There are four modules:
1) Proprecessing: Normalization, Quality Control, FeatureSelection, Batch Correction???
2) BuildNet: WGCNA_based, Bayesian Network based, PICD, GENIE3...
3) NetEnrich: Network traverse, network merge...
4) Plots

20.April
1) Test with wgcna function
2) Test with GENIE3 function
3) Test with BNN function -->> improve
...

21.April
1) gene network enrichment / comparison
2) check sc RNA-seq data for immune cells
3) (SNA): Degree Centrality, Betweeness Centrality,Closeness Centrality, Eigen Centrality, Page Rank...
4) FeatWalk? RandomWalk? Breadth First Search (BFS), Deepth First Search (DFS)...

24.April
1) Random Walk + supervised Random Walk
2) bnlearn
3) PICD...

28.April
1) Network merge by union, intersection and de bruijn mapping
2) Test with first two data set.

29.April
1) multiple Random walk with cycles??? deep walk???
2) De bruijn mapping: hamiltonian path, Eulerian path; Overlap layout consensus mapping
3) Network similarities/matching ranking
4) GMatch4py???
5) Similarity network fusion (SNF)

30.April
1) Test with two simulated datasets
2) Visulization plots
3) graph embedding based method for network merge

08.May
1) modify genie function
2) test with run30: whole / each cluster

10.May
1) test with TF factors for run30 and bulk seq
2) dynamic network function
3) gene-trigger paths test for run30 and run36

13.May-19.May: I was ill during this week...=_=

20.May
1) plotting functions
2) Gnetdata: GNR sepecific object

22.May
1) Learn docker

23.May
1) No progress
2) A wonderful day for me: two cats (Baby and Lumi) came to my life.

24.May
1) update BuildNet with Docker app

28.May
1) Horrible day: I let my two cats out and they are still not coming back.

29.May - 12. June
1) minor progress

13.June
1) Update CORR method
2) test with SJARACNe method

14.June
1) Call function for docker container
2) NetEnrich module

16.June
1) Complete NetEnrich module
2) Test with in-house data

17.June
1) Visualization module

18.June
1) test with in-house dataset

26.June
1) gnetdata save via pickle
2) knn test
3) de bruijin test
4) simulated dataset with merge function
5) Function details

29.June
1) Add other methods to BuildNet module
2) Add functions details
3) Upload to PyPI and readdocio

1.July
1) Update SINCERA

2.July
1) Ensemble classifier
2) BNN classifier
3) Test with old dataset

3.July:
1) Scanpy-based preprocesing
2) affinity network fusion

8.July:
1) Random walk plots
2) poster start

9.July:
1) Visualization!!!
2) Simulated data/published data with Random walk.

5.August
1) After workshop and conference, finally back to work
2) R shiny app for visualization

8.August
1) update the structure of gNetData
2) re-do BuildNet ect and visualization

13.August
1) gene module detection in R
2) publish R shiny

14.August
1) fix the shiny-web bug
2) network with multi-nodeshape

19.August
1) had to switch to side projects...:(

20.August
1) Again! locally publish R shiny:
2) Go to deep learning for GRN

21.August
1) collect one more published GRN method
2) test seq2seq model

22.August
1) Switch to side project

23.August
1) pre-processing module via scanpy

26.August
1) finish work from last week. :D

29.August
1) switch to PNAS and NC papers. :D
2) update NN knowledge

2.Sept
1) re-do GRN plot with smart gene features

3.Sept
1) Optimize Random walk
2) Optimize Rshiny

6.Sept
1) Finish work from last time :(

9.Sept
1) Ideas: Compressive sensing and sparse coding or encoding
2) Re-think about reconstructing GRNs from sparse data. Not must be deep learning.

10.Sept - 24.Sept
1) Ahh! I got lost in my PhD and couldn't come back for this..
2) Not sure if I'll still work on this..

25.Sept
1) fix Docker bug
2) Random walk with re-start

30.Sept
1) try VAE

14.Oct
1) graph embedding python
2) SPT and MST implementation

17.Jul
1) create a pyscnet tutorial
2) update gnetData
3) clean BuildNet module

18.Jul
1) Find consensus GRNs
2) GE for link prediction [![Codacy Badge](https://api.codacy.com/project/badge/Grade/d3c17aac77e14f6bb17b33f875ff7471)](https://app.codacy.com/manual/MingBit/PySCNet?utm_source=github.com&utm_medium=referral&utm_content=MingBit/PySCNet&utm_campaign=Badge_Grade_Dashboard)
[![License](https://img.shields.io/github/license/MingBit/PySCNet)](https://github.com/MingBit/PySCNet/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/MingBit/PySCNet.svg?branch=master)](https://travis-ci.org/MingBit/PySCNet)
[![Documentation Status](https://readthedocs.org/projects/pyscnet/badge/?version=latest)](https://pyscnet.readthedocs.io/en/latest/?badge=latest)

# PySCNet: A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data
There are four modules:
1) **Pro-precessing**: initialize a gnetData object consisting of Expression Matrix, Cell Attributes, Gene Attributes and Network Attributes;
2) **BuildNet**: reconstruct GRNs by various methods implemented in docker;
3) **NetEnrich**: network analysis including consensus network detection, gene module identification and trigger path prediction as well as network fusion;
4) **Visulization**: network illustration.

![Overview](https://github.com/MingBit/PySCNet/blob/master/images/workflow_update.png)

# :tada: :confetti_ball: Create your own GRNs
[Dashboard](https://github.com/MingBit/PySCNet/blob/master/images/pyscnet_dashboard.gif) is available now for creating your own GRNs.
Cell specific GRNs and network analysis results can be saved as a pickle object and upload onto PySCNet-Dashboard.
It provides parameter settings and allows for parameter adjustment and GRNs customization. <br/>
To run the python dashboard: <br/>
`cd PySCNet/pyscnet/dash_pyscnet/` <br/>
`python app.py` 


# Installation
Make sure you have [Docker](https://docs.docker.com/engine/install/ubuntu/) and [graph_tool](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions) manually installed. <br/>
`conda create --name gt -c conda-forge graph-tool python=3.6` <br/>
`conda activate gt`

1) clone from github:
`git clone https://github.com/MingBit/PySCNet`
2) create a new folder and set up:
`mkdir dist | python setup.py sdist`
3) install pyscnet:
`pip install dist/pyscnet-0.0.3.tar.gz`

# Tutorial
Make sure you have [scanpy](https://scanpy.readthedocs.io/en/stable/installation.html) and [stream](https://github.com/pinellolab/STREAM) manually installed. <br/>
`pip install jupyterlab scanpy==1.5.0` <br/>
`conda install -c bioconda stream` <br/>

You might need to re-install anndata: `pip install anndata==0.7.4`

Mouse HSC data preprocessed and analyzed by stream as explained in this 
[tutorial](https://github.com/MingBit/PySCNet/blob/master/tutorial/pyscnet_stream.ipynb). 

open jupyter-notebook with `/miniconda3/envs/gt/bin/./jupyter-notebook ~/PySCNet/tutorial/pyscnet_stream.ipynb `

# TO-DO
1) Add an Auto-ML based pipeline to Pre-Processing module;
2) Collect more GRN methods to BuildNet module;
3) Update network fusion algorithms to NetEnrich module;
5) Test with integrated sc RNA-seq data.

# Cite
- :smile_cat: This work has been presented at [ISMB/ECCB 2019](https://www.iscb.org/ismbeccb2019);
- :paw_prints: Go to [my poster](https://f1000research.com/posters/8-1359);
- :page_with_curl: Reference: *Wu M, Kanev K, Roelli P and Zehn D. PySCNet:
A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data [version 1; not peer reviewed]. F1000Research 2019, 8(ISCB Comm J):1359 (poster) (doi: 10.7490/f1000research.1117280.1)*
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/d3c17aac77e14f6bb17b33f875ff7471)](https://app.codacy.com/manual/MingBit/PySCNet?utm_source=github.com&utm_medium=referral&utm_content=MingBit/PySCNet&utm_campaign=Badge_Grade_Dashboard)
[![License](https://img.shields.io/github/license/MingBit/PySCNet)](https://github.com/MingBit/PySCNet/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/MingBit/PySCNet.svg?branch=master)](https://travis-ci.org/MingBit/PySCNet)
[![Documentation Status](https://readthedocs.org/projects/pyscnet/badge/?version=latest)](https://pyscnet.readthedocs.io/en/latest/?badge=latest)
# PySCNet: A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data
There are four modules:  
1) **Pro-precessing**: initialize a gnetData object consisting of Expression Matrix, Cell Attributes, Gene Attributes and Network Attributes;  
2) **BuildNet**: reconstruct GRNs by various methods implemented in docker;  
3) **NetEnrich**: network analysis including consensus network detection, gene module identification and trigger path prediction as well as network fusion;  
4) **Visulization**: network illustration.  

### GnetData object contains the following parts:  
1) **ExpMatrix**: Raw Gene count matrix;  
2) **CellAttrs**: Cell annotation;  
3) **GeneAttrs**: Gene annotation;  
4) **NetAttrs**: multiple unstructured annotation (e.g. linkage table, graph, gene module, gene centrality).   


# Create your own GRNs
[Dashboard](https://github.com/MingBit/PySCNet/blob/master/images/ShinyApp.gif) is available now for creating your own GRNs.
Once the cells are grouped into several clusters and linkage tables are generated for each/all clusters, you can export the results
as pickle object and uplaod onto Shinyapp. Cell attributes, Gene attributes and Network attributes are illustrated here.
As shown belows, you can set your own thresholds to build each/all cluster-specific GRNs.

# Tutorial
PBMC data preprocessed and analyzed by scanpy as explained in this [tutorial](https://github.com/MingBit/PySCNet/blob/master/tutorial/pyscnet_pbmc.ipynb). 

# Cite
-  This work has been presented at [ISMB/ECCB 2019](https://www.iscb.org/ismbeccb2019);
-  Go to [my poster](https://f1000research.com/posters/8-1359);
-  Reference: *Wu M, Kanev K, Roelli P and Zehn D. PySCNet:
A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data [version 1; not peer reviewed]. F1000Research 2019, 8(ISCB Comm J):1359 (poster) (doi: 10.7490/f1000research.1117280.1)*
# Contact Author
Ming Wu  
Email: ming.wu@tum.de  
Github: https://github.com/MingBit  
