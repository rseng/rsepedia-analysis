# RNASSIST RNA Solutions: Synthesizing Information to Support Transcriptomics

This repo contains the source code of RNASSIST analysis modules developed under the Phase I of the NIAAA ASSIST project. These modules can be run in Jupyter notebooks, Docker containers, as well as workflows managed by the ADE (Active Discovery Engine) with increasing support for provenance.

### Description of the analysis modules:
|#| Analysis Module | Description |
|-|-------|-------------|
|1| Network Analysis | Network construction by WGCNA |
|2| Module Extraction | Network module detection by Louvain Algorithm |
|3| Module Membership Analysis | Genes count per module and module assignment stability check |
|4| Module DE/Diagnostic Correlation | Biological relevance of network module check |
|5| Module Network Embedding | Conversion of network to machine learning-friendly matrix representation |
|6| ML and Critical Gene Identifier | Machine learning to extract features for critical gene identification |

The analysis workflow follows the order of the modules and the modules are inter-related as depicted in the conceptual workflow below.
<p align="center"><img src="https://user-images.githubusercontent.com/12038408/123433360-60762480-d599-11eb-911e-1af52a23df4d.png" width="700" height="650">
</p>


Module ```Critical Gene Validation``` requires a 3rd party license and is thus not included in this repo.

### Where to get data:
Download data from [RNASSIST Dropbox](https://www.dropbox.com/sh/uajkuclelr409e3/AADUigHDIqBXCaaDvcHo9sv-a?dl=0) in the `data` folder. The user guide below assumes that you have downloaded all data files and placed them under the `data` subfolder of this project.

## User Guide

To run RNASSIST software on your machine, you need to have Java, Python and Docker installed. We recommend Java SE Runtime 15.0.2, Python 3.7+ (tested on 3.8.5) and Docker 20.10.6. GNU Make 3.81 is also required to build the Docker images. Below we describe how to set up and run the RNASSIST analysis modules in three different modes.

### 1. How to set up the environment for Jupyter notebooks
Jupyter notebooks for RNASSIST analysis modules are included to allow researchers to test out the analysis code using the Jupyter notebook interface. The `notebooks` folder contains requirements files capturing software dependencies for the three notebooks included. Corresponding requirement file is loaded into each notebook at the beginning of the notebook.

### 2. How to launch containers for each analysis module
Before analysis modules can be launched through standalone containers, the corresponding images need to be loaded. You can either use the included Makefile to generate the corresponding images, or download them from [RNASSIST Dropbox](https://www.dropbox.com/sh/uajkuclelr409e3/AADUigHDIqBXCaaDvcHo9sv-a?dl=0) in the `images/standalone` folder, place them under the `images/standalone` subfolder of this project, and load them using:
```
make load-standalone-images
```
The analysis modules are meant to be launched in sequence in the order listed in the above table and there are configuration files in the `config` folder specifying all input files needed to launch the module and where the module will be generating its output files and plots. Before choosing an analysis module to execute, make sure all the input data specified in the corresponding config file are available.

There is a script called launch.py under the scripts folder that can be used to launch these analysis modules, e.g., to launch `Module Extraction` on the human dataset, use: `python launch.py module_extraction human <path to the data folder>`, where `<path to the data folder>` is the absolute path to the `data` folder under the project root.

### 3. How to run RNASSIST modules in a workflow using ADE

#### Prepare ADE runtime environment
Download `ade_runtime.tgz` from [RNASSIST Dropbox](https://www.dropbox.com/sh/uajkuclelr409e3/AADUigHDIqBXCaaDvcHo9sv-a?dl=0) into the project root folder and unpack it using:
```
tar zxvf ade_runtime.tgz
```
This command will create the following folder structure under the project:
```
ade
├── bin
│   ├── launcher.bat
│   └── launcher.sh
├── create_node_docker_image.py
├── doc
│   ├── README.md
│   ├── action_props.gif
│   ├── connect.gif
│   ├── disconnect.gif
│   ├── doc_props.gif
│   ├── docker_props.png
│   ├── dynamic_props.gif
│   ├── export_data.gif
│   ├── launch.gif
│   ├── new_node.gif
│   ├── persist.gif
│   ├── remove_node.gif
│   ├── scroll_props.gif
│   ├── view_data.gif
│   ├── view_props.gif
│   └── workflow.png
└── repo
    ├── FastInfoset-1.2.16.jar
    ├── ST4-4.0.8.jar
    ├── ade-backend-1.0.0-SNAPSHOT.jar
    ├── ade-frontend-1.0.0-SNAPSHOT.jar
    ├── ade-launcher-1.0.0-SNAPSHOT.jar
...
```


#### Use ADE to run analysis workflow
Use the launch script (`launch.bat` or `launch.sh`) to start up the ADE workflow user interface. There are ready made workflows for both `human` and `mouse` datasets under the `workflows` folder of this repo that can be loaded into the user interface. For this, you need to first download the ADE images for the analysis modules from [RNASSIST Dropbox](https://www.dropbox.com/sh/uajkuclelr409e3/AADUigHDIqBXCaaDvcHo9sv-a?dl=0) in the `images/ade` folder and place them under `~/.ade_image_repo/netrias` on your machine. This is where ADE will be loading images into the runtime environment.

Follow [ADE documentation](./ade/doc/README.md) that provides detailed description on using the ADE user interface.


## Detailed data description of analysis modules

### For all the input/output data below, `Human` means it's for the human example data (Kapoor et al 2019). `Mouse` means it's for the mouse example data (Ferguson et al 2019).
The difference is because the two example datasets (Kapoor and HDID) we used had difference in the availability of the data. For example, `TOM co-expression network` and `gene module assignment by WGCNA hierarchical clustering` for the human data were provided to us but not available for the mouse data so the `Network Analysis` had to be run to construct these two files for the mouse. `subjects' alcohol metadata` was only available for the human data so all the analyses that involve diagnostics were skipped for the mouse data. 

**1. Network Analysis**

Note that for the Kapoor data used in our analysis (aka the human data), the `Network Analysis` module was skipped as the TOM network and the WGCNA module assignment were already published so the example for this module below is for the HDID mouse data. 

<table border="1">
    <thead>
      <tr>
        <th></th>
        <th><sub>File</sub></th>
        <th><sub>Description</sub></th>
        <th><sub>Human</sub></th>
        <th><sub>Mouse</sub></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><sub>Input</sub></td>
        <td><sub>PFC_HDID_norm_exp.txt</sub></td>
        <td><sub>normalized counts from RNA-seq or microarray</sub></td>
        <td><sub</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
    <tbody>
      <tr>
        <td rowspan=2><sub>Output</sub></td>
        <td><sub>tom.csv</sub></td>
        <td><sub>TOM co-expression network</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment by WGCNA hierarchical clustering</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
</table>

**2. Module Extraction**

<table border="1">
    <thead>
      <tr>
        <th></th>
        <th><sub>File</sub></th>
        <th><sub>Description</sub></th>
        <th><sub>Human</sub></th>
        <th><sub>Mouse</sub></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td rowspan=2><sub>Input</sub></td>
        <td><sub>Kapoor_TOM.csv</sub></td>
        <td><sub>TOM co-expression network</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>tom.csv</sub></td>
        <td><sub>TOM co-expression network</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
    <tbody>
      <tr>
        <td rowspan=2><sub>Output</sub></td>
        <td><sub>network_louvain_default.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using its default setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_agg1.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using a different setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
</table>

**3. Module Membership Analysis**

<table border="1">
    <thead>
      <tr>
        <th></th>
        <th><sub>File</sub></th>
        <th><sub>Description</sub></th>
        <th><sub>Human</sub></th>
        <th><sub>Mouse</sub></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td rowspan=4><sub>Input</sub></td>
        <td><sub>kapoor_wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment by WGCNA hierarchical clustering</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment by WGCNA hierarchical clustering</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_default.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using its default setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_agg1.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using a different setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
    <tbody>
      <tr>
        <td rowspan=3><sub>Output</sub></td>
        <td><sub>plot_gene_cnt_each_cluster_wgcna.png</sub></td>
        <td><sub>number of gene per module for WGCNA module assignment</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>plot_gene_cnt_each_cluster_louvain 1.png</sub></td>
        <td><sub>number of gene per module for Louvain module assignment # 1</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>plot_gene_cnt_each_cluster_louvain 2.png</sub></td>
        <td><sub>number of gene per module for Louvain module assignment # 2</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
</table>

**4. Module DE/Diagnostic Correlation**

<table border="1">
    <thead>
      <tr>
        <th></th>
        <th><sub>File</sub></th>
        <th><sub>Description</sub></th>
        <th><sub>Human</sub></th>
        <th><sub>Mouse</sub></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td rowspan=8><sub>Input</sub></td>
        <td><sub>deseq.alc.vs.control.age.rin.batch.gender.PMI. corrected.w.prot.coding.gene.name.xlsx</sub></td>
        <td><sub>differential expression analysis</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>de_data.csv</sub></td>
        <td><sub>differential expression analysis</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>kapoor_expression_Apr5.txt</sub></td>
        <td><sub>normalized counts from RNA-seq or microarray</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>kapoor_wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment by WGCNA hierarchical clustering</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment by WGCNA hierarchical clustering</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>Kapoor2019_coga.inia.detailed.pheno.04.12.17.csv</sub></td>
        <td><sub>subjects' alcohol metadata</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_default.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using its default setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_agg1.csv</sub></td>
        <td><sub>gene module assignment by Louvain algorithm using a different setting</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
    <tbody>
      <tr>
        <td rowspan=4><sub>Output</sub></td>
        <td><sub>expression_meta.csv</sub></td>
        <td><sub>normalized expression data joined with subjects' metadata</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>cluster_DE_perc_xx.png</sub></td>
        <td><sub>DEG distribution across modules</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>plot_sig_perc_xx.png</sub></td>
        <td><sub>% genes in the module that are significant for different alcohol trait group</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>cluster_phenotype_corr_xx.png</sub></td>
        <td><sub>module eigengene and alcohol trait correlation</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
    </tbody>
</table>

**5. Module Network Embedding**

<table border="1">
    <thead>
      <tr>
        <th></th>
        <th><sub>File</sub></th>
        <th><sub>Description</sub></th>
        <th><sub>Human</sub></th>
        <th><sub>Mouse</sub></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td rowspan=7><sub>Input</sub></td>
        <td><sub>Kapoor_TOM.csv</sub></td>
        <td><sub>TOM co-expression network</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>tom.csv</sub></td>
        <td><sub>TOM co-expression network</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>deseq.alc.vs.control.age.rin.batch.gender.PMI. corrected.w.prot.coding.gene.name.xlsx</sub></td>
        <td><sub>differential expression analysis</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>de_data.csv</sub></td>
        <td><sub>differential expression analysis</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>network_louvain_default.csv</sub></td>
        <td><sub>gene module assignment chosen to compare with embedding clusters (user's choice)</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>wgcna_modules.csv</sub></td>
        <td><sub>gene module assignment chosen to compare with embedding clusters (user's choice)</sub></td>
        <td><sub></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>expression_meta.csv</sub></td>
        <td><sub>normalized expression data joined with subjects' metadata</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
    </tbody>
    <tbody>
      <tr>
        <td rowspan=9><sub>Output</sub></td>
        <td><sub>embedding.csv</sub></td>
        <td><sub>network embedding</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>plot_gene_cnt_each_cluster_Network.png</sub></td>
        <td><sub>number of gene per network module (same as the output in <code>Module Membership Analysis</code></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>plot_gene_cnt_each_cluster_epoch=5_alpha=0.1.png</sub></td>
        <td><sub>number of gene per cluster for WGCNA module assignment (compare it with plot_gene_cnt_each_cluster_Network.png)</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>cluster_DE_perc_network.png</sub></td>
        <td><sub>DEG distribution across modules (same as the output in <code>Module DE/Diagnostic Correlation</code></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>cluster_DE_perc_epoch=100_alpha=0.1 embedding.png</sub></td>
        <td><sub>DEG distribution across embedding clusters (compare it with cluster_DE_perc_network.png)</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
      <tr>
        <td><sub>cluster_phenotype_corr_network.png</sub></td>
        <td><sub>module eigengene and alcohol trait correlation (same as the output in <code>Module DE/Diagnostic Correlation</code></sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>cluster_phenotype_corr_embedding.png</sub></td>
        <td><sub>cluster eigengene and alcohol trait correlation (compare it with cluster_phenotype_corr_network.png)</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>alcohol trait correlation network vs embedding.png</sub></td>
        <td><sub>distribution plot to compare cluster_phenotype_corr_network.png and cluster_phenotype_corr_embedding.png</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub></sub></td>
      </tr>
      <tr>
        <td><sub>cluster_jaccard_Network vs epoch=100_alpha=0.1.png</sub></td>
        <td><sub>pairwise jaccard comparison to determine network module and embedding cluster similarity</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
        <td><sub>:heavy_check_mark:</sub></td>
      </tr>
    </tbody>
</table>

**6. ML and Critical Gene Identifier**

<table border="1">
    <thead>
        <tr>
            <th></th>
            <th><sub>File</sub></th>
            <th><sub>Description</sub></th>
            <th><sub>Human</sub></th>
            <th><sub>Mouse</sub></th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=5><sub>Input</sub></td>
            <td><sub>Kapoor_TOM.csv</sub></td>
            <td><sub>TOM co-expression network</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub></sub></td>
        </tr>
        <tr>
            <td><sub>embedding.csv</sub></td>
            <td><sub>network embedding</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>expression_meta.csv</sub></td>
            <td><sub>normalized expression data joined with subjects' metadata</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub></sub></td>
        </tr>
        <tr>
            <td><sub>deseq.alc.vs.control.age.rin.batch.gender.PMI. corrected.w.prot.coding.gene.name.xlsx</sub></td>
            <td><sub>differential expression analysis</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub></sub></td>
        </tr>
        <tr>
            <td><sub>de_data.csv</sub></td>
            <td><sub>differential expression analysis</sub></td>
            <td><sub></sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
    </tbody>
    <tbody>
        <tr>
            <td rowspan=9><sub>Output</sub></td>
            <td><sub>critical_genes.csv</sub></td>
            <td><sub>candidate genes identified by RNASSIST</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>neighbor_genes.csv</sub></td>
            <td><sub>closest DEG neighbors in the co-expression network</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>run_ml_.png</sub></td>
            <td><sub>machine learning accuracy</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>run_ml_top_dims.png</sub></td>
            <td><sub>machine learning accuracy using only the most important dimensions</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>gene_phenotype_corr_for_xx.png</sub></td>
            <td><sub>critical gene/DEG/neighbor gene correlation with alcohol traits</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub></sub></td>
        </tr>
        <tr>
            <td><sub>alcohol trait correlation CG, neighbor & DEG.png</sub></td>
            <td><sub>distribution plot to compare gene_phenotype_corr_for_xx.png</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub></sub></td>
        </tr>
        <tr>
            <td><sub>jaccard_average_Important dim overlap within model repeats.png</sub></td>
            <td><sub>the important dimensions overlap between the repeats of each model</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>jaccard_critical_genes_Critical gene overlap between models.png</sub></td>
            <td><sub>critical gene overlap between each two models</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
        <tr>
            <td><sub>plot_nearby_impact_num_.png</sub></td>
            <td><sub>top 10 critical genes</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
            <td><sub>:heavy_check_mark:</sub></td>
        </tr>
    </tbody>
</table>

## Questions?
Please contact Yi-Pei Chen (ychen@netrias.com) or George Zheng (gzheng@netrias.com)
Last Updated: 10/5/2021

# ADE

The Active Discovery Engine (ADE) is an in-house workflow engine that provides a unified execution environment for various projects at Netrias. Projects are broken down into a set of components, where those components pass data to each other to form a workflow. Each component is isolated and configurable, meaning that re-use in future projects (or alternate versions of the same project) is straight-forward.

Why use ADE over other workflow engines? ADE provides several distinct advantages:

 * Highly robust data parsing infrastructure.
 * Streamlined container image creation and execution.
 * Assistive / predictive technology for building workflows.
 * Global data cache.
 * Execution logging.
 * Triggers, scheduling, and reporting.

The following subsections provide an overview of basic ADE usage via ADE's user interface. Be aware that ADE is **pre-alpha software**: While the backend is stable, the user interface has some bugs (none are blockers).

## Table of Contents

- [ADE](#ade)
  * [Table of Contents](#table-of-contents)
  * [Install Instructions](#install-instructions)
  * [Launch Instructions](#launch-instructions)
  * [Graph Manipulation](#graph-manipulation)
    + [Add Node](#add-node)
    + [Remove Node](#remove-node)
    + [Configure Node](#configure-node)
    + [Connect Nodes](#connect-nodes)
    + [Disconnect Nodes](#disconnect-nodes)
    + [View Data](#view-data)
    + [Export Data](#export-data)
    + [Load/Save](#load-save)
  * [Property Manipulation](#property-manipulation)
    + [Dynamic Properties](#dynamic-properties)
    + [Property Actions](#property-actions)
    + [Property Documentation](#property-documentation)
  * [Containerization](#containerization)
  * [Hotkey Reference](#hotkey-reference)

## Install Instructions

The following instructions are for OS X.

1. Ensure the following prerequisites are installed:
   * Java 13 (later versions may work as well)
   * Docker
1. Unpack the supplied ADE archive into a directory of your choosing.
1. Unpack the supplied image archive into `~/.ade_image_repo`.
1. Ensure `java` executable is accessible via `PATH` environment variable.

## Launch Instructions

![Launch ADE](launch.gif)

The following instructions are for OS X.

Prior to launching, ensure `java` executable is accessible via `PATH` environment variable.

Navigate to the ADE directory and run `./bin/launcher.sh` to launch. Once launched, a browser tab should automatically open and point to [http://localhost:8886](http://localhost:8886).

The browser tab becomes ready to use once a small block of text appears in the upper-right corner. Currently, a bug is present that may cause a delay of ~10s or so before the UI becomes usable.

## Graph Manipulation

![Example workflow](workflow.png)

Workflows are represented as directed acyclic graphs:

 * Nodes represent tasks that transform a set of tables.
 * Edges represent tables flowing between nodes.

Each node has a set of configuration properties that define how many inputs / outputs it has as well as how it transforms data. In the above example, the node `perform_network_embedding` (located near bottom-right) performs the network embedding task in an ASSIST workflow.

Node inputs are outputs are represented as connectors. Connectors span across the top and bottom of each node:

 * Top connectors are for input data (tables).
 * Bottom connectors are for output data (tables).

In the above example, the node `perform_network_embedding` (located near bottom-right) reads in 3 input tables and writes out 1 output table.

Colors represent the execution state of each node:

| Color      | Description                  |
|------------|------------------------------|
| 🟩 (green) | Node is actively executing.  |
| ⬜ (gray)  | Node executed successfully.  |
| 🟥 (red)   | Node executed with error.    |
| ⬛ (black) | Node is in an unknown state. |

In the above example, all nodes are gray, meaning that all nodes in the graph executed successfully.

Graph manipulation instructions are detailed in the sub-sections below. ADE automatically executes the graph as the user manipulates it. Caching is used to ensure that an execution only executes changed nodes and their children, not the entire graph.

### Add Node

![Add new node to graph](new_node.gif)

Press `N` to show the node selection textbox. Type in the node type and press `<ENTER>`.

The new node should appear at the mouse coordinates when `N` was originally pressed. Click-and-drag the node to change its position.

To cancel, press `<ESC>`.

### Remove Node

![Remove existing node from graph](remove_node.gif)

Hover over node and press `D`.

### Configure Node

![Configure node properties](view_props.gif)

Hover over node and press `V`. The tab that opens is horizontally split into 2 panes:

 * Top pane contains the node's configuration properties.
 * Bottom pane contains a preview of the node's input and output data.

As node configurations change, the graph re-executes. Cached executions are used where appropriate to speed up executions.

Invalid configuration inputs automatically revert.

Currently, a user interface bug causes configuration changes made in rapid succession to snap back.

### Connect Nodes

![Connect nodes](connect.gif)

Nodes have connectors (black dots) spanning across their top and bottom:

 * Top connectors are for input tables.
 * Bottom connectors are for output tables.

Hover over a connector and click-and-drag to another connector of opposing type (e.g. input connector to output connector). Connections forming cycles are rejected.

Note that hovering over a connector will pop-up a tooltip displaying the ID of that connector as well as the column names for the table at that connector.

### Disconnect Nodes

![Disconnect nodes](disconnect.gif)

Hover over connection and press `D`.

### View Data

![View node connector's data](view_data.gif)

Hover over connector and press `V`. The tab that opens is visually truncated to to 1000 rows by 25 columns (full data is preserved internally).

On graph re-execution, this tab automatically updates.

### Export Data

![Export node connector's data](export_data.gif)

Hover over connector and press `V`. In the tab that opens, press `E`. Type in the path of the CSV file to save and press `Enter`. 

Although the data displayed on this tab is visually truncated, the full data is exported to the CSV.

To cancel, press `<ESC>`.

### Load/Save

![Load or save workflow](persist.gif)

Press `P`. Type in the path of the YAML file to save and either click `Load` or `Save`. Any subsequent changes to the workflow will automatically persist to that workflow file.

To cancel, click `Cancel`.

To reset the workflow, click `New Graph`

## Property Manipulation

Nodes are configurable via their properties, accessible through ADE's property view ([access instructions](#configure-node)). 

Properties are displayed as a hierarchy tree, where each property in the tree is exposed as one of several UI elements:

| Type         | Description                                                          |
|--------------|----------------------------------------------------------------------|
| Text         | Single-line or multi-line string.                                    |
| Number       | 64-bit integer or float.                                             |
| Boolean      | True or false.                                                       |
| Selection    | Single choice out of several preset values.                         |
| Combinations | List, set, or map that groups together other properties as children. |

In addition, properties may be modified directly in a YAML, displayed at the very end of the tree.

![Scroll node properties](scroll_props.gif)

### Dynamic Properties

A property may have internal logic that causes any of the following behaviours on modification:

 1. A modified value may get reverted.
    * Example: Certain properties may not be negative (e.g. a calendar year can't be negative).
    * Example: Certain properties are length constrained (e.g. a CSV delimiter char must contain exactly 1 char).
    * Example: Certain properties must match some pattern (e.g. an e-mail can't be missing @).
 2. A modified value may get modified.
    * Example: Certain properties may have whitespace trimmed on input (e.g. whitespace at end of e-mail address).  
    * Example: Certain properties may require different encoding (e.g. URL unsafe characters).
 3. A modified value may cause the disappearance, introduction, and modification of other properties.
    * Example: Certain properties may require an entirely different set of surrounding properties based on their value (e.g. a linear regression model requires different tuning parameters than a deep learning model).

![Dynamic node properties](dynamic_props.gif)

### Property Actions

A property may have actions associated with it. These actions are exposed as buttons next to that property's UI element, and when clicked perform high-level tasks that modify that property and potentially other surrounding properties. For example, a property that represents which CSV delimiter to use may have an action that scans through an entire directory of CSVs to determine what the CSV delimiter should be.

![Node property actions](action_props.gif)

### Property Documentation

A property may have short-hand documentation associated with it. The documentation is visible as a tooltip when hovering over a property's label.

Be aware that only internal ADE nodes expose documentation. There is no documentation support for [containerized nodes](#containerization) at the moment.

![Node property documentation](doc_props.gif)

## Containerization

Container images that follow a set of ADE conventions may be wrapped as ADE nodes. Containerization provides the benefit of workflow reproducibility:

 * Reduction in workflow breakage across environments.
 * Increase in reproducible outputs across environments.

Containerized ADE nodes may expose special properties that control resources:

 * `docker_mem` controls the amount of memory (in bytes).
 * `docker_cpu` controls the amount of CPU (in cores).
 * `docker_volume` prefix controls volume mappings (host to container, `:` delimiter).

These special properties are not required to exist. If omitted, volume mappings are not possible / resource allocations are set directly by Docker's runtime.

![Container node properties](docker_props.png)

Streamlined creation of ADE-complaint images is possible through ADE's container creation tool. In most cases, this tool is simple enough that most non-technical users (e.g. scientists) are able to directly create images without supervision.

Note that ADE-complaint container images __are not__ locked-in to ADE. Each container image comes with a command-line interface that operates in a similar fashion as its node counterpart:

 * Input tables are passed in through command-line arguments.
 * Output tables are passed in through command-line arguments.
 * Properties are passed in through command-line arguments.
 * Resource allocations are set manually by the user on the Docker runtime (volume mappings, mem, cpu).

If not enough command-line arguments are supplied by the user, the container outputs the required signature vs the supplied arguments.

## Hotkey Reference

Graph view hotkeys:

**NOTE:** A hotkey may perform different functions depending on what graph element the mouse is hovering on.

| Key       | Context    | Description              |
|-----------|------------|--------------------------|
| `N`       | Graph      | Add node                 |
| `R`       | Node       | Replicate node           |
| `D`       | Node       | Delete node              |
| `D`       | Connection | Delete connection        |
| `.`       | Node       | Toggle inspection mode   |
| `Shift-N` | Node       | Rename node              |
| `V`       | Node       | Open properties          |
| `V`       | Connector  | View data                |
| `Q`       | Connector  | Calculate summary stats (press Q on two connectors and it will compare, Q twice on the same connector will summarize that table) |
| `Shift-Z` | Connector  | Quick line plot          |
| `Shift-X` | Connector  | Quick histogram          |
| `Shift-C` | Connector  | Quick heatmap            |
| `A`       | Node       | View artifacts           |
| `P`       | Graph      | Load or save             |
| `E`       | Graph      | Show environment vars    |
| `C`       | Graph      | Clear cache              |

Property view hotkeys:

| Key                           | Description                    |
|-------------------------------|--------------------------------|
| `ALT-0`                       | Scroll to top                  |
| `ALT-1` to `ALT-9`            | Scroll to corresponding output |
| `CTRL-ALT-1` to `CTRL-ALT+-9` | Scroll to corresponding input  |

Data view hotkeys:

| Key                           | Description                    |
|-------------------------------|--------------------------------|
| `E`                           | Export data as CSV             |
