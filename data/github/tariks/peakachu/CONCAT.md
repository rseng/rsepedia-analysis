> **_NOTE:_**  Peakachu (version>=1.1.2) now supports both [.hic](https://github.com/aidenlab/juicer/wiki/Data) and [.cool](https://cooler.readthedocs.io/en/latest/datamodel.html) formats.

# Introduction
## What is Peakachu
Peakachu is an acronym that standands for Unveil Hi-C Anchors and Peaks. It takes genome-wide contact data as input and returns coordinates of likely interactions such as chromatin loops. A machine learning framework based on sklearn is employed to generate random forest models trained on example interactions predicted by an arbitrary experiment. For example, loops can be predicted in Hi-C data using models trained with the Hi-C profiles of interactions detected via ChIA-PET. Although Hi-C is the namesake of Peakachu, it is designed to accept any genome-wide contact map including those from Micro-C and DNA SPRITE.

## Citation
Salameh, T.J., Wang, X., Song, F. et al. A supervised learning framework for chromatin loop detection in genome-wide contact maps. Nat Commun 11, 3428 (2020). https://doi.org/10.1038/s41467-020-17239-9

Peakachu requires Python3 and several scientific packages to run. It is best to set up a conda environment then install from github. Copy and paste the command snippets below:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n 3dgenome python=3.6 scikit-learn=0.20.2 numpy scipy pandas h5py cooler
source activate 3dgenome
pip install hic-straw
git clone https://github.com/tariks/peakachu
cd peakachu
python setup.py install
```

Peakachu should now be installed as a command-line tool within the new environment. Options for all peakachu commands and sub-commands can be accessed with the -h option. 


```bash
peakachu -h
```

    usage: peakachu [-h] {train,score_chromosome,score_genome,depth,pool} ...
    
    Train Random Forest with Hi-C data and training peaklist.
    
    positional arguments:
      {train,score_chromosome,score_genome,depth,pool}
        train               Train RandomForest model per chromosome
        score_chromosome    Calculate interaction probability per pixel for a
                            chromosome
        score_genome        Calculate interaction probability per pixel for the
                            whole genome
        depth               Print total intra-chromosomal read count
        pool                Print centroid loci from score_x output
    
    optional arguments:
      -h, --help            show this help message and exit


# Example: predicting loops in GM12878 Hi-C

GM12878 is a commonly studied cell-line based on lymphoblasts from an adult individual. The following example will download a cooler file from a public source, train a series of models using ChIA-PET or HiChIP data, then predict loops using the trained models.

## Data preparation

Peakachu requires the contact map to be a .cool file or a .hic file and any training input to be a text file in bedpe format. Example training data can be found at the [training-sets](https://github.com/tariks/peakachu/tree/master/training-sets) subfolder. Cooler files may be found at the [4DN data portal](https://data.4dnucleome.org/).


```bash
wget ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-GM12878-MboI-allreps-filtered.10kb.cool
```

## Train a model and predict loops
It is always a good idea to call the help function immediately before entering a command:


```bash
peakachu train -h
```

    usage: peakachu train [-h] [-r RESOLUTION] [-p PATH] [--balance] [-O OUTPUT]
                          [-w WIDTH] [-b BEDPE]
    
    optional arguments:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp, default 10000
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -O OUTPUT, --output OUTPUT
                            Folder path to store results.
      -w WIDTH, --width WIDTH
                            Number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -b BEDPE, --bedpe BEDPE
                            Path to the bedpe file containing positive training
                            set.



```bash
peakachu train -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O models -b gm12878.mumbach.h3k27ac-hichip.hg19.bedpe
```

This will train one 23 random forest models, each labeled by a chromosome. Every model was trained on all of the interactions from the bedpe files. The purpose of this is to avoid Peakachu to predict loops from the same map it used for training, without overfitting. To use these models, you may either use the score_chromosome function to predict loops in only one chromosome, or the score_genome function when using a trained model to predict loops in a new contact map.


```bash
peakachu score_chromosome -h
```

    usage: peakachu score_chromosome [-h] [-r RESOLUTION] [-p PATH] [--balance]
                                     [-O OUTPUT] [-w WIDTH] [-m MODEL] [-l LOWER]
                                     [-u UPPER]
    
    optional arguments:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp, default 10000
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -O OUTPUT, --output OUTPUT
                            Folder path to store results.
      -w WIDTH, --width WIDTH
                            Number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -m MODEL, --model MODEL
                            Path to pickled model file.
      -l LOWER, --lower LOWER
                            Lower bound of distance between loci in bins (default
                            2).
      -u UPPER, --upper UPPER
                            Upper bound of distance between loci in bins (default
                            300).



```bash
for i in models/*pkl; do peakachu score_chromosome -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O scores -m $i; done
for i in scores/*; do peakachu pool -i $i -t .9 > ${i}.loops.txt; done
```

The pool function serves to select the most significant non-redundant results from per-pixel probabilities calculated by the score functions. It is recommended to try different probability thresholds to achieve the best sensitivity-specificity tradeoff. The output is a standard bedpe file with the 7th and final column containing the predicted probability from the sklearn model, to support further filtering. The results can be visualized in [juicebox](https://github.com/aidenlab/Juicebox) or [higlass](https://docs.higlass.io) by loading as 2D annotations. Here is an example screenshot of predicted GM12878 loops in juicer:
![Predicted loops from model trained on H3K27ac HiChIP interactions](https://github.com/tariks/peakachu/blob/master/example/gm12878-h3k27ac-loops.png)

# Using Peakachu as a standard loop caller

Models for predicting loops in Hi-C have already been prepared using CTCF ChIA-PET and H3K27ac HiChIP as training sets, at a variety of read depths. Simply download the appropriate model file and run the score_genome function.

|   Total intra reads  |  CTCF Model Link | H3K27ac Model  Link |
|----------------------|-------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
|      2 billion       | [CTCF total](https://dl.dropboxusercontent.com/s/enyg2m7ebj8mxsv/down100.ctcf.pkl?dl=0)   | [H3K27ac total](https://dl.dropboxusercontent.com/s/yasl5hu0v510k2v/down100.h3k27ac.pkl?dl=0)   |
|     1.8 billion      | [CTCF 90%](https://dl.dropboxusercontent.com/s/g12hy9f28igh0ng/down90.ctcf.pkl?dl=0)      | [H3K27ac 90%](https://dl.dropboxusercontent.com/s/kdbv52eeilkzqfr/down90.h3k27ac.pkl?dl=0)      |
|     1.6 billion      | [CTCF 80%](https://dl.dropboxusercontent.com/s/n2m4jxxojh0u5ay/down80.ctcf.pkl?dl=0)      | [H3K27ac 80%](https://dl.dropboxusercontent.com/s/45ekayzigeyuown/down80.h3k27ac.pkl?dl=0)      |
|     1.4 billion      | [CTCF 70%](https://dl.dropboxusercontent.com/s/h9vm8z0uysti8xm/down70.ctcf.pkl?dl=0)      | [H3K27ac 70%](https://dl.dropboxusercontent.com/s/mrhe0uayv402vfk/down70.h3k27ac.pkl?dl=0)      |
|     1.2 billion      | [CTCF 60%](https://dl.dropboxusercontent.com/s/cfkfem4w8dhhgwm/down60.ctcf.pkl?dl=0)      | [H3K27ac 60%](https://dl.dropboxusercontent.com/s/0f9xv6ljjlcwnsv/down60.h3k27ac.pkl?dl=0)      |
|      1 billion       | [CTCF 50%](https://dl.dropboxusercontent.com/s/c0b6axxb16p2nd7/down50.ctcf.pkl?dl=0)      | [H3K27ac 50%](https://dl.dropboxusercontent.com/s/3w4befpvu7c7cqe/down50.h3k27ac.pkl?dl=0)      |
|     800 million      | [CTCF 40%](https://dl.dropboxusercontent.com/s/8lvcdjenyoc8ggy/down40.ctcf.pkl?dl=0)      | [H3K27ac 40%](https://dl.dropboxusercontent.com/s/xwlk864nkoafzsy/down40.h3k27ac.pkl?dl=0)      |
|     600 million      | [CTCF 30%](https://dl.dropboxusercontent.com/s/f1383jpzj3addi4/down30.ctcf.pkl?dl=0)      | [H3K27ac 30%](https://dl.dropboxusercontent.com/s/dyvtyqvu3wpq3a5/down30.h3k27ac.pkl?dl=0)      |
|     400 million      | [CTCF 20%](https://dl.dropboxusercontent.com/s/a5nwa1xlg22ud24/down20.ctcf.pkl?dl=0)      | [H3K27ac 20%](https://dl.dropboxusercontent.com/s/qjm84cpw3uzlidp/down20.h3k27ac.pkl?dl=0)      |
|     200 million      | [CTCF 10%](https://dl.dropboxusercontent.com/s/cqi0ws8een9ad4t/down10.ctcf.pkl?dl=0)      | [H3K27ac 10%](https://dl.dropboxusercontent.com/s/q8mlwn4mz6rnumr/down10.h3k27ac.pkl?dl=0)      |
|      30 million      | [CTCF 1.5%](https://dl.dropboxusercontent.com/s/5gxeervadlga1b3/down1.ctcf.pkl?dl=0)      | [H3K27ac 1.5%](https://dl.dropboxusercontent.com/s/uh98lt1rbyauhgn/down1.h3k27ac.pkl?dl=0)      |

To make it clear, let's download another Hi-C dataset from 4DN: https://data.4dnucleome.org/files-processed/4DNFI5IHU27G/@@download/4DNFI5IHU27G.mcool. Peakachu provides a handy function `peakachu depth` to extract the total number of intra-chromosomal pairs from *[cool URI](https://cooler.readthedocs.io/en/latest/concepts.html#uri-string)*:


```bash
peakachu depth -p 4DNFI5IHU27G.mcool:resolutions/10000
```

 592991890

According to the table, for ~600 million intra-reads, we recommend using the 30% models in your prediction. Please refer to our paper to learn the differences between the CTCF and H3K27ac models.


```bash
peakachu score_genome -r 10000 --balance -p 4DNFI5IHU27G.mcool:resolutions/10000 -O scores -m down30.ctcf.pkl
for i in scores/*; do peakachu pool -i $i -t .9 > ${i}.loops.txt; done
```
# Not just Hi-C
Peakachu has been tested on Hi-C, Micrco-C, and DNA SPRITE contact maps with good results. For training sets, ChIA-PET, HiChIP, and PLAC-Seq have been tested. The purpose of this software is ultimately to facilitate the interpretation of results from multiple types of experiments, and the user is encouraged to apply Peakachu's training framework to newer approaches as they become available.
