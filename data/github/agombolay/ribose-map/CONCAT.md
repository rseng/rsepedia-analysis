![Logo](https://github.com/agombolay/Images/blob/master/Logo.png)
# A bioinformatics toolkit for mapping rNMPs in DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

**If you use Ribose-Map, please use the following citation**:  
Gombolay, AL, FO Vannberg, and F Storici. Ribose-Map: a bioinformatics toolkit to map ribonucleotides embedded in genomic DNA. *Nucleic Acids Research*, Volume 47, Issue 1, 10 Jan 2019, Page e5, https://doi.org/10.1093/nar/gky874. 

Ribose-Map is a user-friendly, standardized bioinformatics toolkit for the comprehensive analysis of ribonucleotide (rNMP) sequencing experiments. It can be applied to data generated using any high-throughput ribonucleotide sequencing technique, including ribose-seq, emRiboSeq, RHII-HydEn-seq, Alk-HydEn-seq, and Pu-seq.

## Output of each Module
1. **Alignment**: BAM file of read alignments and log file of alignment statistics
2. **Coordinate**: BED files of rNMP genomic coordinates and TAB files of rNMP counts at each coordinate
3. **Composition**: TXT files of raw counts and normalized percentages for r[A, C, G, U] and barcharts of results
4. **Sequence**: TAB files of raw and normalized frequencies for A, C, G, U/T at/surrounding rNMP sites and plots of results
5. **Distribution**: TAB and BedGraph files of per-nucleotide rNMP coverage normalized for read depth and plots of results
6. **Hotspot**: Sequence logo plots of consensus sequences for top X% most abundant sites of rNMP incorporation

## Hardware and Software Requirements
* Hardware (64-bit computer running Linux or Mac OS X; 4 GB of RAM)  
* Software dependencies and versions are listed in the ribose-map/lib/ribosemap_env.yaml file

## Required Files
1. FASTQ file of NGS rNMP-seq reads (SE or PE)
2. FASTA file of nucleotide sequence of reference genome
3. Chromosome sizes of reference genome (.chrom.sizes)
4. Bowtie2 index files of reference genome (.bt2)
5. [Configuration file](https://github.com/agombolay/ribose-map/blob/master/lib/sample.config)

**Note**: If you have a BED file of any single-nucleotide genomic coordinates (SNPs etc.) and want to input that file directly into the Sequence, Distribution, & Hotspot Modules, create a folder with the path below, save the file in this folder, and create a config file.

"$repository/results/$sample/coordinate$quality" ($variables should be those provided in config)

## Install Software, Timing ~ 30 min

1. **Create conda software environment**:  
* To install conda, please visit [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
* ribosemap_env.yaml is available [here](https://github.com/agombolay/ribose-map/blob/master/lib/ribosemap_env.yaml)
   ```bash
   conda env create --name ribosemap_env --file ribosemap_env.yaml
   ```

2. **Clone Ribose-Map GitHub repository**:  
* To install git, please visit [this link](https://git-scm.com/)
   ```bash
   git clone https://github.com/agombolay/ribose-map.git
   ```
   
## Run Ribose-Map, Timing < 3 hr for most datasets
Before proceeding, close current terminal window and open a new window to refresh the settings  
* If SAMtools gives an error, then install this dependency: conda install -c conda-forge ncurses

1. **Activate environment to access software**:
```bash
conda activate ribosemap_env
```

2. **Run scripts with configuration_file as input**:
```bash
ribose-map/modules/ribosemap {alignment, coordinate, composition, sequence, distribution, hotspot} config
```

3. **Once the analysis is complete, exit environment**:  
```bash
conda deactivate ribosemap_env
```
