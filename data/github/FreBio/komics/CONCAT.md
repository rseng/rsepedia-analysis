### Contents
  * [Introduction](#introduction)
  * [Citation](#citation)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tutorial](#tutorial)
  * [References](#references)


### Introduction
The rKOMICS software suite includes the python3.7 package KOMICS and the R package rKOMICS. These two software packages facilitate the assembly, circularization and downstream analyses of mitochondrial genomes in trypanosomatids. The input of KOMICS is reads in FASTQ format, and the output is maxicircle and circularized minicircles in FASTA format, which can be further processed by rKOMICS.

Please report any issues or questions to frederik.vandenbroeck AT kuleuven.be


### Reading
The following paper includes a **step-by-step tutorial** on how to best use the rKOMICS software suite: Geerts M, Schnaufer A, Van den Broeck F. [rKOMICS: an R package for processing mitochondrial minicircle assemblies in population-scale genome projects](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04384-1). BMC BIOINFORMATICS. 22:468.


### Citation
If you use **KOMICS** please cite: Van den Broeck F, Savill N, Imamura H, Sanders M, Maes I, Cooper S, Mateus D, Jara M, Adaui V, Arevalo J, Llanos-Cuentas A, Garcia L, Cupolillo E, Miles M, Berriman M, Schnaufer A, Cotton J, Dujardin JC. [Ecological divergence and hybridization of Neotropical *Leishmania* parasites](https://www.pnas.org/content/early/2020/09/18/1920136117). PROCEEDINGS OF NATIONAL ACADEMY OF SCIENCES. 117(40): 25159-25168.

If you use **rKOMICS** please cite: Geerts M, Schnaufer A, Van den Broeck F. [rKOMICS: an R package for processing mitochondrial minicircle assemblies in population-scale genome projects](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04384-1). BMC BIOINFORMATICS. 22:468.


### Installation
KOMICS has the following dependencies, which need to be installed first:
  * [MEGAHIT](http://www.metagenomics.wiki/tools/assembly/megahit)
  * [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [VSEARCH](https://github.com/torognes/vsearch)

Once the dependencies are installed, install the latest version of KOMICS using pip3:
```
pip install git+https://github.com/FreBio/komics.git
```

If you are running KOMICS on a supercomputer, you may want to have a look [here](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/python_package_management.html#alternatives-to-conda) on how to setup your local environment and pip install options.


### Usage
```
Usage: komics <command> [options] <required arguments>

To get minimal usage for a command use:
komics command

To get full help for a command use one of:
komics command -h
komics command --help

Available commands:
all:         	 Performs assemble, circularize and polish
assemble:    	 Assembles minicircles using high quality reads from trimfq
circularize: 	 Circularizes minicircles from assemble
polish:      	 Reorientate and filter circular minicircles
```


### Tutorial
If this is the first time you use KOMICS, you may want to follow the tutorial using the data provided on our github page. These sequence reads were generated using whole genome sequencing of *Leishmania peruviana* isolate LCA04. The following files contain the sequence reads that did not align to the *Leishmania braziliensis* M2904 reference genome (see step 1 below).
```
wget https://github.com/FreBio/komics/blob/master/data/LCA04_1.fq.gz
wget https://github.com/FreBio/komics/blob/master/data/LCA04_2.fq.gz
```


#### 1. Extract unaligned reads from BAM file
If you have a BAM file with sequence reads aligned against a nuclear reference genome, you first need to extract the unaligned reads (i.e those reads that likely originate from the mitochondrial genome) from the BAM file using [samtools](http://www.htslib.org), and then convert the BAM file into FASTQ files using [GATK](https://gatk.broadinstitute.org/hc/en-us):
```
samtools view -b -f 4 -o unmapped.reads.bam reads.bam
gatk SamToFastq -I unmapped.reads.bam -F reads1.fq.gz -F2 reads2.fq.gz
```


#### 2. Trim reads for high-quality
Once the sequence reads are in FASTQ format, it is recommended to trim sequences for high quality using e.g. [fastp](https://github.com/OpenGene/fastp):
```
fastp -i LCA04_1.fq.gz -I LCA04_1.fq.gz -o LCA04_trim_1.fq.gz -O LCA04_trim_2.fq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 125
```
You might need to change the setting -l (minimum read length) to a lower value if your reads are shorter, but we recommend to keep -b (maximum read length) to 125bp because Illumina sequencing quality decreases with increasing number of cycles (i.e. longer reads).


#### 3. Assemble the mitochondrial genome
Use `komics assemble` to assemble the mitochondrial maxicircle and minicircles using MEGAHIT:
```
komics assemble --threads 4 --kmin 29 --kmax 119 --kstep 20 LCA04_run1 LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
grep '>' LCA04_run1.maxicircles.fasta
grep -c '>' tmp.LCA04_run1.csb3contigs.fasta
```
The optimal kmer length depends on the complexity of the mitochondrial genome, and it will also be different for maxicircle and minicircles. For **minicircles**, we recommend using high kmer values that are close to the read length (e.g. kmer of 119 for reads that are 125 bp long). For **maxicircles**, we recommend using low kmer values (e.g. 29). In the example above, we follow a kmer sweep strategy, whereby MEGAhit is run using the following k-list: 29,49,69,89,109,119.

**Maxicircle** contigs are identified using BLAST against maxicircle sequences of Leishmania braziliensis, Trypanosoma brucei, T. equiperdum and T. lewisi. The resulting maxicircle contigs can be found in the file LCA04\_run1.maxicircles.fasta including one contig of length 17,202bp that covers the entire coding region. If you are interested in generating a complete circularized assembly of the maxicircle that includes both coding and divergent region, please read Van den Broeck et al. (2018).

**Minicircle** contigs are extracted based on the presence of the Conserved Sequence Block 3 (CSB3), a 12-bp minicircle motif, also known as the universal minicircle sequence, that is highly conserved across all Kinetoplastida species (Ray 1989). By default, KOMICS uses the known CSB-3 motif GGGGTTGGTGTA, GGGGTTGATGTA and their reverse complements to extract contigs of putative minicircle origin. The resulting minicircle contigs can be found in the file tmp.LCA04\_run1.csb3contigs.fasta including a total of 53 minicircle contigs. This seems like a rather low number of minicircles, and is due to the fact that we included low kmer values. Let's do another assembly, this time using only high kmer values:
```
komics assemble --threads 4 --kmin 99 --kmax 119 --kstep 10 LCA04_run2 LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
grep -c '>' tmp.LCA04_run2.csb3contigs.fasta
```
This yields better results as we have now assembled a total of 85 minicircle contigs.


#### 4. Circularize the mitochondrial minicircles
Once we are happy with the number of assembled minicircle contigs, we need to circularize the minicircle contigs. KOMICS uses BLAST as a strategy to identify a sequence that is in common at the start and the end of a given minicircle contig. MEGABLAST is run on the entire set of minicircle contigs with the low complexity filter turned off and allowing a maximum e-value of 10-5. The BLAST output is processed to retain only hits among the same minicircle contig (avoiding artificial dimers) with 100% identity and a minimum 20bp overlap at the start and end of a given contig. Whenever an overlap is found, the contig is classified as circular and the duplicated sequence at the start of the contig is removed.
```
komics circularize LCA04_run2 tmp.LCA04\_run2.csb3contigs.fasta
grep -c '>' tmp.LCA04_run2.circularized.fasta
grep -c '_circularized' tmp.LCA04_run2.circularized.fasta
```
Of the 85 minicircle contigs, 74 (87%) were successfully circularized.


#### 5. Polish the circularized minicircles
Finally, we want to align all minicircles by **(a)** reorienting each minicircle contig based on the CSB3-mer, **(b)** putting the Conserved Sequence Block 1 (CSB1) at the start of each circularized minicircle contig and **(c)** cluster contigs based on a minimum percent identity (e.g. 97%) using VSEARCH.
```
komics polish --minidentity 97 LCA04_run2 tmp.LCA04_run2.circularized.fasta
grep -c '>' LCA04_run2.minicircles.fasta 
```


#### 6. Check the quality and read depth of the assembly
To check the quality of the assembly, we will need to calculate some mapping statistics. One useful metric is the proportion of (near-)perfect alignments of CSB3-containing reads, which serves as a proxy for the total number of minicircles that were initially present within the DNA sample. Note that the CSB-3 12-mer is present within both the minicircles and maxicircles, so it's recommended to use only those CSB3-reads that did not align to the maxicircle.
First, let's create overlapping ends for the circularized minicircles by copying the last 150bp of the circularized contig to the start of each circularized contig. The reason we do this is because **(a)** we want to make sure to include all reads during the alignment process where we use high percent identities to include only those reads that almost perfectly match to the minicircle (see below) and **(b)** 2. read depth at overlapping ends of circular sequences should be approximately half the median read depth of the entire sequence.
```
wget https://github.com/FreBio/komics/raw/master/komics/fasta_extend.py
chmod 755 fasta_extend.py
./fasta_extend.py LCA04_run2.minicircles.extended.fasta LCA04_run2.minicircles.fasta 150
```

Now we can map the reads to the extended minicircles, e.g. using SMALT (or BWA):
```
smalt index -k 5 -s 2 LCA04_run2.minicircles.extended.fasta LCA04_run2.minicircles.extended.fasta
smalt map -f sam -y 0.95 -o run2.mapped.sam LCA04_run2.minicircles.extended.fasta LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
```

And run the following script to get some mapping and depth stats:
```
wget https://github.com/FreBio/komics/raw/master/komics/mapping_stats.sh
chmod 755 mapping_stats.sh
./mapping_stats.sh run2.mapped.sam
```

The following results table gives you a detailed overview of some mapping stats
```
cat run2.mapped.mapping.stats.sh
```

**(a)** Number of reads: 1923482
**(b)** Number of mapped reads: 1763455
**(c)** Number of reads with mapping quality >= 20: 1613847
**(d)** Number of CSB3-containing reads: 212826
**(e)** Number of mapped CSB3-containing reads: 198904
**(f)** Number of CSB3-containing reads with mapping quality >= 20: 163262

Results show that a total of 1,763,455 reads (92%) mapped to the minicircle, of which 84% reads with MQ > 20.
Of a total of 212,826 reads containing CSB3, 94% aligned against the minicircles and 77% with a relatively high mapping quality. This tells us that we have successfully assembled most of the minicircles that were present within the DNA sample.


The following results table gives you a detailed overview of the depth stats per minicircle contig
```
cat run2.mapped.depth.stats.sh
```
Minicircle copy numbers can then be estimated by standardizing median read depths per minicircle to the median genome-wide read depths. For instance, We know that the median genome-wide read depth for our sample is equal to 91. Assuming diploidy, this means that the haploid read depth is equal to ~45. Here, median minicircle read depths range between 218 and 7892. Hence, the minicircle copy numbers range between 2 and 87.


#### 7. Remove intermediate files
Once you are happy with the final set of maxicircles and minicircles, you can remove all intermediate files:
```
rm -r tmp.LCA04_run*
```


### References
Ray. 1989 [Conserved sequence blocks in kinetoplast minicircles from diverse species of trypanosomes](https://dx.doi.org/10.1128%2Fmcb.9.3.1365) Molecular and Cellular Biology.

Van den Broeck et al. 2018 [Mitonuclear Genomics Challenges the Theory of Clonality in Trypanosoma Congolense: Reply to Tibayrenc and Ayala](https://pubmed.ncbi.nlm.nih.gov/30142241/) Molecular Ecology.
