[![DOI](https://joss.theoj.org/papers/10.21105/joss.02743/status.svg)](https://doi.org/10.21105/joss.02743)

## About REGENS :dna:

REGENS (REcombinatory Genome ENumeration of Subpopulations) is an open source Python package :package: that simulates whole genomes from real genomic segments. 
REGENS recombines these segments in a way that simulates completely new individuals while simultaneously preserving the input genomes' linkage disequilibrium (LD) pattern with extremely high fidelity. REGENS can also simulate mono-allelic and epistatic single nucleotide variant (SNV) effects on a continuous or binary phenotype without perturbing the simulated LD pattern. REGENS was measured to be 88.5 times faster and require 6.2 times lower peak RAM on average than a similar algorithm called Triadsim.

## :star2: IMPORTANT NOTICE (PLEASE READ) :star2:

REGENS's simulated genomes are comprised entirely of concatenated segments from the input dataset's real genomes. If your input genomes are not available for public use, then you may not be allowed to publicly release the simulated dataset. Please consult the institutions that provide you access to your input genotype dataset for more information about this matter.

## Instructions to Installing REGENS :hammer_and_wrench:

1. [Install conda](https://docs.conda.io/en/latest/miniconda.html) if you haven't already installed either Anaconda or Miniconda
2. Open your conda terminal. Type "Anaconda" or "Miniconda" into your search bar and open the terminal. It will look like this: <img src="images/conda_terminal.png" width="800" height="400"/>
3. Click the "Anaconda Prompt" app (left) to open the black terminal (right). The terminal's top must say "Anaconda prompt"
4. Enter ```conda create --name regens python=3.7``` in the terminal to create a new environment called regens with python version 3.7
5. Enter ```conda activate regens``` in the terminal to enter your new environment. If that doesn't work, enter ```source activate regens```
6. Once in your regens environment (repeat step 5 if you close and reopen the conda terminal), enter ```pip install regens==0.2.0```
7. Run [this command](https://github.com/EpistasisLab/regens/blob/main/README.md#simulate-genotype-data-computer) to allow regens to download the remaining files. It will write the simulated data into the `examples` folder that it downloads. If you experience permissions issues with this step, [try these remedies](https://github.com/EpistasisLab/regens/blob/main/README.md#remedies-to-known-permission-issues-adhesive_bandage):

## Input :turkey:
REGENS requires the following inputs:
- YOU MUST PROVIDE: real genotype data formatted as a standard (bed, bim, fam) plink _fileset_, ideally containing a minimum of 80 unrelated individuals.
- WE HAVE PROVIDED: a recombination map for every 1000 genomes population. 

The provided recombination maps were created by the [pyrho algorithm](https://github.com/popgenmethods/pyrho) and modified by us to minimize required disk space. [Recombination maps between related populations are highly correlated](https://github.com/EpistasisLab/regens-analysis/blob/master/README.md#about-the-recombination-maps-input-that-we-provided-turkey), so you could pair your input dataset with the recombination map of the most genetically similar 1000 genomes population, as is usually done for SNP imputation. If you wish to make your own recombination map, then it must be formatted as described [here](https://github.com/EpistasisLab/regens/blob/main/README.md#simulate-genotype-data-with-custom-recombination-rate-dataframes-abacus).  

## Output :poultry_leg:
REGENS outputs a standard (bed, bim, fam) plink fileset with the simulated genotype data (and optional phenotype information). 
If plink is not available to you, please consider [bed-reader](https://pypi.org/project/bed-reader/0.1.1/), which reads (bed, bim, fam) plink filesets into the python environment quickly and efficiently. 
 
In phenotype simulation, REGENS also outputs a file containing the R<sup>2</sup> value of the phenotype/genotype correlation and the *inferred* beta coefficients (see [example](https://github.com/EpistasisLab/regens/blob/main/correctness_testing_ACB/ACB_simulated_model_profile.txt)), which will most likely be close to but not equal to the input beta coefficients.

## Simulate genotype data :computer:

The following command uses `ACB.bed`, `ACB.bim`, and `ACB.fam` to simulate 10000 individuals without phenotypes. This command (or any other) will also complete the [final installation step](https://github.com/EpistasisLab/regens/blob/main/README.md#instructions-to-installing-regens-hammer_and_wrench). Windows users should replace all `\`  linebreak characters with `^`.

```shell
python -m regens \
  --in input_files/ACB \
  --out ACB_simulated \
  --simulate_nsamples 10000 \
  --simulate_nbreakpoints 4 \
  --population_code ACB \
  --human_genome_version hg19
```

## Simulate genotype data with custom recombination rate dataframes :abacus:

The following command uses custom recombination rate files as [input for regens](https://github.com/EpistasisLab/regens/blob/main/README.md#input-turkey) instead of the ones provided in the `hg19` and `hg38` folders (though the content in `input_files/hg19_ACB_renamed_as_custom` is just a copy of the content in `hg19/ACB`).  

```shell
python -m regens \
  --in input_files/ACB \
  --out ACB_simulated \
  --simulate_nsamples 10000 \
  --simulate_nbreakpoints 4 \
  --recombination_file_path_prefix input_files/hg19_ACB_renamed_as_custom/custom_chr_
```

Custom recombination rate files are to be named and organized as follows:
- The recombination map must be a single folder (named `hg19_ACB_renamed_as_custom` in the example above) with one gzipped tab seperated dataframe per chromosome.
- Every gzipped tab seperated dataframe must be named as `prefix_chr_1.txt.gz`, then `prefix_chr_2.txt.gz`  all the way through `prefix_chr_22.txt.gz`(`prefix` is named `custom_chr` in the example above).
- the `.txt.gz` files must actually be gzipped (as opposed to a renamed `txt.gz` extension). 
- Each chromosome's recombination map file must contain two tab separated columns named `Position(bp)` and	`Map(cM)`.

The `Position(bp)` column in each chromosome's recombination map is to be formatted as follows:
- The i<sup>th</sup> row of "Position(bp)" contains the genomic position of the left boundary for the i<sup>th</sup> genomic interval.
- The i<sup>th</sup> row of "Position(bp)" is also the genomic position of the right boundary for the (i-1)<sup>th</sup> genomic interval. 
- As such, the last row of "Position(bp)" is only a right boundary, and the first row is only a left boundary. 
- Genomic positions must increase monotonically from top to bottom. 

The `Map(cM)` column in each chromosome's recombination map is to be formatted as follows:
- The i<sup>th</sup> value of "Map(cM)" is the cumulative recombination rate from the first position to the i<sup>th</sup> position in CentiMorgans. 
- In other words, the recombination rate of the interval in between any two rows b and a must equal the Map(cM) value at row b minus the Map(cM) value at row a. 
- As such, the cumulative Map(cM) values must increase monotonically from top to bottom.
- The value of `Map(cM)` in the first row must be 0. 

An example of how this must be formatted is below (remember that there must be one per chromosome, and they must all be gzipped):

```shell
Position(bp)	Map(cM)
16050114	0.0
16058757	0.01366
16071986	0.03912
16072580	0.04013
16073197	0.04079
```

## :apple: Simulate genotype data with phenotype associations :green_apple:

Given at least one set of one or more SNPs, REGENS can simulate a correlation between each set of SNPs and a binary or continuous phenotype.
Different genotype encodings can be applied:

- Normally, if A is the major allele and a is the minor allele, then (AA = 0, Aa = 1, and aa = 2). However, you can _Swap_ the genotype values so that (AA = 2, Aa = 1, and aa = 0).
- You can further transform the values so that they reflect no effect (I), a dominance effect (D), a recessive effect (R), a heterozygous only effect (He), or a homozygous only effect (Ho).

The table below shows how each combination of one step 1 function (columns) and one step 2 function (rows) transforms the original (AA = 0, Aa = 1, and aa = 2) values.

|  Input = {0, 1, 2}   |Identity (I) |   Swap    |
|----------------------|-------------|-----------|
|Identity (I)          |  {0, 1, 2}  | {2, 1, 0} |
|Dominance (D)         |  {0, 2, 2}  | {2, 2, 0} |
|Recessive (R)         |  {0, 0, 2}  | {2, 0, 0} |
|Heterozygous only (He)|  {0, 2, 0}  | {0, 2, 0} |
|Homozygous only (Ho)  |  {2, 0, 2}  | {2, 0, 2} |

### Example 1: a simple additive model :arrow_lower_right:

A full command for REGENS to simulate genomic data with correlated phenotypes would be formatted as follows:

```shell
python -m regens \
  --in input_files/ACB --out ACB_simulated \
  --simulate_nsamples 10000 --simulate_nbreakpoints 4 \
  --phenotype continuous --mean_phenotype 5.75 \
  --population_code ACB --human_genome_version hg19 \
  --causal_SNP_IDs_path input_files/causal_SNP_IDs.txt \
  --noise 0.5 --betas_path input_files/betas.txt
```

This command simulates genotype-phenotype correlations according to the following model.
If we let _y_ be an individual's phenotype, s<sub>i</sub> be the i<sup>th</sup> genotype to influence the value of _y_ such that (AA = 0, Aa = 1, and aa = 2), and _B_ be the bias term. The goal is to simulate the following relationship between genotypes and phenotype:

y = 0.5s<sub>1</sub> + 0.5s<sub>2</sub> + 0.5s<sub>3</sub> + B + &epsilon;

where &epsilon; ~ N(&mu; = 0, &sigma;<sub>&epsilon;</sub> = 0.5E[y]) and E[y] = 5.75.

Notice that the values of &sigma;<sub>&epsilon;</sub> and E[y] are determined by the `--noise` and `--mean_phenotype` arguments.

<!-- h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x -->
<!-- <img src="https://render.githubusercontent.com/render/math?math=y = 0.2s_1 %2B 0.2s_2 %2B 0.2s_3 %2B B %2B \epsilon"> -->
<!-- <img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])"> -->
<!-- <img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75"> -->

The following files, formatted as displayed below, must exist in your working directory.
`input_files/causal_SNP_IDs.txt` contains newline seperated SNP IDs from the input bim file `input_files/ACB.bim`:
```
rs113633859
rs6757623
rs5836360
```
`input_files/betas.txt` contains one (real numbered) beta coefficient for each row in `input_files/causal_SNP_IDs.txt`:
```
0.5
0.5
0.5
```

### Example 2: inclusion of nonlinear single-SNP effects :arrow_heading_down:

```shell
python -m regens \
  --in input_files/ACB --out ACB_simulated \
  --simulate_nbreakpoints 4 --simulate_nsamples 10000 \
  --phenotype continuous --mean_phenotype 5.75 \
  --population_code ACB --human_genome_version hg19 --noise 0.5 \
  --causal_SNP_IDs_path input_files/causal_SNP_IDs.txt \
  --major_minor_assignments_path input_files/major_minor_assignments.txt \
  --SNP_phenotype_map_path input_files/SNP_phenotype_map.txt \
  --betas_path input_files/betas.txt
```

In addition to the notation from the first example, let S<sub>i</sub> = _swap_(s<sub>i</sub>) be the i<sup>th</sup> genotype to influence the value of _y_ such that (AA = 2, Aa = 1, and aa = 0). Also, we recall the definitions for the four nontrivial mapping functions (R, D, He, Ho) defined prior to the first example. The second example models phenotypes as follows:

y = 0.5s<sub>1</sub>+ 0.5R(S<sub>2</sub>) + 0.5He(s<sub>3</sub>) + B + &epsilon;

where &epsilon; ~ N(&mu; = 0, &sigma;<sub>&epsilon;</sub> = 0.5E[y]) and E[y] = 5.75.

<!-- <img src="https://render.githubusercontent.com/render/math?math=y = 0.2R(s_2) %2B 0.2D(s_3) %2B 0.2S_6 %2B B %2B \epsilon"> -->
<!-- <img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])">
<img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75"> -->

Specifying that (AA = 2, Aa = 1, and aa = 0) for one or more alleles is optional and requires `input_files/major_minor_assignments.txt`.

```
0
1
0
```

Specifying the second genotype's recessiveness (AA = 0, Aa = 0, and aa = 2) and third genotype's heterozygosity only (AA = 0, Aa = 2, and aa = 0) is optional and requires `input_files/SNP_phenotype_map.txt`. 

```
regular
recessive
heterozygous_only
```

### Example 3: inclusion of epistatic effects :twisted_rightwards_arrows:

REGENS models epistasis between an arbitrary number of SNPs as the product of transformed genotype values in an individual.

```shell
python -m regens \
  --in input_files/ACB --out ACB_simulated \
  --simulate_nbreakpoints 4 --simulate_nsamples 10000 \
  --phenotype continuous --mean_phenotype 5.75 \
  --population_code ACB --human_genome_version hg19 --noise 0.5 \
  --causal_SNP_IDs_path input_files/causal_SNP_IDs2.txt \
  --major_minor_assignments_path input_files/major_minor_assignments2.txt \
  --SNP_phenotype_map_path input_files/SNP_phenotype_map2.txt \
  --betas_path input_files/betas.txt
```

y = 0.5s<sub>1</sub> + 0.5D(s<sub>2</sub>)s<sub>3</sub>+ 0.2Ho(S<sub>4</sub>)s<sub>5</sub>s<sub>5</sub> + B + &epsilon;

where &epsilon; ~ N(&mu; = 0, &sigma;<sub>&epsilon;</sub> = 0.5E[y]) and E[y] = 5.75.

<!-- <img src="https://render.githubusercontent.com/render/math?math=y = 0.2s_1 %2B 0.2s_2R(s_3) %2B 0.2Ho(S_4)s_5s_5 %2B B %2B \epsilon"> -->
<!-- <img src="https://render.githubusercontent.com/render/math?math=\epsilon ~ N(\mu = 0, \sigma_{\epsilon} = 0.5E[y])">
<img src="https://render.githubusercontent.com/render/math?math=E[y] = 5.75"> -->

Specifying that multiple SNPs interact (or that rs62240045 has a polynomic effect) requires placing all participating SNPs in the same tab seperated line of `input_files/causal_SNP_IDs.txt`

```
rs11852537
rs1867634	rs545673871
rs2066224	rs62240045	rs62240045
```

For each row containing one or more SNP IDs, `input_files/betas.txt` contains a corresponding beta coefficient. (Giving each SNP that participates in a multiplicative interaction its own beta coefficient would be pointless).

```
0.5
0.5
0.5
```

As before, both `input_files/major_minor_assignments.txt` and `input_files/SNP_phenotype_map.txt` are both optional.
If they are not specified, then all genotypes of individual SNPs will have the standard (AA = 0, Aa = 1, and aa = 2) encoding.

For each SNP ID, `input_files/major_minor_assignments.txt` specifies whether or not untransformed genotypes of individual SNPs follow the (AA = 2, Aa = 1, and aa = 0) encoding.
CAUTION: In general, if a SNP appears in two different effects, then it may safely have different major/minor assignments in different effects.
However, if a SNP appears twice in the same effect, then make sure it has the same major/minor assignment within that effect, or that effect may equal 0 depending on the map functions that are used on the SNP. 

```
0	
0	0
1	0	0
```
For each SNP ID, `input_files/SNP_phenotype_map.txt` specifies whether the SNP's effect is regular, recessive, dominant, heterozygous_only, or homozygous_only.

```
regular
dominant	regular
homozygous_only	regular	regular
```

In context to an epistatic interaction, first the _Swap_ function is applied to SNPs specified by `input_files/major_minor_assignments.txt`, then map functions specified by `input_files/SNP_phenotype_map.txt` are applied to their respective SNPs. The transformed genotypes of SNPs in the same row of `input_files/causal_SNP_IDs.txt` are multiplied together and by the corresponding beta coefficient in `input_files/betas.txt`. Each individual's phenotype is then the sum of row products, the bias, and the random noise term.

## REGENS simulates nearly flawless GWAS data :100:

The Triadsim algorithm simulates LD patterns that are almost indistinguishable from those of the input Dataset. REGENS uses Triadsim's method of recombining genomic segments to simulate equally realistic data, and we measured it to be 88.5 times faster (95%CI: [75.1, 105.0]) and require 6.2 times lower peak RAM (95%CI [6.04, 6.33]) than Triadsim on average. The following three figures show that REGENS nearly perfectly replicates the input dataset's LD pattern. 

1. For the 1000 genome project's ACB population, this figure compares (right) every SNP's real maf against it's simulated maf and (left) every SNP pair's real genotype pearson correlation coefficient against its simulated genotype pearson correlation coefficient for SNP pairs less than  200 kilobases apart. 100000 simulated individuals were used (please note that the analogous figures in the correctness testing folders were only made with 10000 samples, so they have slightly more noticable random noise than the figure displayed directly below). 

![Real and simulated R value vs. MAF](images/r_maf_ACB.png)

2. For the 1000 genome project's ACB population, this figure plots SNP pairs' absolute r values against their distance apart (up to 200 kilobases apart) for both real and simulated populations. More specifically, SNP pairs were sorted by their distance apart and seperated into 4000 adjacent bins, so each datapoint plots one bin's average absolute r value against its average position. 100000 simulated individuals were used. 

![Real and simulated R value vs. distance_profile](images/r_dist_ACB.png)

3. These figures compare TSNE plots of the first 10 principal components for real and simulated 1000 genomes subpopulations that are members of the AFR superpopulation. Principal components were computed from all twenty-six 1000 genomes population datasets, and the loadings were used to project the simulated individuals onto the PC space. These results demonstrate that REGENS replicates the the input data's overall population structure in simulated datasets. 10000 simulated individuals were used. 

![TSNE1 vs TSNE2 for 1000 genome African subpopulations](images/tsne.png)

## supplementary details :european_castle:

Thank you concerned reader (for making it this far)!

But our full analysis is in [another repository](https://github.com/EpistasisLab/regens-analysis)!

## Repository structure

### Folders for REGENS to download :inbox_tray:

  * `correctness_testing_ACB`: A directory containing bash scripts to test code correctness on the ACB subpopulation, as well as the output for those tests. Correctness testing part 2 is optional and requires plink version 1.90Beta.
  * `correctness_testing_GBR`: A directory containing bash scripts to test code correctness on the GBR subpopulation, as well as the output for those tests. Correctness testing part 2 is optional and requires plink version 1.90Beta.
  * `examples`: A directory containing bash scripts that run the data simulation examples in the README.
  * `hg19`: for each 1000 genomes project population, contains a folder with one gzipped recombination rate dataframe per hg19 reference human autosome.
  * `hg38`: for each 1000 genomes project population, contains a folder with one gzipped recombination rate dataframe per hg38 reference human autosome.
  * `input_files`: contains examples of regens input that is meant to be provided by the user. The example custom recombination rate information is copied from that of the hg19 mapped ACB population. Also contains input for the Triadsim algorithm. The genetic input labeled as "not_trio" for Triadsim is comprised of ACB population duplicates and is only meant to compare Triadsim's runtime. 
  * `runtime_testing_files`: A directory containing files that were used to compute runtimes, max memory usage values, and improvement ratio bootstrapped confidence intervals.
  * `unit_testing_files`: A directory containing bash scripts to unit test code correctness on the ACB subpopulation, as well as the output for those tests.

### Folders in the repository :file_cabinet:

  * `images`: contains figures that are either displayed or linked to in this github README
  * `paper`: A directory containing the paper's md file, bib file, and figure
  * `thinning_methods`: All code that was used to select 500000 SNPs from the 1000 genomes project's genotype data

### Files :file_folder:

  * `regens.py`: the main file that runs the regens algorithm
  * `regens_library.py`: functions that the regens algorithm uses repeatedly. 
  * `regens_testers.py`: functions used exclusively for correctness testing and unit testing
  * `setup.py` and `_init_.py`: allows regens to be installed with pip
  * `requirements.txt`: lists REGENS' dependencies
  * `regens_tests_info.md`: Installing REGENS also downloads four folders that test REGENS' functionality. The [regens_tests_info.md](https://github.com/EpistasisLab/regens/blob/main/regens_tests_info.md) file explains what they test.   
  
## Remedies to known permission issues :adhesive_bandage:

Try these steps if you had permissions issues with the [final installation step](https://github.com/EpistasisLab/regens/blob/main/README.md#instructions-to-installing-regens-hammer_and_wrench):

1. Right click the "Anaconda Prompt" app (left), then click ```run as administrator```. Reinstalling conda in a different directory may fix this issue permenantly.
2. Your antivirus software might block a file that Anaconda needs (Avast blocked Miniconda's python.exe for us). Try seeing if your antivirus software is blocking anything related to anaconda, and then allow it to stop blocking that file. You could also turn off your antivirus software, though we do not recommend this.  
3. In the worst case, you can download all of the required files with these links:

    1. [input_files](https://ndownloader.figshare.com/files/25515740)
    2. [correctness_testng_ACB, correctness_testng_GBR, examples, runtime_testing, unit_testing_files](https://ndownloader.figshare.com/files/25516322)
    3. [hg19 and hg38](https://ndownloader.figshare.com/articles/13210796/versions/1)
    
Download the three folders containing the aforementioned 8 folders and unzip all folders (only the _folders_, keep the _recombination maps_ in hg19 and hg38 zipped). Then place everything in your working directory and run REGENS from your working directory. You should now be ready to use REGENS in your working directory if you have completed the installation steps. 

## Contributing :thumbsup:
If you find any bugs or have any suggestions/questions, please feel free to [post an issue](https://github.com/EpistasisLab/regens/issues/new)! 
Please refer to our [contribution guide](CONTRIBUTING.md) for more details.
Thanks for your support!

## License
MIT + file LICENSE
## Summary

There is one file named `run_all_tests.sh` in REGENS' repository. If you run REGENS in an anaconda envrionment named `regens` (follow the installation instructions in the README), then running the command `bsub < run_all_tests.sh` will run the majority of tests that regens downloads. There are four testing folders: 

1. correctness_testing_ACB
2. correctness_testing_GBR
3. runtime_testing
4. unit_testing_files

WARNING: the runtime_testing files are not ready to run. They all contain the line `#BSUB -m lambda[i]`, where i is the node. Change this to your node names and make sure that all nodes that you use have identical processors. 

REQUIREMENTS: running `run_all_tests.sh` submits 3 parallel jobs. The correctness_testing jobs take roughly 4.5GB of RAM 2.5 minutes each to complete. The unit_testing jobs take less than 1GB of RAM and 30 seconds to complete. If you want to run one job at a time, then you should run the following files consecutively (paths are all relative to your working directory):

1. correctness_testing_ACB/regens_automated_tests_part1_ACB.sh
2. correctness_testing_GBR/regens_automated_tests_part1_GBR.sh
3. unit_testing_files/unit_tests.sh

Note: all tests take a relatively long time to run, and in general, should not be ran if using regens to simulate a large quantity of data. Tests that `run_all_tests.sh` will run are described below. 

## Testing REGENS' correctness

The tests in `correctness_testing_ACB` and `correctness_testing_GBR` ensure that the most important parts of the regens algorithm work correctly for two different samples. Specifically, it tests four different things:

1. For all possible breakpoint intervals (to which dataset SNPs are not yet assigned) that the expected number of breakpoints in each interval is, on average, equal to the observed number of breakpoints in each interval. It does this for every chromosome, and names the output png file for the ith chromosome `expected_vs_actual_breakpoint_counts_full_range_chr[i].png`. Files corresponding to the ACB and GBR profiles get written into the corresponding correctness_testing folders. 
2. It checks that no breakpoint location is drawn twice for any simulated individual (Each breakpoint can only be drawn once per individual). 
3. It checks that every SNP being used as a breakpoint boundary matches back to the genomic interval in which it resides.
4. It checks that every filled segment of a simulated genome is equivalent to the correct genomic segment of the correct individual. 

Steps 1 and 2 show that breakpoints are drawn from the correct distribution, step 3 shows that drawn breakpoints are mapped to the correct input dataset SNPs to use as boundaries, and step 4 shows that the genotypes inserted into the resulting empty segments are the correct genotypes. If all of these things are true, then the regens algorithm probably works correctly. PySnpTools may contain errors, but this is unlikely. The output of part one is shown in one png file per chromosome. The output of parts 2, 3, and 4 are print statements that get written into the .out file from running `regens_automated_tests_part1_ACB.sh` and `regens_automated_tests_part1_GBR.sh`

If you install plink and cython, then you can also run `regens_automated_tests_optional_part2_LD_getter_GBR.sh.` and `regens_automated_tests_optional_part2_LD_getter_ACB.sh` (both in their respective folders). Doing so will produce the figures `GBR_real_vs_sim_r_val_maf_comparision.png` and `GBR_real_vs_sim_r_vs_distance_profile_comparison.png`. The first figure shows that samples simulated by regens have the same SNP correlations and maf values (with some random noise) as are measured in the real dataset. The second figure shows that the relationship between SNP correlation strength and distance between SNPs is the same for both the real and simulated data. These tests are not included in `run_all_tests.sh` because of the nevessary extra software that needs to be downloaded. 

## Unit tests for REGENS

The tests in `unit_testing_files` confirm, for the first chromosome in the 1000 genomes ACB population only, that actual output is exactly equal to correct reference output at a specific random seed. As such, testing regens with the unit tests is much faster than using the correctness tests, and regens' functionality really doesn't change across populations or chromosomes. These unit tests confirm that the following intermediate output objects equal what they should:

1. correct centimorgans_to_probabilities function output
2. correct choice_with_periodic_replacement function output
3. correct draw_breakpoints function output
4. correct get_samples_fast_breakpoint_interval_minor_allele_counts output (i.e. the total number of minor alleles per segment)
5. correct get_samples_fast_simulated_individual_minor_allele_counts output (i.e. the total number of minor alleles per simulated individual)
6. correct get_samples_fast_SNP_minor_allele_counts output (i.e. the total number of minor alleles per SNP)
7. correct reduce_recomb_rate_info function output
8. correct SNP_positions_to_rcmb_intervals function first output
9. correct SNP_positions_to_rcmb_intervals function second output

Note that, with 4, 5, and 6, the imported genotypes are too large to check for equality directly, so we compare the minor allele counts summed over a single dimension for all three dimensions. It is exceedingly improbable that all of these counts will match perfectly if the numpy arrays themselves do not. All outputs are print statements, so they're written into the .out file from running `unit_tests.sh`

## Efficiency tests for REGENS

The `runtime_testing` folder contains scripts that simply run regens 10 times in the same way. The jobs being run are started by running the following files:

1. runtime_testing/regens1.sh
2. runtime_testing/regens2.sh
3. runtime_testing/regens3.sh
4. runtime_testing/regens4.sh
5. runtime_testing/regens5.sh
6. runtime_testing/regens6.sh
7. runtime_testing/regens7.sh
8. runtime_testing/regens8.sh
9. runtime_testing/regens9.sh
10. runtime_testing/regens10.sh

Each file above contains the specification `#BSUB -m lambda[i]`, which should be removed before use on your system. The `runtime_testing/regens_main.sh` file runs all ten files listed above simultaneouslty. Each job takes requires roughly 10 GB of ram and takes 5 minutes to complete, so these jobs are not included in `run_all_tests.sh`.

## Efficiency tests for Triadsim

The `runtime_testing` folder also contains scripts that run triadsim 10 times in the same way:

1. runtime_testing/triadsim1.sh
2. runtime_testing/triadsim2.sh
3. runtime_testing/triadsim3.sh
4. runtime_testing/triadsim4.sh
5. runtime_testing/triadsim5.sh
6. runtime_testing/triadsim6.sh
7. runtime_testing/triadsim7.sh
8. runtime_testing/triadsim8.sh
9. runtime_testing/triadsim9.sh
10. runtime_testing/triadsim10.sh

You need to install triadsim before attempting to run them. Each replicate requires over 6 hours and 50GB of RAM to run. The `runtime_testing/triadsim_main.sh` file runs all ten files listed above simultaneouslty. Triadsim's runtime distribution is skewed. Roughly 1 or 2 out of every 10 triadsim runs will take roughly double its normal ammount of time or longer to run.  
Thank you for contributing to REGENS! :thumbsup: 

We welcome you to check the [existing issues](https://github.com/EpistasisLab/regens/issues) for bugs or enhancements to work on. 
If you find any bugs or have any questions/ideas for an extension to REGENS, please submit a [new issue](https://github.com/EpistasisLab/regens/issues/new) so we can discuss it.
More detailed guidelines below!

## How to report a bug

When submitting an issue, please make sure to answer these questions:

1. What version of REGENS are you using?
2. What operating system and processor architecture are you using?
3. What did you do?
4. What did you expect to see?
5. What did you see instead?

## How to submit a contribution

1. Create your own fork of this repository
2. Make the changes in your fork on a branch that is NOT master or main branch.
3. If you think the project would benefit from these changes:
    * Make sure you have followed the guidelines above.
    * Submit a pull request.

## Repository structure
Please be sure to familiarize yourself with the repository structure before making any major contributions.

### Folders :file_cabinet:

  * `correctness_testing_ACB`: A directory containing bash scripts to test code correctness on the ACB subpopulation, as well as the output for those tests. Correctness testing part 2 is optional and requires plink version 1.90Beta.
  * `correctness_testing_GBR`: A directory containing bash scripts to test code correctness on the GBR subpopulation, as well as the output for those tests. Correctness testing part 2 is optional and requires plink version 1.90Beta.
  * `examples`: A directory containing bash scripts that run the data simulation examples in the README.
  * `hg19`: for each 1000 genomes project population, contains a folder with one gzipped recombination rate dataframe per hg19 reference human autosome.
  * `hg38`: for each 1000 genomes project population, contains a folder with one gzipped recombination rate dataframe per hg38 reference human autosome.
  * `images`: contains figures that are either displayed or linked to in this github README
  * `input_files`: contains examples of regens input that is meant to be provided by the user. The example custom recombination rate information is copied from that of the hg19 mapped ACB population. Also contains input for the Triadsim algorithm. The genetic input labeled as "not_trio" for Triadsim is comprised of ACB population duplicates and is only meant to compare Triadsim's runtime. 
  * `paper`: A directory containing the paper's md file, bib file, and figure. 
  * `runtime_testing_files`: A directory containing files that were used to compute runtimes, max memory usage values, and improvement ratio bootstrapped confidence intervals.
  * `unit_testing_files`: A directory containing bash scripts to unit test code correctness on the ACB subpopulation, as well as the output for those tests.

### Files :file_folder:

  * `regens.py`: the main file that runs the regens algorithm
  * `regens_library.py`: functions that the regens algorithm uses repeatedly. 
  * `regens_testers.py`: functions used exclusively for correctness testing and unit testing

## Your first contribution

If you haven't contributed to open source code before, check out these friendly tutorials:

 * http://makeapullrequest.com/
 
 * http://www.firsttimersonly.com/
 
 * [How to contribute to an open source project on GitHub](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github).

These guides should provide you with everything you need to start out!
---
title: 'REGENS: an open source Python package for simulating realistic autosomal genotypes'
tags:
  - genomics
  - data simulation
  - bioinformatics
  - python
authors:
  - name: John T. Gregg
    orcid: 0000-0002-2619-3440
    affiliation: 1
  - name: Trang T. Le
    orcid: 0000-0003-3737-6565
    affiliation: 1
  - name: Jason H. Moore
    orcid: 0000-0002-5015-1099
    affiliation: 2
affiliations:
 - name: Department of Biostatistics, Epidemiology and Informatics, University of Pennsylvania, Philadelphia, PA 19087, USA
   index: 1
 - name: Institute for Biomedical Informatics, University of Pennsylvania, Philadelphia, PA 19087, USA
   index: 2
date: 01 October 2020
bibliography: paper.bib
---

# Summary

REcombinatory Genome ENumeration of Subpopulations (REGENS) is an open
source Python package that simulates autosomal genotypes by
concatenating real individuals' genomic segments in a way that preserves
their linkage disequilibrium (LD), which is defined as statistical
associations between alleles at different loci [@slatkin2008linkage].
Recombining segments in a way that preserves LD simulates autosomes that
closely resemble those of the real input population [@source:1]
because real autosomal genotypes can be accurately modeled as genomic
segments from a finite pool of heritable association structures (LD
haplotypes) [@source:3]. REGENS can also simulate mono-allelic and
epistatic single nucleotide variant (SNV) effects of any order without
perturbing the simulated LD pattern. The SNVs involved in an effect can
contribute additively, dominantly, recessively, only if heterozygous, or
only if homozygous. All simulated effects contribute to the value of
either a binary or continuous biological trait (phenotype) with a
specified mean value and a specified amount of random noise.

# Statement of need

The goal of most genome-wide association studies (GWAS) is to identify
associations between single nucleotide variants (SNVs) and a phenotype
to inform researchers and clinicians about potentially causative genetic
factors. Completing this task will require overcoming numerous
challenges such as insufficient sample sizes and over-representation of
European ancestries [@torkamani2018personal]. Computational biologists
build machine learning models that look for genetic associations in such
unconventional datasets, but the majority of genetic associations have
yet to be discovered [@source:4]. Researchers can use simulated datasets
with known ground truths to assess the effectiveness of an algorithm,
such as the power to detect epistatic effects with dimensionality
reduction techniques [@source:11]. The more closely simulated data
matches real-world data, the more accurate such test results will be.
Since humans of different ancestry have different LD patterns
[@source:5], a simulation that can replicate those patterns from a small
number of real samples is desirable. Therefore, intended users of REGENS
are computational biologists who aim to test a statistical learning
model on simulated GWAS data with precise realistic LD patterns.

# Algorithm overview

Two genomic segments are said to be in low LD if alleles are
approximately uncorrelated between the two segments, which is guaranteed
to occur if the boundary separating the segments has a sufficiently high
recombination rate. If two genomic segments from randomly sampled
individuals are concatenated in silico at a boundary with a high
recombination rate (the position of which is referred to as a breakpoint
from here on), then the LD pattern of the resultant in-silico autosomal
genotypes will change minimally [@source:1]. To illustrate this point,
let us let $P(R_i = 1)$ be the probability of *observing* a
recombination event at the $i^{th}$ genomic position. The following
holds: \begin{equation}\label{eq:e_ri}
P(R_i = 1) = 1 \times P(R_i = 1) + 0 \times P(R_i = 0) = E[R_i],
\end{equation}
hence, \begin{equation}\label{eq:e_frac} 
\frac{P(R_i = 1)}{\sum_{i} P(R_i = 1)} = \frac{E[R_i]}{\sum_{i} E[R_i]}.
\end{equation}
Drawing simulated breakpoints from the right hand side of (\autoref{eq:e_frac})
is like drawing differently colored marbles from a jar. Just as the
color composition inferred from drawing (with replacement) a marble from
a jar many times approaches the true distribution of colors, the
sample of simulated segment recombinations learned from drawing breakpoints
for many simulated individuals approaches the input population's empirical
distribution of real recombination events. Genomic segments that only contain
alleles in high LD are rarely separated by breakpoints, which retains
the original LD pattern (\autoref{fig:tsne}).

![Comparison of population whole genomes in 2 dimensional TSNE space.\label{fig:tsne}](tsne.png)

# Differentiating attributes

Many packages have been built to simulate genetic data with different goals
in mind. Genetic Architecture Model Emulator for Testing and Evaluating
Software (GAMETEs) simulates simple and epistatic SNV/phenotype
associations quickly but ignores LD patterns [@source:6]. Genome
Simulation of Linkage and Association (GenomeSIMLA) uses forward time
simulation to produce broadly realistic LD patterns. However, these
patterns do not exactly match those of a particular dataset [@source:7].
Triadsim [@source:1] replicates exact LD patterns, but it requires
(mother, father, kin) trios and takes an average CPU-time of 6.8 hours
and an average peak RAM of 54.6 GB to simulate 10000 trios (20000
unrelated GWAS samples) with 4 breakpoints. REGENS uses the same
recombination principles that Triadsim relies on, but it is 88.5 times
faster (95% CI (75.1, 105.0) via bootstrapping) and requires 6.2 times
lower peak RAM (95% CI (6.04, 6.33) via bootstrapping) on average over
10 replicate simulations (Intel(R) Xeon(R) CPU E5-2690 v4 2.60GHz
processor). REGENS also recombines individuals instead of trios to
simulate GWAS data with small publicly available genomic datasets, such
as those in the 1000 Genomes project. This fact allows REGENS to
accurately simulate the full genetic diversity of the world's population
(representative figures are in the supplementary analysis). Finally, REGENS can simulate
continuous and binary phenotypes that depend on any linear combination
of products of f(SNV) values, where f transforms the standard SNP values
of $\{0,1,2\}$ to represent nonlinear monoallelic effects (such as
dominance). Example implementations of these features are in REGENS'
GitHub repository.

# Supplementary analysis

Figures that demonstrate the similarity between real and simulated populations
for all twenty-six 1000 genomes populations, as well as the methods that were
used to create those figures, are here <https://github.com/EpistasisLab/regens-analysis>

# Inspiration and dependencies

REGENS was inspired by Triadsim's idea to draw simulated breakpoints at
locations with higher recombination rates, as well as well as by
GAMETE's objective of simulating data quickly. REGENS relies on
bed-reader, a spinoff of PySnpTools's core .bed file code [@source:8],
to optimally read re-sampled rows from plink bed files as 8
bit integers and then write the 8 bit integer simulated autosomal
genotypes into new bed files. REGENS also relies on the 1000 genomes
project's whole genomes from 26 distinct sub-populations [@source:9],
and it relies on those populations' corresponding genome-wide
sex-averaged recombination rates inferred by the pyrho algorithm
[@source:10].

# Acknowledgements

We acknowledge contributions from Carl Kadie, who developed PySnpTools,
for implementing its ability to read and write plink bed files as 8 bit
integers. This work was supported by NIH grant LM010098.

# References
