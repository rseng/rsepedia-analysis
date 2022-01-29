<p align="center">
    <img height="150" src="https://uni-muenster.sciebo.de/s/BD2q8HfArBG7juF/download">
</p>

# CopyDetective

CopyDetective is a novel algorithm for detection threshold aware CNV calling in matched whole-exome sequencing data. 

## Requirements

To run CopyDetective, you need R (Version 3.5.2 or higher) and shiny.


## Installation

To install CopyDetective, just download the R-scripts. All required packages are automatically installed.


## Running CopyDetective

CopyDetective is available as a shiny GUI. You can run CopyDetective either by RStudio or by using `shiny::runApp()` in an R console.
On the left, there are two tabs available for analysis:

- *Perform analysis*: actual CNV calling is performed
- *Display results*: display list of CNV calls and output plots for a selected sample

On the right, four output tabs are available:

- *Log*: progress of the analysis
- *Detection thresholds*: individual detection thresholds for every sample, including minimum CNV size and minimum cell fraction for deletions and duplications
- *CNV calls*: merged CNV calls; if filtration by means of quality is selected, the calls are filted based on the selected filtration threshold
- *Plots*: dependent on the selected options, up to four plots are displayed: 1) Raw CNV calls (all), 2) Raw CNV calls (sig), 3) Merged CNV calls, and 4) Filtered CNV calls

For the actual analysis with CopyDetective, several parameters are available.

### Parameters

| Category | Parameter | Details |
| ------ | ------ | ------ |
| Input | Input folder | Path to input files; input files must cover information on heterozygous polymorphisms called for every sample, including VAF and coverage of these polymorphisms in tumor and matching germline sample |
|   | Sample names file | A 2-column txt-file is required; the germline sample names must be defined in the first column, the matching tumor sample names in the seocnd column (no header)|
| Output | Output folder | Path to which output files will be written |
|| Output files | The following can be selected: Raw CNV calls, Merged CNV calls, Filtered CNV calls |
|| Output plots | The following can be selected: Raw CNV calls (all), Raw CNV calls (sig), Merged CNV calls, Filtered CNV calls |
| 1. Quality Analysis | Detection thresholds available? | Select either *yes* or *no*|
||*yes*: Upload detection thresholds file| A txt-file containing the columns: Sample, Window_del, CF_del, SNPs_del, Window_dup, CF_dup and SNPs_dup; an example is available in the *example* folder
||*no*: How shall Quality Analysis be performed? | Select either *Exact* or *Simulation*|
||*no* and *Exact*: VAF for polymorphisms in the control sample | Select either *Exact (variant calling results)* or *Simulated (expected value 0.5)*|
||*no* and *Simulation*: Number of simulations| 10-99999; default: 500|
||*no*: Evaluate Cell Fractions: Fractions to consider| 1-100; default: 5-100|
||*no*: Evaluate Cell Fractions: In steps of...| 1-100; default: 5|
||*no*: Evaluate Window sizes: Percentile to consider for window selection | 0-100; default: 95|
||*no*: Strategy for detection threshold optimization| Select one of the following: Compromize (cell fraction and window size), Force (cell fraction), Force (window size)|
||*no*: Force (cell fraction)| Minimum (possible) cell fraction to consider | 1-100; default: 10 |
||*no*: Force (window size)| Minimum (possible) window size to consider | 100,000-59,128,983; default: 1,000,000 |
| 2. CNV calling | Desired sensitivity |0-1; default: 0.95|
| 3. Merging | Maximum allowed distance between raw calls | 100,000-59128,983; default: 20,000,000|
| 4. Filtration | Perform final filtration? | Select either *yes* or *no*|
|| *yes*: Quality threshold | 0-600; default: 10.76|

## Example

An example, containing all necessary input for two patients can be found in the *Example* folder. Corresponding output files are provided in the subfolder *output*.

## Reporting errors / Requesting features

If you should experience any errors using CopyDetective, or you are looking for some additional features, feel free to open an issue or to contact Sarah Sandmann ( sarah.sandmann@uni-muenster.de ).
