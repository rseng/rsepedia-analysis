# Download and Run TP-MAP

Download the latest version of TP-MAP from the [releases page](https://gitlab.com/ChemBioHub/tpmap/-/releases). Precompiled executable and jar files are provided for Windows and MacOS. To compile the software from source using Apache Maven, follow the instructions in the [INSTALL](https://gitlab.com/ChemBioHub/tpmap/-/blob/master/INSTALL) file.

**Important** TP-MAP requires 64-bit OpenJDK version 11 or higher to run, this can be downloaded from: [https://adoptopenjdk.net](https://adoptopenjdk.net)

**Important** Windows and/or MacOS may issue a security alert the first time TP-MAP is started stating that TP-MAP is an unrecognised app or from an unrecognised developer. Under Windows, click on "More info" and then "Run anyway". Under MacOS, open the System Preferences->Security & Privacy, and click on Allow TP-MAP to run on this computer.

TP-MAP has been tested on Windows 10 and macOS Catalina (10.15).

### To run the JAR file on Windows:

1. Install 64-bit OpenJDK 11 or higher
2. Download TP-MAP-\[RELEASE\]-jar-Windows.jar to your Desktop from the [releases page](https://gitlab.com/ChemBioHub/tpmap/-/releases) (replace \[RELEASE\] with latest version)
3. Run the following command (Windows Key + R):

    `java -Xmx64g -jar %HOMEPATH%\Desktop\TP-MAP-[RELEASE]-jar-Windows.jar`

If the above command doesn't work it may mean that your Java path is set incorrectly, this can be resolved by:

* Entering the full path for your OpenJDK java executable, e.g.

    `"C:\Program Files\AdoptOpenJDK\jdk-11.0.7.10-hotspot\bin\java" -Xmx64g -jar %HOMEPATH%\Desktop\TP-MAP-[RELEASE]-jar-Windows.jar`

* Editing the Path environment variable to ensure that the OpenJDK bin directory appears before any other instances of Java installed on your system:
    1. System->Advanced->Environment Variables..
    2. Under "System variables" select Path and press edit
    * If your AdoptOpenJDK path is not specified, then select "New" and add it to the top, e.g. "C:\Program Files\AdoptOpenJDK\jdk-11.0.7.10-hotspot\bin"
    * Otherwise select your AdoptOpenJDK path and select "Move Up" until it appears at the top

### To run the JAR file on MacOS:

1. Install 64-bit OpenJDK 11 or higher ([https://adoptopenjdk.net](https://adoptopenjdk.net))
2. Download TP-MAP-\[RELEASE\]-MacOS.jar to your desktop from the [releases page](https://gitlab.com/ChemBioHub/tpmap/-/releases) (replace \[RELEASE\] with latest version)
3. Open a terminal and run the following command:

    `java -Xmx64g -jar ~/Desktop/TP-MAP-[RELEASE]-MacOS.jar`

If the above command doesn't work it may mean that your Java path is set incorrectly, this can be resolved by:

* Entering the full path for your OpenJDK java executable, e.g.

    `/Library/Java/JavaVirtualMachines/adoptopenjdk-11.jre/Contents/Home/bin/java -Xmx64g -jar ~/Desktop/TP-MAP-[RELEASE]-MacOS.jar`

* Editing the PATH environment variable to ensure that the OpenJDK bin directory appears before any other instances of Java installed on your system. Run the command:

    `export PATH=/Library/Java/JavaVirtualMachines/adoptopenjdk-11.jre/Contents/Home/bin:$PATH`

  To make this change permanent, this command can be added to ~/.zprofile under ZSH or ~/.bash_profile under BASH
  
### Hardware Requirements

TP-MAP is designed to run on standard desktop computers and has been tested on macOS Catalina (10.15) and Microsoft Windows 10. Compilation may work under other operating systems with support for Java 11.

* Minimum System Requirements: Intel Core i3 1 GHz, 8 GB RAM
* Recommended System Requirements: Intel Core i5 3 GHz, 16 GB RAM

# Quick Start

Execute the TP-MAP application on your system.

The "Import Protein Abundance" button located at the top can be used to load example datasets located in the [datasets](https://gitlab.com/ChemBioHub/tpmap/-/tree/master/datasets) directory. Default options are set so that example datasets can be loaded within 5 minutes on a typical desktop computer.

## Import Panel

![Import panel](https://gitlab.com/ChemBioHub/tpmap/-/blob/master/src/main/java/com/chembiohub/tpmap/doc/import.png)

The import panel contains the following options for importing data.

**Select Experiment Type**

Choose between 1D and 2D thermal profiling data sets.

**Select Input Format**

Choose between either TP-MAP native format or Proteome Discoverer output with a configuration file. Format methods are described below.

**Normalization Method**

Choose between median normalization or no normalization on the data.

**Data file / Configuration file**

Load tab separated text files containing the data. For Proteome Discoverer format an additional configuration file is required to indicate temperature and concentration values associated with each sample and channel.

**2D Bootstrap Analysis**

Select checkbox to perform a bootstrap analysis for 2D data sets. The bootstrap analysis will permute the data for 1,000,000 iterations (by default) to empirically calculate p-values.

**1D Curve Fitting**

Select checkbox to perform curve fitting for 1D data sets. By default the software will attempt 100 times to successfully fit a curve and run the optimization algorithm for up to 10,000 iterations. Try increasing these values if curves are failing to converge.

**Multithreading**

Select checkbox to utilize all available threads/CPUs when performing the 2D bootstrap analysis or 1D curve fitting.

## Input Format

TPMAP can import protein abundance levels either from TPMAP formatted tab-separated text files or from ProteomeDiscoverer output in combination with an additional configuration file in which temperature and concentration values are specified for respective TMT samples and channels.

The native TPMAP formatted tabular text format should be formatted with headers labelled “Accession”, “Description”. If the “Description” is formatted in UniProtKB format, e.g.

> Db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier GN=GeneName PE=ProteinExistence SV=SequenceVersion

Then the description column will be parsed to extract gene name, description, and organism.

Additional columns with headers beginning with “Ref_” specify abundance values.

For 1D experiments, the abundance headers should be of the form “Ref_[Temperature]_[V/T Replicate]”, where [Temperature] is replaced with respective temperature values in decimal format, and [V/T Replicate] is replaced with a string referring to vehicle/treatment  and replicate number.

For 2D experiments, abundance headers should be of the format “Ref_[Temperature]_[Concentration]”, where [Temperature] and [Concentration] are replaced with decimal values referring to temperature and concentration, respectively.

Output from ProteomeDiscoverer can be imported into TPMAP using a configuration file that specifies respective concentration and temperature values for respective TMT samples and channels. The ProteomeDiscoverer output should be exported to a tab-separated text file with columns including “Accession”, “Description” and “Abundances: [Sample], [Channel]” where [Sample] and [Channel] are replaced with respective sample and channel identifiers.

For importing Proteome Discoverer format, an additional tab-separated text configuration file should be provided containing four columns referring to sample, channel, temperature, and either replicate (1D) or concentration (2D). For example, for one sample (F1), ten channels (126-131), two temperatures (37 and 73) and 5 concentrations.


## 2D Thermal Profiling Analysis

![2D Thermal Profiling Analysis](https://gitlab.com/ChemBioHub/tpmap/-/blob/master/src/main/java/com/chembiohub/tpmap/doc/2danalysis.png)

The analysis window contains a left panel displaying properties of the 2D thermal profiling data set, a central panel displaying a list of all proteins and the abundance-dependent thermal stability matrix when a protein is selected, and a right panel including options for down-stream analysis of selected proteins.

The the central table listing all proteins contains the following columns:

**Checkbox**

The checkbox can be used to select proteins for downstream analysis (described below).

**Accession, Gene, Organism, Description**

Descriptors for each protein extracted from the input file.

**Score**

The score ranges from +1 to -1, where +1 indicates that a protein is stabilized and -1 indicates that a protein is destabilized. By default, proteins are ranked by this metric from positive to negative, followed by Mean FC if two proteins have the same score.

**Mean FC**

Mean dose-dependent fold change of all values relative to the base concentration.

**P Value**

The probability of the score based on a bootstrap analysis of the data.

**Effect**

Takes values of "Stabilized", "Destabilized" or "Solubility/Expression", depending on whether the observed effect is thermal stability (positive score), thermal destability (negative score), or if the effect can be explained by changes in solubility or expression (when change is observed at the lowest temperature).

**MD**

The mean difference from the selected protein, calculated after initiating the mean difference analysis (described below).

When a protein is selected, the lower segment displays a table showing all dose-dependent fold changes relative to the base concentration, which is colored with a gradient from green to yellow to red, indicating stability to destability, respectively.

The filter box can be used to search all fields, and has support for common regular expressions.

## 1D Thermal Profiling Analysis

![1D Thermal Profiling Analysis](https://gitlab.com/ChemBioHub/tpmap/-/blob/master/src/main/java/com/chembiohub/tpmap/doc/1danalysis.png)

The analysis window contains a left panel displaying properties of the 1D thermal profiling data set, a central panel displaying a list of all proteins and the temperature-dependent thermal stability matrix as well as melting curves when a protein is selected, and a right panel including options for down-stream analysis of selected proteins.

The the central table listing all proteins contains the following columns:

**Checkbox**

The checkbox can be used to select proteins for downstream analysis (described below).

**Accession, Gene, Organism, Description**

Descriptors for each protein extracted from the input file.

**Tm V1, V2, T1, T2**

Tm values calculated from the dose-response curves for Vehicle 1 and 2, and Treatment 1 and 2, respectively. This is the temperature at which 50% abundance is observed.

**Tm shift V1T1, V2T2, V1V2, Mean Tm shift**

Tm shift calculated as change in Tm between vehicle 1 and treatment 1, vehicle 2 and treatment 2, vehicle 1 and vehicle 2, and mean vehicle and mean treatment, respectively.

**RMSE V1, V2, T1, T2**

Root mean squared error of the curve fit for Vehicle 1 and 2, and treatment 1 and 2, respectively.

**Shift same direction**

Boolean value, true if the shift from vehicle to treatment is the same for both replicates.

**Delta VT > Delta VV**

Boolean value, true if the mean TM shift between vehicle and treatment is greater than the TM shift between two vehicles.

**Score**

Score assessing the quality of the fit and TM shift, calculated as:  

Score = 2 * Rank[Tm] + Rank[RMSE] + Rank[Rep]  

Where RankTm is a vector (TmVT1, TmVT2) containing the scaled rank between [0,1] of increasing thermal melting point shift between vehicle and treatment in the two replicates respectively, RankRMSE is a vector (RMSEV1, RMSEV2, RMSET1, RMSET2) containing the scaled rank between [0,1] in decreasing order of root mean squared error of vehicle replicates 1 and 2, and treatment replicates 1 and 2, respectively, and RankRep is a vector (RepV and RepT) containing the scaled rank between [0,1] of the mean difference between vehicle replicates and treatment replicates, respectively. The resultant score is in the range [0,10] with a higher score prioritising proteins with large Tm shifts, small RMSE of curve fits and small difference between replicates.

## Downstream Analysis

The right panel includes downstream analysis methods, this includes:

**Fold Change**

By default TPMAP uses the 80% percentile of fold change (FC) above and below 1 to determine a threshold at which stability or destability is measured, respectively. The percentile at which the cutoff is determined can be altered using the sliders. The Recalculate button can be used to recalculate scores for proteins once a desired cutoff has been selected.

**String Network Image**

Proteins selected by a checkbox will be uploaded to the String server ([https://www.string-db.org](https://www.string-db.org)) for analysis of protein-protein interactions among selected proteins. _Important_: The accession column must refer to UniProt identifiers for this function work work correctly.

**String Functional Enrichment**

Proteins selected by a checkbox will be uploaded to the String server ([https://www.string-db.org](https://www.string-db.org)) for functional enrichment analysis among selected proteins. _Important_: The accession column must refer to UniProt identifiers for this function work work correctly.

**UniProt Entry**

Load the UniProt pages for selected proteins. _Important_: The accession column must refer to UniProt identifiers for this function work work correctly.

**Mean Difference**

This will calculate the mean fold change difference from the highlighted protein to all other proteins, and rank proteins' from most similar to dissimilar.

**Export Selected Proteins**

Proteins selected by a checkbox can be exported to an Excel file or a tab-separated text file. If the "Expanded export" checkbox is ticked then the list will include the full dose-dependent fold change table for each protein.
