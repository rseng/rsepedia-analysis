# PeacoQC
Peak-based selection of high quality cytometry data

## Introduction
The PeacoQC package provides quality control functions that will check 
for monotonic increasing channels and that will remove outliers and unstable 
events introduced due to e.g. clogs, speed changes etc. during the measurement 
of your sample. It also provides the functionality of visualising the quality
control result of only one sample and the visualisation of the results of 
multiple samples in one experiment.

## Installation
You can install this package using the devtools library.

```{r}
BiocManager::install("ComplexHeatmap")
devtools::install_github("saeyslab/PeacoQC")
```
---
title: "PeacoQC"
author: "Annelies Emmaneel"
package: PeacoQC
output: BiocStyle::html_document
vignette: >
    %\VignetteIndexEntry{PeacoQC_Vignette}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(
collapse=TRUE, 
comment="#>"
)
```

## Introduction

The PeacoQC package provides quality control functions that will check 
for monotonic increasing channels, that will remove outliers and unstable 
events introduced due to e.g. clogs, speed changes etc. during the measurement 
of your sample. It also provides the functionality of visualising the quality
control result of only one sample and the visualisation of the results of 
multiple samples in one experiment.

## Installing PeacoQC

```{r setup}
# install.packages("devtools")
# devtools::install_github("https://github.com/saeyslab/PeacoQC")

library(PeacoQC)

```


## WARNING
Please be aware of the fact that this vignette will create a directory in your
current working directory. Since this package is meant to be run on multiple 
files and the figures are directly created in this directory. 
Please check the folders corresponding to the output_directory variable when 
called for in the chuncks of code.


## Standard pre-processing and quality control pipeline

The following pipeline is recommended to use for pre-processing your data.
The pipeline starts with a flowframe (flow cytometry or mass cytometry data) 
and will first remove margins, then it will compensate and transform while 
ending with the PeacoQC quality control. This will give back a list with the 
cleaned flowframe and the different parameter settings, it will also save the 
flowframe in a new fcs file and it will make a plot of the quality control for 
this sample.

Note that the remove margins functionality and compensation is not necessary 
for mass cytometry data as explained in the code below this example.

```{r, warning=FALSE, fig.show='hide'}
# Specify flowframe path and read your flowframe
fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)


# Determine channels on which quality control should be done
channels <- c(1, 3, 5:14, 18, 21)


# Remove margins
# Make sure you do this before any compensation since the internal parameters
# change and they are neccessary for the RemoveMargins function.
# If this is not possible, you can specify the internal parameters in the 
# channel_specifications parameter.

ff <- RemoveMargins(ff=ff, channels=channels, output="frame")


# Compensate and transform the data
ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
ff <- flowCore::transform(ff, 
    flowCore::estimateLogicle(ff, colnames(flowCore::keyword(ff)$SPILL)))


# Run PeacoQC and save the cleaned flowframe as an fcs file and plot the results
# of this quality control step.
peacoqc_res <- PeacoQC(
    ff=ff, 
    channels=channels, 
    determine_good_cells="all", 
    save_fcs=TRUE, 
    plot=TRUE,
    output_directory = "PeacoQC_Example1")


# Filtered flowframe is stored in peacoqc_res$FinalFF and can be used for 
# further analysis.
ff <- peacoqc_res$FinalFF


```

## Mass cytometry data

If you want to clean mass cytometry data files you should alter some parameters.
The parameter remove_zeros should be set to TRUE. The IT_limit will range 
between 0.6 and 0.65 since some channels will be more sparse than flow cytometry
results. The time_units parameter in the plot function should be also be altered
since mass cytometry data is typically measures for a longer time than the flow 
cytometry data. You should play a bit with the parameter until you don't see
the picketfencing effect anymore in the event rate plot (Top left in the 
overview).

Note that this chunck of code will not be excecuted and that no results will 
appear in your working directory. This is an example of how to work with mass
cytometry data.

```{r, eval = FALSE}
# # Example of how the code could look for mass cytometry data

ff <- read.FCS(file)

# You don't have to remove margin events or compensate the data but you
# should transform it
channels <- c(3, 5, 6:53)

ff <- transform(ff,transformList(colnames(ff)[channels],
                                    arcsinhTransform(a = 0, b = 1/5, c = 0)))

# Make sure the parameters are set correctly and that the remove_zeros variable
# is set to TRUE.
peacoqc_results <- PeacoQC(ff,
    channels=channels,
    IT_limit=0.6,
    remove_zeros=TRUE,
    time_units=50000)

```


## Large dataset

If you have a large dataset and you want to run your preprocessing pipeline
for multiple files, it is recommended to run it first on a couple of files to 
tweak your parameters. If your dataset inlcudes channels that have a sparse 
pattern (e.g. in mass data), the IT_limit should probably be increased in order 
to be more strict. If it seems that the MAD parameter removes too much, you can
also increase this to make PeacoQC less strict.

You can then run all your files with the parameters you chose.

```{r, warning=FALSE}

# Change IT_limit for one compensated and transformed file. 
# (Higher=more strict, lower=less strict)
# The fcs file should not be saved since we are still optimising the parameters
fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Determine channels on which quality control should be done
channels <- c(1, 3, 5:14, 18, 21)

# Remove margins
ff <- RemoveMargins(ff=ff, channels=channels, output="frame")

# Compensate and transform the data
ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
ff <- flowCore::transform(ff, 
    flowCore::estimateLogicle(ff, colnames(flowCore::keyword(ff)$SPILL)))


# Run PeacoQC and save the cleaned flowframe as an fcs file and plot the results
# of this quality control step.
peacoqc_res <- PeacoQC(
    ff=ff, 
    channels=channels, 
    determine_good_cells="all", 
    save_fcs=FALSE, 
    plot=TRUE,
    output_directory = "PeacoQC_Example2",
    IT_limit = 0.65)


```


Note that the next chunck of code will not generate any results. 

```{r, eval = FALSE}

# You can also change the MAD parameter to a lower value 
# (to make it more strict) or to a higher value (to make it less strict). 
# Since the MAD analysis does not remove something, this is not neccesary now.

peacoqc_res <- PeacoQC(
    ff,
    channels,
    determine_good_cells="all",
    save_fcs=FALSE,
    plot=TRUE,
    MAD=8
)


# When the correct parameters are chosen you can run the different files in 
# a for loop

for (file in files){
    ff <- flowCore::read.FCS(file)

    # Remove margins
    ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
    
    # Compensate and transform the data
    ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
    ff <- flowCore::transform(ff, 
                                flowCore::estimateLogicle(ff, 
                                colnames(flowCore::keyword(ff)$SPILL)))
    peacoqc_res <- PeacoQC(
                            ff,
                            channels,
                            determine_good_cells="all",
                            IT_limit=0.6,
                            save_fcs=T,
                            plot=T)
    }

```


## PeacoQCHeatmap

In order to get an overview on how much the quality control removed, 
the PeacoQCHeatmap allows for a visualised representation of the different 
conditions and the different percentages that were removed in different runs.

```{r}
# Find the path to the report that was created by using the PeacoQC function
location <- system.file("extdata", "PeacoQC_report.txt", package="PeacoQC")

# Make heatmap overview of the quality control run
PeacoQCHeatmap(report_location=location)

# Make heatmap with only the runs of the last test
PeacoQCHeatmap(report_location=location, latest_tests=TRUE)

# Make heatmap with row annotation
PeacoQCHeatmap(report_location=location, 
    row_split=c(rep("r1",7), rep("r2", 55)))

```



## PlotPeacoQC without quality control

The PlotPeacoQC function can also be used to only display the peaks of the 
different channels without doing any quality control. It can even be used to 
only display the measurements.

These results will appear in the PeacoQC_plots folder.

```{r, warning=FALSE, fig.show= 'hide'}
# Load in compensated and transformed flowframe

fileName <- system.file("extdata", "111_Comp_Trans.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Plot only the peaks (No quality control)
PlotPeacoQC(ff, channels, display_peaks=TRUE, prefix = "PeacoQC_peaks_")

# Plot only the dots of the file
PlotPeacoQC(ff, channels, display_peaks=FALSE, prefix = "PeacoQC_nopeaks_")


```

















% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{RemoveDoublets}
\alias{RemoveDoublets}
\title{Remove doublet events from flow cytometry data}
\usage{
RemoveDoublets(ff, channel1="FSC-A", channel2="FSC-H", nmad=4,
verbose=FALSE, output="frame")
}
\arguments{
\item{ff}{A flowframe that contains flow cytometry data.}

\item{channel1}{The first channels that will be used to determine the
doublet events. Default is "FSC-A"}

\item{channel2}{The second channels that will be used to determine the
doublet events. Default is "FSC-H"}

\item{nmad}{Bandwidth above the ratio allowed (cells are kept if their
ratio is smaller than the median ratio + \code{nmad} times the median
absolute deviation of the ratios). Default is 4.}

\item{verbose}{If set to TRUE, the median ratio and width will be printed.
Default is FALSE.}

\item{output}{If set to "full", a list with the filtered flowframe and the
indices of the doublet event is returned. If set to "frame", only the
filtered flowframe is returned. The default is "frame".}
}
\value{
This function returns either a filtered flowframe when the
\code{output} parameter is set to "frame" or a list containing the filtered
flowframe and a TRUE/FALSE list indicating the margin events. An extra column
named "Original_ID" is added to the flowframe where the cells are given their
original cell id.
}
\description{
\code{RemoveDoublets} will remove doublet events from the
flowframe based on two channels.
}
\examples{
# Read in data
fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Remove doublets
ff_cleaned <- RemoveDoublets(ff)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{PlotPeacoQC}
\alias{PlotPeacoQC}
\title{Visualise deleted cells of PeacoQC}
\usage{
PlotPeacoQC(ff, channels, output_directory=".", display_cells=2000,
            manual_cells=NULL, title_FR=NULL, display_peaks=TRUE,
            prefix="PeacoQC_", time_unit=100, ...)
}
\arguments{
\item{ff}{A flowframe}

\item{channels}{Indices of names of the channels in the flowframe that have
to be displayed}

\item{output_directory}{Directory where the plots should be generated. Set
to NULL if no plots need to be generated. The default is the working
directory.}

\item{display_cells}{The number of measurements that should be displayed.
(The number of dots that are displayed for every channel) The default is
5000.}

\item{manual_cells}{Give a vector (TRUE/FALSE) with annotations for each cell
to compare the automated QC with. The default is NULL.}

\item{title_FR}{The title that has to be displayed above the flow rate
figure. Default is NULL.}

\item{display_peaks}{If the result of \code{PeacoQC} is given, all the
quality control results will be visualised. If set to TRUE: \code{PeacoQC}
will be run and only the peaks will be displayed without any quality control.
If set to FALSE, no peaks will be displayed and only the events will be
displayed. Default is TRUE.}

\item{prefix}{The prefix that will be given to the generated png file.
Default is "PeacoQC_".}

\item{time_unit}{The number of time units grouped together for visualising
event rate. The default is set to 100, resulting in events per second for
most flow datasets. Suggested to adapt for mass cytometry data.}

\item{...}{Arguments to be given to \code{PeacoQC} if \code{display_peaks}
is set to TRUE.}
}
\value{
This function returns nothing but generates a png file in the
output_directory
}
\description{
\code{PlotPeacoQC} will generate a png file with on overview of
the flow rate and the different selected channels. These will be annotated
based on the measurements that were removed by PeacoQC. It is also possible
to only display the quantiles and median or only the measurements without
any annotation.
}
\examples{

## Plotting the results of PeacoQC

# Read in transformed and compensated data
fileName <- system.file("extdata", "111_Comp_Trans.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Define channels on which the quality control should be done and the
# plots should be made
channels <- c(1, 3, 5:14, 18, 21)

# Run PeacoQC
PeacoQC_res <- PeacoQC(ff,
    channels,
    determine_good_cells="all",
    plot=FALSE,
    save_fcs=TRUE)

# Run PlotPeacoQC
PlotPeacoQC(ff, channels, display_peaks=PeacoQC_res)

## Plot only the peaks (No quality control)
PlotPeacoQC(ff, channels, display_peaks=TRUE)

## Plot only the dots of the file
PlotPeacoQC(ff, channels, display_peaks=FALSE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{RemoveMargins}
\alias{RemoveMargins}
\title{Remove margin events of flow cytometry data}
\usage{
RemoveMargins(ff, channels,
channel_specifications=NULL, output="frame")
}
\arguments{
\item{ff}{A flowframe that contains flow cytometry data.}

\item{channels}{The channel indices or channel names that have to be checked
for margin events}

\item{channel_specifications}{A list of lists with parameter specifications
for certain channels. This parameter should only be used if the values in
the internal parameters description is too strict or wrong for a number or
all channels. This should be one list per channel with first a minRange and
then a maxRange value. This list should have the channel name found back in
\code{colnames(flowCore::exprs(ff))}. If a channel is not listed in this
parameter, its default internal values will be used. The default of this
parameter is NULL.}

\item{output}{If set to "full", a list with the filtered flowframe and the
indices of the margin event is returned. If set to "frame", only the
filtered flowframe is returned. The default is "frame".}
}
\value{
This function returns either a filtered flowframe when the
\code{output} parameter is set to "frame" or a list containing the filtered
flowframe and a TRUE/FALSE list indicating the margin events. An extra column
named "Original_ID" is added to the flowframe where the cells are given their
 original cell id.
}
\description{
\code{RemoveMargins} will remove margin events from the
flowframe based on the internal description of the fcs file.
}
\examples{
# Read in raw data
fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Define channels where the margin events should be removed
channels <- c(1, 3, 5:14, 18, 21)

# Remove margins

ff_cleaned <- RemoveMargins(ff, channels)

# If an internal value is wrong for a channels (e.g. FSC-A)

channel_specifications <- list("FSC-A"=c(-111, 262144))
ff_cleaned <- RemoveMargins(
    ff,
    channels,
    channel_specifications=channel_specifications)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{PeacoQCHeatmap}
\alias{PeacoQCHeatmap}
\title{Make overview heatmap of quality control analysis}
\usage{
PeacoQCHeatmap(report_location, show_values=TRUE, show_row_names=TRUE,
latest_tests=FALSE, title="PeacoQC report", ...)
}
\arguments{
\item{report_location}{The path to the PeacoQC report generated by
\code{PeacoQC}.}

\item{show_values}{If set to TRUE, the percentages of removed values
will be displayed on the heatmap. Default is TRUE.}

\item{show_row_names}{If set to FALSE, the filenames will not be displayed
on the heatmap. Default is TRUE.}

\item{latest_tests}{If this is set to TRUE, only the latest quality control
run will be displayed in the heatmap. Default is FALSE.}

\item{title}{The title that should be given to the heatmap. Default is
"PeacoQC_report".}

\item{...}{Extra parameters to be given to the \code{Heatmap} function
(eg. row_split)}
}
\value{
This function returns nothing but generates a heatmap that can be
saved as pdf or png
}
\description{
\code{PeacoQCHeatmap} will make a heatmap to display all the
results generated by \code{PeacoQC}. It will include the percentages of
measurements that are removed in total, by the IT method and by the MAD
method. It will also show the parameters that were used during the
quality control.
}
\examples{

# Find path to PeacoQC report
location <- system.file("extdata", "PeacoQC_report.txt", package="PeacoQC")

# Make heatmap overview of quality control run
PeacoQCHeatmap(report_location=location)

# Make heatmap with only the runs of the last test
PeacoQCHeatmap(report_location=location, latest_tests=TRUE)

# Make heatmap with row annotation
PeacoQCHeatmap(report_location=location,
              row_split=c(rep("r1",7), rep("r2", 55)))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeacoQC.R
\name{PeacoQC}
\alias{PeacoQC}
\title{Peak-based detection of high quality cytometry data}
\usage{
PeacoQC(ff, channels, determine_good_cells="all",
        plot=20, save_fcs=TRUE, output_directory=".",
        name_directory="PeacoQC_results", report=TRUE,
        events_per_bin=FindEventsPerBin(remove_zeros, ff, channels,
        min_cells, max_bins, step), min_cells=150, max_bins=500, step=500,
        MAD=6, IT_limit=0.6, consecutive_bins=5, remove_zeros=FALSE,
        suffix_fcs="_QC", force_IT=150, peak_removal = (1/3),
        min_nr_bins_peakdetection = 10, ...)
}
\arguments{
\item{ff}{A flowframe or the location of an fcs file. Make sure that the
flowframe is compensated and transformed. If it is mass cytometry data, only
a transformation is necessary.}

\item{channels}{Indices or names of the channels in the flowframe on which
peaks have to be determined.}

\item{determine_good_cells}{If set to FALSE, the algorithm will only
determine peaks. If it is set to "all", the bad measurements will be
filtered out based on the MAD and IT analysis. It can also be put to "MAD"
or "IT" to only use one method of filtering.}

\item{plot}{When PeacoQC removes more than the specified percentage, an
overview plot will be made of all the selected channels and the deleted
measurements. If set to TRUE, the \code{PlotPeacoQC} function is
run to make an overview plot of the deleted measurements, even when
nothing is removed. Default is set to 20. If an increasing or decreasing
trend is found, a figure will also be made except if plot is set to FALSE.}

\item{save_fcs}{If set to TRUE, the cleaned fcs file will be saved in the
\code{output_directory} as: filename_QC.fcs. The _QC name can be altered with
the \code{suffix_fcs} parameter. An extra column named "Original_ID" is added
to this fcs file where the cells are given their original cell id.
Default is TRUE.}

\item{output_directory}{Directory where a new folder will be created that
consists of the generated fcs files, plots and report. If set to NULL,
nothing will be stored.The default folder is the working directory.}

\item{name_directory}{Name of folder that will be generated in
\code{output_directory}. The default is "PeacoQC_results".}

\item{report}{Overview text report that is generated after PeacoQC is run.
If set to FALSE, no report will be generated. The default is TRUE.}

\item{events_per_bin}{Number of events that are put in one bin.
Default is calculated based on the rows in \code{ff}}

\item{min_cells}{The minimum amount of cells (nonzero values) that should be
present in one bin. Lowering this parameter can affect the robustness of the
peak detection. Default is 150.}

\item{max_bins}{The maximum number of bins that can be used in the cleaning
process. If this value is lowered, larger bins will be made. Default is 500.}

\item{step}{The step in events_per_bin to which the parameter is reduced to.
Default is 500.}

\item{MAD}{The MAD parameter. Default is 6. If this is increased, the
algorithm becomes less strict.}

\item{IT_limit}{The IsolationTree parameter. Default is 0.55. If this is
increased, the algorithm becomes less strict.}

\item{consecutive_bins}{If 'good' bins are located between bins that are
removed, they will also be marked as 'bad'. The default is 5.}

\item{remove_zeros}{If this is set to TRUE, the zero values will be removed
before the peak detection step. They will not be indicated as 'bad' value.
This is recommended when cleaning mass cytometry data. Default is FALSE.}

\item{suffix_fcs}{The suffix given to the new fcs files. Default is "_QC".}

\item{force_IT}{If the number of determined bins is less than this number,
the IT analysis will not be performed. Default is 150 bins.}

\item{peak_removal}{During the peak detection step, peaks are only kept if
they are \code{peak_removal} percentage of the maximum height peak. Default is
1/3}

\item{min_nr_bins_peakdetection}{The percentage of number of bins in which
the maximum number of peaks has to be present. Default is 10.}

\item{...}{Options to pass on to the \code{PlotPeacoQC} function
(display_cells, manual_cells, prefix)}
}
\value{
This function returns a \code{list} with a number of items. It will
include "FinalFF" where the transformed, compensated and cleaned flowframe is
stored. It also contains the starting parameters and the information
necessary to give to \code{PlotPeacoQC} if the two functions are run
seperatly. The GoodCells list is also given where 'good' measurements are
indicated as TRUE and the to be removed measurements as FALSE.
}
\description{
\code{PeacoQC} will determine peaks on the channels in the
flowframe. Then it will remove anomalies caused by e.g. clogs, changes in
speed etc. by using an IsolationTree and/or the MAD method.
}
\examples{
# General pipeline for preprocessing and quality control with PeacoQC

# Read in raw fcs file
fileName <- system.file("extdata", "111.fcs", package="PeacoQC")
ff <- flowCore::read.FCS(fileName)

# Define channels where the margin events should be removed
# and on which the quality control should be done
channels <- c(1, 3, 5:14, 18, 21)

ff <- RemoveMargins(ff=ff, channels=channels, output="frame")

# Compensate and transform the data

ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                            colnames(flowCore::keyword(ff)$SPILL)))
#Run PeacoQC
PeacoQC_res <- PeacoQC(ff, channels,
                        determine_good_cells="all",
                        save_fcs=TRUE)

}
