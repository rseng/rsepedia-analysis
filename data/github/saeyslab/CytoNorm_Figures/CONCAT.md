# CytoNorm_Figures
Scripts to replicate the figures reported in the CytoNorm publication. 
---
title: "R Notebook"
output: html_notebook
---

Load the required libraries

```{r}
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(flowWorkspace))
```

Choose the files of interest

```{r}
dir <- "../data/Gates_ExtraControls"

wsp_file <- "manual.wsp"
fcs_files <- list.files(dir, pattern="Gates_P.*fcs$")
```

Choose the gates of interest as final cell labels

```{r}
cell_types <- c("Naive CD4+ T cells" = "CD4+ T cells/CD45RA+ naive T cells",
                "Memory CD4+ T cells" = "CD4+ T cells/CD45RA- mem T cells",
                "TCRgd T cells" = "TCRgd T cells",
                "Naive CD8+ T cells" = "CD8+ T cells/CD45RA+ naive T cells",
                "Memory CD8+ T cells" = "CD8+ T cells/CD45RA- mem T cells",
                "NK cells" = "CD7+CD3-NKcells",
                "B cells" = "B cells",
                "Monocytes" = "lin-CD14+CD11b+CD33+MCs",
                "pDCs" = "pDCs",
                "CD66+CD45-" = "CD66+CD45-")
```

Parse the flowjo workspace

```{r}
wsp <- flowWorkspace::openWorkspace(file.path(dir, wsp_file))
o <- capture.output(
    gates <- suppressMessages(flowWorkspace::parseWorkspace(wsp, 
                                                            "All Samples", 
                                                            sampNloc = "sampleNode"))
)
files_in_wsp <- gates@data@origSampleVector
counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
```


Extract the manual gating for every file

```{r}
result <- list()
for(file in fcs_files){
  print(paste0("Processing ", file))
  
  file_id <- grep(gsub(".*/", "", file), files_in_wsp)
  if(length(file_id) == 0) {stop("File not found. Files available: ",
                                 gsub("_[0-9]*$", "\n", files_in_wsp))}
  
  # Build a matrix with TRUE/FALSE values for every cell / gate combination
  gate_names <- flowWorkspace::getNodes(gates[[file_id]], path = "auto")
  gating_matrix <- matrix(FALSE,
                         nrow = counts[file_id],
                         ncol = length(gate_names),
                         dimnames = list(NULL, gate_names))
  for (gate in gate_names) {
    gating_matrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                                                        gate)
  }
  
  # Build one vector with an individual label for every cell
  manual_vector <- rep("Unknown",nrow(gating_matrix))
  for(cell_type in cell_types){
    manual_vector[gating_matrix[,cell_type]] <- cell_type
  }
  manual_vector <- factor(manual_vector, levels=c("Unknown",cell_types))
  
  # Rename to prettier names than in the FlowJo workspace
  levels(manual_vector) <- c("Unknown", names(cell_types))
  col_ids <- sapply(cell_types, 
                    function(x) { which(colnames(gating_matrix) == x) } )
  colnames(gating_matrix)[col_ids] <- names(cell_types)
  
  result[[file]] <- list(matrix = gating_matrix,
                         manual = manual_vector)
  
}
```

Save results for later use

```{r}
saveRDS(result,
        "../data/Gates_ExtraControls/manual.RDS")
```

```{r}
sessionInfo()
```
---
title: "R Notebook"
output: html_notebook
---

# tSNE comparison 

```{r}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(FlowSOM))
suppressPackageStartupMessages(library(CytoNorm))
suppressPackageStartupMessages(library(flowWorkspace))
suppressPackageStartupMessages(library(CytoML))
```

```{r}
dir <- "D:/Data/Ongoing/Stanford/data/Gates_ExtraControls"
files <- list.files(dir, pattern="Gates_P.*fcs$")
df <- data.frame(File = files,
                 Plate = as.factor(str_match(files, "(PTLG[0-9]*)")[, 2]),
                 Stim = as.factor(str_match(files, "[0-9]_(.*)_Control")[, 2]),
                 Volunteer = as.factor(str_match(files, "_([12]).")[, 2]),
                 stringsAsFactors = FALSE)
rownames(df) <- df$File
df <- df %>% dplyr::filter(Stim == "Unstim" & Volunteer == "2")
df
```

```{r}
set.seed(1)

data <- list("Original"  = AggregateFlowFrames(file.path(dir, df$File), 10000),
             "Normalized" = 
               AggregateFlowFrames(
                 file.path("results_CytoNorm/data/Cluster_QuantileNormLimited_Unstim_and_IFNa_LPS_1", 
                           paste0("Norm_", df$File)),
                 10000))
```


```{r}
cols_of_interest <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                      11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
                      41, 26, 30, 28, 29, 25, 35)
```

```{r}
data <- lapply(data, function(d){transform(d, transformList(colnames(d)[cols_of_interest],
                                                                      cytofTransform))})
```

```{r}
set.seed(1)

dim_red <- list()
for(d in names(data)){
  dim_red[[d]] <- Rtsne::Rtsne(data[[d]]@exprs[, cols_of_interest],
                                perplexity = 30)
}
```

```{r fig.width = 15}
set.seed(1)
rand_order <- sample(nrow(data$Original))

plots <- list()
for(d in names(data)){
  plots[[d]] <- ggplot(data.frame(x = dim_red[[d]]$Y[rand_order,1], 
                                  y = dim_red[[d]]$Y[rand_order,2],
                                  File = factor(data[[d]]@exprs[rand_order,"File"]))) +
    geom_point(aes(x, y, color = File), 
               size = 0.5) + 
    xlab("tSNE 1") + 
    ylab("tSNE 2") +
    theme_minimal() +
    ggtitle(d) +
    theme(legend.position = "none")
}

p <- do.call(gridExtra::grid.arrange, c(plots, list(ncol = 2)))

ggsave("results_CytoNorm/tSNE_before_after.pdf",
       p,
       width = 10, height = 5)

```


---
title: "R Notebook"
output: html_notebook
---

Install libraries if not yet previously installed.
Remove hash sign to actually run the lines.
```{r}
# install.packages("tidyverse")
# install.packages("gridExtra")
# install.packages("pheatmap")

# install.packages("BiocManager")
# BiocManager::install("flowCore")

# install.packages("devtools")
# devtools::install_github("saeyslab/FlowSOM")
# devtools::install_github("saeyslab/CytoNorm")
```

## Setup

Load the required libraries

```{r}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(FlowSOM))
suppressPackageStartupMessages(library(CytoNorm))
```

Directory to save results to

```{r}
results <- "results_CytoNorm"
if (!dir.exists(results)) dir.create(results)
```

Set the path to the data files and extract the metadata from the filenames.
These files can be downloaded from flow repository dataset FR-FCM-Z247.

```{r}
dir <- "D:/Data/Ongoing/Stanford/data/Gates_ExtraControls"
files <- list.files(dir, pattern="Gates_P.*Control.*fcs$")
samples <- data.frame(File = files,
                      Plate = as.factor(str_match(files, "(PTLG[0-9]*)")[, 2]),
                      Stim = as.factor(str_match(files, "[0-9]_(.*)_Control")[, 2]),
                      Volunteer = as.factor(str_match(files, "_([12]).")[, 2]),
                 stringsAsFactors = FALSE)
rownames(samples) <- samples$File
samples
```

Read the marker names from the first file.  
Identify the markers of interest, ordered alphabetically, surface markers first.
```{r}
o <- capture.output(ff <- read.FCS(file.path(dir, samples$File[1])))
marker_names <- get_markers(ff, colnames(ff))

norm_ids <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
              11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
              41, 26, 30, 28, 29, 25, 35)
norm_channels <- colnames(ff)[norm_ids]
norm_markers <- marker_names[norm_ids]

transformList <- transformList(norm_channels, cytofTransform)
transformList.reverse <- transformList(norm_channels, cytofTransform.reverse)

print(norm_markers)
```

Load the manual mapping of the data.  
See ExtractManualGating_CytoNorm.R file to generate the RDS file or download
from flowRepository.
```{r}
manual <- readRDS("D:/Data/Ongoing/Stanford/data/Gates_ExtraControls/manual.RDS")
cell_types <- c("All", levels(manual[[1]]$manual)[-1])
print(cell_types)
```

## Characterize the batch effects

### Density overview

```{r}
# Functions which use global variables, only for use in this notebook!

getDensities_plot <- function(densities, markers){
  
  densities <- lapply(densities, function(x){
    x$File <-  gsub("Norm_", "", gsub(".*/", "", x$File))
    x$Batch <- samples[x$File, "Plate"]
    x$Volunteer_Stim <- paste0(samples[x$File, "Volunteer"],
                                     "_",
                                     samples[x$File, "Stim"])
    x$Volunteer_Stim <- factor(x$Volunteer_Stim,
                               levels = c("1_Unstim", "1_IFNa_LPS",
                                          "2_Unstim", "2_IFNa_LPS"))
    x
  })
  
  plots <- list()
  for(marker in markers){
    p <- ggplot(dplyr::filter(densities$quantiles,
                              Channel == FlowSOM::get_channels(ff, marker))) +
      geom_vline(aes(xintercept = `0.05`), col = "grey", linetype="dashed") +
      geom_vline(aes(xintercept = `0.25`), col = "grey") +
      geom_vline(aes(xintercept = `0.5`), col = "#b30000", lwd = 1) +
      geom_vline(aes(xintercept = `0.75`), col = "grey") +
      geom_vline(aes(xintercept = `0.95`), col = "grey", linetype="dashed") +
      geom_line(aes(x = x, y = y),
                data = dplyr::filter(densities$densities,
                              Channel == FlowSOM::get_channels(ff, marker))) +
      facet_grid(Batch ~ Volunteer_Stim,
                 switch="y") +
      ggtitle(paste0(marker, " expression")) +
      xlab("Control and validation samples included in every plate") +
      ylab("Different plates") +
      xlim(c(-0.5, 6)) +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            strip.text.y = element_text(angle = 180))
    
    plots[[marker]] <- p
  }
  return(plots)
}
```

```{r}
res_file <- file.path(results, "densities_original.RDS")
if(!file.exists(res_file)){
  densities <- getDensities(files = file.path(dir, samples$File),
                            channels = norm_channels,
                            transformList = transformList,
                            quantileValues =  c(0.05, 0.25, 0.5, 0.75, 0.95))
  saveRDS(densities, res_file)
} else {
  densities <- readRDS(res_file)
}
```

```{r}
markers <- c("CD15", "CD66")
plots <- getDensities_plot(densities,
                           markers = c("CD15", "CD66"))

p <- do.call(gridExtra::grid.arrange, c(plots, list(nrow = 1)))

ggsave(p, 
       filename = file.path(results,
                            paste0("Density_original_",
                                   paste(markers, collapse = "_"),".pdf")),
       width = 14, height = 7)
```


```{r}
correlations <- matrix(NA, 
                       nrow = 2,
                       ncol = length(norm_markers),
                       dimnames = list(c("Unstim", "IFNa_LPS"),
                                       norm_markers))
for(channel in norm_channels){
  subset <- densities$quantiles %>% 
    dplyr::filter(Channel == channel) 
  sample_types <- str_match(subset$File, "PTLG[0-9]*_(.*).fcs")[,2]
  medians <- subset$`0.5`
  correlations["Unstim", marker_names[channel]] <- cor(medians[sample_types == "Unstim_Control_1"],
                                                       medians[sample_types == "Unstim_Control_2"])
  correlations["IFNa_LPS", marker_names[channel]] <- cor(medians[sample_types == "IFNa_LPS_Control_1"],
                                                           medians[sample_types == "IFNa_LPS_Control_2"])
}

apply(correlations, 1, mean)
apply(correlations, 1, sd)
```

### Cell type specific densities

```{r fig.width = 14, fig.height = 7 }
cell_types_to_plot <- c("B cells", "Monocytes")
marker <- "HLADR"

res_file <- file.path(results, 
                          paste0("densities_original_",
                                 paste(cell_types_to_plot, collapse = "_"),
                                 ".RDS"))
if (!file.exists(res_file)) {
  
  densities_selected <- list()
  
  for (cell_type in cell_types_to_plot) {
    cell_type_selection <- lapply(manual, function(x){
      x$matrix[, cell_type]
    })
    
    names(cell_type_selection) <- file.path(dir, names(cell_type_selection))
    
    densities_selected[[cell_type]] <- 
      getDensities(files = file.path(dir, samples$File),
                   channels = norm_channels,
                   transformList = transformList,
                   quantileValues =  c(0.05, 0.25, 0.5, 0.75, 0.95),
                   selection = cell_type_selection)
  }

  saveRDS(densities_selected, res_file)
} else {
  densities_selected <- readRDS(res_file)
}
```

```{r}
plots <- list()
for(cell_type in cell_types_to_plot){
  p <- getDensities_plot(densities_selected[[cell_type]],
                         marker)[[1]] + 
       ggtitle(paste0(marker, " expression in ", cell_type)) +
       xlim(c(2.5, 6))
  
  plots[[cell_type]] <- p
}

p <- do.call(gridExtra::grid.arrange, c(plots, list(nrow = 1)))

ggsave(p, 
       filename = file.path(results,
                            paste0("Density_original_",
                                   marker, "_",
                                   paste(cell_types_to_plot, collapse = "_"),
                                   ".pdf")),
       width = 14, height = 7)
              
```


BC densities for reviewer
```{r}
BC_channels <- colnames(ff)[grep("BC", get_markers(ff, colnames(ff)))]
densities_BC <- getDensities(files = file.path(dir, samples$File),
                             channels = BC_channels,
                             transformList = transformList(BC_channels, 
                                                           cytofTransform),
                             quantileValues =  c(0.05, 0.25, 0.5, 0.75, 0.95))

plots <- getDensities_plot(densities_BC,
                           markers = get_markers(ff, BC_channels))

p <- do.call(gridExtra::grid.arrange, c(plots, list(nrow = 2)))

ggsave(p, 
       filename = file.path(results,
                            paste0("Density_original_",
                                   paste(get_markers(ff, BC_channels), 
                                         collapse = "_"),
                                   ".pdf")),
       width = 21, height = 14)

```---
title: "R Notebook"
output: html_notebook
---

# tSNE comparison 

```{r}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(FlowSOM))
suppressPackageStartupMessages(library(CytoNorm))
suppressPackageStartupMessages(library(flowWorkspace))
suppressPackageStartupMessages(library(CytoML))
```

```{r}
set.seed(1)

data_agg <- list("Original"  = AggregateFlowFrames(validation_data$File, 11200),
             "Normalized" = 
               AggregateFlowFrames(
                 file.path(norm_dir, paste0("Norm_", gsub(".*/", "", validation_data$File))),
                 cTotal = 11200))
```


```{r}
cols_of_interest <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                      11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
                      41, 26, 30, 28, 29, 25, 35)
```

```{r}
data_agg_t <- lapply(data_agg, function(d){transform(d, transformList(colnames(d)[cols_of_interest],
                                                                      cytofTransform))})
```

```{r}
set.seed(1)

dim_red <- list()
for(d in names(data_agg)){
  dim_red[[d]] <- Rtsne::Rtsne(data_agg[[d]]@exprs[, cols_of_interest],
                                perplexity = 30)
}
```

```{r fig.width = 15}
set.seed(1)
rand_order <- sample(nrow(data_agg$Original))

plots <- list()
for(d in names(data_agg)){
  plots[[d]] <- ggplot(data.frame(x = dim_red[[d]]$Y[rand_order,1], 
                                  y = dim_red[[d]]$Y[rand_order,2],
                                  File = factor(data_agg[[d]]@exprs[rand_order,"File"]))) +
    geom_point(aes(x, y, color = File), 
               size = 1) + 
    theme_minimal() +
    ggtitle(d) +
    theme(legend.position = "none")
}

p <- do.call(gridExtra::grid.arrange, c(plots, list(ncol = 2)))

ggsave(file.path(results, "tSNE_patients_before_after.pdf"),
       p,
       width = 10, height = 5)

```---
title: "R Notebook"
output: html_notebook
---
```{r}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(flowStats))
```

```{r}
results <- "validation_Gates/GaussNorm"
if(!dir.exists(results)) dir.create(results, recursive = TRUE)

norm_dir <- file.path(results, "Normalized")
if(!dir.exists(norm_dir)) dir.create(norm_dir)


files <- c(list.files("../data/Gates", pattern = "PTLG.*Unstim.*fcs$", full.names = TRUE),
           list.files("../data/Gates_ExtraControls", pattern = "PTLG.*Unstim.*fcs$", full.names = TRUE))
files <- grep(paste0(".*", c("PTLG015", "PTLG016", "PTLG006", "PTLG014", "PTLG023"), "_[123BU]",
                     collapse = "|"),
              files, 
              invert = TRUE,
              value = TRUE) 
files <- grep("Unstim_Control_2",
              files,
              invert = TRUE,
              value = TRUE)

df <- data.frame(File = files,
                 Plate = as.factor(str_match(files, "(PTLG[0-9]*)")[, 2]),
                 Timepoint = as.factor(gsub("Unstim", "Control", str_match(files, "PTLG[0-9]*[_Repeat_]*_([^_]*)")[, 2])),
                 stringsAsFactors = FALSE)
rownames(df) <- df$File

repeat_ctrl_file <- "../data/Gates/Gates_PTLG015_016_Repeats_Unstim_Control.fcs"
df[repeat_ctrl_file, c("Plate","Timepoint")] <- c("PTLG015", "Control") 
head(df)
```

```{r}
train_data <- df %>% 
  dplyr::filter(Timepoint == "Control")

o <- capture.output(ff <- read.FCS(file.path(train_data$File[1])))

norm_markers <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                  11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
                  41, 26, 30, 28, 29, 25, 35)
channels <- colnames(ff)[norm_markers]

transformList <- transformList(colnames(ff)[norm_markers], 
                               cytofTransform)
transformList.reverse <- transformList(colnames(ff)[norm_markers], 
                                       cytofTransform.reverse)
```


Train landmarks on first control file as reference
```{r}
ff <- read.FCS(train_data %>% 
                 dplyr::filter(Plate == "PTLG001") %>% 
                 dplyr::pull(File))

ff <- transform(ff, transformList)

remb.flowset <- flowStats:::remBoundary(flowset = flowSet(ff), 
                                       channel.names = channels)
max.lms <- rep(2, length(channels))
names(max.lms) <- channels

lms <- flowStats:::extract.landmarks(flowset = remb.flowset$res, 
                                     expr.list = 1, 
                                     channel.names = channels, 
                                     max.lms = max.lms, 
                                     peak.density.thr = 0.05, 
                                     peak.distance.thr = 0.05)

base.lms <- flowStats:::extract.base.landmarks(lms$filter, 
                                               channels, 
                                               max.lms)
```

```{r}
validation_data <- df %>% 
  dplyr::filter(Timepoint != "Control")

validation_data$Cohort <- stringr::str_match(validation_data$File, 
                                             "/(Gates[^/]*)/")[,2]

rownames(validation_data) <- gsub(".*/", "", validation_data$File)

validation_data$Plate[validation_data$Plate == "PTLG016"] <- "PTLG015"
```

```{r}
for(i in seq_len(nrow(validation_data))){
  file <- validation_data$File[i]
  ff <- read.FCS(file)
  ff <- transform(ff, transformList)
  gauss_res <- gaussNorm(flowset = flowSet(ff), 
                    channel.names = channels, 
                    base.lms = base.lms)
  ff_n <- gauss_res$flowset[[1]]
  ff_n <- transform(ff_n, transformList.reverse)
  write.FCS(ff_n, 
            filename = file.path(norm_dir, 
                                 paste0("Norm_", gsub(".*/", "", file))))
}
```---
title: "R Notebook"
output: html_notebook
---

This notebook assumes you have run the setup part of
1_data_exploration_of_control_and_validation_samples.Rmd first,
as it makes use of variables created in that notebook.


## Evaluation of the FlowSOM clustering

```{r}
train_files <- samples %>% 
  dplyr::filter(Volunteer == 1) %>% 
  dplyr::filter(Stim == "Unstim") %>%
  dplyr::pull(File)
```

```{r}
res_file <- file.path(results, "prepareFlowSOM.RDS")
if (!file.exists(res_file)) {
  fsom <- prepareFlowSOM(file.path(dir, train_files),
                         norm_channels,
                         nCells = 1000000,
                         FlowSOM.params = list(xdim = 15,
                                               ydim = 15,
                                               nClus = 25,
                                               scale = FALSE),
                         transformList = transformList,
                         plot = TRUE,
                         seed = 1)
  saveRDS(fsom, file.path(results, "prepareFlowSOM.RDS"))
} else {
  fsom <- readRDS(res_file)
}
```

```{r}
# Test different amounts of metaclusters
pdf(file.path(results, "CV_test.pdf"),
    width = 11, height = 7)
cv_test <- testCV(fsom,
                  cluster_values = 3:50,
                  plot = TRUE,
                  verbose = TRUE)
dev.off()

cv_test$cvs$`25`
```

## Run normalization on the samples

### CytoNorm normalization

Experiment setup

```{r}
FlowSOM.params <- list(nCells = 1000000,
                       xdim = 15,
                       ydim = 15,
                       nClus = 25,
                       scale = FALSE)

norm.params <- list("QuantileNormLimited" = list(nQ = 101,
                                                 limit = c(0,8),
                                                 goal = "PTLG021"),
                    
                    "QuantileNorm" = list(nQ = 101,
                                          goal = "PTLG021"),
                    
                    "MinMaxNorm" = list(nQ = 2,
                                        quantile_values = c(0, 1),
                                        goal = "PTLG021"),

                    "QMinMaxNorm" = list(nQ = 2,
                                         quantile_values = c(0.001, 0.999),
                                         goal = "PTLG021"))

exp_setups <- expand.grid(clustering = c("Cluster", ""),
                          norm = names(norm.params),
                          stim = c("Unstim_and_IFNa_LPS"),
                          volunteer = c(1,2))

rownames(exp_setups) <- gsub("^_", "", 
                             apply(exp_setups, 1, paste, collapse = "_"))
exp_setups
```

Setup for testing control samples with different ranges
```{r}
FlowSOM.params <- list(nCells = 1000000,
                       xdim = 15,
                       ydim = 15,
                       nClus = 25,
                       scale = FALSE)

norm.params <- list("QuantileNormLimited" = list(nQ = 101,
                                                 limit = c(0,8),
                                                 goal = "PTLG021"))

exp_setups <- expand.grid(clustering = c("Cluster"),
                          norm = names(norm.params),
                          stim = c("Unstim_and_IFNa_LPS", "Unstim", "IFNa_LPS"),
                          volunteer = c(1,2))

rownames(exp_setups) <- gsub("^_", "", 
                             apply(exp_setups, 1, paste, collapse = "_"))

exp_setups
```

Training 

```{r}
if(!dir.exists(file.path(results, "models")))
  dir.create(file.path(results, "models"))

if(file.exists(file.path(results, "models", "times.rds"))){
  times <- readRDS(file.path(results, "models", "times.rds"))
} else {
  times <- list()
}

for(exp in rownames(exp_setups)){
  if(!dir.exists(file.path(results, "models", exp))){
    dir.create(file.path(results, "models", exp))
  }
  print(paste0("Training ", exp))
  train_data <- samples %>% 
    dplyr::filter(Volunteer == exp_setups[exp, "volunteer"]) %>% 
    dplyr::filter(grepl(gsub("_and_", "|",exp_setups[exp, "stim"]), Stim))
  
  if(exp_setups[exp, "clustering"] == "Cluster"){
    
    t <- system.time(capture.output(
      model <- CytoNorm.train(files = file.path(dir, 
                                                train_data$File),
                              labels = train_data$Plate,
                              channels = norm_channels,
                              transformList = transformList,
                              FlowSOM.params = FlowSOM.params,
                              normMethod.train = QuantileNorm.train,
                              normParams = norm.params[[exp_setups[exp, "norm"]]],
                              seed = 1,
                              outputDir = file.path(results, "models", exp),
                              plot = TRUE)))
  } else {
    t <- system.time(capture.output(
      model <- do.call(QuantileNorm.train,
                       c(list(files = file.path(dir, 
                                                train_data$File),
                              labels = train_data$Plate,
                              channels = norm_channels,
                              transformList = transformList,
                              plot = TRUE),
                         norm.params[[exp_setups[exp, "norm"]]]))))
    
  }
  
  print(t)
  times[[exp]] <- t
  saveRDS(times,
          file = file.path(results,
                           "models",
                           "times.rds"))
  saveRDS(model,
          file = file.path(results,
                           "models",
                           exp,
                           "GatesExtraControls_model.rds"))
  
}
```

Normalization

```{r}
dir.create(file.path(results, "data"))

for(exp in rownames(exp_setups)){
  dir.create(file.path(results, "data", exp))
  print(paste0("Normalizing ", exp))
  model <- readRDS(file.path(results,
                             "models",
                             exp,
                             "GatesExtraControls_model.rds"))
  test_data <- samples %>% 
    dplyr::filter(Volunteer != exp_setups[exp, "volunteer"])
  
  if(exp_setups[exp, "clustering"] == "Cluster"){
    
    t <- system.time(capture.output(
      CytoNorm.normalize(model = model,
                         files = file.path(dir, test_data$File),
                         labels = test_data$Plate,
                         transformList = transformList,
                         transformList.reverse = transformList.reverse,
                         outputDir = file.path(results, "data", exp),
                         normMethod.normalize = QuantileNorm.normalize)))
  } else { 
    
    t <- system.time(capture.output(
      QuantileNorm.normalize(model = model,
                             files = file.path(dir, test_data$File),
                             labels = test_data$Plate,
                             transformList = transformList,
                             transformList.reverse = transformList.reverse,
                             outputDir = file.path(results, "data", exp))))
    
  }
  print(t)
}
```

```{r}
emds <- list()

for(volunteer in c(1,2)){
  for(stim in c("Unstim", "IFNa_LPS")){
    print(paste0("Evaluating original"))
    toEvaluate <- samples %>% 
      dplyr::filter(Volunteer != volunteer) %>% 
      dplyr::filter(Stim == stim)
    
    o <- capture.output(
      emds[[paste0("Original_", volunteer, "_", stim)]] <- 
        emdEvaluation(file.path(dir,
                                toEvaluate$File),
                      transformList = transformList,
                      manual = lapply(manual, function(x){x$manual}),
                      channels = norm_channels))
  }
}

for (exp in rownames(exp_setups)) {
  for(stim in c("Unstim", "IFNa_LPS")){
    print(paste0("Evaluating ", exp))
    toEvaluate <- samples %>% 
      dplyr::filter(Volunteer != exp_setups[exp, "volunteer"]) %>% 
      dplyr::filter(Stim == stim)
    
    o <- capture.output(
      emds[[paste0(exp, "_", stim)]] <- emdEvaluation(file.path(results, "data", exp, 
                                             paste0("Norm_",
                                                    toEvaluate$File)),
                                   transformList = transformList,
                                   manual = lapply(manual, function(x){x$manual}),
                                   channels = norm_channels))
  }
}

saveRDS(emds,
        file = file.path(results, paste0("emds_diff_methods.Rds")))
```

```{r}
emds <- readRDS(file.path(results,
                         "emds_diff_methods.Rds"))
emds_g <- lapply(names(emds), function(state){
  col_name <- paste0("EMD_",state)
  gather(data.frame(Population = rownames(emds[[state]]), 
                    emds[[state]])[-1,], 
         "Channel", col_name, -"Population")  %>% 
    magrittr::set_colnames(c("Population", "Channel", col_name))
}) %>% magrittr::set_names(names(emds))

emds_g <- cbind(emds_g[[1]][,c(1,2)],
                do.call(cbind, lapply(emds_g, function(x) x[,3])))
rownames(emds_g) <- paste0(emds_g[,1], "_", get_markers(ff, emds_g[,2]))


reduction <- list()
for(exp in colnames(emds_g)[-c(1,2)]){
  exp_original <- gsub(".*_([12].*)", "Original_\\1", exp)
  affected <- (emds_g[,exp] > 2) | (emds_g[ ,exp_original] > 2)
  reduction[[exp]] <- data.frame(population_marker = rownames(emds_g)[affected],
                                 reduction = ((emds_g[,exp_original] - emds_g[,exp]) /
                                                emds_g[,exp_original])[affected],
                                 exp = exp,
                                 method = gsub("_[12].*", "", exp),
                                 subset = gsub("^.*_([12].*)", "\\1", exp))
}
reduction <- do.call(rbind, reduction)

reduction_overall <- cbind("mean" = tapply(reduction$reduction, 
                                         reduction$method, 
                                         mean),
                           "sd" = tapply(reduction$reduction,
                                       reduction$method,
                                       sd))
for (volunteer in c(1,2)) {
  for (stim in c("Unstim", "IFNa_LPS")) {
    for(f in c("mean", "sd")){
      subset <- reduction %>% 
                 dplyr::filter(subset == paste0(volunteer, "_", stim)) 
      reduction_subset <- tapply(subset$reduction, 
                                 subset$method, 
                                 eval(parse(text = f)))
     
      
      reduction_overall <- cbind(reduction_overall,
                                 matrix(reduction_subset,
                                        dimnames = list(names(reduction_subset),
                                                        paste0(volunteer, "_", stim, "_", f))))
        
    }

  }
}
                           
print(reduction_overall)
paste0(round(reduction_overall[, grep("mean",colnames(reduction_overall))], 2),
       " (+- ", 
       round(reduction_overall[, grep("sd",colnames(reduction_overall))], 2),
       ")")
pheatmap::pheatmap(reduction_overall[-1, grep("_mean", colnames(reduction_overall))],
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   number_color = "white")
```


```{r}
for(e in unique(reduction$exp)){
 subset <- reduction %>% 
   dplyr::filter(exp == e)
 plot(density(subset$reduction), main = e)
}
```
```{r}
for (volunteer in c(1,2)) {
  for (stim in c("Unstim", "IFNa_LPS")) {
      subset <- reduction %>% 
                 dplyr::filter(subset == paste0(volunteer, "_", stim)) 
      reduction_subset <- tapply(subset$reduction, 
                                 subset$method, 
                                 eval(parse(text = f)))
     
      
      reduction_overall <- cbind(reduction_overall,
                                 matrix(reduction_subset,
                                        dimnames = list(names(reduction_subset),
                                                        paste0(volunteer, "_", stim, "_", f))))
        
  }
  
}

  
```

```{r}
plot_emds <- function(exp_1, exp_2, emds_g, title = NULL, subtitle = NULL){
  emds_g$affected <- rep("Positive", nrow(emds_g))
  emds_g$affected[emds_g[,exp_1] > emds_g[,exp_2]] <- "Negative"
  emds_g$affected[emds_g[,exp_1] < 2 & emds_g[,exp_2] < 2] <- "No"
  
  if(is.null(title)){
    method <- stringr::str_match(exp_1, "(.*)_([1-2].*)")[,2]
    subset <- stringr::str_match(exp_1, "(.*)_([1-2].*)")[,3]
    title <- "Average change: "
    subtitle <- paste0(round(reduction_overall[method, paste0(subset,"_mean")], 4), 
                       " (±", round(reduction_overall[method, paste0(subset,"_sd")],4), ")")
  }
  
  p <- ggplot(emds_g) +
    geom_point(aes_string(x = exp_1, 
                          y = exp_2, 
                          col = "affected")) +
    geom_rect(data = data.frame(xmin = 0, xmax = 2, ymin = 0, ymax = 2),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#555555", alpha = 0.3) +
    scale_x_continuous(limits = c(0, 15)) +
    scale_y_continuous(limits = c(0, 15)) +
    scale_color_manual(values = c("No" = "grey", 
                                  "Positive" = "black", 
                                  "Negative" = "#b30000")) +
    ggplot2::geom_abline(lwd = 1.5) +
    coord_fixed() +
    ggtitle(title, subtitle) +
    theme_minimal() +
    theme(plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none") 
  emds_g$affected <- NULL
  return(p)
}
```

```{r}

plots <- list()
for (test in c("Unstim", "IFNa_LPS")) {
  for (model in c("Unstim_and_IFNa_LPS", "Unstim", "IFNa_LPS")){
    plots[[paste0(model,"_",test)]] <- plot_emds(paste0("Cluster_QuantileNormLimited_", model, "_1_", test),
                                                 "Original_1_Unstim",
                                                 emds_g) +
      xlab("EMD after normalization") +
      ylab("EMD before normalization") +
      ggtitle(paste0("Trained on ", gsub("_and_", " and ", model)), 
              subtitle = paste0("Validated on ", test)) + 
      theme(plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 11))
  }
}

p <- do.call(grid.arrange, c(plots, list(nrow = 2)))
ggsave(file.path(results,
                 "EMD_Unstim_vs_IFNa_LPS.pdf"),
       p,
       width = 9, height = 6)
p
```

```{r}
plot_titles <- c(c("Cluster_QuantileNormLimited" = "Cluster + 101 quantiles",
                   "QuantileNormLimited" = "101 quantiles",
                   "Cluster_QMinMaxNorm" = "Cluster + 2 quantiles",
                   "QMinMaxNorm" = "2 quantiles"))

plots <- list()
for (volunteer in c(1,2)){
  for (test in c("Unstim", "IFNa_LPS")) {
    for(method in c("Cluster_QuantileNormLimited", "QuantileNormLimited",
                    "Cluster_QMinMaxNorm", "QMinMaxNorm")){
      plots[[paste0(method,"_",volunteer,"_",test)]] <- plot_emds(paste0(method, "_Unstim_and_IFNa_LPS_", volunteer, "_", test),
                                                   paste0("Original_",volunteer,"_", test),
                                                   emds_g) +
        xlab("EMD after normalization") +
        ylab("EMD before normalization") +
        ggtitle(plot_titles[method], 
                subtitle = paste0("Validated on donor ", volunteer, " ", test)) + 
        theme(plot.title = element_text(size = 11),
              plot.subtitle = element_text(size = 11))
    }
  }
}
p <- do.call(grid.arrange, c(plots, list(nrow = 4)))
ggsave(file.path(results,
                 "EMD_methods.pdf"),
       p,
       width = 12, height = 12)
```---
title: "R Notebook"
output: html_notebook
---

# Validation

```{r}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(CytoNorm))
suppressPackageStartupMessages(library(flowWorkspace))
suppressPackageStartupMessages(library(CytoML))
```

```{r}
results <- "results_CytoNorm/validation_pregnancy"
if(!dir.exists(results)) dir.create(results, recursive = TRUE)

norm_dir <- file.path(results, "Normalized")
if(!dir.exists(norm_dir)) dir.create(norm_dir)


files <- c(list.files("D:/Data/Ongoing/Stanford/data/Gates", 
                      pattern = "PTLG.*Unstim.*fcs$", 
                      full.names = TRUE),
           list.files("D:/Data/Ongoing/Stanford/data/Gates_ExtraControls", 
                      pattern = "PTLG.*Unstim.*fcs$", 
                      full.names = TRUE))
files <- grep(paste0(".*", c("PTLG015", "PTLG016", "PTLG006", "PTLG014", "PTLG023"), "_[123BU]",
                     collapse = "|"),
              files, 
              invert = TRUE,
              value = TRUE) 
files <- grep("Unstim_Control_2",
              files,
              invert = TRUE,
              value = TRUE)

df <- data.frame(File = files,
                 Plate = as.factor(str_match(files, "(PTLG[0-9]*)")[, 2]),
                 Timepoint = as.factor(gsub("Unstim", "Control", str_match(files, "PTLG[0-9]*[_Repeat_]*_([^_]*)")[, 2])),
                 stringsAsFactors = FALSE)
rownames(df) <- df$File

repeat_ctrl_file <- "D:/Data/Ongoing/Stanford/data/Gates/Gates_PTLG015_016_Repeats_Unstim_Control.fcs"
df[repeat_ctrl_file, c("Plate","Timepoint")] <- c("PTLG015", "Control") 
head(df)
```

```{r}
train_data <- df %>% 
  dplyr::filter(Timepoint == "Control")

o <- capture.output(ff <- read.FCS(file.path(train_data$File[1])))

norm_markers <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                  11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
                  41, 26, 30, 28, 29, 25, 35)

transformList <- transformList(colnames(ff)[norm_markers], 
                               cytofTransform)
transformList.reverse <- transformList(colnames(ff)[norm_markers], 
                                       cytofTransform.reverse)
```

```{r}
model <- CytoNorm.train(files = train_data$File,
                        labels = train_data$Plate,
                        channels = colnames(ff)[norm_markers],
                        transformList = transformList,
                        outputDir = results,
                        FlowSOM.params = list(nCells = 1000000,
                                              xdim = 15,
                                              ydim = 15,
                                              nClus = 5,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 21,
                                          goal = "PTLG001",
                                          limit = c(0,8)),
                        seed = 1,
                        plot = TRUE,
                        verbose = TRUE)
                        
saveRDS(model,
        file.path(results, "CytoNorm_model_190812.RDS"))   
```

```{r}
fsom <- readRDS(file.path(results,"CytoNorm_FlowSOM.RDS"))

pdf(file.path(results, "cv_test.pdf"))
cv_test <- testCV(fsom, cluster_values = c(5,10,15,20))
dev.off()
saveRDS(cv_test, file.path(results, "cv_test.RDS"))

```

```{r}
validation_data <- df %>% 
  dplyr::filter(Timepoint != "Control")

validation_data$Cohort <- stringr::str_match(validation_data$File, 
                                             "/(Gates[^/]*)/")[,2]

rownames(validation_data) <- gsub(".*/", "", validation_data$File)

validation_data$Plate[validation_data$Plate == "PTLG016"] <- "PTLG015"
```

```{r}
CytoNorm.normalize(model = model,
                   files = validation_data$File,
                   labels = validation_data$Plate,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   outputDir = norm_dir,
                   prefix = "Norm_",
                   verbose = TRUE,
                   normMethod.normalize = QuantileNorm.normalize)
```

```{r}
freq_to_extract <- c("Bcells" = "Bcells", 
                     "CD4+Tcells" = "CD4+Tcells",
                     "CD4+Tcells_mem" = "CD4 mem T cells", 
                     "CD4+Tcells_naive" = "CD4 naive T cells",
                     "CD7+NKcells" = "CD7+NKcells",
                     "CD8+Tcells" = "CD8a+Tcells",
                     "CD8+Tcells_mem" = "CD8 mem T cells",
                     "CD8+Tcells_naive" = "CD8 naive T cells",
                     "TCRgd+Tcells" = "TCRgd+Tcells",
                     "cMCs" = "cMCs",
                     "mDCs" = "mDCs",
                     "pDCs" = "pDCs")

mfis_to_extract <- c("CREB", "ERK", "IkB", "MAPKAPK2", "NFkB", "p38", "S6", 
                     "STAT1", "STAT3", "STAT5")
mfis_feature_names <- apply(expand.grid(names(freq_to_extract), 
                                        mfis_to_extract, 
                                        "Unstim"), 
                            1, paste, collapse = "_")
```

```{r}
freq_file <- file.path(results, "validation_data_frequencies.RDS")

if(! file.exists(freq_file)){

  frequencies <- lapply(list(original = "original", 
                             normalized = "normalized", 
                             gaussNorm = "gaussNorm", 
                             adapted = "adapted"),
                        function(x){matrix(NA,
                                           nrow = nrow(validation_data),
                                           ncol = length(freq_to_extract),
                                           dimnames = list(rownames(validation_data),
                                                           names(freq_to_extract)))})
  mfis <- lapply(list(original = "original", 
                             normalized = "normalized", 
                             gaussNorm = "gaussNorm", 
                             adapted = "adapted"),
                        function(x){matrix(NA,
                                           nrow = nrow(validation_data),
                                           ncol = length(freq_to_extract) * 
                                             length(mfis_to_extract),
                                           dimnames = list(rownames(validation_data),
                                                           mfis_feature_names))})
  
  for(file in rownames(validation_data)[73:112]){
    print(file)
    file_path <- validation_data[file, "File"]
    file_versions <- list("original" = file_path,
                          "normalized" = file.path(norm_dir, 
                                                   paste0("Norm_",  
                                                          gsub(".*/", "", file_path))),
                          "gaussNorm" = file.path("validation_Gates/GaussNorm/Normalized", 
                                                  paste0("Norm_", 
                                                         gsub(".*/", "", file_path))))
    for(version in names(file_versions)){
      fcs_path <- file_versions[[version]]
      o <- capture.output(
        gates <- suppressMessages(
          CytoNorm::applyStaticGating(file.path(results,"Gates_PTLG001_Unstim_Control.wsp"),
                                      fcs_path)))
      
      freqs_tmp <- colSums(gates)
      freqs_tmp <- freqs_tmp / freqs_tmp["CD45+CD66-"]
      frequencies[[version]][file, names(freq_to_extract)] <- freqs_tmp[freq_to_extract]
      
      ff <- transform(read.FCS(fcs_path), transformList)
      for(pop in names(freq_to_extract)){
        selection <- gates[, freq_to_extract[pop]]
        if(sum(selection) > 0){
          mfis_tmp <- apply(ff@exprs[selection, FlowSOM::get_channels(ff, mfis_to_extract), drop = FALSE],
                            2, median)
          mfis[[version]][file, paste(pop, mfis_to_extract, "Unstim", sep = "_")] <- mfis_tmp
        }
      }
    }
  }
  saveRDS(frequencies, freq_file)
  saveRDS(mfis, gsub("frequencies", "mfis", freq_file))
} else {
  frequencies <- readRDS(freq_file)
  mfis <- readRDS(gsub("frequencies", "mfis", freq_file))
}
``` 



```{r}
validation_data <- readRDS(file.path(results, "validation_data.RDS"))
load("C:/Users/sofievg/Downloads/Data_orig.Rda")

rownames(validation_data) <- gsub("_Repeat", "", rownames(validation_data))
frequencies <- lapply(frequencies, function(x){
  rownames(x) <- gsub("_Repeat", "", rownames(x))
  x
})
mfis <- lapply(mfis, function(x){
  rownames(x) <- gsub("_Repeat", "", rownames(x))
  x
})

```

```{r}
validation_data$weeks <- NA
validation_data$adaptedValue <- NA

files_meta_order <- paste0("Gates_", featurepatients, "_", names(featuretimes), "_Unstim.fcs")
validation_data[files_meta_order, "weeks"] <- featureweeks
frequencies[["adapted"]][files_meta_order, names(freq_to_extract)] <- freqfeatures[, names(freq_to_extract)]
mfis[["adapted"]][files_meta_order, mfis_feature_names] <- mfifeatures[, mfis_feature_names]

files_meta_order <- paste0("Gates_", v.featurepatients, "_", names(v.featuretimes), "_Unstim.fcs")
validation_data[files_meta_order, "weeks"] <- v.featureweeks
frequencies[["adapted"]][files_meta_order, names(freq_to_extract)] <- v.freqfeatures[, names(freq_to_extract)]
mfis[["adapted"]][files_meta_order, mfis_feature_names] <- v.mfifeatures[, mfis_feature_names]


validation_data$Timepoint <- factor(validation_data$Timepoint,
                                    levels = c("Control", "BL", "1", "2", "3"))
levels(validation_data$Timepoint) <- c("Control", 
                                       "1st\nTrimester", 
                                       "2nd\nTrimester", 
                                       "3rd\nTrimester", 
                                       "Post\nPartum")

```

```{r fig.width = 9, fig.height = 3}
plot_titles <- c("original" = "Original fcs files",
                 "gaussNorm" = "gaussNorm",
                 "normalized" = "CytoNorm",
                 "adapted" = "Manually adapted gates")

plots <- list()
for(type in c("original", "gaussNorm", "normalized", "adapted")){ # 
  p <- ggplot(data.frame(validation_data,
                         frequencies[[type]],
                         mfis[[type]],
                         check.names = FALSE)) +
    geom_boxplot(aes(x = Timepoint, y = `CD4+Tcells_naive_STAT5_Unstim`)) +
    geom_jitter(aes(x = Timepoint, y = `CD4+Tcells_naive_STAT5_Unstim`, col = Cohort), 
                width = 0.2, height = 0) +
    ggtitle(plot_titles[type]) +
    ylab("STAT5 MFI of CD4+ Naive T cells") +
    ylim(c(0,1)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plots[[length(plots)+1]] <- p
}

p_mfi <- do.call(gridExtra::grid.arrange, 
             c(plots,
               list(nrow = 1)))

ggsave(plot = p_mfi, 
       filename = file.path(results, "validation_CD4+Tcells_naive_STAT5_Unstim.pdf"),
       width = 12, height = 3)
```

```{r}
plots <- list()
for(type in c("original", "gaussNorm", "normalized", "adapted")){ # 
  p <- ggplot(data.frame(validation_data,
                         frequencies[[type]],
                         mfis[[type]],
                         check.names = FALSE)) +
    geom_boxplot(aes(x = Timepoint, y = `CD4+Tcells_naive`)) +
    geom_jitter(aes(x = Timepoint, y = `CD4+Tcells_naive`, col = Cohort), 
                width = 0.2, height = 0) +
    ggtitle(plot_titles[type]) +
    ylab("CD4+ Naive T cells frequency") +
    ylim(c(0,1)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plots[[length(plots)+1]] <- p
}

p <- do.call(gridExtra::grid.arrange, 
             c(plots,
               list(nrow = 1)))

ggsave(plot = p, 
       filename = file.path(results, "validation_CD4+Tcells_naive_Unstim.pdf"),
       width = 10, height = 3)
```




```{r fig.height = 3, fig.width = 12}
plots <- list()
for(type in c("original", "gaussNorm", "normalized", "adapted")){
  p <- ggplot(data.frame(validation_data,
                         frequencies[[type]],
                         mfis[[type]],
                         check.names = FALSE) %>% 
                dplyr::filter(Timepoint %in% c("1st\nTrimester", 
                                               "2nd\nTrimester", 
                                               "3rd\nTrimester"))) +
    geom_point(aes_string(x = "weeks", y = "`CD4+Tcells_naive_STAT5_Unstim`", color = "Plate")) +
    geom_line(aes_string(x = "weeks", y = "`CD4+Tcells_naive_STAT5_Unstim`", 
                         group = "Plate", color = "Plate"))+
    ggtitle(plot_titles[type]) +
    ylim(0,1) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plots[[length(plots)+1]] <- p
}

do.call(gridExtra::grid.arrange, 
        c(plots,
          list(nrow = 1)))
```


```{r fig.width=10}
frequencies_longformat <- 
          lapply(names(frequencies), function(type){
            freq <- frequencies[[type]]
            gather(data.frame(File = rownames(freq), 
                              Cohort = validation_data[rownames(freq), "Cohort"],
                              freq), 
                   key = "Population", 
                   value = !!paste0("Frequency_", type),
                   - "File", - "Cohort")
          })

frequencies_longformat <- frequencies_longformat[[1]] %>% 
  left_join(frequencies_longformat[[2]], by = c("File", "Cohort", "Population")) %>% 
  left_join(frequencies_longformat[[3]], by = c("File", "Cohort", "Population")) %>% 
  left_join(frequencies_longformat[[4]], by = c("File", "Cohort", "Population"))
          
          
plots <- list()
for(type in c("original", "gaussNorm", "normalized")){
  plots[[paste0(type, "_cohort")]] <- ggplot(frequencies_longformat) +
    geom_point(aes_string(x = "Frequency_adapted", 
                          y = paste0("Frequency_", type),
                          col = "Cohort")) +
    ggtitle(plot_titles[type]) +
    xlab("Frequency of manually adapted gates")+
    ylab("Frequency of static gates") +
    theme_minimal() +
    theme(legend.position = "none")
}

p_freq <- do.call(gridExtra::grid.arrange, c(plots, list(nrow = 1)))

ggsave(p_freq, filename = file.path(results,"frequencies.pdf"),
       width = 9, height = 3)

```

```{r}
ggsave(plot = do.call(gridExtra::grid.arrange, list(p_freq, p_mfi)),
       filename = file.path(results, "Figure8.pdf"),
       width = 12, height = 6)
       
```

```{r}
for(type in c("original", "gaussNorm", "normalized")){
  p <- 
    ggplot(frequencies_longformat %>% 
             dplyr::filter(Population == "CD4.Tcells_naive_STAT5_Unstim"))+
    geom_point(aes_string(x = "Frequency_adapted", 
                          y = paste0("Frequency_", type),
                          col = "Population")) +
    theme_minimal() 
  print(p)
}
```