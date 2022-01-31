# tradeSeqPaper
Scripts to reproduce analyses of tradeSeq paper.
---
title: "evaluate K"
author: "Koen Van den Berge"
date: "6/22/2019"
output: html_document
---


```{r}

library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(BiocParallel)
palette(wes_palette("Darjeeling1", 10, type="continuous"))
datasetClusters <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasetClustersSlingshot.txt", header=TRUE)
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/20190611_helper.R")

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

data <- readRDS(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_1.rds"))

counts <- t(data$counts)
falseGenes <- data$tde_overall$feature_id[data$tde_overall$differentially_expressed]
nullGenes <- data$tde_overall$feature_id[!data$tde_overall$differentially_expressed]

# get milestones
gid <- data$prior_information$groups_id
gid <- gid[match(colnames(counts),gid$cell_id),]

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g=12)

# quantile normalization
normCounts <- round(FQnorm(counts))

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:3]
## cluster
nClusters <- datasetClusters$nClusters[1]
set.seed(5)
cl <- kmeans(rd, centers = nClusters)$cluster
 rafalib::mypar(mfrow=c(1,2))
 plot(rd, col = wes_palette("Darjeeling1", 10, type="continuous")[cl], pch=16, asp = 1)
 legend("topleft",legend=as.character(1:nClusters),col=wes_palette("Darjeeling1", 10, type="continuous"),pch=16,cex=2/3,bty='n')
 plot(rd, col = brewer.pal(8,"Dark2")[as.numeric(as.factor(gid$group_id))], pch=16, asp = 1)
 legend("topright",paste0("M",1:length(unique(gid$group_id))), col=1:4, pch=16)

#lineages
lin <- getLineages(rd, cl, start.clus=datasetClusters$start[1], end.clus=c(datasetClusters$end1[1], datasetClusters$end2[1]))
plot(rd, col = pal[g], pch=16, asp = 1)
lines(lin,lwd=2)
#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch=16, asp = 1)
lines(crv, lwd=2, col="black")

### tradeSeq: fit smoothers on truth data
trueWeights <- getWeightsBifurcation(data, crv)
trueT <- matrix(truePseudotime, nrow=length(truePseudotime), ncol=2, byrow=FALSE)

aicMat <- evaluateK(counts=counts, pseudotime=trueT, cellWeights=trueWeights, k=3:10, nGenes=250, ncores=2)

```
---
title: "evaluate methods on second generation simulation"
author: "Koen Van den Berge"
date: "5 February 2019"
output: html_document
---

```{r data}
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)

data <- readRDS("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/bifurcating_4.rds")
counts <- t(data$counts)
falseGenes <- data$feature_info$feature_id[!data$feature_info$housekeeping]
# this dataset has no null genes.
set.seed(5)
null1 <- t(apply(counts,1,sample))
dimnames(null1) <- list(paste0("H",1:501),paste0("C",1:2011))
null2 <- t(apply(counts,1,sample))
dimnames(null2) <- list(paste0("H",502:1002),paste0("C",1:2011))
null3 <- t(apply(counts,1,sample))
dimnames(null3) <- list(paste0("H",1003:1503),paste0("C",1:2011))
counts <- rbind(counts,null1,null2,null3)
nullGenes <- rownames(counts)[substr(rownames(counts),1,1)=="H"]

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous[colnames(counts)]
g <- Hmisc::cut2(truePseudotime,g=12)

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:3]
plot(rd, pch=16, asp = 1)
set.seed(9)
cl <- kmeans(rd, centers = 8)$cluster
plot(rd, col = brewer.pal(9,"Set1")[cl], pch=16, asp = 1)
legend("topleft",legend=as.character(1:7),col=brewer.pal(9,"Set1")[1:7],pch=16,cex=2/3,bty='n')
 plot(rd, col = pal[g], pch=16, asp = 1)
#lineages
lin <- getLineages(rd, cl, start.clus=6, end.clus=c(5,2))
plot(rd, col = pal[g], pch=16, asp = 1)
lines(lin,lwd=2)
#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch=16, asp = 1)
lines(crv, lwd=2, col="black")

# milestone ID
gid <- data$prior_information$groups_id
gid <- gid[match(colnames(counts),gid$cell_id),]
plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1)

```

# fit smoothers on raw data

```{r}
cWeights <- slingCurveWeights(crv)
pseudoT <- slingPseudotime(crv, na=FALSE)
gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=FALSE)
```

# Test for differences at end point

```{r}
endRes <- diffEndTest(gamList)
hist(endRes$pvalue)
deGenesEndGam <- rownames(counts)[which(p.adjust(endRes$pvalue,"fdr")<=0.05)]
mean(falseGenes%in%deGenesEndGam) #TPR
mean(!deGenesEndGam%in%falseGenes) #FDR
```

# edgeR analysis on final clusters

```{r}
clF <- as.factor(cl)
design <- model.matrix(~clF)

library(edgeR)
d <- DGEList(counts)
d <- calcNormFactors(d)
d <- estimateDisp(d, design)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef="clF2") #cluster 2 vs cluster 1 (end clusters)
deGenesEdgeR <- rownames(lrt$table)[p.adjust(lrt$table$PValue,"fdr")<=0.05]
mean(falseGenes%in%deGenesEdgeR) #TPR
mean(deGenesEdgeR%in%nullGenes) #FDR
```

# compare GAM end point vs edgeR analysis

```{r}
library(rafalib)
mypar(mfrow=c(1,1))

## compare
edgeR <- p.adjust(lrt$table$PValue,"fdr")<=0.05
tradeR <- p.adjust(endRes$pval,"fdr")<=0.05
vennC <- cbind(edgeR,tradeR)
vennDiagram(vennC)
```

# Test for different expression pattern

```{r}
resPattern <- patternTest(gamList)
hist(resPattern$pval)
deGenesPattern <- rownames(counts)[which(p.adjust(resPattern$pval,"fdr")<=0.05)]
length(deGenesPattern)
mean(falseGenes%in%deGenesPattern) #TPR
mean(!deGenesPattern%in%falseGenes) #FDR

```

# Monocle BEAM analysis

```{r}
### old monocle BEAM analysis
library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
featureInfo <- data.frame(gene_short_name=rownames(counts))
rownames(featureInfo) <- rownames(counts)
fd <- new("AnnotatedDataFrame", featureInfo)
cds <- newCellDataSet(cellData=counts, featureData=fd, expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds)#, max_components = 4, method = 'ICA')
cds <- orderCells(cds)#, num_paths=2)
plot_cell_trajectory(cds, color_by = "State")
# monocle does not find a branching point on 2 components! specify 4.
BEAM_res <- BEAM(cds,  cores = 1)
sum(BEAM_res$qval<0.05)
```

# tradeR downstream of Monocle 2

```{r}
phenoData(cds)$group_id <- gid$group_id
plot_cell_trajectory(cds, color_by = "group_id")
plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1) ; legend("topleft", paste0("M",1:4), col=2:5,pch=16)
# Milestone 2 and 4 are the two branches to be compared.
ptMon <- matrix(phenoData(cds)$Pseudotime, nrow=ncol(counts), ncol=2, byrow=FALSE)
state <- phenoData(cds)$State
cellWeightsMon <- matrix(0, nrow=ncol(counts), ncol=2)
cellWeightsMon[state==1,] <- c(1/2,1/2)
cellWeightsMon[state==3,1] <- 1
cellWeightsMon[state==2,2] <- 1

gamListMon <- fitGAM(counts, pseudotime=ptMon, cellWeights=cellWeightsMon, verbose=TRUE)
resPatternMon <- patternTest(gamListMon)
resEndMon <- diffEndTest(gamListMon)
```



# tradeR on true pseudotime

```{r}
### tradeR on true pseudotime
pst <- matrix(truePseudotime, nrow=ncol(counts),ncol=2, byrow=FALSE)
gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=slingCurveWeights(crv), verbose=FALSE)
# end point Test
waldEndPointResTrueTime <- diffEndTest(gamListTrueTime)
padjWaldTrueTime <- p.adjust(waldEndPointResTrueTime$pvalue,"fdr")
sum(padjWaldTrueTime<=0.05)
# pattern test
patternResTrueTime <- patternTest(gamListTrueTime)
padjPatternTrueTime <- p.adjust(patternResTrueTime$pvalue,"fdr")
sum(padjPatternTrueTime<=0.05, na.rm=TRUE)

```

# GPfates

```{r}
# export
logCpm <- edgeR::cpm(counts, prior.count=.125, log=TRUE)
sampleInfo <- data.frame(global_pseudotime=truePseudotime)
rownames(sampleInfo) <- colnames(counts)
write.table(logCpm, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/simDyntoyLogCpm.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(sampleInfo, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/sampleInfoSimDyntoy.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
# run GPfates python script
GPfatesWeights <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesWeights.txt", header=FALSE)
GPfatesBif <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesBifStats.txt", header=FALSE)
colnames(GPfatesBif) <- c("bif_ll", "amb_ll", "shuff_bif_ll", "shuff_amb_ll", "phi0_corr", "D", "shuff_D")
```

# GPfates weights + tradeR

```{r}
# based on true pseudotime
gamListGPfatesTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=GPfatesWeights, verbose=FALSE)
# end point test
waldEndPointResTrueTimeGPfates <- diffEndTest(gamListGPfatesTrueTime)
# pattern test
patternResTrueTimeGPfates <- patternTest(gamListGPfatesTrueTime)
```


# plot false genes: true pseudotime

```{r}
pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/falseGenesTruePseudotime.pdf")
for(i in 1:length(falseGenes)){
   plotSmoothers(gamListTrueTime[[falseGenes[i]]])
}
dev.off()
```

# plot false genes: estimated pseudotime

```{r}
pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/falseGenesSlingshot.pdf")
for(i in 1:length(falseGenes)){
   plotSmoothers(gamList[[falseGenes[i]]])
}
dev.off()
```

# plot example genes

```{r}
mypar(mfrow=c(1,2))
i=1
plotSmoothers(gamListTrueTime[[falseGenes[i]]])
plotSmoothers(gamListTrueTime[[nullGenes[i]]])
```

# Performance plots

```{r}
########
# FDP-TPR
########
library(iCOBRA)
library(scales)
truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
truth[falseGenes,"status"] <- 1
cols <- hue_pal()(10)
names(cols) <- c("BEAM", "GPfates", "edgeR", "tradeR_slingshot_end", "tradeR_slingshot_pattern", "tradeR_GPfates_end", "tradeR_GPfates_pattern", "tradeR_Monocle2_end", "tradeR_Monocle2_pattern")

### estimated pseudotime
pval <- data.frame( tradeR_slingshot_end=endRes$pval,
                    tradeR_slingshot_pattern=resPattern$pval,
                    BEAM=BEAM_res$pval,
                    edgeR=lrt$table$PValue,
                    tradeR_Monocle2_end=resEndMon$pvalue,
                    tradeR_Monocle2_pattern=resPatternMon$pvalue,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth)
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra, colorscheme=cols[match(sort(unique(cobra@roc$method)),names(cols))])
plot_roc(cobraplot)
plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1))
pEst <- plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1), xaxisrange=c(0,0.5))


### true pseudotime
pval <- data.frame( tradeR_slingshot_end=waldEndPointResTrueTime$pval,
                    tradeR_slingshot_pattern=patternResTrueTime$pval,
                    BEAM=BEAM_res$pval,
                    edgeR=lrt$table$PValue,
                    tradeR_GPfates_end=waldEndPointResTrueTimeGPfates$pval,
                    tradeR_GPfates_pattern=patternResTrueTimeGPfates$pval,
                      row.names=rownames(counts))
score <- data.frame(GPfates=GPfatesBif$D,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth, score=score)
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra, colorscheme=cols[match(sort(unique(cobra@roc$method)),names(cols))])
plot_roc(cobraplot)
plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1), xaxisrange=c(0,0.5))
pTrue <- plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1), xaxisrange=c(0,0.5))
#
library(cowplot)
prow <- plot_grid( pEst + theme(legend.position="none") + xlab("FDP") + ggtitle("Based on estimated pseudotime"),
           pTrue + theme(legend.position="none") + xlab("FDP") + ggtitle("Based on true pseudotime"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_a <- get_legend(pTrue + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_a, ncol = 1, rel_heights = c(1, .2))
p

# png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2sim2_dyntoy_bifurcating_4/performance_sim2Dyntoy.png", width=7,height=6, units="in", res=300)
# p
# dev.off()

#composite one plot
pval <- data.frame( tradeR_slingshot_end=endRes$pval,
                    tradeR_slingshot_pattern=resPattern$pval,
                    BEAM=BEAM_res$pval,
                    edgeR=lrt$table$PValue,
                    tradeR_GPfates_end=waldEndPointResTrueTimeGPfates$pval,
                    tradeR_GPfates_pattern=patternResTrueTimeGPfates$pval,
                    tradeR_Monocle2_end=resEndMon$pvalue,
                    tradeR_Monocle2_pattern=resPatternMon$pvalue,
                      row.names=rownames(counts))
score <- data.frame(GPfates=GPfatesBif$D,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth, score=score)
saveRDS(cobra,"~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/cobra.rds")
cobra <- calculate_adjp(cobra)
 cobra <- calculate_performance(cobra, binary_truth="status")
 cobraplot <- prepare_data_for_plot(cobra, colorscheme=cols[match(sort(unique(cobra@roc$method)),names(cols))])
  plot_roc(cobraplot)
 plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.8,1), xaxisrange=c(0,0.5))
# pAll <- plot_fdrtprcurve(cobraplot, pointsize=3/2, yaxisrange=c(0.8,1), xaxisrange=c(0,0.5))
# png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/performanceSim1All.png", width=7,height=6, units="in", res=300)
# pAll
# dev.off()
```
---
title: "evaluate K"
author: "Koen Van den Berge"
date: "6/22/2019"
output: html_document
---


```{r}
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(BiocParallel)
palette(wes_palette("Darjeeling1", 10, type="continuous"))


  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

dataAll <- readRDS("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/datasets_for_koen.rds")

data <- dataAll[[1]]
counts <- as.matrix(t(data$counts))
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g=12)
pal <- wes_palette("Zissou1", 12, type = "continuous")

# quantile normalization
normCounts <- round(FQnorm(counts))

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch=16, asp = 1, col=pal[g])

library(princurve)
pcc <- principal_curve(rd, smoother="periodic_lowess")
lines(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], col="red", lwd=2)

### tradeSeq
trueT <- matrix(truePseudotime[colnames(counts)], ncol=1)
cWeights <- rep(1,ncol(counts))
infMat <- evaluateK(counts=counts, pseudotime=trueT, cellWeights=cWeights, k=3:10, nGenes=250, ncores=2)

```
---
title: "evaluate methods on second generation simulation: cyclic trajectory"
author: "Koen Van den Berge"
date: "6 February 2019"
output: html_document
---

```{r data}
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)

data=readRDS("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/72.rds")
counts <- t(data$counts)
falseGenes <- data$feature_info$gene_id[!data$feature_info$housekeeping]
nullGenes <- data$feature_info$gene_id[data$feature_info$housekeeping]

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g=12)

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch=16, asp = 1, col=pal[g])

library(princurve)
pcc <- principal_curve(rd, smoother="periodic_lowess")
lines(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], col="red", lwd=2)
```

# fit smoothers on raw data

```{r}
cWeights <- rep(1,ncol(counts))
pseudoT <- matrix(pcc$lambda,nrow=ncol(counts),ncol=1)
gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=TRUE)
```

# Test for association of expression with the trajectory

```{r}
assocTestRes <- associationTest(gamList)
hist(assocTestRes$pvalue)
```

# tradeR on true pseudotime

```{r}
### tradeR on true pseudotime
pst <- matrix(truePseudotime, nrow=ncol(counts), ncol=1, byrow=FALSE)
gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights, verbose=TRUE)
assocTestTrueRes <- associationTest(gamListTrueTime)
hist(assocTestTrueRes$pvalue)
```

# Monocle 3

```{r}
library(monocle)
fd <- data.frame(gene_short_name=rownames(counts))
fd <- new("AnnotatedDataFrame",fd)
rownames(fd) <- rownames(counts)
pd <- data.frame(cellid=colnames(counts))
pd <- new("AnnotatedDataFrame",pd)
rownames(pd) <- colnames(counts)
cds <- newCellDataSet(counts, featureData=fd, phenoData=pd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- preprocessCDS(cds, num_dim = 20)
cds <- reduceDimension(cds, reduction_method = 'UMAP')
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
pr_graph_test <- principalGraphTest(cds, k=3, cores=1)
plot_cell_trajectory(cds,color_by="cellid") # fails to discover cyclic pattern.
```

# mgcv

```{r}
gamP <- getSmootherPvalues(gamList)
```


# plot false genes

```{r}
pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/falseGenesTruePseudotime.pdf")
for(i in 1:length(falseGenes)){
   plotSmoothers(gamListTrueTime[[falseGenes[i]]])
}
dev.off()
```

# plot null genes

```{r}
pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/nullGenesTruePseudotime.pdf")
for(i in 1:length(nullGenes)){
   plotSmoothers(gamListTrueTime[[nullGenes[i]]])
}
dev.off()
```

# plot example genes

```{r}
mypar(mfrow=c(1,2))
i=1
plotSmoothers(gamListTrueTime[[falseGenes[i]]], main="false gene")
plotSmoothers(gamListTrueTime[[nullGenes[i]]], main="null gene")
```

# Performance plots

```{r}
########
# FDP-TPR
########
library(iCOBRA)
library(scales)
truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
truth[falseGenes,"status"] <- 1

### estimated pseudotime
pval <- data.frame( tradeR_slingshot_assoc=assocTestRes$pval,
                    Monocle3=pr_graph_test$pval,
                    gam=gamP[,1],
                    tradeR_slingshot_assoc_trueTime=assocTestTrueRes$pvalue,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth)
saveRDS(cobra,file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/cobraObject.rds")
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra)
plot_roc(cobraplot)
plot_fdrtprcurve(cobraplot, pointsize=3/2, yaxisrange=c(0,1))
```
---
title: "evaluate methods on second generation simulation: cyclic trajectory"
author: "Koen Van den Berge"
date: "6 February 2019"
output: html_document
---

```{r data}
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)

data=readRDS("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/72.rds")
counts <- t(data$counts)
falseGenes <- data$feature_info$gene_id[!data$feature_info$housekeeping]
nullGenes <- data$feature_info$gene_id[data$feature_info$housekeeping]

# we will generate additional null genes by permuting cells of false genes.
set.seed(9)
permNull1 <- t(apply(counts[falseGenes,],1,sample))
dimnames(permNull1) <- list(paste0("pn",1:length(falseGenes)),colnames(counts))
permNull2 <- t(apply(counts[falseGenes,],1,sample))
dimnames(permNull2) <- list(paste0("pn",(length(falseGenes)+1):(2*length(falseGenes))),colnames(counts))
permNull3 <- t(apply(counts[falseGenes,],1,sample))
dimnames(permNull3) <- list(paste0("pn",(2*length(falseGenes)+1):(3*length(falseGenes))),colnames(counts))
counts <- rbind(counts, permNull1, permNull2, permNull3)


pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g=12)

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch=16, asp = 1, col=pal[g])

library(princurve)
pcc <- principal_curve(rd, smoother="periodic_lowess")
lines(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], col="red", lwd=2)
```

# fit smoothers on raw data

```{r}
cWeights <- rep(1,ncol(counts))
pseudoT <- matrix(pcc$lambda,nrow=ncol(counts),ncol=1)
gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=FALSE)
```

# Test for association of expression with the trajectory

```{r}
assocTestRes <- associationTest(gamList)
hist(assocTestRes$pvalue)
```

# tradeR on true pseudotime

```{r}
### tradeR on true pseudotime
pst <- matrix(truePseudotime, nrow=ncol(counts), ncol=1, byrow=FALSE)
gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights, verbose=TRUE)
assocTestTrueRes <- associationTest(gamListTrueTime)
hist(assocTestTrueRes$pvalue)
```

# Monocle 3

```{r}
library(monocle)
fd <- data.frame(gene_short_name=rownames(counts))
fd <- new("AnnotatedDataFrame",fd)
rownames(fd) <- rownames(counts)
pd <- data.frame(cellid=colnames(counts))
pd <- new("AnnotatedDataFrame",pd)
rownames(pd) <- colnames(counts)
cds <- newCellDataSet(counts, featureData=fd, phenoData=pd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- preprocessCDS(cds, num_dim = 20)
cds <- reduceDimension(cds, reduction_method = 'UMAP')
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
pr_graph_test <- principalGraphTest(cds, k=3, cores=1)
plot_cell_trajectory(cds, color_by=NULL)
```

# mgcv

```{r}
hlp=getSmootherPvalues(gamList)
```

# Performance plots

```{r}
########
# FDP-TPR
########
library(iCOBRA)
library(scales)
truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
truth[falseGenes,"status"] <- 1

### estimated pseudotime
pval <- data.frame( tradeR_slingshot_assoc=assocTestRes$pval,
                    Monocle3=pr_graph_test$pval,
                    gam=hlp[,1],
                    tradeR_slingshot_assoc_trueTime=assocTestTrueRes$pvalue,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth)
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra)
plot_roc(cobraplot)
plot_fdrtprcurve(cobraplot, pointsize=3/2, yaxisrange=c(0,1))
```

# tradeR with UMAP

```{r}
umapDims <- t(reducedDimS(cds))
plot(umapDims, pch=16, asp = 1, col=pal[g])

library(princurve)
pccUmap <- principal_curve(umapDims, smoother="periodic_lowess")
lines(x=pccUmap$s[order(pccUmap$lambda),1], y=pccUmap$s[order(pccUmap$lambda),2], col="red", lwd=2)
```

# fit smoothers on UMAP trajectory

```{r}
cWeights <- rep(1,ncol(counts))
pseudoT <- matrix(pccUmap$lambda,nrow=ncol(counts),ncol=1)
gamListUmap <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=FALSE)
```

# Test for association of expression with the trajectory

```{r}
assocTestUmap <- associationTest(gamListUmap)
hist(assocTestUmap$pvalue)
```

# Performance plots

```{r}
### estimated pseudotime
pval <- data.frame( tradeR_slingshot_PCA=assocTestRes$pval,
                    Monocle3=pr_graph_test$pval,
                    #gam=hlp[,1],
                    #tradeR_slingshot_assoc_trueTime=assocTestTrueRes$pvalue,
                    tradeR_slingshot_UMAP=assocTestUmap$pvalue,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth)
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra)
#facetted(cobraplot) <- FALSE
plot_roc(cobraplot)
perUMAP <- plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0,1), xaxisrange=c(0,0.5), stripsize=0)
perUMAP + xlab("FDP") + theme(axis.title.x = element_text(size = rel(1.5)),
                              axis.title.y = element_text(size = rel(1.5)),
                              axis.text.x=element_text(size=rel(1.2)),
                              axis.text.y=element_text(size=rel(1.2)))
```

# plot for supplementary

```{r}
library(ggplot2)

# PCA trajectory
ggPCA <- ggplot(as.data.frame(rd), aes(x=PC1, y=PC2))
ggPCA + geom_point(col="gray") + theme_bw() + geom_path(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2])

# UMAP trajectory
umapDims <- data.frame(UMAP1=umapDims[,1], UMAP2=umapDims[,2])
ggUMAP <- ggplot(as.data.frame(umapDims), aes(x=UMAP1, y=UMAP2))
ggUMAP + geom_point(col="gray") + theme_bw() + geom_path(x=pccUmap$s[order(pccUmap$lambda),1], y=pccUmap$s[order(pccUmap$lambda),2])



library(cowplot)
prow <- plot_grid( ggPCA + geom_point(col="gray") + theme_bw() + geom_path(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2]),
           ggUMAP + geom_point(col="gray") + theme_bw() + geom_path(x=pccUmap$s[order(pccUmap$lambda),1], y=pccUmap$s[order(pccUmap$lambda),2]),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1,
           ncol = 2
           )
prow

plot_grid(prow,
          perUMAP + xlab("FDP") + theme(axis.title.x = element_text(size = rel(1.5)),
                                        axis.title.y = element_text(size = rel(1.5)),
                                        axis.text.x=element_text(size=rel(1.2)),
                                        axis.text.y=element_text(size=rel(1.2))),
          nrow=2, ncol=1, rel_heights=c(2/3,1), labels=c("","c"))

ggsave("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyngen_cycle_72/permutedGenes/trajectoryWithPermGenes.pdf")

```
---
title: "evaluate methods on second generation simulation"
author: "Koen Van den Berge"
date: "5 February 2019"
output: html_notebook
---

```{r data}
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(tidyverse)
library(dyno)
library(dyntoy)
library(patchwork)
RNGversion("3.5.0")


# set.seed(12)
# data <- generate_dataset(
#     model = model_multifurcating(),
#     num_cells = 750,
#     num_features = 5000,
#     differentially_expressed_rate = .2
#   )

data <- readRDS("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/data.rds")
counts <- t(data$counts)
falseGenes <- data$tde_overall$feature_id[data$tde_overall$differentially_expressed]
nullGenes <- data$tde_overall$feature_id[!data$tde_overall$differentially_expressed]

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous[colnames(counts)]
g <- Hmisc::cut2(truePseudotime,g=12)

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:3]
plot(rd, pch=16, asp = 1)
set.seed(9)
cl <- kmeans(rd, centers = 8)$cluster
plot(rd, col = brewer.pal(9,"Set1")[cl], pch=16, asp = 1)
legend("topleft",legend=as.character(1:8),col=brewer.pal(9,"Set1")[1:8],pch=16,cex=2/3,bty='n')
 plot(rd, col = pal[g], pch=16, asp = 1)
#lineages
lin <- getLineages(rd, cl, start.clus=4, end.clus=c(1,2,8))
plot(rd, col = pal[g], pch=16, asp = 1)
lines(lin,lwd=2)
#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch=16, asp = 1)
lines(crv, lwd=2, col="black")

# milestone ID
gid <- data$prior_information$groups_id
gid <- gid[match(colnames(counts),gid$cell_id),]
plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1)

```

# Deriving number of knots

```{r}
cWeights <- slingCurveWeights(crv)
pseudoT <- slingPseudotime(crv, na=FALSE)
icMat <- evaluateK(counts=counts,  pseudotime=pseudoT, cellWeights=cWeights,
                   nGenes=500, k=3:10)
```



# fit smoothers on raw data

```{r}
gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, 
                  nknots=3, verbose=FALSE)
```

# Test for differences at end point

```{r}
endRes <- diffEndTest(gamList)
hist(endRes$pvalue)
deGenesEndGam <- rownames(counts)[which(p.adjust(endRes$pvalue,"fdr")<=0.05)]
mean(falseGenes%in%deGenesEndGam) #TPR
mean(!deGenesEndGam%in%falseGenes) #FDR
```

# edgeR analysis on final clusters

```{r}
clF <- as.factor(cl)
design <- model.matrix(~clF)

library(edgeR)
d <- DGEList(counts)
d <- calcNormFactors(d)
d <- estimateDisp(d, design)
fit <- glmFit(d, design)
L <- matrix(0, nrow=ncol(fit$coefficients), ncol=3)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- c("1vs2", "2vs3", "1vs3")
L[c("clF2"),1] <- c(1)
L[c("clF8"),2] <- c(1)
L[c("clF2","clF8"),3] <- c(1,-1)
lrt <- glmLRT(fit, contrast=L) #omnibus test between end clusters (clusters 4, 6 and 7).
deGenesEdgeR <- rownames(lrt$table)[p.adjust(lrt$table$PValue,"fdr")<=0.05]
mean(falseGenes%in%deGenesEdgeR) #TPR
mean(deGenesEdgeR%in%nullGenes) #FDR
```

# compare GAM end point vs edgeR analysis

```{r}
library(rafalib)
mypar(mfrow=c(1,1))

## compare
edgeR <- p.adjust(lrt$table$PValue,"fdr")<=0.05
tradeSeq <- p.adjust(endRes$pval,"fdr")<=0.05
vennC <- cbind(edgeR,tradeSeq)
vennDiagram(vennC)
```



# Test for different expression pattern

```{r}
resPattern <- patternTest(gamList)
hist(resPattern$pval)
deGenesPattern <- rownames(counts)[which(p.adjust(resPattern$pval,"fdr")<=0.05)]
length(deGenesPattern)
mean(falseGenes%in%deGenesPattern) #TPR
mean(!deGenesPattern%in%falseGenes) #FDR
```

# Monocle BEAM analysis

```{r}
### old monocle BEAM analysis
library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
featureInfo <- data.frame(gene_short_name=rownames(counts))
rownames(featureInfo) <- rownames(counts)
fd <- new("AnnotatedDataFrame", featureInfo)
cds <- newCellDataSet(cellData=counts, featureData=fd, expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state=5)
plot_cell_trajectory(cds, color_by = "State")
# monocle does not find a branching point on 2 components! specify 4.
BEAM1 <- BEAM(cds,  cores = 1, branch_point=1)
BEAM2 <- BEAM(cds,  cores = 1, branch_point=2)
sum(BEAM1$qval<0.05) ; sum(BEAM2$qval<0.05)
library(aggregation)
pvalBeamFisher <- apply(cbind(BEAM1$pval, BEAM2$pval),1,fisher)
```

# tradeSeq downstream of Monocle 2

```{r}
phenoData(cds)$group_id <- gid$group_id
plot_cell_trajectory(cds, color_by = "group_id")
cols <- brewer.pal(9,"Set1")
plot(rd, col = cols[as.numeric(as.factor(gid$group_id))], pch=16, asp = 1) ; legend("topleft", paste0("M",1:7), col=cols,pch=16)
plot_cell_trajectory(cds, color_by = "State")
# Milestone 2 and 4 are the two branches to be compared.
ptMon <- matrix(phenoData(cds)$Pseudotime, nrow=ncol(counts), ncol=3, byrow=FALSE)
state <- phenoData(cds)$State
cellWeightsMon <- matrix(0, nrow=ncol(counts), ncol=3)
cellWeightsMon[state %in% c(2,5),] <- rep(1/3,3)
cellWeightsMon[state==4,3] <- 1 #shortest lineage
cellWeightsMon[state==3,1] <- 1 #lineage 2 as the bottom lineage
cellWeightsMon[state==1,2] <- 1 #lineage 3 as the upper lineage


plot(rd, col = cols[phenoData(cds)$State], pch=16, asp = 1) ; legend("topleft", paste0(1:5), col=cols,pch=16)

gamListMon <- fitGAM(counts[1:50,], pseudotime=ptMon, cellWeights=cellWeightsMon,
                     verbose=FALSE, nknots=3)
plotSmoothers(gamListMon[[1]])


gamListMon <- fitGAM(counts, pseudotime=ptMon, cellWeights=cellWeightsMon, 
                     verbose=FALSE, nknots=3)
resPatternMon <- patternTest(gamListMon)
resEndMon <- diffEndTest(gamListMon)
```

# GPfates

```{r}
# # export
# logCpm <- edgeR::cpm(counts, prior.count=.125, log=TRUE)
# sampleInfo <- data.frame(global_pseudotime=truePseudotime)
# rownames(sampleInfo) <- colnames(counts)
# write.table(logCpm, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/simDyntoyLogCpm.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
# write.table(sampleInfo, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/sampleInfoSimDyntoy.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
# # run GPfates python script
# system("python3 /Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/20190326_runGPfates_simDyntoy_multifurcating4.py")
# # import output
# GPfatesWeights <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/GPfatesWeights.txt", header=FALSE)
# GPfatesBif <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/GPfatesBifStats.txt", header=FALSE)
# colnames(GPfatesBif) <- c("bif_ll", "amb_ll", "shuff_bif_ll", "shuff_amb_ll", "phi0_corr", "D", "shuff_D")
```

# Performance plots

```{r}
########
# FDP-TPR
########
library(iCOBRA)
library(scales)
truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
truth[falseGenes,"status"] <- 1
cols <- hue_pal()(10)
names(cols) <- c("BEAM", "GPfates", "edgeR", "tradeSeq_slingshot_end", "tradeSeq_slingshot_pattern", "tradeSeq_GPfates_end", "tradeSeq_GPfates_pattern", "tradeSeq_Monocle2_end", "tradeSeq_Monocle2_pattern")

### estimated pseudotime
pval <- data.frame( tradeSeq_slingshot_end=endRes$pval,
                    tradeSeq_slingshot_pattern=resPattern$pval,
                    BEAM=pvalBeamFisher,
                    edgeR=lrt$table$PValue,
                    tradeSeq_Monocle2_end=resEndMon$pvalue,
                    tradeSeq_Monocle2_pattern=resPatternMon$pvalue,
                      row.names=rownames(counts))
cobra <- COBRAData(pval=pval, truth=truth)
saveRDS(cobra, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/cobra.rds")
cobra <- calculate_adjp(cobra)
cobra <- calculate_performance(cobra, binary_truth="status")
cobraplot <- prepare_data_for_plot(cobra, colorscheme=cols[match(sort(unique(cobra@roc$method)),names(cols))])
plot_roc(cobraplot)
plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1))
pEst <- plot_fdrtprcurve(cobraplot, pointsize=0, yaxisrange=c(0.6,1), xaxisrange=c(0,0.5))

```




# OLD

# edgeR using basis functions

## 3 knots

```{r}
gamListIk <- fitGAM(counts[1:100,], pseudotime=pseudoT, cellWeights=cWeights, 
                  nknots=3, verbose=FALSE)

d <- DGEList(counts)
d <- calcNormFactors(d)
designSmooth <- predict(gamListIk[[1]], type="lpmatrix")[,-1]
d <- estimateDisp(d, designSmooth)
fit <- glmFit(d, designSmooth)
# every lineage ends at a knot so we may take the estimated coefficient as a mean.
plotGeneCount(crv, counts, gene=rownames(counts)[1], models=gamListIk)
# lineage 1 has 3 knots, lineage 2 has 2 knots, lineage 3 has 3 knots.

L <- matrix(0, nrow=ncol(fit$coefficients), ncol=3)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- c("3v1", "2v1", "3v2")
L[ c("s(t3):l3.3", "s(t1):l1.3") ,"3v1"] <- c(1,-1)
L[ c("s(t2):l2.2", "s(t1):l1.3") ,"2v1"] <- c(1,-1)
L[ c("s(t3):l3.3", "s(t2):l2.2") ,"3v2"] <- c(1,-1)
lrtSmooth <- glmLRT(fit, contrast=L) #omnibus test between end clusters (clusters 4, 6 and 7).
deGenesEdgeRSmooth <- rownames(lrtSmooth$table)[p.adjust(lrtSmooth$table$PValue,"fdr")<=0.05]
mean(falseGenes%in%deGenesEdgeRSmooth) #TPR
mean(deGenesEdgeRSmooth%in%nullGenes) #FDR

mean(p.adjust(lrtSmooth$table$PValue,"fdr")<=0.05) # everything is significant?
```

## 6 knots

```{r}
gamListIk6 <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, 
                  nknots=6, verbose=FALSE)
endRes6 <- diffEndTest(gamListIk6)
deGenesEndGam6 <- rownames(counts)[which(p.adjust(endRes6$pvalue,"fdr")<=0.05)]
mean(falseGenes%in%deGenesEndGam6) #TPR
mean(!deGenesEndGam6%in%falseGenes) #FDR
```


```{r}
d <- DGEList(counts)
d <- calcNormFactors(d)
designSmooth <- predict(gamListIk6[[1]], type="lpmatrix")[,-1]
d <- estimateDisp(d, designSmooth)
fit <- glmFit(d, designSmooth)
# every lineage ends at a knot so we may take the estimated coefficient as a mean.
plotGeneCount(crv, counts, gene=rownames(counts)[1], models=gamListIk6)
# lineage 1 has 3 knots, lineage 2 has 2 knots, lineage 3 has 3 knots.
L <- matrix(0, nrow=ncol(fit$coefficients), ncol=3)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- c("3v1", "2v1", "3v2")
L[ c("s(t3):l3.6", "s(t1):l1.6") ,"3v1"] <- c(1,-1)
L[ c("s(t2):l2.5", "s(t1):l1.6") ,"2v1"] <- c(1,-1)
L[ c("s(t3):l3.6", "s(t2):l2.5") ,"3v2"] <- c(1,-1)
lrtSmooth <- glmLRT(fit, contrast=L) #omnibus test between end clusters (clusters 4, 6 and 7).
deGenesEdgeRSmooth <- rownames(lrtSmooth$table)[p.adjust(lrtSmooth$table$PValue,"fdr")<=0.05]
mean(falseGenes%in%deGenesEdgeRSmooth) #TPR
mean(deGenesEdgeRSmooth%in%nullGenes) #FDR

```


## 10 knots

```{r}
gamListIk10 <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, 
                  nknots=10, verbose=FALSE)
endRes10 <- diffEndTest(gamListIk10)
deGenesEndGam10 <- rownames(counts)[which(p.adjust(endRes10$pvalue,"fdr")<=0.05)]
mean(falseGenes%in%deGenesEndGam10) #TPR
mean(!deGenesEndGam10%in%falseGenes) #FDR
```

```{r}
d <- DGEList(counts)
d <- calcNormFactors(d)
designSmooth <- predict(gamListIk10[[1]], type="lpmatrix")[,-1]
d <- estimateDisp(d, designSmooth)
fit <- glmFit(d, designSmooth)
# every lineage ends at a knot so we may take the estimated coefficient as a mean.
plotGeneCount(crv, counts, gene=rownames(counts)[1], models=gamListIk6)
# lineage 1 has 3 knots, lineage 2 has 2 knots, lineage 3 has 3 knots.
L <- matrix(0, nrow=ncol(fit$coefficients), ncol=3)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- c("3v1", "2v1", "3v2")
L[ c("s(t3):l3.6", "s(t1):l1.6") ,"3v1"] <- c(1,-1)
L[ c("s(t2):l2.5", "s(t1):l1.6") ,"2v1"] <- c(1,-1)
L[ c("s(t3):l3.6", "s(t2):l2.5") ,"3v2"] <- c(1,-1)
lrtSmooth10 <- glmLRT(fit, contrast=L) #omnibus test between end clusters (clusters 4, 6 and 7).
deGenesEdgeRSmooth10 <- rownames(lrtSmooth10$table)[p.adjust(lrtSmooth10$table$PValue,"fdr")<=0.05]
mean(falseGenes%in%deGenesEdgeRSmooth10) #TPR
mean(deGenesEdgeRSmooth10%in%nullGenes) #FDR

```

# compare GAM end point vs edgeR analysis and edgeR smooth

```{r}
library(rafalib)
mypar(mfrow=c(1,1))

## compare
edgeR <- p.adjust(lrt$table$PValue,"fdr")<=0.05
edgeRSmooth <- p.adjust(lrtSmooth$table$PValue,"fdr")<=0.05
tradeSeq <- p.adjust(endRes$pval,"fdr")<=0.05
vennC <- cbind(edgeR, edgeRSmooth, tradeSeq)
vennDiagram(vennC)
```


---
title: "clustering comparison"
author: "Koen Van den Berge"
date: "11/8/2019"
output: html_document
---


```{r}
library(monocle) # Load Monocle
RNGversion("3.5.0")
library(rafalib)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(clusterExperiment)
library(cluster)
library(scales)
library(tradeSeq)
gcolpal <- c(brewer.pal(8, "Dark2")[-c(2, 3, 5)],
             brewer.pal(12, "Paired")[c(1, 2, 8, 10, 9)],
             brewer.pal(12, "Set3")[c(7, 8, 12)],
             brewer.pal(8, "Pastel2")[8], brewer.pal(11, "BrBG")[11],
             brewer.pal(11, "PiYG")[1], "cyan", "darkblue", "darkorchid2",
             "brown1", "springgreen1", "deepskyblue4", "darkolivegreen",
             "antiquewhite2")
cell_type_color <- c("Basophils" = "#E088B8",
                     "Dendritic cells" = "#46C7EF",
                     "Eosinophls" = "#EFAD1E",
                     "Erythrocyte" = "#8CB3DF",
                     "Monocytes" = "#53C0AD",
                     "Multipotent progenitors" = "#4EB859",
                     "GMP" = "#D097C4",
                     "Megakaryocytes" = "#ACC436",
                     "Neutrophils" = "#F5918A",
                     'NA' = '#000080')


download.file(
  "https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",
  destfile = "./se_paul.rda")
load("./se_paul.rda")
rd <- reducedDim(se)
set.seed(97)
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
plot(rd, col = brewer.pal(9, "Set1")[cl], xlab = "UMAP1", ylab = "UMAP2")
lines(lin, lwd = 2)
crv <- getCurves(lin)
plot(rd, col = brewer.pal(9, "Set1")[cl], main = "color by cluster", xlab = "UMAP1", ylab = "UMAP2")
lines(crv, lwd = 2)
plot(rd, col = cell_type_color[colData(se)$cell_type2], main = "color by cell type", xlab = "UMAP1", ylab = "UMAP2", pch = 16)
lines(crv, lwd = 2)
counts <- as.matrix(assays(se)$counts)


gamListPaul <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=6)

# pattern test
patternResPaul <- patternTest(gamListPaul)
sum(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05, na.rm = TRUE)
patternGenes <- rownames(counts)[which(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05)]


#### fitted values
#RSEC on PCA
resRSEC <- clusterExpressionPatterns(gamListPaul, nPoints=100, genes=patternGenes)
clRSEC <- primaryCluster(resRSEC$rsec)

#RSEC directly
resRSECNoDR <- clusterExpressionPatterns(gamListPaul, nPoints=100, genes=patternGenes,
                                         reduceMethod="none")

## same number of clusters with PAM.
set.seed(882)
yhatScaled <- resRSEC$yhatScaled
resPam <- cluster::pam(x=yhatScaled, k=length(unique(clRSEC))-1)
clPam <- resPam$clustering

## the silhouette value is biased to k-means / PAM since it's based on the same distance.
## a non-parametric approach would be to use bootstrapping.
pt <- slingPseudotime(crv,na=FALSE)
cw <- slingCurveWeights(crv)

ariRSEC <- c()
ariRSECNoDR <- c()
ariPam <- c()
for(ii in 1:6){
  set.seed(ii)
  ## bootstrap cells
  bootCells <- sample(1:ncol(counts), replace=TRUE)
  ptBoot <- pt[bootCells,]
  cwBoot <- cw[bootCells,]
  countsBoot <- counts[,bootCells]
  
  ## refit tradeSeq
  glBoot <- fitGAM(countsBoot[patternGenes,], pseudotime=ptBoot, cellWeights=cwBoot, 
                   verbose=TRUE, nknots=6)
  
  ## cluster using RSEC
  rsecBoot <- clusterExpressionPatterns(glBoot, nPoints=100, genes=patternGenes)
  
  ## cluster using RSEC without PCA
  rsecBootNoDR <- clusterExpressionPatterns(glBoot, nPoints=100, genes=patternGenes,
                                            reduceMethod="none")
  
  ## cluster using PAM
  yhatBoot <- rsecBoot$yhatScaled
  pamBoot <- cluster::pam(yhatBoot, k=length(unique(clRSEC))-1)
  
  ## calculate ARI with full data clustering of respective clustering method.
  ariRSEC[ii] <- mclust::adjustedRandIndex(primaryCluster(rsecBoot$rsec),
                            primaryCluster(resRSEC$rsec))
  
  ariRSECNoDR[ii] <- mclust::adjustedRandIndex(primaryCluster(rsecBootNoDR$rsec),
                            primaryCluster(resRSECNoDR$rsec))
  
  ariPam[ii] <- mclust::adjustedRandIndex(resPam$clustering,
                            pamBoot$clustering)
}

ariRSEC
ariRSECNoDR
ariPam
```

```{r}
ariRSEC <- c(0.1522267, 0.1252056, 0.2309406, 0.2118156, 0.2353261, 0.2095161)
ariRSECNoDR <- c(0.1745011, 0.1634443, 0.2390643, 0.2544524, 0.1995562, 0.2817630)
ariPam <- c(0.1562204, 0.1651667, 0.1519890, 0.1616606, 0.1524510, 0.1511518)

par(bty='l')
df <- data.frame(ari=c(ariRSEC, ariRSECNoDR, ariPam),
                 method=factor(rep(c("RSEC", "RSEC_noDR", "PAM"), each=6),
                               levels=c("PAM", "RSEC_noDR", "RSEC")))
boxplot(ari~method, data=df,
        xlab="Clustering methods", ylab="Adjusted Rand Index")
stripchart(ari~method, data=df, vertical=TRUE, method="jitter", col=c("bisque", "coral", "darkcyan"), pch=19, add=TRUE, cex=1.4)
```

---
title: "evaluate K Paul et al. data"
author: "Koen Van den Berge"
date: "6/22/2019"
output: html_document
---

```{r}
#setwd("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/case/paul/")

cell_type_color <- c("Basophils" = "#E088B8",
                    "Dendritic cells" = "#46C7EF",
                    "Eosinophls" = "#EFAD1E",
                    "Erythrocyte" = "#8CB3DF",
                    "Monocytes" = "#53C0AD",
                    "Multipotent progenitors" = "#4EB859",
                    "GMP" = "#D097C4",
                    "Megakaryocytes" = "#ACC436",
                    "Neutrophils" = "#F5918A",
                    'NA' = '#000080')

### no internet:  use data stored in tradeSeq package
#data(se,package="tradeSeq")
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)

download.file(
  "https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",
  destfile = "./se_paul.rda")
load("./se_paul.rda")
rd <- reducedDim(se)
set.seed(97)
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
library(slingshot)
lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
plot(rd, col = brewer.pal(9, "Set1")[cl], xlab = "UMAP1", ylab = "UMAP2")
lines(lin, lwd = 2)
crv <- getCurves(lin)
plot(rd, col = brewer.pal(9, "Set1")[cl], main = "color by cluster", xlab = "UMAP1", ylab = "UMAP2")
lines(crv, lwd = 2)
plot(rd, col = cell_type_color[colData(se)$cell_type2], main = "color by cell type", xlab = "UMAP1", ylab = "UMAP2", pch = 16)
lines(crv, lwd = 2)
counts <- as.matrix(assays(se)$counts)


######## look at AIC to get K
rafalib::mypar()
library(tradeSeq)
infMat <- evaluateK(counts=counts, pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), k=3:20, nGenes=250, ncores=2)

# with standard errors
aicAvg <- colMeans(infMat$AIC, na.rm=TRUE)
seAic <- sqrt(matrixStats::colVars(infMat$AIC))/sqrt(nrow(infMat$AIC))
plot(x=3:20, y=aicAvg, type='b', ylim=c(3700, 4200), xlab="Number of knots", ylab="Average AIC")
lines(x=3:20, y=aicAvg + seAic, lty=2)
lines(x=3:20, y=aicAvg - seAic, lty=2)
legend("topright",c("average", "average +/- SE"), lty=c(1,2), bty='n')

# plot BIC and AIC for 16 random genes
rafalib::mypar(mfrow=c(4,4)); for(ii in 1:16) plot(infMat$BIC[ii,], type='b')
rafalib::mypar(mfrow=c(4,4)); for(ii in 1:16) plot(infMat$AIC[ii,], type='b')


## we used only 250 genes. reproducible with different seed?
infMat2 <- evaluateK(counts=counts, pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), k=3:20, nGenes=250, ncores=2, seed=2)

```
---
title: "evaluate K Fletcher et al. data"
author: "Koen Van den Berge"
date: "6/22/2019"
output: html_document
---

```{r}
load("~/PhD_Data/singleCell/fletcher/ShareWithKelly/E4c2b_slingshot_wsforkelly.RData")
library(slingshot)
library(rgl)
library(rafalib)
mypar()
library(RColorBrewer)
library(mgcv)
library(tradeSeq)

origData <- t(pcax$x %*% t(pcax$rotation))
colpal <- cc
rd <- X[, 1:5]
lin <- getLineages(rd, clusterLabels = clus.labels, start.clus = "1", end.clus = "4")
crv <- getCurves(lin)
plot(X[, 1:2], col = colpal[as.factor(clus.labels)], pch = 16)
lines(crv)

suppressPackageStartupMessages(library(SummarizedExperiment))
load("~/Downloads/GSE95601_oeHBCdiff_Cufflinks_eSet.rda")
counts <- assayData(Cufflinks_eSet)$counts_table
counts <- counts[!apply(counts, 1, function(row) any(is.na(row))), ]
counts <- counts[, colnames(counts) %in% colnames(origData)]
keep <- rowSums(edgeR::cpm(counts) > 5) >= 15
countsFiltered <- counts[keep, ]
countsFiltered <- countsFiltered[-grep(rownames(countsFiltered), pattern="ERCC"),]

# batch from Github
library(scone)
library(edgeR)
load("~/p63-HBC-diff/output/clust/oeHBCdiff/oeHBCdiff_scone.Rda")
batch <- droplevels(colData(scone_out)$batch[match(colnames(countsFiltered), rownames(colData(scone_out)))])
run <- droplevels(phenoData(Cufflinks_eSet)[[3]])[colnames(exprs(Cufflinks_eSet)) %in% colnames(origData)]

library(zinbwave)
library(BiocParallel)
library(doParallel)
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())

core <- SummarizedExperiment(countsFiltered, colData = data.frame(clusLabel = clus.labels, batch = batch))
#zinb_c <- zinbFit(core, X = '~ clusLabel + batch', commondispersion = TRUE)

#save(zinb_c,file="~/zinbFletcherTradeSeq.rda")
load("~/zinbFletcherTradeSeq.rda")
weights <- computeObservationalWeights(zinb_c, countsFiltered)

######## look at IC to get K
library(tradeSeq)
k <- 3:30
infMat <- evaluateK(counts=countsFiltered, pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights, k=k, nGenes=1000, ncores=2)

infMat2 <- evaluateK(counts=countsFiltered, pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights, k=k, nGenes=1000, ncores=2,
                     seed=1)


# plot BIC and AIC for 16 random genes
rafalib::mypar(mfrow=c(4,4)); for(ii in 1:16) plot(infMat$BIC[ii,], type='b')
rafalib::mypar(mfrow=c(4,4)); for(ii in 1:16) plot(infMat$AIC[ii,], type='b')


### verify AIC calculation
# hlp <- fitGAM(counts=countsFiltered[1:100,], pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights)
# m <- hlp[[1]]
# summ <- summary(m)
# # how AIC is calculated:
# summ$family$aic
# ll <- summ$sp.criterion #is this the REML?
#
# 2*ll + 2*sum(summ$edf[1]) ; m$aic #does not correspond
#
# getAIC <- function(model){
#   summ <- summary(model)
#   ll <- summ$sp.criterion #REML
#   aic <- 2*ll + 2*summ$edf[1]
#   return(aic)
# }
#
# aicIk <- unlist(lapply(hlp, getAIC))
# aicGam <- unlist(lapply(hlp, function(x) x$aic))
# plot(aicIk, aicGam) #ok
#
# getBIC <- function(model){
#   summ <- summary(model)
#   ll <- summ$sp.criterion #REML
#   n <- nrow(model$model) #sample size
#   bic <- 2*ll + summ$edf[1]*log(n)
#   return(bic)
# }
#
# bicIk <- unlist(lapply(hlp, getBIC))
# plot(aicIk, bicIk)

```
---
title: 'Case study: OE dataset'
author: "Koen Van den Berge"
date: "`r Sys.Date()`"
---

# pca data

```{r loadPCAData}
load("~/PhD_Data/singleCell/fletcher/ShareWithKelly/E4c2b_slingshot_wsforkelly.RData")
library(slingshot)
library(rgl)
library(rafalib)
mypar()
library(RColorBrewer)
library(mgcv)
library(tradeSeq)

origData <- t(pcax$x %*% t(pcax$rotation))
colpal <- cc
rd <- X[, 1:5]
lin <- getLineages(rd, clusterLabels = clus.labels, start.clus = "1", end.clus = "4")
crv <- getCurves(lin)
plot(X[, 1:2], col = colpal[as.factor(clus.labels)], pch = 16)
lines(crv)
# 3D plot
rgl::plot3d(X[, 1:3], t = "p", col = colpal[as.factor(clus.labels)], alpha = 0.3, pch = 19, cex = 2, size = 8, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", aspect = "iso", box = FALSE, axes = FALSE)
# rgl::plot3d(X[,1:3], t='p', col=c("white","black")[hlp+1],alpha=0.3, pch = 19, cex = 2, size=8, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
rgl::axes3d(tick = FALSE)
rgl::par3d(windowRect = c(20, 30, 800, 800))
for (i in seq_along(curves)) {
  rgl::plot3d(crv@curves[[i]]$s[order(crv@curves[[i]]$lambda), 1:3], type = "l", add = TRUE, lwd = 4, col = colpal[which.max(tail(lin@lineages[[i]], 1) == levels(clus.labels))])
}
#rgl::plot3d(20, 90, 0, col = "black", add = TRUE)


# rgl.postscript("~/fletcher3d.pdf", fmt="pdf")
```

In the plot below, trajectory 1 (neuronal trajectory) is the green curve, trajectory 2 is the yellow curve, and trajectory 3 (sustentacular trajectory) is the brown curve.

```{r plot1, echo=TRUE, fig.cap="Fletcher 3D PCA", include=identical(knitr:::pandoc_to(), 'html')}
knitr::include_graphics("/Users/koenvandenberge/fletcher3d.pdf")
#![Fletcher 3D PCA](/Users/koenvandenberge/fletcher3d.pdf)
```

```{r cufflinksCountData}
suppressPackageStartupMessages(library(SummarizedExperiment))
# download object from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95601
load("~/Downloads/GSE95601_oeHBCdiff_Cufflinks_eSet.rda")
counts <- assayData(Cufflinks_eSet)$counts_table
counts <- counts[!apply(counts, 1, function(row) any(is.na(row))), ]
counts <- counts[, colnames(counts) %in% colnames(origData)]
keep <- rowSums(edgeR::cpm(counts) > 5) >= 15
countsFiltered <- counts[keep, ]
countsFiltered <- countsFiltered[-grep(rownames(countsFiltered), pattern="ERCC"),]
```

```{r EDA}
# batch from Github
library(scone)
library(edgeR)
load("~/p63-HBC-diff/output/clust/oeHBCdiff/oeHBCdiff_scone.Rda")
batch <- droplevels(colData(scone_out)$batch[match(colnames(countsFiltered), rownames(colData(scone_out)))])
run <- droplevels(phenoData(Cufflinks_eSet)[[3]])[colnames(exprs(Cufflinks_eSet)) %in% colnames(origData)]
# table(batch,run)
```

# fit smoothers on unaligned data

## get ZI weights from zinbwave

```{r zinbwave}
library(zinbwave)
library(BiocParallel)
library(doParallel)
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())

core <- SummarizedExperiment(countsFiltered, colData = data.frame(clusLabel = clus.labels, batch = batch))
#zinb_c <- zinbFit(core, X = '~ clusLabel + batch', commondispersion = TRUE)
#save(zinb_c,file="~/zinbFletcherTradeSeq.rda")
load("~/zinbFletcherTradeSeq.rda")
weights <- computeObservationalWeights(zinb_c, countsFiltered)

```


```{r fitGAM}
#gamList <- tradeSeq::fitGAM(countsFiltered, U=model.matrix(~-1+batch), pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights, nknots=6)
#save(gamList,file="~/gamListOE_tradeSeq_6k.rda")
load("~/gamListOE_tradeSeq_6k.rda")
## check convergence
converged <- unlist(lapply(gamList, function(m){
  if(class(m)=="try-error"){
    return(FALSE)
  } else return(m$converged)
}))
mean(converged)
sum(converged) ; length(gamList)
```


# BETWEEN LINEAGE COMPARISONS

## Comparing the expression at the differentiated cell type

```{r}
endTestGam <- diffEndTest(gamList, global = TRUE, pairwise = TRUE)
endOmnibusPval <- endTestGam$pvalue
sum(is.na(endOmnibusPval)) # genes we could not fit or test.
endOmnibusPadj <- p.adjust(endOmnibusPval, "fdr")
sum(endOmnibusPadj <= 0.05, na.rm = TRUE)
mean(endOmnibusPadj <= 0.05, na.rm = TRUE)
deGenesEndGam <- which(endOmnibusPadj <= 0.05)
hist(endOmnibusPval)
```

# edgeR analysis

In the edgeR analysis, I compare the gene expression of the differentiated cell types (i.e. final cluster of a trajectory) to each other, between the different trajectories.

```{r}
library(edgeR)
d <- DGEList(countsFiltered)
d <- calcNormFactors(d)
design <- model.matrix(~clus.labels+batch)
d$weights <- weights
d <- estimateDisp(d,design)
fit <- glmFit(d,design)
L <- matrix(0, nrow = ncol(fit$coefficients), ncol = 3)
rownames(L) <- colnames(fit$coefficients)
# trajectory 1 vs. 2
L[c("clus.labels12", "clus.labels15"), 1] <- c(1, -1)
# trajectory 3 vs 1
L[c("clus.labels4", "clus.labels12"), 2] <- c(1, -1)
# trajectory 3 vs 2
L[c("clus.labels4", "clus.labels15"), 3] <- c(1, -1)
lrt <- zinbwave::glmWeightedF(fit, contrast = L)
deGenesEdgeR <- which(p.adjust(lrt$table$PValue, "fdr") <= 0.05)
length(deGenesEdgeR)
mean(deGenesEdgeR %in% deGenesEndGam)
hist(lrt$table$PValue)
```

# compare analyses

We retrieve 85% of the genes that edgeR finds, but also obtain a bunch of other genes.

```{r}
mypar(mfrow=c(1,1))

## compare
edgeRDE <- p.adjust(lrt$table$PValue, "fdr") <= 0.05
gamDE <- endOmnibusPadj <= 0.05
vennC <- cbind(edgeRDE, gamDE)
vennDiagram(vennC, main = "end point comparison across all trajectories")
```

# plot unique GAM genes

```{r}
uniqueGamEndId <- which(gamDE == TRUE & edgeRDE == FALSE)

i <- 0
while (i < 10) {
  i <- i + 1
  plotSmoothers(gamList[[uniqueGamEndId[i]]])
}
```

# plot unique edgeR genes

```{r}
uniqueEdgeREndId <- which(gamDE == FALSE & edgeRDE == TRUE)

i <- 0
while (i < 10) {
  i <- i + 1
  plotSmoothers(gamList[[uniqueEdgeREndId[i]]])
}
```

# are unique GAM genes relevant?

```{r}
uniqGamEndFletcher <- rownames(countsFiltered)[uniqueGamEndId]
write.table(uniqGamEndFletcher, file = "~/uniqGamEndFletcher.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# submit to http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp for top 20 gene sets
overlapUniqEndGam <- readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/case/fletcher/overlapUniqGamEndFletcher_6k")
overlapSets <- overlapUniqEndGam[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x, split = "\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit, "[[", 1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit, "[[", 2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit, "[[", 4)))
pvalUniqGam <- unname(unlist(lapply(overlapSetsSplit, "[[", 6)))
qval <- unname(unlist(lapply(overlapSetsSplit, "[[", 7)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x = gsNames, pattern = "_", replacement = " ")
gsNames[-1] <- unname(sapply(gsNames[-1], function(x) substr(x, 4, nchar(x))))
tabUniqGam <- data.frame(
  geneSet = gsNames[-1],
  overlap = as.numeric(genesInOverlap[-1]),
  genesInSet = as.numeric(genesInSet[-1]),
  qvalue = qval[-1]
)
library(xtable)
xtable(tabUniqGam)

## figure with overlap as scale
pval <- pvalUniqGam
tab <- tabUniqGam
library(scales)
genesInComp <- 1576 #nr of genes assessed by MSigDB
genesInUniverse <- 45956 #according to MSigDB
background <- genesInComp / genesInUniverse
tabIk <- tab[nrow(tab):1,]
setOverlap <- as.numeric(as.character(tabIk$overlap))/as.numeric(as.character(tabIk$genesInSet))

pal <- colorRampPalette(c("red","yellow"))(20)
hlpBar <- barplot(setOverlap, horiz=TRUE, xlim=c(0,0.2), col=alpha(pal[order(as.numeric(as.character(tabIk$qvalue)), decreasing=FALSE)],.2))
barplot(rep(background,20), horiz=TRUE, add=TRUE, col=alpha("grey",.3))
abline(v=background, col="red", lty=2, lwd=2)
text(y=hlpBar, x=setOverlap, label=paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet)), pos=4, cex=1)
text(y=hlpBar, x=rep(0,20), label=as.character(tabIk$geneSet), pos=4, cex=.8, font=2)

## figure with p-value as scale
library(scales)
genesInComp <- 1576 #nr of genes assessed by MSigDB
genesInUniverse <- 45956 #according to MSigDB
background <- genesInComp / genesInUniverse
tabIk <- tab[nrow(tab):1,]
pvalScale <- rev(-log10(as.numeric(pval[-1])))

pal <- colorRampPalette(c("red","yellow"))(20)
hlpBar <- barplot(pvalScale, horiz=TRUE, xlim=c(0,66), col=alpha(pal[order(as.numeric(as.character(tabIk$qvalue)), decreasing=FALSE)],.2),
xlab="-log10 p-value")
#barplot(rep(background,20), horiz=TRUE, add=TRUE, col=alpha("grey",.3))
#abline(v=background, col="red", lty=2, lwd=2)
text(y=hlpBar, x=pvalScale, label=paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet)), pos=4, cex=1)
text(y=hlpBar, x=rep(0,20), label=as.character(tabIk$geneSet), pos=4, cex=1.1, font=2)

```



```{r patternTest}
resPat <- patternTest(gamList, global = TRUE, pairwise = TRUE)
o <- order(resPat$waldStat, decreasing = TRUE)

i <- 0
while (i < 10) {
  i <- i + 1
  plotSmoothers(gamList[[o[i]]], main = i)
}

# for paper
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/case/fletcher/plotSmoothersFletcher.R")
library(scales)
png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/topGenesFletcherPattern_6k.png", width = 6, height = 7, units = "in", res = 300)
mypar(mfrow = c(3, 2))
i <- 0
while (i < 6) {
  i <- i + 1
  plotSmoothersIk(gamList[[o[i]]], main = rownames(resPat)[o][i])
}
dev.off()

# Stage-wise testing
library(stageR)
pScreen <- resPat$pvalue
names(pScreen) <- rownames(resPat)
pConfirmation <- cbind(resPat$pvalue_1vs2, resPat$pvalue_1vs3, resPat$pvalue_2vs3)
dimnames(pConfirmation) <- list(rownames(resPat), c("1v2", "1v3", "2v3"))
stageObj <- stageR(pScreen, pConfirmation, pScreenAdjusted = FALSE)
stageObj <- stageWiseAdjustment(stageObj, alpha = .05, method = "holm", allowNA = TRUE)
res <- getResults(stageObj)
colSums(res)
allCompGenes <- names(which(rowSums(res) == 4))
write.table(allCompGenes, file = "~/allCompGenes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# submit these in http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp

overlap <- readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/overlapCellCycleOE_6k")
overlapSets <- overlap[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x, split = "\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit, "[[", 1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit, "[[", 2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit, "[[", 4)))
qval <- unname(unlist(lapply(overlapSetsSplit, "[[", 7)))
pvalPatternOE <- unname(unlist(lapply(overlapSetsSplit, "[[", 6)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x = gsNames, pattern = "_", replacement = " ")
gsNames[-1] <- unname(sapply(gsNames[-1], function(x) substr(x, 4, nchar(x))))
tabPatternOE <- data.frame(
  geneSet = gsNames[-1],
  overlap = genesInOverlap[-1],
  genesInSet = genesInSet[-1],
  qvalue = qval[-1]
)
 library(xtable)
 xtable(tabPatternOE)

# cluster fitted values of genes significant for pattern test and create a heatmap.
# or plot gene clusters their smoothing profiles
log(fitted(gamList[[1]]) + .1)
```


# within lineage comparisons

### cell cycle patterns within lineages: tradeSeq associationTest

```{r}
assocRes <- associationTest(gamList, global = TRUE, lineages = TRUE)
deAssoc1 <- rownames(assocRes)[p.adjust(assocRes$pvalue_1, "fdr") <= 0.05]
deAssoc1 <- deAssoc1[!is.na(deAssoc1)]
length(deAssoc1)

o = order(assocRes$waldStat_1, decreasing=TRUE)
write.table(head(assocRes[o,"waldStat_1",drop=FALSE], 2000), file="~/topGenesNeuronalLineage.txt", quote=FALSE, col.names=FALSE)
```

#### cell cycle gene set

```{r}
cellGenes <- read.table("~/Downloads/GO_term_summary_20181001_091947.txt", header = TRUE, sep = "\t", row.names = NULL)
cellGenes$MGI.Gene.Marker.ID <- as.character(cellGenes$MGI.Gene.Marker.ID)
cellGenes <- cellGenes[cellGenes$MGI.Gene.Marker.ID %in% rownames(countsFiltered), ]
# out of a total of 2889 genes present in filtered dataset
sum(cellGenes$MGI.Gene.Marker.ID %in% deAssoc1)

# make enrichment plot
# how many cell cycle genes were discovered along ordering of significance?
o <- order(assocRes$waldStat_1, decreasing = TRUE)
sumFoundCC <- sapply(1:nrow(assocRes), function(ii) {
  sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(assocRes)[o[1:ii]], na.rm = TRUE)
})
# how many do we expect by chance?
sumFoundRandom <- (nrow(cellGenes) / nrow(assocRes)) * (1:nrow(assocRes))
# absolute
plot(x = 1:nrow(assocRes), y = sumFoundCC, type = "l")
lines(x = 1:nrow(assocRes), y = (2889 / nrow(assocRes)) * (1:nrow(assocRes)), type = "l", col = "steelblue")
expRandom <- (2889 / nrow(smootherStats)) * (1:nrow(smootherStats))
# relative
plot(x = 1:nrow(assocRes), y = sumFoundCC / sumFoundRandom, type = "l", ylab = "# genes found with tradeSeq / # genes found under random selection")
abline(h = 1, lty = 2, col = "red")


png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/cellCycleFletcher.png", width = 7, height = 7, units = "in", res = 200)
plot(x = 1:nrow(assocRes), y = sumFoundCC / sumFoundRandom, type = "l", ylab = "# genes found with tradeSeq / # genes found under random selection", bty = "n", xlab = "Gene list ordered according to significance", col = "darkgray", lwd = 1.5)
abline(h = 1, lty = 2, col = "red")
dev.off()
```


#### Heatmap of top genes

```{r}
o1 <- order(assocRes[deAssoc1, "waldStat_1"], decreasing = TRUE)
top200Assoc <- deAssoc1[o1[1:200]]

# heatmap for neuronal lineage.
source("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeq/R/utils.R")
df <- .getPredictRangeDf(gamList[[1]], lineageId = 1, nPoints = 100)
y <- do.call(rbind, lapply(gamList[top200Assoc], predict, newdata = df, type = "link"))
yScaled <- t(scale(t(y)))
pst <- slingPseudotime(crv)
pst <- cbind(pst, clus.labels)
pst <- pst[!is.na(pst[, 1]), ]
pst <- pst[order(pst[, 1], decreasing = FALSE), ]

df$ct <- NA
df$cols <- NA
df$ct[df$t1 < 80] <- "HBC"
df$cols[df$t1 < 80] <- colpal[1]
# df$ct[df$t1>=80 & df$t1<115] <- expression(paste(Delta, "HBC1"))
df$ct[df$t1 >= 80 & df$t1 < 115] <- "HBC1"
df$cols[df$t1 >= 80 & df$t1 < 115] <- colpal[3]
# df$ct[df$t1>=90 & df$t1<100] <- expression(paste(Delta, "HBC2"))
df$ct[df$t1 >= 90 & df$t1 < 100] <- "HBC2"
df$cols[df$t1 >= 90 & df$t1 < 100] <- colpal[6]
df$ct[df$t1 >= 115 & df$t1 < 220] <- "GBC"
df$cols[df$t1 >= 115 & df$t1 < 220] <- colpal[17]
df$ct[df$t1 >= 220 & df$t1 < 250] <- "INP1"
df$cols[df$t1 >= 220 & df$t1 < 250] <- colpal[24]
df$ct[df$t1 >= 250 & df$t1 < 310] <- "INP2"
df$cols[df$t1 >= 250 & df$t1 < 310] <- colpal[4]
df$ct[df$t1 >= 310 & df$t1 < 325] <- "INP3"
df$cols[df$t1 >= 310 & df$t1 < 325] <- colpal[19]
df$ct[df$t1 >= 325 & df$t1 < 385] <- "iOSN"
df$cols[df$t1 >= 325 & df$t1 < 385] <- colpal[13]
df$ct[df$t1 >= 385] <- "mOSN"
df$cols[df$t1 >= 385] <- colpal[8]

library(pheatmap)
annoCol <- data.frame(celltype = df$ct)
rownames(annoCol) <- colnames(yScaled)
dfUniq <- unique(df[, c("ct", "cols")])
colhlp <- dfUniq$cols
names(colhlp) <- dfUniq$ct
annColors <- list(celltype = colhlp)
heatRaw <- pheatmap(yScaled, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, cutree_rows = 5, clustering.method = "ward.D", annotation_col = annoCol, annotation_colors = annColors, annotation_names_col = FALSE)

## reorder big clusters.
library(vegan)
set.seed(7)
origHClust <- heatRaw$tree_row
hlpClust <- function(hc, mat) {
  reorder(hc, wts = runif(200))
}
pheatmap(yScaled, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, cutree_rows = 5, clustering.method = "ward.D", annotation_col = annoCol, annotation_colors = annColors, annotation_names_col = FALSE, clustering_callback = hlpClust)
```

# Check for HBC markers

```{r}
startRes <- startVsEndTest(gamList)
top250Genes <- rownames(head(startRes[order(startRes$waldStat, decreasing = TRUE), ], 250))
write.table(top250Genes, file = "~/startTestOETop250_tradeSeq.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# submit to http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp for top 20 gene sets


overlap <- readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/overlapStartTestOE_tradeSeq_6k")
overlapSets <- overlap[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x, split = "\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit, "[[", 1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit, "[[", 2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit, "[[", 4)))
pvalStart <- unname(unlist(lapply(overlapSetsSplit, "[[", 6)))
qval <- unname(unlist(lapply(overlapSetsSplit, "[[", 7)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x = gsNames, pattern = "_", replacement = " ")
gsNames[-1] <- unname(sapply(gsNames[-1], function(x) substr(x, 4, nchar(x))))
tabStart <- data.frame(
  geneSet = gsNames[-1],
  overlap = genesInOverlap[-1],
  genesInSet = genesInSet[-1],
  qvalue = qval[-1]
)
# library(xtable)
# xtable(tabStart)
```

# figure of gene set enrichment

```{r}
library(scales)
rafalib::mypar(mar=c(2.5,4,1.5,1),mfrow=c(1,3))
pal <- colorRampPalette(c("red","yellow"))(20)

## startVsEndTest: top 250 genes
tabIk <- tabStart[nrow(tabStart):1,]
pvalScale <- rev(-log10(as.numeric(pvalStart[-1])))
hlpBar <- barplot(pvalScale, horiz=TRUE, xlim=c(0,22), col=alpha(pal[order(as.numeric(as.character(tabIk$qvalue)), decreasing=FALSE)],.2),
xlab="-log10 p-value")
#text(y=hlpBar, x=pvalScale, label=paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet)), pos=4, cex=1)
txts <- paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet))
for(ii in 1:20) mtext(txts[ii], side=2, at=hlpBar[ii,], las=1, cex=2/3, line=.2)
text(y=hlpBar, x=rep(0,20), label=as.character(tabIk$geneSet), pos=4, cex=1.1, font=2)
mtext("a", side=3, at=-2.5, font=2, cex=3/2)

##  patternTest: genes significant in all three comparisons
tabIk <- tabPatternOE[nrow(tabPatternOE):1,]
pvalScale <- rev(-log10(as.numeric(pvalPatternOE[-1])))
hlpBar <- barplot(pvalScale, horiz=TRUE, xlim=c(0,30), col=alpha(pal[order(as.numeric(as.character(tabIk$qvalue)), decreasing=FALSE)],.2),
xlab="-log10 p-value")
#text(y=hlpBar, x=pvalScale, label=paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet)), pos=4, cex=1)
txts <- paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet))
for(ii in 1:20) mtext(txts[ii], side=2, at=hlpBar[ii,], las=1, cex=2/3, line=.2)
text(y=hlpBar, x=rep(0,20), label=as.character(tabIk$geneSet), pos=4, cex=1.1, font=2)
mtext("b", side=3, at=-3, font=2, cex=3/2)

##  genes unique in tradeSeq diffEnd (vs edgeR)
tabIk <- tabUniqGam[nrow(tabUniqGam):1,]
pvalScale <- rev(-log10(as.numeric(pvalUniqGam[-1])))
hlpBar <- barplot(pvalScale, horiz=TRUE, xlim=c(0,66), col=alpha(pal[order(as.numeric(as.character(tabIk$qvalue)), decreasing=FALSE)],.2),
xlab="-log10 p-value")
#text(y=hlpBar, x=pvalScale, label=paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet)), pos=4, cex=1)
txts <- paste0(as.character(tabIk$overlap),"/",as.character(tabIk$genesInSet))
for(ii in 1:20) mtext(txts[ii], side=2, at=hlpBar[ii,], las=1, cex=2/3, line=.2)
text(y=hlpBar, x=rep(0,20), label=as.character(tabIk$geneSet), pos=4, cex=1.1, font=2)
mtext("c", side=3, at=-9, font=2, cex=3/2)

```



#### early DE test


```{r}
early13 <- earlyDETest(gamList, knots = c(1, 3), nPoints = 50, global = TRUE, pairwise = TRUE)
library(stageR)
pScreen <- early13$pvalue
names(pScreen) <- rownames(countsFiltered)
pConfirmation <- cbind(early13$pvalue_1vs2, early13$pvalue_1vs3, early13$pvalue_2vs3)
rownames(pConfirmation) <- rownames(countsFiltered)
colnames(pConfirmation) <- c("1v2", "1v3", "2v3")
stageRObj <- stageR(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted = FALSE)
stageRObj <- stageWiseAdjustment(stageRObj, method = "holm", alpha = 0.05, allowNA = TRUE)
res <- getResults(stageRObj)
sigAll13 <- names(which(rowSums(res) == 4))
# most of these genes seem to be very relevant.
oo <- order(early13$waldStat, decreasing=TRUE)
head(early13[oo,])

library(UpSetR)
resDf <- as.data.frame(res)
colnames(resDf) <- c("Global", "Neur. vs. Microv.", "Neur. vs. Sust.", "Microv. vs. Sust.")
upset(resDf[, -1], order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1))

plotSmoothersIk2 <- function(m, nPoints = 100, legendPos="topright", ...) {
  data <- m$model
  y <- data$y
  # construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  # cols <- c("#E7298A", "#FF7F00", "#1F78B4")
  cols <- c("#FF7F00", "#1F78B4", "#E7298A")
  plot(
    x = timeAll, y = log(y + 1), col = alpha(cols[col], 2 / 3), pch = 16, cex = 2 / 3,
    ylab = "log(count + 1)", xlab = "Pseudotime", ...
  )

  # predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
    # average batch
    df$U[] <- 0
    df$U[,"batchGBC09B"] <- 1
    yhat <- predict(m, newdata = df, type = "response")
    lines(x = df[, paste0("t", jj)], y = log(yhat+1), col = cols[jj], lwd = 2)
  }
  # knots
  abline(v = gamList[[1]]$smooth[[1]]$xp[1], lty = 2, col = "black", lwd = 1.5)
  abline(v = gamList[[1]]$smooth[[1]]$xp[3], lty = 2, col = "black", lwd = 1.5)
  legend(legendPos, c("Neuronal", "Microvillous", "Sustentacular"),
    col = cols,
    lty = 1, lwd = 2, bty = "n", cex = 4 / 5
  )
}

plotSmoothersIk2(gamList[["Sox11"]])


control=gam.control()
control$maxit=1000
hlp = tradeSeq::fitGAM(countsFiltered[9600:9700,], U=model.matrix(~-1+batch, contrasts.arg=list(batch=contr.sum)), pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights[1300:1400,], nknots=6, control=control)

plotSmoothersIk2(hlp[["Sox11"]])


### look at TFs
sigAll <- rownames(res)[res[,1]==1]
# download list of TFs from http://www.tfcheckpoint.org/index.php/browse
tfCheck <- read.table("~/Downloads/export.txt", sep="\t", header=TRUE)
tfAll <- as.character(tfCheck[,1])
sigAllTFs = sigAll[toupper(sigAll) %in% tfAll]
write.table(sigAllTFs, file="~/sigAllTFs.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# submit for GSEA
overlapTF <- readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/case/fletcher/overlapEarlyDETFs")
overlapSets <- overlapTF[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x, split = "\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit, "[[", 1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit, "[[", 2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit, "[[", 4)))
pvalTF <- unname(unlist(lapply(overlapSetsSplit, "[[", 6)))
pvalTF <- as.numeric(pvalTF[-1])
qval <- unname(unlist(lapply(overlapSetsSplit, "[[", 7)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x = gsNames, pattern = "_", replacement = " ")
gsNames[-1] <- unname(sapply(gsNames[-1], function(x) substr(x, 4, nchar(x))))
tabTF <- data.frame(
  geneSet = gsNames[-1],
  overlap = as.numeric(genesInOverlap[-1]),
  genesInSet = as.numeric(genesInSet[-1]),
  qvalue = qval[-1]
)
library(xtable)
xtable(tabTF)


plot(x=tabTF$overlap/tabTF$genesInSet, y=-log10(pvalTF), pch=16)

## plot epithelial TFs
# download epithelial cell differentiation gene set at http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=GO_EPITHELIAL_CELL_DIFFERENTIATION
epiGenes <- as.character(read.table("~/Downloads/geneset.txt", skip=2)[,1])
epiSigTFs <- sigAllTFs[toupper(sigAllTFs) %in% epiGenes]
mypar(mfrow=c(2,2))
for(ii in 1:length(epiSigTFs)) plotSmoothersIk2(gamList[[epiSigTFs[ii]]])
earlyEpiTF <- early13[epiSigTFs,]
oo13 <- order(earlyEpiTF[,"waldStat"], decreasing=TRUE)

pos=c("topleft", "topright", "topright", "topleft", rep("topright",4))
mypar(mfrow=c(2,2), bty='l')
for(ii in c(1,3:5)) plotSmoothersIk2(gamList[[rownames(earlyEpiTF[oo13,]) [ii] ]], legendPos=pos[ii], ylim=c(0,10))


###### ADDITIONAL RESULTS
### cross-reference with TFs from Fletcher et al.
# these are TFs that are DE between clusters within a lineage
deTFs <- openxlsx::read.xlsx("~/Downloads/mmc4.xlsx")[,1] # supplementary table 3 from fletcher et al. paper
# which ones are also early DE between lineages?
earlyTFs <- deTFs[deTFs %in% sigAll]
# plot earlyTFs
pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/earlyTFs_OE.pdf", width = 9, height = 7)
rafalib::mypar(mfrow = c(4, 4), bty = "l")
for(ii in 1:length(earlyTFs)) plotSmoothersIk2(gamList[[earlyTFs[ii]]], main=earlyTFs[ii])
dev.off()

# plot significant genes
library(scales)
png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/earlyDE13SignGenesInAll.png", width = 9, height = 7, units = "in", res = 200)
rafalib::mypar(mfrow = c(3, 3), bty = "l")
for (ii in 1:9) {
  plotSmoothersIk2(gamList[[sigAll13[ii]]], ylim = c(0, 10), main = sigAll13[ii])
}
dev.off()

pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/earlyDE13SignGenesInAll_allGenes.pdf", width = 9, height = 7)
rafalib::mypar(mfrow = c(4, 4), bty = "l")
for (ii in 1:length(sigAll13)) {
  plotSmoothersIk2(gamList[[sigAll13[ii]]], ylim = c(0, 10), main = sigAll13[ii])
}
dev.off()

### custom background GSEA
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
bm <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), # this is what we want to extract
            filters="mgi_symbol", # this determines the filter
            values=tfAll, # this are the values to filter for (get only the genes in our list)
            mart=ensembl)
ensemblTF <- bm[,1]
# for DAVID
write.table(ensemblTF, file="~/tfAllEnsembl.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


```

# plot for paper v2

```{r}
plotSmoothersIk <- function(m, nPoints = 100, ...) {
  data <- m$model
  y <- data$y
  # construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  # cols <- c("#E7298A", "#FF7F00", "#1F78B4")
  cols <- c("#FF7F00", "#1F78B4", "#E7298A")
  plot(
    x = timeAll, y = log(y + 1), col = alpha(cols[col], 2 / 3), pch = 16, cex = 2 / 3,
    ylab = "log(count + 1)", xlab = "Pseudotime", ...
  )

  # predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
    yhat <- predict(m, newdata = df, type = "response")
    lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = cols[jj], lwd = 2)
  }
  # knots
  abline(v = gamList[[1]]$smooth[[1]]$xp[2], lty = 2, col = "black", lwd = 1.5)
  abline(v = gamList[[1]]$smooth[[1]]$xp[4], lty = 2, col = "black", lwd = 1.5)
  legend("topleft", c("Neuronal", "Microvillous", "Sustentacular"),
    col = cols,
    lty = 1, lwd = 2, bty = "n", cex = 4 / 5
  )
}

## heatmap
png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/FigCaseFletcher_heatmap.png", width = 7, height = 6, units = "in", res = 200)
## reorder big clusters.
set.seed(7)
origHClust <- heatRaw$tree_row
hlpClust <- function(hc, mat) {
  stats::reorder(origHClust, wts = runif(200))
}
pheatmap(yScaled, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, cutree_rows = 5, clustering.method = "ward.D", annotation_col = annoCol, annotation_colors = annColors, annotation_names_col = FALSE, clustering_callback = hlpClust)
dev.off()

### early DE genes
library(scales)
png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/FigCaseFletcher_earlyDE_TFs.png", width = 7, height = 6, units = "in", res = 200)
rafalib::mypar(mfrow = c(2, 2), bty='l')
for(ii in c(1,3:5)) plotSmoothersIk2(gamList[[rownames(earlyEpiTF[oo13,]) [ii] ]], legendPos=pos[ii], ylim=c(0,10), main=rownames(earlyEpiTF[oo13,])[ii])
dev.off()

```


# Session information

```{r}
sessionInfo()
```
---
title: "10X ik"
author: "Koen Van den Berge"
date: "11/7/2019"
output: html_document
---

```{r}
library(Matrix)
library(aroma.light)
countsFull <- readMM("~/Downloads/GSE128889_RAW/GSM3717977_SCmurinep12_matrix.mtx.gz")
ft <- read.table("~/Downloads/GSE128889_RAW/GSM3717977_SCmurinep12_genes.tsv.gz")
bc <- read.table("~/Downloads/GSE128889_RAW/GSM3717977_SCmurinep12_barcodes.tsv.gz")
rownames(countsFull) <- as.character(ft[,2])
```

# Preprocessing, filtering and trajectory inference

```{r}
keep <- rowSums(countsFull > 1) >= 400
counts <- countsFull[keep,]
rownames(counts) <- ft[keep,2]
countsFQ <- normalizeQuantileRank(as.matrix(counts))
rownames(countsFQ) <- ft[keep,2]

pca <- prcomp(t(log1p(countsFQ)), scale. = FALSE)
rd <- pca$x[,1:8]

plot(rd[,1:2], pch=16, cex=1/2)

plotByGene <- function(rd, geneCount, ng=10, main=NULL, ...){
  pal <- wesanderson::wes_palette("Zissou1", n=ng, type="continuous")
  gg <- Hmisc::cut2(geneCount, g=ng)
  plot(rd, pch=16, cex=1/2, col=pal[gg], main=main, ...)
}

# progenitors
plotByGene(rd[, 1:2], countsFQ["Dpp4",], ng=2, main="Dpp4")
plotByGene(rd[, 1:2], countsFQ["Wnt2",], ng=4, main="Wnt2")

# group 2 cells
plotByGene(rd[, 1:2], countsFQ["Icam1",], ng=4, main="Icam1")
plotByGene(rd[, 1:2], countsFQ["Dlk1",], ng=4, main="Pref1") #Pref1 is also called Dlk1

# group 3 cells
plotByGene(rd[, 1:2], countsFQ["Clec11a",], ng=4, main="Clec11a")

# group 4 cells
plotByGene(rd[, 1:2], countsFull["Wnt6",], ng=2, main="Wnt6")


# cluster
set.seed(9)
cl <- kmeans(rd, centers = 10)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pal <- sample(color, 12)
plot(rd[,1:2], pch=16, cex=1/2, col=pal[cl$cluster])
legend("bottomleft", legend=1:10, pch=16, col=pal, bty='n')

# cluster 3 and 7 are progenitor cells
# cluster 1, 5 and 10 are group 2 cells
# cluster 2 and 6 are group 3 cells
keepCells <- which(cl$cluster %in% c(1,2,3,5,7,10))
countsFQ2 <- countsFQ[,keepCells]
rd2 <- prcomp(t(log1p(countsFQ2)), scale. = FALSE)

### umap
library(umap)
rdUmap <- umap(rd2$x[,1:20])
plot(rdUmap$layout, pch=16, cex=1/3)

rafalib::mypar(mfrow=c(3,2))
# progenitors
plotByGene(rdUmap$layout, countsFQ2["Dpp4",], ng=2, main="Prog: Dpp4")
plotByGene(rdUmap$layout, countsFQ2["Wnt2",], ng=4, main="Prog: Wnt2")

# group 2 cells
plotByGene(rdUmap$layout, countsFQ2["Icam1",], ng=4, main="G2: Icam1")
plotByGene(rdUmap$layout, countsFQ2["Dlk1",], ng=4, main="G2: Pref1") #Pref1 is also called Dlk1

# group 3 cells
plotByGene(rdUmap$layout, countsFQ2["Clec11a",], ng=4, main="G3: Clec11a")

# group 4 cells
plotByGene(rdUmap$layout, countsFull["Wnt6",keepCells], ng=2, main="G4: Wnt6")

# cluster
set.seed(12)
pal <- wesanderson::wes_palette("Darjeeling1", n=6, type="continuous")
cl <- kmeans(rdUmap$layout, centers = 6)
plot(rdUmap$layout, pch=16, cex=1/2, col=pal[cl$cluster])
legend("topleft", legend=1:6, pch=16, col=pal[1:6], bty='n')

# slingshot
library(slingshot)
lin <- getLineages(rdUmap$layout, clusterLabels=cl$cluster, 
                   start.clus=1, end.clus=c(5,2))
plot(rdUmap$layout, pch=16, cex=1/2, col=pal[cl$cluster])
lines(lin, lwd=2)
crv <- getCurves(lin)
plot(rdUmap$layout, pch=16, cex=1/2, col=pal[cl$cluster],
     xlab="UMAP1", ylab="UMAP2", bty='l', asp=1)
lines(crv, lwd=2, col="black")
saveRDS(crv, file="crvUmap.rds")
```

# Evaluate optimal number of knots

8 knots seem appropriate.

```{r}
library(tradeSeq)
set.seed(3)
aicMat1 <- evaluateK(countsFQ2, k=3:10, nGenes=200, sds=crv)
set.seed(4)
aicMat2 <- evaluateK(countsFQ2, k=3:10, nGenes=200, sds=crv)
```

# Model fitting and inference

```{r}
sce <- fitGAM(countsFQ2, sds=crv, nknots=8)
saveRDS(sce, file="~/data/sceAdipose.rds")
```

# Progenitor cell markers other than Dpp4 and Wnt2 (startVsEndTest)

```{r}
seRes <- startVsEndTest(sce)
ose <- order(seRes$waldStat, decreasing=TRUE)
which(rownames(seRes)[ose] == "Wnt2") # ranked 263
which(rownames(seRes)[ose] == "Dpp4") # ranked 450

head(seRes[ose,], n=16)

rafalib::mypar(mfrow=c(3,2))
sapply(1:6, function(ii){
  gene <- rownames(seRes)[ose][ii]
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, asp=1))
})
sapply(7:12, function(ii){
  gene <- rownames(seRes)[ose][ii]
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, asp=1))
})


#figure for paper
rafalib::mypar(mfrow=c(3,2))
genes <- c("Dpp4", "Wnt2", 
           "Pi16", "Akr1c18",
           "Fn1", "Fbn1")
sapply(genes, function(gene){
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, bty='n',
                   xlab="UMAP1", ylab="UMAP2"))
})
```


# Differentiated group 2 and group 3 markers (diffEndTest)

Note that this is only restricted to the end points and is therefore limited.

```{r}
deRes <- diffEndTest(sce)
ode <- order(deRes$waldStat, decreasing=TRUE)


head(seRes[ode,], n=10)

rafalib::mypar(mfrow=c(3,2))
sapply(1:6, function(ii){
  gene <- rownames(deRes)[ode][ii]
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, asp=1))
})
```

# earlyDETest for entire cluster of group 2 and 3 cells

```{r}
plotGeneCount(crv, countsFQ2, models=sce, clusters = cl$cluster)
resedt <- earlyDETest(sce, knots=c(3,7))
oedt <- order(resedt$waldStat, decreasing=TRUE)

head(resedt[oedt,], n=10)

rafalib::mypar(mfrow=c(3,2))
sapply(1:6, function(ii){
  gene <- rownames(resedt)[oedt][ii]
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, asp=1))
})


#figure for paper
rafalib::mypar(mfrow=c(1,2))
genes <- c("Mgp", "Meox2")
sapply(genes, function(gene){
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, bty='n',
                   xlab="UMAP1", ylab="UMAP2", asp=1))
})

#figure for paper
rafalib::mypar(mfrow=c(1,2))
genes <- c("H19", "Col14a1")
sapply(genes, function(gene){
  print(plotByGene(rdUmap$layout, countsFQ2[gene,], ng=4, main=gene, bty='n',
                   xlab="UMAP1", ylab="UMAP2", asp=1))
})
```

---
title: "10X case study"
author: "Koen Van den Berge"
date: "10/15/2019"
output: html_document
---

```{r}
library(tidyverse)
library(dyndimred)
library(slingshot)
library(Seurat)
library(tradeSeq)
library(S4Vectors)
library(SingleCellExperiment)
library(tradeSeq)
```


```{r}
library(Seurat)
seu_trajectory <- readRDS("~/Downloads/adipose_differentiation-merick/seu_trajectory.rds")
feature_info <- read_tsv("~/Downloads/adipose_differentiation-merick/feature_info.tsv")
feature_mapper <- function(x) {feature_info %>% dplyr::slice(base::match(symbol, x)) %>% dplyr::pull(feature_id)}

```


# MDS

```{r}
set.seed(1)

# MDS, we use the landmark mds because it is much faster and memory efficient
lmds <- dyndimred::dimred_landmark_mds(Matrix::t(seu_trajectory@assays$spliced@scale.data))
colnames(lmds) <- paste0("lmds_", seq_len(ncol(lmds)))
lmds_object <- CreateDimReducObject(lmds, key = "lmds_", assay = "spliced")
seu_trajectory@reductions$lmds <- lmds_object
DimPlot(seu_trajectory, reduction = "lmds",pt.size = 0.5, label = TRUE, repel =TRUE)
```
# Clustering

```{r}
set.seed(1)
seu_trajectory <- FindNeighbors(seu_trajectory, verbose = FALSE) %>% 
  FindClusters(resolution = 0.25, verbose = FALSE)
```

# Trajectory inference

```{r}
dimred <- seu_trajectory@reductions$lmds@cell.embeddings
clustering <- seu_trajectory@meta.data$spliced_snn_res.0.25


start_cluster_id <- tibble(
  expression = seu_trajectory@assays$spliced@counts[feature_mapper("Dpp4"), ],
  cluster_id = clustering
) %>% 
  group_by(cluster_id) %>% 
  summarise(expression = mean(expression)) %>% 
  arrange(desc(expression)) %>% 
  pull(cluster_id) %>% 
  dplyr::first() %>% 
  as.character()

set.seed(1)
lineages <- getLineages(dimred, clustering, start.clus = start_cluster_id)
curves <- getCurves(lineages)

plot(dimred, col = RColorBrewer::brewer.pal(9,"Set1")[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = 'black')
```


# EDA

DPP4+ cells are multipotent mesenchymal progenitors according to paper.
However, this is not evident from the current dimension reduction or trajectory.

```{r}
plotGeneCount(curves, counts=counts, gene=feature_mapper("Dpp4"))
plotGeneCount(curves, counts=counts, gene=feature_mapper("Wnt2"))
```


```{r}
library(aroma.light) ; library(umap)
countsFQ <- normalizeQuantileRank(counts)
rdFQ <- umap(t(countsFQ))

plot(rdFQ$layout, pch=16, cex=1/3)
```







# Differential expression

## Evaluate optimal number of knots

```{r}
counts <- seu_trajectory@assays$spliced@counts
keep <- rowSums(counts > 1) >= 100
table(keep) #we keep 3280 genes for analysis.
counts <- as.matrix(counts[keep,])
#aicMat <- evaluateK(counts, sds=curves, k=3:10, nGenes=250, seed=3)
#aicMat2 <- evaluateK(counts, sds=curves, k=3:10, nGenes=250, seed=4)
```

## Fit NB-GAM

```{r}
# sce <- fitGAM(counts = as.matrix(counts), sds = curves,
#               nknots = 7)
# saveRDS(sce, file="20191018_sceTradeSeq.rds")
sce <- readRDS("20191018_sceTradeSeq.rds")
```



## Inference

```{r}

```

