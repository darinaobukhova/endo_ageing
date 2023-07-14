#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

# Loading libraries

library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer) 
library(DESeq2)
library(pheatmap)
library(ggdendro)
library(sva)
library(rafalib)
library(pvclust)

# Loading files

load(corrected_data)
samplesheet <- read.csv2(samplesheet, header=T)
samplesheet <- samplesheet[order(match(samplesheet[,2], colnames(corrected_data_fc))), ]
print(all(samplesheet[,2] == colnames(corrected_data_fc)))

# Prepare for variance stabilizing transformation

dds_fc <- DESeqDataSetFromMatrix(countData = as.matrix(corrected_data_fc_filtered), 
                                 colData = samplesheet, 
                                 design = ~Receptivity_Score + Group)

dds_fc$Group <- relevel(dds_fc$Group, ref="YMA")

dds <- DESeq(dds_fc)

# Variancestabilizing transformation

vst <- varianceStabilizingTransformation(dds)
vsd <- assay(vst)

# Perform sample clustering
# Sample-to-sample distances

sampleDists <- dist(t(vsd), method="euclidean")
sampleDistsMat <- as.matrix(sampleDists)

clusterSample <- hclust(sampleDists, method="average")
samplemod <- as.dendrogram(clusterSample)

sampleDendrogramData <- segment(dendro_data(samplemod, type = "rectangle"))

pdf(file=paste0(outputdir, "/samples_hclust_vst_average.pdf"), height=7, width=9)
plot(clusterSample, xlab=NULL)
dev.off()

clusterSample_c <- hclust(sampleDists, method="complete")
samplemod <- as.dendrogram(clusterSample_c)

sampleDendrogramData <- segment(dendro_data(samplemod, type = "rectangle"))

pdf(file=paste0(outputdir, "/samples_hclust_vst_complete.pdf"), height=7, width=9)
plot(clusterSample_c, xlab=NULL)
dev.off()

#Sample clustering on log2 cpm values

cpm <- t( round(t(corrected_data_fc_filtered) / colSums(corrected_data_fc_filtered) * 1e6))
logcpm <- log2(cpm + 1)

d <- dist(t(logcpm), method="euclidean")

# Compute sample correlations
sample_cor <- cor( logcpm )
round(sample_cor,4)

# Transform the scale from correlations
cor_distance <- -(sample_cor - 1)/2
round(cor_distance,4)

# Convert it to a distance object
d2 <- as.dist(cor_distance)

pdf(file=paste0(outputdir, "/samples_hclust_logcpm2.pdf"), height=12, width=12)
mypar(1,2,mar=c(6,4,2,1))
h <- hclust(d, method="average")
plot( as.dendrogram(h) , las=1, main="D=Euclidean")
points(1:ncol(corrected_data_fc_filtered) ,rep(0,ncol(corrected_data_fc_filtered)), pch= 16, cex=2, col=ifelse(samplesheet$Group == 'YMA', "darkorchid4", "darkorange2")[h$order])
h2 <- hclust(d2, method="average")
plot( as.dendrogram(h2) , las=1, main="D=Correlation")
points(1:ncol(corrected_data_fc_filtered) ,rep(0,ncol(corrected_data_fc_filtered)), pch= 16, cex=2, col=ifelse(samplesheet$Group == 'YMA', "darkorchid4", "darkorange2")[h$order])
dev.off()

# Clustering with bootstrapping

pvc <- pvclust(data = corrected_data_fc_filtered, method.dist = "correlation", method.hclust = "average", nboot=1000, parallel = T)

pdf(file=paste0(outputdir, "/samp_clustering_bootstrapping_average_corr.pdf"), height=12, width=12)
plot(pvc,las=2,hang = -0.5)
pvrect(pvc, alpha = 0.95)
points(1:ncol(corrected_data_fc_filtered) ,rep(0,ncol(corrected_data_fc_filtered)), pch= 16, cex=2, col=ifelse(samplesheet$Group == 'YMA', "darkorchid4", "darkorange2")[pvc$hclust$order])
dev.off()
