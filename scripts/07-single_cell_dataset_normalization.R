#!/usr/bin/env Rscript

###########################################################################################
#' Author: Darina Obukhova
#' Lab: Cellular Genomic Medicine, Maastricht University Medical Center (MUMC+)
#' Script purpose: Perform single cell dataset (GEO acc code: GSE111976) for later plotting

#' Input file: Single cell C1 (Fluidigm) dataset
#' Output file: sce object (.RData)

########################################################################################### 

options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)

# Arguments
config <- args[1]
source(as.character(config))

# Loading libraries

library(transformGamPoi)
library(scran)
library(SingleCellExperiment)
library(scRNAseq)

# Single cell data
cat("Loading single cell dataset...")
C1 <- read.csv(C1_data, header = T)
C1_labels <- read.csv(C1_labels, header = T)

#Removing X in column names
colnames(C1) <- sub("^X", "", colnames(C1))
rownames(C1) <- C1[, 1]

# Remove first column
C1 <- C1[, -1]

# Gene level filtering

cat("Perfoming gene level filtering:")
# Remove genes with 0 counts in all the cells
C1_filt <- C1[rowSums(C1[]) > 0, ]
sprintf("%i genes with 0 counts across all the cells were removed", nrow(C1) - nrow(C1_filt))
sprintf("Now the expression matrix has %i genes", nrow(C1_filt))

# Keep genes that are 'detectable' (have at least 10 counts) in at least 5% of cells
keep_genes <- rowSums(C1_filt >= 0) >= round(0.05 * ncol(C1_filt))

C1_filt2 <- C1_filt[keep_genes,]
sprintf("%i genes expressed in less than 10 cells were removed", nrow(C1_filt) - nrow(C1_filt2))
sprintf("Now the expression matrix has %i genes", nrow(C1_filt2))

# C1 data as matrix
C1_filt_mat <- as.matrix(C1_filt2, drop=F)

# C1 data as sce object for latter convenience

C1_sce <- SingleCellExperiment::SingleCellExperiment(list(counts=C1_filt_mat))

assay(C1_sce, "shifted_log_transform") <- transformGamPoi::shifted_log_transform(C1_sce, 
                                                                                 size_factors = "deconvolution")

save(C1_sce, file=paste0(outputdir, "C1_transformed.RData"))

