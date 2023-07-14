#!/usr/bin/env Rscript

#########################################################################################################################
#' Author: Darina Obukhova
#' Lab: Cellular Genomic Medicine, Maastricht University Medical Center (MUMC+)

#' Script purpose: Create gene expression signature matrices from endometrium single-cell data subsets corresponding to the midsecretory
#' phase of the menstrual cycle (GEO acc code: GSE111976), generated in the paper
#' "Single cell RNA-seq analysis on human endometrium across the natural menstrual cycle" (PMID: 32929266). It is the
#' first and most computation intensive part of the Dampened Weighted Least Squares (DWLS) method. Recommended to be
#' carried out on the HPC.

#' Input files: (1) C1 scRNA-seq dataset (Fluidigm) from GSE111976 
#' (2) labels for C1 scRNA-seq dataset
#' (3) 10X scRNA-seq dataset (10X Genomics) from GSE111976
#' (4) labels for 10X scRNA-seq dataset 
#' The paths to them are provided in the config file supplied from the command line

#' Output files: Two signature matrices of C1 and 10X datasets, respectively (.RData)

########################################################################################################################

options(scipen=999)
set.seed(746)

args <- commandArgs(trailingOnly=TRUE)

# Arguments

config <- args[1]
source(as.character(config))

# Loading packages

library(dplyr)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(DWLS)
library(ggpubr)
library(Seurat)

dir.create(file.path(outputdir, "results", Sys.Date()), showWarnings = F)
setwd(file.path(outputdir, "results", Sys.Date()))

# Single cell data

cat("Loading single cell dataset...")
C1 <- read.csv(C1_data, header = T)
C1_labels <- read.csv(C1_labels, header = T)

tenx <- readRDS(tenx_data)
tenx_labels <- read.csv(tenx_labels, header = T)

#Filtering of single cell (C1 data)

cat("Pre-processing single cell dataset (C1)")

#Removing X in column names

colnames(C1) <- sub("^X", "", colnames(C1))
rownames(C1) <- C1[, 1]

#Remove first column

C1 <- C1[, -1]

#Subsetting cells, corresponding to the mid-secretory phase of the menstrual cycle

C1_midsec_labels <- C1_labels[C1_labels$day == '21' | C1_labels$day == '22' | C1_labels$day == '24', ]
C1_midsec <- C1[,colnames(C1) %in% C1_midsec_labels$cell_name]

#Gene level filtering

cat("Perfoming gene level filtering:")

#Remove genes with 0 counts in all the cells
C1_filt <- C1_midsec[rowSums(C1[]) > 0, ]
sprintf("%i genes with 0 counts across all the cells were removed", nrow(C1_midsec) - nrow(C1_filt))
sprintf("Now the expression matrix has %i genes", nrow(C1_filt))

#Keep genes that are 'detectable' (have at least 10 counts) in at least 5% of cells

keep_genes <- rowSums(C1_filt >= 0) >= round(0.05 * ncol(C1_filt))

C1_filt2 <- C1_filt[keep_genes,]
sprintf("%i genes expressed in less than 10 cells were removed", nrow(C1_filt) - nrow(C1_filt2))
sprintf("Now the expression matrix has %i genes", nrow(C1_filt2))

#Some editing of the labels

C1_midsec_labels$cell_type <- sub(' ', '_', C1_midsec_labels$cell_type)

#Filtering of single cell (10X data)

cat("Pre-processing single cell dataset (10x)...")

seurt <- Seurat::CreateSeuratObject(counts = tenx, min.cells = round(0.05 * length(tenx@Dimnames[[2]])), min.features = 350)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurt@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(seurt@assays$RNA@data[mito.genes, ])/Matrix::colSums(seurt@assays$RNA@data)
seurt[["percent.mito"]] <- PercentageFeatureSet(seurt, pattern = "^MT-")

seurt_filt <- subset(seurt, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mito < 0.15)
seurt_filt <- NormalizeData(object = seurt, normalization.method = "LogNormalize")

#Subsetting cells, corresponding to the mid-secretory phase of the menstrual cycle

tenx_midsec_labels <- tenx_labels[tenx_labels$day == '20' | tenx_labels$day == '22' | tenx_labels$day == '23', ]
seurt_midsec <- seurt_filt[,colnames(seurt_filt) %in% tenx_midsec_labels$cell_name]

#Some editing of the labels

tenx_midsec_labels$cell_type <- gsub(' ', '_', tenx_midsec_labels$cell_type)

#Build signature from single-cell data (C1)

cat("Building signature matrix from C1 data...")

signature_C1_midsec <- buildSignatureMatrixMAST(scdata = C1_filt2, id = C1_midsec_labels$cell_type, path = getwd(),
diff.cutoff = 0.5, pval.cutoff = 0.01)

save(signature_C1_midsec, file= "signature_C1_midsec.RData")

sprintf("Signature matrix from C1 data has %i genes", dim(signature_C1_midsec)[1])

#Build signature from single-cell data (10X)

cat("Building signature matrix from 10x data...")

signature_tenx_midsec <- buildSignatureMatrixUsingSeurat(scdata = seurt_midsec@assays$RNA@counts, 
id = tenx_midsec_labels[tenx_midsec_labels$X %in% colnames(seurt_midsec@assays$RNA@counts),]$cell_type,
path = getwd(), diff.cutoff = 0.5, pval.cutoff = 0.01)

save(signature_tenx_midsec, file = "signature_tenx_midsec.RData")

sprintf("Signature matrix from 10x data has %i genes", dim(signature_tenx_midsec)[1])

