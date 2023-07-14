#!/usr/bin/env Rscript

#########################################################################################################################
#' Author: Darina Obukhova
#' Lab: Cellular Genomic Medicine, Maastricht University Medical Center (MUMC+)

#' Script purpose: To perform gene expression deconvolution of bulk RNA-seq profiles of young and advanced maternal age
#' endometrium to infer their cell-type composition using Dampened Weighted Least Squares (DWLS) method, published in
#' the paper "Accurate estimation of cell-type composition from gene expression data" (PMID: 31278265) and performing
#' well in benchmarking studies in terms of the accurate detection across diverse cell types. The method uses WLS with
#' weights dampened by a certain dampening constant. 

#' Input files: (1) bulk RNA-seq matrix
#' (2) signature matrices for single-cell datasets

#' Output files: tables with cellular proportions of all major endometrial cell types in bulk data (.csv)
#########################################################################################################################

options(scipen=999)
set.seed(8)

# Arguments

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

# Loading libraries

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(DWLS)
library(scales)

# Loading data

cat("Loading data...")
load(data_for_deconvolution)
load(signature_c1)
load(signature_tenx)
samplesheet <- read.csv2(samplesheet, header = T)

# Important to have data object as a double type for later processing

# Order columns alphabetically

corrected_data_fc_filtered <- corrected_data_fc_filtered[, order(colnames(corrected_data_fc_filtered), decreasing = T)]
nms_bulk <- rownames(corrected_data_fc_filtered)
corrected_data_fc_filtered_db <- as.double(corrected_data_fc_filtered)
names(corrected_data_fc_filtered_db) <- nms_bulk

# Deconvolution

dir.create(file.path(outputdir, "results", gsub("-", "_",Sys.Date())))
setwd(file.path(outputdir, "results", gsub("-", "_",Sys.Date())))

# For C1 dataset

cat("Deconvolving cellular types in the bulk dataset using 10x dataset as the reference...")

# Using 2 methods of DWLS package: dampened weighted least squares and a support vector regression

# Raw data

allProport_DWLS <- NULL
allProport_SVR<- NULL

for(j in 1:(ncol(corrected_data_fc_filtered))){
  sign <- signature_tenx_midsec
  bulk <- corrected_data_fc_filtered[,j]
  names(bulk) <- rownames(corrected_data_fc_filtered)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS <- cbind(allProport_DWLS, proport_DWLS)
  allProport_SVR <- cbind(allProport_SVR, proport_SVR)
}

#Modifying for plotting purposes later on

allProport_SVR <- as.data.frame(t(allProport_SVR))
rownames(allProport_SVR) <- samplesheet$SampleID
write.csv(allProport_SVR, file="cellprop_tenx_midsec_svr_raw.csv")
allProport_SVR$Group <- samplesheet$Group
allProport_SVR <- allProport_SVR[order(allProport_SVR$Group, decreasing=T),]
allProport_DWLS <- as.data.frame(t(allProport_DWLS))
rownames(allProport_DWLS) <- samplesheet$SampleID
write.csv2(allProport_DWLS, file="cellprop_tenx_midsec_dwls_raw.csv")
allProport_DWLS$Group <- samplesheet$Group
allProport_DWLS <- allProport_DWLS[order(allProport_DWLS$Group, decreasing=T),]

save(allProport_SVR, allProport_DWLS, file = "deconv_res_tenxref.RData")

cat("Performing statistics tests...")

comp_dwls_c1_raw <- apply(allProport_DWLS[, -ncol(allProport_DWLS)], 2, function(x) wilcox.test(x[1:12], x[12:24], exact = FALSE, 
                                                                                                paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()

comp_svr_c1_raw <- apply(allProport_SVR[, -ncol(allProport_SVR)], 2, function(x) wilcox.test(x[1:12], x[12:24], exact = FALSE, 
                                                                                             paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()

write.csv(comp_dwls_c1_raw, file="prop_DWLS_tenx_raw.csv")
write.csv(comp_svr_c1_raw, file="prop_svr_tenx_raw.csv")

