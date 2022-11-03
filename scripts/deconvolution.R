#!/usr/bin/env Rscript

options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)

#Arguments
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

cat("Loading libraries...")
library(dplyr); library(ggplot2); library(ggpubr); library(data.table); library(DWLS);library(scales)

####################
####Loading data####
####################

cat("Loading data...")

load(matrices)
load(signatures)

#####################
####Deconvolution####
#####################

#For C1 dataset
cat("Deconvolving cellular types in the bulk dataset using C1 dataset as the reference...")

#Using 2 methods of DWLS package: dampened weighted least squares and a support vector regression

#Raw data

allProport_DWLS <- NULL
allProport_SVR  <- NULL

for(j in 1:(ncol(expr_matrix_filtered2))){
  sign <- signature_C1_epaggr
  bulk <- expr_matrix_filtered2[,j]
  names(bulk) <- rownames(expr_matrix_filtered2)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS <- cbind(allProport_DWLS, proport_DWLS)
  allProport_SVR <- cbind(allProport_SVR, proport_SVR)
}

#Normalized data

allProport_DWLS_norm_tmm <- NULL
allProport_SVR_norm_tmm  <- NULL

for(j in 1:(ncol(normalized_expr_matrix_filtered2_sampnames))){
  sign <- signature_C1_epaggr
  bulk <- normalized_expr_matrix_filtered2_sampnames[,j]
  names(bulk) <- rownames(normalized_expr_matrix_filtered2_sampnames)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS_norm_tmm <- cbind(allProport_DWLS_norm_tmm, proport_DWLS)
  allProport_SVR_norm_tmm <- cbind(allProport_SVR_norm_tmm, proport_SVR)
}

#Modifying for plotting purposes later on

allProport_SVR <- as.data.frame(t(allProport_SVR))
rownames(allProport_SVR) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_SVR, file=paste0(outputdir, "cellprop_C1_svr_raw_epaggr.csv"))
allProport_SVR$Group <- colnames(expr_matrix_filtered2)
allProport_DWLS <- as.data.frame(t(allProport_DWLS))
rownames(allProport_DWLS) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_DWLS, file=paste0(outputdir, "cellprop_C1_dwls_raw_epaggr.csv"))
allProport_DWLS$Group <- colnames(expr_matrix_filtered2)

allProport_SVR_norm_tmm <- as.data.frame(t(allProport_SVR_norm_tmm))
rownames(allProport_SVR_norm_tmm) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_SVR_norm_tmm, file=paste0(outputdir, "cellprop_C1_svr_norm_tmm_epaggr.csv"))
allProport_SVR_norm_tmm$Group <- colnames(normalized_expr_matrix_filtered2)
allProport_DWLS_norm_tmm <- as.data.frame(t(allProport_DWLS_norm_tmm))
rownames(allProport_DWLS_norm_tmm) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_DWLS_norm_tmm, file=paste0(outputdir, "cellprop_C1_dwls_norm_tmm_epaggr.csv"))
allProport_DWLS_norm_tmm$Group <- colnames(normalized_expr_matrix_filtered2)

cat("Finished plotting the results of C1 as the reference..")

save(allProport_SVR, allProport_SVR_norm_tmm, allProport_DWLS, allProport_DWLS_norm_tmm, 
     file = paste0(outputdir, "deconv_res_C1ref_tmm_eppagr.RData"))

cat("Deconvolving using 10x dataset as the reference..")

#Raw data

allProport_DWLS_tenx <- NULL
allProport_SVR_tenx  <- NULL

for(j in 1:(ncol(expr_matrix_filtered2))){
  sign <- signature_tenx_epaggr
  bulk <- expr_matrix_filtered2[,j]
  names(bulk) <- rownames(expr_matrix_filtered2)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS_tenx <- cbind(allProport_DWLS_tenx, proport_DWLS)
  allProport_SVR_tenx <- cbind(allProport_SVR_tenx, proport_SVR)
}

#Log normalized data

allProport_DWLS_tenx_norm_tmm <- NULL
allProport_SVR_tenx_norm_tmm  <- NULL

for(j in 1:(ncol(normalized_expr_matrix_filtered2_sampnames))){
  sign <- signature_tenx_epaggr
  bulk <- normalized_expr_matrix_filtered2_sampnames[,j]
  names(bulk) <- rownames(normalized_expr_matrix_filtered2_sampnames)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS_tenx_norm_tmm <- cbind(allProport_DWLS_tenx_norm_tmm, proport_DWLS)
  allProport_SVR_tenx_norm_tmm <- cbind(allProport_SVR_tenx_norm_tmm, proport_SVR)
}

#Modifying for plotting purposes later on

allProport_SVR_tenx <- as.data.frame(t(allProport_SVR_tenx))
rownames(allProport_SVR_tenx) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_SVR_tenx, file=paste0(outputdir, "cellprop_tenx_svr_raw_epaggr.csv"))
allProport_SVR_tenx$Group <- colnames(expr_matrix_filtered2)
allProport_DWLS_tenx <- as.data.frame(t(allProport_DWLS_tenx))
rownames(allProport_DWLS_tenx) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_DWLS_tenx, file=paste0(outputdir, "cellprop_tenx_dwls_raw_epaggr.csv"))
allProport_DWLS_tenx$Group <- colnames(expr_matrix_filtered2)

allProport_SVR_tenx_norm_tmm <- as.data.frame(t(allProport_SVR_tenx_norm_tmm))
rownames(allProport_SVR_tenx_norm_tmm) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_SVR_tenx_norm_tmm, file=paste0(outputdir, "cellprop_tenx_svr_norm_tmm_epaagr.csv"))
allProport_SVR_tenx_norm_tmm$Group <- colnames(normalized_expr_matrix_filtered2)
allProport_DWLS_tenx_norm_tmm <- as.data.frame(t(allProport_DWLS_tenx_norm_tmm))
rownames(allProport_DWLS_tenx_norm_tmm) <- colnames(normalized_expr_matrix_filtered2_sampnames)
write.csv(allProport_DWLS_tenx_norm_tmm, file=paste0(outputdir, "cellprop_tenx_dwls_norm_tmm_epaggr.csv"))
allProport_DWLS_tenx_norm_tmm$Group <- colnames(normalized_expr_matrix_filtered2)

save(allProport_SVR_tenx, allProport_SVR_tenx_norm_tmm, allProport_DWLS_tenx, allProport_DWLS_tenx_norm_tmm, 
     file = paste0(outputdir, "deconv_res_tenxref_tmm_epaggr.RData"))

#Stats

cat("Performing statistics tests...")

comp_dwls_c1_raw <- apply(allProport_DWLS[, -ncol(allProport_DWLS)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                        paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()

comp_svr_c1_raw <- apply(allProport_SVR[, -ncol(allProport_SVR)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                     paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()
comp_dwls_c1_norm <- apply(allProport_DWLS_norm_tmm[, -ncol(allProport_DWLS_norm_tmm)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                 paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()

comp_svr_c1_norm <- apply(allProport_SVR_norm_tmm[, -ncol(allProport_SVR_norm_tmm)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                              paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()
comp_dwls_tenx_raw <- apply(allProport_DWLS_tenx[, -ncol(allProport_DWLS_tenx)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()

comp_svr_tenx_raw <- apply(allProport_SVR_tenx[, -ncol(allProport_SVR_tenx)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                             paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()
comp_dwls_tenx_norm <- apply(allProport_DWLS_tenx_norm_tmm[, -ncol(allProport_DWLS_tenx_norm_tmm)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                        paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()

comp_svr_tenx_norm <- apply(allProport_SVR_tenx_norm_tmm[, -ncol(allProport_SVR_tenx_norm_tmm)], 2, function(x) wilcox.test(x[1:12], x[12:23], exact = FALSE, 
                                                                                                     paired=FALSE, p.adjust = "BH")$p.value) %>% as.data.frame()

write.csv(comp_dwls_c1_raw, file=paste0(outputdir, "prop_DWLS_C1_raw_epaggr.csv"))
write.csv(comp_svr_c1_raw, file=paste0(outputdir, "prop_svr_C1_raw_epaggr.csv"))
write.csv(comp_dwls_c1_norm, file=paste0(outputdir, "prop_DWLS_C1_norm_tmm_epaggr.csv"))
write.csv(comp_svr_c1_norm, file=paste0(outputdir, "prop_svr_C1_norm_tmm_epaggr.csv"))
write.csv(comp_dwls_tenx_raw, file=paste0(outputdir, "prop_DWLS_tenx_raw_epaggr.csv"))
write.csv(comp_svr_tenx_raw, file=paste0(outputdir, "prop_svr_tenx_raw_epaggr.csv"))
write.csv(comp_dwls_tenx_norm, file=paste0(outputdir, "prop_DWLS_tenx_norm_tmm.csv"))
write.csv(comp_svr_tenx_norm, file=paste0(outputdir, "prop_svr_tenx_norm_tmm.csv"))
