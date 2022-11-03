#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

library(ggplot2); library(dplyr); library(RColorBrewer); library(PCAtools); library(VennDiagram);
library(rafalib); library(sva); library(biomaRt)

#################################################################
##Generating a single expression matrix of FeatureCounts method##
################################################################

setwd(fc_dir)
dir.create(matrices)

filenames = list.files(pattern = ".tsv$")
count_files = lapply(filenames, FUN = function(x) {read.table(x, sep = "\t", header = T)[,7]})
expr_matrix_fc <- do.call("cbind", count_files)
colnames(expr_matrix_fc) <- sub("\\..*", "", filenames)
gene_id <- read.table(filenames[[1]], sep = "\t", header = T)[,1]
rownames(expr_matrix_fc) <- gene_id
rownames(expr_matrix_fc) <- sub("\\..*", "", rownames(expr_matrix_fc))

write.csv2(expr_matrix_fc, file = paste0(matrices, "/raw_expr_matrix_fc.csv"))

setwd(htseq_dir)

filenames = list.files()
count_files_htseq = lapply(filenames, FUN = function(x) {read.table(x, sep = "\t", header = T)[,2]})
expr_matrix_htseq <- do.call("cbind", count_files_htseq)
colnames(expr_matrix_htseq) <- sub("\\..*", "", filenames)
#Deleting last 5 rows as they present summary
n <- dim(expr_matrix_htseq)[1]
expr_matrix_htseq <- expr_matrix_htseq[1:(n-5),]
rownames(expr_matrix_htseq) <- head(read.table(filenames[[1]], sep = "\t", header = F)$V1, -6)
rownames(expr_matrix_htseq) <- sub("\\..*", "", rownames(expr_matrix_htseq))

write.csv2(expr_matrix_htseq, file = paste0(matrices, "/raw_expr_matrix_htseq.csv"))


