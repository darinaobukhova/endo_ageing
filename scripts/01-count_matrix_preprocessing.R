#!/usr/bin/env Rscript

options(scipen=999)
args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

# Loading packages

library(dplyr)
library(tidyverse)
library(biomaRt)
library(sva)


# Generating a single expression matrix from all files output from FeatureCounts

setwd(fc_dir)
dir.create(matrices)

filenames <- list.files(pattern = ".tsv$")
count_files <- lapply(filenames, FUN = function(x) {read.table(x, sep = "\t", header = T)[,7]})
expr_matrix_fc <- do.call("cbind", count_files)
colnames(expr_matrix_fc) <- sub("\\..*", "", filenames)
gene_id <- read.table(filenames[[1]], sep = "\t", header = T)[,1]
rownames(expr_matrix_fc) <- gene_id
rownames(expr_matrix_fc) <- sub("\\..*", "", rownames(expr_matrix_fc))

write.csv2(expr_matrix_fc, file = paste0(matrices, "/raw_expr_matrix_fc.csv"))

# Performing batch correction for a sequencing run taking into condsideration their receptivity status

samplesheet <- read.csv2(samplesheet, header = T)

samplesheet <- samplesheet[order(match(samplesheet[,2], colnames(expr_matrix_fc))), ]
print(all(samplesheet[,2] == colnames(expr_matrix_fc)))

batch <- samplesheet[,5]
group <- samplesheet[,3]

corrected_data_fc <- sva::ComBat_seq(as.matrix(expr_matrix_fc), batch=batch, group=group)

# Remove rows where values are >=5 in >75% samples

corrected_data_fc_filtered <- corrected_data_fc[rowSums(corrected_data_fc >= 5) >= round(ncol(corrected_data_fc)/100 * 75),]

write.csv2(corrected_data_fc_filter, file = paste0(matrices, "/corrected_data_fc_filter.csv"))

# Save data 

save(samplesheet, expr_matrix_fc, corrected_data_fc_filtered, file=paste0(cor_data, "/corrected_countdata_fc_allsamp.RData"))
