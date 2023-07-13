#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

set.seed(10)

setwd(matrices)

# Load Libraries

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(rafalib)


# Load objects

load(expr_matrices)

# PCA on the uncorrected data

dir.create(output_dir_eda)

# Remove rows where values are >=5 in >75% samples

expr_matrix_fc_filtered <- expr_matrix_fc_filtered[rowSums(expr_matrix_fc_filtered >= 5) >= round(ncol(expr_matrix_fc_filtered)/100 * 75),]

expr_matrix_fc_filtered <- expr_matrix_fc_filtered + 1
expr_matrix_fc_filtered_log <- log2(expr_matrix_fc_filtered)

pcDat <- prcomp(as.matrix(t(expr_matrix_fc_filtered_log)), scale = F, center = T)
var_explained <- pcDat$sdev^2/sum(pcDat$sdev^2)

pdf(file=paste0(output_dir_eda, "/fc_uncorrected_allsamp_PCA_labeled.pdf"), height=7, width=9)

# Trim rownames for a cleaner appearance on a plot
 
theme_set(theme_bw())
pcDat$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(samplesheet[,3]), shape = factor(samplesheet[,5]))) + 
  ggrepel::geom_text_repel(aes(label = rownames(pcDat$x))) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank()) + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100, 1), "%"), 
       title = "Uncorrected data") +
  scale_color_manual(values=ifelse(samplesheet[,3] == "YMA", "darkorchid4", "darkorange2")) +
  guides(color = guide_legend(title="Group",override.aes = list(size=4, linetype=0), shape = guide_legend(title="Sequencing Run"))) +
  stat_ellipse(level=0.95, linetype=2)

dev.off()

# PCA on the batch-corrected data 

corrected_data_fc_filtered <- corrected_data_fc_filtered + 1
corrected_data_fc_filtered_log <- log2(corrected_data_fc_filtered)

pcDat <- prcomp(as.matrix(t(corrected_data_fc_filtered_log)), scale = F, center = T)
var_explained <- pcDat$sdev^2/sum(pcDat$sdev^2)

pdf(file=paste0(output_dir_eda, "/fc_corrected_allsamp_PCA_labeled.pdf"), height=7, width=9)

#Trim rownames for a cleaner appearance on a plot 

theme_set(theme_bw())

pcDat$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(samplesheet[,3]), shape = factor(samplesheet[,5]))) + 
  ggrepel::geom_text_repel(aes(label = rownames(pcDat$x))) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank()) + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100, 1), "%"),
       title = "Batch corrected data") +
  scale_color_manual(values=ifelse(samplesheet[,3] == "YMA", "darkorchid4", "darkorange2")) +
  guides(color = guide_legend(title="Group",override.aes = list(size=4, linetype=0), shape = guide_legend(title="Sequencing Run"))) +
  stat_ellipse(level=0.95, linetype=2)

dev.off()

