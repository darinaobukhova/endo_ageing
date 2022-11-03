#!/usr/bin/env Rscript

options(scipen=999)
options(digits=2)

args <- commandArgs(trailingOnly=TRUE)

#Arguments
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

cat("Loading libraries...")
library(dplyr); library(ggplot2); library(ggpubr); library(tidyr); library(tidyverse); library(cowplot);
library(ggsci); library(scales)

####################
####Loading data####
####################

cat("Loading data...")
load(deconv_C1)
load(deconv_tenx)

################
####Plotting####
################

cat("Plotting...")

#C1 as the reference
deconv_plots_c1_log_svr <- lapply(colnames(allProport_SVR[, -ncol(allProport_SVR)]), function(var) {
  ggpar(ggviolin(allProport_SVR, x = "Group", y = var, fill = "Group", 
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"))
        
  )})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_C1_SVR_raw_epaggr.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_c1_log_svr, ncol = 3)

dev.off()

allProport_DWLS <- round(allProport_DWLS[,-ncol(allProport_DWLS)], digits = 5)

deconv_plots_c1_log_dwls <- lapply(colnames(allProport_DWLS[, -ncol(allProport_DWLS)]), function(var) {
  ggpar(ggviolin(allProport_DWLS, x = "Group", y = var, fill = "Group", 
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.8, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"), axis.text = element_text(size=12),
                axis.title = element_text(size=12))
        
  )})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_C1_midsec_DWLS_raw.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_c1_log_dwls, ncol = 3)

dev.off()

#10x as the reference
deconv_plots_tenx_log_svr <- lapply(colnames(allProport_SVR_tenx[, -ncol(allProport_SVR_tenx)]), function(var) {
  ggpar(ggviolin(allProport_SVR_tenx, x = "Group", y = var, fill = "Group", 
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.5, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"))
        
  )})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_10x_SVR_raw_epaggr.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_tenx_log_svr, ncol = 3)

dev.off()

deconv_plots_tenx_log_dwls <- lapply(colnames(allProport_DWLS_tenx[, -ncol(allProport_DWLS_tenx)]), function(var) {
  ggpar(ggviolin(allProport_DWLS_tenx, x = "Group", y = var, fill = "Group", scale = "width",
                 palette = c("darkorchid4", "darkorange2"), add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.8, width = 0.8, title = var) +
          geom_jitter(alpha = 0.2, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"))
        
  )})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_10x_midsec_DWLS_raw.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_tenx_log_dwls, ncol = 3)

dev.off()

#Stacked barchats

allProport_DWLS_c1_bc <- gather(allProport_DWLS, 1:6, key="Cell_type", value="Proportion")
allProport_DWLS_c1_bc$Group <- factor(allProport_DWLS_c1_bc$Group, levels=c("YMA", "AMA"))

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_C1_midsec_DWLS_raw_cellprop.pdf"), height=7, width=7)
ggplot(allProport_DWLS_c1_bc, aes(fill=Cell_type, y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_aaas() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

allProport_SVR_c1_bc <- gather(allProport_SVR, 1:5, key="Cell_type", value="Proportion")

svr_c1 <- ggplot(allProport_SVR_c1_bc, aes(fill=Cell_type, y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_aaas() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

allProport_DWLS_tenx_bc <- gather(allProport_DWLS_tenx, 1:8, key="Cell_type", value="Proportion") 
allProport_DWLS_tenx_bc$Group <- factor(allProport_DWLS_tenx_bc$Group, levels=c("YMA", "AMA"))

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_10x_midsec_DWLS_raw_cellprop.pdf"), height=7, width=7)
ggplot(allProport_DWLS_tenx_bc, aes(fill=Cell_type, y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_aaas() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

allProport_SVR_tenx_bc <- gather(allProport_SVR_tenx, 1:6, key="Cell_type", value="Proportion")

svr_tenx <- ggplot(allProport_SVR_tenx_bc, aes(fill=Cell_type, y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_aaas() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

united_plot <- plot_grid(dwls_c1, dwls_tenx, nrow=1, labels="AUTO")

ggsave(filename = paste0(outputdir,"CellProport_barcharts_DWLS_raw_epaggr.pdf"), 
       plot = united_plot, width = 25, height = 15, units = "cm")

united_plot2 <- plot_grid(svr_c1, svr_tenx, nrow=1, labels="AUTO")
ggsave(filename = paste0(outputdir,"CellProport_barcharts_SVR_raw_epaggr.pdf"), 
       plot = united_plot2, width = 25, height = 15, units = "cm")
