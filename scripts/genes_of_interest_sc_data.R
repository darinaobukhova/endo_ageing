#########################################################################################################################
#' Author: Darina Obukhova
#' Lab: Cellular Genomic Medicine, Maastricht University Medical Center (MUMC+)

#' Script purpose: Visualize expression of 57 DEGs discovered in bulk RNA-seq aging data using single cell RNA-seq data
#' of the endometrium across the menstrual cycle (GEO acc code: GSE111976), generated in the paper
#' "Single cell RNA-seq analysis on human endometrium across the natural menstrual cycle" (PMID: 32929266)
#' and see the dynamics of gene expression across the menstrual cycle stages and cell types

#' Input files: (1) DEG names from bulk RNA-seq analysis
#' (2) normalized matrix of single-cell counts; normalized using shifted log transformation from transforGamPoi R package 
#' (3) labels for single-cell data 
#' The paths to them are provided in the config file supplied from the command line

#' Output files: Various plots (in .pdf)

########################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

#' Libraries

library(ggplot2); library(cowplot); library(reshape2); library(tidyverse)

#' Loading the data and subsetting 

genes_of_interest <- read.csv2(DEG_bulk)
load(C1_dataset_transformed)
C1_labels <- read.csv(C1_labels, header = T)
gene_names <- genes_of_interest$gene_s
features <- rownames(C1_sce@assays@data@listData$shifted_log_transform)[which(rownames(C1_sce@assays@data@listData$shifted_log_transform) %in% gene_names)]

#' Remove Day9 as it doesn't have some of the cell types
#' indetified

C1_labels_order <- C1_labels[!C1_labels$day == "9",] 
C1_labels_order <- C1_labels_order[order(C1_labels_order$day),]

#' Subsetting for genes of interest
C1_sce_subset <- C1_sce[,colnames(C1_sce@assays@data@listData$shifted_log_transform) %in% C1_labels_order$X]
C1_sce_subset <- C1_sce_subset[rownames(C1_sce@assays@data@listData$shifted_log_transform) %in% features,]
idx <- match(colnames(C1_sce_subset@assays@data@listData$shifted_log_transform), C1_labels_order$cell_name)
C1_sce_subset_ord <- C1_sce_subset@assays@data@listData$shifted_log_transform[,order(idx)]
print(all(colnames(C1_sce_subset_ord) == C1_labels_order$cell_name))

#' Transforming the C1 subsetted df for futher plotting

C1_sce_subset_ord <- as.data.frame(t(as.matrix(C1_sce_subset_ord)))
C1_sce_subset_ord$Cell <- C1_labels_order$cell_name
C1_sce_subset_ord$Cell_type <- C1_labels_order$cell_type
C1_sce_subset_ord$Day <- C1_labels_order$day
C1_sce_subset_ord_r <- reshape2::melt(C1_sce_subset_ord, id.vars = c("Cell", "Day", "Cell_type"),
                                    measure.vars = features, variable.name = "Feat", value.name = "Expr")

pdf(file=paste0(output_dir, "/sc_markers_1_allsamp.pdf"), height=10, width=8)
a <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[1:12],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, size = 1.5, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        panel.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Day after onset of menstrual bleeding", y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
a
dev.off()

pdf(file=paste0(output_dir, "/sc_markers_2_allsamp.pdf"), height=10, width=8)
b <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[13:24],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Day after onset of menstrual bleeding", y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
b
dev.off()

pdf(file=paste0(output_dir, "/sc_markers_3_allsamp.pdf"), height=10, width=8)
c <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[25:37],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Day after onset of menstrual bleeding", y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
c
dev.off()

pdf(file=paste0(output_dir, "/sc_markers_4_allsamp.pdf"), height=10, width=8)
d <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[38:50],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Day after onset of menstrual bleeding", y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
d
dev.off()

#' Minor modifications such as the removal of y axis label in 2 plots
#' and removal of legend annotation in all plots to plot and adjust it separately

a <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[1:12],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, size = 1.5, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        panel.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0),
        legend.position = "none")+
  labs(x=NULL, y=NULL, color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
a

b <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[13:24],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        panel.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0),
        legend.position = "none")+
  labs(x=NULL, y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
b

c <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[25:37],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        panel.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0),
        legend.position = "none")+
  labs(x="Day after onset of menstrual bleeding", y=NULL, color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
c

d <- ggplot(C1_sce_subset_ord_r[C1_sce_subset_ord_r$Feat %in% features[38:50],], aes(factor(Day), Expr)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        panel.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0), 
        legend.position = "none")+
  labs(x="Day after onset of menstrual bleeding", y="Expression", color="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
d

united_plot <- plot_grid(a, b, c, d, ncol=2, labels=c("a", "b", "c", "d"))

ggsave(filename = paste0(output_dir,"/sc_markers_unitedplot_allsamp.pdf"), 
       plot = united_plot, width = 40, height = 30 , units = "cm")

#Genes of interest during WOI

WOI_genes <- c("RIMKLB", "GALNT12", "GAST", "CDA", "LYPD3", "HPSE", "STC1", "C2CD4A", "C2CD4B")

C1_sce_subset_ord_woi <- C1_sce_subset_ord[,colnames(C1_sce_subset_ord) %in% WOI_genes]
C1_sce_subset_ord_woi$Cell <- C1_labels_order$cell_name
C1_sce_subset_ord_woi$Cell_type <- C1_labels_order$cell_type
C1_sce_subset_ord_woi$Day <- C1_labels_order$day
C1_sce_subset_ord_woi_r <- reshape2::melt(C1_sce_subset_ord_woi, id.vars = c("Cell", "Day", "Cell_type"),
                                      measure.vars = WOI_genes, variable.name = "Feat", value.name = "Expr")
C1_sce_subset_ord_woi_r$Expr <- as.numeric(C1_sce_subset_ord_woi_r$Expr)

pdf(file=paste0(output_dir, "/genes_woi_inclviolinplots.pdf"), height=10, width=8)
genes_woi <- ggplot(C1_sce_subset_ord_woi_r, aes(factor(Day), Expr)) +
  geom_violin(scale = "width", trim = F)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right", breaks=seq(0, 40, 10))+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  labs(x="Time after onset of menstrual bleeding (Day)", y="Expression", colour="Cell type")+
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        panel.spacing.y = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_text(vjust = 5))+
        #legend.position = "none")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
genes_woi
dev.off()

genes_cilia <- c("PPP1R1B", "DNAI1", "AKAP14", "C11orf97", "CAPSL")
C1_sce_subset_ord_cilia <- C1_sce_subset_ord[,colnames(C1_sce_subset_ord) %in% genes_cilia]
C1_sce_subset_ord_cilia$Cell <- C1_labels_order$cell_name
C1_sce_subset_ord_cilia$Cell_type <- C1_labels_order$cell_type
C1_sce_subset_ord_cilia$Day <- C1_labels_order$day
C1_sce_subset_ord_cilia_r <- reshape2::melt(C1_sce_subset_ord_cilia, id.vars = c("Cell", "Day", "Cell_type"),
                                          measure.vars = genes_cilia, variable.name = "Feat", value.name = "Expr")
C1_sce_subset_ord_cilia_r$Expr <- as.numeric(C1_sce_subset_ord_cilia_r$Expr)

pdf(file=paste0(output_dir, "/genes_cilia_inclvioliplots.pdf"), height=10, width=8)
genes_cilia <- ggplot(C1_sce_subset_ord_cilia_r, aes(factor(Day), Expr)) +
  geom_violin(scale = "width", trim = F)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right", breaks=seq(0, 30, 10))+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        panel.spacing.y = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Time after onset of menstrual bleeding (Day)", y="Expression", colour="Cell type")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
genes_cilia
dev.off()

#' Remove color legend and y axis label in the plot on WOI genes for
#' a better readability of the combined plot

genes_woi <- ggplot(C1_sce_subset_ord_woi_r, aes(factor(Day), Expr)) +
  geom_violin(scale = "width", trim = F)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0, 0), position="right", breaks=seq(0, 40, 10))+
  geom_jitter(aes(colour=Cell_type), alpha = 0.7, position = position_dodge2(width=0.75))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  labs(x="Time after onset of menstrual bleeding (Day)", y=NULL, colour="Cell type")+
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        panel.spacing.y = unit(0.2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_text(vjust = 5),
        legend.position = "none")+
  scale_color_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "cyan", "#8B008B", "gray27"))
genes_woi


united_plot <- plot_grid(genes_woi, genes_cilia, ncol=2, labels=c("a", "b"))

ggsave(filename = paste0(output_dir,"/selected_sc_markers_unitedplot_allsamp_inclviolinplots.pdf"), 
       plot = united_plot, width = 30, height = 20, units = "cm")
