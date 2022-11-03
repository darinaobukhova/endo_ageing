#Vizualization 

genes_of_interest <- readxl::read_xlsx("/Users/darinaobukhova/projects/Endometrial_Ageing/reanalysis_070622/EA_protein_genes_logFC1.xlsx")

gene_names <- genes_of_interest$gene_name
features <- rownames(C1_sce@assays@data@listData$shifted_log_transform)[which(rownames(C1_sce@assays@data@listData$shifted_log_transform) %in% gene_names)][41:50]

#Subsetting for genes of interest
C1_sce_subset <- C1_sce[rownames(C1_sce@assays@data@listData$shifted_log_transform) %in% features]

#Order labels according to the day
C1_labels_order <- C1_labels[order(C1_labels$day),]

idx <- match(colnames(C1_sce_subset@assays@data@listData$shifted_log_transform), C1_labels_order$cell_name)
C1_sce_subset_ord <- C1_sce_subset@assays@data@listData$shifted_log_transform[,order(idx)]
print(all(colnames(C1_sce_subset_ord) == C1_labels_order$cell_name))

#Renaming columns as day
#colnames(C1_sce_subset_ord) <- C1_labels_order$day
C1_sce_subset_ord <- as.data.frame(t(as.matrix(C1_sce_subset_ord)))
C1_sce_subset_ord$Cell <- C1_labels_order$cell_name
C1_sce_subset_ord$Cell_type <- C1_labels_order$cell_type
C1_sce_subset_ord$Day <- C1_labels_order$day
C1_sce_subset_ord_r <- reshape2::melt(C1_sce_subset_ord, id.vars = c("Cell", "Day", "Cell_type"),
                                    measure.vars = features, variable.name = "Feat", value.name = "Expr")

pdf(file=paste0(outputdir, "/scdata_fifthten_celltypes.pdf"), height=10, width=8)
ggplot(C1_sce_subset_ord_r, aes(factor(Day), Expr)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE) +
  scale_y_continuous(expand = c(0, 0), position="right")+
  geom_jitter(aes(colour=Cell_type), alpha = 0.5,position = position_jitter(seed = 1, width = 0.2))+
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0))+
  labs(x="Day after onset of menstrual bleeding", y="Expression level", fill="Cell type")

dev.off()

c <- ggplot(C1_sce_subset_ord_r, aes(factor(Feat), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Day), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Feature on x-axis") + xlab("Feature") + ylab("Expression Level")
