#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

library(ggplot2); library(dplyr); library(RColorBrewer);
library(rafalib); library(sva); library(biomaRt); library(DESeq2); library(pheatmap); 

#########################
####Loading files########
#########################

load(corrected_data_fc)

samplesheet <- read.csv2(samplesheet, header = T)

samplesheet <- samplesheet[order(match(samplesheet[,2], colnames(corrected_data_fc))), ]
print(all(samplesheet[,2] == colnames(corrected_data_fc)))

#####################
####Filtering########
#####################

corrected_data_fc_filtered <- corrected_data_fc[rowSums(corrected_data_fc[]) > 0, ]
#Remove rows where values are >=5 in >75% samples
corrected_data_fc_filtered <- corrected_data_fc_filtered[rowSums(corrected_data_fc_filtered  >= 5) >= round(ncol(corrected_data_fc_filtered)/100 * 75),]

dds_fc <- DESeqDataSetFromMatrix(countData = as.matrix(corrected_data_fc_filtered), 
                                colData = samplesheet, 
                                design = ~Receptivity_Score + Group)

dds_fc$Group <- relevel(dds_fc$Group, ref="YMA")

dds <- DESeq(dds_fc)

res_groups <- as.data.frame(results(dds,pAdjustMethod = "fdr"))

res_sort <- res_groups[order(res_groups$padj), ]

res_p_sig <- res_sort %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange), desc(padj)) #significant p-adj values 

res_p_sig2 <- res_sort %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), desc(padj))

sig <- res_p_sig[res_p_sig$log2FoldChange > 1.0 | res_p_sig$log2FoldChange < -1.0, ]

write.csv2(counts(dds, normalized=T), file=paste0(output_dir_deg, "/expr_matrix_normalized_ct.csv"))

##############################
####Getting gene names########
#############################

mart <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')

getinfo <- c("ensembl_gene_id", "hgnc_symbol")

gene_names_fc <- biomaRt::getBM(attributes=getinfo,
                                filters="ensembl_gene_id", 
                                values=names(corrected_data_fc_filtered), 
                                mart=mart)

res_p_sig <- cbind(gene_s = gene_names_fc$hgnc_symbol[match(rownames(res_p_sig), 
                                                      gene_names_fc$refseq_mrna)],res_p_sig)

sig <- cbind(gene_s = gene_names_fc$hgnc_symbol[match(rownames(sig), 
                                                gene_names_fc$refseq_mrna)], sig)

sig_upreg <- sig[sig$log2FoldChange > 1.0, ]

sig_downreg <- sig[sig$log2FoldChange < -1.0, ]

dir.create(output_dir_deg)


save(corrected_data_fc_filtered, file=paste0(cor_data_dir, "/expr_matrix_filt_allsamp.RData"))
save(gene_names_fc, file=paste0(output_dir_deg, "/datProbes.RData"))
write.csv2(res_p_sig, paste0(output_dir_deg, "/allsamp_allgenes_pval005_receptivityadj.csv"))
write.csv2(res_p_sig2, paste0(output_dir_deg, "/allsamp_allgenes_pval001_receptivityadj.csv"))
write.csv2(sig, paste0(output_dir_deg, "/allsamp_signDEGs_receptivityadj.csv"))
write.csv2(sig_upreg, file = paste0(output_dir_deg, "/allsamp_upreg_sig_DEGs_adj.csv"), quote = FALSE)
write.csv2(sig_downreg, file = paste0(output_dir_deg, "/allsamp_downreg_sig_DEGs_adj.csv"), quote = FALSE)

vst = varianceStabilizingTransformation(dds)
vsd=assay(vst)

save(vsd, file=paste0(output_dir_deg, "/allsamp_vst.RData"))

matrix <- vsd[rownames(res_p_sig),]

heat <- t(scale(t(matrix)))

#Trying different agglomeration methods for clustering

pdf(file=paste0(output_dir_deg, "/allsamp_heatmap_for_all_sign_genes_fc_receptadj_average.pdf"), height=7, width=9)

annotation_col = data.frame(Group = samplesheet[,3]) 
rownames(annotation_col) <- colnames(vsd)
Group <- c("darkorchid4", "darkorange2")
names(Group) <- c("YMA", "AMA")
anno_colors <- list(Group = Group)

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
               scale = "none", clustering_method = "average", 
               clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
               color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()

pdf(file=paste0(output_dir_deg, "/allsamp_heatmap_for_all_sign_genes_fc_receptadj_complete.pdf"), height=7, width=9)

annotation_col = data.frame(Group = samplesheet[,3]) 
rownames(annotation_col) <- colnames(vsd)
Group <- c("darkorchid4", "darkorange2")
names(Group) <- c("YMA", "AMA")
anno_colors <- list(Group = Group)

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
               scale = "none", clustering_method = "complete", 
               clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
               color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()


matrix <- vsd[rownames(sig),]

heat <- t(scale(t(matrix)))


pdf(file=paste0(output_dir_deg, "/allsamp_heatmap_for_DEGs_fc_receptadj_average.pdf"), height=7, width=9)

annotation_col = data.frame(Group = samplesheet[,3]) 
rownames(annotation_col) <- colnames(vsd)


Group <- c("darkorchid4", "darkorange2")
names(Group) <- c("YMA", "AMA")
anno_colors <- list(Group = Group)

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
         scale = "none", clustering_method = "average", 
         clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()

pdf(file=paste0(output_dir_deg, "/allsamp_heatmap_for_DEGs_fc_receptadj_complete.pdf"), height=7, width=9)

annotation_col = data.frame(Group = samplesheet[,3]) 
rownames(annotation_col) <- colnames(vsd)


Group <- c("darkorchid4", "darkorange2")
names(Group) <- c("YMA", "AMA")
anno_colors <- list(Group = Group)

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
         scale = "none", clustering_method = "complete", 
         clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()

top10 <- vsd[rownames(head(sig, 10)), ]
bottom10 <- vsd[rownames(tail(sig, 10)), ]

top20 <- rbind(top10, bottom10)
write.csv2(rownames(top20), paste0(output_dir_deg, "/83excl_top_20_DEGs_receptivityadj.csv"))

gene_s <- gene_names_fc$hgnc_symbol[match(rownames(top20), gene_names_fc$ensembl_gene_id)]
rownames(top20) <- gene_s

heat2 <- t(scale(t(top20)))

pdf(file=paste0(output_dir_deg, "/allsamp_heatmap_for_20_top_degs_receptivityadj_complete.pdf"), height=7, width=9)

pheatmap(heat2, annotation_col=annotation_col, annotation_colors = anno_colors, 
         scale = "none", clustering_method = "complete", 
         clusering_distance_cols = "euclidean", show_rownames = T, show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()


