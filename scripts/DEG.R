#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

#' Loading libraties

library(sva); library(biomaRt); library(DESeq2); library(pheatmap);
library(ggplot2); library(dplyr); library(RColorBrewer); library(rafalib); 

#' Loading files

load(corrected_data)

samplesheet <- read.csv2(samplesheet, header = T)
samplesheet <- samplesheet[order(match(samplesheet[,2], colnames(corrected_data_fc))), ]
print(all(samplesheet[,2] == colnames(corrected_data_fc)))

#' Filtering

corrected_data_fc_filtered <- corrected_data_fc[rowSums(corrected_data_fc[]) > 0, ]
#' Remove rows where values are >=5 in >75% samples
corrected_data_fc_filtered <- corrected_data_fc_filtered[rowSums(corrected_data_fc_filtered  >= 5) >= round(ncol(corrected_data_fc_filtered)/100 * 75),]

#' Performing DEG analysis

dds_fc <- DESeqDataSetFromMatrix(countData = as.matrix(corrected_data_fc_filtered), 
                                 colData = samplesheet, 
                                 design = ~Receptivity_Score + Group)

dds_fc$Group <- relevel(dds_fc$Group, ref="YMA")

dds <- DESeq(dds_fc,
             test = "Wald",
             fitType = "parametric",
             sfType = "ratio",
             modelMatrixType = "standard")

res_groups <- as.data.frame(results(dds,pAdjustMethod = "fdr"))

res_sort <- res_groups[order(res_groups$padj), ]

#' significant p-adjusted values
#' with different thresholds

res_p_sig <- res_sort %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange), desc(padj)) 

res_p_sig2 <- res_sort %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), desc(padj))

sig <- res_p_sig[res_p_sig$log2FoldChange > 1.0 | res_p_sig$log2FoldChange < -1.0, ]

#'Extract and save matrix of normalized counts

write.csv2(counts(dds, normalized=T), file=paste0(output_dir_deg, "/expr_matrix_normalized_ct.csv"))

#'Getting gene names

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


save(corrected_data_fc_filtered, gene_names_fc, file=paste0(cor_data_dir, "/expr_matrix_filt_allsamp.RData"))
write.csv2(res_p_sig, paste0(output_dir_deg, "/allsamp_allgenes_pval005_receptivityadj.csv"))
write.csv2(res_p_sig2, paste0(output_dir_deg, "/allsamp_allgenes_pval001_receptivityadj.csv"))
write.csv2(sig, paste0(output_dir_deg, "/allsamp_signDEGs_receptivityadj.csv"))
write.csv2(sig_upreg, file = paste0(output_dir_deg, "/allsamp_upreg_sig_DEGs_adj.csv"), quote = FALSE)
write.csv2(sig_downreg, file = paste0(output_dir_deg, "/allsamp_downreg_sig_DEGs_adj.csv"), quote = FALSE)

#' VST for vizualization

vst = varianceStabilizingTransformation(dds)
vsd=assay(vst)

save(vsd, file=paste0(output_dir_deg, "/allsamp_vst.RData"))

#' Trying different agglomeration methods for clustering on the matrix of all statistically significant genes
matrix <- vsd[rownames(res_p_sig),]
heat <- t(scale(t(matrix)))

#' Plotting
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

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
               scale = "none", clustering_method = "complete", 
               clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
               color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()

#' Trying different agglomeration methods for clustering on the matrix of DEGS 
#' logFC threshold of |1|

matrix <- vsd[rownames(sig),]
heat <- t(scale(t(matrix)))

#'Plotting

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

pheatmap(heat, annotation=annotation_col, annotation_colors = anno_colors, 
         scale = "none", clustering_method = "complete", 
         clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()

