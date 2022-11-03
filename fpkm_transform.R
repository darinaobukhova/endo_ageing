load("/Users/darinaobukhova/projects/Endometrial_Ageing/reanalysis_070622/corrected_data/filtered/expr_matrices_filtered.RData")
expr_matrix_filtered2 <- expr_matrix_filtered[rownames(expr_matrix_filtered) %in% gene_names_fc$ensembl_gene_id,]
sel_genes <- all_genes[all_genes@elementMetadata@listData$gene_id %in% rownames(expr_matrix_filtered2),]
idx <- match(sel_genes@elementMetadata@listData$gene_id, rownames(expr_matrix_filtered2))
sel_genes <- sel_genes[order(idx),]
expr_matrix_filtered_sel <- expr_matrix_filtered2[rownames(expr_matrix_filtered2) %in% sel_genes@elementMetadata@listData$gene_id,]
dds <- DESeqDataSetFromMatrix(expr_matrix_filtered_sel, samplesheet, design = ~Group)
rowRanges(dds) <- sel_genes
normalized_expr_matrix_filtered2 <- fpkm(dds)
gene_names_fc_sel <- gene_names_fc[gene_names_fc$ensembl_gene_id %in% rownames(normalized_expr_matrix_filtered2),]
idx <- match(gene_names_fc_sel$ensembl_gene_id, rownames(normalized_expr_matrix_filtered2))
gene_names_fc_sel <- gene_names_fc[order(idx),]
rownames(normalized_expr_matrix_filtered2) <- gene_names_fc_sel$hgnc_symbol
fpkm_matrix <- normalized_expr_matrix_filtered2

save(fpkm_matrix, file="/Users/darinaobukhova/projects/Endometrial_Ageing/reanalysis_070622/corrected_data/filtered/fpkm_matrix.RData")

 