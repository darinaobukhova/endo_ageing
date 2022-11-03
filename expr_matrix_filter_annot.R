#!/usr/bin/env Rscript

options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)

#Arguments
config <- args[1]
source(as.character(config))

#Loading libraties
library(biomaRt); library(dplyr); library(GenomicRanges); library(GenomicFeatures)

#Loading data
load(corrected_data_fc)

#################
####Filtering####
#################

cat("Filtering the count matrix...")
expr_matrix_filtered <- corrected_data_fc[rowSums(corrected_data_fc[]) > 0, ]
expr_matrix_filtered <- expr_matrix_filtered[rowSums(expr_matrix_filtered >= 5) >= round(ncol(expr_matrix_filtered)/100 * 75),]
sprintf("%i rows were removed and filtered matrix has now %i rows", nrow(corrected_data_fc) - nrow(expr_matrix_filtered), nrow(expr_matrix_filtered))


##################
####Annotating####
#################

cat("Annotating the count matrix...")

#Gene names 

mart <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', mirror = "useast")

getinfo <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", 
             "end_position", "strand", "gene_biotype", "percentage_gene_gc_content", "transcript_start",
             "transcript_end")

gene_names_fc <- biomaRt::getBM(attributes=getinfo,
                                filters="ensembl_gene_id", 
                                values=rownames(expr_matrix_filtered), 
                                mart=mart)

#Transcript length
gene_names_fc <- gene_names_fc[!duplicated(gene_names_fc[, c("ensembl_gene_id")]),]
gene_names_fc$transcript_length <- gene_names_fc$transcript_end - gene_names_fc$transcript_start
save(gene_names_fc, file = paste0(outputdir, "probes.RData"))

#GRangesListObject with transcript length for later use
len_for_fpkm <- data.frame(chr=gene_names_fc$chromosome_name, start=gene_names_fc$transcript_start, 
                           end=gene_names_fc$transcript_end, gene_symbol = gene_names_fc$hgnc_symbol)
len_for_fpkm_gr <- makeGRangesListFromDataFrame(len_for_fpkm, names.field = 'gene_symbol')

#Matching the order
idx <- match(gene_names_fc$ensembl_gene_id, rownames(expr_matrix_filtered))
gene_names_fc <- gene_names_fc[order(idx),]

#Excluding unannotated genes and filtering the matrix

unannot <- rownames(expr_matrix_filtered)[!rownames(expr_matrix_filtered) %in% gene_names_fc$ensembl_gene_id]
sprintf("There are %i genes which were not annotated with biomart", length(unannot))
expr_matrix_filtered2 <- expr_matrix_filtered[rownames(expr_matrix_filtered) %in% gene_names_fc$ensembl_gene_id,]

#Gene length
txdb <- makeTxDbFromGFF(gtf_file)
all_genes <- genes(txdb)
pres_genes <- all_genes[all_genes@ranges@NAMES %in% rownames(expr_matrix_filtered2),]

pres_genes_lengths <- width(pres_genes)
names(pres_genes_lengths) <- pres_genes@ranges@NAMES
pres_genes_lengths <- as.data.frame(pres_genes_lengths)
colnames(pres_genes_lengths) <- c("gene_length")
gene_names_pres <- gene_names_fc[gene_names_fc$ensembl_gene_id %in% rownames(pres_genes_lengths),]
pres_genes_lengths$gene_Sym <- gene_names_pres$hgnc_symbol

#rownames(expr_matrix_filtered2) <- gene_names_fc$hgnc_symbol
sprintf("Annotated matrix has %i rows", nrow(expr_matrix_filtered2))
#Saving intermediates
save(expr_matrix_filtered, expr_matrix_filtered2, gene_names_fc, all_genes, pres_genes, len_for_fpkm, len_for_fpkm_gr, pres_genes_lengths, file = paste0(outputdir, "expr_matrices_filtered.RData"))
sprintf("The filtered matrices are stored at %s", outputdir)
