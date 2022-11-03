#!/usr/bin/env Rscript

options(scipen=999)

#Arguments
args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

library(dplyr); library(biomaRt); library(DESeq2); library(ggplot2); library(DWLS); library(ggpubr);
library(Seurat);

dir.create(file.path(outputdir, "results", gsub("-", "_",Sys.Date())), showWarnings = F)
setwd(file.path(outputdir, "results", gsub("-", "_",Sys.Date())))

####################
####Loading data####
####################

#Bulk data
cat("Loading bulk dataset...")
load(expr_matrix_filtered)
samplesheet <- read.csv2(samplesheet, header = T)

#Single cell data
cat("Loading single cell datasets (C1 and 10x)...")
C1 <- read.csv(C1_data, header = T)
C1_labels <- read.csv(C1_labels, header = T)
tenx <- readRDS(tenx_data)
tenx_labels <- read.csv(tenx_labels, header = T)

#######################
####Data processing####
#######################

cat("Pre-processing bulk data...")

#Putting conditions instead of sample names

#samplesheet <- samplesheet[!rownames(samplesheet) %in% c("RNA83"), ]
idx <- match(colnames(corrected_data_fc_filtered), samplesheet$SampleID)
corrected_data_fc_filtered <- corrected_data_fc_filtered[, order(idx)]
#Check the equivalence
print(all(samplesheet$SampleID == colnames(corrected_data_fc_filtered)))
colnames(corrected_data_fc_filtered) <- samplesheet$Group

#Naming genes
mart <- biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')

getinfo <- c("ensembl_gene_id", "hgnc_symbol")

gene_names <- biomaRt::getBM(attributes=getinfo,
                                filters="ensembl_gene_id", 
                                values=rownames(corrected_data_fc_filtered), 
                                mart=mart)

corrected_data_fc_filtered <- corrected_data_fc_filtered[rownames(corrected_data_fc_filtered) %in% gene_names$ensembl_gene_id,]
idx2 <- match(rownames(corrected_data_fc_filtered), gene_names$ensembl_gene_id)
corrected_data_fc_filtered <- corrected_data_fc_filtered[order(idx2),]
print(all(rownames(corrected_data_fc_filtered) == gene_names$ensembl_gene_id))
rownames(corrected_data_fc_filtered) <- gene_names$hgnc_symbol

#Important to have data object as a double type for later processing

#Order columns alphabetically
corrected_data_fc_filtered <- corrected_data_fc_filtered[, order(colnames(corrected_data_fc_filtered), decreasing = T)]
nms_bulk <- rownames(corrected_data_fc_filtered)
corrected_data_fc_filtered_db <- as.double(corrected_data_fc_filtered)
names(corrected_data_fc_filtered_db) <- nms_bulk

cat("Finished pre-processing bulk data...")

#Filtering of single cell (C1 data)
cat("Pre-processing single cell dataset (C1)")

#Adjusting the column names
colnames(C1) <- sub("^X", "", colnames(C1))
rownames(C1) <- C1[, 1]
C1 <- C1[, -1]

#Subsetting the mid-secretory cells
C1_midsec_labels <- C1_labels[C1_labels$day == '21' | C1_labels$day == '22' | C1_labels$day == '24', ]
C1_midsec <- C1[,colnames(C1) %in% C1_midsec_labels$cell_name]

#Gene level filtering
#Remove genes with 0 counts in all the cells
C1_midsec <- C1_midsec[rowSums(C1[]) > 0, ]
sprintf("%i genes with 0 counts across all the cells were removed", nrow(C1) - nrow(C1_midsec))
sprintf("Now the C1 expression matrix has %i genes", nrow(C1_midsec))

#Keep genes that are 'detectable' (have at least 10 counts) in at least 5% of cells
keep_genes <- rowSums(C1_midsec >= 0) >= round(0.05 * ncol(C1_midsec))

C1_midsec_filt <- C1_midsec[keep_genes,]
sprintf("%i genes expressed in less than 10 cells were removed", nrow(C1_midsec) - nrow(C1_midsec_filt))
sprintf("Now the C1 expression matrix has %i genes", nrow(C1_midsec_filt))

#Collapsing cell labels for ciliated and unciliated epithelia in "Epithelial cells" type
#C1_labels[C1_labels$cell_type == "Ciliated" | C1_labels$cell_type == "Unciliated epithelia",]$cell_type <- c("Epithelial cells")

#Some editing of the labels
C1_midsec_labels$cell_type <- sub(' ', '_', C1_midsec_labels$cell_type)

#Filtering of single cell (10X data)

cat("Pre-processing single cell dataset (10x)...")

seurt <- Seurat::CreateSeuratObject(counts = tenx, min.cells = round(0.05 * length(tenx@Dimnames[[2]])), min.features = 350)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurt@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(seurt@assays$RNA@data[mito.genes, ])/Matrix::colSums(seurt@assays$RNA@data)
seurt[["percent.mito"]] <- PercentageFeatureSet(seurt, pattern = "^MT-")

seurt_filt <- subset(seurt, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mito < 0.15)
seurt_filt <- NormalizeData(object = seurt, normalization.method = "LogNormalize")

#Collapsing cell labels for ciliated and unciliated epithelia 1,2 in "Epithelial cells" type
#tenx_labels[tenx_labels$cell_type == "Ciliated" | tenx_labels$cell_type == "Unciliated epithelia 1" | tenx_labels$cell_type == "Unciliated epithelia 2",]$cell_type <- c("Epithelial cells")

#Some editing of the labels
tenx_labels$cell_type <- gsub(' ', '_', tenx_labels$cell_type)

#Subsetting the mid-secretory cells
tenx_midsec_labels <- tenx_labels[tenx_labels$day == '20' | tenx_labels$day == '22' | tenx_labels$day == '23', ]
seurt_midsec <- seurt_filt[,colnames(seurt_filt) %in% tenx_midsec_labels$cell_name]

save(corrected_data_fc_filtered, corrected_data_fc_filtered_db, C1_midsec_filt, C1_midsec_labels, seurt_midsec, tenx_midsec_labels, file = "data_for_deconv.RData")

#Build signature from single-cell data (C1)

cat("Building signature matrix from C1 data...")

#signature_C1_midsec <- buildSignatureMatrixMAST(scdata = C1_midsec_filt, id = C1_midsec_labels$cell_type, path = getwd(),
                                                #diff.cutoff = 0.5, pval.cutoff = 0.01)

#save(signature_C1_midsec, file= "signature_C1_midsec.RData")

#sprintf("Signature matrix from C1 data has %i genes", dim(signature_C1_midsec)[1])

#Build signature from single-cell data (10X)

cat("Building signature matrix from 10x data...")

signature_tenx_midsec <- buildSignatureMatrixMAST(scdata = as.matrix(seurt_midsec@assays$RNA@counts), 
                                                  id = tenx_midsec_labels[tenx_midsec_labels$X %in% colnames(seurt_midsec@assays$RNA@counts),]$cell_type,
                                                  path = getwd(), diff.cutoff = 0.5, pval.cutoff = 0.01)

save(signature_tenx_midsec, file = "signature_tenx_midsec.RData")

sprintf("Signature matrix from 10x data has %i genes", dim(signature_tenx_midsec)[1])
