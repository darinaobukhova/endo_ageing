# Endometrial ageing 

Scripts and notebooks for the analyses used in our paper _Aging affects ciliated cells development in human endometrial epithelium_.

In this project, we compared the endometrial transcriptome of young and advanced age females undergoing hormonal replacement therapy (HRT)
before frozen embryo transfer, followed by immunohistology analysis and single-cell based deconvolution. RNA was extracted from the biopsies and 
and single-end sequenced them (75bp) on a NextSeq500 (Illumina).

The _align_count_workflow directory_ contains a snakemake workflow for mapping to the reference genome and counting, the _scripts_ directory contains scripts used for data analysis and 
_figures_ directory contains scripts and notebooks used for vizualisation.

# Data analysis and vizualization 

1. [Combining expression count matrix for each sample in one joint expression count matrix and processing:](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/01-count_matrix_preprocessing.R)
    - removal of rows with all 0s and values >=5% in >75% samples
    - retrieval of gene names in the place of Ensembl identifiers
    - correcting for sequencing run using Combat_seq() from sva package
2. [Principal component analysis](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/02-PCA.R) and [clustering analysis](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/02-sample_clustering.R)
3. [Differential gene expression analysis using DESeq2 package](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/03-DEG_analysis.R)
4. [Gene ontology analysis](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/04-GO_analysis.R)
5. [Preparation of single cell datasets for deconvolution analysis with DWLS method](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/05-deconvolution_preparation.R)
6. [Deconvolution analysis with DWLS method](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/06-deconvolution.R)
7. [Normalization of single-cell dataset for later visualization](https://github.com/darinaobukhova/endo_ageing/blob/main/scripts/07-single_cell_dataset_normalization.R)
