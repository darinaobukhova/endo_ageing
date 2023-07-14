# Endometrial ageing 

Scripts and notebooks for the analyses used in our paper _Aging affects ciliated cells development in human endometrial epithelium_.

In this project, we compared the endometrial transcriptome of young and advanced age females undergoing hormonal replacement therapy (HRT)
before frozen embryo transfer, followed by immunohistology analysis and single-cell based deconvolution. RNA was extracted from the biopsies and 
and single-end sequenced them (75bp) on a NextSeq500 (Illumina).

The _align_count_workflow directory_ contains a snakemake workflow for mapping to the reference genome and counting, and the _scripts_ directory contains scripts used for data analysis and vizualisation.

# Data analysis and vizualization 

1. [Combining expression count matrix for each sample in one joint expression count matrix and processing:]

