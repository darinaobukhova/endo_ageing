#!/usr/bin/env Rscript
options(scipen=999)
set.seed(8)

#Arguments
args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

cat("Loading libraries...")
library(dplyr);library(ggplot2);library(ggpubr);library(data.table);library(DWLS);library(scales)

####################
####Loading data####
####################

cat("Loading data...")
load(data_for_deconvolution)
load(signature_tenx)
samplesheet <- read.csv2(samplesheet, header = T)

#####################
####Deconvolution####
#####################

dir.create(file.path(outputdir, "results", gsub("-", "_",Sys.Date())))
setwd(file.path(outputdir, "results", gsub("-", "_",Sys.Date())))

#For C1 dataset
cat("Deconvolving cellular types in the bulk dataset using 10x dataset as the reference...")

#Using 2 methods of DWLS package: dampened weighted least squares and a support vector regression

#Raw data

allProport_DWLS <- NULL
allProport_SVR<- NULL

for(j in 1:(ncol(corrected_data_fc_filtered))){
  sign <- signature_tenx_midsec
  bulk <- corrected_data_fc_filtered[,j]
  names(bulk) <- rownames(corrected_data_fc_filtered)
  genes <- intersect(rownames(sign), names(bulk))
  B <- bulk[genes]
  S <- sign[genes,]
  proport_DWLS <- solveDampenedWLS(S,B)
  proport_SVR  <- solveSVR(S,B)
  
  allProport_DWLS <- cbind(allProport_DWLS, proport_DWLS)
  allProport_SVR <- cbind(allProport_SVR, proport_SVR)
}

#Modifying for plotting purposes later on

allProport_SVR <- as.data.frame(t(allProport_SVR))
rownames(allProport_SVR) <- samplesheet$SampleID
write.csv(allProport_SVR, file="cellprop_tenx_midsec_svr_raw.csv")
allProport_SVR$Group <- samplesheet$Group
allProport_SVR <- allProport_SVR[order(allProport_SVR$Group, decreasing=T),]
allProport_DWLS <- as.data.frame(t(allProport_DWLS))
rownames(allProport_DWLS) <- samplesheet$SampleID
write.csv2(allProport_DWLS, file="cellprop_tenx_midsec_dwls_raw.csv")
allProport_DWLS$Group <- samplesheet$Group
allProport_DWLS <- allProport_DWLS[order(allProport_DWLS$Group, decreasing=T),]

save(allProport_SVR, allProport_DWLS, file = "deconv_res_tenxref.RData")

cat("Performing statistics tests...")

comp_dwls_c1_raw <- apply(allProport_DWLS[, -ncol(allProport_DWLS)], 2, function(x) wilcox.test(x[1:12], x[12:24], exact = FALSE, 
                                                                                                paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()

comp_svr_c1_raw <- apply(allProport_SVR[, -ncol(allProport_SVR)], 2, function(x) wilcox.test(x[1:12], x[12:24], exact = FALSE, 
                                                                                             paired=FALSE, p.adjust = "bonferroni")$p.value) %>% as.data.frame()

write.csv(comp_dwls_c1_raw, file="prop_DWLS_tenx_raw.csv")
write.csv(comp_svr_c1_raw, file="prop_svr_tenx_raw.csv")

allProport_DWLS2 <- round(allProport_DWLS[,-ncol(allProport_DWLS)], digits = 5)

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

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_tenx_midsec_DWLS_raw.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_c1_log_dwls, ncol = 3)

dev.off()

allProport_DWLS2$Patient <- factor(rownames(allProport_DWLS2), levels=unique(rownames(allProport_DWLS2)))
allProport_DWLS$Age <- c("26", "23", "25", "26", "27", "26", "23", "26", "25", "24", "24", "27", 
                         "47", "47", "48", "50", "48", "48", "47", "47", "49", "47", "50", "47")
allProport_DWLS2 <- allProport_DWLS[order(allProport_DWLS$Age),]

pdf(file=paste0(outputdir, "/tenx_DWLS_prop_ageordered.pdf"), height=7, width=9)
theme_set(theme_bw())
ggplot(allProport_DWLS2, aes(x = Patient,y=`Ciliated`, color = as.factor(Group))) +
  geom_point(size = 2.5, color=ifelse(allProport_DWLS2$Group == "YMA", "darkorchid4", "darkorange2")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank())
dev.off()

allProport_DWLS_c1_bc <- gather(allProport_DWLS, 1:6, key="Cell_type", value="Proportion")
allProport_DWLS_c1_bc$Group <- factor(allProport_DWLS_c1_bc$Group, levels=c("YMA", "AMA"))

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_tenx_midsec_DWLS_raw_cellprop.pdf"), height=7, width=7)
ggplot(allProport_DWLS_c1_bc, aes(fill=Cell_type, y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "dodgerblue3", "#8B008B", "gray27"))+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


