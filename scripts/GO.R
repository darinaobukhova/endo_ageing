#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
config <- args[1]
source(as.character(config))

#########################
####Loading libraries####
#########################

library(RColorBrewer); library(enrichplot); library(clusterProfiler.dplyr); library(DOSE); 
library(gprofiler2); library(ggplot2)

load(paste0(output_dir_deg, "/datProbes.RData"))

#Enrichment analysis of DEGs


GO_res <- gost(list(gene_names_fc$hgnc_symbol),
               multi_query = F, evcodes = T)

GO_res_enr = GO_res$result[,c("query", "source", "term_id",
                              "term_name", "p_value", "query_size", 
                              "intersection_size", "term_size", 
                              "effective_domain_size", "intersection")]

GO_res_enr$GeneRatio = paste0(GO_res_enr$intersection_size,  "/", GO_res_enr$query_size)

GO_res_enr$BgRatio = paste0(GO_res_enr$term_size, "/", GO_res_enr$effective_domain_size)

names(GO_res_enr) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                      "query_size", "Count", "term_size", "effective_domain_size", 
                      "geneID", "GeneRatio", "BgRatio")

GO_res_enr$geneID = gsub(",", "/", GO_res_enr$geneID)


GO_res_enr_ord = GO_res_enr[order(GO_res_enr$p.adjust), ] %>% 
  filter(p.adjust <= 0.05) %>% 
  filter(Count > 1) %>% 
  arrange(desc(Count), p.adjust)

GO_res_enr_ord$Description <- factor(GO_res_enr_ord$Description, 
                                     levels = unique(GO_res_enr_ord$Description[order(GO_res_enr_ord$Count, GO_res_enr_ord$p.adjust)]))

write.csv(GO_res_enr_ord, file = paste0(outputdir, "/allgenes_GO.csv"), quote = FALSE)

pdf(file=paste0(output_dir_deg, "/", Sys.Date(), "GO_enrichment.pdf"), width = 16, height=9)


ggplot(data = GO_res_enr_ord, aes(x=Cluster, y = Description,
                                  color = p.adjust, size = Count)) +
  geom_point() + 
  scale_size_continuous(range = c(1,5)) +
  scale_colour_gradient(low = "darkblue", high="darkred") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14, family = "Helvetica"),
        axis.text.x = element_text(size = 14),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("") + 
  xlab("") + 
  guides(fill=guide_legend("q-value")) +
  ggtitle("GO enrichment analysis") 

dev.off()
