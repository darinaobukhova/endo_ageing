#!/usr/bin/env Rscript

library(ggplot2)
library(ggsankey)

set.seed(15)

#Creating dataframe 


t1 <- c(rep("5", 5), rep("8", 5), rep("16", 5))
t2 <- c(rep("Full ROI", 15))
t3 <- c(rep("Blood_vessel", 3), rep("Follicles", 3), rep("Abnormal_follicles", 3), rep("Follicles_and_granulosa cells", 3), rep("Stromal_cells", 3))

t1 <- sample(x = c("5", "8", "16") , size = 15, replace=TRUE)
t2 <- sample(x = c("Full ROI"), size = 15, replace=TRUE)
t3 <- sample(x = c("Blood vessel", "Follicles", "Abnormal follicles", "Follicles_and_granulosa cells", "Stromal cells") , size = 15, replace=TRUE)

d <- data.frame(cbind(t1,t2,t3))
names(d) <- c('Patient', 'ROI',  'Cell type')

df <- d %>%
  make_long(Patient, ROI, Cell_type)
df

pdf(file=paste0(getwd(), "/ovarianST_sankey.pdf"), height=7, width=12)
ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node), 
               label = node)) +
  geom_sankey(flow.alpha = 0.75, node.color = "black", show.legend = FALSE) +
  geom_sankey_label(size=2.5, color = "black", fill="white", hjust=-0.07) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank())+
  scale_fill_viridis_d(option = "B", alpha = 0.7) +
  theme_sankey(base_size = 16) + labs(fill = 'Nodes')

dev.off()