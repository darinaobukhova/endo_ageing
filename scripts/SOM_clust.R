library(kohonen); library(reshape2)

set.seed(7)

samplesheet <- samplesheet[order(match(samplesheet[,2], colnames(corrected_data_fc))), ]
print(all(samplesheet[,2] == colnames(corrected_data_fc)))

#Remove rows with 0s only and where values are >=5 in >75% samples
corrected_data_fc_filtered <- corrected_data_fc[rowSums(corrected_data_fc[]) > 0, ]
corrected_data_fc_filtered <- corrected_data_fc_filtered[rowSums(corrected_data_fc_filtered  >= 5) >= round(ncol(corrected_data_fc_filtered)/100 * 75),]

#Scaling data and converting to matrix
som_data <- as.matrix(t(scale(t(corrected_data_fc_filtered))))
pcDat = prcomp(as.matrix(t(corrected_data_fc_filtered)), scale = T)
pca_scores <- data.frame(pcDat$x)

#Dim data
dim <- dim(som_data) |> dplyr::glimpse()

#Creating SOM (Vesanto(2005) suggests optimal size as 5*sqrt(k))
#Best results achieved with Gaussian kernel 

xdim = 5 * round(sqrt(dim[1] * dim[2]))/2
ydim = xdim

som_map <- som_data |> 
  kohonen::som(grid=somgrid(4, 6, "hexagonal", "gaussian"))

plot(som_map, type ="mapping")

som_cluster <- cutree(hclust(dist(som_map$codes[[1]])), 2)

plot(som_map, type="mapping", bgcol = som_cluster, main = "Clusters") 
add.cluster.boundaries(som_map, som_cluster) 

som_clusterKey <- data.frame(som_cluster)
som_clusterKey$unit.classif <- c(1:24)

theme_set(theme_bw())

ggplot(pca_scores, aes(PC1, PC2, colour=factor(samplesheet[,3]), shape=factor(som_cluster))) +
  ggrepel::geom_text_repel(aes(label = rownames(pcDat$x))) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_blank()) +
  scale_color_manual(values=ifelse(samplesheet[,3] == "YMA", "darkorchid4", "darkorange2")) +
  guides(color = guide_legend(title="Group",override.aes = list(size=4, linetype=0), shape = guide_legend(title="SOM cluster"))) +
  stat_ellipse(level=0.95, linetype=2)

annotation_col = data.frame(Group = samplesheet[,3]) 
rownames(annotation_col) <- colnames(corrected_data_fc)
Group <- c("darkorchid4", "darkorange2")
names(Group) <- c("YMA", "AMA")
anno_colors <- list(Group = Group)

pheatmap(som_data, annotation=annotation_col, annotation_colors = anno_colors, 
         scale = "none", clustering_method = "average", 
         clusering_distance_cols = "euclidean", show_rownames = F, show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red")) ((50)), border_color = F)

dev.off()


