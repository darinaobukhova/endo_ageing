#Garzia-Alonso dataset characteristics

library(rcartocolor)

secretory_donorE1 <- read.table("/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/Garzia_Alonso_10x/annotations/MRC_Endo8625698_cells.tsv",
                                header=T, sep="\t")
secretory_donorE1_extra <- read.table("/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/Garzia_Alonso_10x/annotations/MRC_Endo8625699_cells.tsv",
                                      header=T, sep="\t")
secretory_donorE2 <- read.table("/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/Garzia_Alonso_10x/annotations/MRC_Endo8712024_cells.tsv",
                                header=T, sep="\t")
secretory_donorE2_extra <- read.table("/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/Garzia_Alonso_10x/annotations/MRC_Endo8712032_cells.tsv",
                                      header=T, sep="\t")


head(secretory_donorE1, 10)
head(secretory_donorE1_extra, 10)
head(secretory_donorE2, 10)
head(secretory_donorE2_extra, 10)

#Filter out the cells with "Other" and "filtered" labels in the first part of dataset
secretory_donorE1 <- as.data.frame(secretory_donorE1)
secretory_donorE1_subset <- secretory_donorE1[!(secretory_donorE1$general_celltypes == "Other" | secretory_donorE1$general_celltypes == "filtered"),]
secretory_donorE1_subset$general_celltypes <- sub(' ', '_', secretory_donorE1_subset$general_celltypes)

#Filter out the cells with "Other" and "filtered" labels in the second part of dataset
secretory_donorE1_extra <- as.data.frame(secretory_donorE1_extra)
secretory_donorE1_extra_subset <- secretory_donorE1_extra[!(secretory_donorE1_extra$general_celltypes == "Other" | secretory_donorE1_extra$general_celltypes == "filtered"),]
secretory_donorE1_extra_subset$general_celltypes <- sub(' ', '_', secretory_donorE1_extra_subset$general_celltypes)
full_sec_data_labels <- rbind(secretory_donorE1_subset, secretory_donorE1_extra_subset)

full_sec_data_labels_df <- data.frame(table(full_sec_data_labels$general_celltypes))
full_sec_data_labels_df$Day <- c("25")

full_sec_E1_sec_prop <- full_sec_data_labels_df %>% 
  pivot_wider(names_from=Var1, values_from = Freq) %>% 
  mutate(total = sum(c_across(Arterial_Endothelia:Venular_Endothelia))) %>%
  ungroup() %>% 
  mutate(across(Arterial_Endothelia:Venular_Endothelia, ~ ./total))

write.table(full_sec_E1_sec_prop, file=paste0(outputdir, "/E1_latesec_prop.csv"), sep=",", row.names = F)

#Filter out the cells with "Other" and "filtered" labels in the first part of dataset
secretory_donorE2 <- as.data.frame(secretory_donorE2)
secretory_donorE2_subset <- secretory_donorE2[!(secretory_donorE2$general_celltypes == "Other" | secretory_donorE2$general_celltypes == "filtered"),]
secretory_donorE2_subset$general_celltypes <- sub(' ', '_', secretory_donorE2_subset$general_celltypes)

#Filter out the cells with "Other" and "filtered" labels in the second part of dataset
secretory_donorE2_extra <- as.data.frame(secretory_donorE2_extra)
secretory_donorE2_extra_subset <- secretory_donorE2_extra[!(secretory_donorE2_extra$general_celltypes == "Other" | secretory_donorE2_extra$general_celltypes == "filtered"),]
secretory_donorE2_extra_subset$general_celltypes <- sub(' ', '_', secretory_donorE2_extra_subset$general_celltypes)
full_sec_E2_data_labels <- rbind(secretory_donorE2_subset, secretory_donorE2_extra_subset)

full_sec_E2_data_labels_df <- data.frame(table(full_sec_E2_data_labels$general_celltypes))
full_sec_E2_data_labels_df$Day <- c("20")

full_sec_E2_sec_prop <- full_sec_E2_data_labels_df %>% 
  pivot_wider(names_from=Var1, values_from = Freq) %>% 
  mutate(total = sum(c_across(Arterial_Endothelia:Venular_Endothelia))) %>%
  ungroup() %>% 
  mutate(across(Arterial_Endothelia:Venular_Endothelia, ~ ./total))

write.table(full_sec_E2_sec_prop, file=paste0(outputdir, "/E2_latesec_prop.csv"), sep=",", row.names = F)

safe_pal <- carto_pal(12, "Safe")
pdf(file=paste0(outputdir, "/E1_Day25_prop.pdf"), height=8, width=8)
ggplot(full_sec_data_labels_df, aes(fill=Var1, y=Freq, x=as.factor(Day))) + 
  geom_col(width = 0.8) + 
  theme_bw() +
  labs(x = "Day", y="Number of cells", fill="Cell type") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=14)) +
  scale_fill_manual(values=safe_pal)
dev.off()

pdf(file=paste0(outputdir, "/E2_Day20_prop.pdf"), height=8, width=8)
ggplot(full_sec_E2_data_labels_df, aes(fill=Var1, y=Freq, x=as.factor(Day))) + 
  geom_col(width = 0.8) + 
  theme_bw() +
  labs(x = "Day", y="Number of cells", fill="Cell type") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=14)) +
  scale_fill_manual(values=safe_pal)
dev.off()


