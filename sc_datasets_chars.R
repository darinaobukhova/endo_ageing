
C1_labels_n <- data.frame(table(C1_midsec_labels$day, C1_midsec_labels$cell_type))
C1_labels_midsec_prop <- C1_labels_n %>% 
  pivot_wider(names_from=Var2, values_from = Freq) %>% 
  rowwise(Var1) %>% 
  mutate(total = sum(c_across(Ciliated:Unciliated_epithelia))) %>%
  ungroup() %>% 
  mutate(across(Ciliated:Unciliated_epithelia, ~ ./total))


pdf(file="/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/C1/secretoryphase_prop.pdf", height=6, width=8)
ggplot(C1_labels_n, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_col() + 
  theme_bw() +
  labs(title="C1 dataset", x="Day after onset of menstrual bleeding", y="Number of cells", fill="Cell type") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=14)) + 
  scale_fill_aaas()
dev.off()

tenx_labels_n <- data.frame(table(tenx_midsec_labels$day, tenx_midsec_labels$cell_type))
tenx_labels_midsec_prop <- tenx_labels_n %>% 
  pivot_wider(names_from=Var2, values_from = Freq) %>% 
  rowwise(Var1) %>% 
  mutate(total = sum(c_across(Ciliated:Unciliated_epithelia_2))) %>%
  ungroup() %>% 
  mutate(across(Ciliated:Unciliated_epithelia_2, ~ ./total))

pdf(file="/Users/darinaobukhova/projects/Endometrial_Ageing/deconvolution/data/10x/secretoryphase_prop.pdf", height=6, width=8)
ggplot(tenx_labels_n, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_col() + 
  theme_bw() +
  labs(title="tenx dataset", x="Day after onset of menstrual bleeding", y="Number of cells", fill="Cell type") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=14)) + 
  scale_fill_aaas()
dev.off()