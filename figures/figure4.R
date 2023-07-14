#########################################################################################################################
######################################################## FIGURE 4 #######################################################
#########################################################################################################################

options(scipen=999)
set.seed(8)

# Arguments

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
source(as.character(config))

# Load data

load(deconvolved_data)
samplesheet <- read.csv2(samplesheet, header = T)

if (!dir.exists(outputdir))
    {dir.create(outputdir)}

# Plotting violin plots for C1 dataset (a)

palette <- c("gray27", "mediumseagreen", "darkmagenta", "indianred3", "mediumslateblue", "dodgerblue3")

deconv_plots_c1_dwls <- lapply(colnames(allProport_DWLS_C1[1:6]), function(var) {
  idx <- match(var, colnames(allProport_DWLS_C1))
  ggpar(ggviolin(allProport_DWLS_C1, x = "Group", y = var, fill = palette[idx], 
                 palette = palette, add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.8, width = 0.8) +
          geom_jitter(alpha = 0.4, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                plot.title = element_text(face="bold"), axis.text = element_text(size=12),
                axis.title = element_text(size=12))
  )
})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_C1_midsec_DWLS_raw.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_c1_dwls, ncol = 3)

dev.off()

# Plotting violin plots for 10x dataset (c)

palette2 <- c("gray35", "darkmagenta", "gray50", "indianred3", "mediumslateblue", "green4", "mediumseagreen", "dodgerblue3")

deconv_plots_10x_dwls <- lapply(colnames(allProport_DWLS_10x[, -ncol(allProport_DWLS_10x)]), function(var) {
  idx <- match(var, colnames(allProport_DWLS_10x))
  ggpar(ggviolin(allProport_DWLS_10x, x = "Group", y = var, fill = palette2[idx], 
                 palette = palette2, add = "boxplot", add.params = list(fill = "white"),
                 xlab = FALSE, ylab = FALSE, alpha = 0.8, width = 0.8,
                 scale = "width") +
          geom_jitter(alpha = 0.4, position = position_jitter(seed = 1, width = 0.2)) + 
          scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                axis.text.x = element_blank(), axis.title.x = element_blank())
  )
})

theme_set(theme_bw())

pdf(file=paste0(outputdir, "/EndoAgeing_deconv_10x_midsec_DWLS_raw.pdf"), height=10, width=15)

gridExtra::grid.arrange(grobs = deconv_plots_10x_dwls, ncol = 3)

dev.off()

# Barcharts (b,d)

allProport_DWLS_c1_bc <- gather(allProport_DWLS_C1, 1:6, key="Cell_type", value="Proportion")
allProport_DWLS_c1_bc$Group <- factor(allProport_DWLS_c1_bc$Group, levels=c("YMA", "AMA"))


# SVG works better for showing in Word and PowerPoint presentations 

c1_dwls_upd_bar_upd <- ggplot(allProport_DWLS_c1_bc, aes(fill=gsub("_", " ", Cell_type), y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "dodgerblue3", "#8B008B", "gray27"))+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 30, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.spacing.y = unit(0.5, "cm")) +
  labs(x = NULL, y = "Proportion of cells (%)", fill = NULL) +
  guides(fill = guide_legend(nrow=2))
c1_dwls_upd_bar_upd
ggsave(paste0(outputdir, "/EndoAgeing_deconv_C1_midsec_DWLS_raw_cellprop_upd.svg"), 
       device = "svg", height = 18, width = 24, units = "cm")
#dev.off()

# Same for 10x

allProport_DWLS_10x_bc <- gather(allProport_DWLS, 1:8, key="Cell_type", value="Proportion")
allProport_DWLS_10x_bc$Group <- factor(allProport_DWLS_10x_bc$Group, levels=c("YMA", "AMA"))

tenx_dwls_upd_bar_upd <- ggplot(allProport_DWLS_10x_bc, aes(fill=gsub("_", " ", Cell_type), y=Proportion, x=Group)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_manual(values=c("indianred3", "mediumseagreen", "mediumslateblue", "dodgerblue3", "green4", "#8B008B", "gray50", "gray35"))+
  scale_y_continuous(limits = c(0, 1)) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 30, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom", 
        legend.spacing.y = unit(0.5, "cm")) +
  labs(x = NULL, y = "Proportion of cells (%)", fill =  NULL)+
  guides(fill = guide_legend(nrow=3))
tenx_dwls_upd_bar_upd
ggsave(paste0(outputdir, "/EndoAgeing_deconv_tenx_midsec_DWLS_raw_cellprop_upd.svg"), 
       device = "svg", height = 18, width = 24, units = "cm")

# Grouped by age (e)

dwls_c1_grouped <- allProport_DWLS_C1 %>%  group_by(Age)
theme_set(theme_bw())
age_ord <- ggplot(dwls_c1_grouped, aes(x = Age, y = Ciliated, color = as.factor(Group))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(y = "Proportion of ciliated cells (%)", color = "Group") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_blank())
age_ord
ggsave(paste0(outputdir, "/C1_DWLS_prop_ageordered.pdf"), 
       device = "pdf", height = 8, width = 6, units = "cm")


