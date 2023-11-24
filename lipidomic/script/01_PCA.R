library(FactoMineR)
library(ggplot2)

# load data
# read lipids concentration data (value = retention time of lipid/retention time of internalstandard)
lipid_df <- readxl::read_xlsx("./data/lipids_abundance.xlsx")
# remove NA rows
lipid_df <- as.data.frame(na.omit(lipid_df))
# give rownames and rename colnames
rownames(lipid_df) <- lipid_df$Metabolite
lipid_mtx <- lipid_df[,colnames(lipid_df)[!colnames(lipid_df) 
                                  %in% c("Metabolite", "VIP_score")]]
colnames(lipid_mtx) <- gsub("(.*)_Group_.*", "\\1", colnames(lipid_mtx))

# ---------------------------------
# boxplot to check if data is normalized between different groups
# ---------------------------------
lipid_mtx_melt <- reshape2::melt(lipid_mtx)
colnames(lipid_mtx_melt) <- c("variable", "value")
# add group information
lipid_mtx_melt$group <- gsub("(.*)_.*", "\\1", lipid_mtx_melt$variable)

# boxplot
pdf(file = "plots/expression_boxplots.pdf", height = 5, width = 6)
ggplot(lipid_mtx_melt, aes(x = variable, y = log10(value), fill = group)) +
  geom_boxplot(color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
dev.off()

# ---------------------------------
################### PCA analysis  
# ---------------------------------
# analysis main components
lipid_mtx_t <- as.data.frame(t(lipid_mtx))
normalized_counts_pca <- FactoMineR::PCA(lipid_mtx_t, graph = FALSE)

# check main components
normalized_counts_pca$ind$coord
dim(normalized_counts_pca$var$coord)
summary(normalized_counts_pca)

# extract the information we want from our ‘normalized_counts_pca’ object
pca_df <- data.frame("pca1" = normalized_counts_pca$ind$coord[, 1], 
                     "pca2" = normalized_counts_pca$ind$coord[, 2])
pca_df$group <- gsub("(.*)_.*.","\\1",rownames(pca_df))
pca_df$group <- factor(pca_df$group, 
                       levels = c("HT115", "ECHDC1", "MRPS5", "MRPS5.ECHDC1"))

# make plots
pdf(file = "./plots/lipidomics_PCA_without_legend.pdf", height = 4.5, width = 4.5)
ggplot(data = pca_df, aes(x = pca1, y = pca2)) +
  geom_point(aes(color = group), size = 6) + 
  stat_ellipse(geom="polygon", 
               aes(fill = group, color = group), 
               alpha = 0.2, 
               show.legend = FALSE, 
               level = 0.95)  +
  scale_color_manual(values=c("#606060", "#077e97", "#a00000", "#c06000"))+
  xlab("PC1 (46.07%)") + 
  ylab("PC1 (16.50%)") +
  theme_bw() +
  theme(legend.position = "none") 
dev.off()


