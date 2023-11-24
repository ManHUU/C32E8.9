library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# -----------------------------
# lipids class enrichment analysis
# -----------------------------

# load DE lipids file
df_data <- readRDS("./data/R_data/DE_lipids_limma_voom.RDS")

# Lipid Class Custom annotation sets
# Add column for Lipid Class to each top table 
# add lipid class
df_data$Lipid_Class <- gsub("(.*)\\(.*\\)", "\\1", rownames(df_data))
# add lipid identification
df_data$indentification <- rownames(df_data)
TERM2GENE <- df_data[,c("Lipid_Class", "indentification")]

# Prepare Background genes
lipids.background <- df_data$indentification

# prepare different accumulated lipids set
ajust_p <- df_data[,grepl("(.*)_adjustP", colnames(df_data))]
groups_adjust_p <- colnames(ajust_p)

log_fc <- df_data[,grepl("(.*)_logFC", colnames(df_data))]
groups_log_fc <- colnames(log_fc)

df_lipids <- mapply(x = groups_log_fc,
                    y = groups_adjust_p,
                    FUN = function(x, y){
                      group <- gsub("(.*)_logFC", "\\1", x)
                      df_gene_up <- rownames(df_data[df_data[, x] >= 0.5 & 
                                                       df_data[, y] <= 0.05 & 
                                                       !is.na(df_data[, y]), c(x, y)])
                      df_gene_down <- rownames(df_data[df_data[, x] <= -0.5 & 
                                                         df_data[, y] <= 0.05 & 
                                                         !is.na(df_data[, y]), c(x, y)])
                      group <- list("up" = df_gene_up, "down" = df_gene_down)
                      
                    }, SIMPLIFY = FALSE)

df_lipids <- unlist(df_lipids, recursive = F)
names(df_lipids) <- gsub("(logFC.)", "", names(df_lipids))


# -----------------------------
# # perform the lipids enrichment (over-representation)
# -----------------------------
lipids_groups <- names(df_lipids)[!names(df_lipids) %in% c("MRPS5_vs_HT115_up", "MRPS5_vs_HT115_down")]

mapply(group = lipids_groups, FUN = function(group){
  lipids <- df_lipids[[group]]
  enrich <- clusterProfiler::enricher(gene = lipids,
                                      TERM2GENE = TERM2GENE,
                                      universe = lipids.background,
                                      maxGSSize = 800)
  
  # define file name all following files
  file_name <- paste0("OR_enrich","_", group)
  # write results table and save in reports
  write.csv(enrich@result, file = paste0("./reports/over-representation/", file_name, ".csv"))
  
  # make plots
  if(sum(enrich@result$qvalue < 0.05) > 0) {
    pdf(file = paste0("./plots/over_representation/dotplot/", file_name, "_dotplot",".pdf"))
    print(enrichplot::dotplot(enrich, showCategory = 30))
    dev.off()}
  
  }, SIMPLIFY = FALSE)

group <- "MRPS5.ECHDC1_vs_MRPS5_down"
mapply(group = lipids_groups, FUN = function(group){
  lipids <- df_lipids[[group]]
  enrich <- clusterProfiler::enricher(gene = lipids,
                                      TERM2GENE = TERM2GENE,
                                      universe = lipids.background,
                                      maxGSSize = 800)
  
  # define file name all following files
  file_name <- paste0("OR_enrich","_", group)
  # write results table and save in reports
  write.csv(enrich@result, file = paste0("./reports/over-representation/", file_name, ".csv"))
  
  # make plots
  if(sum(enrich@result$qvalue < 0.05) > 0) {
    # make plots
    result <- enrich@result[enrich@result$p.adjust < 0.2, ]
    plot_data1 <- result[order(result$pvalue, decreasing = FALSE),]
    plot_data1$Description <- factor(plot_data1$Description, levels = rev(plot_data1$Description))
    
    pdf(file = paste0("./plots/over_representation/dotplot/", file_name,".pdf"))
    print(ggplot(data = plot_data1, mapping = aes(x = -log10(pvalue), y = Description)) +
            geom_point(mapping = aes(size = Count ,color = -log10(p.adjust))) +
            scale_color_gradient(low = "grey", high ="deepskyblue4", breaks = c(4,8,12)) +
            labs(y = "", title = paste0("gsea ", group)) + 
            theme(axis.text = element_text(size = 12), 
                  plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
            theme_bw())
    dev.off()
    
    }
  
}, SIMPLIFY = FALSE)


# -----------------------------
# # perform the lipids enrichment (gsea)
# -----------------------------
log_fc <- df_data[,grepl("(.*)_logFC", colnames(df_data))]
groups_log_fc <- colnames(log_fc)


group <- groups_log_fc[1]


mapply(group = groups_log_fc, FUN = function(group){
  lipids_group <- df_data[,group]
  # names input
  names(lipids_group) <- rownames(df_data)
  # order in desending
  lipids_group <- lipids_group[order(-lipids_group)]

  enrich <- clusterProfiler::GSEA(geneList = lipids_group,
                                  TERM2GENE = TERM2GENE,
                                  minGSSize = 4, 
                                  maxGSSize = 1000,
                                  pvalueCutoff = 1, 
                                  verbose = F, 
                                  seed = T)
  
  # define file name all following files
  group_info <- gsub("(.*)_logFC", "\\1", group)
  file_name <- paste0("gsea_enrich","_", group_info)
  # write results table and save in reports
  write.csv(enrich@result, file = paste0("./reports/gsea/", file_name, ".csv"))
  
  # make plots
  result <- enrich@result
  plot_data1 <- result[order(result$NES, decreasing = TRUE),]
  plot_data1$Description <- factor(plot_data1$Description, levels = rev(plot_data1$Description))
  
  colnames(plot_data)
  pdf(file = paste0("./plots/gsea/dotplot/", group_info,".pdf"))
  print(ggplot(data = plot_data1, mapping = aes(x = NES, y = Description)) +
          geom_point(mapping = aes(size = setSize,color = -log10(p.adjust))) +
          scale_color_gradient(low = "grey", high ="deepskyblue4", breaks = c(1.3,2,4,6)) +
          labs(y = "", title = paste0("gsea ", group)) + 
          theme(axis.text = element_text(size = 12), 
                plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
          theme_bw())
  dev.off()
  }, SIMPLIFY = FALSE)


