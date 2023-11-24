# enrichment analysis
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Ce.eg.db)

# ----------------------------------------------------
# Prepare data for plotting
# ----------------------------------------------------
df_data <- readRDS("./transcriptomic/data/DE_analysis/DE_genes_DESeq_with_geneName.RDS")
groups_adjust_p <- colnames(df_data[,grepl("(.*)_adjustP", colnames(df_data))])
groups_log_fc <- colnames(df_data[,grepl("(.*)_logFC", colnames(df_data))])

df_genes <- mapply(log_fc = groups_log_fc,
                   adjust_p = groups_adjust_p, FUN = function(log_fc, adjust_p){
                     group <- gsub("(.*)_logFC", "\\1", log_fc)
                     df_gene_up <- df_data[rownames(df_data[df_data[, log_fc] >= 0.5
                                                            & df_data[, adjust_p] <= 0.05 
                                                            & !is.na(df_data[, adjust_p]),
                                                            c(log_fc, adjust_p)]), 'gene_name']
                     df_gene_down <- df_data[rownames(df_data[df_data[log_fc] <= -0.5 & 
                                                                df_data[, adjust_p] <= 0.05 & 
                                                                !is.na(df_data[, adjust_p]), 
                                                              c(log_fc, adjust_p)]), 'gene_name']
                     group <- list("up" = df_gene_up, "down" = df_gene_down)
                   }, SIMPLIFY = FALSE)

df_genes <- unlist(df_genes, recursive = F)
names(df_genes) <- gsub("(logFC.)", "", names(df_genes))


# ----------------------------------------------------
# Barplot of number of up and down genes per condition
# ----------------------------------------------------
DE.stat.df <-
  data.frame(count = unlist(lapply(df_genes, length)),
             comp = names(df_genes),
             group = gsub("^(.*_vs_.*)_.*$", "\\1", names(df_genes)),
             direction = gsub("^.*_vs_.*_(.*)$", "\\1", names(df_genes)))

# keep only comparisons of interest
DE.stat.df <- DE.stat.df[!DE.stat.df$group %in% 'MC_vs_ctr',]

# multiply down ones with -1 for plotting
DE.stat.df$count[DE.stat.df$direction == "down"] <- DE.stat.df$count[DE.stat.df$direction == "down"] * -1

# order the comparisons
DE.stat.df$group <- factor(DE.stat.df$group, levels = c("ech_vs_ctr", "mrp_vs_ctr", "MC_vs_mrp"))
DE.stat.barplot <-
  ggplot(DE.stat.df, aes(x = group, y = count, fill = direction)) +
  geom_bar(stat = 'identity', color = "black") +
  scale_fill_manual(values = c('up' = 'red',
                               'down' = 'blue')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

pdf(file = paste0("./plots/transcriptomic/DeSeq2/barplot/DE.count.barplot.pdf"),
    width = 3, height = 2.5)
print(DE.stat.barplot)
dev.off()


# ----------------------------------------------------
# Over-representation enrichment
# ----------------------------------------------------
# anotation <- c("BP", "CC", "MF")

anotation <- c("BP")
for(group_info in names(df_genes)){
  print(group_info)
  # input genelist for enrichGO
  group_genes <- df_genes[[group_info]]
  # loop enrichGO on "BP", "CC", "MF"
  for(ano in anotation) {
    print(ano)
    go <- enrichGO(group_genes,
                   ont = ano,
                   OrgDb = org.Ce.eg.db,
                   keyType = "SYMBOL",
                   pvalueCutoff = 1)
    
    file_name <- paste0("enrichGO","_", group_info, "_", ano)
    
    write.csv(go@result, file = paste0("./reports/Over_Representation_enrich_GO/", file_name, ".csv"))
    # remove not significant enriched terms and make plots
    if(sum(go@result$p.adjust < 0.05) > 0){
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_GO/dotplot/", file_name, "_", "dotplot", ".pdf"),
          width = 10, height = 10)
      print(dotplot(go, showCategory = 30) +
              ggtitle(paste0("dotplot", "_", file_name)))
      dev.off()
      
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_GO/upsetR/", file_name, "_", "upsetR", ".pdf"),
          width = 10, height = 10)
      print(upsetplot(go) +
              ggtitle(paste0("upsetR ", "_", file_name)))
      dev.off()
      ### plot network ### 
      go <- pairwise_termsim(go)
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_GO/enrichment_map/", file_name, "_",
                        "enrichment_map", ".pdf"), width = 10, height = 10)
      print(emapplot(go, layout = 'circle') +
              ggtitle(paste0("enrichment_map", "_", file_name)))
      
      dev.off() 
    }
  }
}


# ----------------------------------------------------
# EnrichKEGG
# ----------------------------------------------------
group_info <- names(df_genes)

mapply(group_info = names(df_genes), FUN = function(group_info){
  group_genes <- df_genes[[group_info]]
  # convert geneID to ENTRIZID
  gene_list <- bitr(group_genes, from  = "SYMBOL", 
                    to = c("ENTREZID", "WORMBASE","UNIPROT"),
                    OrgDb = org.Ce.eg.db)
  # genelist input for enrichKEGG
  gene_list_input <- unique(gene_list[, "ENTREZID"])
  kegg <- clusterProfiler::enrichKEGG(gene_list_input, 
                                      organism = 'cel', 
                                      keyType = 'ncbi-geneid', 
                                      pvalueCutoff = 1, 
                                      # package version: clusterProfiler_4.0.5, KEGG.db_2.8.0
                                      use_internal_data = TRUE)
  # remove no gene mapped comparison
  if(is.null(kegg)){
    print("No gene can be mapped")
  } else{
    # define file name for all plots and csv.file
    file_name <- paste0("enrichKEGG","_", group_info)
    write.csv(kegg@result, file = paste0("./reports/Over_Representation_enrich_KEGG/", 
                                         file_name, ".csv"))
    # remove not significant enriched terms and make plots
    if(sum(kegg@result$p.adjust < 0.05) > 0){
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_KEGG/dotplot/", file_name, "_", "dotplot", ".pdf"),
          width = 10, height = 10)
      print(dotplot(kegg, showCategory = 30) +
              ggtitle(paste0("dotplot", "_", file_name)))
      dev.off()
      
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_KEGG/upsetR/", file_name, "_", "upsetR", ".pdf"),
          width = 10, height = 10)
      print(upsetplot(kegg) +
              ggtitle(paste0("upsetR ", "_", file_name)))
      dev.off()
      
      ### plot network ### 
      kegg <- pairwise_termsim(kegg)
      pdf(file = paste0("./plots/transcriptomic/DeSeq2/enrich_KEGG/enrichment_map/", file_name, "_",
                        "enrichment_map", ".pdf"), width = 10, height = 10)
      print(emapplot(kegg) +
              ggtitle(paste0("enrichment_map", "_", file_name)))
      
      dev.off() 
    }
  }
}, SIMPLIFY = FALSE)


