# Use DESeq2_1.32.0 downloaded 2023.5.1 to do differential analysis
# You need two things: cout and sample_information (in same order)

# Read count matrix generated from STAR
count <- readRDS("./data/R_data/counts.RDS")
# Create sample information table, which we will name coldata,
# note: make sure  count matrix
# and column data, they are consistent in terms of sample order
coldata <- data.frame("condition" = as.factor(gsub("(.*)-.*", "\\1",
                                                   colnames(count))))
rownames(coldata) <- colnames(count)
# Construct DESEQDataSet Object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~ condition)

# Differential expression analysis, The standard differential expression
# analysis steps are wrapped into a single function, DESeq.
dds <- DESeq2::DESeq(dds)
saveRDS(dds, "./data/DE_analysis/DeSeq2_dds.RDS")
# lists the coefficients
DESeq2::resultsNames(dds)


dds <- readRDS("./data/DE_analysis/DeSeq2_dds.RDS")
# Generate result tables, which extracts a results table with log2 fold changes,
# p values and adjusted p values.
# Method1 to get results for all comparisons you are interested in
res_all <- mapply(x = c("ech", "mrp", "MC", "MC"),
       y = c(rep("ctr", 3), "mrp"),
       FUN = function(x, y) {
         res <- DESeq2::results(dds, contrast = c("condition", x, y))
         colnames(res) <- gsub("log2FoldChange", "logFC", colnames(res))
         colnames(res) <- gsub("pvalue", "Pvalue", colnames(res))
         colnames(res) <- gsub("padj", "adjustP", colnames(res))
         colnames(res) <- paste0(x, "_vs_", y, "_", colnames(res))
         return(res)
         }, SIMPLIFY = FALSE)

res_all_1 <- as.data.frame(do.call(cbind, unname(res_all)))
saveRDS(res_all_1, "./data/DE_analysis/DE_genes_DESeq.RDS")


# Data transformations and visualization
# Count data transformations
# Extracting transformed values, use vst here, and the blind = FALSE, as
# we expect differences in counts which are explainable by the
# experimental design, and want use transformed data for downstream analysis
## 3 methods for data transformation
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "./data/R_data/DeSeq2_vsd.RDS")
vsd <- readRDS("./data/R_data/DeSeq2_vsd.RDS")
# check cst transformed data
head(assay(vsd))
dim(assay(vsd))
# make plots
plotPCA(vsd, intgroup = c("condition")) +
  geom_text(label = colnames(assay(vsd)),
            nudge_x = 0.5,
            nudge_y = 1,
            check_overlap = TRUE)



