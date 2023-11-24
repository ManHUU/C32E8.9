# use edgeR to get cpm
library(edgeR)

# read counts
counts <- readRDS("./data/R_data/counts.RDS")
samples <- colnames(counts)

# generate DGElist
dge_list <- edgeR::DGEList(counts = counts)
# add group information to DGElist
dge_list$samples$group <- gsub("(.*)-.*", "\\1", samples)
# add design matrix to DGElist
group_ <- dge_list$samples$group
dge_list$design <- model.matrix(~ 0 + group_)

# filtering low expressed genes
# get rid of genes which did not occur frequently enough, using default cut-offs
genes_keep <- edgeR::filterByExpr(dge_list, design = dge_list$design)
dge_list$counts <- dge_list$counts[genes_keep, ]

# compare edger results calcNormfactors before/after filter
# after filter
dge_list <- edgeR::calcNormFactors(dge_list)
# get cpm normalized data, for visualization/other analysis
cpm <- edgeR::cpm(dge_list)
saveRDS(cpm, file = "./data/R_data/cpm.RDS")
