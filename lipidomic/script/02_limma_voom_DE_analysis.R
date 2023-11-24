library(limma)
library(edgeR)

# use limma voom to calculate DE lipids
# read DE lipids data
lipid_df <- readxl::read_xlsx("./data/lipids_abundance.xlsx")
# remove NA rows
lipid_df <- as.data.frame(na.omit(lipid_df))

# give row names and rename colnames
rownames(lipid_df) <- lipid_df$Metabolite
lipid_mtx <- lipid_df[,colnames(lipid_df)[!colnames(lipid_df) 
                                          %in% c("Metabolite", "VIP_score")]]
colnames(lipid_mtx) <- gsub("(.*)_Group_.*", "\\1", colnames(lipid_mtx))

# use edgeR to normalize aboundance based on library size of sample
# create DGEList object
exp <- edgeR::DGEList(counts = lipid_mtx)
# add group information to DGElist
exp$samples$group <- gsub("(.*)_.*.","\\1",colnames(lipid_mtx))
# add design matrix to DGElist
group_ <- exp$samples$group
exp$design <- model.matrix(~ 0 + group_)

# use voom check data and prepare input data for fitting model
voom_d <- voom(counts = exp,  # put DEGlist object here, results are different when using counts and DEGlist object
               design = exp$design,
               # check if there are too many low expressed with hight variance genes (if needs filtering again)
               plot = TRUE)

# Fitting linear models in limma
fit_limma <- lmFit(object = voom_d,
                   design = exp$design)

# Specify which groups to compare, and get DE genes table
de_lipids_limma <- mapply(condition1 = c("ECHDC1", "MRPS5", "MRPS5.ECHDC1", "MRPS5.ECHDC1"),
                          condition2 = c(rep("HT115", 3), "MRPS5"), FUN = function(condition1, condition2) {
                            print(condition1)
                            print(condition2)
                            coparison <- paste0("group_", condition1, " - ", "group_", condition2)
                            contr <- makeContrasts(contrasts = coparison, levels = colnames(coef(fit_limma)))
                            tmp <- contrasts.fit(fit_limma, contr)
                            # calcute p value
                            tmp <- eBayes(tmp)
                            # get differential expressed (DE) genes table
                            top_table <- topTable(tmp, sort.by = "none", number=Inf)
                            colnames(top_table) <- gsub("logFC", "logFC", colnames(top_table))
                            colnames(top_table) <- gsub(" P.Value", "Pvalue", colnames(top_table))
                            colnames(top_table) <- gsub("adj.P.Val", "adjustP", colnames(top_table))
                            # add comparison information to colnames of DE genes table
                            colnames(top_table) <- paste0(condition1, "_vs_", condition2, "_", colnames(top_table))
                            # return the differential expressed (DE) genes table
                            return(top_table)
                            }, SIMPLIFY = FALSE)

# cbind all DE genes table together and save as RDS
de_lipids_limma <- do.call(cbind, unname(de_lipids_limma))
saveRDS(de_lipids_limma, file = "./data/R_data/DE_lipids_limma_voom.RDS")

