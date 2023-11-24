# prepare plotting data
## x is length of lipids
## y is the seturation of lipids
## color is the logFC of lipids
## size is the P.value of lipids

library(ggplot2)
library(stringr)
# read Defferential lipids dataframe
de_lipid <- readRDS("./DE_lipids_limma_voom.RDS")
# select interesting comparasion
de_lipid <- de_lipid[,c("MRPS5.ECHDC1_vs_MRPS5_logFC", "MRPS5.ECHDC1_vs_MRPS5_adjustP")]
# add class, length of lipids and seturation to dataframe
de_lipid$length <- gsub(".*\\((.*):.*\\)", "\\1", rownames(de_lipid))
de_lipid$class <- gsub("(.*)\\(.*\\)", "\\1", rownames(de_lipid))

de_lipid$class[grepl(".*-.*", de_lipid$length)] <- gsub("(.*\\(.*)-.*\\)", "\\1", rownames(de_lipid))[grepl(".*-.*", de_lipid$length)]
de_lipid$class <- str_replace(de_lipid$class, "\\(", " ")

de_lipid$length[grepl(".*-.*", de_lipid$length)] <- gsub(".*-(.*)", "\\1", de_lipid$length[grepl(".*-.*", de_lipid$length)])
de_lipid$length[grepl("(d|t).*", de_lipid$length)] <- gsub("(d|t)(.*)", "\\2", de_lipid$length[grepl("(d|t).*", de_lipid$length)])
de_lipid$length <- as.numeric(de_lipid$lengt)

de_lipid$seturation <- gsub(".*:(.*)\\)", "\\1", rownames(de_lipid))
de_lipid$seturation <- as.numeric(de_lipid$seturation)

# -----------------------------
# make seteuration scatter plot
# -----------------------------
# change column names
file_dir <- paste0("/Users/manhu/Desktop/Projects/echdc_1/lipidomic/plots/scatter_plot/double_check/")

# select interested lipid classes
classes <- c("TG", "Cer", "2-acyl LPC", "PG", "PC", "PE")

# loop over lipid classes and make plots
mapply(class = classes, FUN = function(class){
  print(class)
  mapply(group = c( "all", "sub"), FUN = function(group){
    # expression of lipids class
    # if group = all, use all lipids, else use filtered lipids
    # if(group == "all") {
    #   data_plot <- de_lipid[de_lipid$class == class,]
    # } else {
    #   data_plot <- de_lipid[de_lipid$class == class & abs(de_lipid$MRPS5.ECHDC1_vs_MRPS5_logFC) > 0.5, ]
    # }
    data_plot <- de_lipid[de_lipid$class == class,]
    if(group == "sub") {
      data_plot <- data_plot[abs(data_plot$MRPS5.ECHDC1_vs_MRPS5_logFC) > 0.5  & data_plot$MRPS5.ECHDC1_vs_MRPS5_adjustP < 0.05, ]
    }
    
    colnames(data_plot)[grepl("_logFC", colnames(data_plot))] <- "logFC"
    colnames(data_plot)[grepl("adjustP", colnames(data_plot))] <- "adjustP"
    
    pdf(paste0(file_dir, class, group, "_seturation.pdf"), height = 3.5, width = 5)
    print(ggplot(data = data_plot, mapping = aes(x = length, y = seturation)) + 
            geom_point(aes(fill = logFC, size = -log10(adjustP)), colour = "dimgrey", shape = 21) +
            scale_fill_gradient2(low= "#435040",mid = "white", high= "#78343c") +
            scale_x_continuous(breaks = scales::pretty_breaks())+
            theme_bw())
    dev.off()
  })
  })


