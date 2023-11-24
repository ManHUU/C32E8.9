# Parameters
n_clust <- 3
fastq_dir <- "./data/fastp_output/"
fasta_ref <- "./resources/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
gtf_ref <- "./resources/Caenorhabditis_elegans.WBcel235.105.gtf"
star_genome <- "./data/STAR_genome"

# raw file names
fq_files <- list.files(fastq_dir, full.names = T, recursive = T, pattern = "fq\\.gz$") 
names(fq_files) <- gsub(".fq.gz$", "", gsub(".*\\/", "", fq_files))
samples <- unique(gsub("_cleaned.*", "", names(fq_files)))

# ---------------------------------
# Read mapping using STAR
# ---------------------------------
# Track Progress
# check which samples are already done
done_samples <- c()
if(file.exists("./logs/STAR_Mapping.log")) {
  mv_tmp_log <- Log_data("mv ./data/STAR_output_tmp/(.*)__. ./data/STAR_output", 
                           Log_path = "./logs/STAR_Mapping.log")
  done_samples <- mv_tmp_log$Sample[mv_tmp_log$Exit.Status == 0]
}

# Map remaining samples
map <- parallel::mclapply(mc.cores = 2, 
                             X = samples[!samples %in% done_samples], 
                             FUN = function(sample){
                               align_Command <- paste0("./tools/STAR-2.7.10b/bin/MacOSX_x86_64/STAR ",
                          "--runThreadN ", n_clust,
                          " --genomeDir ", star_genome,
                          " --readFilesCommand gunzip -c",
                          " --outSAMtype BAM SortedByCoordinate",
                          " --quantMode TranscriptomeSAM GeneCounts",
                          " --readFilesIn ", paste(fq_files[grepl(sample, names(fq_files))], collapse = " "),
                          " --outFileNamePrefix", " ./data/STAR_output_tmp/", sample, "__")

  System(align_Command, Log_name = "./logs/STAR_Mapping.log", Append = TRUE)
  
  move_Command <- paste0("mv ./data/STAR_output_tmp/", sample, "__* ./data/STAR_output")
  System(move_Command, Log_name = "./logs/STAR_Mapping.log", Append = TRUE)
})

# Check Log for All Jobs
mapping_log <- Log_data(".*tmp/(.*-\\d)__.*", Log_path = "./logs/STAR_Mapping.log")

# ---------------------------------
# Index bam files
# ---------------------------------
bam_files <- list.files("./data/STAR_output/",
                        full.names = T, recursive = T,
                        pattern = ".*Coord.out.bam$") 

bam_index_done <- list.files("./data/STAR_output/",
                             full.names = T, recursive = T,
                             pattern = ".*Coord.out.bam.bai$") 


adding_index <- parallel::mclapply(mc.cores = 3,
                                   X = bam_files[!bam_files %in% bam_index_done],
                                   FUN = function(bam_file){
                                     index_command <- paste0("samtools ",
                                                             "index ",
                                                             bam_file)
                                     System(index_command, 
                                            Log_name = "./logs/bam_add_index.log", 
                                            Append = TRUE)
                                     })

index_info <- Log_data("samtools index ./data/STAR_output//.*__.*", 
                       Log_path = "./logs/bam_add_index.log")
mean(index_info$Exit.Status)


# ---------------------------------
# Read Counts
# ---------------------------------
# list the counts files output by STAR, which ended with ReadsPerGene.out.tab
count_files <- list.files("./transcriptomic/data/STAR_output/",
                          full.names = T, recursive = T,
                          pattern = ".*ReadsPerGene.out.tab$") 
file <- count_files[1]
counts <- mapply(file = count_files, FUN = function(file){
  # read counts table, remove the first 4 comments rows (mapping information),
  # select columns 1 (geneID), and 4 (reverse stranded counts) 
  # !!!!!!!!!!! be careful of which strand to select, it's based on the RNA-seq method, can use RseQC to infer the experiments information
  sample_counts <- as.data.frame(read.table(file, sep = "\t")[-c(1:4),c(1,4)])
  # give counts table rownames as geneID
  rownames(sample_counts) <- sample_counts[,1]
  # give counts table colnames as 'geneID' and Sample-ID
  colnames(sample_counts) <- c("geneID", gsub(".*//(.*)__ReadsPerGene.out.tab", "\\1",file))
  # return column 2 (unstranded counts)
  return(sample_counts[, 2, drop=FALSE])
}, SIMPLIFY = FALSE)

# cbind all tables in counts in one table
counts <- do.call(cbind, counts)

# save count matrix
saveRDS(counts, file = "./transcriptomic/data/R_data/counts.RDS")
