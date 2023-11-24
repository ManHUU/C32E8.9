#' @Author Iman M. Hu
#' @Date 2023.3.28
#' @Description Perform quality control for fastq files
#' @Last_modified 2023.3.28
 

## define global parameters 
# parameter the number of threads to use
n_clust <- 6 
# define fastq sample directory
fastq_dir <- "./data/raw/"
# Read Fastq File names
fq_files <- list.files(fastq_dir, full.names = T, recursive = T, pattern = "\\.fastq\\.gz$") 
names(fq_files) <- gsub(".fastq.gz$", "", gsub(".*\\/", "", fq_files))
samples <- unique(gsub("_L.*", "", names(fq_files)))


# ---------------------------------
# Perform QC for each sample
# ---------------------------------
# Create Directory for Quality Control
System("mkdir ./data/QC", Log_name = "./logs/QC.log", Append = FALSE)

# Track Progress (Needed bz of Large Sample Size)
done_samples <- list.dirs("./Data/QC/", recursive = F)
html <- list.files(done_samples, pattern = "\\.html$")
done_file <- gsub("_fastqc.html", "", html)
length(done_file)
length(fq_files)
remaining_files <- names(fq_files)[!names(fq_files) %in% done_file]
length(remaining_files)
remaining_files_path <- fq_files[remaining_files]

parallel::mclapply(mc.cores = n_clust, X = remaining_files_path, function(file_path){
  # get the name of the file_path
  file_name <- gsub(".fastq.gz$", "", gsub(".*\\/", "", file_path))
  # print which file is currently being processed
  print(paste0("Performing QC for sequencing reads of sample: ", file_name))
  # create output directory
  dir.create(paste0("./Data/QC/", file_name))
  # create the command to run fastqc of file, and output in specified folder named in file_name
  fq_command <- paste0("fastqc ", file_path," -o ./data/QC/", file_name)
  # run the fastqc command, and output to the log file (running information)
  System(fq_command, Log_name = "./logs/QC.log", Append = TRUE)
})


QC_Log_data <- Log_data("fastqc ./data/raw//(.*).fastq.gz -o.*", 
                        Log_path = "./logs/QC.log")
mean(QC_Log_data$Exit.Status)


# ---------------------------------
# MultiQC Report
# ---------------------------------
# Run command in terminal to generate multiQC aggregating all fastqc reports
# multiqc ./QC --filename 'MultiQC-report' --outdir ../reports 


# ---------------------------------
# Perform fastp clean-up 
# ---------------------------------
# Track Progress
# check which samples are already done
done_samples <- c()
if(file.exists("./logs/fastp_clean.log")) {
  fastp_log <- Log_data(".*15 -i ./data/raw/(.*)_L002_R1.fastq.gz.*", 
                          Log_path = "./logs/fastp_clean.log")
  done_samples <- fastp_log$Sample[fastp_log$Exit.Status == 0]
}

# Map remaining samples
fastp <- parallel::mclapply(mc.cores = 3, 
                          X = samples[!samples %in% done_samples], 
                          FUN = function(sample){
                            fastp_Command <- paste0("fastp ",
                                                    "--qualified_quality_phred 15",
                                                    " --unqualified_percent_limit 4", 
                                                    " --n_base_limit 5",
                                                    " --length_required 15",
                                                    " -i ", "./data/raw/", sample, "_L002_R1.fastq.gz",
                                                    " -I ", "./data/raw/", sample, "_L002_R2.fastq.gz",
                                                    " -o ", "./data/fastp_output/", sample, "_cleaned_R1.fq.gz",
                                                    " -O ", "./data/fastp_output/", sample, "_cleaned_R2.fq.gz",
                                                    " --html ", "./reports/fastp_reports/", sample, "_fastp_report.html",
                                                    " --json ", "./reports/fastp_reports/", sample, "_fastp_report.json")
                            
                            System(fastp_Command, Log_name = "./logs/fastp_clean.log", Append = TRUE)
                          })

# Check Log for All Jobs
fastp_log <- Log_data(".*15 -i ./data/raw/(.*)_L002_R1.fastq.gz.*", Log_path = "./logs/fastp_clean.log")

