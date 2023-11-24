# ---------------------------------
# LOG function: Custom System Command to capture: Exit Status, Output & Errors
# ---------------------------------
#' Log the detailed information of executing a shell command. 
#' 
#' @param Command  string: shell command to give to the system 
#' @param Log_name string: the file path and name of the log file to which to write the information of the command
#' @param Append boolean: if true append to log file, otherwises overwrite
#' @return List object. Including the running system time, exit status, standard out/error. 
#' Function will write the information to a log file "Log_name"
#' @examples
#' System("mkdir ./Data/QC", Log_name = "./Logs/QC.log", Append = FALSE)

System <- function(Command, Log_name, Append){
  # Create temporary files to hold command outputs
  ## create 3 names of temporary files for the script we wanna run, and 
  ## store the output of the Command
  stderrFile <- tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
  stdoutFile <- tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
  scriptFile <- tempfile(pattern="R_robust.system_script", fileext=as.character(Sys.getpid()))
  write(paste0(Command, " 2> ", shQuote(stderrFile), " 1> ", shQuote(stdoutFile)), file = scriptFile)

  # Create a list to hold output
  log_info <- list()
  # run command capturing error messages
  log_info$StartSysTime <- Sys.time()
  log_info$exitStatus = system(paste0("sh ", scriptFile))
  log_info$EndSysTime <- Sys.time()
  log_info$stdout = readLines(stdoutFile)
  log_info$stderr = readLines(stderrFile)
  
  # remove temporary files
  unlink(c(stdoutFile, stderrFile, scriptFile))

  # Write to Log File
  write(c(format(log_info$StartSysTime, "%a %b %d %X %Y"), "Command:"), 
        file = Log_name, sep = "\n\n", append = Append)
  write(Command, file = Log_name, sep = "\n\n", append = TRUE)
  write(c("\n", "Run_Time:"), file = Log_name, sep = "\n\n", append = TRUE)
  write(paste0(difftime(log_info$EndSysTime, log_info$StartSysTime, units = "secs")[[1]], ' Secs'), 
        file = Log_name, sep = "\n\n", append = TRUE)
  write(c("\n", "Exit_Status:"), file = Log_name, sep = "\n\n", append = TRUE)
  write(log_info$exitStatus, file = Log_name, sep = "\n\n", append = TRUE)
  write(c("\n", "Std_Out:"), file = Log_name, sep = "\n\n", append = TRUE)
  write(log_info$stdout, file = Log_name, sep = "\n\n", append = TRUE)
  write(c("\n", "Std_Err:"), file = Log_name, sep = "\n\n",append = TRUE)
  write(log_info$stderr, file = Log_name, sep = "\n\n", append = TRUE)
  return(log_info)
}


# ---------------------------------
# Custom Function to retrieve Log File Information
# ---------------------------------
#' Retrive information from Log File. 
#' 
#' @param Sample_regex  string: The regexp of the pattern we want to search to modify the command 
#' @param Log_path string: The path of log file to retrive information from
#' @return data.frame. Including the command, exit status, running time
#' @examples
#' Log_data('fastqc ./data/raw//(.+).fastq.gz -o.+', Log_path = "./logs/QC.log")
Log_data <- function(Sample_regex, Log_path){
  # get information from specified log file and store it in a dataframe
  Command_Stats <- cbind(data.frame(gsub(Sample_regex, "\\1", 
                                         system(paste0("awk 'f{print;f=0} /Command:/{f=1}' ", 
                                                       Log_path), intern = TRUE))),
                         # get the Exit_Status from the log file
                         data.frame(system(paste0("awk 'f{print;f=0} /Exit_Status/{f=1}' ", 
                                                  Log_path), intern = TRUE)),
                         # get the Run_Time of each command from the log file
                         data.frame(system(paste0("awk 'f{print;f=0} /Run_Time/{f=1}' ", 
                                                  Log_path, " | sed 's/ Secs//g'"), intern = TRUE)))
  # rename the columns of the dataframe
  colnames(Command_Stats) <- c("Sample", "Exit.Status", "Run.Time")
  # convert exit.status to numeric
  Command_Stats$Exit.Status <- as.numeric(as.character(Command_Stats$Exit.Status))
  # convert run time to numeric
  Command_Stats$Run.Time <- as.numeric(as.character(Command_Stats$Run.Time))
  return(Command_Stats)
}


# ---------------------------------
# Custom Function to retrieve STAR mapping information
# ---------------------------------
#' Retrive information from STAR Log.final.out 
#' 
#' @param STAR_output_Path  string: path specifying the directory containing the STAR outputs
#' @return data.frame. Including the sample, infor_2_extract (Uniquely_mapped, Mapped_multiple_loci...)
#' @examples
#' Mapping_data(STAR_output_Path = "./data/STAR_output/")

Mapping_data <- function(STAR_output_Path){
  # specify the information to extract from STAR Log.final.out 
  infor_2_extract <- list("Uniquely_mapped" = " *Uniquely mapped reads number |\t/",
                          "Mapped_multiple_loci" = " *Number of reads mapped to multiple loci |\t/",
                          "Mapped_too_many_loci" = " *Number of reads mapped to too many loci |\t/",
                          "Unmapped_too_short" = " *Number of reads unmapped: too short |\t/",
                          "Unmapped_other" = " *Number of reads unmapped: other |\t/", 
                          "Avg_input_length" = " *Average input read length |\t/",
                          "Avg_mapped_length" = " *Average mapped length |\t/",
                          "Uniquely_mapped_%" = " *Uniquely mapped reads % |\t/",
                          "Multi_mapped_%" = " *% of reads mapped to multiple loci |\t/")
  # get the file name in the output directory
  Log_file_names <-list.files(STAR_output_Path, "Log.final.out$", recursive = FALSE)
  # extract sample name and infor_2_extract 
  map_infor_all <- mapply(log_file_name = Log_file_names, FUN = function(log_file_name){
    # get file name, then extract sample name
    sample_name <- data.frame("Sample" = gsub("(.*)__Log.final.out", "\\1", log_file_name))
    # use 'sed' get unique mapped reads
    sample_info <- mapply(infor = infor_2_extract, FUN = function(infor){
      system(paste0("sed -n ",
                    "'s/",infor,
                    "/","p'", " ./data/STAR_output/",log_file_name), intern = TRUE)}, 
      SIMPLIFY = FALSE)
    sample_info <- do.call(cbind, sample_info)
    map_infor <- cbind(sample_name, sample_info)
    return(map_infor)
  }, SIMPLIFY = FALSE)
  map_infor_all <- do.call(rbind, map_infor_all)
  return(map_infor_all)
}

  
# ---------------------------------
# Custom Function to retrieve fastp cleaning and quality check  information
# ---------------------------------
#' Retrieve information from fastp_reports (output reports of fastp cleaning and quality checking)
#' @param fastp_reports_path  string: path specifying the directory containing the fastp report outputs
#' @return data.frame. Including the sample, infor_2_extract (reads Before'after filtering...)
#' @examples
#' Extract_Fastp_info(Path = "./reports/fastp_reports/")


Extract_Fastp_info <- function(fastp_reports_path){
  # get json files name in specified directory
  file_name <- list.files(fastp_reports_path, "json$", recursive = FALSE)
  ## specify information to extract
  fastp_info_2_extract <- list("Before_filter" = "\"before_filtering",
                               "After_filter" = "\"after_filtering",
                               "insert_size" = "\"insert_size",
                               "duplication" = "\"duplication"
  )
  fastp_info <- mapply(file = file_name, FUN = function(file){
    # get sample name
    sample_name <- gsub("(.*)_fastp_report.json", "\\1", file)
    # get before/after filtering reads
    ### get reads before/after filtering (use awk to get the first line after specified pattern (awk 'f{print;f=0} /pattern/{f=1}')
    fastp_reads <- mapply(info = fastp_info_2_extract, FUN = function(info){
      # use system awk command to search get the line which contains specified information
      total_reads <- system(paste0("awk 'f{print;f=0} /", info, 
                                   "/{f=1}' ", 
                                   "./reports/fastp_reports/",
                                   file), intern = TRUE)
      # use gsub to get information 
      total_reads <- gsub(".*:(.*),", "\\1", total_reads)
    }, SIMPLIFY = FALSE)
    
    fastp_reads <- as.data.frame(do.call(cbind, fastp_reads))
    fastp_reads$insert_size <- gsub(" (\\d*)", "\\1", fastp_reads$insert_size)
    fastp_reads$duplication <- gsub(".*: (\\d*)", "\\1", fastp_reads$duplication)
    
    # get "low_quality_reads" 
    low_quality_reads <- system(paste0("awk 'c&&!--c;/",
                                       ".*\"filtering_result",
                                       "/{c=2}' ", 
                                       "./reports/fastp_reports/", file), intern = TRUE)
    fastp_reads$low_quality_reads <-  gsub(".*\":.(.*),", "\\1", low_quality_reads)
    # "too_many_N_reads"
    too_many_N_reads <- system(paste0("awk 'c&&!--c;/",
                                      ".*\"filtering_result",
                                      "/{c=3}' ", 
                                      "./reports/fastp_reports/", file), intern = TRUE)
    fastp_reads$too_many_N_reads <- gsub(".*\":.(.*),", "\\1", too_many_N_reads)
    # "too_short_reads"
    too_short_reads <- system(paste0("awk 'c&&!--c;/",
                                     ".*\"filtering_result",
                                     "/{c=4}' ", 
                                     "./reports/fastp_reports/", file), intern = TRUE)
    fastp_reads$too_short_reads <- gsub(".*\":.(.*),", "\\1", too_short_reads)
    
    sample_read <- cbind("Sample" = sample_name, fastp_reads)
  }, SIMPLIFY = FALSE)
  
  fastp_all <- do.call(rbind, fastp_info)
  return(fastp_all)
}



